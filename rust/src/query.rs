//! Memory-mapped YACHT run index and the two query phases used by `yacht run`.
//!
//! Phase A (`query`): look up every sample hash; record which rows / shared sets were hit.
//! Phase B (`exclusive_counts`): given the set N of "nontrivial" rows (rows sharing at least one hash with
//! the sample, extended by organism name as the legacy multisearch-based code does), count for every row
//! g in N the hashes that occur in g and in no other row of N (`excl`), and how many of those are in the
//! sample (`matches`). This is exactly what `hypothesis_recovery_src.get_exclusive_hashes` computes.

use crate::format::Header;
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

pub struct Index {
    mmap: Mmap,
    pub hdr: Header,
}

/// Result of phase A.
#[derive(Debug, Clone, Default)]
pub struct QueryResult {
    /// (row, number of sample hashes found only in this row), sorted by row
    pub single_hits: Vec<(u32, u32)>,
    /// (shared set id, number of sample hashes whose posting set is this set), sorted by set id
    pub set_hits: Vec<(u32, u32)>,
    /// number of sample hashes found in the index
    pub n_found: u64,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct Entry32 {
    key: u32,
    val: u32,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct Entry64 {
    key: u64,
    val: u32,
    _pad: u32,
}

impl Index {
    pub fn open(path: &Path) -> Result<Index, String> {
        if cfg!(target_endian = "big") {
            return Err("YACHT index files are only supported on little-endian platforms".into());
        }
        let file = File::open(path).map_err(|e| format!("cannot open index {}: {}", path.display(), e))?;
        // SAFETY: the file is opened read-only; index files are replaced atomically (rename), never modified in place.
        let mmap = unsafe { Mmap::map(&file) }.map_err(|e| format!("cannot map index {}: {}", path.display(), e))?;
        let hdr = Header::from_bytes(&mmap)?;
        hdr.validate(mmap.len() as u64)?;
        Ok(Index { mmap, hdr })
    }

    fn slice<T>(&self, off: u64, n: u64) -> &[T] {
        // offsets/lengths were bounds- and alignment-checked in Header::validate; mmap is page aligned
        let ptr = self.mmap[off as usize..].as_ptr();
        debug_assert_eq!(ptr as usize % std::mem::align_of::<T>(), 0);
        unsafe { std::slice::from_raw_parts(ptr as *const T, n as usize) }
    }

    pub fn sizes(&self) -> &[u32] {
        self.slice(self.hdr.sizes_off, self.hdr.n_rows)
    }

    pub fn unique_counts(&self) -> &[u32] {
        self.slice(self.hdr.unique_off, self.hdr.n_rows)
    }


    fn set_offsets(&self) -> &[u64] {
        self.slice(self.hdr.set_offsets_off, self.hdr.n_sets + 1)
    }

    fn set_members(&self) -> &[u32] {
        self.slice(self.hdr.set_members_off, self.hdr.n_members)
    }

    fn set_weights(&self) -> &[u32] {
        self.slice(self.hdr.set_weights_off, self.hdr.n_sets)
    }

    pub fn name(&self, row: usize) -> Result<String, String> {
        let offs: &[u64] = self.slice(self.hdr.name_offsets_off, self.hdr.n_rows + 1);
        let bytes: &[u8] = self.slice(self.hdr.name_bytes_off, self.hdr.name_bytes_len);
        let (a, b) = (offs[row] as usize, offs[row + 1] as usize);
        if a > b || b > bytes.len() {
            return Err("corrupt YACHT index (name table)".into());
        }
        String::from_utf8(bytes[a..b].to_vec()).map_err(|_| "corrupt YACHT index (name is not UTF-8)".into())
    }

    /// Look up one hash; returns the stored value (row id, or n_rows + set id).
    #[inline]
    fn lookup<E: Copy, D: Copy + Into<u64>>(
        entries: &[E],
        dir: &[D],
        key_shift: u64,
        key_mask: u64,
        h: u64,
        key_of: impl Fn(&E) -> u64,
        val_of: impl Fn(&E) -> u32,
    ) -> Option<u32> {
        let b = crate::format::bucket_of(h, key_shift) as usize;
        if b + 1 >= dir.len() {
            return None;
        }
        let (lo, hi) = (dir[b].into() as usize, dir[b + 1].into() as usize);
        let bucket = &entries[lo..hi];
        let key = h & key_mask;
        bucket
            .binary_search_by(|e| key_of(e).cmp(&key))
            .ok()
            .map(|i| val_of(&bucket[i]))
    }

    fn lookup_all<E: Copy, D: Copy + Into<u64>>(
        &self,
        hashes: &[u64],
        entries: &[E],
        dir: &[D],
        key_of: impl Fn(&E) -> u64 + Copy,
        val_of: impl Fn(&E) -> u32 + Copy,
        mut push: impl FnMut(u32),
    ) {
        let key_shift = self.hdr.key_shift;
        let key_mask: u64 = if key_shift >= 64 { u64::MAX } else { (1u64 << key_shift) - 1 };
        for &h in hashes {
            if let Some(v) = Self::lookup(entries, dir, key_shift, key_mask, h, key_of, val_of) {
                push(v);
            }
        }
    }

    /// Phase A. `hashes` need not be sorted or unique (duplicates are ignored).
    pub fn query(&self, hashes: &[u64]) -> QueryResult {
        let owned;
        let hs: &[u64] = if hashes.windows(2).all(|w| w[0] < w[1]) {
            hashes
        } else {
            let mut v = hashes.to_vec();
            v.sort_unstable();
            v.dedup();
            owned = v;
            &owned
        };
        let n_rows = self.hdr.n_rows as u32;
        let mut single: Vec<u32> = Vec::new();
        let mut sets: Vec<u32> = Vec::new();
        let mut n_found = 0u64;
        let push = |v: u32| {
            n_found += 1;
            if v < n_rows {
                single.push(v);
            } else {
                sets.push(v - n_rows);
            }
        };
        let nb = self.hdr.n_buckets() + 1;
        match (self.hdr.key64(), self.hdr.dir64()) {
            (true, true) => {
                let e: &[Entry64] = self.slice(self.hdr.entries_off, self.hdr.n_keys);
                let d: &[u64] = self.slice(self.hdr.dir_off, nb);
                self.lookup_all(hs, e, d, |e| e.key, |e| e.val, push)
            }
            (true, false) => {
                let e: &[Entry64] = self.slice(self.hdr.entries_off, self.hdr.n_keys);
                let d: &[u32] = self.slice(self.hdr.dir_off, nb);
                self.lookup_all(hs, e, d, |e| e.key, |e| e.val, push)
            }
            (false, true) => {
                let e: &[Entry32] = self.slice(self.hdr.entries_off, self.hdr.n_keys);
                let d: &[u64] = self.slice(self.hdr.dir_off, nb);
                self.lookup_all(hs, e, d, |e| e.key as u64, |e| e.val, push)
            }
            (false, false) => {
                let e: &[Entry32] = self.slice(self.hdr.entries_off, self.hdr.n_keys);
                let d: &[u32] = self.slice(self.hdr.dir_off, nb);
                self.lookup_all(hs, e, d, |e| e.key as u64, |e| e.val, push)
            }
        }
        QueryResult {
            single_hits: run_length(single),
            set_hits: run_length(sets),
            n_found,
        }
    }

    /// Rows that share at least one hash with the sample (sorted, unique).
    pub fn marked_rows(&self, q: &QueryResult) -> Vec<u32> {
        let mut rows: Vec<u32> = q.single_hits.iter().map(|&(r, _)| r).collect();
        let offs = self.set_offsets();
        let members = self.set_members();
        for &(s, _) in &q.set_hits {
            rows.extend_from_slice(&members[offs[s as usize] as usize..offs[s as usize + 1] as usize]);
        }
        rows.sort_unstable();
        rows.dedup();
        rows
    }

    /// Phase B. `in_n[row]` marks the nontrivial rows. Returns (excl, matches), each of length n_rows
    /// (zero for rows outside N).
    pub fn exclusive_counts(&self, q: &QueryResult, in_n: &[bool]) -> Result<(Vec<u64>, Vec<u64>), String> {
        let n_rows = self.hdr.n_rows as usize;
        if in_n.len() != n_rows {
            return Err(format!("mask has length {}, expected {}", in_n.len(), n_rows));
        }
        let mut excl = vec![0u64; n_rows];
        let mut matches = vec![0u64; n_rows];
        let unique = self.unique_counts();
        for r in 0..n_rows {
            if in_n[r] {
                excl[r] = unique[r] as u64;
            }
        }
        for &(r, c) in &q.single_hits {
            if in_n[r as usize] {
                matches[r as usize] += c as u64;
            }
        }
        let offs = self.set_offsets();
        let members = self.set_members();
        let weights = self.set_weights();
        let mut hits = q.set_hits.iter().peekable();
        for s in 0..self.hdr.n_sets as usize {
            let mut owner = u32::MAX;
            let mut count = 0u32;
            for &m in &members[offs[s] as usize..offs[s + 1] as usize] {
                if in_n[m as usize] {
                    count += 1;
                    if count > 1 {
                        break;
                    }
                    owner = m;
                }
            }
            let hit = match hits.peek() {
                Some(&&(hs, c)) if hs as usize == s => {
                    hits.next();
                    c
                }
                _ => 0,
            };
            if count == 1 {
                excl[owner as usize] += weights[s] as u64;
                matches[owner as usize] += hit as u64;
            }
        }
        Ok((excl, matches))
    }
}

fn run_length(mut xs: Vec<u32>) -> Vec<(u32, u32)> {
    xs.sort_unstable();
    let mut out: Vec<(u32, u32)> = Vec::new();
    for x in xs {
        match out.last_mut() {
            Some((v, c)) if *v == x => *c += 1,
            _ => out.push((x, 1)),
        }
    }
    out
}
