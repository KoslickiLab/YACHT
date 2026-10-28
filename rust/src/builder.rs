//! Build a YACHT run index from the genome signatures of a trained database.

use crate::format::{align_up, choose_layout, Header, FLAG_DIR64, FLAG_KEY64, FLAG_MIXED_SEED, FORMAT_VERSION, HEADER_SIZE};
use crate::sig::load_genome_sketch;
use rustc_hash::FxHashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Default)]
pub struct BuildStats {
    pub n_rows: u64,
    pub n_hashes: u64,
    pub n_keys: u64,
    pub n_unique_keys: u64,
    pub n_sets: u64,
    pub n_members: u64,
    pub file_len: u64,
}

fn write_u32s<W: Write>(w: &mut W, xs: &[u32]) -> std::io::Result<()> {
    for x in xs {
        w.write_all(&x.to_le_bytes())?;
    }
    Ok(())
}

fn write_u64s<W: Write>(w: &mut W, xs: &[u64]) -> std::io::Result<()> {
    for x in xs {
        w.write_all(&x.to_le_bytes())?;
    }
    Ok(())
}

fn pad_to<W: Write + Seek>(w: &mut W, off: u64) -> std::io::Result<()> {
    let pos = w.stream_position()?;
    assert!(pos <= off, "section overlap while writing index");
    if pos < off {
        w.write_all(&vec![0u8; (off - pos) as usize])?;
    }
    Ok(())
}

/// Build the index for the genomes in `sig_paths` (one per manifest row, in manifest order) and write it
/// atomically to `out_path`. `digest` identifies the manifest the index belongs to.
pub fn build_index(
    sig_paths: &[PathBuf],
    ksize: u64,
    scaled: u64,
    digest: [u8; 32],
    out_path: &Path,
) -> Result<BuildStats, String> {
    build_index_opts(sig_paths, ksize, scaled, digest, out_path, false)
}

/// As build_index; `force_dir64` selects the 64-bit bucket directory even for small indexes (for testing).
pub fn build_index_opts(
    sig_paths: &[PathBuf],
    ksize: u64,
    scaled: u64,
    digest: [u8; 32],
    out_path: &Path,
    force_dir64: bool,
) -> Result<BuildStats, String> {
    let n_rows = sig_paths.len();
    if n_rows as u64 >= (u32::MAX / 2) as u64 {
        return Err("too many genomes for a YACHT index".into());
    }

    // 1. load every genome sketch (sorted, de-duplicated hashes) into one flat array
    let mut all: Vec<u64> = Vec::new();
    let mut starts: Vec<usize> = Vec::with_capacity(n_rows + 1);
    let mut names: Vec<String> = Vec::with_capacity(n_rows);
    let mut sizes: Vec<u32> = Vec::with_capacity(n_rows);
    let mut seeds: Vec<u64> = Vec::new();
    for path in sig_paths {
        let g = load_genome_sketch(path, ksize)?;
        if !seeds.contains(&g.seed) {
            seeds.push(g.seed);
        }
        starts.push(all.len());
        sizes.push(u32::try_from(g.hashes.len()).map_err(|_| "genome sketch too large")?);
        all.extend_from_slice(&g.hashes);
        names.push(g.name);
    }
    starts.push(all.len());
    let n_hashes = all.len() as u64;
    let max_hash = all.iter().copied().max().unwrap_or(0);

    let (bucket_bits, key_shift, key64) = choose_layout(n_hashes, max_hash);
    let n_buckets = 1u64 << bucket_bits;
    let key_mask: u64 = if key_shift >= 64 { u64::MAX } else { (1u64 << key_shift) - 1 };

    // 2. fixed-size sections before the entries
    let mut name_offsets: Vec<u64> = Vec::with_capacity(n_rows + 1);
    let mut name_bytes: Vec<u8> = Vec::new();
    for n in &names {
        name_offsets.push(name_bytes.len() as u64);
        name_bytes.extend_from_slice(n.as_bytes());
    }
    name_offsets.push(name_bytes.len() as u64);
    drop(names);

    let mut h = Header {
        version: FORMAT_VERSION,
        flags: (if key64 { FLAG_KEY64 } else { 0 })
            | (if force_dir64 || n_hashes >= u32::MAX as u64 { FLAG_DIR64 } else { 0 })
            | (if seeds.len() > 1 { FLAG_MIXED_SEED } else { 0 }),
        seed: seeds.first().copied().unwrap_or(42),
        ksize,
        scaled,
        n_rows: n_rows as u64,
        bucket_bits,
        key_shift,
        digest,
        ..Default::default()
    };
    h.sizes_off = align_up(HEADER_SIZE as u64);
    h.unique_off = align_up(h.sizes_off + 4 * n_rows as u64);
    h.name_offsets_off = align_up(h.unique_off + 4 * n_rows as u64);
    h.name_bytes_off = align_up(h.name_offsets_off + 8 * (n_rows as u64 + 1));
    h.name_bytes_len = name_bytes.len() as u64;
    h.dir_off = align_up(h.name_bytes_off + h.name_bytes_len);
    h.entries_off = align_up(h.dir_off + h.dir_entry_size() * (n_buckets + 1));

    let parent = out_path.parent().unwrap_or_else(|| Path::new("."));
    let tmp_path = parent.join(format!(
        ".{}.tmp{}",
        out_path.file_name().and_then(|s| s.to_str()).unwrap_or("index"),
        std::process::id()
    ));
    let io = |e: std::io::Error| format!("cannot write index file {}: {}", tmp_path.display(), e);
    let file = File::create(&tmp_path).map_err(io)?;
    let mut w = BufWriter::with_capacity(1 << 22, file);
    w.seek(SeekFrom::Start(h.entries_off)).map_err(io)?;

    // 3. stream the (key, val) entries in global hash order, partition by partition
    let bits = crate::format::bit_width(max_hash).max(1);
    let pbits = (crate::format::bit_width(n_hashes) as i64 - 21).clamp(0, bits.min(16) as i64) as u64;
    let pshift = bits - pbits;
    let n_parts = 1u64 << pbits;
    let mut cursor: Vec<usize> = starts[..n_rows].to_vec();
    let mut dir: Vec<u64> = vec![0; n_buckets as usize + 1];
    let mut unique: Vec<u32> = vec![0; n_rows];
    let mut set_ids: FxHashMap<Box<[u32]>, u32> = FxHashMap::default();
    let mut set_offsets: Vec<u64> = vec![0];
    let mut set_members: Vec<u32> = Vec::new();
    let mut set_weights: Vec<u32> = Vec::new();
    let mut pairs: Vec<u128> = Vec::new();
    let mut group: Vec<u32> = Vec::new();
    let mut n_keys: u64 = 0;
    let mut n_unique_keys: u64 = 0;

    for p in 0..n_parts {
        pairs.clear();
        for r in 0..n_rows {
            let seg = &all[cursor[r]..starts[r + 1]];
            let take = if p + 1 == n_parts {
                seg.len()
            } else {
                seg.partition_point(|&x| (x >> pshift) <= p)
            };
            for &x in &seg[..take] {
                pairs.push(((x as u128) << 32) | r as u128);
            }
            cursor[r] += take;
        }
        pairs.sort_unstable();
        let mut i = 0;
        while i < pairs.len() {
            let hash = (pairs[i] >> 32) as u64;
            group.clear();
            while i < pairs.len() && (pairs[i] >> 32) as u64 == hash {
                group.push(pairs[i] as u32);
                i += 1;
            }
            let val: u32 = if group.len() == 1 {
                unique[group[0] as usize] += 1;
                n_unique_keys += 1;
                group[0]
            } else {
                let id = match set_ids.get(group.as_slice()) {
                    Some(&id) => id,
                    None => {
                        let id = set_weights.len() as u32;
                        set_ids.insert(group.clone().into_boxed_slice(), id);
                        set_members.extend_from_slice(&group);
                        set_offsets.push(set_members.len() as u64);
                        set_weights.push(0);
                        id
                    }
                };
                set_weights[id as usize] += 1;
                n_rows as u32 + id
            };
            if (n_rows as u64 + set_weights.len() as u64) >= u32::MAX as u64 {
                return Err("too many distinct shared hash sets for a YACHT index".into());
            }
            let key = hash & key_mask;
            if key64 {
                w.write_all(&key.to_le_bytes()).map_err(io)?;
                w.write_all(&val.to_le_bytes()).map_err(io)?;
                w.write_all(&0u32.to_le_bytes()).map_err(io)?;
            } else {
                w.write_all(&(key as u32).to_le_bytes()).map_err(io)?;
                w.write_all(&val.to_le_bytes()).map_err(io)?;
            }
            dir[crate::format::bucket_of(hash, key_shift) as usize + 1] += 1;
            n_keys += 1;
        }
    }
    drop(all);
    drop(set_ids);
    for b in 0..n_buckets as usize {
        dir[b + 1] += dir[b];
    }

    // 4. shared posting sets
    h.n_keys = n_keys;
    h.n_sets = set_weights.len() as u64;
    h.n_members = set_members.len() as u64;
    h.set_offsets_off = align_up(h.entries_off + h.entry_size() * n_keys);
    pad_to(&mut w, h.set_offsets_off).map_err(io)?;
    write_u64s(&mut w, &set_offsets).map_err(io)?;
    h.set_members_off = align_up(h.set_offsets_off + 8 * (h.n_sets + 1));
    pad_to(&mut w, h.set_members_off).map_err(io)?;
    write_u32s(&mut w, &set_members).map_err(io)?;
    h.set_weights_off = align_up(h.set_members_off + 4 * h.n_members);
    pad_to(&mut w, h.set_weights_off).map_err(io)?;
    write_u32s(&mut w, &set_weights).map_err(io)?;
    h.file_len = align_up(h.set_weights_off + 4 * h.n_sets);
    pad_to(&mut w, h.file_len).map_err(io)?;

    // 5. header and the sections before the entries
    w.seek(SeekFrom::Start(0)).map_err(io)?;
    w.write_all(&h.to_bytes()).map_err(io)?;
    pad_to(&mut w, h.sizes_off).map_err(io)?;
    write_u32s(&mut w, &sizes).map_err(io)?;
    pad_to(&mut w, h.unique_off).map_err(io)?;
    write_u32s(&mut w, &unique).map_err(io)?;
    pad_to(&mut w, h.name_offsets_off).map_err(io)?;
    write_u64s(&mut w, &name_offsets).map_err(io)?;
    pad_to(&mut w, h.name_bytes_off).map_err(io)?;
    w.write_all(&name_bytes).map_err(io)?;
    pad_to(&mut w, h.dir_off).map_err(io)?;
    if h.dir64() {
        write_u64s(&mut w, &dir).map_err(io)?;
    } else {
        let dir32: Vec<u32> = dir.iter().map(|&x| x as u32).collect();
        write_u32s(&mut w, &dir32).map_err(io)?;
    }
    let file = w.into_inner().map_err(|e| io(e.into_error()))?;
    file.sync_all().map_err(io)?;
    drop(file);
    fs::rename(&tmp_path, out_path).map_err(|e| {
        let _ = fs::remove_file(&tmp_path);
        format!("cannot move index into place at {}: {}", out_path.display(), e)
    })?;

    Ok(BuildStats {
        n_rows: n_rows as u64,
        n_hashes,
        n_keys,
        n_unique_keys,
        n_sets: h.n_sets,
        n_members: h.n_members,
        file_len: h.file_len,
    })
}
