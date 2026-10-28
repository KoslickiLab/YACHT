//! On-disk layout of the YACHT run index (`<prefix>_index.yidx`).
//!
//! All integers are little-endian. Every section starts on a 64-byte boundary.
//!
//! ```text
//! [header: 256 bytes; includes ksize, scaled and the hash seed of the genome sketches]
//! SIZES        u32[n_rows]          number of hashes in each genome sketch (manifest row order)
//! UNIQUE       u32[n_rows]          number of hashes found in this row and in no other row
//! NAME_OFFSETS u64[n_rows + 1]      byte offsets into NAME_BYTES
//! NAME_BYTES   u8[...]              UTF-8 signature names (what `sourmash scripts multisearch` reports as match_name)
//! DIR          u32|u64[2^bucket_bits+1] bucket directory: entries of bucket b are ENTRIES[DIR[b]..DIR[b+1]]
//!                                   (u64 only if there are >= 2^32 entries)
//! ENTRIES      (key, val)[n_keys]   key = hash & ((1 << key_shift) - 1), sorted within bucket (bucket = hash >> key_shift)
//!                                   val < n_rows: the single row containing the hash
//!                                   val >= n_rows: shared posting set id (val - n_rows)
//! SET_OFFSETS  u64[n_sets + 1]      offsets into SET_MEMBERS
//! SET_MEMBERS  u32[n_members]       sorted row ids of each distinct shared posting set
//! SET_WEIGHTS  u32[n_sets]          number of hashes whose posting set is exactly this set
//! ```
//! Entries are 8 bytes (u32 key, u32 val) when `key_shift <= 32`, otherwise 16 bytes (u64 key, u32 val, u32 pad).

pub const MAGIC: &[u8; 8] = b"YACHTIDX";
pub const FORMAT_VERSION: u32 = 1;
pub const HEADER_SIZE: usize = 256;
pub const ALIGN: u64 = 64;
pub const FLAG_KEY64: u32 = 1;
pub const FLAG_DIR64: u32 = 2;
/// the genome sketches do not all use the same hash seed
pub const FLAG_MIXED_SEED: u32 = 4;

#[derive(Debug, Clone, Default)]
pub struct Header {
    pub version: u32,
    pub flags: u32,
    pub ksize: u64,
    pub scaled: u64,
    pub n_rows: u64,
    pub n_keys: u64,
    pub n_sets: u64,
    pub n_members: u64,
    pub bucket_bits: u64,
    pub key_shift: u64,
    pub digest: [u8; 32],
    pub sizes_off: u64,
    pub unique_off: u64,
    pub name_offsets_off: u64,
    pub name_bytes_off: u64,
    pub name_bytes_len: u64,
    pub dir_off: u64,
    pub entries_off: u64,
    pub set_offsets_off: u64,
    pub set_members_off: u64,
    pub set_weights_off: u64,
    pub file_len: u64,
    /// hash seed of the genome sketches (sourmash default: 42)
    pub seed: u64,
}

pub fn align_up(x: u64) -> u64 {
    x.div_ceil(ALIGN) * ALIGN
}

impl Header {
    pub fn key64(&self) -> bool {
        self.flags & FLAG_KEY64 != 0
    }

    pub fn dir64(&self) -> bool {
        self.flags & FLAG_DIR64 != 0
    }

    pub fn dir_entry_size(&self) -> u64 {
        if self.dir64() {
            8
        } else {
            4
        }
    }

    pub fn entry_size(&self) -> u64 {
        if self.key64() {
            16
        } else {
            8
        }
    }

    pub fn n_buckets(&self) -> u64 {
        1u64 << self.bucket_bits
    }

    pub fn to_bytes(&self) -> [u8; HEADER_SIZE] {
        let mut b = [0u8; HEADER_SIZE];
        b[0..8].copy_from_slice(MAGIC);
        b[8..12].copy_from_slice(&self.version.to_le_bytes());
        b[12..16].copy_from_slice(&self.flags.to_le_bytes());
        let fields = [
            self.ksize,
            self.scaled,
            self.n_rows,
            self.n_keys,
            self.n_sets,
            self.n_members,
            self.bucket_bits,
            self.key_shift,
        ];
        for (i, v) in fields.iter().enumerate() {
            b[16 + 8 * i..24 + 8 * i].copy_from_slice(&v.to_le_bytes());
        }
        b[80..112].copy_from_slice(&self.digest);
        let offs = [
            self.sizes_off,
            self.unique_off,
            self.name_offsets_off,
            self.name_bytes_off,
            self.name_bytes_len,
            self.dir_off,
            self.entries_off,
            self.set_offsets_off,
            self.set_members_off,
            self.set_weights_off,
            self.file_len,
            self.seed,
        ];
        for (i, v) in offs.iter().enumerate() {
            b[112 + 8 * i..120 + 8 * i].copy_from_slice(&v.to_le_bytes());
        }
        b
    }

    pub fn from_bytes(b: &[u8]) -> Result<Header, String> {
        if b.len() < HEADER_SIZE {
            return Err("file is too small to be a YACHT index".into());
        }
        if &b[0..8] != MAGIC {
            return Err("not a YACHT index file (bad magic)".into());
        }
        let u32_at = |o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
        let u64_at = |o: usize| u64::from_le_bytes(b[o..o + 8].try_into().unwrap());
        let mut h = Header {
            version: u32_at(8),
            flags: u32_at(12),
            ..Default::default()
        };
        if h.version != FORMAT_VERSION {
            return Err(format!(
                "unsupported YACHT index format version {} (this YACHT supports version {})",
                h.version, FORMAT_VERSION
            ));
        }
        h.ksize = u64_at(16);
        h.scaled = u64_at(24);
        h.n_rows = u64_at(32);
        h.n_keys = u64_at(40);
        h.n_sets = u64_at(48);
        h.n_members = u64_at(56);
        h.bucket_bits = u64_at(64);
        h.key_shift = u64_at(72);
        h.digest.copy_from_slice(&b[80..112]);
        h.sizes_off = u64_at(112);
        h.unique_off = u64_at(120);
        h.name_offsets_off = u64_at(128);
        h.name_bytes_off = u64_at(136);
        h.name_bytes_len = u64_at(144);
        h.dir_off = u64_at(152);
        h.entries_off = u64_at(160);
        h.set_offsets_off = u64_at(168);
        h.set_members_off = u64_at(176);
        h.set_weights_off = u64_at(184);
        h.file_len = u64_at(192);
        h.seed = u64_at(200);
        Ok(h)
    }

    /// Check that every section lies inside a file of length `len` and is suitably aligned.
    pub fn validate(&self, len: u64) -> Result<(), String> {
        if self.file_len != len {
            return Err(format!(
                "index file is truncated or corrupt (expected {} bytes, found {})",
                self.file_len, len
            ));
        }
        if self.bucket_bits > 40 || self.key_shift > 64 || self.n_rows >= u32::MAX as u64 {
            return Err("corrupt YACHT index header".into());
        }
        if self.n_rows + self.n_sets >= u32::MAX as u64 {
            return Err("corrupt YACHT index header (too many rows/sets)".into());
        }
        let sections = [
            (self.sizes_off, 4 * self.n_rows),
            (self.unique_off, 4 * self.n_rows),
            (self.name_offsets_off, 8 * (self.n_rows + 1)),
            (self.name_bytes_off, self.name_bytes_len),
            (self.dir_off, self.dir_entry_size() * (self.n_buckets() + 1)),
            (self.entries_off, self.entry_size() * self.n_keys),
            (self.set_offsets_off, 8 * (self.n_sets + 1)),
            (self.set_members_off, 4 * self.n_members),
            (self.set_weights_off, 4 * self.n_sets),
        ];
        for (off, size) in sections {
            if off % ALIGN != 0 || off < HEADER_SIZE as u64 || off.checked_add(size).map_or(true, |end| end > len) {
                return Err("corrupt YACHT index (section out of bounds)".into());
            }
        }
        Ok(())
    }
}

/// Number of significant bits needed to represent `x` (0 for x == 0).
pub fn bit_width(x: u64) -> u64 {
    64 - x.leading_zeros() as u64
}

/// Choose (bucket_bits, key_shift, key64) for `n_items` hashes whose maximum value is `max_hash`.
/// Aims for ~32 entries per bucket; uses more buckets (up to 8x) if that makes 32-bit keys possible.
pub fn choose_layout(n_items: u64, max_hash: u64) -> (u64, u64, bool) {
    let bits = bit_width(max_hash).max(1);
    let want = (bit_width(n_items.saturating_sub(1)) as i64 - 5).clamp(0, 28) as u64;
    let need32 = bits.saturating_sub(32);
    let bucket_bits = if need32 <= want || (need32 <= want + 3 && need32 <= 28) {
        want.max(need32)
    } else {
        want
    }
    .min(bits);
    let key_shift = bits - bucket_bits;
    (bucket_bits, key_shift, key_shift > 32)
}

/// Bucket of hash `h` for a layout with the given key shift.
#[inline]
pub fn bucket_of(h: u64, key_shift: u64) -> u64 {
    if key_shift >= 64 {
        0
    } else {
        h >> key_shift
    }
}
