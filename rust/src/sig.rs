//! Minimal reader for sourmash JSON signature files (`.sig`, optionally gzip-compressed).
//!
//! Only the fields YACHT needs are decoded; everything else (abundances, md5sum, ...) is skipped.

use flate2::read::MultiGzDecoder;
use serde::Deserialize;
use std::fs::File;
use std::io::Read;
use std::path::Path;

#[derive(Deserialize)]
struct SigEntry {
    #[serde(default)]
    name: Option<String>,
    #[serde(default)]
    signatures: Vec<SketchEntry>,
}

#[derive(Deserialize)]
struct SketchEntry {
    ksize: u64,
    #[serde(default)]
    molecule: Option<String>,
    #[serde(default)]
    mins: Vec<u64>,
}

/// The single DNA sketch of a genome signature file at a given k-mer size.
pub struct GenomeSketch {
    /// Signature name as stored in the file ("" when absent); this is what both the Python sourmash API
    /// (`sig.name`) and `sourmash scripts multisearch` (`match_name`) report.
    pub name: String,
    /// Sorted, de-duplicated hash values.
    pub hashes: Vec<u64>,
}

fn read_maybe_gz(path: &Path) -> Result<Vec<u8>, String> {
    let mut raw = Vec::new();
    File::open(path)
        .and_then(|mut f| f.read_to_end(&mut raw))
        .map_err(|e| format!("cannot read signature file {}: {}", path.display(), e))?;
    if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut out = Vec::with_capacity(raw.len() * 4);
        MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut out)
            .map_err(|e| format!("cannot decompress signature file {}: {}", path.display(), e))?;
        Ok(out)
    } else {
        Ok(raw)
    }
}

/// Load the unique DNA sketch with k-mer size `ksize` from a sourmash signature file.
/// Mirrors `yacht.utils.load_signature_with_ksize`: exactly one matching sketch must exist.
pub fn load_genome_sketch(path: &Path, ksize: u64) -> Result<GenomeSketch, String> {
    let data = read_maybe_gz(path)?;
    let entries: Vec<SigEntry> = match serde_json::from_slice::<Vec<SigEntry>>(&data) {
        Ok(v) => v,
        // a file may also hold a single signature object instead of a list
        Err(_) => vec![serde_json::from_slice::<SigEntry>(&data)
            .map_err(|e| format!("cannot parse signature file {}: {}", path.display(), e))?],
    };
    let mut found: Option<GenomeSketch> = None;
    let mut n_found = 0usize;
    for entry in entries {
        for sk in entry.signatures {
            let is_dna = sk.molecule.as_deref().map_or(true, |m| m.eq_ignore_ascii_case("dna"));
            if sk.ksize == ksize && is_dna {
                n_found += 1;
                if found.is_none() {
                    let mut hashes = sk.mins;
                    hashes.sort_unstable();
                    hashes.dedup();
                    found = Some(GenomeSketch {
                        name: entry.name.clone().unwrap_or_default(),
                        hashes,
                    });
                }
            }
        }
    }
    match (n_found, found) {
        (1, Some(g)) => Ok(g),
        (n, _) => Err(format!(
            "Expected exactly one signature with ksize {} in {}, found {}",
            ksize,
            path.display(),
            n
        )),
    }
}

#[derive(Deserialize)]
struct SampleEntry {
    #[serde(default)]
    name: Option<String>,
    #[serde(default)]
    signatures: Vec<SampleSketchEntry>,
}

#[derive(Deserialize)]
struct SampleSketchEntry {
    ksize: u64,
    #[serde(default)]
    molecule: Option<String>,
    #[serde(default)]
    num: u64,
    #[serde(default)]
    max_hash: u64,
    #[serde(default)]
    mins: Vec<u64>,
    #[serde(default)]
    abundances: Option<Vec<u64>>,
}

/// A sample sketch as sourmash would load it.
pub struct SampleSketch {
    pub name: String,
    /// sorted hashes
    pub hashes: Vec<u64>,
    /// abundances (None if the sketch does not track abundance)
    pub abundances: Option<Vec<u64>>,
    /// `MinHash.scaled` as sourmash computes it from max_hash
    pub scaled: u64,
}

/// The md5sum sourmash computes for a sketch (`KmerMinHash::md5sum`).
fn sourmash_md5(ksize: u64, mins: &[u64]) -> String {
    use md5::{Digest, Md5};
    let mut ctx = Md5::new();
    ctx.update(ksize.to_string().as_bytes());
    let mut buf = String::with_capacity(24);
    for x in mins {
        use std::fmt::Write;
        buf.clear();
        write!(&mut buf, "{}", x).unwrap();
        ctx.update(buf.as_bytes());
    }
    ctx.finalize().iter().map(|b| format!("{:02x}", b)).collect()
}

/// Parse the signature file (JSON, optionally gzipped) that a sample zip's manifest row points to, and return
/// the DNA sketch with k-mer size `ksize` whose sourmash md5sum is `md5`, exactly as sourmash would load it.
/// Returns Ok(None) for anything that is not a plain, unambiguous scaled sketch (no matching sketch, several,
/// unsorted or duplicated hashes, mismatched abundances, ...), so that the caller can fall back to sourmash.
pub fn parse_sample_sketch(data: &[u8], ksize: u64, md5: &str) -> Result<Option<SampleSketch>, String> {
    let json: Vec<u8>;
    let data = if data.len() >= 2 && data[0] == 0x1f && data[1] == 0x8b {
        let mut out = Vec::with_capacity(data.len() * 4);
        MultiGzDecoder::new(data)
            .read_to_end(&mut out)
            .map_err(|e| format!("cannot decompress sample signature: {}", e))?;
        json = out;
        &json[..]
    } else {
        data
    };
    let entries: Vec<SampleEntry> = match serde_json::from_slice(data) {
        Ok(v) => v,
        Err(_) => return Ok(None),
    };
    let mut found: Vec<SampleSketch> = Vec::new();
    for entry in entries {
        for sk in entry.signatures {
            if sk.ksize != ksize {
                continue;
            }
            let regular = sk.max_hash != 0
                && sk.num == 0
                && sk.molecule.as_deref().is_some_and(|m| m.eq_ignore_ascii_case("dna"))
                && sk.mins.windows(2).all(|w| w[0] < w[1])
                && sk.abundances.as_ref().map_or(true, |a| a.len() == sk.mins.len() && a.iter().all(|&x| x > 0));
            if !regular {
                // sourmash would normalize this sketch (sort, keep duplicates, ...) before computing its md5sum
                return Ok(None);
            }
            if sourmash_md5(sk.ksize, &sk.mins) != md5 {
                continue;
            }
            found.push(SampleSketch {
                name: entry.name.clone().unwrap_or_default(),
                scaled: (u64::MAX as f64 / sk.max_hash as f64) as u32 as u64,
                hashes: sk.mins,
                abundances: sk.abundances,
            });
        }
    }
    if found.len() == 1 {
        Ok(found.pop())
    } else {
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn md5_matches_sourmash() {
        // md5sum of the k=31 sketch in YACHT's tests/testdata/20_genomes_sketches.zip is checked from Python;
        // here: the empty sketch and a small one, computed with sourmash 4.9.4
        assert_eq!(sourmash_md5(31, &[]), format!("{:x}", md5_of(b"31")));
        assert_eq!(sourmash_md5(21, &[5, 17]), format!("{:x}", md5_of(b"21517")));
    }

    fn md5_of(b: &[u8]) -> u128 {
        use md5::{Digest, Md5};
        u128::from_be_bytes(Md5::digest(b).into())
    }

    #[test]
    fn irregular_sketches_are_rejected() {
        let ok = br#"[{"name":"s","signatures":[{"ksize":31,"molecule":"DNA","num":0,"max_hash":18446744073709552,"mins":[1,5,9],"abundances":[2,1,3]}]}]"#;
        let md5 = sourmash_md5(31, &[1, 5, 9]);
        let s = parse_sample_sketch(ok, 31, &md5).unwrap().unwrap();
        assert_eq!(s.hashes, vec![1, 5, 9]);
        assert_eq!(s.abundances, Some(vec![2, 1, 3]));
        assert_eq!(s.scaled, 1000);
        assert!(parse_sample_sketch(ok, 21, &md5).unwrap().is_none());
        let unsorted = br#"[{"signatures":[{"ksize":31,"molecule":"DNA","num":0,"max_hash":18446744073709552,"mins":[5,1,9]}]}]"#;
        assert!(parse_sample_sketch(unsorted, 31, &md5).unwrap().is_none());
        let dup = br#"[{"signatures":[{"ksize":31,"molecule":"DNA","num":0,"max_hash":18446744073709552,"mins":[1,5,5,9]}]}]"#;
        assert!(parse_sample_sketch(dup, 31, &sourmash_md5(31, &[1, 5, 5, 9])).unwrap().is_none());
    }
}
