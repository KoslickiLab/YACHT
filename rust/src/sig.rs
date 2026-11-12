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
    #[serde(default = "default_seed")]
    seed: u64,
    #[serde(default)]
    mins: Vec<u64>,
}

fn default_seed() -> u64 {
    42
}

/// The single DNA sketch of a genome signature file at a given k-mer size.
pub struct GenomeSketch {
    /// Signature name as stored in the file ("" when absent); this is what both the Python sourmash API
    /// (`sig.name`) and `sourmash scripts multisearch` (`match_name`) report.
    pub name: String,
    /// Sorted, de-duplicated hash values.
    pub hashes: Vec<u64>,
    /// Hash seed of the sketch.
    pub seed: u64,
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
                        seed: sk.seed,
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
    #[serde(default = "default_seed")]
    seed: u64,
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
    /// hash seed
    pub seed: u64,
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

/// Outcome of reading a sample sketch with [`parse_sample_sketch`].
pub enum SampleParse {
    /// The sketch, exactly as sourmash would load it.
    Sketch(SampleSketch),
    /// Anything else that sourmash can still load (no unique matching sketch, unsorted hashes, num sketches, ...).
    NeedsSourmash,
    /// A scaled sketch that is not a valid FracMinHash sketch; the reason, see [`fracminhash_problem`].
    Invalid(String),
}

/// Why a scaled sketch is not a valid FracMinHash sketch, or None if it is one. A FracMinHash sketch is a set
/// of distinct hash values, each at most `max_hash`, with an abundance of at least 1 for every hash if
/// abundances are tracked; sourmash's sketchers only write such sketches. (sourmash's loader keeps repeated
/// hashes, so that `len(minhash)` and `minhash.hashes` disagree, and pairs up mismatched abundance lists.)
pub fn fracminhash_problem(mins: &[u64], abundances: Option<&[u64]>, max_hash: u64) -> Option<String> {
    if let Some(a) = abundances {
        if a.len() != mins.len() {
            return Some(format!("it has {} abundances for {} hashes", a.len(), mins.len()));
        }
        if a.contains(&0) {
            return Some("it has hashes with abundance 0".into());
        }
    }
    if mins.iter().any(|&h| h > max_hash) {
        return Some(format!("it has hash values above its max_hash ({})", max_hash));
    }
    if !mins.windows(2).all(|w| w[0] < w[1]) {
        let mut sorted = mins.to_vec();
        sorted.sort_unstable();
        if sorted.windows(2).any(|w| w[0] == w[1]) {
            return Some("it has repeated hash values".into());
        }
    }
    None
}

/// Parse the signature file (JSON, optionally gzipped) that a sample zip's manifest row points to, and return
/// the DNA sketch with k-mer size `ksize` whose sourmash md5sum is `md5`, exactly as sourmash would load it.
/// A scaled DNA sketch with this k-mer size that is not a valid FracMinHash sketch gives `Invalid`; anything
/// else that is not a plain, unambiguous scaled sketch (no matching sketch, several, unsorted hashes, ...)
/// gives `NeedsSourmash`, so that the caller can fall back to sourmash.
pub fn parse_sample_sketch(data: &[u8], ksize: u64, md5: &str) -> Result<SampleParse, String> {
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
        Err(_) => return Ok(SampleParse::NeedsSourmash),
    };
    let mut found: Vec<SampleSketch> = Vec::new();
    for entry in entries {
        for sk in entry.signatures {
            if sk.ksize != ksize {
                continue;
            }
            let scaled_dna = sk.max_hash != 0
                && sk.num == 0
                && sk.molecule.as_deref().is_some_and(|m| m.eq_ignore_ascii_case("dna"));
            if !scaled_dna {
                return Ok(SampleParse::NeedsSourmash);
            }
            if let Some(problem) = fracminhash_problem(&sk.mins, sk.abundances.as_deref(), sk.max_hash) {
                return Ok(SampleParse::Invalid(problem));
            }
            if !sk.mins.windows(2).all(|w| w[0] < w[1]) {
                // sourmash sorts the hashes before computing the md5sum
                return Ok(SampleParse::NeedsSourmash);
            }
            if sourmash_md5(sk.ksize, &sk.mins) != md5 {
                continue;
            }
            found.push(SampleSketch {
                name: entry.name.clone().unwrap_or_default(),
                scaled: (u64::MAX as f64 / sk.max_hash as f64) as u32 as u64,
                seed: sk.seed,
                hashes: sk.mins,
                abundances: sk.abundances,
            });
        }
    }
    match (found.pop(), found.is_empty()) {
        (Some(sketch), true) => Ok(SampleParse::Sketch(sketch)),
        _ => Ok(SampleParse::NeedsSourmash),
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

    fn parse(data: &[u8], ksize: u64, md5: &str) -> SampleParse {
        parse_sample_sketch(data, ksize, md5).unwrap()
    }

    fn sketch_json(mins: &str, abundances: &str) -> Vec<u8> {
        format!(
            r#"[{{"name":"s","signatures":[{{"ksize":31,"molecule":"DNA","num":0,"max_hash":18446744073709552,"mins":{}{}}}]}}]"#,
            mins, abundances
        )
        .into_bytes()
    }

    #[test]
    fn irregular_sketches_need_sourmash() {
        let ok = sketch_json("[1,5,9]", r#","abundances":[2,1,3]"#);
        let md5 = sourmash_md5(31, &[1, 5, 9]);
        match parse(&ok, 31, &md5) {
            SampleParse::Sketch(s) => {
                assert_eq!(s.hashes, vec![1, 5, 9]);
                assert_eq!(s.abundances, Some(vec![2, 1, 3]));
                assert_eq!(s.scaled, 1000);
            }
            _ => panic!("expected a sketch"),
        }
        assert!(matches!(parse(&ok, 21, &md5), SampleParse::NeedsSourmash));
        // unsorted but distinct hashes are a valid set (sourmash sorts them)
        assert!(matches!(parse(&sketch_json("[5,1,9]", ""), 31, &md5), SampleParse::NeedsSourmash));
        let num = br#"[{"signatures":[{"ksize":31,"molecule":"DNA","num":500,"max_hash":0,"mins":[1,5,5]}]}]"#;
        assert!(matches!(parse(num, 31, &md5), SampleParse::NeedsSourmash));
    }

    #[test]
    fn invalid_fracminhash_sketches_are_reported() {
        let md5 = sourmash_md5(31, &[1, 5, 9]);
        for (mins, abundances, reason) in [
            ("[1,5,5,9]", "", "repeated hash values"),
            ("[9,5,1,5]", "", "repeated hash values"),
            ("[1,5,18446744073709553]", "", "above its max_hash"),
            ("[1,5,9]", r#","abundances":[2,1]"#, "2 abundances for 3 hashes"),
            ("[1,5,9]", r#","abundances":[2,0,3]"#, "abundance 0"),
        ] {
            match parse(&sketch_json(mins, abundances), 31, &md5) {
                SampleParse::Invalid(problem) => assert!(problem.contains(reason), "{}: {}", mins, problem),
                _ => panic!("expected an invalid sketch: {} {}", mins, abundances),
            }
        }
        assert_eq!(fracminhash_problem(&[1, 2], None, 2), None);
        assert_eq!(fracminhash_problem(&[], Some(&[]), 2), None);
    }
}
