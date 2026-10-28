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
