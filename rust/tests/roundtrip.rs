//! Build random indexes and compare the query results with a brute-force implementation of the
//! legacy YACHT semantics (exclusive hashes relative to the nontrivial genomes only).

use _yacht_core::builder::{build_index, build_index_opts};
use _yacht_core::query::Index;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::PathBuf;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // splitmix64
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn write_sig(path: &PathBuf, name: Option<&str>, ksize: u64, hashes: &[u64], gz: bool) {
    let mut sorted = hashes.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mins: Vec<String> = sorted.iter().map(|h| h.to_string()).collect();
    let name_field = match name {
        Some(n) => format!("\"name\":{},", serde_json::to_string(n).unwrap()),
        None => String::new(),
    };
    let json = format!(
        "[{{\"class\":\"sourmash_signature\",\"email\":\"\",\"hash_function\":\"0.murmur64\",\"filename\":\"x.fa\",{}\"license\":\"CC0\",\"signatures\":[{{\"num\":0,\"ksize\":21,\"seed\":42,\"max_hash\":18446744073709552,\"mins\":[1,2,3],\"md5sum\":\"x\",\"molecule\":\"DNA\"}},{{\"num\":0,\"ksize\":{},\"seed\":42,\"max_hash\":18446744073709552,\"mins\":[{}],\"abundances\":[{}],\"md5sum\":\"x\",\"molecule\":\"DNA\"}}],\"version\":0.4}}]",
        name_field,
        ksize,
        mins.join(","),
        vec!["1"; mins.len()].join(",")
    );
    let mut f = std::fs::File::create(path).unwrap();
    if gz {
        let mut e = flate2::write::GzEncoder::new(f, flate2::Compression::fast());
        e.write_all(json.as_bytes()).unwrap();
        e.finish().unwrap();
    } else {
        f.write_all(json.as_bytes()).unwrap();
    }
}

fn brute(rows: &[Vec<u64>], sample: &HashSet<u64>, in_n: &[bool]) -> (Vec<u64>, Vec<u64>) {
    let mut cnt: HashMap<u64, u32> = HashMap::new();
    for (r, hs) in rows.iter().enumerate() {
        if in_n[r] {
            for h in hs.iter().collect::<HashSet<_>>() {
                *cnt.entry(*h).or_default() += 1;
            }
        }
    }
    let mut excl = vec![0u64; rows.len()];
    let mut matches = vec![0u64; rows.len()];
    for (r, hs) in rows.iter().enumerate() {
        if !in_n[r] {
            continue;
        }
        for h in hs.iter().collect::<HashSet<_>>() {
            if cnt[h] == 1 {
                excl[r] += 1;
                if sample.contains(h) {
                    matches[r] += 1;
                }
            }
        }
    }
    (excl, matches)
}

#[test]
fn random_roundtrips() {
    let dir = std::env::temp_dir().join(format!("yacht_core_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut rng = Rng(12345);
    for case in 0..300u64 {
        let n_rows = 1 + rng.below(40) as usize;
        // small hash universe => lots of sharing; occasionally huge hashes to exercise key64 layouts
        let universe = match case % 4 {
            0 => 50,
            1 => 2000,
            2 => 1 << 40,
            _ => u64::MAX,
        };
        let pool: Vec<u64> = (0..(n_rows * 30).max(10)).map(|_| rng.below(universe)).collect();
        let mut rows: Vec<Vec<u64>> = Vec::new();
        let mut paths = Vec::new();
        for r in 0..n_rows {
            let k = rng.below(60) as usize;
            let mut hs: Vec<u64> = (0..k).map(|_| pool[rng.below(pool.len() as u64) as usize]).collect();
            if r > 0 && rng.below(10) == 0 {
                hs = rows[r - 1].clone(); // duplicate genome
            }
            let p = dir.join(format!("c{}_r{}.sig", case, r));
            let name = match rng.below(5) {
                0 => None,
                1 => Some("dup name".to_string()),
                _ => Some(format!("genome, \"{}\" é", r)),
            };
            write_sig(&p, name.as_deref(), 31, &hs, rng.below(2) == 0);
            hs.sort_unstable();
            hs.dedup();
            rows.push(hs);
            paths.push(p);
        }
        let out = dir.join(format!("c{}.yidx", case));
        let digest = [case as u8; 32];
        let stats = build_index_opts(&paths, 31, 1000, digest, &out, case % 3 == 0).unwrap();
        assert_eq!(Index::open(&out).unwrap().hdr.dir64(), case % 3 == 0);
        assert_eq!(stats.n_rows as usize, n_rows);
        let idx = Index::open(&out).unwrap();
        assert_eq!(idx.hdr.digest, digest);
        for (r, hs) in rows.iter().enumerate() {
            assert_eq!(idx.sizes()[r] as usize, hs.len());
        }
        for _ in 0..5 {
            let sample: Vec<u64> = (0..rng.below(80)).map(|_| {
                if rng.below(3) == 0 { rng.below(universe) } else { pool[rng.below(pool.len() as u64) as usize] }
            }).collect();
            let sample_set: HashSet<u64> = sample.iter().copied().collect();
            let q = idx.query(&sample);
            let marked = idx.marked_rows(&q);
            let expect_marked: Vec<u32> = (0..n_rows as u32)
                .filter(|&r| rows[r as usize].iter().any(|h| sample_set.contains(h)))
                .collect();
            assert_eq!(marked, expect_marked, "case {}", case);
            // N = marked rows plus a few random extra rows (models the by-name extension)
            let mut in_n = vec![false; n_rows];
            for &r in &marked {
                in_n[r as usize] = true;
            }
            if rng.below(2) == 0 {
                for _ in 0..rng.below(4) {
                    in_n[rng.below(n_rows as u64) as usize] = true;
                }
            }
            let got = idx.exclusive_counts(&q, &in_n).unwrap();
            let want = brute(&rows, &sample_set, &in_n);
            assert_eq!(got, want, "case {}", case);

            // SuperYACHT: per-row sample abundances of the sample hashes in the row (sorted)
            let mut sorted: Vec<u64> = sample_set.iter().copied().collect();
            sorted.sort_unstable();
            let abunds: Vec<u64> = sorted.iter().map(|_| 1 + rng.below(7)).collect();
            let qh = idx.query_with_hits(&sorted).unwrap();
            assert_eq!(idx.exclusive_counts(&qh, &in_n).unwrap(), want, "case {}", case);
            let sel: Vec<u32> = (0..n_rows as u32).filter(|&r| in_n[r as usize]).collect();
            let (offs, vals) = idx.coverage_values(&qh, &abunds, &sel).unwrap();
            assert_eq!(offs.len(), sel.len() + 1);
            for (i, &r) in sel.iter().enumerate() {
                let mut expect: Vec<u64> = rows[r as usize]
                    .iter()
                    .filter_map(|h| sorted.binary_search(h).ok().map(|p| abunds[p]))
                    .collect();
                expect.sort_unstable();
                assert_eq!(&vals[offs[i] as usize..offs[i + 1] as usize], &expect[..], "case {} row {}", case, r);
            }
        }
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn coverage_values_needs_sorted_hashes_and_hits() {
    let dir = std::env::temp_dir().join(format!("yacht_core_cov_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("g.sig");
    write_sig(&p, Some("g"), 31, &[5, 9, 12], false);
    let out = dir.join("i.yidx");
    build_index_opts(&[p], 31, 1000, [0u8; 32], &out, false).unwrap();
    let idx = Index::open(&out).unwrap();
    assert!(idx.query_with_hits(&[9, 5]).is_err());
    let q = idx.query(&[5, 9]); // no hits recorded
    assert!(idx.coverage_values(&q, &[1, 2], &[0]).is_err());
    let q = idx.query_with_hits(&[5, 9, 40]).unwrap();
    assert!(idx.coverage_values(&q, &[1, 2, 3], &[1]).is_err()); // row out of range
    assert_eq!(idx.coverage_values(&q, &[4, 2, 3], &[0]).unwrap(), (vec![0, 2], vec![2, 4]));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn rejects_ambiguous_ksize() {
    let dir = std::env::temp_dir().join(format!("yacht_core_test_amb_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("a.sig");
    write_sig(&p, Some("a"), 21, &[5, 6, 7], false); // only k=21 sketches -> two of them, none at 31
    let err = build_index(&[p], 31, 1000, [0; 32], &dir.join("a.yidx")).err().unwrap();
    assert!(err.contains("Expected exactly one signature with ksize 31"), "{}", err);
    std::fs::remove_dir_all(&dir).ok();
}
