# The YACHT run index (YACHT ≥ 2.0)

`yacht run` used to compare the sample against every reference genome on every run. It ran
`sourmash scripts multisearch` over all genome signatures, then re-read every overlapping genome
signature twice in Python to find the k-mers exclusive to each genome. For the GTDB databases,
that meant parsing tens of thousands of JSON files for every sample.

Since version 2.0, `yacht train` writes a binary index of the selected genomes, `<prefix>_index.yidx`,
next to `<prefix>_config.json`. `yacht run` memory-maps this index and computes the same results
from it, single-threaded, in a few seconds. The work is done by a small Rust extension module,
`yacht._yacht_core` (sources in `rust/`).

## What `yacht run` computes

For a sample with sketch hashes `S` and the genomes of the processed manifest (rows `g`, with
sketch hashes `H_g`):

1. **Overlapping genomes.** These are the rows with `|H_g ∩ S| > 0`. This is exactly what
   `multisearch -t 0` reported. The genomes analysed (`N`) are all manifest rows whose
   `organism_name` equals the name of an overlapping genome. This join by name is kept as before:
   the names are round-tripped through CSV and `pandas.read_csv` exactly like the multisearch
   output was.
2. **Exclusive k-mers.** For every `g` in `N`:
   - `E_g` is the set of hashes in `H_g` that occur in no other genome of `N`;
   - `num_exclusive_kmers_to_genome = |E_g|`;
   - `num_matches = |E_g ∩ S|`.

   Exclusivity is relative to `N`, not to the whole database, as before.
3. **Hypothesis tests** for every genome and every `min_coverage` value. These use the same scipy
   functions as before, vectorized over genomes, and give bit-for-bit identical numbers (see
   `hyp_test_vectorized`).

## Index contents

Everything is stored in manifest row order:

- the number of hashes of each genome sketch (checked against the manifest);
- the number of hashes found in only that genome;
- the genome signature names (what multisearch reported as `match_name`);
- a hash table mapping each distinct hash to either the single row containing it, or a *shared
  posting set*;
- the distinct shared posting sets, i.e. the sets of rows that share a hash. Each set is stored
  once, together with the number of hashes that have exactly that set.

With this, step 1 is a lookup of every sample hash. Step 2 is one pass over the shared posting
sets:

- a hash is exclusive to `g` in `N` exactly when its posting set intersected with `N` is `{g}`;
- so `|E_g|` is the number of hashes unique to `g` in the whole database, plus the weights of the
  shared sets whose intersection with `N` is `{g}`;
- `|E_g ∩ S|` is the same sum restricted to the hashes found in the sample.

Duplicate genomes (identical md5sums or identical sketches) are handled naturally, because postings
are per manifest row.

## File format (version 1)

Little-endian; every section starts on a 64-byte boundary.

| Section | Type | Content |
|---|---|---|
| header | 256 bytes | magic `YACHTIDX`, format version, flags, ksize, scaled, counts, SHA-256 of the manifest md5sum column, section offsets |
| SIZES | `u32[n_rows]` | hashes per genome sketch |
| UNIQUE | `u32[n_rows]` | hashes found only in this genome |
| NAME_OFFSETS / NAME_BYTES | `u64[n_rows+1]` / UTF-8 | genome signature names |
| DIR | `u32` (or `u64`) `[2^bucket_bits + 1]` | bucket directory (bucket = high bits of the hash) |
| ENTRIES | `(u32 key, u32 val)` or `(u64 key, u32 val, pad)` | low bits of each distinct hash, sorted within its bucket; `val < n_rows` is a row, otherwise a shared set id |
| SET_OFFSETS / SET_MEMBERS / SET_WEIGHTS | `u64` / `u32` / `u32` | shared posting sets |

For the GTDB r214 representative genomes (k=31, scaled=1000), the index has about 250 million
hashes and takes about 2 GB:

- the 0.95 ANI database (78,861 genomes) is 2.0 GB;
- the 0.9995 ANI database (85,205 genomes) is 2.2 GB.

It is built in under a minute with ~3 GB of RAM. The file is written to a temporary name and
renamed into place, so concurrent readers never see a partial file.

Because the index is memory-mapped, the operating system keeps a single copy in its page cache that
is shared by all concurrently running `yacht run` processes.

## Staleness and legacy databases

The header stores ksize, scaled and a SHA-256 digest of the manifest's md5sum column (in row order).
`yacht run` checks these, and the per-genome sketch sizes, against the config and manifest.

If no valid index exists, `yacht run` prints instructions and exits with status 1. This is the case
for databases created with YACHT < 2.0, or databases whose manifest changed. The index can be
created in place with:

```bash
python -m yacht.yacht_index --json <prefix>_config.json [--force]
```

`yacht download pretrained_ref_db` does this automatically when a downloaded database has no
index. It uses an index already present in the downloaded folder as is, so databases zipped
together with their index work without rebuilding.

## Output writing

`result_all.txt` is written by `yacht.utils.write_dataframe_tsv`, which produces the same bytes as
`DataFrame.to_csv(sep="\t", index=False)`:

- floats are printed like Python's `repr` (shortest round-trip digits, ties to even, via Ryu);
- text fields use `csv.QUOTE_MINIMAL` rules.

It falls back to pandas for any column type it does not handle. `result.xlsx` is still written by
pandas/openpyxl.
