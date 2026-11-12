# Retrieving the genomes behind a YACHT result

YACHT reports each genome by the name of its sketch (`organism_name` in `result.xlsx` and `result_all.txt`). In the
pretrained databases, the first word of that name is a stable identifier of the genome:

| identifier looks like | what it is | where the genome is |
|---|---|---|
| `GCA_000005845.2`, `GCF_000005845.2` | versioned INSDC/NCBI assembly accession (`GCA_` GenBank, `GCF_` RefSeq) | NCBI (below) |
| `MOTU40_133241`, `SPIREOTU_03021551`, `MGYG000535630`, `GEMOTU_00001`, ... | GlobDB r232 genome ID (only in `globdb-r232.k31` and `yacht-union-globdb-2026.k31`) | GlobDB r232 release (below) |

Each pretrained database also contains `<prefix>_lineages.csv` (taxonomy of every genome) and
`<prefix>_removed_genomes.tsv`. The latter lists the genomes that `yacht train` left out because they are equivalent
(ANI at or above the database's threshold) to a genome in the database. When YACHT reports genome X, the genomes listed
next to X in that file are indistinguishable from X at this ANI threshold.

## Genomes with an NCBI assembly accession

The accession includes its version (`.2`). NCBI keeps every version, so the accession always identifies the exact
sequence that was sketched.

With the [NCBI Datasets command-line tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
(`conda install -c conda-forge ncbi-datasets-cli`):

```bash
# one genome
datasets download genome accession GCF_000005845.2 --include genome --filename GCF_000005845.2.zip
# many genomes: one accession per line in accessions.txt
datasets download genome accession --inputfile accessions.txt --include genome --filename genomes.zip
unzip genomes.zip   # ncbi_dataset/data/<accession>/<accession>_<assembly name>_genomic.fna
```

Without the tool, every assembly is on the NCBI FTP site under a path derived from its accession, e.g.
`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/` holds the folder
`GCF_000005845.2_ASM584v2/` with `GCF_000005845.2_ASM584v2_genomic.fna.gz`. The web page of an assembly is
`https://www.ncbi.nlm.nih.gov/datasets/genome/<accession>/`.

The bacterial and archaeal genomes of the GTDB-based databases are also in GTDB's own release archive
(`genomic_files_reps/gtdb_genomes_reps_r232.tar.gz` at https://data.gtdb.ecogenomic.org/releases/release232/232.0/,
192 GB), e.g. `gtdb_genomes_reps_r232/database/GCF/000/020/485/GCF_000020485.1_genomic.fna.gz`.

## Genomes with a GlobDB ID

GlobDB IDs are assigned by [GlobDB](https://globdb.org) and are stable within a GlobDB release; the YACHT databases use
release **r232**. The prefix of the ID names the source catalogue (e.g. `MOTU40_` mOTUs 4.0, `SPIREOTU_` SPIRE,
`MGYG` MGnify, `GEMOTU_` GEM; the list is in `globdb_r232_dataset_list.tsv`). The files of the release are on the
University of Vienna file share: https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/

| file | content |
|---|---|
| `globdb_r232_genome_fasta.tar.gz` (292 GB, md5 `c95ef3b3278b9f035e2089bd729b0b48`, see `globdb_r232_md5sum.txt`) | all 346,233 genomes, as `globdb_r232_genome_fasta/R232CHUNK<nnn>/<GlobDB ID>.fa.gz` |
| `globdb_r232_taxonomy.tsv.gz` | GlobDB ID and GTDB-style taxonomy (`d__...;s__...`) |
| `globdb_r232_dictionaries.tar.gz` | per source catalogue, the original genome IDs and their GlobDB IDs (see its `README`: most files have two columns without a header, original ID then GlobDB ID; the SPIRE and GEM files have a header and also the species cluster ID) |
| `globdb_r232_dataset_list.tsv` | the source catalogues, their papers and how many genomes each contributes |

The genome archive is a single compressed tar file, so a genome cannot be fetched on its own, but it can be extracted
while the archive is downloaded, without storing the archive (this still reads the whole 292 GB stream):

```bash
# one or a few genomes
curl -s https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/globdb_r232_genome_fasta.tar.gz \
  | tar -xzf - --wildcards '*/MOTU40_133241.fa.gz' '*/SPIREOTU_03021551.fa.gz'
# many genomes: one pattern per line, e.g. created with sed 's|.*|*/&.fa.gz|' ids.txt > patterns.txt
curl -s https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/globdb_r232_genome_fasta.tar.gz \
  | tar -xzf - --wildcards -T patterns.txt
```

If you need genomes repeatedly, download the archive once, check its md5 against `globdb_r232_md5sum.txt`, and extract
from the local copy in the same way.

To find a genome in its source catalogue (for example to get its sample metadata), look up its GlobDB ID in the
dictionary file of that catalogue. Catalogues whose IDs were already standardized (GTDB, NGDC, MGnify, CRBC, HRGM2,
MRGM) have no dictionary: their GlobDB ID is the catalogue's own ID.

In GlobDB, the GTDB genomes are named by their accession without version (e.g. `GCA_047199055`). The YACHT databases
name them with the versioned accession of GTDB R232 instead (e.g. `GCA_047199055.1`), so they can be retrieved from
NCBI as described above.

Please cite GlobDB (Speth et al., Bioinformatics Advances, 2025, https://doi.org/10.1093/bioadv/vbaf280) and the
source catalogues of the genomes you use (https://globdb.org/citations).
