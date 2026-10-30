"""
Data-card texts of the YACHT 2.0 pretrained databases. Writes db_sources.json (used by package_db.py).
Placeholders filled by package_db.py: {prefix} {name} {ani} {n_genomes} {n_trained} {n_removed} {ksize} {scaled}
{yacht_version} {n_no_lineage} {n_no_taxid} {ani_warning} {coarseness} {n_paired}.
"""
import json
import os

REPO_DOCS = "https://github.com/KoslickiLab/YACHT/blob/main/docs"

COMMON_TOP = """# YACHT pretrained database `{prefix}`

{title}

| | |
|---|---|
| genomes in the database | {n_genomes} |
| genome sketches it was trained on | {n_trained} ({n_removed} left out as equivalent to another genome at ANI >= {ani}) |
| ANI threshold | {ani} |
| k-mer size / scaled | {ksize} / {scaled} |
| trained with | YACHT {yacht_version} (`yacht train --ksize {ksize} --ani_thresh {ani}`) |
| license | {license} |
{ani_warning}
## Use

```bash
yacht download pretrained_ref_db --name {name} --ani_thresh {ani} --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer {ksize} --scaled {scaled} --outfile sample.sig.zip
yacht run --json {prefix}_pretrained/{prefix}_config.json --sample_file sample.sig.zip --outdir ./
```

A sample sketched with a smaller `--scaled` than {scaled} also works (it is downsampled to {scaled}).

## Sources

{sources}

## Genome identifiers

The `organism_name` of each genome in YACHT's output is its sketch name. Its first word is the genome's identifier:
{identifiers}
How to retrieve a genome from its identifier: {retrieve}

## Taxonomy

{taxonomy}

- `{prefix}_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  {n_no_lineage} genomes have no lineage.
- `{prefix}_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. {n_no_taxid} genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= {ani}). For each taxonomic rank:
the taxa of the {n_trained} genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the {n_paired} left-out genomes that have a genome of this database at ANI >= {ani}:

{coarseness}

## Files

| file | content |
|---|---|
| `{prefix}_config.json`, `{prefix}_processed_manifest.tsv`, `{prefix}_index.yidx` | the database `yacht run` uses |
| `{prefix}_lineages.csv`, `{prefix}_genome_to_taxid.tsv` | taxonomy (see above) |
| `{prefix}_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= {ani}; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= {ani} (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.
{notes}"""

# shown below the summary table of the ANI 0.80 databases ({containment} = 0.80^k)
ANI_080_WARNING = """
> **ANI 0.80 is a very coarse setting, much coarser than genus level.** At k = {ksize}, ANI 0.80 corresponds to a
> k-mer containment of only {containment}: two genomes are equivalent when about 1 in 1,000 hashes of the smaller one
> is shared (about 4 of the 4,000 hashes of a 4 Mbp genome at scaled 1000, a single hash of a 30 kb virus at scaled
> 50). Overlaps this small also occur between genomes of different families, orders and phyla. The table "How coarse
> is this database?" below shows how many taxa keep no genome of their own and how often a left-out genome is
> represented by a genome of another taxon. Use this database only for a coarse overview; for genus- or species-level
> results, use ANI 0.95 or higher.
"""

PROKARYOTE_NOTE = ""
VIRUS_NOTE = """
## Notes

- Viral genomes are small: at scaled {scaled}, a 30 kb genome has about {virus_hashes} hashes. The fewer hashes a genome
  has, the less evidence YACHT has for or against its presence.
"""
GENBANK_VIRUS_NOTE = """
- **Not every GenBank "viral" assembly is a single viral genome.** 92 assemblies of this set are larger than 1 Mbp,
  while the largest known viral genomes are about 2.5 Mbp: e.g. `GCA_024181985.1` "uncultured virus" (334 Mbp, a
  collection of metagenomic viral contigs) and `GCA_038430685.1` "Enterococcus phage Efs 1-1-1" (3 Mbp, with its
  host's genome). `yacht train` keeps the larger of two equivalent genomes, so such an assembly can represent many
  viruses (at ANI 0.95, `GCA_024181985.1` alone represents 5,042), and an assembly that contains a bacterial genome is
  reported present in samples with that bacterium. `ncbi-viruses-refseq-2025.01.k31-sc50` (the curated RefSeq genomes
  only) and `ictv-vmr-msl40.k31-sc50` avoid this.
"""
VMR_HOST_NOTE = """
- **A few ICTV exemplars are sketched with their host's genome.** For some exemplars (prophages, gene transfer
  agents) the VMR points to the GenBank assembly of the host that carries the virus, e.g. `GCA_000091085.2` "Chlamydia
  phage phiCPAR39" (the 1.25 Mbp *Chlamydia pneumoniae* AR39 genome) and `GCA_002951655.1` "Escherichia coli strain
  13P477T plasmid p1" (5.1 Mbp). Such an exemplar is reported present in samples that contain its host."""
EUK_NOTE = """
## Notes

- The assemblies NCBI flags as contaminated are left out (see Sources). Other GenBank eukaryotic assemblies can still
  contain contaminating (e.g. bacterial) sequence. With this eukaryote-only
  database, a genome reported present with many matches in a bacteria-rich sample may be detected through such
  sequence. In the union databases, k-mers shared with the bacterial genomes are not exclusive to the eukaryotic genome,
  so they do not count as evidence for it.
"""
UNION_NOTE = """
## Notes

- **Viruses** are the NCBI RefSeq viral genomes and the ICTV exemplars (VMR MSL40): curated genomes only. The GenBank
  viral assemblies are not included, because some of them are collections of many viral genomes or contain their host's
  genome, and `yacht train` keeps the larger of two equivalent genomes. Viruses are included at scaled 1000, where a
  30 kb genome has about 30 hashes and many small viruses have none (those are left out: they can never be detected).
  For viruses, the scaled-50 viral databases (`ncbi-viruses-refseq-2025.01.k31-sc50`, `ictv-vmr-msl40.k31-sc50`, with
  samples sketched at `--scaled 50`) are more sensitive.
- A few ICTV exemplars are sketches of the host assembly that carries the virus (prophages, gene transfer agents, e.g.
  `GCA_000091085.2` "Chlamydia phage phiCPAR39", the 1.25 Mbp *Chlamydia pneumoniae* AR39 genome). Since `yacht train`
  keeps the larger of two equivalent genomes, such an exemplar can represent its host's GTDB genome (a handful of
  genomes; see `{prefix}_removed_genomes.tsv`).
- **Eukaryotes:** the GenBank fungal and protozoan assemblies NCBI flags as contaminated are left out (listed in
  `fungi.ncbi-contaminated.csv` and `protozoa.ncbi-contaminated.csv`). Other eukaryotic assemblies can still contain
  bacterial sequence, and since `yacht train` keeps the larger of two equivalent genomes, such an assembly can represent
  a bacterial genome: `{prefix}_removed_genomes.tsv` lists, for each left-out genome, the genome that represents it.
- **Prophages:** `yacht train` treats two genomes as equivalent when the smaller one is contained in the larger one at
  ANI >= {ani}. A phage whose genome is (as a prophage) part of a larger genome in the database is therefore represented by
  that genome (see `{prefix}_removed_genomes.tsv`).
- Taxonomies differ between the parts of the database: GTDB for bacteria and archaea, NCBI for eukaryotes and RefSeq
  viruses, ICTV MSL40 for the VMR exemplars.
"""

GTDB_R232 = {
    "name": "GTDB R11-RS232 (release 232.0, 15 April 2026), species representative genomes (199,923)",
    "url": "https://data.gtdb.ecogenomic.org/releases/release232/232.0/",
    "files": "genomic_files_reps/gtdb_genomes_reps_r232.tar.gz (md5 a63a586186871a429bf12ea5bc7be0db), bac120/ar53 metadata and taxonomy",
    "license": "CC BY-SA 4.0",
    "cite": "Parks et al., GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy, Nucleic Acids Res. 50:D785 (2022), https://doi.org/10.1093/nar/gkab776",
    "how": "Sketched with `sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000` (sourmash_plugin_branchwater), named `<accession> <ncbi_organism_name>` from the GTDB metadata (the naming of sourmash's GTDB databases); k=31 used for training.",
}
GTDB_RS226 = {
    "name": "GTDB R10-RS226 species representative genomes, prepared by the sourmash team (143,384 of 143,614 genomes)",
    "url": "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/ (gtdb-rs226-reps.k31.sig.zip, gtdb-rs226-reps.lineages.csv)",
    "files": "",
    "license": "GTDB data: CC BY-SA 4.0; sourmash sketches: CC0",
    "cite": "Parks et al., Nucleic Acids Res. 50:D785 (2022), https://doi.org/10.1093/nar/gkab776; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830",
    "how": "Used as published (k=31, scaled=1000).",
}
GLOBDB_R232 = {
    "name": "GlobDB release 232 (June 2026): 346,233 species-representative genomes, the GTDB R232 species representatives plus species-level MAGs from 25 further catalogues",
    "url": "https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/",
    "files": "globdb_r232_genome_fasta.tar.gz (md5 c95ef3b3278b9f035e2089bd729b0b48), globdb_r232_taxonomy.tsv.gz, globdb_r232_dictionaries.tar.gz",
    "license": "CC BY-SA 4.0 (GlobDB propagates the licenses of the underlying data sources)",
    "cite": "Speth et al., GlobDB: a comprehensive species-dereplicated microbial genome resource, Bioinformatics Advances (2025), https://doi.org/10.1093/bioadv/vbaf280, and the source catalogues listed at https://globdb.org/citations",
    "how": "The GTDB R232 part reuses the GTDB R232 sketches above (same genomes, versioned accessions from the GTDB metadata); the other genomes were sketched with `sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000`, named `<GlobDB ID> <GlobDB species>`.",
}
FUNGI = {
    "name": "GenBank fungal genomes, 8 April 2025 (21,135 genomes), prepared by the sourmash team",
    "url": "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/ (genbank-20250408-fungi-k31.zip, lineages.fungi.csv)",
    "files": "",
    "license": "NCBI GenBank (no NCBI restrictions on use); sourmash sketches: CC0",
    "cite": "NCBI GenBank, https://www.ncbi.nlm.nih.gov/genbank/; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830",
    "how": "Used as published (k=31, scaled=1000), without the assemblies NCBI flags as contaminated (`excluded_from_refseq` "
           "contains \"contaminated\" in the NCBI GenBank assembly summary, current or historical, downloaded 30 October 2026), "
           "left out with `sourmash sig extract --picklist ...:exclude`; they are listed in `fungi.ncbi-contaminated.csv` in this zip.",
}
PROTOZOA = dict(FUNGI, name="GenBank protozoan genomes, 8 April 2025 (2,392 genomes), prepared by the sourmash team",
                url="https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/ (genbank-20250408-protozoa-k31.zip, lineages.protozoa.csv)",
                how=FUNGI["how"].replace("fungi.ncbi-contaminated.csv", "protozoa.ncbi-contaminated.csv"))
VMR = {
    "name": "ICTV Virus Metadata Resource, MSL40 v1 (2025): 17,765 exemplar virus genomes, prepared by the sourmash team (ICTV has since released MSL41, July 2026)",
    "url": "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/2025-ictv-vmr-msl40/ (vmr_MSL40_v1.dna-k31.sig.zip, vmr_MSL40_v1.lineages.csv)",
    "files": "",
    "license": "ICTV VMR: CC BY 4.0 (https://doi.org/10.5281/zenodo.15042310); sourmash sketches: CC0",
    "cite": "ICTV Virus Metadata Resource, https://ictv.global/vmr, https://doi.org/10.5281/zenodo.15042310",
    "how": "Used as published (k=31, scaled=50); for the unions downsampled to scaled 1000 with `sourmash sig downsample`.",
}
NCBI_VIRUSES = {
    "name": "NCBI Virus genome assemblies (NCBI:txid10239), January 2025 (230,973 assemblies), prepared by the sourmash team",
    "url": "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/ncbi-viruses-2025.01/ (ncbi-viruses.dna.k=31.scaled=50.sig.zip, ncbi-viruses.lineages.csv)",
    "files": "",
    "license": "NCBI (no NCBI restrictions on use); sourmash sketches: CC0",
    "cite": "NCBI Virus, https://www.ncbi.nlm.nih.gov/labs/virus/; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830",
    "how": "Used as published (k=31, scaled=50).",
}
NCBI_VIRUSES_REFSEQ = dict(
    NCBI_VIRUSES,
    name="NCBI RefSeq viral genomes of the NCBI Virus assemblies of January 2025 (16,280 of 230,973 assemblies, the `GCF_` accessions), prepared by the sourmash team",
    how="The RefSeq (`GCF_`) sketches, selected with `sourmash sig extract --picklist` (k=31, scaled=50); for the unions "
        "downsampled to scaled 1000 with `sourmash sig downsample`.",
)

INSDC_IDS = (
    "a versioned INSDC/NCBI assembly accession (`GCA_...` for GenBank, `GCF_...` for RefSeq), e.g. `GCF_000005845.2`."
)
INSDC_RETRIEVE = (
    "`datasets download genome accession <accession> --include genome` (NCBI Datasets), or "
    "https://www.ncbi.nlm.nih.gov/datasets/genome/<accession>/. Details: "
    f"{REPO_DOCS}/retrieving_genomes.md"
)
GLOBDB_IDS = (
    "a versioned INSDC/NCBI assembly accession (`GCA_...`/`GCF_...`) for the GTDB R232 genomes, and the GlobDB genome ID "
    "(e.g. `MOTU40_133241`, `SPIREOTU_03021551`, `MGYG000535630`; the prefix names the source catalogue) for the "
    "other genomes. GlobDB IDs are stable within a GlobDB release; `globdb_r232_dictionaries.tar.gz` maps them to the "
    "IDs of the source catalogues."
)
GLOBDB_RETRIEVE = (
    "for an accession, NCBI Datasets (`datasets download genome accession <accession>`); for a GlobDB ID, the genome is "
    "`globdb_r232_genome_fasta/<chunk>/<GlobDB ID>.fa.gz` in `globdb_r232_genome_fasta.tar.gz` of the GlobDB r232 "
    f"release (it can be extracted from the download stream without storing the whole archive). Details: {REPO_DOCS}/retrieving_genomes.md"
)


def sources_md(*sources):
    out = []
    for s in sources:
        lines = [f"- **{s['name']}**", f"  - {s['url']}"]
        if s["files"]:
            lines.append(f"  - files: {s['files']}")
        lines += [f"  - license: {s['license']}", f"  - please cite: {s['cite']}", f"  - preparation: {s['how']}"]
        out.append("\n".join(lines))
    return "\n".join(out)


def card(title, license, sources, identifiers, retrieve, taxonomy, notes):
    text = COMMON_TOP.replace("{notes}", notes)
    for k, v in (("{title}", title), ("{license}", license), ("{sources}", sources_md(*sources)),
                 ("{identifiers}", identifiers), ("{retrieve}", retrieve), ("{taxonomy}", taxonomy)):
        text = text.replace(k, v)
    return text


GTDB_TAX = "GTDB R232 taxonomy (domain to species, `d__...;s__...`), from the GTDB metadata; NCBI taxids from the GTDB metadata (`ncbi_taxid`)."
CARDS = {
    "gtdb-rs232-reps": card(
        "GTDB R232 species representative genomes (bacteria and archaea).", "CC BY-SA 4.0 (as GTDB)", [GTDB_R232],
        INSDC_IDS, INSDC_RETRIEVE, GTDB_TAX, PROKARYOTE_NOTE),
    "gtdb-rs226-reps": card(
        "GTDB R226 species representative genomes (bacteria and archaea).", "CC BY-SA 4.0 (as GTDB)", [GTDB_RS226],
        INSDC_IDS, INSDC_RETRIEVE, "GTDB R226 taxonomy (`gtdb-rs226-reps.lineages.csv` of the sourmash team). No NCBI taxids are included.",
        PROKARYOTE_NOTE),
    "globdb-r232": card(
        "GlobDB r232: the most comprehensive species-level set of bacterial and archaeal genomes (GTDB R232 plus MAG catalogues).",
        "CC BY-SA 4.0 (as GlobDB and GTDB)", [GTDB_R232, GLOBDB_R232], GLOBDB_IDS, GLOBDB_RETRIEVE,
        "GlobDB r232 taxonomy (GTDB R232 taxonomy, assigned with GTDB-Tk for the non-GTDB genomes). NCBI taxids only for the GTDB R232 genomes (from the GTDB metadata).",
        PROKARYOTE_NOTE),
    "genbank-2025.04-fungi": card(
        "GenBank fungal genomes (April 2025).", "NCBI data; sketches CC0", [FUNGI], INSDC_IDS.replace("(`GCA_...` for GenBank, `GCF_...` for RefSeq)", "(`GCA_...`)"),
        INSDC_RETRIEVE, "NCBI taxonomy and taxids (`lineages.fungi.csv` of the sourmash team).", EUK_NOTE),
    "genbank-2025.04-protozoa": card(
        "GenBank protozoan genomes (April 2025).", "NCBI data; sketches CC0", [PROTOZOA], INSDC_IDS.replace("(`GCA_...` for GenBank, `GCF_...` for RefSeq)", "(`GCA_...`)"),
        INSDC_RETRIEVE, "NCBI taxonomy and taxids (`lineages.protozoa.csv` of the sourmash team).", EUK_NOTE),
    "ictv-vmr-msl40": card(
        "ICTV exemplar virus genomes (Virus Metadata Resource MSL40), at scaled 50.", "CC BY 4.0 (as the ICTV VMR)", [VMR],
        INSDC_IDS, INSDC_RETRIEVE,
        "ICTV taxonomy (realm to species, MSL40; `vmr_MSL40_v1.lineages.csv` of the sourmash team). NCBI taxids of the assemblies from `ncbi-viruses.lineages.csv` (NCBI Virus, January 2025, same versioned accession).",
        VIRUS_NOTE.replace("{virus_hashes}", "600") + VMR_HOST_NOTE),
    "ncbi-viruses-2025.01": card(
        "NCBI viral genome assemblies (January 2025), all GenBank and RefSeq assemblies, at scaled 50. See the note on "
        "GenBank assemblies that are not single viral genomes; `ncbi-viruses-refseq-2025.01.k31-sc50` holds the curated "
        "RefSeq genomes only.", "NCBI data; sketches CC0", [NCBI_VIRUSES], INSDC_IDS,
        INSDC_RETRIEVE, "NCBI taxonomy and taxids (`ncbi-viruses.lineages.csv` of the sourmash team).",
        VIRUS_NOTE.replace("{virus_hashes}", "600") + GENBANK_VIRUS_NOTE),
    "ncbi-viruses-refseq-2025.01": card(
        "NCBI RefSeq viral genomes (January 2025), at scaled 50: the curated subset of the NCBI viral assemblies.",
        "NCBI data; sketches CC0", [NCBI_VIRUSES_REFSEQ], INSDC_IDS.replace("(`GCA_...` for GenBank, `GCF_...` for RefSeq)", "(`GCF_...`, RefSeq)"),
        INSDC_RETRIEVE, "NCBI taxonomy and taxids (`ncbi-viruses.lineages.csv` of the sourmash team).",
        VIRUS_NOTE.replace("{virus_hashes}", "600")),
    "yacht-union-insdc-2026": card(
        "Bacteria, archaea, fungi, protozoa and viruses, with INSDC/NCBI assembly accessions only: GTDB R232 species "
        "representatives, GenBank fungi and protozoa (April 2025, without the assemblies NCBI flags as contaminated), "
        "NCBI RefSeq viral genomes (January 2025) and the ICTV exemplar viruses (VMR MSL40). Every genome can be "
        "retrieved from NCBI by its accession; this is the easiest database to go from a result to the genome.",
        "CC BY-SA 4.0 (as GTDB)", [GTDB_R232, FUNGI, PROTOZOA, NCBI_VIRUSES_REFSEQ, VMR], INSDC_IDS, INSDC_RETRIEVE,
        "GTDB R232 taxonomy for bacteria and archaea (`d__...`), NCBI taxonomy for eukaryotes and RefSeq viruses, ICTV "
        "MSL40 taxonomy for the VMR exemplars (in the same columns, `superkingdom` = `Viruses`). NCBI taxids for all "
        "genomes that have one (VMR exemplars: from the NCBI viral lineages of the same accession).", UNION_NOTE),
    "yacht-union-globdb-2026": card(
        "Bacteria, archaea, fungi, protozoa and viruses, the most comprehensive YACHT database: GlobDB r232 (GTDB R232 "
        "species representatives plus MAG catalogues), GenBank fungi and protozoa (April 2025, without the assemblies "
        "NCBI flags as contaminated), NCBI RefSeq viral genomes (January 2025) and the ICTV exemplar viruses (VMR MSL40). "
        "Genomes with a GlobDB ID are retrieved from the GlobDB release (see below), which is more involved than for the "
        "accession-only database `yacht-union-insdc-2026.k31`.",
        "CC BY-SA 4.0 (as GlobDB and GTDB)", [GTDB_R232, GLOBDB_R232, FUNGI, PROTOZOA, NCBI_VIRUSES_REFSEQ, VMR], GLOBDB_IDS,
        GLOBDB_RETRIEVE,
        "GlobDB r232 (GTDB R232) taxonomy for bacteria and archaea, NCBI taxonomy for eukaryotes and RefSeq viruses, ICTV "
        "MSL40 taxonomy for the VMR exemplars (in the same columns, `superkingdom` = `Viruses`). NCBI taxids for all "
        "genomes that have one, i.e. not for the GlobDB-only MAGs.", UNION_NOTE),
}

if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "db_sources.json")
    with open(out, "w") as f:
        json.dump({**{k: {"readme": v} for k, v in CARDS.items()}, "ani_0.80_warning": ANI_080_WARNING}, f, indent=1)
    print(out)
