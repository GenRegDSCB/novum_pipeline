# data-prep — pipeline input bootstrap

Resolve + download bacterial assemblies (NCBI) — **CDS plus the non-coding
`rna_from_genomic` set, merged into one target FASTA** — and build miRNA query
FASTAs (miRBase), then surgically write the `queries:`/`targets:` blocks of
`Config/config.yaml` (every other param and comment is preserved). It does **not**
run Snakemake — you keep manual control.

## Usage

    python data-prep/prepare_inputs.py --pairs pairs.txt              # full run
    python data-prep/prepare_inputs.py --pairs pairs.txt --dry-run    # preview, no writes/downloads
    python data-prep/prepare_inputs.py Veillonella_parvula:hsa-miR-200b-3p   # inline pair(s)

Input lines are **taxon-first**; the separator is a tab, a `:`, or spaces, and
case is flexible:

    Veillonella_parvula     hsa-miR-200b-3p
    Dysosmobacter_welbionis:hsa-mir-200c-3p
    fusobacterium_nucleatum hsa-mir-515-5p

The same taxon on multiple lines collects its miRNAs into one query FASTA.
A bare miRNA name (no species prefix) is tried as `hsa-`, then `mmu-`
(`--prefix-chain` to change). Unresolved miRNAs / taxa are skipped and recorded
in `prep_report.tsv` rather than aborting the run.

### Pinning one assembly

The first field may be an **NCBI assembly accession** instead of a name, which
fetches that exact assembly rather than whichever one the resolver ranks first
for the species — the way to reach a paper's strain when it is not the species
reference:

    GCF_000153625.3         hsa-mir-515-5p    # F. nucleatum ATCC 10953
    GCF_000020225.1:hsa-miR-200b-3p           # A. muciniphila ATCC BAA-835

Both forms mix freely in one file. No mode flag is needed: an accession is
recognised by its shape — `GCF_`/`GCA_`, nine digits, a version suffix — which no
organism name or taxid can take. GenBank (`GCA_`) accessions are accepted too,
since a literature strain does not always have a RefSeq assembly.

Two consequences worth knowing:

- **A pin fails closed.** A typo resolves to nothing and is reported as
  `bad_accession`, distinct from `no_assembly`; it never falls back to a name
  search, which would quietly fetch some *other* assembly into the pinned tree.
- **The accession becomes the sample key**, so it names the query FASTA
  (`Data/Raw/RNAs/GCF_000153625.3.fa`), the `queries:`/`targets:` keys, and the
  `Data/Results/GCF_000153625.3/` output directory. Pinning buys exactness and
  costs readable paths; re-pinning to a new assembly version renames the tree.

## Options

    --config PATH        config to edit          [Config/config.yaml]
    --rnas-dir DIR       query FASTA output       [Data/Raw/RNAs]
    --genomes-dir DIR    genome output            [Data/Raw/genomes]
    --mirnas-dir DIR     miRBase mature.fa cache  [Data/Raw/miRNAs]
    --prefix-chain CSV   species fallback order   [hsa,mmu]
    --report PATH        per-(taxon,miRNA) TSV    [data-prep/prep_report.tsv]
    --force              re-download even if outputs exist
    --dry-run            resolve + print the blocks; no genome download, no config edit

## Target FASTAs

Each taxon gets **one** merged file in `--genomes-dir`:

    <accession>_<assembly>_cds_rna_from_genomic.fna    CDS + ncRNA/rRNA/tRNA
    <accession>_<assembly>_cds_from_genomic.fna        CDS only (RNA unavailable)

The two halves come from different NCBI services: CDS from the Datasets v2 API,
which does **not** publish `rna_from_genomic.fna` for prokaryotes, and the RNA
set from the FTP mirror at a URL derived from the accession and assembly name.
Neither half is kept separately — they are streamed into one file.

Only the **merged** name counts as a cached download. A CDS-only file is a
degraded result, so a later run retries it; that is also what migrates a
pre-merge `genomes/` directory without needing `--force`.

## Notes

- A taxon reaches `config.yaml` only if its genome resolved **and** ≥1 miRNA
  resolved (keeps the written config runnable).
- Assembly selection prefers a RefSeq *reference* genome, then *representative*,
  then the latest/most-complete assembly — so long-tail taxa with no designated
  reference still resolve. An accession in the taxon field skips this ranking
  entirely and takes no selection filters, so a `GCA_` pin is not rejected.
- `GENOME_STATUS` separates four outcomes that all look like "no genome":
  `lookup_error` (the NCBI call never completed), `no_assembly` (NCBI has no
  genome for that organism), `bad_accession` (a pin naming an assembly that does
  not exist), and `download_failed` (resolved, but the payload did not arrive —
  what an unannotated `GCA_` pin with no CDS lands as, and what every accession
  in a CDS batch of up to 200 lands as when that batch's download fails or its
  NCBI package is corrupt; the run carries on and still writes `config.yaml`
  and `prep_report.tsv`).
- A trailing `# comment` on a pairs line is stripped, so a pinned line can say
  which strain an opaque accession is.
- `prep_report.tsv` carries an `RNA_STATUS` column: `merged`, `cached`,
  `missing` (NCBI has no RNA file — 404), `error` (fetch failed; the URL is
  printed), or `-`. A missing RNA set never drops a taxon — it degrades to
  CDS-only and is still written to `config.yaml`.
- `NCBI_API_KEY` / `NCBI_API_EMAIL` (optional) raise the NCBI rate limit.

See `docs/superpowers/specs/2026-06-30-data-prep-design.md` for the full design,
and `docs/superpowers/specs/2026-07-23-cds-rna-merge-design.md` for the CDS+RNA merge.
