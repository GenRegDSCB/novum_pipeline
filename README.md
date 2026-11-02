# Novum Pipeline

Novum Pipeline predicts which bacterial genes a set of miRNAs may bind. It runs two RNA–RNA interaction tools, **RNAhybrid** and **IntaRNA**, over the same miRNAs and genomes, and reports each tool's predictions plus the **consensus**: the miRNA–gene pairs both tools predict.

A [Snakemake](https://snakemake.github.io/) workflow takes you from FASTA files to annotated tables and readable reports with one command.

Start with [Installation](#installation) and the [Quick start](#quick-start). [How it works](#how-it-works) describes each step.

## Installation

There are two ways to install. The **simple way** is Docker: one image holds Snakemake and every tool, and runs the same on Linux, WSL and Macs. The **direct way** installs everything with conda, and is the one the config editor needs. Commands assume bash.

### Simple way: Docker

1. **Install Docker:** [Docker Desktop](https://docs.docker.com/desktop/) on a Mac, [Docker Engine](https://docs.docker.com/engine/install/) on Linux. In WSL, use Docker Engine inside WSL or Docker Desktop with WSL integration.

2. **Pull the image:**

   ```bash
   docker pull --platform linux/amd64 albertosola/novum_pipeline:v1.0.0
   ```

   The image is built for Intel/AMD processors, because RNAhybrid has no ARM build. Apple Silicon Macs run it through Rosetta automatically. It is slower than on Intel, as the direct way is on those Macs.

3. **Make a project folder.** A clone of this repository works as one. Alternatively, start from an empty folder and let the image write a starter `Config/config.yaml`:

   ```bash
   mkdir my_project && cd my_project
   docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$PWD":/work albertosola/novum_pipeline:v1.0.0 init
   ```

   The container sees only the project folder, so keep your FASTAs inside it, for example under `Data/Raw/`.

4. **Run**, from the project folder:

   ```bash
   docker run --rm -it --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$PWD":/work albertosola/novum_pipeline:v1.0.0 --cores 8
   ```

   Everything after the image name goes to Snakemake. Add `-n` for a dry run.

   To run one sample, name its four final files ([Outputs](#outputs)). With `results_dir: Data/Results/`, which is what `init` writes:

   ```bash
   docker run --rm -it --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$PWD":/work albertosola/novum_pipeline:v1.0.0 --cores 8 \
     Data/Results/<sample>/{w_calibration,wo_calibration}/rnahybrid_enhanced.txt \
     Data/Results/<sample>/w_accessibility/intarna_enhanced.txt \
     Data/Results/<sample>/consensus/consensus_enhanced.txt
   ```

   To run the [input bootstrap](#preparing-inputs-automatically) instead, put `prepare-inputs` after the image name, e.g. `… v1.0.0 prepare-inputs --pairs pairs.txt`.

   For long runs, replace `--rm -it` with `-d --name novum`, then follow the run with `docker logs -f novum`.

**CPUs and memory.**
- **Linux and WSL:** a container can use every core and all memory, within WSL's own limits in `.wslconfig`.
- **Docker Desktop on a Mac:** it uses all cores and half the memory by default. For large panels, raise the memory under Settings → Resources; the RNAhybrid arm has needed about 12 GB.
- Pass `--cores N` rather than `--cores all`.

**Clusters** without Docker usually have Apptainer, which runs the same image. Pull it once on a node with internet access, then run it from the project folder:

```bash
apptainer pull novum.sif docker://albertosola/novum_pipeline:v1.0.0
apptainer run --cleanenv --bind "$PWD":/work --pwd /work novum.sif --cores 8
```

`--cleanenv` keeps your shell's variables, such as `PYTHONPATH`, out of the pipeline's environments.

**Versions.** Image `vX.Y.Z` is built from git tag `vX.Y.Z`, and `latest` is the newest.

### Direct way: conda

You install **conda** and **Snakemake** yourself. On the first run, Snakemake installs everything else (RNAhybrid, IntaRNA, the Python and R post-processing) into one conda environment per step.

#### Conda

Skip this if you already have conda. [Miniforge](https://github.com/conda-forge/miniforge) is conda set up for conda-forge:

```bash
curl -LO "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
source ~/miniforge3/bin/activate
conda config --set channel_priority strict
```

Miniconda and Anaconda also work once you accept Anaconda's terms for their `defaults` channels. Otherwise Snakemake cannot build the environments, even though they do not use `defaults`:

```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

#### Snakemake

```bash
conda env update -f Workflow/Envs/environment.yaml --prune
conda activate snakemake-modern
```

The first command creates the environment, or resets an existing one to the file. To update it later, run it without `--prune`, which would also remove the config editor and `pytest`. It pins Snakemake 9.17.2, as the full `snakemake` package, which also installs pandas for the unit tests; `snakemake-minimal` does not.

#### Apple Silicon Macs

RNAhybrid and libfaketime have no `osx-arm64` conda builds, so on M-series Macs every pipeline environment must be built as `osx-64` (Intel) and run under Rosetta 2. Skip this on Linux and Intel Macs.

1. Install Rosetta 2, once per machine:

   ```bash
   sudo softwareupdate --install-rosetta --agree-to-license
   ```

2. Export this in every shell that runs Snakemake (or add it to `~/.zshrc`):

   ```bash
   export CONDA_SUBDIR=osx-64
   ```

   If you already ran the pipeline without it, delete `.snakemake/conda/` so the environments are rebuilt.

The first run builds the Intel environments and takes longer. `setup-macos.sh` tries to automate both steps, but on a real M-series Mac it failed: it runs the Rosetta installer without `sudo`, and the `subdir` setting it writes did not reach the environments Snakemake builds. Use the two steps above.

#### Config editor (optional)

A small desktop app can edit `Config/config.yaml` and start runs. To use it, add its GUI libraries to the same environment:

```bash
conda install -c conda-forge pyqt pyqtwebengine
pip install pywebview
```

### Results made before v1.0.0

This applies to both ways. Snakemake ties every result to the exact environment that made it. The environments are pinned since v1.0.0, so a results folder made before it is rerun from the start, calibration included. A folder made by the other install route is rerun too.

To keep such a folder as it is, add `--rerun-triggers mtime params input code` at the end of the command, after any output files. Placed before them, the flag takes their names as trigger values and fails.

## Quick start

With `snakemake-modern` active:

1. **Set your inputs.** In `Config/config.yaml`, list your miRNA FASTAs under `queries:` and the target CDS FASTAs under `targets:`. The keys pair them, so both blocks must use the same keys:

   ```yaml
   queries:
     my_sample: Data/Raw/RNAs/my_mirnas.fa
   targets:
     my_sample: Data/Raw/genomes/my_genome_cds.fna
   ```

   `data-prep/` can download these files and write both blocks for you; see [Preparing inputs automatically](#preparing-inputs-automatically). Every other setting has a working default.

2. **Run:**

   ```bash
   snakemake --use-conda --cores 8
   ```

   The first run is slower while Snakemake builds the conda environments.

3. **Read the results** in `<results_dir>/<sample>/`. `results_dir` defaults to `Data/Results`. See [Outputs](#outputs).

**With Docker**, follow the same steps from your project folder. In step 2, use the `docker run` line from [Simple way: Docker](#simple-way-docker) in place of `snakemake --use-conda`. The environments are already in the image, so the first run is not slower.

### Config editor

```bash
python UI/launcher.py
```

The editor edits `Config/config.yaml` and starts a run with one click. You can launch it from any directory. Snakemake's output goes to `Data/Results/.pipeline.log`; when a run fails or is cancelled, the app shows the last 80 lines.

Saving rewrites the whole file. Comments, and any key the editor does not know, are lost.

### Tests

The unit tests cover the pipeline's Python logic and the editor's config serializer. Add `pytest` to the environment, then run:

```bash
conda install pytest
pytest data-prep/tests Workflow/Scripts/tests UI/tests
```

In the bare environment 16 tests fail, because they call the IntaRNA, RNAcalibrate and libfaketime binaries. The UI tests are skipped unless `node` is on `PATH`. With the rule environments' binaries and `node` on `PATH`, all 487 pass.

With Docker, you can run the two Python suites against the image's own binaries. Build the test stage from a clone:

```bash
docker build --platform linux/amd64 --target test -t novum:test .
mkdir -p Data/Raw
docker run --rm -v "$PWD/Data/Raw":/work/Data/Raw:ro novum:test
```

All 413 pass. Two of them need the E. coli genome under `Data/Raw/genomes/`, and are skipped without it.

## Outputs

Each sample gets its own folder, `<results_dir>/<sample>/`. With every variant enabled:

```text
<sample>/
├── w_calibration/         RNAhybrid, calibrated p-values
│   ├── rnacalibrate.json
│   ├── rnahybrid_output.tsv
│   ├── tidy_output.csv
│   ├── rnahybrid_annotated.csv
│   └── rnahybrid_enhanced.txt
├── wo_calibration/        RNAhybrid, uncalibrated (same files, no JSON)
├── w_accessibility/       IntaRNA with accessibility
│   ├── intarna_output.csv
│   ├── intarna_tidy.csv
│   ├── intarna_annotated.csv
│   └── intarna_enhanced.txt
├── wo_accessibility/      IntaRNA without accessibility (same files)
├── consensus/
│   ├── consensus_annotated.csv
│   └── consensus_enhanced.txt
└── plots_<type>.pdf
```

A run creates only the variant folders it is configured for. `consensus/` is always created. `plots_<type>.pdf` appears only when `plots.type` is set.

Within each folder the files build on each other:

1. the raw tool output;
2. a tidy table;
3. `*_annotated.csv`, which adds gene metadata from the FASTA headers;
4. `*_enhanced.txt`, a readable report per hit with the duplex drawn out.

## Configuration

All settings live in [`Config/config.yaml`](Config/config.yaml):

| Block | Controls |
|---|---|
| `queries:`, `targets:` | Inputs, paired by key. The only block you must fill in. |
| `threads`, `max_suboptimal_hits`, `seed` | Settings both tools share. See [Shared settings](#shared-settings) and [Seed](#seed). |
| `rnacalibrate:` | Calibration of RNAhybrid's p-values. See [Calibration](#calibration-rnahybrid-arm). |
| `rnahybrid:` | RNAhybrid: energy cutoff, loop limits, p-value threshold, uncalibrated null model. |
| `intarna:` | IntaRNA: accessibility, prediction mode, model, energy table, and the `helix`, `seed`, `accessibility` and `output` sub-blocks. |
| `plots:` | Optional PDF plots. See [Plots](#plots). |
| `results_dir` | Where results go (default `Data/Results`). |

The older single `query: <path>` form still works: when `queries:` is absent, that one FASTA is run against every target.

### Preparing inputs automatically

`data-prep/prepare_inputs.py` builds the inputs from `(taxon, miRNA)` pairs:

1. It downloads each taxon's CDS and non-coding RNA from NCBI, and the miRNA sequences from miRBase.
2. It rewrites the `queries:` and `targets:` blocks of `Config/config.yaml` and leaves every other line untouched.

It does not run Snakemake.

```bash
python data-prep/prepare_inputs.py --pairs pairs.txt                  # one "<taxon> <miRNA>" per line
python data-prep/prepare_inputs.py --pairs pairs.txt --dry-run        # preview: no genome download, no config edit
python data-prep/prepare_inputs.py Veillonella_parvula:hsa-miR-200b-3p
python data-prep/prepare_inputs.py GCF_000153625.3:hsa-mir-515-5p     # pin one exact assembly
```

Pinning an assembly:

- An NCBI assembly accession in place of a taxon name pins that exact assembly.
- A mistyped accession is reported as `bad_accession`; it is never swapped for a name search.
- The accession becomes the sample key, so it also names the query FASTA and the results folder.

A taxon reaches the config only if its genome and at least one of its miRNAs resolved. Anything unresolved is listed in `data-prep/prep_report.tsv` instead of stopping the run. Full details are in [`data-prep/README.md`](data-prep/README.md).

## Options in detail

The defaults suit a first run. These sections are for tuning.

### Shared settings

`threads`, `max_suboptimal_hits` and `seed` sit at the top level, so both tools read the same value.

- **`threads`** sets the CPU count. IntaRNA runs one multi-threaded process per sample (`--threads`). RNAhybrid splits the targets into chunks and runs them through GNU Parallel.
- **`max_suboptimal_hits`** caps hits per pair:
  - A number `N` keeps at most `N` hits per pair on both arms.
  - `null` keeps every hit below the arm's energy cutoff. If RNAhybrid also has no cutoff, it then keeps only its best hit per pair.
  - IntaRNA reports at most 1000 hits per pair, so the arms can differ only for a pair with more than 1000 qualifying hits.

### Energy cutoffs

Each arm has its own hybridization-energy cutoff, because the two tools report energy on different scales:

```yaml
rnahybrid:
  max_hybrid_energy: -18      # RNAhybrid -e, on its mfe
intarna:
  max_hybrid_energy: -12.9    # on IntaRNA's E_hybrid
```

−18 kcal/mol is this project's own cutoff, not a value taken from a publication.

**How −12.9 is derived:**

- IntaRNA's `E_hybrid` includes a duplex-initiation penalty (+4.1 kcal/mol under Turner04) that RNAhybrid never adds.
- On identical duplexes, IntaRNA's number is higher by a roughly constant +5 kcal/mol.
- The pipeline uses +5.10, so −18 on RNAhybrid's scale is −12.9 on IntaRNA's.
- It warns when the two cutoffs are more than 1 kcal/mol away from that offset. Setting both to the same number triggers this warning.

**The match holds only without accessibility.** With accessibility on, IntaRNA often picks a different site in the gene, a more accessible one with weaker hybridization, so −12.9 is a stricter gate. On identical miRNA–gene pairs:

- 94.8% pass RNAhybrid's −18;
- 90.6% pass IntaRNA's −12.9 without accessibility;
- 68.0% pass IntaRNA's −12.9 with accessibility.

**Where the IntaRNA cutoff is applied:**

- Without accessibility, IntaRNA applies it itself (`--outMaxE`), and `tidy_intarna` applies it again afterwards.
- With accessibility, only `tidy_intarna` applies it. IntaRNA's own `--outMaxE` default of 0 still applies there, so a site whose accessibility penalty `ED1 + ED2` exceeds `|E_hybrid|` is never reported.

### Accessibility (IntaRNA arm)

`intarna.accessibility_variant` works like `calibration_variant`:

- `"on"` runs IntaRNA with accessibility (`--acc=C`) and writes to `w_accessibility/`.
- `"off"` runs it without accessibility (`--acc=N`) and writes to `wo_accessibility/`.
- `"both"` runs both.

Quote the values: YAML reads a bare `on`/`off` as true/false. The pipeline converts those back, but quoting avoids the ambiguity.

The two variants differ in four things: `--acc`, where the energy cutoff is applied, the search depth below, and the `Pu` floors. Every other `intarna` key applies to both. `helix.min_unpaired_probability` has no effect without accessibility, and the editor greys it out then.

**Accessibility columns.** IntaRNA reports accessibility as `ED1`/`ED2` (energy penalties) and `Pu1`/`Pu2` (unpaired probabilities). The `1` columns describe the gene and the `2` columns the miRNA.

**`Pu` floors.** Accessibility filters results only if you set one of two floors (both `null` by default):

- `intarna.min_target_unpaired_probability`: minimum `Pu1` (gene side).
- `intarna.min_query_unpaired_probability`: minimum `Pu2` (miRNA side).

How the floors behave:

- `tidy_intarna` applies them before choosing each pair's best hits. So a floor picks which site represents a pair, rather than only dropping pairs whose best site is inaccessible.
- Without accessibility (`--acc=N`) no `Pu` is computed. The floors are then not applied, and the pipeline warns.
- Published thresholds do not carry over directly, because `Pu` depends on the accessibility window. With this config's 80/40 target window, `Pu1` reads about 3.6× higher than at IntaRNA's default window of 150/100.
- `intarna.output.min_unpaired_probability` (`--outMinPu`) is a stricter, different bar: every interacting position on both sides must clear it. Setting it together with a floor applies both, and the pipeline warns.

**Search depth.** `intarna.accessibility_search_depth` sets the minimum number of hits IntaRNA reports per pair with accessibility on. The shipped config uses 150; the default if unset is 20.

- It exists because, with accessibility on, IntaRNA ranks hits by total energy while the cutoff applies to `E_hybrid`. The hit with the best `E_hybrid` can therefore rank below the top few.
- It never lowers the number of hits `max_suboptimal_hits` asks for.
- It has no effect when `max_suboptimal_hits` is `null`.

**Output columns.** `intarna.output.columns` is passed to IntaRNA as `--outCsvCols`. Keep at least:

- `id1,id2,start1,end1,start2,end2,E,E_hybrid` (needed by `tidy_intarna`);
- `subseqDP,hybridDP` (needed by `enhance_intarna`);
- `Pu1`/`Pu2` when you set a floor.

The seed columns are dropped automatically when no seed is set. The config editor always writes the default list.

### Calibration (RNAhybrid arm)

`rnacalibrate.calibration_variant` chooses which RNAhybrid variants a run produces:

- `on`: RNAcalibrate fits a null distribution `(xi, theta)`, and RNAhybrid computes p-values against it. Output goes to `w_calibration/`.
- `off`: no calibration. RNAhybrid uses `rnahybrid.distribution` if it is set. Otherwise it uses its built-in model for `rnahybrid.species` (default `3utr_human`). Output goes to `wo_calibration/`.
- `both`: both, in one run. This is the shipped setting.

**Length classes.** Calibration fits a separate null for each miRNA and each target-length class:

- `rnacalibrate.length_anchors` lists the class centres. The shipped list is `76,150,300,600,900,1500,3000`.
- Each gene joins the nearest class on a log scale.
- Classes up to 1500 nt are fitted directly. Longer classes take `theta` from the miRNA's fitted curve.
- `length_anchors: null` gives one fit per miRNA over the whole genome instead.

`rnacalibrate.json` records every fit. RNAhybrid then scores each target chunk with the parameters of its miRNA and length class. `rnahybrid.species` and `rnahybrid.distribution` do not affect the calibrated variant.

**`rnacalibrate.randomize_targets`** maps to RNAcalibrate's `-s`.

- With `true` (shipped), RNAcalibrate draws random sequences that match each target file's dinucleotide composition.
- With `false` (the default if unset), it uses the real target sequences as the random database. This is rarely what you want.

**`rnacalibrate.rng_seed`** makes the fits reproducible: two runs on identical inputs give identical fits.

- RNAcalibrate has no seed flag, so the pipeline fixes its clock with libfaketime. If it cannot confirm the clock is fixed, the run stops with an error.
- The shipped config uses `rng_seed: 1`. `null` seeds from the wall clock instead.
- The seed matters only with `randomize_targets: true`. Otherwise the pipeline warns and ignores it.

### Seed

The top-level `seed: "x,y"` (positions on the miRNA) forces a seed on both arms:

- RNAhybrid gets `-f x,y`.
- IntaRNA gets `--seedBP=y−x+1` (must be 2–20), `--seedQRange=x-y` and `--seedMaxUP=0`.

The shipped config has `seed: null`: IntaRNA gets `--noSeed`, and RNAhybrid gets no `-f`.

The `intarna.seed.*` keys apply only when a seed is set: `max_unpaired_bases` (`--seedMaxUP`, default 0), `target_range` (`--seedTRange`), `max_energy`, `max_hybrid_energy`, `min_unpaired_probability`, `forbid_gu`, `forbid_gu_at_ends` and `report_best_only`.

### Helix block (IntaRNA arm)

`intarna.model: B` selects IntaRNA's helix-block model, which builds an interaction from short helices joined by loops.

- The six `intarna.helix.*` keys tune that model. They have no effect under the other models (`X`, `S`, `P`). The pipeline warns if you set them anyway, and the editor greys them out.
- Model `B` changes predictions even with no helix keys set. On a let-7a test duplex it gave −27.6 kcal/mol, against −32.2 under `X`.
- Under model `B`, a seed must fit inside one helix: IntaRNA stops with an error when the seed has more base pairs than `helix.max_bp` (default 10). The pipeline checks this before any job starts.
- A seed `x,y` has `y − x + 1` base pairs, so `seed: 2,7` needs `helix.max_bp` of at least 6.

### Recipe: make IntaRNA behave like RNAhybrid

For a like-for-like comparison (`--noSeed --acc=N --mode=M --intLoopMax=30`):

```yaml
seed: null
intarna:
  accessibility_variant: "off"
  prediction_mode: M
  max_loop_size: 30
```

Without accessibility, IntaRNA's total energy equals `E_hybrid`, so its energy cutoff compares directly with RNAhybrid's.

### Plots

Plots are optional and drawn from the RNAhybrid tables only. With `plots.type: null`, no plot job runs and no R environment is built.

To re-render after changing the `plots` block, without re-running the rest:

```bash
snakemake --use-conda --cores 8 --forcerun build_plots
```

- **`type`**: `all`, or a comma-separated subset of the panels below. Panels are drawn in this order: `pvalue_distribution`, `position_pvalue`, `position_energy`, `volcano`, `pvalue_ecdf`, `calibration_delta`, `per_mirna`, `position_density`, `seed_class`. `calibration_delta` is skipped unless both calibration variants ran. The type is part of the file name (commas become hyphens), so different selections never overwrite each other.
- **`basesize`**: base font size.
- **`pvalue_threshold`**: keeps only points at or below it in `position_pvalue`, `position_energy` and `volcano`, and splits `position_density` into significant and non-significant hits.
- **`per_mirna_top_n`**: how many miRNAs the `per_mirna` panel shows, ranked by hit count (default 12).
- **`locus`**, **`gene`**, **`protein`**: known hits to highlight.
  - An entry is either `<target>` (any miRNA) or `<target>,<miRNA>` (that miRNA only).
  - Quote entries that contain spaces or punctuation, e.g. `"30S ribosomal protein S2,hsa-miR-X"`.
  - `pvalue_distribution` marks each one with a dashed line at its lowest p-value.

With both calibration variants, panels show the two side by side.

## How it works

Two prediction arms run from the same inputs, and the consensus combines them. The Snakefile defines the exact rule graph; the scripts live in `Workflow/Scripts/`.

```text
                    ┌ (RNAcalibrate) ── RNAhybrid ── tidy ── annotate ── enhance
                    │                                         │
queries + targets ──┤                                         ├── intersect ── enhance (consensus)
                    │                                         │
                    └ IntaRNA ────────────────────── tidy ── annotate ── enhance
```

### RNAhybrid arm

1. **rnacalibrate** (calibrated variant only): fits `(xi, theta)` per miRNA and length class.
2. **rnahybrid**: runs RNAhybrid on target chunks in parallel through GNU Parallel.
3. **tidy_rnahybrid**: turns the raw output into a table. `Position` is the first paired base divided by `Gene_length`, the same convention as the IntaRNA arm.
4. **annotate_rnahybrid**: adds gene metadata parsed from the target FASTA headers (`annotate.py`, shared with IntaRNA).
5. **enhance_rnahybrid**: writes the readable report with the alignment drawn out.
6. **build_plots** (optional): renders the PDF.

### IntaRNA arm

1. **intarna**: one multi-threaded IntaRNA process per sample and variant.
2. **tidy_intarna**: applies the energy cutoff and the `Pu` floors, then keeps each pair's best hits by `E_hybrid`.
3. **annotate_intarna**: adds gene metadata, as on the RNAhybrid side.
4. **enhance_intarna**: writes the readable report with energy, accessibility and the duplex.

### Consensus

1. **intersect**: joins the two arms on `(miRNA, Gene)`.
   - Each arm contributes its best hit per pair: lowest `Energy` for RNAhybrid (ties broken by `P_value`), lowest `E_hybrid` for IntaRNA (ties broken by `E`).
   - Positions do not constrain the join. A consensus row means both tools found a qualifying site in that gene, not necessarily the same site.
   - `Site_offset_nt` gives the distance between the two sites in nucleotides; it is positive when IntaRNA's site is downstream.
   - When an arm ran both variants, the consensus uses its `w_*` variant (calibrated, or with accessibility).
2. **enhance_consensus**: writes the readable report with both duplexes.

## Citation

If you use this workflow in a publication, please cite this repository and the external tools it relies on:

- **RNAhybrid**: Rehmsmeier, M., Steffen, P., Höchsmann, M., and Giegerich, R. (2004). Fast and effective prediction of microRNA/target duplexes. *RNA*, 10(10), 1507–1517. https://doi.org/10.1261/rna.5248604
- **IntaRNA**:
  - Mann, M., Wright, P. R., and Backofen, R. (2017). IntaRNA 2.0: enhanced and customizable prediction of RNA–RNA interactions. *Nucleic Acids Research*, 45(W1), W435–W439. https://doi.org/10.1093/nar/gkx279
  - Busch, A., Richter, A. S., and Backofen, R. (2008). IntaRNA: efficient prediction of bacterial sRNA targets incorporating target site accessibility and seed regions. *Bioinformatics*, 24(24), 2849–2856. https://doi.org/10.1093/bioinformatics/btn544
- **Snakemake**: Mölder, F., Jablonski, K. P., Letcher, B., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. https://doi.org/10.12688/f1000research.29032.2
- **GNU Parallel**: Tange, O. (2011). GNU Parallel: The command-line power tool. *;login: The USENIX Magazine*, 36(1), 42–47.
