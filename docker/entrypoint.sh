#!/usr/bin/env bash
#----- Image entry point: `init` and `prepare-inputs` are ours, everything else goes to Snakemake -----#
set -euo pipefail
app=/opt/novum_pipeline

case "${1:-}" in
    init)
        if [ -e Config/config.yaml ]; then
            echo "init: Config/config.yaml already exists; left untouched." >&2
            exit 1
        fi
        mkdir -p Config
        # The template is the repo's own config, with placeholder inputs and results in Data/Results/
        sed 's|^results_dir:.*|results_dir: Data/Results/|' "$app/Config/config.yaml" > Config/config.yaml
        PYTHONPATH="$app/data-prep" python -c 'from config_writer import update_path_blocks as u
u("Config/config.yaml", {"my_sample": "Data/Raw/RNAs/my_mirnas.fa"}, {"my_sample": "Data/Raw/genomes/my_genome_cds.fna"})'
        echo "init: wrote Config/config.yaml (results go to Data/Results/). Replace my_sample under queries: and targets: with your inputs."
        ;;
    prepare-inputs)
        shift
        exec python "$app/data-prep/prepare_inputs.py" "$@"
        ;;
    *)
        exec snakemake --snakefile "$app/Snakefile" --use-conda --conda-prefix /opt/conda-envs "$@"
        ;;
esac
