#!/usr/bin/env bash
#----- Build time only: Snakemake creates every rule env at the prefix the image runs with. -----#
#----- A one-record fixture with plots on makes the DAG reach all four env files -----#
set -euo pipefail
app=/opt/novum_pipeline
tmp=$(mktemp -d)

mkdir -p "$tmp/Config"
cp "$app/Config/config.yaml" "$tmp/Config/config.yaml"
cp "$app/docker/build-fixture/query.fa" "$app/docker/build-fixture/target.fna" "$tmp/"
PYTHONPATH="$app/data-prep" python - "$tmp/Config/config.yaml" <<'EOF'
import sys
from config_writer import update_path_blocks
update_path_blocks(sys.argv[1], {"fixture": "query.fa"}, {"fixture": "target.fna"})
EOF

"$app/docker/entrypoint.sh" -d "$tmp" --conda-create-envs-only --cores 1 --config 'plots={type: all}'
rm -rf "$tmp"
