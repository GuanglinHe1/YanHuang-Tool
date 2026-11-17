#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$PROJECT_ROOT"
mkdir -p temp/1-preprocess log
export TMPDIR="$PROJECT_ROOT/temp/1-preprocess"
source script/load_config.sh
load_config conf/1-preprocess.yaml
mkdir -p "$paths_temp" log
export TMPDIR="$PROJECT_ROOT/$paths_temp"
args=()
while IFS= read -r key; do
    case "$key" in
        paths_*) args+=("--${key#paths_}" "${!key}") ;;
        runtime_*) args+=("--${key#runtime_}" "${!key}") ;;
        analysis_*) args+=("--${key#analysis_}" "${!key}") ;;
    esac
done < <(compgen -A variable | LC_ALL=C sort)
"$tools_rscript" R/1-1-preprocess.R "${args[@]}" 2>&1 | tee log/1-preprocess.log
