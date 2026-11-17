#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$PROJECT_ROOT"
mkdir -p temp/3-distance log
export TMPDIR="$PROJECT_ROOT/temp/3-distance"
source script/load_config.sh
load_config conf/3-distance.yaml
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
"$tools_rscript" R/3-1-distance.R "${args[@]}" 2>&1 | tee log/3-distance.log
