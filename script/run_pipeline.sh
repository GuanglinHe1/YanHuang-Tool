#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
for step in "$root"/pipe/[1-6]-*.sh; do
    bash "$step"
done
