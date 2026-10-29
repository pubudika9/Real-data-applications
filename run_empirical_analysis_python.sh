#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

python3 empirical_ccemg_abc_hybrid.py \
  --bootstrap-reps "${BOOT_REPS:-9999}" \
  --factor-count "${FACTOR_COUNT:-2}" \
  --seed "${EMPIRICAL_SEED:-20261001}" \
  --output-dir "${OUTPUT_DIR:-${script_dir}/results}"
