#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

BOOT_REPS="${BOOT_REPS:-9999}"
CONFIDENCE_LEVEL="${CONFIDENCE_LEVEL:-0.95}"
FACTOR_COUNT="${FACTOR_COUNT:-2}"
EMPIRICAL_SEED="${EMPIRICAL_SEED:-20261001}"
OUTPUT_DIR="${OUTPUT_DIR:-${script_dir}/results}"

export BOOT_REPS CONFIDENCE_LEVEL FACTOR_COUNT EMPIRICAL_SEED OUTPUT_DIR

Rscript empirical_ccemg_abc_hybrid.R
