#!/usr/bin/env bash
#
# Compatibility wrapper for historical createAF3query_withSMILES.sh behavior.
# Uses the unified implementation in createAF3query.sh with SMILES defaults.
#
set -euo pipefail

AF3_BASE_DIR="${AF3_BASE_DIR:-/scratch/groups/ogozani/alphafold3}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec "${SCRIPT_DIR}/createAF3query.sh" \
  --output-dir "${AF3_BASE_DIR}/jobs/human_test_set" \
  --append-csv "${AF3_BASE_DIR}/folding_jobs_nsd2i.csv" \
  "$@"

