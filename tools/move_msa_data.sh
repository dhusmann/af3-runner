#!/usr/bin/env bash
set -euo pipefail          # safer bash: exit on errors, undefined vars, and pipeline failures

for topdir in */ ; do                       # iterate over every item in the current directory
  if [[ -d "${topdir}output" ]]; then       # step 1: does it have an output/ sub-dir?
    mkdir -p "${topdir}output_msa"          # step 2: create output_msa/ if it doesn’t exist
    # step 3: find and move *_data.json
    # for dry run:
    # find "${topdir}output" -maxdepth 2 -type f -name '*_data.json' -print
    # for wet run:
    find "${topdir}output" -maxdepth 2 -type f -name '*_data.json' -exec mv -t "${topdir}output_msa" {} +
  fi
done
