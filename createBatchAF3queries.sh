#!/usr/bin/env bash
# run_AF3_queries.sh  ─ iterate over a file of FASTA names and launch createAF3query.sh
# Usage: ./run_AF3_queries.sh [list_file]
# If no list_file is given, default to arabidopsis_kmt_list.txt in the current directory.

set -euo pipefail

LIST_FILE="${1:-arabidopsis_kmt_list.txt}"

if [[ ! -f $LIST_FILE ]]; then
  echo "Error: list file '$LIST_FILE' not found." >&2
  exit 1
fi

while IFS= read -r protein_fa || [[ -n $protein_fa ]]; do
  # skip empty or whitespace-only lines
  [[ -z $protein_fa ]] && continue

  echo "▶  Running createAF3query.sh on $protein_fa"
  createAF3query.sh "$protein_fa" atEEF1A1_noM.fa --ptm 2 43 me1 --lig SAH
  #createAF3query.sh "$protein_fa" atEEF1A1_noM.fa --ptm 2 186 me1 --lig SAH
  #createAF3query.sh "$protein_fa" atEEF1A1_noM.fa --ptm 2 395 me1 --lig SAH
  #createAF3query.sh "$protein_fa" atEEF1A1_noM.fa --ptm 2 226 me1 --lig SAH
  #createAF3query.sh "$protein_fa" atEEF1A1_noM.fa --ptm 2 35 me1 --lig SAH
  #createAF3query.sh "$protein_fa"  atEEF1A1_noM.fa --ptm 2 178 me1 --lig SAH
  #createAF3query.sh "$protein_fa"  atEEF1A1_noM.fa --ptm 2 379 me1 --lig SAH
  #createAF3query.sh "$protein_fa"  atEEF1A1_noM.fa --ptm 2 383 me1 --lig SAH
done < "$LIST_FILE"
