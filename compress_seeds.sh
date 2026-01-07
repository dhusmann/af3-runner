#!/usr/bin/env bash
# Usage:
#   ./compress_seeds.sh /scratch/groups/ogozani/alphafold3/output/human_test_set
#
# Tunables via env:
#   THREADS (default 8)            # pigz threads
#   COMPRESSION_LEVEL (default 7)  # 1-9 (7 is fast/near-max)
#   COMPRESSOR (auto|pigz|gzip|none; default auto)
#   PARALLEL_DIRS (default 1)      # process N runs at once (use >1 if gzip/none)

set -euo pipefail

ROOT="${1:-/scratch/groups/ogozani/alphafold3/output/human_test_set}"
THREADS="${THREADS:-4}"
COMPRESSION_LEVEL="${COMPRESSION_LEVEL:-7}"
COMPRESSOR="${COMPRESSOR:-auto}"
PARALLEL_DIRS="${PARALLEL_DIRS:-1}"

# Detect compressor if auto
if [[ "$COMPRESSOR" == "auto" ]]; then
  if command -v pigz >/dev/null 2>&1; then
    COMPRESSOR="pigz"
  elif command -v gzip >/dev/null 2>&1; then
    COMPRESSOR="gzip"
  else
    COMPRESSOR="none"
  fi
fi

echo "Root: $ROOT"
echo "Compressor: $COMPRESSOR  (level -$COMPRESSION_LEVEL)"
[[ "$COMPRESSOR" == "pigz" ]] && echo "pigz threads: $THREADS"
echo "Parallel directories: $PARALLEL_DIRS"

compress_one_dir() {
  local d="$1"
  echo ">>> Processing: $d"

  # Any seed-* dirs?
  shopt -s nullglob
  local seeds=("$d"/seed-*)
  shopt -u nullglob
  if [[ ${#seeds[@]} -eq 0 ]]; then
    echo "    No seed-* directories; skipping."
    return 0
  fi

  # Skip if archive exists already
  if [[ -f "$d/seeds.tar.gz" || -f "$d/seeds.tar" ]]; then
    echo "    Archive already exists; skipping."
    return 0
  fi

  # Create archive (keep originals)
  pushd "$d" >/dev/null
  case "$COMPRESSOR" in
    pigz)
      # stream tar -> pigz (avoid tar -I "pigz ...")
      tar -cf - seed-*/ | pigz -p "$THREADS" -"$COMPRESSION_LEVEL" > seeds.tar.gz
      tar -tzf seeds.tar.gz >/dev/null || { echo "    Verify failed; removing"; rm -f seeds.tar.gz; popd >/dev/null; return 1; }
      echo "    Created seeds.tar.gz ($(du -h seeds.tar.gz | awk '{print $1}')); originals kept."
      ;;
    gzip)
      tar -cf - seed-*/ | gzip -"$COMPRESSION_LEVEL" > seeds.tar.gz
      tar -tzf seeds.tar.gz >/dev/null || { echo "    Verify failed; removing"; rm -f seeds.tar.gz; popd >/dev/null; return 1; }
      echo "    Created seeds.tar.gz ($(du -h seeds.tar.gz | awk '{print $1}')); originals kept."
      ;;
    none)
      tar -cf seeds.tar seed-*/
      tar -tf seeds.tar >/dev/null || { echo "    Verify failed; removing"; rm -f seeds.tar; popd >/dev/null; return 1; }
      echo "    Created seeds.tar ($(du -h seeds.tar | awk '{print $1}')); originals kept."
      ;;
  esac
  popd >/dev/null
}

export -f compress_one_dir
export THREADS COMPRESSION_LEVEL COMPRESSOR

# Process immediate subdirectories of ROOT
find "$ROOT" -mindepth 1 -maxdepth 1 -type d -print0 \
| xargs -0 -I{} -P "$PARALLEL_DIRS" bash -c 'compress_one_dir "$@"' _ {}
echo "Done."
