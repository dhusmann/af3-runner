#!/usr/bin/env bash
#
# createAF3query.sh
#
# Unified AlphaFold3 job creation script supporting:
# - Multi-chain (protein/DNA/RNA) + ligand jobs
# - CCD ligands (--lig SAH,GTP or with counts SAH:2)
# - SMILES ligands via .sml files (--lig mycompound.sml or mycompound.sml:2)
# - PTMs:
#   - Single site:        --ptm <fasta_index> <position> <type>
#   - All lysines:        --ptm ALL <type>            (applies to last protein)
#   - All lysines (idx):  --ptm <fasta_index> ALL <type>
#   - Each lysine:        --ptm <fasta_index> EACH <type>  (creates one job per lysine)
#
# Defaults are cluster-specific but configurable via env vars or flags:
# - AF3_BASE_DIR  (default: /scratch/groups/ogozani/alphafold3)
# - BASE_INPUT_DIR (default: $AF3_BASE_DIR/jobs/inputs)
# - OUTPUT_DIR     (default: $AF3_BASE_DIR/jobs)
# - APPEND_CSV     (default: $AF3_BASE_DIR/folding_jobs.csv)
#
set -euo pipefail

AF3_BASE_DIR="${AF3_BASE_DIR:-/scratch/groups/ogozani/alphafold3}"
BASE_INPUT_DIR="${BASE_INPUT_DIR:-$AF3_BASE_DIR/jobs/inputs}"
OUTPUT_DIR="${OUTPUT_DIR:-$AF3_BASE_DIR/jobs}"
APPEND_CSV="${APPEND_CSV:-$AF3_BASE_DIR/folding_jobs.csv}"
APPEND_ENABLED=1
FORCE=0
DRY_RUN=0

MODEL_SEEDS_JSON='[1, 2, 8, 42, 88]'

usage() {
  cat <<'USAGE'
Usage:
  createAF3query.sh <p1.fa> [p2.fa ...] [--ptm ...] [--lig ...] [options]

PTMs:
  --ptm <idx> <pos> <type>        Apply PTM at position <pos> on FASTA #<idx>
  --ptm ALL <type>                Apply PTM to all lysines in the last protein
  --ptm <idx> ALL <type>          Apply PTM to all lysines in protein FASTA #<idx>
  --ptm <idx> EACH <type>         Create separate jobs for each lysine in protein #<idx>

Ligands:
  --lig <item[,item...]>          Items are CCD codes (e.g., SAH,GTP) or .sml files
                                 Counts supported via :N (e.g., SAH:2,mycompound.sml:2)

Options:
  --base-input-dir <path>         Where FASTA/.sml inputs live (default: $AF3_BASE_DIR/jobs/inputs)
  --output-dir <path>             Where job folders are created (default: $AF3_BASE_DIR/jobs)
  --append-csv <path>             CSV to append job names to (default: $AF3_BASE_DIR/folding_jobs.csv)
  --no-append                     Do not append job names to CSV
  --force                         Overwrite existing job directories
  -n, --dry-run                   Print what would be created, but do not write anything
  -h, --help                      Show this help

Examples:
  createAF3query.sh proteinA.fa proteinB.fa --ptm 2 43 me1 --lig SAH
  createAF3query.sh proteinA.fa --ptm ALL me1 --lig SAH
  createAF3query.sh proteinA.fa proteinB.fa --ptm 2 EACH me1 --lig SAH
  createAF3query.sh proteinA.fa proteinB.fa --lig SAH:2,GTP,mycompound.sml
USAGE
}

die() { echo "Error: $*" >&2; exit 1; }

json_escape() {
  local s="$1"
  if command -v python3 >/dev/null 2>&1; then
    python3 -c 'import json,sys; print(json.dumps(sys.argv[1]))' "$s"
  elif command -v python >/dev/null 2>&1; then
    python -c 'import json,sys; print(json.dumps(sys.argv[1]))' "$s"
  elif command -v jq >/dev/null 2>&1; then
    printf '%s' "$s" | jq -R .
  else
    die "Need python3/python/jq to escape SMILES strings"
  fi
}

ensure_csv_header() {
  local csv_path="$1"
  local header

  if [[ ! -f "$csv_path" || ! -s "$csv_path" ]]; then
    mkdir -p "$(dirname "$csv_path")"
    printf 'input_folder_name\n' > "$csv_path"
    return 0
  fi

  header="$(head -n 1 "$csv_path" | tr -d '\r')"
  if [[ "$header" != *input_folder_name* && "$header" != *folder* ]]; then
    # Prepend header (rare, but makes the file compatible with batch_reuse_msa.py / submit_dist.sh).
    local tmp
    tmp="$(mktemp)"
    {
      printf 'input_folder_name\n'
      cat "$csv_path"
    } > "$tmp"
    mv "$tmp" "$csv_path"
  fi
}

append_job_to_csv() {
  local job_name="$1"
  if [[ "$APPEND_ENABLED" -ne 1 ]]; then
    return 0
  fi
  ensure_csv_header "$APPEND_CSV"

  # Avoid duplicates (fast path: exact line match).
  if grep -Fxq "$job_name" "$APPEND_CSV" 2>/dev/null; then
    return 0
  fi
  printf '%s\n' "$job_name" >> "$APPEND_CSV"
}

# --- PTM CCD Code Mapping ---
declare -A CCD_MAP=(
  [me1]=MLZ
  [me2]=MLY
  [me3]=M3L
  [ac]=ALY
)

# --- Argument parsing ---
FASTA_FILES=()
PTM_ARGS=()
EACH_PTM_ARGS=()
LIGANDS_STR=""

while (($#)); do
  case "$1" in
    --base-input-dir)
      [[ $# -ge 2 ]] || die "--base-input-dir requires a path"
      BASE_INPUT_DIR="$2"
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || die "--output-dir requires a path"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --append-csv)
      [[ $# -ge 2 ]] || die "--append-csv requires a path"
      APPEND_CSV="$2"
      shift 2
      ;;
    --no-append)
      APPEND_ENABLED=0
      shift
      ;;
    --force)
      FORCE=1
      shift
      ;;
    -n|--dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --ptm)
      [[ $# -ge 3 ]] || die "--ptm requires arguments"
      if [[ "$2" == "ALL" ]]; then
        [[ $# -ge 3 ]] || die "--ptm ALL requires a PTM type"
        PTM_ARGS+=("ALL $3")
        shift 3
      elif [[ $# -ge 4 && "$3" == "ALL" ]]; then
        PTM_ARGS+=("$2 ALL $4")
        shift 4
      elif [[ $# -ge 4 && "$3" == "EACH" ]]; then
        EACH_PTM_ARGS+=("$2 $4")
        shift 4
      else
        [[ $# -ge 4 ]] || die "--ptm requires 3 arguments: <idx> <pos> <type>"
        PTM_ARGS+=("$2 $3 $4")
        shift 4
      fi
      ;;
    --lig)
      [[ $# -ge 2 ]] || die "--lig requires a comma-separated list of CCD codes or .sml files"
      LIGANDS_STR="$2"
      shift 2
      ;;
    -*)
      die "Unsupported flag: $1 (use --help)"
      ;;
    *)
      FASTA_FILES+=("$1")
      shift
      ;;
  esac
done

if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
  usage
  exit 1
fi

# --- Load sequences and detect molecule types ---
SEQUENCES=()
MOLECULE_TYPES=()
CLEAN_NAMES=()

for file in "${FASTA_FILES[@]}"; do
  local_path="$file"
  if [[ "$file" != */* ]]; then
    local_path="${BASE_INPUT_DIR}/${file}"
  fi
  [[ -f "$local_path" ]] || die "FASTA file not found: $local_path"

  seq="$(grep -v '^>' "$local_path" | tr -d '\n\r' | tr '[:lower:]' '[:upper:]')"
  [[ -n "$seq" ]] || die "Empty sequence in FASTA: $local_path"
  SEQUENCES+=("$seq")

  name="$(basename "${file%.fa}")"
  name="${name#h}"
  CLEAN_NAMES+=("$name")

  if [[ "$seq" =~ ^[GATC]+$ ]]; then
    MOLECULE_TYPES+=("dna")
  elif [[ "$seq" =~ ^[GAUC]+$ ]]; then
    MOLECULE_TYPES+=("rna")
  else
    MOLECULE_TYPES+=("protein")
  fi
done

# --- Build stoichiometric name part for sequences ---
declare -A MOLECULE_COUNTS=()
for name in "${CLEAN_NAMES[@]}"; do
  MOLECULE_COUNTS["$name"]=$(( ${MOLECULE_COUNTS[$name]:-0} + 1 ))
done

MOLECULE_NAME_PARTS=()
mapfile -t UNIQUE_CLEAN_NAMES < <(printf "%s\n" "${CLEAN_NAMES[@]}" | awk '!a[$0]++')
for name in "${UNIQUE_CLEAN_NAMES[@]}"; do
  count="${MOLECULE_COUNTS[$name]}"
  if ((count > 1)); then
    MOLECULE_NAME_PARTS+=("${count}x${name}")
  else
    MOLECULE_NAME_PARTS+=("$name")
  fi
done
MOLECULE_NAME_PART="$(IFS=-; echo "${MOLECULE_NAME_PARTS[*]}")"

# --- PTMs ---
declare -A PTM_MAP=()   # 1-indexed key -> "json,json,"
declare -A PTM_NAME_SUFFIX=()  # clean_name -> "_K43me1_KALLme1"

for ptm_arg in "${PTM_ARGS[@]+"${PTM_ARGS[@]}"}"; do
  if [[ "$ptm_arg" =~ ^ALL[[:space:]](.+)$ ]]; then
    ptm_type="${BASH_REMATCH[1]}"

    last_protein_idx=-1
    for i in "${!MOLECULE_TYPES[@]}"; do
      if [[ "${MOLECULE_TYPES[$i]}" == "protein" ]]; then
        last_protein_idx=$((i + 1))
      fi
    done
    (( last_protein_idx != -1 )) || die "No protein found for --ptm ALL"

    ccd_code="${CCD_MAP[$ptm_type]:-}"
    [[ -n "$ccd_code" ]] || die "Unknown PTM type '$ptm_type' (known: ${!CCD_MAP[*]})"

    sequence="${SEQUENCES[$((last_protein_idx - 1))]}"
    for ((pos=0; pos<${#sequence}; pos++)); do
      if [[ "${sequence:$pos:1}" == "K" ]]; then
        ptm_pos=$((pos + 1))
        PTM_MAP["$last_protein_idx"]+="{\"ptmType\": \"$ccd_code\", \"ptmPosition\": $ptm_pos},"
      fi
    done
    target_name="${CLEAN_NAMES[$((last_protein_idx - 1))]}"
    if (( ${MOLECULE_COUNTS[$target_name]:-0} > 1 )); then
      echo "Warning: PTM targets '$target_name' which appears ${MOLECULE_COUNTS[$target_name]} times; job naming may be ambiguous." >&2
    fi
    PTM_NAME_SUFFIX["$target_name"]+="_KALL${ptm_type}"
    continue
  fi

  if [[ "$ptm_arg" =~ ^([0-9]+)[[:space:]]ALL[[:space:]](.+)$ ]]; then
    ptm_file_idx="${BASH_REMATCH[1]}"
    ptm_type="${BASH_REMATCH[2]}"

    (( ptm_file_idx >= 1 && ptm_file_idx <= ${#SEQUENCES[@]} )) || die "Invalid FASTA index for --ptm: $ptm_file_idx"
    sequence_idx=$((ptm_file_idx - 1))
    [[ "${MOLECULE_TYPES[$sequence_idx]}" == "protein" ]] || die "PTMs can only be applied to proteins (idx $ptm_file_idx is ${MOLECULE_TYPES[$sequence_idx]})"

    ccd_code="${CCD_MAP[$ptm_type]:-}"
    [[ -n "$ccd_code" ]] || die "Unknown PTM type '$ptm_type' (known: ${!CCD_MAP[*]})"

    sequence="${SEQUENCES[$sequence_idx]}"
    lysine_count=0
    for ((pos=0; pos<${#sequence}; pos++)); do
      if [[ "${sequence:$pos:1}" == "K" ]]; then
        ptm_pos=$((pos + 1))
        PTM_MAP["$ptm_file_idx"]+="{\"ptmType\": \"$ccd_code\", \"ptmPosition\": $ptm_pos},"
        lysine_count=$((lysine_count + 1))
      fi
    done
    if (( lysine_count == 0 )); then
      echo "Warning: No lysine residues found in protein at index $ptm_file_idx" >&2
    fi
    target_name="${CLEAN_NAMES[$sequence_idx]}"
    if (( ${MOLECULE_COUNTS[$target_name]:-0} > 1 )); then
      echo "Warning: PTM targets '$target_name' which appears ${MOLECULE_COUNTS[$target_name]} times; job naming may be ambiguous." >&2
    fi
    PTM_NAME_SUFFIX["$target_name"]+="_KALL${ptm_type}"
    continue
  fi

  # Regular PTM case: "<idx> <pos> <type>"
  read -r ptm_file_idx ptm_pos ptm_type <<<"$ptm_arg"
  (( ptm_file_idx >= 1 && ptm_file_idx <= ${#SEQUENCES[@]} )) || die "Invalid FASTA index for --ptm: $ptm_file_idx"
  (( ptm_pos >= 1 )) || die "Invalid PTM position for --ptm: $ptm_pos"
  sequence_idx=$((ptm_file_idx - 1))
  [[ "${MOLECULE_TYPES[$sequence_idx]}" == "protein" ]] || die "PTMs can only be applied to proteins (idx $ptm_file_idx is ${MOLECULE_TYPES[$sequence_idx]})"

  sequence="${SEQUENCES[$sequence_idx]}"
  (( ptm_pos <= ${#sequence} )) || die "PTM position $ptm_pos out of range for FASTA index $ptm_file_idx (len=${#sequence})"

  residue="${sequence:$((ptm_pos - 1)):1}"
  target_name="${CLEAN_NAMES[$sequence_idx]}"
  if (( ${MOLECULE_COUNTS[$target_name]:-0} > 1 )); then
    echo "Warning: PTM targets '$target_name' which appears ${MOLECULE_COUNTS[$target_name]} times; job naming may be ambiguous." >&2
  fi
  PTM_NAME_SUFFIX["$target_name"]+="_${residue}${ptm_pos}${ptm_type}"

  ccd_code="${CCD_MAP[$ptm_type]:-}"
  [[ -n "$ccd_code" ]] || die "Unknown PTM type '$ptm_type' (known: ${!CCD_MAP[*]})"
  PTM_MAP["$ptm_file_idx"]+="{\"ptmType\": \"$ccd_code\", \"ptmPosition\": $ptm_pos},"
done

build_molecule_name_with_ptm_suffixes() {
  local extra_target_name="${1:-}"
  local extra_suffix="${2:-}"
  local parts=()

  for i in "${!UNIQUE_CLEAN_NAMES[@]}"; do
    local name="${UNIQUE_CLEAN_NAMES[$i]}"
    local base="${MOLECULE_NAME_PARTS[$i]}"
    local suffix="${PTM_NAME_SUFFIX[$name]:-}"

    if [[ -n "$extra_target_name" && "$name" == "$extra_target_name" ]]; then
      suffix+="$extra_suffix"
    fi

    parts+=("${base}${suffix}")
  done

  (IFS=-; echo "${parts[*]}")
}

FINAL_JOB_NAME="$(build_molecule_name_with_ptm_suffixes)"

# --- EACH PTMs (one job per lysine) ---
EACH_JOBS=()
HAS_EACH_PTM=false
for each_ptm_arg in "${EACH_PTM_ARGS[@]+"${EACH_PTM_ARGS[@]}"}"; do
  HAS_EACH_PTM=true
  read -r ptm_file_idx ptm_type <<<"$each_ptm_arg"
  (( ptm_file_idx >= 1 && ptm_file_idx <= ${#SEQUENCES[@]} )) || die "Invalid FASTA index for --ptm EACH: $ptm_file_idx"
  sequence_idx=$((ptm_file_idx - 1))
  [[ "${MOLECULE_TYPES[$sequence_idx]}" == "protein" ]] || die "PTMs can only be applied to proteins (idx $ptm_file_idx is ${MOLECULE_TYPES[$sequence_idx]})"

  target_name="${CLEAN_NAMES[$sequence_idx]}"
  if (( ${MOLECULE_COUNTS[$target_name]:-0} > 1 )); then
    echo "Warning: EACH PTM targets '$target_name' which appears ${MOLECULE_COUNTS[$target_name]} times; job naming may be ambiguous." >&2
  fi

  ccd_code="${CCD_MAP[$ptm_type]:-}"
  [[ -n "$ccd_code" ]] || die "Unknown PTM type '$ptm_type' (known: ${!CCD_MAP[*]})"

  sequence="${SEQUENCES[$sequence_idx]}"
  lysine_count=0
  for ((pos=0; pos<${#sequence}; pos++)); do
    if [[ "${sequence:$pos:1}" == "K" ]]; then
      ptm_pos=$((pos + 1))
      EACH_JOBS+=("$ptm_file_idx:$ptm_pos:$ptm_type:$ccd_code")
      lysine_count=$((lysine_count + 1))
    fi
  done
  if (( lysine_count == 0 )); then
    echo "Warning: No lysine residues found in protein at index $ptm_file_idx for EACH PTM" >&2
  else
    echo "Info: Found $lysine_count lysines in protein index $ptm_file_idx; creating $lysine_count jobs." >&2
  fi
done

# --- Ligands (CCD + SMILES via .sml files) ---
LIGAND_ENTRIES=()                 # per instance: "id|ccd|smiles"
declare -A LIGAND_SMILES_JSON=()  # id -> json string literal
LIGAND_NAME_PART_STR=""

if [[ -n "$LIGANDS_STR" ]]; then
  IFS=',' read -r -a LIGAND_ARRAY <<<"$LIGANDS_STR"
  for item in "${LIGAND_ARRAY[@]}"; do
    count=1
    code="$item"
    if [[ "$item" == *":"* ]]; then
      code="${item%:*}"
      count="${item#*:}"
      [[ "$count" =~ ^[1-9][0-9]*$ ]] || die "Invalid ligand count in '$item' (must be positive int)"
    fi

    if [[ "$code" == *.sml ]]; then
      smiles_path="$code"
      if [[ "$code" != */* ]]; then
        smiles_path="${BASE_INPUT_DIR}/${code}"
      fi
      [[ -f "$smiles_path" ]] || die "SMILES file not found: $smiles_path"
      smiles_string="$(cat "$smiles_path")"
      smiles_string="$(printf '%s' "$smiles_string" | tr -d '\r\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
      [[ -n "$smiles_string" ]] || die "Empty SMILES in file: $smiles_path"

      ligand_id="$(basename "${code%.sml}")"
      LIGAND_SMILES_JSON["$ligand_id"]="$(json_escape "$smiles_string")"
      for ((j=0; j<count; j++)); do
        LIGAND_ENTRIES+=("$ligand_id|smiles")
      done
    else
      ligand_id="$code"
      for ((j=0; j<count; j++)); do
        LIGAND_ENTRIES+=("$ligand_id|ccd")
      done
    fi
  done

	  # Build stoichiometric ligand name part (preserve order of first appearance)
	  declare -A LIGAND_COUNTS=()
	  declare -A LIGAND_SEEN=()
	  UNIQUE_LIGANDS=()
	  for entry in "${LIGAND_ENTRIES[@]+"${LIGAND_ENTRIES[@]}"}"; do
	    ligand_id="${entry%|*}"
	    LIGAND_COUNTS["$ligand_id"]=$(( ${LIGAND_COUNTS[$ligand_id]:-0} + 1 ))
	    if [[ -z "${LIGAND_SEEN[$ligand_id]:-}" ]]; then
	      UNIQUE_LIGANDS+=("$ligand_id")
	      LIGAND_SEEN["$ligand_id"]=1
	    fi
	  done

  LIGAND_NAME_PARTS=()
  for ligand_id in "${UNIQUE_LIGANDS[@]}"; do
    count="${LIGAND_COUNTS[$ligand_id]}"
    if ((count > 1)); then
      LIGAND_NAME_PARTS+=("${count}x${ligand_id}")
    else
      LIGAND_NAME_PARTS+=("$ligand_id")
    fi
  done

  LIGAND_NAME_PART_STR="$(IFS=-; echo "${LIGAND_NAME_PARTS[*]}")"
  FINAL_JOB_NAME+="-$LIGAND_NAME_PART_STR"
fi

create_job() {
  local job_name="$1"
  local each_ptm_file_idx="${2:-}"
  local each_ptm_pos="${3:-}"
  local each_ptm_type="${4:-}"
  local each_ccd_code="${5:-}"

  local job_dir="${OUTPUT_DIR}/${job_name}"

  if [[ -d "$job_dir" && -f "$job_dir/alphafold_input.json" && "$FORCE" -ne 1 ]]; then
    echo "Skipping existing job (use --force to overwrite): $job_dir" >&2
    return 2
  fi

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY RUN] Would create: $job_dir"
    return 0
  fi

  mkdir -p "$job_dir"

  local json_sequences_entries=""
  local chain_id_ascii=65 # 'A'

  for i in "${!SEQUENCES[@]}"; do
    local chain_id
    chain_id="$(printf "\\$(printf '%03o' "$chain_id_ascii")")"
    local sequence="${SEQUENCES[$i]}"
    local molecule_type="${MOLECULE_TYPES[$i]}"
    local molecule_idx=$((i + 1))

    [[ -n "$json_sequences_entries" ]] && json_sequences_entries+=","

    local modifications_json=""
    if [[ "$molecule_type" == "protein" ]]; then
      local mods=""
      if [[ -n "${PTM_MAP[$molecule_idx]:-}" ]]; then
        mods="${PTM_MAP[$molecule_idx]}"
      fi
      if [[ -n "$each_ptm_file_idx" && "$molecule_idx" -eq "$each_ptm_file_idx" ]]; then
        local each_ptm_json="{\"ptmType\": \"$each_ccd_code\", \"ptmPosition\": $each_ptm_pos}"
        mods+="$each_ptm_json,"
      fi
      if [[ -n "$mods" ]]; then
        mods="$(printf '%s' "$mods" | sed 's/,$//')"
        modifications_json=",\"modifications\": [ $mods ]"
      fi
    fi

    local molecule_entry
    molecule_entry="$(printf '\n    {\n      "%s": {\n        "id": "%s",\n        "sequence": "%s"%s\n      }\n    }' \
      "$molecule_type" "$chain_id" "$sequence" "$modifications_json")"
    json_sequences_entries+="$molecule_entry"
    ((chain_id_ascii++))
  done

  for entry in "${LIGAND_ENTRIES[@]+"${LIGAND_ENTRIES[@]}"}"; do
    local ligand_id="${entry%|*}"
    local ligand_type="${entry#*|}"
    local chain_id
    chain_id="$(printf "\\$(printf '%03o' "$chain_id_ascii")")"

    [[ -n "$json_sequences_entries" ]] && json_sequences_entries+=","

    local ligand_entry=""
    if [[ "$ligand_type" == "smiles" ]]; then
      local escaped_smiles="${LIGAND_SMILES_JSON[$ligand_id]}"
      [[ -n "$escaped_smiles" ]] || die "Internal error: missing SMILES for ligand '$ligand_id'"
      ligand_entry="$(printf '\n    {\n      "ligand": {\n        "id": "%s",\n        "smiles": %s\n      }\n    }' "$chain_id" "$escaped_smiles")"
    else
      ligand_entry="$(printf '\n    {\n      "ligand": {\n        "id": "%s",\n        "ccdCodes": ["%s"]\n      }\n    }' "$chain_id" "$ligand_id")"
    fi

    json_sequences_entries+="$ligand_entry"
    ((chain_id_ascii++))
  done

  cat > "$job_dir/alphafold_input.json" <<EOF
{
  "name": "$job_name",
  "modelSeeds": $MODEL_SEEDS_JSON,
  "sequences": [${json_sequences_entries}
  ],
  "dialect": "alphafold3",
  "version": 1
}
EOF

  append_job_to_csv "$job_name"
  echo "✅ Success! AlphaFold3 input created at: $job_dir"
}

if [[ "$HAS_EACH_PTM" == true ]]; then
  for each_job in "${EACH_JOBS[@]+"${EACH_JOBS[@]}"}"; do
    IFS=':' read -r ptm_file_idx ptm_pos ptm_type ccd_code <<<"$each_job"

    ligand_suffix=""
    if [[ -n "${LIGAND_NAME_PART_STR:-}" ]]; then
      ligand_suffix="-$LIGAND_NAME_PART_STR"
    fi

    target_name="${CLEAN_NAMES[$((ptm_file_idx - 1))]}"
    each_job_name="$(build_molecule_name_with_ptm_suffixes "$target_name" "_K${ptm_pos}${ptm_type}")${ligand_suffix}"
    create_job "$each_job_name" "$ptm_file_idx" "$ptm_pos" "$ptm_type" "$ccd_code"
  done
else
  create_job "$FINAL_JOB_NAME"
fi
