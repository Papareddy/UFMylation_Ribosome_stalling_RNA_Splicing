#!/usr/bin/env bash
set -euo pipefail

# download_genomes.sh
#
# Robust download + decompression + indexing + validation of Ensembl "dna.primary_assembly" FASTA files.
# - Atomic downloads (writes *.part then mv)
# - Retries + resume for flaky networks
# - Optional disk-space preflight
# - Validation of required contigs (catches partial downloads)
# - MT vs M conflict handling (optional renamed FASTA)
# - Quick self-test mode: validate only, no downloads
#
# Usage:
#   ./download_genomes.sh
#   ./download_genomes.sh --test
#   ./download_genomes.sh --force
#   ./download_genomes.sh --make-m-alias
#
# Notes:
# - This script uses Ensembl release 115 URLs (HTTPS).
# - Run from repo root (Figure4_splicing_analysis).

# -------------------- config --------------------
ENSEMBL_RELEASE="115"
BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta"

HUMAN_URL="${BASE}/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
MOUSE_URL="${BASE}/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

HUMAN_DIR="data/human"
MOUSE_DIR="data/mouse"

HUMAN_GZ="${HUMAN_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
HUMAN_FA="${HUMAN_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

MOUSE_GZ="${MOUSE_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
MOUSE_FA="${MOUSE_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"

# Conservative free-space requirements (GB). Includes overhead for temp files.
NEED_GB_HUMAN=8
NEED_GB_MOUSE=6

# -------------------- args --------------------
MODE_TEST=0
MODE_FORCE=0
MAKE_M_ALIAS=0

for arg in "$@"; do
  case "$arg" in
    --test) MODE_TEST=1 ;;
    --force) MODE_FORCE=1 ;;
    --make-m-alias) MAKE_M_ALIAS=1 ;;
    -h|--help)
      cat <<EOF
Usage:
  ./download_genomes.sh [--test] [--force] [--make-m-alias]

  --test         Validate existing FASTA+FAI only (no download).
  --force        Re-download and rebuild everything even if files exist.
  --make-m-alias Create an additional human FASTA with MT header renamed to M.
EOF
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $arg" >&2
      exit 1
      ;;
  esac
done

# -------------------- helpers --------------------
msg() { echo "[INFO] $*"; }
warn() { echo "[WARN] $*" >&2; }
die() { echo -e "[ERROR] $*\n" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
}

avail_gb() {
  # Available GB on filesystem containing path
  local path="$1"
  local kb
  kb=$(df -Pk "$path" | awk 'NR==2{print $4}')
  echo $(( kb / 1024 / 1024 ))
}

need_space_gb() {
  local path="$1" need="$2"
  local have
  have="$(avail_gb "$path")"
  if (( have < need )); then
    die "Not enough free space at '$path'. Need ~${need}GB, have ~${have}GB.
Hint: run this on a work/scratch filesystem (e.g., GPFS), not /tmp. Or delete old data."
  fi
}

download_gz() {
  local url="$1" out="$2"
  local tmp="${out}.part"

  msg "Downloading: $(basename "$out")"
  # -f fail on HTTP errors, -L follow redirects
  # --retry* survive transient issues
  # -C - resume if tmp exists
  curl -fL --retry 10 --retry-delay 2 --retry-connrefused -C - \
    -o "$tmp" "$url"
  mv -f "$tmp" "$out"
}

inflate_fa() {
  local gz="$1" fa="$2"
  local tmp="${fa}.part"
  msg "Decompressing: $(basename "$gz") -> $(basename "$fa")"
  gunzip -c "$gz" > "$tmp"
  mv -f "$tmp" "$fa"
}

index_fa() {
  local fa="$1"
  msg "Indexing with samtools faidx: $(basename "$fa")"
  samtools faidx "$fa"
}

has_contig() {
  local fai="$1" contig="$2"
  awk -v c="$contig" '$1==c{found=1} END{exit found?0:1}' "$fai"
}

validate_mouse() {
  local fai="$1"
  for c in {1..19} X Y; do
    has_contig "$fai" "$c" || return 1
  done
  return 0
}

validate_human() {
  local fai="$1"
  for c in {1..22} X Y; do
    has_contig "$fai" "$c" || return 1
  done
  return 0
}

print_manual_fix() {
  local species="$1" url="$2" fa="$3" gz="$4"
  cat >&2 <<EOF

[MANUAL FIX] $species
1) Remove potentially broken files:
   rm -f "$gz" "$fa" "${fa}.fai"

2) Re-download (resume + retries):
   curl -fL --retry 10 --retry-delay 2 --retry-connrefused -C - -o "$gz" "$url"

3) Decompress to a real .fa:
   gunzip -c "$gz" > "$fa"

4) Index:
   samtools faidx "$fa"

5) Quick contig checks:
   cut -f1 "${fa}.fai" | head
   cut -f1 "${fa}.fai" | wc -l

EOF
}

quick_summary() {
  local fai="$1"
  msg "Index summary for $(basename "$fai"): contigs=$(wc -l < "$fai")"
  msg "First 10 contigs:"
  cut -f1 "$fai" | head -10 | sed 's/^/  - /'
}

make_m_alias_fasta() {
  local fa="$1"
  local out="${fa%.fa}.M.fa"

  if ! has_contig "${fa}.fai" "MT"; then
    warn "No MT contig found in $(basename "$fa"); not creating M alias."
    return 0
  fi
  if has_contig "${fa}.fai" "M"; then
    msg "Contig M already exists in $(basename "$fa"); not creating alias."
    return 0
  fi

  msg "Creating human FASTA with MT renamed to M: $(basename "$out")"
  local tmp="${out}.part"
  awk '
    BEGIN{OFS=""}
    /^>MT( |$)/ { sub(/^>MT/,">M"); print; next }
    { print }
  ' "$fa" > "$tmp"
  mv -f "$tmp" "$out"
  samtools faidx "$out"
  msg "Created: $out and ${out}.fai"
}

# -------------------- preflight --------------------
need_cmd curl
need_cmd gunzip
need_cmd samtools
need_space_gb "." 2  # minimal sanity check for current dir

mkdir -p "$HUMAN_DIR" "$MOUSE_DIR"

# -------------------- test-only mode --------------------
if (( MODE_TEST == 1 )); then
  msg "Running in --test mode (validate only, no downloads)."

  [[ -s "$MOUSE_FA" && -s "${MOUSE_FA}.fai" ]] || die "Mouse FASTA/index missing: $MOUSE_FA and/or ${MOUSE_FA}.fai"
  [[ -s "$HUMAN_FA" && -s "${HUMAN_FA}.fai" ]] || die "Human FASTA/index missing: $HUMAN_FA and/or ${HUMAN_FA}.fai"

  validate_mouse "${MOUSE_FA}.fai" || die "Mouse validation failed: missing one or more of 1..19, X, Y"
  validate_human "${HUMAN_FA}.fai" || die "Human validation failed: missing one or more of 1..22, X, Y"

  quick_summary "${MOUSE_FA}.fai"
  quick_summary "${HUMAN_FA}.fai"

  if has_contig "${HUMAN_FA}.fai" "MT" && ! has_contig "${HUMAN_FA}.fai" "M"; then
    warn "Human FASTA uses MT (mitochondrion). If your pipeline expects 'M', use --make-m-alias or update pipeline to use MT."
  fi

  msg "[DONE] --test passed."
  exit 0
fi

# -------------------- main flow --------------------
msg "=== MOUSE protects against partial downloads ==="
need_space_gb "$MOUSE_DIR" "$NEED_GB_MOUSE"

if (( MODE_FORCE == 1 )); then
  msg "--force set: removing existing mouse files."
  rm -f "$MOUSE_GZ" "$MOUSE_FA" "${MOUSE_FA}.fai"
fi

if [[ ! -s "$MOUSE_GZ" ]]; then
  download_gz "$MOUSE_URL" "$MOUSE_GZ" || {
    print_manual_fix "Mouse" "$MOUSE_URL" "$MOUSE_FA" "$MOUSE_GZ"
    die "Mouse download failed."
  }
else
  msg "Mouse .gz exists, skipping download."
fi

if [[ ! -s "$MOUSE_FA" ]]; then
  inflate_fa "$MOUSE_GZ" "$MOUSE_FA" || {
    print_manual_fix "Mouse" "$MOUSE_URL" "$MOUSE_FA" "$MOUSE_GZ"
    die "Mouse decompression failed."
  }
fi

if [[ ! -s "${MOUSE_FA}.fai" ]]; then
  index_fa "$MOUSE_FA" || {
    print_manual_fix "Mouse" "$MOUSE_URL" "$MOUSE_FA" "$MOUSE_GZ"
    die "Mouse indexing failed."
  }
fi

if ! validate_mouse "${MOUSE_FA}.fai"; then
  warn "Mouse FASTA index looks incomplete (missing one or more of 1..19, X, Y)."
  warn "Likely partial/corrupt download or decompression issue."
  rm -f "$MOUSE_GZ" "$MOUSE_FA" "${MOUSE_FA}.fai"
  print_manual_fix "Mouse" "$MOUSE_URL" "$MOUSE_FA" "$MOUSE_GZ"
  die "Mouse validation failed."
fi

msg "Mouse FASTA + index OK."
quick_summary "${MOUSE_FA}.fai"

msg "=== HUMAN protects against partial downloads ==="
need_space_gb "$HUMAN_DIR" "$NEED_GB_HUMAN"

if (( MODE_FORCE == 1 )); then
  msg "--force set: removing existing human files."
  rm -f "$HUMAN_GZ" "$HUMAN_FA" "${HUMAN_FA}.fai"
fi

if [[ ! -s "$HUMAN_GZ" ]]; then
  download_gz "$HUMAN_URL" "$HUMAN_GZ" || {
    print_manual_fix "Human" "$HUMAN_URL" "$HUMAN_FA" "$HUMAN_GZ"
    die "Human download failed."
  }
else
  msg "Human .gz exists, skipping download."
fi

if [[ ! -s "$HUMAN_FA" ]]; then
  inflate_fa "$HUMAN_GZ" "$HUMAN_FA" || {
    print_manual_fix "Human" "$HUMAN_URL" "$HUMAN_FA" "$HUMAN_GZ"
    die "Human decompression failed."
  }
fi

if [[ ! -s "${HUMAN_FA}.fai" ]]; then
  index_fa "$HUMAN_FA" || {
    print_manual_fix "Human" "$HUMAN_URL" "$HUMAN_FA" "$HUMAN_GZ"
    die "Human indexing failed."
  }
fi

if ! validate_human "${HUMAN_FA}.fai"; then
  warn "Human FASTA index looks incomplete (missing one or more of 1..22, X, Y)."
  warn "Likely partial/corrupt download or decompression issue."
  rm -f "$HUMAN_GZ" "$HUMAN_FA" "${HUMAN_FA}.fai"
  print_manual_fix "Human" "$HUMAN_URL" "$HUMAN_FA" "$HUMAN_GZ"
  die "Human validation failed."
fi

msg "Human FASTA + index OK."
quick_summary "${HUMAN_FA}.fai"

# MT vs M conflict info + optional alias FASTA
if has_contig "${HUMAN_FA}.fai" "MT" && ! has_contig "${HUMAN_FA}.fai" "M"; then
  warn "Human FASTA uses MT (mitochondrion). If your pipeline expects 'M', either:"
  warn "  (a) update pipeline to use 'MT', or"
  warn "  (b) run: ./download_genomes.sh --make-m-alias"
fi

if (( MAKE_M_ALIAS == 1 )); then
  make_m_alias_fasta "$HUMAN_FA"
fi

msg "[DONE] All required genomes are downloaded, validated, and indexed."
