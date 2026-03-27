#!/usr/bin/env bash
# Extract_SEnt_and_BGexons.sh
# Purpose:
#   (A) From a pre-filtered coordinates TSV of skipped exons (SE), write:
#       - BED6 of skipped exons
#       - metadata TSV
#       - strand-aware nucleotide FASTA for skipped exons
#       - optional deduplicated nucleotide FASTA
#   (B) Optionally (one-time, cached) build a genome-wide exon background FASTA from a GTF:
#       - BED6 of exons
#       - strand-aware nucleotide FASTA for exons
#       - optional deduplicated nucleotide FASTA
#   The script intentionally STOPS after AA best-ORF generation + optional AA dedup/similarity filtering (no Pfam scan).
#
# Usage:
#   bash Extract_SEnt_and_BGexons.sh \
#        --lost /path/to/lost.tsv \
#        --preserved /path/to/preserved.tsv \
#        -g /path/to/GRCh38.genome.fa \
#        -o /path/to/output_dir \
#        [-u] \
#        [--event_types "SE,MXE" ] \
#        [--no_aa] [--aa_dedup] [--aa_filter_align] [--identity 0.9] \
#        [--bg_gtf /path/to/gencode.gtf] [--bg_cache_dir /path/to/cache]
#
# Notes:
# - Input TSV must have either the plain columns or .WT-suffixed versions:
#   ID or ID...1.WT, GeneID, geneSymbol, chr or chr.WT, strand or strand.WT, exonStart_0base or exonStart_0base.WT, exonEnd or exonEnd.WT
# - Genome exon background is derived from GTF “exon” features (all exons). It is cached and reused on reruns.
# - Requires: awk, sed, grep, coreutils, samtools; bedtools recommended.
# ------------------------ defaults ------------------------
LOST_TSV=""
PRESERVED_TSV=""
GENOME_FA=""
OUTDIR=""
DEDUP=0
BG_GTF=""
BG_CACHE_DIR=""
EVENT_TYPES="SE"   # default: process skipped exons only; optional comma-separated list, e.g. "SE,MXE"; empty means keep all
DIRECTION=0      # when 1, split each input by dPSI_num.WT sign into favored (>0) and unfavoured (<0)

# AA / filtering controls
MAKE_AA=1
AA_DEDUP=1   # default ON: remove identical AA sequences
AA_FILTER_ALIGN=0
IDENTITY=0.9

# ------------------------ args ----------------------------
usage(){
  echo "Usage: $0 --lost lost_coords.tsv --preserved preserved_coords.tsv -g genome.fa -o outdir [-u] [--event_types \"SE,MXE\"] [--direction] [--no_aa] [--aa_dedup] [--aa_filter_align] [--identity 0.9] [--bg_gtf gencode.gtf] [--bg_cache_dir cache_dir]" >&2
  echo
  echo "Notes:"
  echo "  - --aa_dedup is ON by default (identical AA sequences are removed unless you override with --aa_dedup 0)."
  echo "  - --event_types defaults to SE (skipped exons only); pass --event_types \"\" (empty) to keep all event types."
  exit 1
}

# ------------------------ args (short + long, any order) ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g)
      GENOME_FA="$2"; shift 2 ;;
    -o)
      OUTDIR="$2"; shift 2 ;;
    -u)
      DEDUP=1; shift 1 ;;

    --lost)
      LOST_TSV="$2"; shift 2 ;;
    --preserved)
      PRESERVED_TSV="$2"; shift 2 ;;

    --bg_gtf)
      BG_GTF="$2"; shift 2 ;;
    --bg_cache_dir)
      BG_CACHE_DIR="$2"; shift 2 ;;

    --event_types)
      EVENT_TYPES="$2"; shift 2 ;;
    --direction)
      DIRECTION=1; shift 1 ;;

    --no_aa)
      MAKE_AA=0; shift 1 ;;
    --aa_dedup)
      AA_DEDUP=1; shift 1 ;;
    --aa_filter_align)
      AA_FILTER_ALIGN=1; shift 1 ;;
    --identity)
      IDENTITY="$2"; shift 2 ;;

    --help|-h)
      usage ;;

    *)
      echo "[ERR] Unknown argument: $1" >&2
      usage ;;
  esac
done

if [[ -z "$LOST_TSV" || -z "$PRESERVED_TSV" || -z "$GENOME_FA" || -z "$OUTDIR" ]]; then
  echo "[ERR] Missing required arguments." >&2
  echo "      LOST_TSV=$LOST_TSV" >&2
  echo "      PRESERVED_TSV=$PRESERVED_TSV" >&2
  echo "      GENOME_FA=$GENOME_FA" >&2
  echo "      OUTDIR=$OUTDIR" >&2
  usage
fi

mkdir -p "$OUTDIR"
if [[ ! -f "$GENOME_FA" ]]; then echo "[ERR] Genome FASTA not found: $GENOME_FA" >&2; exit 2; fi

# Ensure FASTA is indexed for samtools fallback
if [[ ! -f "${GENOME_FA}.fai" ]]; then
  echo "[INFO] Indexing FASTA with samtools faidx..."
  samtools faidx "$GENOME_FA"
fi

# ------------------ normalize chromosome labels to FASTA style ------------------
normalize_bed_contigs_to_fasta(){
  local bed="$1"
  local fa="$2"

  local FA_HDR
  FA_HDR=$(grep -m1 '^>' "$fa" | sed 's/^>//; s/ .*//') || FA_HDR=""

  if [[ -z "$FA_HDR" ]]; then
    echo "[WARN] Could not read FASTA header; skipping contig normalization for $bed" >&2
    return 0
  fi

  if [[ "$FA_HDR" == chr* ]]; then
    # FASTA uses UCSC-style (chr1, chr2, chrM)
    local HAS_CHR_IN_BED
    HAS_CHR_IN_BED=$(awk 'BEGIN{FS="\t"} NR==1{print ($1 ~ /^chr/ ? 1 : 0); exit}' "$bed")
    if [[ "$HAS_CHR_IN_BED" -eq 0 ]]; then
      echo "[INFO] Normalizing BED contigs: adding 'chr' prefix to match FASTA (and MT->chrM): $bed"
      awk 'BEGIN{OFS="\t"} {c=$1; if(c!~ /^chr/){ if(c=="MT"||c=="M") c="chrM"; else c="chr" c;} $1=c; print }' \
        "$bed" > "$bed.tmp" && mv "$bed.tmp" "$bed"
    fi
  else
    # FASTA uses Ensembl-style (1,2,…,MT)
    local HAS_CHR_IN_BED
    HAS_CHR_IN_BED=$(awk 'BEGIN{FS="\t"} NR==1{print ($1 ~ /^chr/ ? 1 : 0); exit}' "$bed")
    if [[ "$HAS_CHR_IN_BED" -eq 1 ]]; then
      echo "[INFO] Normalizing BED contigs: removing 'chr' prefix to match FASTA (and chrM->MT): $bed"
      awk 'BEGIN{OFS="\t"} {c=$1; if(c=="chrM") c="MT"; else sub(/^chr/,"",c); $1=c; print }' \
        "$bed" > "$bed.tmp" && mv "$bed.tmp" "$bed"
    fi
  fi
}

# ----------- Process lost and preserved datasets ------------
process_coords_tsv(){
  local TSV="$1"
  local TAG="$2"   # lost | preserved
  local MODE="$3"  # all | favored | unfavoured

  if [[ ! -f "$TSV" ]]; then
    echo "[ERR] Coordinates TSV not found: $TSV" >&2
    exit 2
  fi

  local PREFIX="$TAG"
  if [[ -n "$MODE" && "$MODE" != "all" ]]; then
    PREFIX="${TAG}.${MODE}"
  fi
  local COORDS_BED="$OUTDIR/${PREFIX}.coords.bed"
  local META_TSV="$OUTDIR/${PREFIX}.meta.tsv"
  local FA_OUT_BEDTOOLS="$OUTDIR/${PREFIX}.nt.fa"
  local FA_OUT_FIDX="$OUTDIR/${PREFIX}.nt.fa"

  echo "[INFO] Building BED from coordinates TSV ($TAG): $TSV"
  if [[ -n "$EVENT_TYPES" ]]; then echo "[INFO] Filtering $TAG TSV by EventType: $EVENT_TYPES"; fi
  if [[ "$MODE" == "favored" ]]; then echo "[INFO] Direction split ON: keeping dPSI_num.WT > 0 ($TAG)"; fi
  if [[ "$MODE" == "unfavoured" ]]; then echo "[INFO] Direction split ON: keeping dPSI_num.WT < 0 ($TAG)"; fi

  # Accept either plain headers or .WT-suffixed headers, and filter by EventType if requested
  awk -F'\t' -v keep_types="$EVENT_TYPES" -v mode="$MODE" '
    BEGIN{
      n_keep=0;
      if(keep_types != ""){
        split(keep_types, tmp, ",");
        for(i in tmp){
          gsub(/^[ \t]+|[ \t]+$/, "", tmp[i]);
          if(tmp[i] != ""){ keep[tmp[i]]=1; n_keep++; }
        }
      }
    }
    NR==1{for(i=1;i<=NF;i++){h[$i]=i} next}
    {
      eti=(h["EventType"]?h["EventType"]:(h["EventType.WT"]?h["EventType.WT"]:0));
      ev=(eti?$(eti):"");
      if(n_keep>0 && !(ev in keep)) next;

      # optional direction split based on dPSI_num.WT
      if(mode != "" && mode != "all"){
        dpi=(h["dPSI_num.WT"]?h["dPSI_num.WT"]:(h["dPSI_num"]?h["dPSI_num"]:0));
        if(dpi==0){
          # if direction split requested but dPSI column missing, skip row
          next;
        }
        dpsi=$(dpi)+0.0;
        if(mode=="favored" && !(dpsi>0)) next;
        if(mode=="unfavoured" && !(dpsi<0)) next;
      }

      ci=(h["chr"]?h["chr"]:h["chr.WT"]);
      s0i=(h["exonStart_0base"]?h["exonStart_0base"]:h["exonStart_0base.WT"]);
      e1i=(h["exonEnd"]?h["exonEnd"]:h["exonEnd.WT"]);
      si=(h["strand"]?h["strand"]:h["strand.WT"]);
      idi=(h["ID"]?h["ID"]:h["ID...1.WT"]);
      geni=(h["geneSymbol"]?h["geneSymbol"]:h["geneSymbol.WT"]);

      chr=$(ci); s0=$(s0i); e1=$(e1i); strand=$(si); id=$(idi); gene=$(geni);
      if(chr=="" || s0=="" || e1=="" || strand=="" || id=="" || gene=="") next;
      gsub(/\"/,"",id); gsub(/\"/,"",gene);
      name=gene"|"id;
      print chr"\t"s0"\t"e1"\t"name"\t0\t"strand;
    }' "$TSV" > "$COORDS_BED"

  # Metadata TSV (minimal)
  (
    echo -e "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd";
    awk -F'\t' -v keep_types="$EVENT_TYPES" -v mode="$MODE" '
      BEGIN{
        n_keep=0;
        if(keep_types != ""){
          split(keep_types, tmp, ",");
          for(i in tmp){
            gsub(/^[ \t]+|[ \t]+$/, "", tmp[i]);
            if(tmp[i] != ""){ keep[tmp[i]]=1; n_keep++; }
          }
        }
      }
      NR==1{for(i=1;i<=NF;i++){h[$i]=i} next}
      {
        eti=(h["EventType"]?h["EventType"]:(h["EventType.WT"]?h["EventType.WT"]:0));
        ev=(eti?$(eti):"");
        if(n_keep>0 && !(ev in keep)) next;

        # optional direction split based on dPSI_num.WT
        if(mode != "" && mode != "all"){
          dpi=(h["dPSI_num.WT"]?h["dPSI_num.WT"]:(h["dPSI_num"]?h["dPSI_num"]:0));
          if(dpi==0) next;
          dpsi=$(dpi)+0.0;
          if(mode=="favored" && !(dpsi>0)) next;
          if(mode=="unfavoured" && !(dpsi<0)) next;
        }

        idi=(h["ID"]?h["ID"]:h["ID...1.WT"]);
        gi=(h["GeneID"]?h["GeneID"]:h["GeneID.WT"]);
        gsi=(h["geneSymbol"]?h["geneSymbol"]:h["geneSymbol.WT"]);
        chi=(h["chr"]?h["chr"]:h["chr.WT"]);
        si=(h["strand"]?h["strand"]:h["strand.WT"]);
        s0i=(h["exonStart_0base"]?h["exonStart_0base"]:h["exonStart_0base.WT"]);
        e1i=(h["exonEnd"]?h["exonEnd"]:h["exonEnd.WT"]);

        ID=$(idi); GeneID=$(gi); geneSymbol=$(gsi); chr=$(chi); strand=$(si); s0=$(s0i); e1=$(e1i);
        gsub(/\"/,"",ID); gsub(/\"/,"",GeneID); gsub(/\"/,"",geneSymbol);
        if(chr=="" || s0=="" || e1=="" || strand=="" || ID=="" || geneSymbol=="") next;
        print ID"\t"GeneID"\t"geneSymbol"\t"chr"\t"strand"\t"s0"\t"e1;
      }' "$TSV"
  ) > "$META_TSV"

  local countSig
  countSig=$(wc -l < "$COORDS_BED"); countSig=$((countSig))
  echo "[INFO] $TAG coordinate rows: $countSig (BED lines)"

  # Normalize contigs
  normalize_bed_contigs_to_fasta "$COORDS_BED" "$GENOME_FA"

  # Step 2: strand-aware NT FASTA
  if command -v bedtools >/dev/null 2>&1; then
    echo "[INFO] bedtools found; extracting strand-aware FASTA… ($TAG)"
    bedtools getfasta -fi "$GENOME_FA" -bed "$COORDS_BED" -s -name -fo "$FA_OUT_BEDTOOLS"
    echo "[OK] Wrote nucleotide FASTA: $FA_OUT_BEDTOOLS"
  else
    echo "[WARN] bedtools not found; using samtools faidx + Python for strand handling. ($TAG)"
    python - "$GENOME_FA" "$COORDS_BED" "$FA_OUT_FIDX" <<'PY'
import sys
import subprocess
comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
fa, bed, out = sys.argv[1:4]

def fetch(chrom, start1, end1):
    region = f"{chrom}:{start1}-{end1}"
    x = subprocess.check_output(["samtools","faidx", fa, region], text=True)
    return ''.join(l.strip() for l in x.splitlines() if not l.startswith(">"))

def rc(seq):
    return seq.translate(comp)[::-1]

with open(out, 'w') as fo:
    with open(bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, s0, e1, name, score, strand = line.rstrip().split('\t')
            s0 = int(s0); e1 = int(e1)
            start1 = s0 + 1
            end1 = e1
            seq = fetch(chrom, start1, end1)
            if strand == '-':
                seq = rc(seq)
            fo.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fo.write(seq[i:i+60] + "\n")
print(f"[OK] Wrote nucleotide FASTA: {out}")
PY
  fi

  # Step 2b: optional NT dedup
  local FA_NT_ORIG="$FA_OUT_BEDTOOLS"
  if [[ ! -s "$FA_NT_ORIG" ]]; then FA_NT_ORIG="$FA_OUT_FIDX"; fi
  local FA_NT="$FA_NT_ORIG"

  if [[ "$DEDUP" -eq 1 ]]; then
    local NT_DEDUP="$OUTDIR/${PREFIX}.nt.dedup.fa"
    echo "[INFO] Deduplicating nucleotide FASTA… ($TAG)"
    awk 'BEGIN{RS=">"; ORS=""} NR>1{n=split($0,a,"\n"); h=a[1]; s=""; for(i=2;i<=n;i++) s=s a[i]; if(!seen[s]++){print ">" h "\n" s "\n"}}' \
      "$FA_NT_ORIG" > "$NT_DEDUP"
    echo "[OK] Wrote deduplicated FASTA: $NT_DEDUP"
    FA_NT="$NT_DEDUP"
  fi

  # Step 3: best-ORF AA
  local AA_OUT="$OUTDIR/${PREFIX}.bestORF.aa.fa"
  local AA_OUT_DEDUP="$OUTDIR/${PREFIX}.bestORF.aa.dedup.fa"
  local AA_OUT_FILTER_ALIGN="$OUTDIR/${PREFIX}.bestORF.aa.filtered.align${IDENTITY}.fa"

  if [[ "$MAKE_AA" -eq 1 ]]; then
    echo "[INFO] Translating NT FASTA -> best-ORF AA FASTA (3-frame)… ($TAG)"
    python - "$FA_NT" "$AA_OUT" <<'PY'
import sys
nt_fa, aa_out = sys.argv[1:3]

codon_table = {
  'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
  'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
  'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
  'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
  'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
  'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
  'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
  'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
  'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
  'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
  'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
  'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
  'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
  'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
  'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
  'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def read_fasta(path):
  h=None; seq=[]
  with open(path) as f:
    for line in f:
      line=line.rstrip('\n')
      if not line: continue
      if line.startswith('>'):
        if h is not None:
          yield h, ''.join(seq)
        h=line[1:].strip(); seq=[]
      else:
        seq.append(line.strip())
    if h is not None:
      yield h, ''.join(seq)

def translate_frame(s, frame):
  s=s.upper(); pep=[]
  for i in range(frame, len(s)-2, 3):
    pep.append(codon_table.get(s[i:i+3], 'X'))
  return ''.join(pep)

def best_orf(nt):
  best=""; best_frame=0
  for fr in (0,1,2):
    pep = translate_frame(nt, fr)
    for seg in pep.split('*'):
      if len(seg) > len(best):
        best = seg
        best_frame = fr
  return best, best_frame

with open(aa_out, 'w') as out:
  for h, nt in read_fasta(nt_fa):
    pep, fr = best_orf(nt)
    out.write(f">{h}|best_frame={fr}|aa_len={len(pep)}\n")
    for i in range(0, len(pep), 60):
      out.write(pep[i:i+60] + "\n")

print(f"[OK] Wrote AA best-ORF FASTA: {aa_out}")
PY

    # Step 4: optional AA dedup
    local AA_TO_FILTER="$AA_OUT"
    if [[ "$AA_DEDUP" -eq 1 ]]; then
      echo "[INFO] Deduplicating AA FASTA… ($TAG)"
      awk 'BEGIN{RS=">"; ORS=""} NR>1{n=split($0,a,"\n"); h=a[1]; s=""; for(i=2;i<=n;i++) s=s a[i]; if(!seen[s]++){print ">" h "\n" s "\n"}}' \
        "$AA_OUT" > "$AA_OUT_DEDUP"
      echo "[OK] Wrote AA deduplicated FASTA: $AA_OUT_DEDUP"
      AA_TO_FILTER="$AA_OUT_DEDUP"
    fi

    # Step 4: optional AA similarity filtering
    if [[ "$AA_FILTER_ALIGN" -eq 1 ]]; then
      echo "[INFO] AA similarity filtering by alignment at identity >= $IDENTITY … ($TAG)"
      python - "$AA_TO_FILTER" "$AA_OUT_FILTER_ALIGN" "$IDENTITY" <<'PY'
import sys
fa_in, fa_out, ident = sys.argv[1], sys.argv[2], float(sys.argv[3])

try:
  from Bio import pairwise2
except Exception:
  sys.stderr.write("[ERR] Biopython is required for --aa_filter_align (pip install biopython)\n")
  raise

def read_fa(p):
  h=None; s=[]
  with open(p) as f:
    for line in f:
      line=line.rstrip('\n')
      if not line: continue
      if line.startswith('>'):
        if h is not None:
          yield h, ''.join(s)
        h=line[1:].strip(); s=[]
      else:
        s.append(line.strip())
    if h is not None:
      yield h, ''.join(s)

def identity(a, b):
  al = pairwise2.align.globalms(a, b, 2, -1, -4, -1, one_alignment_only=True)[0]
  aa, bb = al.seqA, al.seqB
  match=0; denom=0
  for x,y in zip(aa, bb):
    if x=='-' or y=='-':
      continue
    denom += 1
    if x==y:
      match += 1
  return (match/denom) if denom else 0.0

kept=[]
for h, seq in read_fa(fa_in):
  if not seq:
    continue
  drop=False
  for _, kseq in kept:
    if identity(seq, kseq) >= ident:
      drop=True
      break
  if not drop:
    kept.append((h, seq))

with open(fa_out, 'w') as out:
  for h, seq in kept:
    out.write(f">{h}\n")
    for i in range(0, len(seq), 60):
      out.write(seq[i:i+60] + "\n")

print(f"[OK] Wrote AA alignment-filtered FASTA: {fa_out} (kept {len(kept)})")
PY
    fi
  fi

  echo "[DONE] $TAG outputs:"
  echo "  BED      : $COORDS_BED"
  echo "  META     : $META_TSV"
  echo "  NT FASTA : ${FA_OUT_BEDTOOLS:-$FA_OUT_FIDX}"
  [[ "$DEDUP" -eq 1 ]] && echo "  NT dedup : $OUTDIR/${PREFIX}.nt.dedup.fa"
  [[ "$MAKE_AA" -eq 1 ]] && echo "  AA bestORF: $AA_OUT"
  [[ "$MAKE_AA" -eq 1 && "$AA_DEDUP" -eq 1 ]] && echo "  AA dedup : $AA_OUT_DEDUP"
  [[ "$MAKE_AA" -eq 1 && "$AA_FILTER_ALIGN" -eq 1 ]] && echo "  AA align : $AA_OUT_FILTER_ALIGN"
}

if [[ "$DIRECTION" -eq 1 ]]; then
  process_coords_tsv "$LOST_TSV" "lost" "favored"
  process_coords_tsv "$LOST_TSV" "lost" "unfavoured"
  process_coords_tsv "$PRESERVED_TSV" "preserved" "favored"
  process_coords_tsv "$PRESERVED_TSV" "preserved" "unfavoured"
else
  process_coords_tsv "$LOST_TSV" "lost" "all"
  process_coords_tsv "$PRESERVED_TSV" "preserved" "all"
fi

# ------------------ optional: genome-wide exon background (cached) ------------------
# If --bg_gtf is provided, build a BED6 of all exons from the GTF and extract a strand-aware NT FASTA.
# Results are cached and reused on reruns.
if [[ -n "$BG_GTF" ]]; then
  if [[ ! -f "$BG_GTF" ]]; then
    echo "[ERR] --bg_gtf file not found: $BG_GTF" >&2
    exit 3
  fi

  if [[ -z "$BG_CACHE_DIR" ]]; then
    BG_CACHE_DIR="$OUTDIR/cache_bg_exons"
  fi
  mkdir -p "$BG_CACHE_DIR"

  # Cache key based on filenames + sizes (fast, good enough for HPC workflows)
  FA_KEY=$(basename "$GENOME_FA")_$(stat -c%s "$GENOME_FA" 2>/dev/null || stat -f%z "$GENOME_FA")
  GTF_KEY=$(basename "$BG_GTF")_$(stat -c%s "$BG_GTF" 2>/dev/null || stat -f%z "$BG_GTF")
  CACHE_TAG="bgexons__${FA_KEY}__${GTF_KEY}"

  BG_BED="$BG_CACHE_DIR/${CACHE_TAG}.exons.bed"
  BG_NT_FA="$BG_CACHE_DIR/${CACHE_TAG}.exons.nt.fa"
  BG_NT_DEDUP_FA="$BG_CACHE_DIR/${CACHE_TAG}.exons.nt.dedup.fa"

  if [[ -s "$BG_NT_FA" ]]; then
    echo "[OK] Using cached background exon NT FASTA: $BG_NT_FA"
  else
    echo "[INFO] Building genome-wide exon background BED from GTF: $BG_GTF"

    # Extract all GTF exon features; convert to BED6 (0-based start, end-exclusive), keep strand.
    # Name uses gene_id|transcript_id|exon_number when available.
    if [[ "$BG_GTF" == *.gz ]]; then
      gunzip -c "$BG_GTF" | awk 'BEGIN{FS="\t"; OFS="\t"}
      $0 !~ /^#/ && $3=="exon" {
        chrom=$1; start0=$4-1; end=$5; strand=$7;
        attrs=$9;
        gene=""; tx=""; exonnum="";

        # gene_id
        if (match(attrs, /gene_id "[^"]+"/)) {
          s = substr(attrs, RSTART, RLENGTH);
          sub(/^gene_id "/, "", s);
          sub(/"$/, "", s);
          gene = s;
        }
        # transcript_id
        if (match(attrs, /transcript_id "[^"]+"/)) {
          s = substr(attrs, RSTART, RLENGTH);
          sub(/^transcript_id "/, "", s);
          sub(/"$/, "", s);
          tx = s;
        }
        # exon_number
        if (match(attrs, /exon_number "[^"]+"/)) {
          s = substr(attrs, RSTART, RLENGTH);
          sub(/^exon_number "/, "", s);
          sub(/"$/, "", s);
          exonnum = s;
        }

        name=gene;
        if(tx!="") name=name"|"tx;
        if(exonnum!="") name=name"|exon"exonnum;
        if(name=="") name="exon";
        print chrom, start0, end, name, 0, strand;
      }' > "$BG_BED"
    else
      awk 'BEGIN{FS="\t"; OFS="\t"}
      $0 !~ /^#/ && $3=="exon" {
        chrom=$1; start0=$4-1; end=$5; strand=$7;
        attrs=$9;
        gene=""; tx=""; exonnum="";

        # gene_id
        if (match(attrs, /gene_id "[^"]+"/)) {
          s = substr(attrs, RSTART, RLENGTH);
          sub(/^gene_id "/, "", s);
          sub(/"$/, "", s);
          gene = s;
        }
        # transcript_id
        if (match(attrs, /transcript_id "[^"]+"/)) {
          s = substr(attrs, RSTART, RLENGTH);
          sub(/^transcript_id "/, "", s);
          sub(/"$/, "", s);
          tx = s;
        }
        # exon_number
        if (match(attrs, /exon_number "[^"]+"/)) {
          s = substr(attrs, RSTART, RLENGTH);
          sub(/^exon_number "/, "", s);
          sub(/"$/, "", s);
          exonnum = s;
        }

        name=gene;
        if(tx!="") name=name"|"tx;
        if(exonnum!="") name=name"|exon"exonnum;
        if(name=="") name="exon";
        print chrom, start0, end, name, 0, strand;
      }' "$BG_GTF" > "$BG_BED"
    fi

    if [[ ! -s "$BG_BED" ]]; then
      echo "[ERR] Background exon BED is empty after parsing GTF. Check GTF format and awk compatibility: $BG_GTF" >&2
      exit 3
    fi

    # Normalize contigs to match genome FASTA
    normalize_bed_contigs_to_fasta "$BG_BED" "$GENOME_FA"

    echo "[INFO] Extracting strand-aware background exon NT FASTA…"
    if command -v bedtools >/dev/null 2>&1; then
      bedtools getfasta -fi "$GENOME_FA" -bed "$BG_BED" -s -name -fo "$BG_NT_FA"
    else
      echo "[WARN] bedtools not found; using samtools faidx + Python for background exon extraction." >&2
      python - "$GENOME_FA" "$BG_BED" "$BG_NT_FA" <<'PY'
import sys
import subprocess
comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
fa, bed, out = sys.argv[1:4]

def fetch(chrom, start1, end1):
    region = f"{chrom}:{start1}-{end1}"
    x = subprocess.check_output(["samtools","faidx", fa, region], text=True)
    return ''.join(l.strip() for l in x.splitlines() if not l.startswith(">"))

def rc(seq):
    return seq.translate(comp)[::-1]

with open(out, 'w') as fo:
    with open(bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, s0, e1, name, score, strand = line.rstrip().split('\t')
            s0 = int(s0); e1 = int(e1)
            start1 = s0 + 1
            end1 = e1
            seq = fetch(chrom, start1, end1)
            if strand == '-':
                seq = rc(seq)
            fo.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fo.write(seq[i:i+60] + "\n")
print(f"[OK] Wrote background exon NT FASTA: {out}")
PY
    fi

    echo "[OK] Cached background exon BED: $BG_BED"
    echo "[OK] Cached background exon NT FASTA: $BG_NT_FA"
  fi

  if [[ "$DEDUP" -eq 1 ]]; then
    if [[ -s "$BG_NT_DEDUP_FA" ]]; then
      echo "[OK] Using cached background exon dedup FASTA: $BG_NT_DEDUP_FA"
    else
      echo "[INFO] Deduplicating background exon NT FASTA…"
      awk 'BEGIN{RS=">"; ORS=""} NR>1{n=split($0,a,"\n"); h=a[1]; s=""; for(i=2;i<=n;i++) s=s a[i]; if(!seen[s]++){print ">" h "\n" s "\n"}}' \
        "$BG_NT_FA" > "$BG_NT_DEDUP_FA"
      echo "[OK] Wrote background exon deduplicated FASTA: $BG_NT_DEDUP_FA"
    fi
  fi

  # ------------------ optional: background exon AA best-ORF + filtering (cached) ------------------
  if [[ "$MAKE_AA" -eq 1 ]]; then
    BG_AA_FA="$BG_CACHE_DIR/${CACHE_TAG}.exons.bestORF.aa.fa"
    BG_AA_DEDUP_FA="$BG_CACHE_DIR/${CACHE_TAG}.exons.bestORF.aa.dedup.fa"
    BG_AA_ALIGN_FA="$BG_CACHE_DIR/${CACHE_TAG}.exons.bestORF.aa.filtered.align${IDENTITY}.fa"

    if [[ -s "$BG_AA_FA" ]]; then
      echo "[OK] Using cached background exon AA FASTA: $BG_AA_FA"
    else
      echo "[INFO] Translating background exon NT FASTA -> best-ORF AA FASTA…"
      python - "$BG_NT_FA" "$BG_AA_FA" <<'PY'
import sys

nt_fa, aa_out = sys.argv[1:3]

codon_table = {
  'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
  'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
  'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
  'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
  'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
  'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
  'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
  'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
  'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
  'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
  'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
  'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
  'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
  'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
  'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
  'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def read_fasta(path):
  h=None; seq=[]
  with open(path) as f:
    for line in f:
      line=line.rstrip('\n')
      if not line: continue
      if line.startswith('>'):
        if h is not None:
          yield h, ''.join(seq)
        h=line[1:].strip(); seq=[]
      else:
        seq.append(line.strip())
    if h is not None:
      yield h, ''.join(seq)

def translate_frame(s, frame):
  s=s.upper(); pep=[]
  for i in range(frame, len(s)-2, 3):
    cod=s[i:i+3]
    pep.append(codon_table.get(cod, 'X'))
  return ''.join(pep)

def best_orf(nt):
  best=""; best_frame=0
  for fr in (0,1,2):
    pep = translate_frame(nt, fr)
    for seg in pep.split('*'):
      if len(seg) > len(best):
        best = seg
        best_frame = fr
  return best, best_frame

with open(aa_out, 'w') as out:
  for h, nt in read_fasta(nt_fa):
    pep, fr = best_orf(nt)
    out.write(f">{h}|best_frame={fr}|aa_len={len(pep)}\n")
    for i in range(0, len(pep), 60):
      out.write(pep[i:i+60] + "\n")

print(f"[OK] Wrote background exon AA best-ORF FASTA: {aa_out}")
PY
    fi

    BG_AA_TO_FILTER="$BG_AA_FA"
    if [[ "$AA_DEDUP" -eq 1 ]]; then
      if [[ -s "$BG_AA_DEDUP_FA" ]]; then
        echo "[OK] Using cached background exon AA dedup FASTA: $BG_AA_DEDUP_FA"
      else
        echo "[INFO] Deduplicating background exon AA FASTA…"
        awk 'BEGIN{RS=">"; ORS=""} NR>1{n=split($0,a,"\n"); h=a[1]; s=""; for(i=2;i<=n;i++) s=s a[i]; if(!seen[s]++){print ">" h "\n" s "\n"}}' \
          "$BG_AA_FA" > "$BG_AA_DEDUP_FA"
        echo "[OK] Wrote background exon AA dedup FASTA: $BG_AA_DEDUP_FA"
      fi
      BG_AA_TO_FILTER="$BG_AA_DEDUP_FA"
    fi

    if [[ "$AA_FILTER_ALIGN" -eq 1 ]]; then
      if [[ -s "$BG_AA_ALIGN_FA" ]]; then
        echo "[OK] Using cached background exon AA alignment-filtered FASTA: $BG_AA_ALIGN_FA"
      else
        echo "[INFO] AA similarity filtering (background exons) at identity >= $IDENTITY …"
        python - "$BG_AA_TO_FILTER" "$BG_AA_ALIGN_FA" "$IDENTITY" <<'PY'
import sys

fa_in, fa_out, ident = sys.argv[1], sys.argv[2], float(sys.argv[3])

try:
  from Bio import pairwise2
except Exception:
  sys.stderr.write("[ERR] Biopython is required for --aa_filter_align (pip install biopython)\n")
  raise

def read_fa(p):
  h=None; s=[]
  with open(p) as f:
    for line in f:
      line=line.rstrip('\n')
      if not line: continue
      if line.startswith('>'):
        if h is not None:
          yield h, ''.join(s)
        h=line[1:].strip(); s=[]
      else:
        s.append(line.strip())
    if h is not None:
      yield h, ''.join(s)

def identity(a, b):
  al = pairwise2.align.globalms(a, b, 2, -1, -4, -1, one_alignment_only=True)[0]
  aa, bb = al.seqA, al.seqB
  match=0; denom=0
  for x,y in zip(aa, bb):
    if x=='-' or y=='-':
      continue
    denom += 1
    if x==y:
      match += 1
  return (match/denom) if denom else 0.0

kept=[]
for h, seq in read_fa(fa_in):
  if not seq:
    continue
  drop=False
  for _, kseq in kept:
    if identity(seq, kseq) >= ident:
      drop=True
      break
  if not drop:
    kept.append((h, seq))

with open(fa_out, 'w') as out:
  for h, seq in kept:
    out.write(f">{h}\n")
    for i in range(0, len(seq), 60):
      out.write(seq[i:i+60] + "\n")

print(f"[OK] Wrote background exon AA alignment-filtered FASTA: {fa_out} (kept {len(kept)})")
PY
      fi
    fi

    echo "[DONE] Background exon AA cache ready."
    echo "  BG AA  : $BG_AA_FA"
    [[ "$AA_DEDUP" -eq 1 ]] && echo "  BG AA dedup: $BG_AA_DEDUP_FA"
    [[ "$AA_FILTER_ALIGN" -eq 1 ]] && echo "  BG AA align: $BG_AA_ALIGN_FA"
  fi

  echo "[DONE] Background exon cache ready."
  echo "  BG BED : $BG_BED"
  echo "  BG NT  : $BG_NT_FA"
  [[ "$DEDUP" -eq 1 ]] && echo "  BG NT dedup: $BG_NT_DEDUP_FA"
fi

echo "[DONE] Finished coordinate extraction for lost + preserved."
[[ -n "$BG_GTF" ]] && echo "[DONE] Background exon cache dir: $BG_CACHE_DIR"

