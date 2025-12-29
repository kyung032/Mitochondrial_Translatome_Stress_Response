#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# === Configuration (Adjust according to environment) ===
BATCH="batch_15"
THREADS=12

BASE="$HOME/riboseq_atlas/mapping"
# Input: Original Trimmed (101nt included)
TRIM_DIR="$BASE/fastq/trimmed/${BATCH}/ribo"
# Output: Mapping results
OUT="$BASE/aligned/ribo/${BATCH}"
WORK="/mnt/disks/ssd0/work/${BATCH}/ribo"
STARIDX="$BASE/ref/star_index"

# Automatically detect all Batch 15 samples
RIBO_RUNS=($(ls "$TRIM_DIR"/*.trimmed.fastq.gz | xargs -n 1 basename | cut -d'.' -f1 | sort -u))

need(){ command -v "$1" >/dev/null 2>&1 || { echo "❌ missing: $1"; exit 1; }; }
need STAR; need samtools

echo "=== ${BATCH} RIBO: Valley Cut (Max Length 50nt) & Map ==="
echo "Target Samples: ${RIBO_RUNS[*]}"

for s in "${RIBO_RUNS[@]}"; do
  INPUT_FQ="$TRIM_DIR/${s}.trimmed.fastq.gz"
  # Clean temporary file keeping only <= 50nt (Using SSD)
  CLEAN_FQ="/mnt/disks/ssd0/tmp_${BATCH}_${s}_le50.fastq.gz"

  echo ">>> Processing $s"

  # [Step 1] Valley Cut: Reads longer than 50nt are considered noise and removed
  echo "  [Filter] Keeping reads <= 50nt..."
  zcat "$INPUT_FQ" | \
  awk '{
      header=$0; getline seq; getline plus; getline qual;
      # Output only if length is <= 50
      if (length(seq) <= 50) {
          print header"\n"seq"\n"plus"\n"qual
      }
  }' | gzip > "$CLEAN_FQ"

  # [Step 2] STAR Mapping
  OD="$OUT/$s"; WD="$WORK/$s"
  mkdir -p "$OD" "$WD"
  
  BAM="$OD/$s.bam"
  FINAL="$OD/Log.final.out"
  
  rm -rf "${WD:?}/"* 2>/dev/null || true
  TMPDIR="$WD/_STARtmp.$(date +%s).$RANDOM"

  echo "  [STAR] Mapping clean reads..."
  STAR --runThreadN "$THREADS" \
    --genomeDir "$STARIDX" \
    --readFilesIn "$CLEAN_FQ" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$WD/" \
    --outTmpDir "$TMPDIR" \
    --outSAMtype BAM Unsorted \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 2 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --alignIntronMax 1 \
    > "$OD/$s.star.log" 2>&1

  # [Step 3] Sorting and Finalizing
  echo "  [Sort] $s"
  samtools sort -@ "$THREADS" -m 2G -o "$BAM" "$WD/Aligned.out.bam"
  samtools index -@ "$THREADS" "$BAM"
  
  [[ -s "$WD/Log.final.out" ]] && mv -f "$WD/Log.final.out" "$FINAL" || true
  
  # Remove temporary files (Free up space)
  rm -f "$CLEAN_FQ"
  rm -rf "${WD:?}/"* 2>/dev/null || true

  echo "--- Stats for $s ---"
  grep "Uniquely mapped reads %" "$FINAL" || true
  echo
done

echo "✅ DONE: All Batch 15 Samples Processed."
