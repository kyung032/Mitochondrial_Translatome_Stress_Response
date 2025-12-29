#!/usr/bin/env bash
set -euo pipefail

# === Config ===
BATCH="batch_15"
TRIM_DIR="$HOME/riboseq_atlas/mapping/fastq/trimmed/${BATCH}/ribo"
QC_DIR="$HOME/riboseq_atlas/mapping/qc/adapter_verify/${BATCH}"
mkdir -p "$QC_DIR"

echo "=== Batch 15 Trimming QC Report ==="
echo "Input Directory: $TRIM_DIR"
echo "Output Directory: $QC_DIR"
echo "---------------------------------------------------------------------------------------"
printf "%-15s %-15s %-15s %-15s %-10s\n" "Sample" "RPF(25-35nt)" "Noise(101nt)" "Signal:Noise" "Status"

# Find files
FILES=("$TRIM_DIR"/*.trimmed.fastq.gz)
if [ ! -e "${FILES[0]}" ]; then
    echo "❌ No trimmed files found in $TRIM_DIR"
    exit 1
fi

for f in "${FILES[@]}"; do
    SRR=$(basename "$f" | cut -d'.' -f1)
    DIST_FILE="$QC_DIR/${SRR}.length_dist.tsv"

    # Calc Distribution (Sample 2M reads for speed)
    zcat "$f" | head -n 8000000 | \
    awk 'NR%4==2{l[length($0)]++} END{for(k in l) print k, l[k]}' | sort -n > "$DIST_FILE"

    # Extract Stats
    STATS=$(awk '
    {
        if ($1 >= 25 && $1 <= 35) signal += $2;
        if ($1 == 101) noise += $2;
    }
    END {
        print (signal?signal:0), (noise?noise:0)
    }' "$DIST_FILE")

    SIG=$(echo "$STATS" | awk '{print $1}')
    NOISE=$(echo "$STATS" | awk '{print $2}')
    
    # Calc Ratio
    if [ "$NOISE" -eq 0 ]; then
        RATIO="Inf"
    else
        RATIO=$(awk -v s="$SIG" -v n="$NOISE" 'BEGIN {printf "%.2f", s/n}')
    fi

    # Judge
    if [ "$SIG" -gt "$NOISE" ] && [ "$SIG" -gt 1000 ]; then
        STATUS="✅ PASS"
    else
        STATUS="⚠️ CHECK"
    fi

    printf "%-15s %-15s %-15s %-15s %-10s\n" "$SRR" "$SIG" "$NOISE" "$RATIO" "$STATUS"
done

echo "---------------------------------------------------------------------------------------"
echo "QC Done. Detailed distribution files saved in $QC_DIR"
