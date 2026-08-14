#!/bin/bash

RESULTS_DIR="$1"
if [[ -z "$RESULTS_DIR" ]]; then
    echo "Usage: $0 <results_dir>" >&2
    exit 1
fi

RESULT_NAME=$(basename "$RESULTS_DIR")
OUT_DIR="$(dirname "$0")/$RESULT_NAME"
mkdir -p "$OUT_DIR"
OUT_FILE="$OUT_DIR/libsizes_treecells_nodes.txt"

> "$OUT_FILE"

for dir in "$RESULTS_DIR"/*/; do
    dir=${dir%/}
    file=$(find "$dir" -maxdepth 2 -name '*output_segments_all.tsv_final_cn_profiles.tsv' | head -n 1)

    if [[ -n "$file" ]]; then
        count=$(grep -oE '(AT|SA)[A-Za-z0-9]+-A[A-Za-z0-9]+-R[0-9]+-C[0-9]+' "$file" | sort -u | wc -l)
        adjusted_count=$((count))
        echo $(basename "$dir"), $adjusted_count

        echo "$(basename "$dir"), $adjusted_count" >> "$OUT_FILE"
    else
        echo "$(basename "$dir"), 0, No Tree file found. Pipe likely failed" >> "$OUT_FILE"
    fi
done