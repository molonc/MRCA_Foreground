#!/bin/bash

RESULTS_DIR="$1"
if [[ -z "$RESULTS_DIR" ]]; then
    echo "Usage: $0 <results_dir>" >&2
    exit 1
fi

RESULT_NAME=$(basename "$RESULTS_DIR")
OUT_DIR="$(dirname "$0")/$RESULT_NAME"
mkdir -p "$OUT_DIR"
OUT_FILE="$OUT_DIR/libsizes_FGcells.txt"

> "$OUT_FILE"

# Loop through each subdirectory
for dir in "$RESULTS_DIR"/*/; do
    # Remove the trailing slash
    dir=${dir%/}

    # Find the matching file in the subdirectory
    file=$(find "$dir" -maxdepth 1 -name '*_breakpoint_file.csv.gz' | head -n 1)

    if [[ -n "$file" ]]; then
        # Count unique values in the first column, minus 1
        count=$(cat "$file" | cut -d',' -f1 | sort | uniq | wc -l)
        adjusted_count=$((count - 1)) # Remove the column header
        echo "$(basename "$dir"), $adjusted_count" >> "$OUT_FILE"
    else
        echo "$(basename "$dir"), 0,  No Brekapoit file found. Pipe Likely Failed" >> "$OUT_FILE"
    fi
done
