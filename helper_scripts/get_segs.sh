#!/bin/bash

RESULTS_DIR="$1"
if [[ -z "$RESULTS_DIR" ]]; then
    echo "Usage: $0 <results_dir>" >&2
    exit 1
fi

RESULT_NAME=$(basename "$RESULTS_DIR")
OUT_DIR="$(dirname "$0")/$RESULT_NAME"
mkdir -p "$OUT_DIR"

OUT="$OUT_DIR/all_magic_segs.csv.gz"
TMP=$(mktemp)

# Remove output if it exists
rm -f "$OUT"

# Loop over directories in the results dir
for dir in "$RESULTS_DIR"/*/; do
    find "$dir" -maxdepth 2 -type f -name "*_seg.csv.gz" | while read file; do
        if [ ! -s "$TMP" ]; then
            echo "Processing $file (including header)"
            gzip -cd "$file" > "$TMP"
        else
            echo "Processing $file"
            gzip -cd "$file" | tail -n +2 >> "$TMP"
        fi
    done
done

gzip -c "$TMP" > "$OUT"
rm "$TMP"

echo "All segments concatenated into $OUT"
