#!/bin/bash

RESULTS_DIR="$1"
if [[ -z "$RESULTS_DIR" ]]; then
    echo "Usage: $0 <results_dir>" >&2
    exit 1
fi

RESULT_NAME=$(basename "$RESULTS_DIR")
OUT_DIR="$(dirname "$0")/$RESULT_NAME"

FG="$OUT_DIR/libsizes_FGcells.txt"
HQ="$OUT_DIR/libsizes_HQ.txt"
TREE="$OUT_DIR/libsizes_treecells.txt"
OUT="$OUT_DIR/libsizes_merged.txt"

for f in "$FG" "$HQ" "$TREE"; do
    if [[ ! -f "$f" ]]; then
        echo "Missing required file: $f" >&2
        echo "Run get_all_summary_results.sh first." >&2
        exit 1
    fi
done

awk -F', ' -v fg_f="$FG" -v hq_f="$HQ" -v tree_f="$TREE" '
    FILENAME == fg_f   { fg[$1]=$2;   keys[$1]=1; next }
    FILENAME == hq_f   { hq[$1]=$2;   keys[$1]=1; next }
    FILENAME == tree_f { tree[$1]=$2; keys[$1]=1 }
    END {
        print "UID, CELLS_W_FG, MEDICC_INPUT, MEDICC_OUTPUT"
        for (k in keys)
            print k", "fg[k]+0", "hq[k]+0", "tree[k]+0
    }
' "$FG" "$HQ" "$TREE" | sort -t',' -k1,1 > "$OUT"

echo "Merged results saved to: $OUT"
