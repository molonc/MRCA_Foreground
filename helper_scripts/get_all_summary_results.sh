#!/bin/bash

RESULTS_DIR="$1"
if [[ -z "$RESULTS_DIR" ]]; then
    echo "Usage: $0 <results_dir>" >&2
    exit 1
fi

SCRIPT_DIR="$(dirname "$0")"

echo "=== Running all summary scripts for: $RESULTS_DIR ==="

bash "$SCRIPT_DIR/get_cells_in_tree_nodes.sh" "$RESULTS_DIR"
bash "$SCRIPT_DIR/get_cells_in_tree.sh" "$RESULTS_DIR"
bash "$SCRIPT_DIR/get_libsize_FGcells.sh" "$RESULTS_DIR"
bash "$SCRIPT_DIR/get_libsize_HQCells.sh" "$RESULTS_DIR"
#bash "$SCRIPT_DIR/get_reads.sh" "$RESULTS_DIR"
#bash "$SCRIPT_DIR/get_segs.sh" "$RESULTS_DIR"
bash merge_libsize_summary.sh
echo "=== Done. Results saved to: $SCRIPT_DIR/$(basename "$RESULTS_DIR")/ ==="
