#!/bin/bash
set -e

SCRIPT_DIR=$(dirname "$0")
PROJECT_DIR="$SCRIPT_DIR/.."

mkdir -p "$PROJECT_DIR/data/iupred"
mkdir -p "$PROJECT_DIR/data/tmp_fastas"

for pdb in "$PROJECT_DIR"/data/models/*.pdb; do
    fname=$(basename "$pdb" .pdb)
    fasta_file="$PROJECT_DIR/data/tmp_fastas/${fname}.fasta"
    out_file="$PROJECT_DIR/data/iupred/${fname}.txt"

    python3 "$SCRIPT_DIR/extract_seq_from_pdb.py" "$pdb" "$fasta_file"

    python3 /cs/labs/dina/tsori/af3_example/RedefineSubunit/iupred3/iupred3.py "$fasta_file" long > "$out_file"

    echo "Finished $fname"
done
