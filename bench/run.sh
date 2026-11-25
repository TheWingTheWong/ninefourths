#!/bin/sh
# Regenerates bench/cost_table.jsonl exactly as committed (about 30 minutes on 6 cores).
set -e
cd "$(dirname "$0")/.."
cargo build --release --bin cost
B=./target/release/cost
{
  MODE=faithful TAUS=strassen,2.75,2.7,2.6,2.5,2.45,2.4 $B --json
  MODE=search TAUS=strassen,2.75,2.7,2.6,2.5,2.45 CMAX=60 DMAX=2500 $B --json
  # at 2.4 the optimum sits beyond d = 2500, so the d range is wider
  MODE=search TAUS=2.4 CMAX=60 DMAX=8000 $B --json
} > bench/cost_table.jsonl
python3 scripts/table.py bench/cost_table.jsonl > bench/cost_table.md
python3 scripts/chart.py bench/cost_table.jsonl assets/cost.svg
