#!/usr/bin/env python3
"""CI check: recompute the target-2.807 rows with a small search range that
contains the optimum and compare with bench/cost_table.jsonl."""
import json, os, subprocess, sys

def run(**kw):
    env = dict(os.environ, **kw)
    out = subprocess.run(["./target/release/cost", "--json"], env=env, capture_output=True, text=True, check=True).stdout
    return [json.loads(l) for l in out.splitlines() if l.strip()]

table = [json.loads(l) for l in open("bench/cost_table.jsonl")]
# 2.807: both rows on a grid containing the optimum; 2.7: the searched row on c = 6 only
got = run(TAUS="strassen", CMAX="6", DMAX="250") + run(TAUS="2.7", CMIN="6", CMAX="6", DMAX="400", MODE="search")
ok = True
for g in got:
    w = next(x for x in table if x["kind"] == g["kind"] and abs(x["tau"] - g["tau"]) < 1e-5)
    for key in ("c", "d", "k"):
        if g[key] != w[key]:
            print(f"{g['kind']}: {key} {g[key]} != {w[key]}"); ok = False
    for key in ("h", "log10_V", "log10_R"):
        if abs(g[key] - w[key]) > 1e-9 * abs(w[key]):
            print(f"{g['kind']}: {key} {g[key]} != {w[key]}"); ok = False
print("bench rows reproduced" if ok else "MISMATCH")
sys.exit(0 if ok else 1)
