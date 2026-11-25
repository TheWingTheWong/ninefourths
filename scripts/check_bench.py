#!/usr/bin/env python3
"""CI check: recompute the target-2.807 rows with a small search range that
contains the optimum and compare with bench/cost_table.jsonl."""
import json, os, subprocess, sys

env = dict(os.environ, TAUS="strassen", CMAX="6", DMAX="250")
out = subprocess.run(["./target/release/cost", "--json"], env=env, capture_output=True, text=True, check=True).stdout
got = [json.loads(l) for l in out.splitlines() if l.strip()]
want = [json.loads(l) for l in open("bench/cost_table.jsonl") if abs(json.loads(l)["tau"] - 2.807355) < 1e-5]
ok = True
for g in got:
    w = next(x for x in want if x["kind"] == g["kind"])
    for key in ("c", "d", "k"):
        if g[key] != w[key]:
            print(f"{g['kind']}: {key} {g[key]} != {w[key]}"); ok = False
    for key in ("h", "log10_V", "log10_R"):
        if abs(g[key] - w[key]) > 1e-9 * abs(w[key]):
            print(f"{g['kind']}: {key} {g[key]} != {w[key]}"); ok = False
print("bench rows reproduced" if ok else "MISMATCH")
sys.exit(0 if ok else 1)
