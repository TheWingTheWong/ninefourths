#!/usr/bin/env python3
"""Markdown table from bench/cost_table.jsonl."""
import json, math, sys

rows = [json.loads(l) for l in open(sys.argv[1])]

def p10(x):
    # 10^x with x itself possibly huge -> "10^(7.33×10^15)" or "10^1,078,625"
    if x < 1e7:
        return f"10^{round(x):,}"
    e = int(math.floor(math.log10(x)))
    return f"10^({x / 10**e:.2f}×10^{e})"

def hfmt(h):
    return f"{int(h):,}" if h < 1e6 else f"{h / 10**int(math.log10(h)):.2f}×10^{int(math.log10(h))}"

print("| target exponent | parameters | c | d | k | h | mod-ε side v (exponent) | ε-degree | exact side V | multiplications R | exact exponent | Strassen at side V |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for tau in sorted({r["tau"] for r in rows}, reverse=True):
    for kind in ("faithful", "searched"):
        for r in rows:
            if r["tau"] == tau and r["kind"] == kind:
                t = "log2 7 = 2.807" if abs(tau - math.log2(7)) < 1e-4 else f"{tau:g}"
                strassen = math.ceil(r["log10_V"] / math.log10(2)) * math.log10(7)
                print(f"| {t} | {'Pass' if kind == 'faithful' else 'searched'} | {r['c']} | {r['d']:,} | {r['k']:,} | {hfmt(r['h'])} | {p10(r['log10_v'])} ({r['exp_mod_eps']:.3f}) | {p10(r['log10_deg'])} | **{p10(r['log10_V'])}** | {p10(r['log10_R'])} | {r['exp_exact']:.6f} | {p10(strassen)} |")
