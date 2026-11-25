#!/usr/bin/env python3
"""Cost chart (SVG) from bench/cost_table.jsonl: how many decimal digits the
side of the identity has, per target exponent, on a log scale."""
import json, math, sys

rows = [json.loads(l) for l in open(sys.argv[1])]
dst = sys.argv[2]
series = [
    ("Exact identity, Pass's parameter rules", "#2a78d6", lambda r: r["kind"] == "faithful", "log10_V"),
    ("Exact identity, searched parameters", "#eb6834", lambda r: r["kind"] == "searched", "log10_V"),
    ("Before removing eps (mod-eps), searched", "#1baf7a", lambda r: r["kind"] == "searched", "log10_v"),
]
W, H = 860, 470
ml, mr, mt, mb = 78, 250, 58, 58
pw, ph = W - ml - mr, H - mt - mb
taus = sorted({r["tau"] for r in rows})
xmin, xmax = min(taus) - 0.02, max(taus) + 0.02
ys = [math.log10(r[k]) for _, _, f, k in series for r in rows if f(r)]
ymin, ymax = math.floor(min(ys)), math.ceil(max(ys))
X = lambda t: ml + (xmax - t) / (xmax - xmin) * pw  # higher exponent (easier) on the left
Y = lambda y: mt + (ymax - y) / (ymax - ymin) * ph
o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">',
     f'<rect width="100%" height="100%" fill="#fcfcfb"/>',
     f'<text x="{ml}" y="26" font-size="17" font-weight="600" fill="#1f1f1e">How big is the matrix-multiplication identity in Pass\'s ω ≤ 9/4 proof?</text>',
     f'<text x="{ml}" y="45" font-size="12.5" fill="#5f5e5a">Side of the smallest identity whose exponent is below the target, as a number of decimal digits (log scale)</text>']
for e in range(ymin, ymax + 1):
    y = Y(e)
    o.append(f'<line x1="{ml}" x2="{ml+pw}" y1="{y:.1f}" y2="{y:.1f}" stroke="#e6e5e0" stroke-width="1"/>')
    o.append(f'<text x="{ml-8}" y="{y+4:.1f}" font-size="12" fill="#5f5e5a" text-anchor="end">10^{e}</text>')
o.append(f'<text x="18" y="{mt+ph/2}" font-size="12" fill="#5f5e5a" transform="rotate(-90 18 {mt+ph/2})" text-anchor="middle">digits in the side</text>')
for t in taus:
    lab = "2.807 (Strassen)" if abs(t - math.log2(7)) < 1e-4 else f"{t:g}"
    o.append(f'<text x="{X(t):.1f}" y="{mt+ph+20}" font-size="12" fill="#5f5e5a" text-anchor="middle">{lab}</text>')
o.append(f'<text x="{ml+pw/2}" y="{H-12}" font-size="12" fill="#5f5e5a" text-anchor="middle">target exponent  (identity must use fewer than side^target multiplications)</text>')
o.append(f'<line x1="{ml}" x2="{ml+pw}" y1="{mt+ph}" y2="{mt+ph}" stroke="#8a8984" stroke-width="1"/>')
for name, col, f, k in series:
    pts = sorted(((r["tau"], math.log10(r[k]), r) for r in rows if f(r)), key=lambda p: -p[0])
    if not pts:
        continue
    d = " ".join(f"{'M' if i == 0 else 'L'}{X(t):.1f},{Y(y):.1f}" for i, (t, y, _) in enumerate(pts))
    o.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="2"/>')
    for t, y, r in pts:
        o.append(f'<circle cx="{X(t):.1f}" cy="{Y(y):.1f}" r="4.5" fill="{col}" stroke="#fcfcfb" stroke-width="2"><title>{name}, target {t:.4f}: side 10^({r[k]:.4g}); c={r["c"]} d={r["d"]} k={r["k"]} h={r["h"]}</title></circle>')
    t, y, r = pts[-1]
    o.append(f'<text x="{ml+pw+10}" y="{Y(y)+4:.1f}" font-size="12" fill="#1f1f1e"><tspan fill="{col}">●</tspan> {name}</text>')
o.append("</svg>")
open(dst, "w").write("\n".join(o))
