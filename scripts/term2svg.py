#!/usr/bin/env python3
"""Render captured terminal text to an SVG (monospace, PASS green / FAIL red)."""
import html, sys

src, dst = sys.argv[1], sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else ""
lines = open(src).read().rstrip("\n").split("\n")
cw, lh, pad = 7.8, 18, 16
width = int(max(len(l) for l in lines) * cw + 2 * pad)
top = 34 if title else pad
height = int(len(lines) * lh + top + pad)
out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
       f'<rect width="100%" height="100%" rx="8" fill="#1a1a19"/>',
       '<style>text{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px;white-space:pre}</style>']
if title:
    out.append(f'<circle cx="18" cy="17" r="5" fill="#d95926"/><circle cx="34" cy="17" r="5" fill="#c98500"/><circle cx="50" cy="17" r="5" fill="#199e70"/>')
    out.append(f'<text x="68" y="21" fill="#c3c2b7">{html.escape(title)}</text>')
for i, l in enumerate(lines):
    y = top + (i + 1) * lh - 4
    if l.startswith("PASS") or l.startswith("FAIL"):
        col = "#3ecf8e" if l.startswith("PASS") else "#ff6b6b"
        out.append(f'<text x="{pad}" y="{y}"><tspan fill="{col}">{l[:4]}</tspan><tspan fill="#f0efec">{html.escape(l[4:])}</tspan></text>')
    else:
        col = "#3ecf8e" if "ALL PASS" in l else "#c3c2b7"
        out.append(f'<text x="{pad}" y="{y}" fill="{col}">{html.escape(l)}</text>')
out.append("</svg>")
open(dst, "w").write("\n".join(out))
