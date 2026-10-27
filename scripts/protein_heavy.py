#!/usr/bin/env python3
"""Receptor (PDB ou PDBQT) -> PDB so com atomos pesados de residuos (ATOM), TER entre cadeias.
O pdb2gmx -ignh recoloca os H do campo de forca.

Uso: protein_heavy.py RECEPTOR OUT_PDB
"""
import sys

AD4 = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "O": "O",
       "SA": "S", "S": "S", "HD": "H", "HS": "H", "H": "H"}

src, out = sys.argv[1], sys.argv[2]
pdbqt = src.lower().endswith(".pdbqt")
lines, last, n = [], None, 0
for l in open(src, errors="replace"):
    l = l.rstrip("\r\n")
    if not l.startswith("ATOM"):
        continue
    name = l[12:16].strip()
    el = AD4.get(l[77:79].strip(), name[0]) if pdbqt else (l[76:78].strip() or name.lstrip("0123456789")[0])
    if el.upper() == "H":
        continue
    ch = l[21]
    if last is not None and ch != last:
        lines.append("TER")
    last = ch
    n += 1
    nm = name if len(name) == 4 else " " + name
    lines.append(f"ATOM  {n % 100000:>5} {nm:<4} {l[17:20]:>3} {ch}{l[22:26]}    {l[30:54]}  1.00  0.00          {el:>2}")
lines += ["TER", "END"]
open(out, "w").write("\n".join(lines) + "\n")
print(f"{n} atomos pesados -> {out}")
