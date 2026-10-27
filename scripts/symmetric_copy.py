#!/usr/bin/env python3
"""Copia o ligante para o sitio equivalente de um homodimero (operacao de simetria A<->B).

A operacao vem do ajuste (Kabsch) dos CA comuns: [A; B] -> [B; A], que para um dimero
simetrico e uma rotacao de ~180 graus. Ela e aplicada a todos os atomos do ligante (GRO, nm).
Audita: RMSD do ajuste, angulo, menor distancia da copia a proteina/cofatores/1o ligante e
os residuos a 4 A de cada copia (que devem se corresponder com as cadeias trocadas).

Choques: com --fix-out, cada residuo com atomo pesado a menos de --clash A da copia recebe a cadeia
lateral (CB em diante) da imagem simetrica do residuo equivalente da outra cadeia, ou seja, o 2o sitio
fica igual ao sitio onde a pose foi obtida (tipicamente uma Lys/Arg com rotamero diferente no cristal).
O receptor corrigido (so atomos pesados) e escrito em --fix-out, para o pdb2gmx.

Uso: symmetric_copy.py RECEPTOR(pdb|pdbqt) A:B LIG_IN.gro LIG_OUT.gro --resid N
         [--fix-out RECEPTOR_FIXED.pdb] [--system COMPLEX.gro]
"""
import argparse
import sys

import numpy as np

AD4_H = {"HD", "HS", "H"}


def read_receptor(path):
    pdbqt = path.lower().endswith(".pdbqt")
    atoms = []
    for l in open(path, errors="replace"):
        if not l.startswith("ATOM"):
            continue
        name = l[12:16].strip()
        is_h = (l[77:79].strip() in AD4_H) if pdbqt else ((l[76:78].strip() or name[0]) == "H")
        el = l[76:78].strip() if not pdbqt else ""
        atoms.append(dict(name=name, chain=l[21], resi=int(l[22:26]), resn=l[17:20], h=is_h, el=el,
                          xyz=np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return atoms


def read_gro(path):
    L = open(path).read().splitlines()
    n = int(L[1])
    out = []
    for l in L[2:2 + n]:
        out.append(dict(line=l, name=l[10:15].strip(), resn=l[5:10].strip(),
                        xyz=np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]) * 10.0))
    return out, L[2 + n]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def is_heavy_name(n):
    return not n.lstrip("0123456789").upper().startswith("H")


def site(lig_xyz, rec, cut=4.0):
    heavy = [a for a in rec if not a["h"]]
    X = np.array([a["xyz"] for a in heavy])
    d = np.linalg.norm(X[:, None] - lig_xyz[None], axis=2).min(axis=1)
    return sorted({(a["chain"], a["resi"], a["resn"]) for a, di in zip(heavy, d) if di <= cut})


def swap_chain(c, c1, c2):
    return c2 if c == c1 else c1


def clashing(lig_xyz, rec, cut):
    bad = set()
    for a in rec:
        if a["h"] or a["name"] in ("N", "CA", "C", "O"):
            continue
        if np.linalg.norm(lig_xyz - a["xyz"], axis=1).min() < cut:
            bad.add((a["chain"], a["resi"]))
    return sorted(bad)


def symmetrize_sidechain(rec, chain, resi, partner, R, t):
    src = {a["name"]: a for a in rec if a["chain"] == partner and a["resi"] == resi}
    n = 0
    for a in rec:
        if a["chain"] == chain and a["resi"] == resi and a["name"] not in ("N", "CA", "C", "O") and a["name"] in src:
            if src[a["name"]]["resn"] != a["resn"]:
                sys.exit(f"{chain}{resi} e {partner}{resi} sao residuos diferentes")
            a["xyz"] = R @ src[a["name"]]["xyz"] + t
            n += 1
    return n


def write_receptor(rec, path):
    lines, last, n = [], None, 0
    for a in rec:
        if a["h"]:
            continue
        if last is not None and a["chain"] != last:
            lines.append("TER")
        last = a["chain"]
        n += 1
        nm = a["name"] if len(a["name"]) == 4 else " " + a["name"]
        el = a.get("el") or a["name"].lstrip("0123456789")[0]
        x = a["xyz"]
        lines.append(f"ATOM  {n % 100000:>5} {nm:<4} {a['resn']:>3} {a['chain']}{a['resi']:>4}    "
                     f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00          {el:>2}")
    lines += ["TER", "END"]
    open(path, "w").write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("receptor"); ap.add_argument("pair"); ap.add_argument("lig_in"); ap.add_argument("lig_out")
    ap.add_argument("--resid", type=int, default=1)
    ap.add_argument("--system", help="GRO com proteina+cofatores+1o ligante, para checar choques")
    ap.add_argument("--min-dist", type=float, default=1.5, help="distancia minima pesado-pesado aceita (A)")
    ap.add_argument("--clash", type=float, default=2.2, help="corte para trocar a cadeia lateral (A)")
    ap.add_argument("--fix-out", help="escreve o receptor (pesados) com as cadeias laterais simetrizadas")
    a = ap.parse_args()
    c1, c2 = a.pair.split(":")

    rec = read_receptor(a.receptor)
    ca = {(x["chain"], x["resi"]): x["xyz"] for x in rec if x["name"] == "CA"}
    common = sorted({r for c, r in ca if c == c1} & {r for c, r in ca if c == c2})
    if len(common) < 10:
        sys.exit(f"poucos CA em comum entre {c1} e {c2}: {len(common)}")
    P = np.array([ca[(c1, r)] for r in common] + [ca[(c2, r)] for r in common])
    Q = np.array([ca[(c2, r)] for r in common] + [ca[(c1, r)] for r in common])
    R, t = kabsch(P, Q)
    res = np.linalg.norm(P @ R.T + t - Q, axis=1)
    rms = float(np.sqrt((res ** 2).mean()))
    ang = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
    print(f"simetria {c1}<->{c2}: {len(common)} CA, RMSD {rms:.2f} A, rotacao {ang:.1f} graus")
    if rms > 2.0 or abs(ang - 180) > 10:
        sys.exit("o dimero nao parece simetrico o bastante para copiar o ligante (RMSD > 2 A ou rotacao longe de 180)")

    lig, box = read_gro(a.lig_in)
    X = np.array([x["xyz"] for x in lig])
    Y = X @ R.T + t
    with open(a.lig_out, "w") as fh:
        fh.write(f"copia simetrica ({c1}<->{c2})\n{len(lig)}\n")
        for x, y in zip(lig, Y / 10.0):
            l = x["line"]
            fh.write(f"{a.resid:>5}{l[5:20]}{y[0]:8.3f}{y[1]:8.3f}{y[2]:8.3f}\n")
        fh.write(box + "\n")

    heavy = np.array([is_heavy_name(x["name"]) for x in lig])

    if a.fix_out:
        fixed = set()
        for _ in range(5):
            bad = clashing(Y[heavy], rec, a.clash)
            if not bad:
                break
            for c, r in bad:
                if (c, r) in fixed:
                    continue
                n = symmetrize_sidechain(rec, c, r, swap_chain(c, c1, c2), R, t)
                fixed.add((c, r))
                print(f"cadeia lateral de {c}{r} trocada pela imagem simetrica de {swap_chain(c, c1, c2)}{r} ({n} atomos)")
        bad = clashing(Y[heavy], rec, a.clash)
        if bad:
            print("ainda perto da copia: " + " ".join(f"{c}{r}" for c, r in bad))
        write_receptor(rec, a.fix_out)

    s1, s2 = site(X[heavy], rec), site(Y[heavy], rec)
    swap = {c1: c2, c2: c1}
    mapped = {(swap[c], r) for c, r, _ in s1}
    got = {(c, r) for c, r, _ in s2}
    print(f"sitio original ({len(s1)} res): " + " ".join(f"{n}{r}{c}" for c, r, n in s1))
    print(f"sitio da copia ({len(s2)} res): " + " ".join(f"{n}{r}{c}" for c, r, n in s2))
    print(f"correspondencia dos residuos com cadeias trocadas: {len(mapped & got)}/{len(mapped | got)}")

    if a.system:
        sysat, _ = read_gro(a.system)
        S = np.array([x["xyz"] for x in sysat])
        sh = np.array([is_heavy_name(x["name"]) for x in sysat])
        dh = np.linalg.norm(S[sh][:, None] - Y[heavy][None], axis=2).min()
        da = np.linalg.norm(S[:, None] - Y[None], axis=2).min()
        dl = np.linalg.norm(X[:, None] - Y[None], axis=2).min()
        print(f"menor distancia da copia: pesado-pesado {dh:.2f} A, qualquer {da:.2f} A, ao 1o ligante {dl:.2f} A")
        if dh < a.min_dist:
            sys.exit(f"choque: copia a {dh:.2f} A de um atomo pesado do sistema (< {a.min_dist} A)")


if __name__ == "__main__":
    main()
