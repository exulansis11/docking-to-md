#!/usr/bin/env python3
"""Coloca uma molecula do ATB (ITP) nas coordenadas de uma pose/estrutura de origem (PDB ou PDBQT).

- Atomos pesados: isomorfismo de grafo (elemento + ligacoes) entre o ITP e a origem
  (ligacoes da origem por distancia). Funciona com nomes diferentes (FAD do cristal x FAD do ATB).
  Entre os isomorfismos possiveis (simetria), escolhe o que mais preserva os nomes.
- H do ITP (no united-atom so os polares): se a origem tem um H com o mesmo nome ligado ao
  mesmo atomo, usa essa coordenada (H polar da pose do AutoDock); senao, ajuste de Kabsch do
  fragmento local do PDB de referencia do ATB, ou geometria ideal se nao houver PDB.
- Escreve .gro (nm) na ordem do ITP, uma copia do ITP com #ifdef POSRES e o posre (pesados).

A origem pode ser a pose do AutoDock (PDBQT: tipos AD4 -> elementos) ou um PDB (ex. cofator do cristal).

Uso: place_molecule.py ITP SOURCE OUT_PREFIX [--ref ATB_PDB] [--resid N] [--select CHAIN:RESNUM] [--gro-only]
"""
import argparse
import re
import sys
from pathlib import Path

import networkx as nx
import numpy as np
from networkx.algorithms import isomorphism


def element_of_type(t):
    t = t.upper()
    if t.startswith("CL"):
        return "Cl"
    if t.startswith("BR"):
        return "Br"
    return t[0]


def read_itp(path):
    atoms, bonds, sect, moltype = [], [], None, None
    for raw in open(path, errors="replace"):
        s = raw.split(";")[0].strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith("["):
            sect = s.strip("[] ").lower()
            continue
        f = s.split()
        if sect == "moleculetype" and moltype is None:
            moltype = f[0]
        elif sect == "atoms":
            atoms.append(dict(nr=int(f[0]), type=f[1], resname=f[3], name=f[4], q=float(f[6]),
                              el=element_of_type(f[1])))
        elif sect == "bonds":
            bonds.append((int(f[0]), int(f[1])))
    return moltype, atoms, bonds


AD4_ELEMENT = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "O": "O",
               "SA": "S", "S": "S", "HD": "H", "HS": "H", "H": "H", "P": "P", "F": "F", "Cl": "Cl",
               "CL": "Cl", "Br": "Br", "BR": "Br", "I": "I"}


def read_pdb(path, select=None):
    out = []
    pdbqt = str(path).lower().endswith(".pdbqt")
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if select and (l[21], l[22:26].strip()) != select:
            continue
        name = l[12:16].strip()
        if pdbqt:
            el = AD4_ELEMENT.get(l[77:79].strip(), name.lstrip("0123456789")[0])
        else:
            el = l[76:78].strip() or name.lstrip("0123456789")[0]
            el = el.capitalize() if el.upper() in ("CL", "BR") else el.upper()
        out.append(dict(name=name, el=el, xyz=np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def dist_graph(atoms, heavy_only=True):
    g = nx.Graph()
    idx = [i for i, a in enumerate(atoms) if not heavy_only or a["el"] != "H"]
    for i in idx:
        g.add_node(i, el=atoms[i]["el"])
    X = np.array([atoms[i]["xyz"] for i in idx])
    D = np.linalg.norm(X[:, None] - X[None], axis=2)
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            cut = 1.25 if "H" in (atoms[idx[a]]["el"], atoms[idx[b]]["el"]) else 1.95
            if D[a, b] < cut:
                g.add_edge(idx[a], idx[b])
    return g


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def ideal_h(parent, nbrs, nh, k, ref_far):
    """posicao do k-esimo de nh H num atomo pai com vizinhos pesados nbrs (geometria aproximada)."""
    b = 1.0
    u = [(n - parent) / np.linalg.norm(n - parent) for n in nbrs]
    if len(u) >= 2:
        d = -np.sum(u, axis=0)
        if nh == 1 or np.linalg.norm(d) < 1e-3:
            return parent + b * d / np.linalg.norm(d)
        perp = np.cross(u[0], u[1]); perp /= np.linalg.norm(perp)
        d /= np.linalg.norm(d)
        return parent + b * (np.cos(0.955) * d + (1 if k == 0 else -1) * np.sin(0.955) * perp)
    # um vizinho: H a ~109.5 graus, girando em torno do eixo pai-vizinho
    axis = u[0]
    ref = ref_far - parent if ref_far is not None else np.array([1.0, 0, 0])
    ref = ref - np.dot(ref, axis) * axis
    if np.linalg.norm(ref) < 1e-3:
        ref = np.cross(axis, [0, 0, 1.0])
    ref /= np.linalg.norm(ref)
    third = np.cross(axis, ref)
    phi = np.pi + 2 * np.pi * k / nh
    ang = np.deg2rad(109.5)
    return parent + b * (np.cos(ang) * axis + np.sin(ang) * (np.cos(phi) * ref + np.sin(phi) * third))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("itp"); ap.add_argument("source"); ap.add_argument("out")
    ap.add_argument("--ref"); ap.add_argument("--select"); ap.add_argument("--resid", type=int, default=1)
    ap.add_argument("--gro-only", action="store_true", help="nao escreve ITP/posre (copias extras da mesma molecula)")
    a = ap.parse_args()
    sel = tuple(a.select.split(":")) if a.select else None

    moltype, T, bonds = read_itp(a.itp)
    S = read_pdb(a.source, sel)
    if not S:
        sys.exit("origem vazia")
    tg = nx.Graph()
    for t in T:
        tg.add_node(t["nr"], el=t["el"])
    tg.add_edges_from(bonds)
    th = tg.subgraph([t["nr"] for t in T if t["el"] != "H"]).copy()
    sh = dist_graph(S)
    if th.number_of_nodes() != sh.number_of_nodes() or th.number_of_edges() != sh.number_of_edges():
        sys.exit(f"[{moltype}] pesados/ligacoes diferem: ITP {th.number_of_nodes()}/{th.number_of_edges()} "
                 f"origem {sh.number_of_nodes()}/{sh.number_of_edges()}")
    gm = isomorphism.GraphMatcher(th, sh, node_match=lambda x, y: x["el"] == y["el"])
    byname = {t["nr"]: t["name"] for t in T}
    best, best_score = None, -1
    for n, m in enumerate(gm.isomorphisms_iter()):
        score = sum(byname[i] == S[j]["name"] for i, j in m.items())
        if score > best_score:
            best, best_score = m, score
        if score == len(m) or n > 20000:
            break
    if best is None:
        sys.exit(f"[{moltype}] sem isomorfismo entre ITP e origem")

    X = {i: S[j]["xyz"] for i, j in best.items()}
    # H polares da origem (mesmo nome e ligados ao atomo correspondente)
    sall = dist_graph(S, heavy_only=False)
    ref = read_pdb(a.ref) if a.ref else None
    refxyz = {r["name"]: r["xyz"] for r in ref} if ref else {}
    from_src = from_ref = ideal = 0
    hs = [t for t in T if t["el"] == "H"]
    for t in hs:
        parent = next(iter(tg[t["nr"]]))
        sp = best[parent]
        cand = [j for j in sall[sp] if S[j]["el"] == "H" and S[j]["name"] == t["name"]]
        if cand:
            X[t["nr"]] = S[cand[0]]["xyz"]; from_src += 1
            continue
        heavy_n = [n for n in th[parent]]
        frame = [parent] + heavy_n
        if len(frame) < 3:
            frame += [n2 for n in heavy_n for n2 in th[n] if n2 not in frame]
        if ref and all(byname[f] in refxyz for f in frame) and t["name"] in refxyz:
            P = np.array([refxyz[byname[f]] for f in frame]); Q = np.array([X[f] for f in frame])
            R, tr = kabsch(P, Q)
            X[t["nr"]] = R @ refxyz[t["name"]] + tr; from_ref += 1
        else:
            sib = [h for h in tg[parent] if tg.nodes[h]["el"] == "H"]
            far = [n2 for n in heavy_n for n2 in th[n] if n2 != parent]
            X[t["nr"]] = ideal_h(X[parent], [X[n] for n in heavy_n], len(sib), sib.index(t["nr"]),
                                 X[far[0]] if far else None); ideal += 1

    out = Path(a.out)
    with open(f"{out}.gro", "w") as fh:
        fh.write(f"{moltype}\n{len(T)}\n")
        for t in T:
            x = X[t["nr"]] / 10
            fh.write(f"{a.resid:>5}{t['resname'][:5]:<5}{t['name'][:5]:>5}{t['nr'] % 100000:>5}"
                     f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}\n")
        fh.write("   0.00000   0.00000   0.00000\n")
    if a.gro_only:
        q = sum(t["q"] for t in T)
        print(f"{moltype}\tatoms={len(T)}\theavy={th.number_of_nodes()}\tH_from_source={from_src}\t"
              f"H_from_ATBref={from_ref}\tH_ideal={ideal}\tnames_kept={best_score}/{len(best)}\tcharge={q:.3f}")
        return
    with open(f"{out}_posre.itp", "w") as fh:
        fh.write(f"; restricao de posicao dos atomos pesados de {moltype}\n[ position_restraints ]\n")
        for t in T:
            if t["el"] != "H":
                fh.write(f"{t['nr']:>6}     1  1000  1000  1000\n")
    itp_txt = open(a.itp, errors="replace").read().replace("\r", "")
    # troca qualquer bloco POSRES que ja venha no ITP pelo nosso
    itp_txt = re.sub(r"#ifdef POSRES[^\n]*\n.*?#endif[^\n]*\n?", "", itp_txt, flags=re.S)
    itp_txt = itp_txt.rstrip() + f'\n\n#ifdef POSRES\n#include "{out.name}_posre.itp"\n#endif\n'
    Path(f"{out}.itp").write_text(itp_txt)

    q = sum(t["q"] for t in T)
    print(f"{moltype}\tatoms={len(T)}\theavy={th.number_of_nodes()}\tH_from_source={from_src}\t"
          f"H_from_ATBref={from_ref}\tH_ideal={ideal}\tnames_kept={best_score}/{len(best)}\tcharge={q:.3f}")


if __name__ == "__main__":
    main()
