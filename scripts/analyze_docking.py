#!/usr/bin/env python3
"""Le os DLG do AutoDock4, monta a tabela de clusters e escolhe a pose de cada ligante.

Regra de escolha: entre os clusters cujo menor dG esta a ate DELTA kcal/mol do minimo global,
o mais populoso (empate -> menor energia); a pose e a de menor energia desse cluster.
O erro padrao da funcao do AD4 e ~2 kcal/mol: clusters a menos de DELTA do minimo sao tratados
como energeticamente equivalentes e o tamanho do cluster desempata.

Saidas em <out>/04_analysis/: docking_summary.tsv, <L>_clusters.tsv, <L>_best.pdbqt,
<L>_cluster_reps.pdbqt (multi-model, um representante por cluster).

Uso: analyze_docking.py OUT RECEPTOR_PDBQT DELTA REFERENCE_DIR|- LIG1 [LIG2 ...]
"""
import csv
import re
import sys
from pathlib import Path

import numpy as np


def parse_dlg(path):
    runs, cur = {}, None
    clusters, members = [], {}
    sect = None
    for raw in open(path, errors="replace"):
        line = raw.rstrip("\r\n")
        if line.startswith("DOCKED: MODEL"):
            cur = {"atoms": []}
        elif cur is not None and line.startswith("DOCKED: USER    Run ="):
            cur["run"] = int(line.split("=")[1])
        elif cur is not None and "Estimated Free Energy of Binding" in line:
            cur["dG"] = float(line.split("=")[1].split()[0])
        elif cur is not None and "Estimated Inhibition Constant" in line:
            m = re.search(r"=\s*([\d.]+)\s*(\w+)", line)
            cur["Ki"] = f"{m.group(1)} {m.group(2)}" if m else "NA"
        elif cur is not None and line.startswith("DOCKED: ATOM"):
            cur["atoms"].append(line[8:])
        elif cur is not None and line.startswith("DOCKED: ENDMDL"):
            cur.setdefault("Ki", "NA")
            runs[cur["run"]] = cur
            cur = None
        elif "CLUSTERING HISTOGRAM" in line:
            sect = "hist"
        elif "RMSD TABLE" in line:
            sect = "rmsd"
        elif sect == "hist":
            m = re.match(r"\s*(\d+)\s*\|\s*([-\d.]+)\s*\|\s*(\d+)\s*\|\s*([-\d.]+)\s*\|\s*(\d+)\s*\|", line)
            if m:
                clusters.append(dict(rank=int(m.group(1)), best=float(m.group(2)), run=int(m.group(3)),
                                     mean=float(m.group(4)), n=int(m.group(5))))
        elif sect == "rmsd":
            m = re.match(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)\s+RANKING", line)
            if m:
                members.setdefault(int(m.group(1)), []).append(int(m.group(3)))
    return runs, clusters, members


def xyz(atom_lines):
    return {l[12:16].strip(): np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in atom_lines}


def rmsd(a, b):
    common = [k for k in a if k in b and not k.startswith("H")]
    return float(np.sqrt(np.mean([np.sum((a[k] - b[k]) ** 2) for k in common]))) if common else float("nan")


def contacts(lig, rec_lines, cut=4.0):
    L = np.array(list(lig.values()))
    res = set()
    for l in rec_lines:
        p = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if np.min(np.linalg.norm(L - p, axis=1)) <= cut:
            res.add((l[21], int(l[22:26]), l[17:20]))
    return sorted(res, key=lambda r: (r[0], r[1]))


def main():
    out, rec_pdbqt, delta, refdir = Path(sys.argv[1]), sys.argv[2], float(sys.argv[3]), sys.argv[4]
    ligs = sys.argv[5:]
    ana = out / "04_analysis"
    ana.mkdir(exist_ok=True)
    rec = [l for l in open(rec_pdbqt) if l.startswith(("ATOM", "HETATM"))]

    summary = []
    for L in ligs:
        runs, clusters, members = parse_dlg(out / "03_docking" / L / f"{L}.dlg")
        if not clusters:
            sys.exit(f"{L}: DLG sem histograma de clusters")
        nruns = len(runs)
        gmin = min(c["best"] for c in clusters)
        cand = [c for c in clusters if c["best"] - gmin <= delta + 1e-9]
        sel = sorted(cand, key=lambda c: (-c["n"], c["best"]))[0]
        pose = runs[sel["run"]]

        with open(ana / f"{L}_clusters.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["cluster_rank", "lowest_dG_kcal_mol", "run_of_lowest", "mean_dG_kcal_mol", "n_members",
                        "pct", "Ki_of_lowest", "runs"])
            for c in clusters:
                w.writerow([c["rank"], c["best"], c["run"], c["mean"], c["n"], f"{100 * c['n'] / nruns:.0f}",
                            runs[c["run"]]["Ki"], ",".join(map(str, members.get(c["rank"], [])))])

        (ana / f"{L}_best.pdbqt").write_text(
            f"REMARK selected run {sel['run']} cluster {sel['rank']} dG {pose['dG']} Ki {pose['Ki']}\n"
            + "".join(a + "\n" for a in pose["atoms"]) + "TER\n")
        with open(ana / f"{L}_cluster_reps.pdbqt", "w") as fh:
            for c in clusters:
                fh.write(f"MODEL {c['rank']}\nREMARK cluster {c['rank']} run {c['run']} dG {c['best']} n {c['n']}\n")
                fh.writelines(a + "\n" for a in runs[c["run"]]["atoms"])
                fh.write("ENDMDL\n")

        new = xyz(pose["atoms"])
        ref_rmsd = "NA"
        if refdir != "-":
            f = Path(refdir) / f"{L}.pdbqt"
            if f.exists():
                old = xyz([l for l in open(f, errors="replace") if l.startswith(("ATOM", "HETATM"))])
                ref_rmsd = f"{rmsd(new, old):.2f}"
        top = sorted(clusters, key=lambda c: c["best"])[0]
        summary.append(dict(
            ligand=L, n_runs=nruns, n_clusters=len(clusters),
            selected_cluster=sel["rank"], selected_run=sel["run"],
            dG_kcal_mol=pose["dG"], Ki=pose["Ki"], mean_dG_cluster=sel["mean"],
            cluster_size=sel["n"], cluster_pct=f"{100 * sel['n'] / nruns:.0f}",
            global_min_dG=gmin, global_min_run=top["run"], global_min_cluster_size=top["n"],
            rmsd_vs_reference_A=ref_rmsd,
            residues_4A=" ".join(f"{r[2]}{r[1]}{r[0]}" for r in contacts(new, rec)),
        ))

    with open(ana / "docking_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    print(f"{'ligante':14s} {'dG':>7s} {'Ki':>12s} {'cluster':>9s} {'run':>4s} {'min.global':>11s} {'RMSD ref':>8s}")
    for s in summary:
        print(f"{s['ligand']:14s} {s['dG_kcal_mol']:7.2f} {s['Ki']:>12s} {s['cluster_size']:>4d}/{s['n_runs']:<4d} "
              f"{s['selected_run']:4d} {s['global_min_dG']:11.2f} {s['rmsd_vs_reference_A']:>8s}")


if __name__ == "__main__":
    main()
