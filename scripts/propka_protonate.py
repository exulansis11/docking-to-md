#!/usr/bin/env python3
"""Estados de protonacao pelo PROPKA (pKa no microambiente da proteina) para um dado pH.

Le o resumo do .pka e renomeia os residuos do PDB (so pesados) com os nomes que o pdb2gmx do
GROMOS 54A7 entende (via aminoacids.r2b), sem perguntas interativas:
  ASP  pKa > pH -> ASPH (neutro)      senao ASP (-1)
  GLU  pKa > pH -> GLUH (neutro)      senao GLU (-1)
  LYS  pKa < pH -> LYSN (neutro)      senao LYS (+1)
  ARG  pKa < pH -> ARGN (neutro)      senao ARG (+1)
  HIS  pKa > pH -> HISH (+1)          senao neutra: tautomero escolhido pelo pdb2gmx (rede de H-bond),
                                       lido do log de uma rodada automatica (HISD/HISE)
CYS/TYR desprotonadas e termini neutros nao tem bloco no GROMOS 54A7: so aparecem como aviso.
Cisteinas em ponte (pKa 99.99 no PROPKA) sao ignoradas.

Uso: propka_protonate.py PDB_PESADOS ARQUIVO.pka PH LOG_PDB2GMX_AUTO PDB_SAIDA RELATORIO.tsv
"""
import re
import sys


def read_pka(path):
    out, on = [], False
    for l in open(path, errors="replace"):
        if "SUMMARY OF THIS PREDICTION" in l:
            on = True
            continue
        if on:
            if l.startswith("-") or not l.strip():
                if out:
                    break
                continue
            f = l.split()
            if len(f) >= 5 and re.fullmatch(r"-?\d+(\.\d+)?", f[3]):
                out.append(dict(res=f[0], num=int(f[1]), chain=f[2], pka=float(f[3]), model=float(f[4])))
    if not out:
        sys.exit(f"resumo nao encontrado em {path}")
    return out


def his_auto(log, pdb):
    """{(cadeia, resnum): HISD|HISE|HISH} da rodada automatica do pdb2gmx.
    O log lista 'Will use HISx for residue N' na ordem das His do PDB, sem dizer a cadeia;
    por isso as escolhas sao associadas, em ordem, as His do PDB de entrada."""
    his = []
    for l in open(pdb):
        if l.startswith("ATOM") and l[17:20] == "HIS" and l[12:16].strip() == "CA":
            his.append((l[21], int(l[22:26])))
    picks = [(m.group(1), int(m.group(2))) for m in
             (re.search(r"Will use (HIS[DEH]) for residue (\d+)", l) for l in open(log, errors="replace")) if m]
    if len(picks) != len(his) or any(p[1] != h[1] for p, h in zip(picks, his)):
        sys.exit(f"log do pdb2gmx ({len(picks)} His) nao bate com o PDB ({len(his)} His)")
    return {h: p[0] for h, p in zip(his, picks)}


def main():
    pdb, pka_file, ph, log, out_pdb, report = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6]
    pka = read_pka(pka_file)
    auto = his_auto(log, pdb)
    new, rows, warns = {}, [], []
    for r in pka:
        key, name, k = (r["chain"], r["num"]), r["res"], r["pka"]
        state = None
        if k >= 99:
            continue                      # cisteina em ponte dissulfeto
        if name == "ASP":
            state = "ASPH" if k > ph else "ASP"
        elif name == "GLU":
            state = "GLUH" if k > ph else "GLU"
        elif name == "LYS":
            state = "LYSN" if k < ph else "LYS"
        elif name == "ARG":
            state = "ARGN" if k < ph else "ARG"
        elif name == "HIS":
            state = "HISH" if k > ph else auto.get(key, "HISE").replace("HISH", "HISE")
        elif name == "CYS" and k < ph:
            warns.append(f"CYS{r['num']}{r['chain']} pKa {k:.2f} < pH: GROMOS 54A7 nao tem tiolato; mantida neutra")
        elif name == "TYR" and k < ph:
            warns.append(f"TYR{r['num']}{r['chain']} pKa {k:.2f} < pH: sem tirosinato no GROMOS 54A7; mantida neutra")
        elif name == "N+" and k < ph:
            warns.append(f"N-terminal {r['num']}{r['chain']} pKa {k:.2f} < pH: neutro pelo PROPKA; o pdb2gmx usa NH3+ (mude com -ter)")
        elif name == "C-" and k > ph:
            warns.append(f"C-terminal {r['num']}{r['chain']} pKa {k:.2f} > pH: neutro pelo PROPKA; o pdb2gmx usa COO- (mude com -ter)")
        if state:
            new[key] = state
            default = {"ASP": "ASP", "GLU": "GLU", "LYS": "LYS", "ARG": "ARG"}.get(name, auto.get(key, "HIS?"))
            rows.append((r["chain"], r["num"], name, k, r["model"], state, "" if state == default else "ALTERADO"))

    n = 0
    with open(out_pdb, "w") as fh:
        for l in open(pdb):
            if l.startswith(("ATOM", "HETATM")):
                key = (l[21], int(l[22:26]))
                if key in new and l[17:20] == new[key][:3]:
                    l = l[:17] + f"{new[key]:<4}" + l[21:]
                    n += 1
            fh.write(l)

    with open(report, "w") as fh:
        fh.write(f"# pH {ph}; PROPKA: {pka_file}\nchain\tresnum\tres\tpKa\tmodel_pKa\testado\tvs_padrao_pdb2gmx\n")
        for row in rows:
            fh.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.2f}" for x in row) + "\n")
        for w in warns:
            fh.write(f"# AVISO: {w}\n")

    changed = [r for r in rows if r[6]]
    print(f"PROPKA pH {ph}: {len(rows)} residuos ionizaveis; {len(changed)} diferem do padrao do pdb2gmx")
    for c, num, name, k, _, st, _ in changed:
        print(f"    {name}{num}{c}: pKa {k:.2f} -> {st}")
    for w in warns:
        print(f"    AVISO: {w}")


if __name__ == "__main__":
    main()
