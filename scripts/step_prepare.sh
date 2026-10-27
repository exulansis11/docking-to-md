#!/usr/bin/env bash
# prepare: receptor e ligantes -> PDBQT (H polares, cargas de Gasteiger, tipos AD4, torsoes)

step_prepare() {
    local D="$OUT/01_prepared"
    mkdir -p "$D"; cd "$D"

    # arquivos do Windows (CRLF) confundem o MolKit
    local ext="${RECEPTOR##*.}"
    copy_lf "$RECEPTOR" "${REC}_in.$ext"
    # Atencao: o MolKit so liga H com nomes no padrao PDB v2/ADT (ex. 1HH1). PDBs com H renomeados
    # (HH11...) ficam com H soltos e cargas erradas: confira a carga total no prepare_audit.tsv.
    # shellcheck disable=SC2086
    "$PYTHONSH" "$UTIL/prepare_receptor4.py" -r "${REC}_in.$ext" -o "$REC.pdbqt" $RECEPTOR_PREP_FLAGS \
        > prepare_receptor.log 2>&1 || { tail -20 prepare_receptor.log; die "prepare_receptor4 falhou"; }

    local L
    for L in "${LIGANDS[@]}"; do
        copy_lf "${LIG_PDB[$L]}" "${L}_atb.pdb"
        # sem reprotonar: a protonacao e a do PDB (ex. ATB). Padrao: Gasteiger, H apolares fundidos,
        # amidas fixas. Rodar na propria pasta (o MolKit antigo descarta diretorios do nome).
        "$PYTHONSH" "$UTIL/prepare_ligand4.py" -l "${L}_atb.pdb" -o "$L.pdbqt" > "prepare_$L.log" 2>&1 \
            || { tail -20 "prepare_$L.log"; die "prepare_ligand4 falhou: $L"; }
    done

    {
        printf "molecule\tatoms\tgasteiger_charge\texpected_charge\tTORSDOF\tatom_types\n"
        local M
        for M in "$REC" "${LIGANDS[@]}"; do
            LC_ALL=C awk -v m="$M" -v e="${LIG_CHARGE[$M]:-NA}" '
                /^(ATOM|HETATM)/ {n++; q+=substr($0,71,6); t[$NF]=1}
                /^TORSDOF/ {tor=$2}
                END {s=""; for (k in t) s=s k " "; printf "%s\t%d\t%.3f\t%s\t%s\t%s\n", m, n, q, e, (tor==""?"-":tor), s}' "$M.pdbqt"
        done
    } | tee prepare_audit.tsv

    # carga Gasteiger dos ligantes deve bater com a esperada (se informada)
    for L in "${LIGANDS[@]}"; do
        [ -n "${LIG_CHARGE[$L]:-}" ] || continue
        LC_ALL=C awk -v e="${LIG_CHARGE[$L]}" -v l="$L" -F'\t' '$1==l {d=$3-e; if (d<0) d=-d; if (d>0.1) {print "carga de " l " = " $3 ", esperada " e; exit 1}}' prepare_audit.tsv \
            || die "carga do ligante nao confere"
    done
    cd - >/dev/null
}
