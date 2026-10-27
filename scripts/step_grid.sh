#!/usr/bin/env bash
# grid: GPF (prepare_gpf4) + AutoGrid. Um conjunto de mapas cobre os tipos de todos os ligantes.

grid_center() {
    # GRID_CENTER explicito, ou centroide de residuos do receptor, ou de um ligante de referencia
    if [ -n "$GRID_CENTER" ]; then
        printf '%s' "$GRID_CENTER"; return
    fi
    "$PYTHON" - "$OUT/01_prepared/$REC.pdbqt" "$GRID_CENTER_RESIDUES" "$GRID_CENTER_LIGAND" <<'PY'
import sys
rec, residues, ref = sys.argv[1], sys.argv[2].split(), sys.argv[3]
xyz = []
if residues:
    want = {(r.split(":")[0], r.split(":")[1]) for r in residues}
    for l in open(rec):
        if l.startswith(("ATOM", "HETATM")) and (l[21], l[22:26].strip()) in want:
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
elif ref:
    for l in open(ref, errors="replace"):
        if l.startswith(("ATOM", "HETATM")):
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
if not xyz:
    sys.exit("defina GRID_CENTER, GRID_CENTER_RESIDUES ou GRID_CENTER_LIGAND")
c = [sum(v[i] for v in xyz) / len(xyz) for i in range(3)]
print("%.3f,%.3f,%.3f" % tuple(c))
PY
}

step_grid() {
    local D="$OUT/02_grid"
    mkdir -p "$D"; cd "$D"
    cp "$OUT/01_prepared/$REC.pdbqt" .

    local types center
    types=$(for L in "${LIGANDS[@]}"; do awk '/^(ATOM|HETATM)/{print $NF}' "$OUT/01_prepared/$L.pdbqt"; done | sort -u | paste -sd,)
    center=$(grid_center)
    log "caixa: centro $center A, $GRID_NPTS pontos, espacamento $GRID_SPACING A; tipos $types"

    "$PYTHONSH" "$UTIL/prepare_gpf4.py" -r "$REC.pdbqt" -o "$REC.gpf" \
        -p ligand_types="$types" -p npts="$GRID_NPTS" -p gridcenter="$center" -p spacing="$GRID_SPACING" \
        > prepare_gpf.log 2>&1 || { tail prepare_gpf.log; die "prepare_gpf4 falhou"; }

    "$AUTOGRID4" -p "$REC.gpf" -l "$REC.glg"
    grep -q "Successful Completion" "$REC.glg" || die "AutoGrid falhou (ver $D/$REC.glg)"
    log "AutoGrid ok"
    cd - >/dev/null
}
