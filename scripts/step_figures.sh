#!/usr/bin/env bash
# figures: imagens das poses no receptor (PyMOL)

step_figures() {
    activate_env
    need_cmd "$PYMOL"
    local D="$OUT/07_figures" L args=()
    mkdir -p "$D"
    for L in "${LIGANDS[@]}"; do args+=("$L=$OUT/04_analysis/${L}_best.pdbqt"); done
    # homodimero: copia simetrica montada na etapa md
    if [ -n "$LIG_SYMMETRY" ]; then
        for L in "${LIGANDS[@]}"; do
            [ -s "$OUT/05_md/$L/lig2.gro" ] && args+=("${L}_sitio2=$OUT/05_md/$L/lig2.gro")
        done
    fi
    "$PYMOL" -cq "$HERE/scripts/figures.py" -- "$D" "$OUT/02_grid/$REC.pdbqt" "${COFACTOR_PDB:--}" "${args[@]}" \
        > "$D/pymol.log" 2>&1 || { tail "$D/pymol.log"; die "PyMOL falhou"; }
    log "figuras: $(ls "$D"/*.png | wc -l) PNG em $D"
}
