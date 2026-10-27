#!/usr/bin/env bash
# analyze: clusters, escolha da pose e tabela final

step_analyze() {
    "$PYTHON" "$HERE/scripts/analyze_docking.py" "$OUT" "$OUT/02_grid/$REC.pdbqt" "$SELECT_DELTA" \
        "${REFERENCE_POSES:--}" "${LIGANDS[@]}"
    log "tabela: $OUT/04_analysis/docking_summary.tsv"
}
