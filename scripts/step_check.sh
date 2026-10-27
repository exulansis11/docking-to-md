#!/usr/bin/env bash
# check: software e arquivos de entrada

step_check() {
    activate_env
    local ok=1

    if [ -n "$MGLTOOLS" ] && [ -x "$PYTHONSH" ]; then
        local t; t=$(mktemp -d)
        echo 'import MolKit, AutoDockTools; print "ok"' > "$t/t.py"
        if "$PYTHONSH" "$t/t.py" >/dev/null 2>&1; then log "MGLTools: ok ($MGLTOOLS)"; else warn "MGLTools: pythonsh nao importa MolKit"; ok=0; fi
        rm -rf "$t"
        for s in prepare_receptor4 prepare_ligand4 prepare_gpf4 prepare_dpf4; do
            [ -f "$UTIL/$s.py" ] || { warn "falta $UTIL/$s.py"; ok=0; }
        done
    else
        warn "MGLTOOLS nao definido ou sem bin/pythonsh (necessario para prepare/grid/dock)"; ok=0
    fi

    local c
    for c in "$AUTOGRID4" "$AUTODOCK4" "$PYTHON"; do
        command -v "$c" >/dev/null 2>&1 && log "$c: $(command -v "$c")" || { warn "nao encontrado: $c"; ok=0; }
    done
    "$PYTHON" -c "import numpy, networkx" 2>/dev/null && log "python: numpy e networkx ok" \
        || { warn "python precisa de numpy e networkx"; ok=0; }
    command -v "$GMX" >/dev/null 2>&1 && log "GROMACS: $("$GMX" --version 2>/dev/null | awk '/GROMACS version/{print $3; exit}')" \
        || warn "GROMACS ($GMX) nao encontrado: so necessario para a etapa md"
    if [ "$PROTONATION" = propka ]; then
        command -v "$PROPKA" >/dev/null 2>&1 && log "PROPKA: $("$PROPKA" --version 2>&1 | head -1) (pH $PH)" \
            || { warn "PROPKA ($PROPKA) nao encontrado: necessario com PROTONATION=propka"; ok=0; }
    fi
    log "MD: NaCl $SALT_CONC M, protonacao $PROTONATION, agua $WATER, campo de forca $FF_NAME"
    command -v "$PYMOL" >/dev/null 2>&1 || warn "PyMOL ($PYMOL) nao encontrado: so necessario para a etapa figures"

    [ -s "$RECEPTOR" ] && log "receptor: $RECEPTOR" || { warn "receptor ausente: $RECEPTOR"; ok=0; }
    local L
    for L in "${LIGANDS[@]}"; do
        [ -s "${LIG_PDB[$L]:-}" ] || { warn "$L: PDB ausente (${LIG_PDB[$L]:-nao definido})"; ok=0; }
        [ -s "${LIG_ITP[$L]:-}" ] || warn "$L: ITP ausente (so necessario para a etapa md)"
    done
    local entry name
    for entry in "${COFACTORS[@]}"; do
        name=${entry%%:*}
        [ -s "${COFACTOR_ITP[$name]:-}" ] || warn "cofator $name: ITP ausente"
    done
    [ ${#COFACTORS[@]} -eq 0 ] || [ -s "$COFACTOR_PDB" ] || warn "COFACTOR_PDB ausente"
    [ -z "$FF_DIR" ] || [ -s "$FF_DIR/forcefield.itp" ] || warn "FF_DIR sem forcefield.itp: $FF_DIR"
    for m in empr em cg mdpr298 mdnvt298; do [ -s "$MDP_DIR/$m.mdp" ] || warn "falta $MDP_DIR/$m.mdp"; done

    [ $ok -eq 1 ] && log "check: ok" || die "check falhou (ver avisos acima)"
}
