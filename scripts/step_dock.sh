#!/usr/bin/env bash
# dock: DPF (prepare_dpf4) + AutoDock4.
# O padrao do prepare_dpf4 ja e o Lamarckian GA (GALS). NAO usar -L: e so busca local.

dock_one() {
    local L=$1 D="$OUT/03_docking/$1"
    mkdir -p "$D"; cd "$D"
    cp "$OUT/01_prepared/$L.pdbqt" .
    local f
    for f in "$OUT/02_grid/$REC".*map "$OUT/02_grid/$REC".maps.* "$OUT/02_grid/$REC.pdbqt"; do ln -sf "$f" .; done

    "$PYTHONSH" "$UTIL/prepare_dpf4.py" -l "$L.pdbqt" -r "$REC.pdbqt" -o "$L.dpf" \
        -p ga_run="$GA_RUN" -p ga_pop_size="$GA_POP" -p ga_num_evals="$GA_EVALS" \
        -p ga_num_generations="$GA_GENS" -p rmstol="$RMSTOL" > prepare_dpf.log 2>&1 \
        || { tail prepare_dpf.log; die "prepare_dpf4 falhou: $L"; }

    # sementes fixas: o mesmo config reproduz o mesmo resultado
    sed -i "s/^seed .*/seed $SEED1 $SEED2                  # sementes fixas/" "$L.dpf"
    grep -q '^analysis' "$L.dpf" || echo "analysis" >> "$L.dpf"
    # parametros do GA antes de set_ga, e nada de do_local_only
    awk '/^set_ga/{s=NR} /^ga_num_evals/{e=NR} /^ga_num_generations/{g=NR} /^ga_pop_size/{p=NR} /^do_local_only/{bad=1}
         END{ if(!s || !e || !g || e>s || g>s || p>s || bad){exit 1} }' "$L.dpf" || die "DPF invalido: $L"

    log "[$L] autodock4 iniciado"
    "$AUTODOCK4" -p "$L.dpf" -l "$L.dlg" > "$L.run.log" 2>&1 || true
    grep -q "Successful Completion" "$L.dlg" && log "[$L] ok" || { log "[$L] FALHOU (ver $D/$L.dlg)"; return 1; }
}

step_dock() {
    local L pids=() names=() fail=0 i
    for L in "${LIGANDS[@]}"; do
        while [ "$(jobs -rp | wc -l)" -ge "$DOCK_JOBS" ]; do sleep 5; done
        ( dock_one "$L" ) &
        pids+=($!); names+=("$L")
    done
    for i in "${!pids[@]}"; do wait "${pids[$i]}" || { warn "docking falhou: ${names[$i]}"; fail=1; }; done
    [ $fail -eq 0 ] || die "docking com falhas"
}
