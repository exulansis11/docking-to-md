#!/usr/bin/env bash
# hpc: pacote para o cluster (PBS). Uma pasta e um PBS por ligante, campo de forca um nivel acima.
# Saida: <OUT>/06_hpc/<PROJECT>_md/ e <PROJECT>_md.tar.gz

step_hpc() {
    local H="$OUT/06_hpc" NAME="${PROJECT}_md"
    local PKG="$H/$NAME"
    rm -rf "$PKG"; mkdir -p "$PKG"
    cp -rL "$FF_DIR" "$PKG/$FF_NAME.ff"

    local L P S D job moltype prod_ps systems=()
    prod_ps=$((PROD_NS * 1000))
    for L in "${LIGANDS[@]}"; do
        P=$(md_prefix "$L"); S="$OUT/05_md/$L"; D="$PKG/$P"
        need_file "$S/system.gro"
        mkdir -p "$D"; systems+=("$P")
        cp "$S/system.gro" "$D/${P}_04_solvated.gro"
        cp "$S/index.ndx" "$D/${P}_04_index.ndx"
        cp "$S"/*.itp "$D/"
        sed "s|#include \"$FF_NAME.ff/|#include \"../$FF_NAME.ff/|" "$S/topol.top" > "$D/${P}_04_topology.top"
        local m
        for m in empr em cg mdpr298 mdnvt298; do
            # MDPs as vezes vem em Latin-1: converter para UTF-8
            iconv -f utf-8 -t utf-8 "$MDP_DIR/$m.mdp" >/dev/null 2>&1 && cp "$MDP_DIR/$m.mdp" "$D/" \
                || iconv -f latin1 -t utf-8 "$MDP_DIR/$m.mdp" > "$D/$m.mdp"
        done
        moltype=$(itp_moltype "$S/$L.itp")
        job=$(printf '%s' "$L" | tr '[:lower:]' '[:upper:]' | cut -c1-9)-${PROD_NS}NS
        sed -e "s|@PREFIX@|$P|g" -e "s|@JOBNAME@|${job:0:15}|g" -e "s|@LIG@|$L|g" -e "s|@MOLTYPE@|$moltype|g" \
            -e "s|@PROJECT@|$PROJECT|g" -e "s|@QUEUE@|$PBS_QUEUE|g" -e "s|@NODES@|$PBS_NODES|g" -e "s|@PPN@|$PBS_PPN|g" \
            -e "s|@MEM@|$PBS_MEM|g" -e "s|@WALLTIME@|$PBS_WALLTIME|g" -e "s|@MODULE@|$HPC_MODULE|g" \
            -e "s|@GMX@|$HPC_GMX|g" -e "s|@MPIRUN@|$HPC_MPIRUN|g" -e "s|@FF_NAME@|$FF_NAME|g" \
            -e "s|@MAXWARN_EM@|$MAXWARN_EM|g" -e "s|@MAXWARN_MD@|$MAXWARN_MD|g" \
            -e "s|@FMAX_EMPR@|$FMAX_TARGET_EMPR|g" -e "s|@FMAX_EM@|$FMAX_TARGET_EM|g" -e "s|@FMAX_CG@|$FMAX_TARGET_CG|g" \
            -e "s|@FMAX_ABORT@|$FMAX_ABORT|g" -e "s|@PROD_NS@|$PROD_NS|g" -e "s|@PROD_PS@|$prod_ps|g" \
            -e "s|@MDPR_PS@|$MDPR_PS|g" -e "s|@SAFE_HOURS@|$PBS_SAFE_HOURS|g" \
            "$HERE/templates/run_md.pbs.in" > "$D/${P}_run_md.pbs"
        chmod +x "$D/${P}_run_md.pbs"
        bash -n "$D/${P}_run_md.pbs" || die "PBS com erro de sintaxe: $P"
    done

    sed "s|@SYSTEMS@|${systems[*]}|" "$HERE/templates/submit_all.sh.in" > "$PKG/submit_all.sh"
    chmod +x "$PKG/submit_all.sh"
    sed -e "s|@PROJECT@|$PROJECT|g" -e "s|@SYSTEMS@|${systems[*]}|g" -e "s|@PROD_NS@|$PROD_NS|g" \
        -e "s|@FF_NAME@|$FF_NAME|g" -e "s|@NAME@|$NAME|g" "$HERE/templates/README_HPC.txt.in" > "$PKG/README_HPC.txt"

    ( cd "$PKG" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS )
    ( cd "$H" && tar czf "$NAME.tar.gz" "$NAME" && sha256sum "$NAME.tar.gz" > "$NAME.tar.gz.sha256" )
    log "pacote: $H/$NAME.tar.gz ($(du -h "$H/$NAME.tar.gz" | cut -f1))"
}
