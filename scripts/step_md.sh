#!/usr/bin/env bash
# md: monta receptor + cofatores + ligante (pose escolhida; 2 copias num homodimero) em GROMACS, solvata, neutraliza, grompp e EM local.
#
# Ligantes e cofatores usam ITPs do ATB (all-atom recomendado; o united-atom e "deprecated" no ATB V4).
# Atomos pesados casados por isomorfismo de grafo (funciona com nomes diferentes, ex. FAD do cristal x ATB).
# H polares do ligante vem da pose; os demais do PDB do ATB (geometria de referencia).

run_pdb2gmx() {
    # run_pdb2gmx <pasta> <pdb so com pesados>: protein.gro, protein.top, itps e molecules.txt na pasta
    local P=$1 heavy=$2
    mkdir -p "$P"; cd "$P"
    ln -sfn "$FF_DIR" "$FF_NAME.ff"
    [ "$heavy" = "$P/protein_heavy.pdb" ] || cp "$heavy" protein_heavy.pdb
    # shellcheck disable=SC2086
    "$GMX" pdb2gmx -f protein_heavy.pdb -o protein.gro -p protein.top -i posre.itp \
        -ff "$FF_NAME" -water "$WATER" $PDB2GMX_FLAGS > pdb2gmx.log 2>&1 \
        || { tail -30 pdb2gmx.log; die "pdb2gmx falhou (ver $P/pdb2gmx.log)"; }
    # protonacao pelo PROPKA no pH do config: renomeia ASP/GLU/LYS/ARG/HIS e refaz o pdb2gmx
    if [ "$PROTONATION" = propka ]; then
        need_cmd "$PROPKA"
        mv pdb2gmx.log pdb2gmx_auto.log
        "$PROPKA" protein_heavy.pdb > propka.log 2>&1 || { tail propka.log; die "PROPKA falhou (ver $P/propka.log)"; }
        "$PYTHON" "$HERE/scripts/propka_protonate.py" protein_heavy.pdb protein_heavy.pka "$PH" \
            pdb2gmx_auto.log protein_heavy_pka.pdb protonation.tsv | sed 's/^/    /'
        # shellcheck disable=SC2086
        "$GMX" pdb2gmx -f protein_heavy_pka.pdb -o protein.gro -p protein.top -i posre.itp \
            -ff "$FF_NAME" -water "$WATER" $PDB2GMX_FLAGS > pdb2gmx.log 2>&1 \
            || { tail -30 pdb2gmx.log; die "pdb2gmx (protonacao PROPKA) falhou (ver $P/pdb2gmx.log)"; }
    fi
    grep -E 'Linking|Total charge' pdb2gmx.log | sed 's/^/    /'
    # cadeia unica: o pdb2gmx poe o moleculetype dentro do .top -> extrair para protein.itp
    if ! grep -qE '#include "[^"]*Protein[^"]*\.itp"' protein.top; then
        awk '/\[ *moleculetype *\]/{f=1} /; Include water topology/{f=0} f' protein.top > protein.itp
    fi
    # nomes dos moleculetypes da proteina, na ordem do protein.top
    awk '/\[ *molecules *\]/{f=1;next} f && !/^;/ && NF==2 && $1!="SOL"{print $1}' protein.top > molecules.txt
    rm -f \#*
    cd - >/dev/null
}

md_protein() {
    local P="$OUT/05_md/protein"
    mkdir -p "$P"
    "$PYTHON" "$HERE/scripts/protein_heavy.py" "$RECEPTOR" "$P/protein_heavy.pdb"
    [ -n "$LIG_SYMMETRY" ] || run_pdb2gmx "$P" "$P/protein_heavy.pdb"
}

merge_gro() {
    # merge_gro <saida> <titulo> <gro...>
    "$PYTHON" - "$@" <<'PY'
import sys
out, title, parts = sys.argv[1], sys.argv[2], sys.argv[3:]
body, box = [], None
for p in parts:
    lines = open(p).read().splitlines()
    n = int(lines[1])
    body += lines[2:2 + n]
    box = box or lines[2 + n]
open(out, "w").write(f"{title}\n{len(body)}\n" + "\n".join(body) + "\n" + box + "\n")
PY
}

md_one() {
    local L=$1
    local P="$OUT/05_md/protein" S="$OUT/05_md/$L"
    rm -rf "$S"; mkdir -p "$S"; cd "$S"
    ln -sfn "$FF_DIR" "$FF_NAME.ff"
    local ncof=${#COFACTORS[@]}

    # ligante na pose escolhida
    "$PYTHON" "$HERE/scripts/place_molecule.py" "${LIG_ITP[$L]}" "$OUT/04_analysis/${L}_best.pdbqt" lig \
        --ref "${LIG_PDB[$L]}" --resid "$((ncof + 1))" | tee -a build.log
    mv lig.itp "$L.itp"; mv lig_posre.itp "${L}_posre.itp"; sed -i "s/lig_posre.itp/${L}_posre.itp/" "$L.itp"

    # homodimero: 2a copia no sitio simetrico; cadeias laterais que colidem com ela recebem a imagem
    # simetrica da outra cadeia (o 2o sitio fica igual ao do docking); pdb2gmx proprio deste sistema
    if [ -n "$LIG_SYMMETRY" ]; then
        "$PYTHON" "$HERE/scripts/symmetric_copy.py" "$P/protein_heavy.pdb" "$LIG_SYMMETRY" lig.gro lig2.gro \
            --resid "$((ncof + 2))" --fix-out protein_heavy_sym.pdb | tee -a build.log
        run_pdb2gmx "$S/protein" "$S/protein_heavy_sym.pdb"
        P="$S/protein"
    fi

    cp "$P"/*.itp .
    local entry name chain resnum moltype i=0 gros=("$P/protein.gro") includes=() molecules=()
    while read -r m; do molecules+=("$m 1"); done < "$P/molecules.txt"
    # itps da proteina (com posre) incluidos pelo pdb2gmx
    while read -r inc; do includes+=("$inc"); done < <(grep -oE '#include "[^"]*Protein[^"]*\.itp"' "$P/protein.top" | sed 's/#include "\(.*\)"/\1/')
    [ ${#includes[@]} -gt 0 ] || includes=("protein.itp")

    # cofatores
    declare -A seen=()
    for entry in "${COFACTORS[@]}"; do
        IFS=: read -r name chain resnum <<<"$entry"
        i=$((i + 1))
        local args=("${COFACTOR_ITP[$name]}" "$COFACTOR_PDB" "cof$i" --select "$chain:$resnum" --resid "$i")
        [ -n "${COFACTOR_REF[$name]:-}" ] && args+=(--ref "${COFACTOR_REF[$name]}")
        if [ -n "${seen[$name]:-}" ]; then
            "$PYTHON" "$HERE/scripts/place_molecule.py" "${args[@]}" --gro-only | tee -a build.log
        else
            "$PYTHON" "$HERE/scripts/place_molecule.py" "${args[@]}" | tee -a build.log
            mv "cof$i.itp" "$name.itp"; mv "cof${i}_posre.itp" "${name}_posre.itp"
            sed -i "s/cof${i}_posre.itp/${name}_posre.itp/" "$name.itp"
            includes+=("$name.itp"); seen[$name]=1
        fi
        gros+=("cof$i.gro")
        molecules+=("$(itp_moltype "$name.itp") 1")
    done

    includes+=("$L.itp"); gros+=("lig.gro")
    moltype=$(itp_moltype "$L.itp")
    molecules+=("$moltype 1")
    if [ -n "$LIG_SYMMETRY" ]; then
        # auditoria final da copia contra proteina (com H do campo de forca) + cofatores + 1o ligante
        merge_gro check.gro "checagem" "${gros[@]}"
        "$PYTHON" "$HERE/scripts/symmetric_copy.py" "$S/protein_heavy_sym.pdb" "$LIG_SYMMETRY" lig.gro lig2.gro \
            --resid "$((ncof + 2))" --system check.gro | grep -E 'menor distancia|correspondencia' | tee -a build.log \
            || die "copia simetrica com choque: $L"
        rm -f check.gro
        gros+=("lig2.gro")
        molecules+=("$moltype 1")
    fi

    # coordenadas
    merge_gro complex.gro "complexo + $L" "${gros[@]}"

    # topologia
    {
        echo "; $PROJECT: receptor + cofatores + $L ($moltype)"
        echo "#include \"$FF_NAME.ff/forcefield.itp\""
        echo
        for inc in "${includes[@]}"; do echo "#include \"$inc\""; done
        echo
        echo "#include \"$FF_NAME.ff/$WATER.itp\""
        echo "#ifdef POSRES_WATER"
        echo "[ position_restraints ]"
        echo "   1    1       1000       1000       1000"
        echo "#endif"
        echo "#include \"$FF_NAME.ff/ions.itp\""
        echo
        echo "[ system ]"
        echo "$PROJECT + $L"
        echo
        echo "[ molecules ]"
        printf '%s\n' "${molecules[@]}"
    } > topol.top

    "$GMX" editconf -f complex.gro -o box.gro -bt "$BOX_TYPE" -d "$BOX_D" -c > editconf.log 2>&1
    "$GMX" solvate -cp box.gro -cs "$WATER_GRO" -o solv.gro -p topol.top > solvate.log 2>&1
    # avisos esperados aqui: GROMOS + PME e carga ainda nao neutralizada
    printf 'integrator=steep\nnsteps=0\ncutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\n' > ions.mdp
    "$GMX" grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 > grompp_ions.log 2>&1 \
        || { grep -B2 -A10 'ERROR' grompp_ions.log; die "grompp (ions) falhou: $L"; }
    echo SOL | "$GMX" genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral -conc "$SALT_CONC" > genion.log 2>&1
    printf 'q\n' | "$GMX" make_ndx -f system.gro -o index.ndx > /dev/null 2>&1
    "$GMX" editconf -f complex.gro -o complex.pdb > /dev/null 2>&1   # complexo sem agua, para visualizar

    "$GMX" grompp -f "$MDP_DIR/empr.mdp" -c system.gro -r system.gro -p topol.top -o empr.tpr -maxwarn "$MAXWARN_EM" \
        > grompp_empr.log 2>&1 || { grep -B2 -A10 'ERROR' grompp_empr.log; die "grompp (empr) falhou: $L"; }
    log "[$L] $(sed -n 2p system.gro | tr -d ' ') atomos; $(grep -E '^(SOL|NA|CL) ' topol.top | tr -s ' ' | tr '\n' ' ')"

    if [ "$LOCAL_EM" = yes ]; then
        local nt=(); [ "$LOCAL_NT" -gt 0 ] && nt=(-nt "$LOCAL_NT")
        "$GMX" mdrun -deffnm empr "${nt[@]}" > mdrun_empr.log 2>&1
        "$GMX" grompp -f "$MDP_DIR/em.mdp" -c empr.gro -p topol.top -o em.tpr -maxwarn "$MAXWARN_EM" > grompp_em.log 2>&1
        "$GMX" mdrun -deffnm em "${nt[@]}" > mdrun_em.log 2>&1
        log "[$L] EM local: $(grep -E 'Maximum force' em.log | tail -1 | sed 's/  */ /g')"
        rm -f step*.pdb
    fi
    rm -f \#* mdout.mdp
    cd - >/dev/null
}

step_md() {
    activate_env
    [ -s "$FF_DIR/forcefield.itp" ] || die "FF_DIR invalido: $FF_DIR"
    need_cmd "$GMX"
    local L
    for L in "${LIGANDS[@]}"; do need_file "${LIG_ITP[$L]}"; need_file "$OUT/04_analysis/${L}_best.pdbqt"; done
    md_protein
    for L in "${LIGANDS[@]}"; do md_one "$L"; done
}
