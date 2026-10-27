#!/usr/bin/env bash
# Funcoes comuns e valores padrao. Tudo que esta aqui pode ser sobrescrito no config.

log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
warn() { printf '[%s] AVISO: %s\n' "$(date '+%F %T')" "$*" >&2; }
die()  { printf '[%s] ERRO: %s\n' "$(date '+%F %T')" "$*" >&2; exit 1; }

need_file() { [ -s "$1" ] || die "arquivo ausente ou vazio: $1"; }
need_cmd()  { command -v "$1" >/dev/null 2>&1 || die "comando nao encontrado: $1"; }

# copia removendo CR (arquivos vindos do Windows)
copy_lf() { tr -d '\r' < "$1" > "$2"; }

lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }

load_config() {
    # ---- padroes ----
    PROJECT=project
    SYSTEM_TAG=sys
    OUT=./run
    LIGANDS=()
    declare -gA LIG_PDB=() LIG_ITP=() LIG_CHARGE=()

    MGLTOOLS=""
    AUTOGRID4=autogrid4
    AUTODOCK4=autodock4
    PYTHON=python3
    GMX=gmx
    PYMOL=pymol

    RECEPTOR=""
    RECEPTOR_PREP_FLAGS="-A checkhydrogens -U nphs_lps_waters_nonstdres"

    GRID_CENTER=""              # "x,y,z" em A
    GRID_CENTER_RESIDUES=""     # ou "B:467 B:397 A:62" (centroide dos residuos no receptor)
    GRID_CENTER_LIGAND=""       # ou um PDB de ligante de referencia (centroide)
    GRID_NPTS="60,60,60"
    GRID_SPACING=0.375

    GA_RUN=100
    GA_POP=150
    GA_EVALS=2500000
    GA_GENS=27000
    RMSTOL=2.0
    SEED1=20261026
    SEED2=4242
    DOCK_JOBS=4                 # ligantes em paralelo

    SELECT_DELTA=1.0            # kcal/mol: entre clusters ate DELTA do minimo global, o mais populoso
    REFERENCE_POSES=""          # opcional: pasta com <Ligante>.pdbqt para comparar RMSD
    LIG_SYMMETRY=""             # homodimero: "A:B" copia o ligante para o sitio simetrico na MD

    COFACTOR_PDB=""             # PDB com os cofatores nas coordenadas do receptor (opcional)
    COFACTORS=()                # "NOME:CADEIA:RESNUM", ex. "FAD:A:601"
    declare -gA COFACTOR_ITP=() COFACTOR_REF=()

    FF_DIR=""                   # ex. .../gromos54a7_atb.ff
    WATER=spc
    WATER_GRO=spc216.gro
    BOX_TYPE=dodecahedron
    BOX_D=1.2
    SALT_CONC=0.15              # mol/L de NaCl, alem da neutralizacao
    PDB2GMX_FLAGS="-ignh -merge no"
    PROTONATION=pdb2gmx         # pdb2gmx (padrao do GROMACS, ~pH 7) | propka (pKa no microambiente, no pH abaixo)
    PH=7.0
    PROPKA=propka3
    MDP_DIR=""                  # padrao: mdp/ do repositorio
    LOCAL_EM=yes
    LOCAL_NT=0                  # 0 = deixa o mdrun decidir

    PBS_QUEUE=par16
    PBS_NODES=1
    PBS_PPN=16
    PBS_MEM=60gb
    PBS_WALLTIME=240:00:00
    PBS_SAFE_HOURS=230
    HPC_MODULE="gromacs/2021.2-intel-2021.3.0"
    HPC_GMX=gmx_d
    HPC_MPIRUN="mpirun -np 16"
    PROD_NS=100
    MDPR_PS=1000
    FMAX_TARGET_EMPR=100
    FMAX_TARGET_EM=50
    FMAX_TARGET_CG=10
    FMAX_ABORT=10000
    MAXWARN_EM=1
    MAXWARN_MD=2

    # ---- config do usuario ----
    local cfgdir; cfgdir=$(dirname "$1")
    pushd "$cfgdir" >/dev/null
    # shellcheck disable=SC1090
    source "$1"
    popd >/dev/null

    # caminhos relativos ao config
    OUT=$(abspath_from "$cfgdir" "$OUT")
    RECEPTOR=$(abspath_from "$cfgdir" "$RECEPTOR")
    [ -n "$COFACTOR_PDB" ] && COFACTOR_PDB=$(abspath_from "$cfgdir" "$COFACTOR_PDB")
    [ -n "$FF_DIR" ] && FF_DIR=$(abspath_from "$cfgdir" "$FF_DIR")
    [ -n "$GRID_CENTER_LIGAND" ] && GRID_CENTER_LIGAND=$(abspath_from "$cfgdir" "$GRID_CENTER_LIGAND")
    [ -n "$REFERENCE_POSES" ] && REFERENCE_POSES=$(abspath_from "$cfgdir" "$REFERENCE_POSES")
    MDP_DIR=${MDP_DIR:-$HERE/mdp}
    MDP_DIR=$(abspath_from "$cfgdir" "$MDP_DIR")
    local k
    for k in "${!LIG_PDB[@]}"; do LIG_PDB[$k]=$(abspath_from "$cfgdir" "${LIG_PDB[$k]}"); done
    for k in "${!LIG_ITP[@]}"; do LIG_ITP[$k]=$(abspath_from "$cfgdir" "${LIG_ITP[$k]}"); done
    for k in "${!COFACTOR_ITP[@]}"; do COFACTOR_ITP[$k]=$(abspath_from "$cfgdir" "${COFACTOR_ITP[$k]}"); done
    for k in "${!COFACTOR_REF[@]}"; do COFACTOR_REF[$k]=$(abspath_from "$cfgdir" "${COFACTOR_REF[$k]}"); done

    [ ${#LIGANDS[@]} -gt 0 ] || die "config: LIGANDS vazio"
    REC=receptor
    PYTHONSH=${MGLTOOLS:+$MGLTOOLS/bin/pythonsh}
    UTIL=${MGLTOOLS:+$MGLTOOLS/MGLToolsPckgs/AutoDockTools/Utilities24}
    FF_NAME=$(basename "${FF_DIR:-none.ff}" .ff)
}

abspath_from() {
    # abspath_from <dir base> <caminho>
    local base=$1 p=$2
    [ -z "$p" ] && { printf ''; return; }
    case "$p" in
        /*) printf '%s' "$p" ;;
        *)  printf '%s' "$(cd "$base" && pwd)/$p" ;;
    esac
}

# prefixo dos arquivos de MD de um ligante: <tag>_<ligante em minusculas>
md_prefix() { printf '%s_%s' "$SYSTEM_TAG" "$(lower "$1")"; }

# primeira linha de [ moleculetype ] de um ITP
itp_moltype() { awk '/\[ *moleculetype *\]/{f=1;next} f && !/^[;#]/ && NF{print $1; exit}' "$1"; }

# ativa o GROMACS/Python se estiverem num ambiente conda indicado no config
activate_env() {
    [ -n "${CONDA_ENV:-}" ] || return 0
    local base="${CONDA_BASE:-}" c
    if [ -z "$base" ] && [ -n "${CONDA_EXE:-}" ]; then base=$(dirname "$(dirname "$CONDA_EXE")"); fi
    if [ -z "$base" ] && command -v conda >/dev/null 2>&1; then base=$(conda info --base 2>/dev/null); fi
    for c in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" /opt/conda; do
        [ -z "$base" ] && [ -f "$c/etc/profile.d/conda.sh" ] && base=$c
    done
    if [ -z "$base" ] || [ ! -f "$base/etc/profile.d/conda.sh" ]; then
        warn "conda nao encontrado (defina CONDA_BASE); seguindo com o PATH atual"
        return 0
    fi
    set +u
    # shellcheck disable=SC1091
    source "$base/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV" || warn "nao consegui ativar o ambiente $CONDA_ENV"
    set -u
}
