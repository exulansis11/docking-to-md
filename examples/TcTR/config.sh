# docking-to-md -- exemplo: tripanotiona redutase de T. cruzi (dimero, PDB 1BZL) + 4 ligantes (LIG1-LIG4)
# Caminhos relativos sao resolvidos a partir da pasta deste arquivo.
# Arrays associativos: use atribuicoes NOME[chave]=valor (nao "declare -A", que criaria variavel local).

PROJECT=TcTR
SYSTEM_TAG=1bzl                     # prefixo das pastas de MD: 1bzl_<ligante>
OUT=./run                           # saida (de preferencia fora de um disco lento/instavel)
# CONDA_ENV=docking-to-md          # ambiente com gmx, numpy, networkx, pymol (opcional)

# ---- software ----
MGLTOOLS=/caminho/mgltools_x86_64Linux2_1.5.7
AUTOGRID4=autogrid4
AUTODOCK4=autodock4
PYTHON=python3
GMX=gmx
PYMOL=pymol

# Entradas (nao incluidas no repositorio) organizadas em entrada/ ao lado deste arquivo.

# ---- receptor: TcTR (dimero, PDB 1BZL), H polares no padrao ADT, sem cofatores ----
# (o FAD fica a ~9 A do sitio)
RECEPTOR=entrada/receptor/1bzl_receptor.pdbqt

# ---- ligantes: PDB e ITP all-atom do ATB (mesmo molid/revisao), carga esperada ----
# Os nomes reais foram substituidos por LIG1-LIG4 (compostos ainda em estudo).
LIGANDS=(LIG1 LIG2 LIG3 LIG4)
for L in "${LIGANDS[@]}"; do
    LIG_PDB[$L]=entrada/ligantes/${L}_allatom_optimised_geometry.pdb
    LIG_ITP[$L]=entrada/ligantes/${L}_GROMACS_G54A7FF_allatom.itp
done
LIG_CHARGE[LIG1]=1
LIG_CHARGE[LIG2]=1
LIG_CHARGE[LIG3]=0
LIG_CHARGE[LIG4]=0

# ---- caixa: sitio da tripanotiona, na interface do dimero ----
GRID_CENTER="25.183,9.590,-0.750"
GRID_NPTS="84,80,72"
GRID_SPACING=0.375
# alternativas: GRID_CENTER_RESIDUES="B:467 B:397 B:463 A:62"  ou  GRID_CENTER_LIGAND=ligante_cristal.pdb

# ---- Lamarckian GA ----
GA_RUN=100
GA_POP=150
GA_EVALS=2500000
GA_GENS=27000
RMSTOL=2.0
SEED1=20261026
SEED2=4242
DOCK_JOBS=4
SELECT_DELTA=1.0

# ---- homodimero: a TcTR tem 2 sitios da tripanotiona equivalentes (interface A/B).
# O docking e feito num deles; na MD o ligante e copiado para o outro pela simetria do dimero.
LIG_SYMMETRY="A:B"

# ---- cofatores (MD): FAD oxidado -2, ATB molid 2053771, nos dois sitios ----
COFACTOR_PDB=entrada/cofatores/FAD_cristal_alinhado.pdb
COFACTORS=("FAD:A:601" "FAD:B:602")
COFACTOR_ITP[FAD]=entrada/cofatores/9Y2Q_GROMACS_G54A7FF_allatom.itp
COFACTOR_REF[FAD]=entrada/cofatores/9Y2Q_allatom_optimised_geometry.pdb

# ---- MD (protocolo TcTR do laboratorio) ----
# pH 7,5 e 298 K (dados de TSA); NaCl 500 mM; protonacao pelo PROPKA 3.1; GROMOS 54A7 do ATB.
# A temperatura (298 K) esta nos MDPs (mdpr298.mdp e mdnvt298.mdp: ref-t e gen-temp).
FF_DIR=entrada/gromos54a7_atb.ff            # GROMOS 54A7 do ATB (versao atual!)
WATER=spc
BOX_TYPE=dodecahedron
BOX_D=1.2
SALT_CONC=0.5                               # mol/L de NaCl, alem da neutralizacao
PROTONATION=propka                          # pKa no microambiente (PROPKA) ...
PH=7.5                                      # ... no pH do experimento
PROPKA=propka3
LOCAL_EM=yes
LOCAL_NT=12

# ---- HPC (CENAPAD, PBS) ----
PBS_QUEUE=par16
PBS_PPN=16
PBS_MEM=60gb
PBS_WALLTIME=240:00:00
PBS_SAFE_HOURS=230
HPC_MODULE="gromacs/2021.2-intel-2021.3.0"
HPC_GMX=gmx_d
HPC_MPIRUN="mpirun -np 16"
PROD_NS=100                         # producao (ns), depois de 1 ns de MDPR
MDPR_PS=1000                        # = nsteps*dt do mdpr298.mdp
