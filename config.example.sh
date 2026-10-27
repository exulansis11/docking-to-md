# docking-to-md -- modelo de configuracao. Copie, edite e rode:
#   ./pipeline.sh meu_config.sh all
# Caminhos relativos sao resolvidos a partir da pasta deste arquivo.
# Arrays associativos: use NOME[chave]=valor (nao "declare -A").
# Valores omitidos usam os padroes de scripts/lib.sh.

PROJECT=MeuAlvo
SYSTEM_TAG=alvo                     # pastas de MD: alvo_<ligante>
OUT=./run                           # saida
# CONDA_ENV=docking-to-md           # opcional: ativa este ambiente antes de gmx/pymol
# CONDA_BASE=$HOME/miniconda3

# ---- software ----
MGLTOOLS=/caminho/mgltools_x86_64Linux2_1.5.7   # com bin/pythonsh (rode o install.sh antes)
AUTOGRID4=autogrid4
AUTODOCK4=autodock4
PYTHON=python3                      # com numpy e networkx
GMX=gmx
PYMOL=pymol

# ---- receptor: PDB ou PDBQT, so a proteina (sem agua, ligantes e cofatores) ----
RECEPTOR=entrada/receptor.pdb

# ---- ligantes: PDB e ITP all-atom do ATB (mesmo molid e revisao) ----
LIGANDS=(LIG1 LIG2)
LIG_PDB[LIG1]=entrada/LIG1_allatom_optimised_geometry.pdb
LIG_ITP[LIG1]=entrada/LIG1_GROMACS_G54A7FF_allatom.itp
LIG_CHARGE[LIG1]=0                  # opcional: confere a carga de Gasteiger
LIG_PDB[LIG2]=entrada/LIG2_allatom_optimised_geometry.pdb
LIG_ITP[LIG2]=entrada/LIG2_GROMACS_G54A7FF_allatom.itp

# ---- caixa: escolha UMA forma de definir o centro ----
GRID_CENTER="10.0,20.0,30.0"                 # x,y,z em A
# GRID_CENTER_RESIDUES="A:145 A:41 A:166"    # centroide dos residuos (cadeia:numero)
# GRID_CENTER_LIGAND=entrada/ligante_cristal.pdb
GRID_NPTS="60,60,60"
GRID_SPACING=0.375

# ---- Lamarckian GA ----
GA_RUN=100
GA_POP=150
GA_EVALS=2500000
GA_GENS=27000
RMSTOL=2.0
SEED1=20261026
SEED2=4242
DOCK_JOBS=4                         # ligantes em paralelo (1 nucleo cada)
SELECT_DELTA=1.0                    # kcal/mol (ver README)

# ---- homodimero (opcional): copia o ligante para o sitio simetrico na MD ----
# LIG_SYMMETRY="A:B"

# ---- cofatores para a MD (opcional): coordenadas do cristal alinhadas ao receptor ----
# COFACTOR_PDB=entrada/cofatores.pdb
# COFACTORS=("FAD:A:601")           # NOME:CADEIA:RESNUM (uma entrada por copia)
# COFACTOR_ITP[FAD]=entrada/FAD_GROMACS_G54A7FF_allatom.itp
# COFACTOR_REF[FAD]=entrada/FAD_allatom_optimised_geometry.pdb

# ---- MD ----
FF_DIR=/caminho/gromos54a7_atb.ff   # versao ATUAL do ATB (atb_gromacs/5/gromos54a7_atb.ff.tar.gz)
WATER=spc
BOX_TYPE=dodecahedron
BOX_D=1.2
SALT_CONC=0.15                      # mol/L de NaCl, alem da neutralizacao
PDB2GMX_FLAGS="-ignh -merge no"
PROTONATION=pdb2gmx                 # ou propka: estados de ASP/GLU/LYS/ARG/HIS pelo pKa no microambiente
# PH=7.0                            # pH usado com PROTONATION=propka
# PROPKA=propka3
# MDP_DIR=meus_mdps                 # padrao: mdp/ do repositorio
LOCAL_EM=yes

# ---- HPC (PBS) ----
PBS_QUEUE=par16
PBS_NODES=1
PBS_PPN=16
PBS_MEM=60gb
PBS_WALLTIME=240:00:00
PBS_SAFE_HOURS=230                  # o mdrun para antes do limite e o PBS pode ser ressubmetido
HPC_MODULE="gromacs/2021.2-intel-2021.3.0"
HPC_GMX=gmx_d
HPC_MPIRUN="mpirun -np 16"
PROD_NS=100                         # producao (ns), depois de 1 ns de MDPR
MDPR_PS=1000                        # = nsteps*dt do mdpr298.mdp
