# docking-to-md

Do **docking com AutoDock4** até o **pacote de dinâmica molecular pronto para o cluster** (GROMACS + PBS), num único comando: prepara receptor e ligantes, roda o AutoGrid e o AutoDock, escolhe a melhor pose de cada ligante, monta proteína + cofatores + ligante com topologias do ATB (GROMOS 54A7), solvata, neutraliza, minimiza localmente e gera um PBS por ligante, empacotado em `tar.gz`, só para enviar e disparar.

*English summary at the end.*

```
receptor + ligantes (ATB)
   │ prepare   MGLTools: H polares, Gasteiger, torsões → PDBQT
   │ grid      uma caixa para todos os ligantes → AutoGrid
   │ dock      Lamarckian GA, sementes fixas, ligantes em paralelo → AutoDock4
   │ analyze   clusters, escolha da pose, ΔG, Ki, contatos → tabela
   │ md        pdb2gmx + cofatores + ligante na pose (+ cópia simétrica em homodímeros) → caixa, água, íons, grompp, EM
   │ hpc       um PBS por ligante (empr → em → cg → MDPR → produção) → tar.gz
   ▼ figures   PyMOL: visão geral e sítio de cada ligante
```

## Instalação

```bash
git clone https://github.com/exulansis11/docking-to-md.git
cd docking-to-md
conda env create -f environment.yml      # autodock4, autogrid4, gromacs, numpy, networkx, pymol
conda activate docking-to-md
```

**MGLTools 1.5.7** (os `prepare_*4.py`, em Python 2) não está no conda. Baixe em
<https://ccsb.scripps.edu/mgltools/downloads/>, extraia e rode o `install.sh` dentro da pasta; depois aponte `MGLTOOLS=` no config para ela. Teste: `bin/pythonsh -c "import MolKit"`.

**Arquivos do ATB** (<https://atb.uq.edu.au>, requer login para baixar):

- de cada ligante (e de cada cofator): o **ITP all-atom** (`*_GROMACS_G54A7FF_allatom.itp`) e o **PDB all-atom otimizado** (`*_allatom_optimised_geometry.pdb`), do mesmo molid e da mesma revisão;
- o campo de força **GROMOS 54A7 do ATB, versão atual**: página *Forcefield Files* → "Modified GROMOS Force Field files … for Gromacs 4.5.x_5.x.x" (`atb_gromacs/5/gromos54a7_atb.ff.tar.gz`).

## Uso

```bash
cp config.example.sh meu_config.sh       # edite caminhos, ligantes, caixa, cluster
./pipeline.sh meu_config.sh check        # confere software e entradas
./pipeline.sh meu_config.sh all          # tudo
./pipeline.sh meu_config.sh dock analyze # ou só algumas etapas
```

`examples/TcTR/config.sh` é um exemplo completo: tripanotiona redutase de *T. cruzi* (dímero, PDB 1BZL) com 2 FAD e 4 ligantes (LIG1–LIG4). Os arquivos de entrada não estão no repositório; o config mostra como organizá-los em `entrada/`.

### Saídas (`OUT=` no config)

```
01_prepared/   receptor.pdbqt, <L>.pdbqt, prepare_audit.tsv (átomos, carga de Gasteiger, TORSDOF, tipos)
02_grid/       receptor.gpf, mapas, receptor.glg
03_docking/<L>/  <L>.dpf, <L>.dlg
04_analysis/   docking_summary.tsv, <L>_clusters.tsv, <L>_best.pdbqt, <L>_cluster_reps.pdbqt
05_md/protein/ pdb2gmx (uma vez, compartilhado)
05_md/<L>/     topol.top, system.gro (solvatado e neutro), index.ndx, complex.pdb, *.itp, empr/em (EM local)
06_hpc/        <PROJECT>_md/ e <PROJECT>_md.tar.gz (+ .sha256)
07_figures/    00_visao_geral.png, NN_<L>_sitio.png
pipeline.log
```

### No cluster

```bash
tar xzf <PROJECT>_md.tar.gz && cd <PROJECT>_md
bash submit_all.sh            # qsub de todos os sistemas
```

Cada `<sistema>_run_md.pbs` roda empr → em → cg → MDPR (POSRES) → produção, com QC de força máxima, checkpoint a cada 15 min e `-maxh` antes do limite do PBS. Se o job acabar antes do fim da produção, **ressubmeta o mesmo PBS**: etapas prontas são puladas e a produção continua do checkpoint. O `grompp` roda no próprio cluster (um `.tpr` de uma versão mais nova do GROMACS não abre numa mais antiga).

## Como cada decisão foi tomada

**Escolha da pose.** Entre os clusters (RMSD 2 Å, do próprio AutoDock) cujo menor ΔG está a até `SELECT_DELTA` (1,0 kcal/mol) do mínimo global, fica o **mais populoso**; empate → menor energia. A pose é a de menor energia desse cluster. O erro padrão da função do AD4 é ~2 kcal/mol, então um mínimo global isolado (cluster de 1 pose) não vence um cluster bem povoado quase tão bom. O Ki do AD4 serve para ordenar ligantes, não como afinidade.

**Reprodutibilidade.** `seed` fixo no DPF: o mesmo config gera o mesmo resultado.

**Pose → topologia do ATB.** Os átomos pesados da pose são casados com o ITP por **isomorfismo de grafo** (elemento + ligações), então funciona mesmo quando os nomes diferem (ex.: FAD do cristal × FAD do ATB). Os H polares vêm da pose; os apolares, do PDB do ATB por ajuste local (Kabsch). A ordem e os nomes saem iguais aos do `[ atoms ]` do ITP.

**Proteína.** O receptor sai sem H e o `pdb2gmx -ignh` recoloca os H do campo de força. Protonação das His por rede de H-bond (padrão do `pdb2gmx`); para escolher, acrescente `-his` em `PDB2GMX_FLAGS`.

**Homodímeros (`LIG_SYMMETRY="A:B"`).** Quando o alvo tem dois sítios equivalentes (ex.: na interface de um homodímero), o docking é feito num deles e, na MD, o ligante é **copiado para o outro** pela operação de simetria do dímero: ajuste dos Cα comuns [A; B] → [B; A], que num dímero simétrico é uma rotação de ~180° (o pipeline confere: RMSD < 2 Å e ângulo perto de 180°). Se algum resíduo do 2º sítio colidir com a cópia (< 2,2 Å, tipicamente uma Lys/Arg com rotâmero diferente no cristal), a cadeia lateral dele recebe a imagem simétrica da cadeia lateral equivalente do 1º sítio, para que a cópia fique no mesmo ambiente em que a pose foi obtida. Tudo é registrado em `05_md/<L>/build.log`: simetria, resíduos alterados, correspondência dos resíduos de contato entre os dois sítios e a menor distância da cópia ao resto do sistema. A alternativa mais conservadora é fazer um novo docking no 2º sítio.

**Restrição de posição.** `-DPOSRES` restringe os átomos pesados da proteína, dos cofatores e do ligante (1000 kJ mol⁻¹ nm⁻²).

## Armadilhas que o pipeline já trata (e que custaram tempo)

- **`prepare_dpf4.py -L` não é "Lamarckian"**: é só busca local. O padrão já é o LGA. O pipeline também confere que `ga_num_evals` e afins vêm antes do `set_ga`.
- **Nomes de H no receptor.** O MolKit só liga H no padrão PDB v2/ADT (`1HH1`, `2HH2`…). PDBs com H renomeados (`HH11`…) ficam com H soltos e cargas de Gasteiger erradas, sem erro aparente. Confira a carga total em `prepare_audit.tsv`.
- **CRLF.** Arquivos vindos do Windows são convertidos antes do MolKit.
- **Versão do campo de força do ATB.** ITPs recentes usam tipos (`COpt`, `HOpt`, `SOpt`, `C_H3`, `HCH3`…) que não existem em cópias antigas do `gromos54a7_atb.ff`, e tipos antigos (`OEOpt`) somem nas novas. Use ITPs e campo de força da mesma época.
- **United-atom é "deprecated" no ATB V4.** O pipeline usa ITPs all-atom para ligantes e cofatores.
- **QC com all-atom.** As ligações X–H quárticas do GROMOS são muito rígidas: minimizações sem constraints param com Fmax ~10³ nelas (ex.: N⁺–H em ponte salina). Por isso o PBS mostra o alvo clássico (EMPR < 100, EM < 50, CG < 10) como PASS/WARN e só interrompe acima de `FMAX_ABORT` (10⁴, choque de verdade). Na MD as ligações ficam restritas.
- **`-maxwarn` auditado**: 1 nas minimizações (aviso GROMOS/integrador), 2 na MD (+ barostato Berendsen dos MDPs incluídos). O PBS imprime quantos avisos cada `grompp` gerou.
- **MDPs em Latin-1** são convertidos para UTF-8 no pacote.

## Validação

O exemplo da TcTR, rodado do zero com `./pipeline.sh examples/TcTR/config.sh all` (ligantes anonimizados como LIG1–LIG4, compostos ainda em estudo; os arquivos de entrada não estão no repositório):

| Ligante | ΔG (kcal/mol) | Ki estimado | Cluster | RMSD vs. execução anterior | Átomos (solvatado, 2 cópias) |
|---|---|---|---|---|---|
| LIG3 | −9,22 | 174 nM | 22/100 | 0,00 Å | 146 930 |
| LIG2 | −8,60 | 496 nM | 60/100 | 0,00 Å | 146 935 |
| LIG1 | −7,75 | 2,10 µM | 32/100 | 0,00 Å | 146 955 |
| LIG4 | −7,01 | 7,30 µM | 24/100 | 0,00 Å | 146 914 |

- 18 min no total numa máquina de 12 núcleos (docking ~11 min com os 4 ligantes em paralelo; montagem + EM local ~5 min; figuras ~2 min).
- Com as mesmas sementes, as poses são idênticas às de uma execução independente feita etapa por etapa (RMSD 0,00 Å). Aplicada aos DLGs de uma execução anterior (sem sementes fixas), a regra de escolha reproduz as poses escolhidas manualmente.
- O PBS gerado foi executado localmente com etapas encurtadas: 05 → 09 completo e, na ressubmissão, todas as etapas prontas puladas.

## Limitações

- Receptor rígido; o docking usa uma caixa por execução (em homodímeros, o 2º sítio recebe a cópia simétrica).
- Cofatores não entram no docking, só na MD (coloque-os no receptor se quiser que influenciem o docking; o `prepare_receptor4` remove resíduos não padrão por padrão, ver `RECEPTOR_PREP_FLAGS`).
- O campo de força é GROMOS 54A7/ATB; outros exigiriam trocar a etapa `md`.
- Os MDPs em `mdp/` são os de um laboratório (SPC, PME 1,0 nm, v-rescale, Berendsen). Revise antes de produção longa: Berendsen não gera o ensemble correto (C-rescale é o padrão atual), e o `mdnvt298.mdp` grava a trajetória de todos os átomos a cada 10 ps.
- O PBS supõe `module load` e `mpirun`; para SLURM, adapte `templates/run_md.pbs.in`.

## Citação

Se usar, cite as ferramentas: AutoDock4 (Morris et al., *J. Comput. Chem.* 2009, 30, 2785), ATB (Malde et al., *JCTC* 2011, 7, 4026; Stroet et al., *JCTC* 2018, 14, 5834), GROMOS 54A7 (Schmid et al., *Eur. Biophys. J.* 2011, 40, 843) e GROMACS (Abraham et al., *SoftwareX* 2015, 1–2, 19).

## Licença

MIT.

---

## English summary

**docking-to-md** runs the whole path from **AutoDock4 docking** to a **ready-to-submit GROMACS MD package for a PBS cluster** with one config file and one command: MGLTools preparation (polar H, Gasteiger charges, torsions), AutoGrid + Lamarckian GA AutoDock4 with fixed seeds, cluster-based pose selection (most populated cluster within 1 kcal/mol of the global minimum), mapping of the pose onto ATB all-atom topologies by graph isomorphism (works with renamed atoms, e.g. crystal cofactors), `pdb2gmx` with GROMOS 54A7, cofactors, solvation, ions, `grompp` and local EM, and one restartable PBS script per ligand (EM → CG → restrained MD → production with checkpoint/`-maxh`), packed as `tar.gz`. It also guards against several silent pitfalls (the `-L` flag of `prepare_dpf4.py`, unbonded renamed hydrogens in MolKit, mismatched ATB force-field versions, stiff X–H bonds in all-atom minimization QC).

```bash
conda env create -f environment.yml && conda activate docking-to-md
./pipeline.sh config.example.sh all
```
