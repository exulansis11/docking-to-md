#!/usr/bin/env bash
# docking-to-md: do docking com AutoDock4 ate o pacote de dinamica molecular para HPC (GROMACS + PBS).
#
# Uso:
#   ./pipeline.sh <config.sh> [etapa ...]
#
# Etapas (na ordem; "all" roda todas):
#   check     confere software e arquivos de entrada
#   prepare   receptor e ligantes -> PDBQT (MGLTools: H polares, Gasteiger, torsoes)
#   grid      GPF + AutoGrid (uma caixa para todos os ligantes)
#   dock      DPF (Lamarckian GA) + AutoDock4, ligantes em paralelo
#   analyze   clusters, escolha da pose, tabela (dG, Ki, tamanho do cluster), contatos
#   md        pdb2gmx, cofatores e ligante (ITPs do ATB), topol.top, caixa, agua, ions, grompp [+ EM local]
#   hpc       um PBS por ligante + submit_all.sh + tar.gz
#   figures   imagens (PyMOL): visao geral e sitio de cada ligante
#
# Exemplos:
#   ./pipeline.sh examples/TcTR/config.sh all
#   ./pipeline.sh meu_config.sh dock analyze
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/scripts/lib.sh"

[ $# -ge 1 ] || { sed -n '2,20p' "$0"; exit 1; }
CONFIG=$(readlink -f "$1"); shift
[ -f "$CONFIG" ] || die "config nao encontrado: $CONFIG"
load_config "$CONFIG"

STEPS=("$@"); [ ${#STEPS[@]} -eq 0 ] && STEPS=(all)
[ "${STEPS[0]}" = "all" ] && STEPS=(check prepare grid dock analyze md hpc figures)

mkdir -p "$OUT"
exec > >(tee -a "$OUT/pipeline.log") 2>&1
log "config: $CONFIG"
log "saida:  $OUT"
log "etapas: ${STEPS[*]}"

for s in "${STEPS[@]}"; do
    case "$s" in
        check|prepare|grid|dock|analyze|md|hpc|figures)
            log "===== $s ====="
            source "$HERE/scripts/step_$s.sh"
            "step_$s"
            ;;
        *) die "etapa desconhecida: $s" ;;
    esac
done
log "concluido: ${STEPS[*]}"
