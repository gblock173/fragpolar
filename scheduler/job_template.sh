#!/bin/bash
# Generic batch-job template (PBS/Torque syntax shown; adapt the directives to your scheduler, e.g. SLURM #SBATCH).
# Usage examples (see README "Run order"):
#   single job : qsub -q screen -v SCRIPT=pipelines/01_group_effects/ge1_aggregate.py,ARGS="full 8" scheduler/job_template.sh
#   array job  : qsub -q screen -J 1-10 -v SCRIPT=pipelines/05_direct_fragments/df1_extract.py,ARGS='$ARRAY_INDEX' scheduler/job_template.sh
#                (per-sample scripts accept a 1-based task index of the sample table, or a sample identifier)
#PBS -N fragpolar
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -l walltime=06:00:00
#PBS -j oe
set -euo pipefail
cd "${PBS_O_WORKDIR:-$(pwd)}"                      # repository root
# activate the environment created from environment.yml, e.g.:  conda activate fragpolar
export FRAGPOLAR_CONFIG="${FRAGPOLAR_CONFIG:-$(pwd)/config/config.yaml}"
export OPENBLAS_NUM_THREADS=${THREADS:-1} OMP_NUM_THREADS=${THREADS:-1} MKL_NUM_THREADS=${THREADS:-1} PYTHONDONTWRITEBYTECODE=1
ARRAY_INDEX="${PBS_ARRAY_INDEX:-${SLURM_ARRAY_TASK_ID:-${TASK:-}}}"
ARGS="${ARGS//\$ARRAY_INDEX/$ARRAY_INDEX}"
python -u "$SCRIPT" $ARGS
