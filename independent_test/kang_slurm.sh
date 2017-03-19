#!/bin/bash
#SBATCH -n 8
#SBATCH -c 1
#SBATCH -J Connectivity
#SBATCH -p general
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH -o SBATCH_Connectivity_out.txt
#SBATCH -e SBATCH_Connectivity_err.txt

module load NWS

SQDIR="/ysm-gpfs/project/dj333/hongyu/Project/Pathway_Connectiviy/independent_test/SQ_Files_${SLURM_JOB_ID}"
mkdir "$SQDIR"
exec > "$SQDIR/PBS_script.log"
echo "$(date +'%F %T') Batch script starting in earnest (pid: $$)."
echo "$(date +'%F %T') About to execute /ysm-gpfs/apps/software/SimpleQueue/3.1/SQDedDriver.py using task file: kang2.sh"

python "/ysm-gpfs/apps/software/SimpleQueue/3.1/SQDedDriver.py" \
  --logFile="$SQDIR/SQ.log" \
  --pnwss --wrapperVerbose \
  "kang2.sh"
RETURNCODE=$?
echo "$(date +'%F %T') Writing exited file."
touch "$SQDIR/exited"
echo "$(date +'%F %T') Batch script exiting, $RETURNCODE."

