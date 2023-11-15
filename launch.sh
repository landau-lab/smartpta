#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=nextflow
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jzinno@nygenome.org
#SBATCH --output=darkshore-log_%j.out



module load nextflow/22.10.4


nextflow main.nf -with-report report-nextflow-log.html -with-dag flowchart.html -with-timeline timeline.html -resume