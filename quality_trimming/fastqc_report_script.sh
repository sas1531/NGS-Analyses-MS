#!/bin/bash
#SBATCH --job-name=fastqc_job # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=first.last@nyumc.org # Where to send mail
#SBATCH --ntasks=1 # Run on mulitple CPU
#SBATCH --mem=4gb # Job memory request
#SBATCH --time=02:00:00 # Time limit hrs:min:sec
#SBATCH --output=fastqc_%j.log # Standard output and error log
#SBATCH --p=cpu_short # Which partition you want to run it on

# This script is used to generate an FASTQC report from a fastqc.gz files

module load fastqc

fastqc -o /gpfs/scratch/sas1531/ /gpfs/scratch/sas1531/SRR1523657_1.fastq.gz /gpfs/scratch/sas1531/SRR1523657_2.fastq.gz /gpfs/scratch/sas1531/SRX747060_1.fastq.gz /gpfs/scratch/sas1531/SRX747060_2.fastq.gz /gpfs/scratch/sas1531/ERR218285_1.fastq.gz /gpfs/scratch/sas1531/ERR218285_2.fastq.gz /gpfs/scratch/sas1531/ERR218285_3.fastq.gz

### Submit job using sbatch
# sbatch fastqc_report_script.sh

### View queue
# squeue -u sas1531

### Export html files using filezilla
