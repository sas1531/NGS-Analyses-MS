#!/bin/bash
#SBATCH --job-name=sra_job # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=first.last@nyumc.org # Where to send mail
#SBATCH --ntasks=4 # Run on mulitple CPU
#SBATCH --mem=4gb # Job memory request
#SBATCH --time=06:00:00 # Time limit hrs:min:sec
#SBATCH --output=sra_%j.log # Standard output and error log
#SBATCH --partition=cpu_short # Which partition you want to run it on

# This script is used to download fastq data directly from SRA ncbi, split the data into two fastq files and compress the results files to a gzip (fastq.gz)

module load sratoolkit/2.9.1

fastq-dump --split-files SRR7992453 --gzip -O /gpfs/scratch/sas1531/
fastq-dump --split-files SRX747060 --gzip -O /gpfs/scratch/sas1531/
fastq-dump --split-files ERR218285 --gzip -O /gpfs/scratch/sas1531/

rm -r ~/ncbi # fastq-dump creates a temp dir but doesn't remove it... plays merry havoc with home storage space

### Submit job using sbatch
# sbatch sra_download_script.sh

### View queue
# squeue -u sas1531
