#!/bin/bash
#SBATCH --job-name=trimgalore_test # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Shaleigh.Smith@nyumc.org # Where to send mail
#SBATCH --ntasks=4 # Run on multiple CPU
#SBATCH --mem=4gb # Job memory request
#SBATCH --time=6:00:00 # Time limit hrs:min:sec
#SBATCH --output=trim_galore_%j.log # Standard output and error log
#SBATCH --partition=cpu_short # Which partition you want to run it on

# This script is used to trim sequence reads using trimgalore

module load trimgalore/0.5.0
module load python/cpu/2.7.15-ES ### CutAdapt is hidden in here
module load fastqc

# Trim sequences with PHRED quality score above 30 and remove adapters (trim galore will automatically detect these), run FASTQC on output:

trim_galore --q 20 --phred33 -o /gpfs/scratch/sas1531/ngs2_coursework/ --fastqc /gpfs/scratch/sas1531/ngs2_coursework/ERR218285_3.fastq.gz

trim_galore --q 20 --phred33 --paired -o /gpfs/scratch/sas1531/ngs2_coursework/ --fastqc /gpfs/scratch/sas1531/ngs2_coursework/SRR1523657_1.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRR1523657_2.fastq.gz

trim_galore --q 20 --phred33 --small_rna --paired -o /gpfs/scratch/sas1531/ngs2_coursework/ --fastqc /gpfs/scratch/sas1531/ngs2_coursework/SRX747060_1.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRX747060_2.fastq.gz

# Don't need to trim these:
# ERR218285_1.fastq.gz
# ERR218285_2.fastq.gz

### For command help
# trim_galore --help

### Submit job using sbatch
# sbatch trimgalore_script.sh

### View queue
# squeue -u sas1531

### Export html files using filezilla
