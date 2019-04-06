#!/bin/bash
#SBATCH --job-name=trimmomatic_job_test # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Shaleigh.Smith@nyumc.org # Where to send mail
#SBATCH --ntasks=1 # Run on a single CPU
#SBATCH --mem=4gb # Job memory request
#SBATCH --time=06:00:00 # Time limit hrs:min:sec
#SBATCH --output=trimmomatic_%j.log # Standard output and error log
#SBATCH --partition=cpu_short # Which partition you want to run it on

# This script is used to trim sequence reads using trimmomatic

module load trimmomatic/0.36
module load fastqc

#Single End:
java -jar /gpfs/share/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 /gpfs/scratch/sas1531/ngs2_coursework/ERR218285_3.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/ERR218285_3_trimmomatic_out.fastq.gz ILLUMINACLIP:/gpfs/scratch/sas1531/ngs2_coursework/adapter.fasta:2:30:10 SLIDINGWINDOW:4:20

#Paired End:
java -jar /gpfs/share/apps/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 /gpfs/scratch/sas1531/ngs2_coursework/SRR1523657_1.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRR1523657_2.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRR1523657_1_trimmomatic_paired_out.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRR1523657_1_trimmomatic_unpaired_out.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRR1523657_2_trimmomatic_paired_out.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRR1523657_2_trimmomatic_unpaired_out.fastq.gz ILLUMINACLIP:/gpfs/scratch/sas1531/ngs2_coursework/adapter.fasta:2:30:10 SLIDINGWINDOW:4:20

java -jar /gpfs/share/apps/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 /gpfs/scratch/sas1531/ngs2_coursework/SRX747060_1.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRX747060_2.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRX747060_1_trimmomatic_paired_out.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRX747060_1_trimmomatic_unpaired_out.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRX747060_2_trimmomatic_paired_out.fastq.gz /gpfs/scratch/sas1531/ngs2_coursework/SRX747060_2_trimmomatic_unpaired_out.fastq.gz ILLUMINACLIP:/gpfs/scratch/sas1531/ngs2_coursework/adapter.fasta:2:30:10 SLIDINGWINDOW:4:20

# Don't need to trim these:
# ERR218285_1.fastq.gz
# ERR218285_2.fastq.gz

### Put each output file through fastqc to generate report
fastqc -o /gpfs/scratch/sas1531/ngs2_coursework /gpfs/scratch/sas1531/ngs2_coursework/*trimmomatic*.gz

### Manual
# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

### Submit job using sbatch
# sbatch fastqc_report_script.sh

### View queue
# squeue -u sas1531

### Export html files using filezilla
