#!/bin/bash
#SBATCH --job-name=ngs6  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shaleigh.smith@nyulangone.org # Where to send mail
#SBATCH --ntasks=8 # Run on a multiple CPU
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=12:00:00 # Time limit hrs:min:sec
#SBATCH --output=/gpfs/scratch/sas1531/ngs6_coursework/test/ngs6_%j.log # Standard output and error log
#SBATCH -p cpu_short

# Load modules
module load sratoolkit/2.9.1
module load fastqc/0.11.7
module load trimgalore/0.5.0
module load python/cpu/2.7.15-ES ### CutAdapt is hidden in here
module load bowtie2/2.3.4.1
module load samtools/1.9

# Download datasets
fastq-dump ${1} --gzip -O /gpfs/scratch/sas1531/ngs6_coursework/

# Remove fastq-dump directory
rm -r ~/ncbi

# Rename files
mv /gpfs/scratch/sas1531/ngs6_coursework/${1}*fastq.gz /gpfs/scratch/sas1531/ngs6_coursework/${2}.fastq.gz

# Trim
# These are single end reads (not paired) sequenced using RIP-Seq
trim_galore --q 30 \
--phred33 \
-o /gpfs/scratch/sas1531/ngs6_coursework/ \
--fastqc ./${2}.fastq.gz

# Build bowtie2 index (done once)
#bowtie2-build -f /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa hg38_index

# Align to hg38
bowtie2 -q \
--end-to-end \
--very-sensitive \
--no-unal \
-x /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
-U ./${2}_trimmed.fq.gz \
-S ./${2}.sam


# Convert to bam and sort
samtools view -S -b ${2}.sam > ${2}.bam
samtools sort ${2}.bam -o ${2}_sorted.bam
samtools index ${2}_sorted.bam
