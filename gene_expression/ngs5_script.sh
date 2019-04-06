#!/bin/bash
#SBATCH --job-name=ngs5  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shaleigh.smith@nyulangone.org # Where to send mail
#SBATCH --ntasks=4 # Run on a single CPU
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --ntasks=4
#SBATCH --output=/gpfs/scratch/sas1531/ngs5_coursework/ngs5_%j.log # Standard output and error log
#SBATCH -p cpu_short

module load fastqc/0.11.7
module load trimgalore/0.5.0
module load python/cpu/2.7.15-ES ### CutAdapt is hidden in here
module load bbmap/38.25
module load samtools/1.9
module load subread/1.6.3


trim_galore --q 30 --phred33 --paired -o ./ --fastqc ./assignment/${1}.fastq.gz ./assignment/${2}.fastq.gz


bbmap.sh \
-Xmx26G \
ref=/gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
in=./${1}_val_1.fq.gz \
in2=./${2}_val_2.fq.gz \
outm=./${1}_R2.sam \
minid=0.90 \
ambiguous=random \
nodisk \

samtools view -S -b ${1}_R2.sam > ${1}_R2.bam
samtools sort ${1}_R2.bam -o ${1}_R2_sorted.bam

featureCounts -s 2 -p -B -C -P --ignoreDup --primary -a /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf -g gene_id -o ${1}_R2_feature_counts ./first_output/${1}_R2_sorted.bam
