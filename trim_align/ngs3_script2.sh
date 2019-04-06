#!/bin/bash
#SBATCH --job-name=ngs3_2  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shaleigh.smith@nyulangone.org # Where to send mail
#SBATCH --ntasks=4 # Run on multiple CPU
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=ngs3_%j.log # Standard output and error log
#SBATCH -p cpu_medium

### Script for NGS Coursework 3

### Load Modules
module load sratoolkit/2.9.1
module load fastqc/0.11.7
module load trimgalore/0.5.0
module load python/cpu/2.7.15-ES ### CutAdapt is hidden in here
module load bbmap/38.25
module load samtools/1.3
module load bedtools/2.26.0

### Download datasets
fastq-dump --split-files SRR1523666 --gzip -O /gpfs/scratch/sas1531/ngs3_coursework/
rm -r ~/ncbi # fastq-dump creates a temp dir that needs to be removed

### Run fastQC on datasets
fastqc -o /gpfs/scratch/sas1531/ngs3_coursework/ /gpfs/scratch/sas1531/ngs3_coursework/SRR1523666_1.fastq.gz /gpfs/scratch/sas1531/ngs3_coursework/SRR1523666_2.fastq.gz

### Trim datasets and run fastQC again
trim_galore --q 20 --phred33 --paired -o /gpfs/scratch/sas1531/ngs3_coursework/clean --fastqc /gpfs/scratch/sas1531/ngs3_coursework/SRR1523666_1.fastq.gz /gpfs/scratch/sas1531/ngs3_coursework/SRR1523666_2.fastq.gz

######### SRR1523666
### Align dataset against the human genome
bbmap.sh -Xmx26G ref=/gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa in=/gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_1_val_1.fq.gz in2=/gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_2_val_2.fq.gz outm=/gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sam minid=0.90 ambiguous=random nodisk

### Parse alignment to generate sorted and indexed bam files
samtools view -b -o /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.bam /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sam
samtools sort -o /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sorted.bam /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.bam
samtools index /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sorted.bam

### Parse for forward strand and generate sorted and indexed bam files
samtools view -b -f99 /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sorted.bam > /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.for2.bam
samtools view -b -f147 /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sorted.bam > /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.for1.bam
samtools merge -f /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out_sorted.forward.bam /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.for1.bam /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.for2.bam
samtools index /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out_sorted.forward.bam

### Parse for reverse strand and generate sorted and indexed bam files
samtools view -b -f83 /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sorted.bam > /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.rev2.bam
samtools view -b -f163 /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sorted.bam > /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.rev1.bam
samtools merge -f /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out_sorted.reverse.bam /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.rev1.bam /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.rev2.bam
samtools index /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out_sorted.reverse.bam

### Parse and generate bedgraphs for Gvis
samtools view -b /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out_sorted.forward.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out_sorted.forward.bedgraph
samtools view -b /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out_sorted.reverse.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out_sorted.reverse.bedgraph
samtools view -b /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sorted.bam | genomeCoverageBed -ibam stdin -bg -split -g /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa > /gpfs/scratch/sas1531/ngs3_coursework/clean/SRR1523666_out.sorted.bedgraph

### Samtools Notes
# Explanantion for sam flags: https://broadinstitute.github.io/picard/explain-flags.html
# Explanantion of sam paired flags: https://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/

### Submit job using sbatch
# sbatch ngs3_script2.sh

### View queue
# squeue -u sas1531
