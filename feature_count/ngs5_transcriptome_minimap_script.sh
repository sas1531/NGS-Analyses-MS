#!/bin/bash
#SBATCH --job-name=minimap2_transcriptome # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shaleigh.smith@nyulangone.org # Where to send mail
#SBATCH --ntasks=4 # Run on mulitple CPU
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=2:00:00 # Time limit hrs:min:sec
#SBATCH --output=/gpfs/scratch/sas1531/ngs5_coursework/minimap2_transcriptome _%j.log # Standard output and error log
#SBATCH -p cpu_short

### This script is to align nanopore data against the human transcriptome with minimap2 and generate transcript counts

module load minimap2/2.15
module load samtools/1.9
module load subread/1.6.3

# align nanpore data against human transcriptome with minimap2
# /gpfs/scratch/sas1531/ngs5_coursework/Homo_sapiens.GRCh38.all.cds.ncrna.fa

# Align with minimap2
# -ax: long reads with cigar
# -k14: kmer size
# splice: enable splice alignment
# -uf: find splice sites on transcript strand
# --secondary: Whether to output secondary alignments
minimap2 -ax splice -uf -k14 --secondary=no /gpfs/scratch/sas1531/ngs5_coursework/Homo_sapiens.GRCh38.all.cds.ncrna.fa /gpfs/scratch/sas1531/ngs5_coursework/nanopore.fastq > /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.sam

# Get stats, convert sam to bam, sort
# -S: input sam
# -b: output bam
samtools view -S -b /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.sam > /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.bam
samtools sort /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.bam -o /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.sorted.bam

### Select first 4 columns of sam, this will be imported into R
samtools view nanopore_transcriptome_secondary.sam | cut -f 1,2,3,4 > truncated_transcript_secondary.sam
