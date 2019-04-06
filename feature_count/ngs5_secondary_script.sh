#!/bin/bash
#SBATCH --job-name=feature_secondary # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shaleigh.smith@nyulangone.org # Where to send mail
#SBATCH --ntasks=4 # Run on mulitple CPU
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --output=/gpfs/scratch/sas1531/ngs5_coursework/feature_secondary_%j.log # Standard output and error log
#SBATCH -p cpu_short

### This script is to align nanopore data against the human transcriptome with minimap2 and generate transcript counts

module load minimap2/2.15
module load samtools/1.3
module load subread/1.6.3

# align nanpore data against human transcriptome with minimap2
# /gpfs/scratch/sas1531/ngs5_coursework/Homo_sapiens.GRCh38.all.cds.ncrna.fa

# Align with minimap2
# -ax: long reads with cigar
# -k14: kmer size
# splice: enable splice alignment
# -uf: find splice sites on transcript strand
# --secondary: Whether to output secondary alignments [yes]
minimap2 -ax splice -uf -k14 --secondary=no /gpfs/scratch/sas1531/ngs5_coursework/Homo_sapiens.GRCh38.all.cds.ncrna.fa /gpfs/scratch/sas1531/ngs5_coursework/nanopore.fastq > /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.sam
minimap2 -ax splice -uf -k14 --secondary=no /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa /gpfs/scratch/sas1531/ngs5_coursework/nanopore.fastq > /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome_secondary.sam
# Secondary? --secondary=no

# Get stats, convert sam to bam, sort
# -S: input sam
# -b: output bam
samtools view -S -b /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.sam > /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.bam
samtools sort /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.bam -o /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.sorted.bam

samtools view -S -b /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome_secondary.sam > /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome_secondary.bam
samtools sort /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome_secondary.bam -o /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome_secondary.sorted.bam


# Generate transcript counts (this is NOT paired)
# -a: annotation file
# -s: strand-specific read counting - 1
# -L: used when counting long reads
# -M: multi-mapping reads/fragments will be counted
# -O: reads (or fragments if -p is specified) will be al- lowed to be assigned to more than one matched meta-feature
featureCounts -M -s 1 -L -a /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf -g transcript_id -o /gpfs/scratch/sas1531/ngs5_coursework/feature_counts_transcriptome_secondary /gpfs/scratch/sas1531/ngs5_coursework/nanopore_transcriptome_secondary.sorted.bam
featureCounts -M -s 1 -L -a /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf -g gene_id -o /gpfs/scratch/sas1531/ngs5_coursework/feature_counts_genome_secondary /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome_secondary.sorted.bam


### Altering and counting the sam
samtools view nanopore_transcriptome_secondary.sam | cut -f 1,3 > truncated_transcript_secondary.sam

# Convert transcript bam to bed
bedtools bamtobed -i nanopore_transcriptome_secondary.sorted.bam > nanopore_transcriptome_secondary.bed
