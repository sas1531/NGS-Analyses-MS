#!/bin/bash
#SBATCH --job-name=minimap2_genome # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=shaleigh.smith@nyulangone.org # Where to send mail
#SBATCH --ntasks=4 # Run on mulitple CPU
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=2:00:00 # Time limit hrs:min:sec
#SBATCH --output=/gpfs/scratch/sas1531/ngs5_coursework/mimimap2_genome_%j.log # Standard output and error log
#SBATCH -p cpu_short

### This script is to align nanopore data against the human genome with minimap2 and generate gene counts

module load minimap2/2.15
module load samtools/1.9
module load subread/1.6.3

# Align with minimap2
# -ax: long reads with cigar
# -k14: kmer size
# --cs=long: output the cigar string (not necessary but could be interesting)
# splice: enable splice alignment
# -uf: find splice sites on transcript strand
# --secondary: Whether to output secondary alignments
minimap2 -ax splice -uf -k14 --secondary=no /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa /gpfs/scratch/sas1531/ngs5_coursework/nanopore.fastq > /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome.sam

# Get stats, convert sam to bam, sort
# -S: input sam
# -b: output bam
# --secondary: Whether to output secondary alignments [yes]
samtools view -S -b /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome.sam > /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome.bam
samtools sort /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome.bam -o /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome.sorted.bam

# Generate gene counts (this is NOT paired)
# -a: annotation file
# -s: strand-specific read counting - 1
# -L: used when counting long reads
# -M: multi-mapping reads/fragments will be counted
# -g: specifies the attribute to group
# -O: reads (or fragments if -p is specified) will be al- lowed to be assigned to more than one matched meta-feature
featureCounts -s 1 -L -a /gpfs/scratch/sas1531/hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf -g gene_id -o /gpfs/scratch/sas1531/ngs5_coursework/feature_counts_genome /gpfs/scratch/sas1531/ngs5_coursework/nanopore_genome.sorted.bam
