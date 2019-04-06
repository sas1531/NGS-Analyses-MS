### Command line code for  manipulation

module load samtools/1.3
module load bedtools/2.27.1

### Altering and counting the sam
samtools view nanopore_transcriptome_secondary.sam | cut -f 1,2,3,4 > truncated_transcript_secondary.sam


#### Don't need this stuff:
### Alter gtf
gtf2bed < genes.gtf > genes.bed
cut -f 1,2,3,4,8 genes.bed > genes_simple.bed


# Convert transcript bam to bed
bedtools bamtobed -i nanopore_transcriptome.sorted.bam > nanopore_transcriptome.bed
