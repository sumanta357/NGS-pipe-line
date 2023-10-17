#!/bin/bash

# Define input and output directories
INPUT_DIR="input_reads"
OUTPUT_DIR="output_results"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Step 1: Quality control using FastQC
mkdir -p $OUTPUT_DIR/fastqc_reports
for file in $INPUT_DIR/*.fastq.gz; do
    fastqc -o $OUTPUT_DIR/fastqc_reports $file
done

# Step 2: Generate MultiQC report
multiqc $OUTPUT_DIR/fastqc_reports -o $OUTPUT_DIR/multiqc_report

# Step 3: Quality trimming using Trimmomatic
mkdir -p $OUTPUT_DIR/trimmed_reads
for file in $INPUT_DIR/*.fastq.gz; do
    filename=$(basename "$file")
    base="${filename%.fastq.gz}"
    trimmomatic SE -phred33 $INPUT_DIR/$base.fastq.gz \
        $OUTPUT_DIR/trimmed_reads/${base}_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Step 4: Build Bowtie2 index (only needs to be done once)
BOWTIE2_INDEX="path/to/bowtie2_index"
REFERENCE_GENOME="path/to/reference_genome.fasta"
if [ ! -f "$BOWTIE2_INDEX.1.bt2" ]; then
    bowtie2-build $REFERENCE_GENOME $BOWTIE2_INDEX
fi

# Step 5: Align reads to the reference genome using Bowtie2
mkdir -p $OUTPUT_DIR/alignment
for file in $OUTPUT_DIR/trimmed_reads/*_trimmed.fastq.gz; do
    filename=$(basename "$file")
    base="${filename%_trimmed.fastq.gz}"
    bowtie2 -x $BOWTIE2_INDEX -U $OUTPUT_DIR/trimmed_reads/${base}_trimmed.fastq.gz -S $OUTPUT_DIR/alignment/${base}.sam
done

# Step 6: Convert SAM to BAM and sort using Samtools
mkdir -p $OUTPUT_DIR/bam_files
for file in $OUTPUT_DIR/alignment/*.sam; do
    filename=$(basename "$file")
    base="${filename%.sam}"
    samtools view -b -o $OUTPUT_DIR/bam_files/${base}.bam $OUTPUT_DIR/alignment/${base}.sam
    samtools sort -o $OUTPUT_DIR/bam_files/${base}.sorted.bam $OUTPUT_DIR/bam_files/${base}.bam
done

# Step 7: Index the sorted BAM files
for file in $OUTPUT_DIR/bam_files/*.sorted.bam; do
    samtools index $file
done

# Step 8: Perform gene quantification using featureCounts
mkdir -p $OUTPUT_DIR/quantification
GTF_FILE="path/to/annotation.gtf"
for file in $OUTPUT_DIR/bam_files/*.sorted.bam; do
    filename=$(basename "$file")
    base="${filename%.sorted.bam}"
    featureCounts -a $GTF_FILE -o $OUTPUT_DIR/quantification/${base}_counts.txt $file
done

# Clean up intermediate files if needed
# rm -r $OUTPUT_DIR/fastqc_reports
# rm -r $OUTPUT_DIR/trimmed_reads
# rm -r $OUTPUT_DIR/alignment
# rm -r $OUTPUT_DIR/bam_files

echo "RNA-seq pipeline for single-end data completed!"
