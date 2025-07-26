# Computational Practice: Indexing a Reference Genome with HISAT2

## Lab Overview
A key aspect of measuring gene expression by sequencing (RNA-seq) is aligning reads in a FASTQ file to a reference genome. Before you can do the sequence alignment, you must index the reference genome FASTA file. In essence, this step will convert the genome into a quick look-up data structure so the alignment algorithm can map millions of reads against a giant genome very quickly. Once the FASTQ reads are mapped to the reference genome, they will be stored in SAM or BAM alignment files. Then one can use a GFF/GTF coordinate file to count the number of cDNA reads that overlap known gene coordinates (or find new ones!). See Figure 1 for an example of this workflow.

![HISAT2 Workflow](hisat2.png)

In this lab, we will prepare the human genome for mapping using HISAT2 software.

## Lab Objectives:
* Install HISAT2 on the Palmetto Cluster using conda/bioconda
* Download the human genome (FASTA) and gene coordinates file (GTF) from ENSEMBL build 114
* Index the reference genome using HISAT2
* Submit a SLURM batch job for genome indexing

## Prerequisites
* Access to Palmetto2 cluster
* Basic knowledge of Linux command line
* Understanding of FASTA and GTF file formats

## Task A: Install HISAT2 on Palmetto using Conda

### Step 1. Get an interactive node.
Access Palmetto2 Cluster via ondemand and launch a Jupyter Notebook code

### Step 2. Set up Conda Environment
```bash
# Load anaconda module
module load anaconda3/2023.09-0

# Create a new conda environment for RNA-seq analysis
conda create -n rnaseq_env python=3.9
source activate rnaseq_env

# Install HISAT2 and related tools
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda subread

# Verify installation
hisat2 --version
hisat2-build --help
```

### Step 3. Create Working Directory
```bash
# Navigate to scratch space
cd /scratch/$USER

# Create project directory
mkdir -p hisat2_indexing
cd hisat2_indexing

# Create subdirectories
mkdir -p data/genome data/annotation results logs
```

---

## Task B: Download Human Genome and GTF Files from ENSEMBL Build 114

### Step 1. Download Human Genome (Primary Assembly)
```bash
# Navigate to genome data directory
cd /scratch/$USER/hisat2_indexing/data/genome

# Download human genome GRCh38 primary assembly (build 114)
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Verify download
ls -lh Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Uncompress the genome file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Check file integrity
head -n 20 Homo_sapiens.GRCh38.dna.primary_assembly.fa
tail -n 20 Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Get file statistics
wc -l Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Step 2. Download Gene Annotation (GTF) File
```bash
# Navigate to annotation directory
cd /scratch/$USER/hisat2_indexing/data/annotation

# Download GTF annotation file (build 114)
wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz

# Verify download
ls -lh Homo_sapiens.GRCh38.114.gtf.gz

# Uncompress the GTF file
gunzip Homo_sapiens.GRCh38.114.gtf.gz

# Examine GTF structure
head -n 20 Homo_sapiens.GRCh38.114.gtf
tail -n 20 Homo_sapiens.GRCh38.114.gtf

# Count genes and transcripts
grep -c "gene_id" Homo_sapiens.GRCh38.114.gtf
grep -c "transcript_id" Homo_sapiens.GRCh38.114.gtf
```

### Step 3. Download Additional Reference Files (Optional)
```bash
# Download cDNA sequences (for validation)
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Download known splice sites (for improved alignment)
wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
```

---

## Task C: Index the Human Genome using SLURM Batch Job

### Step 1. Create SLURM Batch Script
```bash
# Navigate to main project directory
cd /scratch/$USER/hisat2_indexing

# Create SLURM script using nano editor
nano hisat2_index.slurm
```

### Step 2. SLURM Script Content
```bash
#!/bin/bash
#SBATCH -J HISAT2_INDEX
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=64GB
#SBATCH -t 24:00:00
#SBATCH -o logs/hisat2_index_%j.out
#SBATCH -e logs/hisat2_index_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@clemson.edu

# Load required modules
module load anaconda3/2023.09-0
source activate rnaseq_env

# Set variables
GENOME_DIR="/scratch/$USER/hisat2_indexing/data/genome"
ANNOTATION_DIR="/scratch/$USER/hisat2_indexing/data/annotation"
RESULTS_DIR="/scratch/$USER/hisat2_indexing/results"
GENOME_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="Homo_sapiens.GRCh38.114.gtf"
INDEX_BASE="HG38_114"

# Change to results directory
cd $RESULTS_DIR

# Print job information
echo "Job started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Working directory: $(pwd)"
echo "Genome file: $GENOME_DIR/$GENOME_FILE"
echo "GTF file: $ANNOTATION_DIR/$GTF_FILE"

# Extract splice sites and exons from GTF file
echo "Extracting splice sites and exons..."
hisat2_extract_splice_sites.py $ANNOTATION_DIR/$GTF_FILE > ${INDEX_BASE}_splice_sites.txt
hisat2_extract_exons.py $ANNOTATION_DIR/$GTF_FILE > ${INDEX_BASE}_exons.txt

# Build HISAT2 index with splice sites and exons
echo "Building HISAT2 index..."
hisat2-build \
    --ss ${INDEX_BASE}_splice_sites.txt \
    --exon ${INDEX_BASE}_exons.txt \
    -p 8 \
    $GENOME_DIR/$GENOME_FILE \
    $INDEX_BASE

# Verify index creation
echo "Index files created:"
ls -lh ${INDEX_BASE}*

# Print completion information
echo "Job completed: $(date)"
echo "Total runtime: $SECONDS seconds"

# Generate summary report
echo "=== HISAT2 Index Summary ===" > ${INDEX_BASE}_summary.txt
echo "Date: $(date)" >> ${INDEX_BASE}_summary.txt
echo "Genome: $GENOME_FILE" >> ${INDEX_BASE}_summary.txt
echo "Annotation: $GTF_FILE" >> ${INDEX_BASE}_summary.txt
echo "Index base name: $INDEX_BASE" >> ${INDEX_BASE}_summary.txt
echo "Number of splice sites: $(wc -l < ${INDEX_BASE}_splice_sites.txt)" >> ${INDEX_BASE}_summary.txt
echo "Number of exons: $(wc -l < ${INDEX_BASE}_exons.txt)" >> ${INDEX_BASE}_summary.txt
echo "Index files:" >> ${INDEX_BASE}_summary.txt
ls -lh ${INDEX_BASE}*.ht2 >> ${INDEX_BASE}_summary.txt
```

### Step 3. Submit SLURM Job
```bash
# Make sure script is executable
chmod +x hisat2_index.slurm

# Submit job to SLURM scheduler
sbatch hisat2_index.slurm

# Check job status
squeue -u $USER

# Monitor job progress
tail -f logs/hisat2_index_*.out
```

## Task D: Verify Index Creation and Quality Control

### 1. Check Index Files
```bash
# Navigate to results directory
cd /scratch/$USER/hisat2_indexing/results

# List all index files
ls -lh HG38_114*

# Expected files:
# HG38_114.1.ht2
# HG38_114.2.ht2
# HG38_114.3.ht2
# HG38_114.4.ht2
# HG38_114.5.ht2
# HG38_114.6.ht2
# HG38_114.7.ht2
# HG38_114.8.ht2
```

The indexed genome is now ready for use with HISAT2 for RNA-seq read alignment!
