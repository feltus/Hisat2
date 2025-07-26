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

---

## Step A: Install HISAT2 on Palmetto using Conda

### 1. Access Palmetto2 Cluster
```bash
# Log into Palmetto2
ssh your_username@ondemand.rcd.clemson.edu
```

### 2. Request Interactive Computing Resources
```bash
# Request an interactive node
srun -p interactive --pty bash
```

### 3. Set up Conda Environment
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

### 4. Create Working Directory
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

## Step B: Download Human Genome and GTF Files from ENSEMBL Build 114

### 1. Download Human Genome (Primary Assembly)
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

### 2. Download Gene Annotation (GTF) File
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

### 3. Download Additional Reference Files (Optional)
```bash
# Download cDNA sequences (for validation)
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Download known splice sites (for improved alignment)
wget https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
```

---

## Step C: Index the Human Genome using SLURM Batch Job

### 1. Create SLURM Batch Script
```bash
# Navigate to main project directory
cd /scratch/$USER/hisat2_indexing

# Create SLURM script using nano editor
nano hisat2_index.slurm
```

### 2. SLURM Script Content
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

### 3. Submit SLURM Job
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

### 4. Alternative: Quick Test Index (for testing)
```bash
# For quick testing, create a smaller index with limited chromosomes
# Create test script
nano hisat2_test_index.slurm
```

```bash
#!/bin/bash
#SBATCH -J HISAT2_TEST
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16GB
#SBATCH -t 2:00:00
#SBATCH -o logs/hisat2_test_%j.out
#SBATCH -e logs/hisat2_test_%j.err

# Load modules
module load anaconda3/2023.09-0
source activate rnaseq_env

cd /scratch/$USER/hisat2_indexing/results

# Extract chromosome 22 only for testing
grep "^>22\|^22" /scratch/$USER/hisat2_indexing/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa > chr22_test.fa

# Build test index
hisat2-build -p 4 chr22_test.fa HG38_chr22_test

# List results
ls -lh HG38_chr22_test*
```

---

## Step D: Verify Index Creation and Quality Control

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

### 2. Test Index with Sample Alignment
```bash
# Test the index with a small sample (if you have FASTQ files)
# This is just an example - you'll need actual FASTQ files
hisat2 -x HG38_114 -U sample.fastq -S test_alignment.sam
```

### 3. Generate Final Report
```bash
# Create comprehensive report
cat > final_report.txt << EOF
=== HISAT2 Human Genome Index Report ===
Date: $(date)
User: $USER
Working Directory: $(pwd)

Genome Information:
- Source: ENSEMBL Release 114
- Assembly: GRCh38
- File: Homo_sapiens.GRCh38.dna.primary_assembly.fa

Annotation Information:
- GTF File: Homo_sapiens.GRCh38.114.gtf
- Splice Sites: $(wc -l < HG38_114_splice_sites.txt) sites
- Exons: $(wc -l < HG38_114_exons.txt) exons

Index Files:
$(ls -lh HG38_114*.ht2)

Total Index Size: $(du -sh HG38_114*.ht2 | tail -n 1)

Job Information:
- SLURM Job ID: Check logs/hisat2_index_*.out
- Runtime: Check job completion time in logs
- Resources Used: 8 cores, 64GB RAM

Next Steps:
1. Use this index for RNA-seq alignment with hisat2
2. Convert alignments to BAM format with samtools
3. Count reads per gene with featureCounts
EOF

echo "Index creation complete! Check final_report.txt for summary."
```

---

## Key Differences from Previous Version

### 1. **Installation Method**
- **Old**: Manual download and installation
- **New**: Conda/bioconda package management (more reliable, includes dependencies)

### 2. **ENSEMBL Version**
- **Old**: Release 101
- **New**: Release 114 (current as of 2024)

### 3. **Batch System**
- **Old**: PBS (Portable Batch System)
- **New**: SLURM (Simple Linux Utility for Resource Management)

### 4. **Enhanced Features**
- Added splice site and exon extraction for better alignment
- Improved resource allocation
- Better error handling and logging
- Comprehensive quality control steps

### 5. **File Organization**
- Better directory structure
- Separate data, results, and logs directories
- More comprehensive documentation

---

## Troubleshooting

### Common Issues and Solutions

1. **Memory Issues**
   - Increase memory allocation in SLURM script
   - Use `--mem=128GB` for large genomes

2. **Time Limits**
   - Genome indexing can take 4-12 hours
   - Adjust walltime accordingly

3. **Module Loading**
   - Ensure anaconda module is available
   - Check module availability with `module avail`

4. **File Permissions**
   - Ensure files are readable: `chmod +r filename`
   - Check disk space: `df -h /scratch`

### Performance Tips

1. **Parallel Processing**
   - Use `-p 8` or higher for faster indexing
   - Balance cores vs. memory usage

2. **Storage Location**
   - Use `/scratch` for large files
   - Clean up intermediate files after completion

3. **Index Reuse**
   - Share index files among project members
   - Document index version and parameters

---

## Next Steps

After completing this lab, you should be able to:
1. Use the created index for RNA-seq alignment
2. Process multiple samples in batch
3. Integrate with downstream analysis pipelines
4. Understand the relationship between genome indexing and alignment performance

The indexed genome is now ready for use with HISAT2 for RNA-seq read alignment!
