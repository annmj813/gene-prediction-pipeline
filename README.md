# gene-prediction-pipeline
This project focuses on improving the gene annotation of Cryptococcus neoformans var. grubii H99, a major fungal pathogen, using the BRAKER3 pipeline. By integrating RNA-seq evidence with ab initio predictions, we generated an updated gene model and compared it against the existing annotation from FungiDB.

The initial step for genome annotation is to download the genome and SRA files. 
## Downloading Genome from FungiDB

<pre> ```wget https://fungidb.org/common/downloads/Current_Release/CneoformansH99/fasta/data/FungiDB-68_CneoformansH99_Genome.fasta ``` </pre>

## Downloading SRA Files
To download SRA files, 
<pre>#!/bin/bash

# Load the SRA Toolkit module 
module load sra-tools

# List of SRA accession numbers 
SRA_LIST="sra_accessions.txt"

# Directory to store the downloaded data
DOWNLOAD_DIR="/path/to/your/output_directory"

# Loop through each accession and download + convert
while read -r SRA_ID; do
    echo "Processing $SRA_ID..."
    
    # Download the SRA file
    prefetch "$SRA_ID"
    
    # Convert to FASTQ using fasterq-dump
    fasterq-dump "$SRA_ID" --split-files -e 6 -O "$DOWNLOAD_DIR"
    
    echo "$SRA_ID processing done."
done < "$SRA_LIST"

echo "All downloads and conversions complete."</pre>

## Quality check by FastQC

After downloading all the datas, the next is to check the quality of the SRA files, for that we are doing FastQC. 

<pre> ``` #!/bin/bash

# List of SRA accession numbers
SRA_LIST="sra_accessions.txt"

# Output directories
FASTQ_DIR="/path/to/fastq_output"
FASTQC_DIR="/path/to/fastqc_output"
mkdir -p "$FASTQ_DIR" "$FASTQC_DIR"

# Load modules 
module load sra-tools
module load fastqc

# Process each SRA ID
while read -r SRA_ID; do
    echo "Processing $SRA_ID..."

    # Step 1: Convert to FASTQ
    fasterq-dump "$SRA_ID" --split-files -e 6 -O "$FASTQ_DIR"

    # Step 2: Run FastQC on all resulting FASTQ files
    for fq in "$FASTQ_DIR"/"$SRA_ID"_*.fastq; do
        fastqc "$fq" -o "$FASTQC_DIR"
    done

    echo "$SRA_ID done."
done < "$SRA_LIST"

echo "All SRA samples processed with FastQC." ``` </pre>

### Trimming of SRA files
After analising the quality of the SRA files, trimming should be done if neccessary.

<pre> ``` #!/bin/bash

# List of SRA accession numbers
SRA_LIST="sra_accessions.txt"

# Output directories
FASTQ_DIR="/path/to/fastq"
TRIM_DIR="/path/to/trimmed"
ADAPTERS="/path/to/Trimmomatic/adapters/TruSeq3-PE.fa"

# Create output folders
mkdir -p "$FASTQ_DIR" "$TRIM_DIR"

# Load required tools (adjust for your HPC environment)
module load sra-tools
module load trimmomatic

# Loop through each SRA sample
while read -r SRA_ID; do
    echo "Processing $SRA_ID..."

    # Step 1: Download and convert SRA to FASTQ
    fasterq-dump "$SRA_ID" --split-files -e 6 -O "$FASTQ_DIR"

    # Define input/output file names
    R1="$FASTQ_DIR/${SRA_ID}_1.fastq"
    R2="$FASTQ_DIR/${SRA_ID}_2.fastq"
    P1="$TRIM_DIR/${SRA_ID}_1_paired.fq.gz"
    U1="$TRIM_DIR/${SRA_ID}_1_unpaired.fq.gz"
    P2="$TRIM_DIR/${SRA_ID}_2_paired.fq.gz"
    U2="$TRIM_DIR/${SRA_ID}_2_unpaired.fq.gz"

    # Step 2: Run Trimmomatic for paired-end reads
    trimmomatic PE -threads 6 \
        "$R1" "$R2" \
        "$P1" "$U1" "$P2" "$U2" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "Trimming done for $SRA_ID"
done < "$SRA_LIST"

echo "All samples trimmed." ``` </pre>


## Softmasking the genome
Before gene prediction, the genome is soft-masked to identify and mask repetitive elements using RepeatModeler and RepeatMasker. This helps BRAKER3 avoid falsely predicting genes in repetitive regions. To keep the process organized and trackable, it's recommended to write the softmasking commands into a job script (e.g., using `nano`), and submit it to your HPC scheduler.

<pre> ``` # soft masking
#!/bin/bash
#SBATCH --job-name=softmask_cneo
#SBATCH --output=softmask_%j.out
#SBATCH --error=softmask_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=100G

set -e

DB=FungiDB68_Cneo_H99_repeatdb
GENOME=FungiDB-68_CneoformansH99_Genome.fasta
SIF=dfam-tetools-latest.sif
WORKDIR=$(pwd)

echo "Starting softmasking in $WORKDIR"

# Mount current working directory into container at /mnt and work from there
singularity exec -B ${WORKDIR}:/mnt $SIF bash -c "
  cd /mnt && \
  BuildDatabase -name ${DB} ${GENOME} && \
  RepeatModeler -database ${DB} -threads 32 -LTRStruct && \
  RepeatMasker -pa 32 -lib ${DB}-families.fa -xsmall ${GENOME} && \
  mv ${GENOME}.masked ${GENOME%.fasta}.softmasked.fasta
"

echo "Softmasking complete: Output = ${GENOME%.fasta}.softmasked.fasta" ``` </pre>

## 🎯 RNA-Seq Alignment Using HISAT2

Once the genome is soft-masked and the RNA-seq reads are trimmed, the next step is to align those reads to the reference genome using **HISAT2**, a fast and memory-efficient aligner for spliced reads.

This step helps prepare the alignment file (in BAM format) that will be used later for training the gene prediction model (BRAKER3).

Below is a job script to run HISAT2 and convert the output to BAM using `samtools`:

<pre>  
#!/bin/bash -l
#SBATCH -J hisat2_align
#SBATCH -o slurm_%j.out
#SBATCH -e slurm_%j.err
#SBATCH -t 12:00:00
#SBATCH -c 8

# Load required modules
module load hisat2
module load samtools

# Define file names and prefixes
GENOME=genome.cleaned.fa
IDX=genome_hisat2_index
NAME=rnaseq_alignment

# Build index if it doesn't already exist
if [ ! -e ${IDX}.1.ht2 ]; then
    echo "Building HISAT2 index..."
    hisat2-build ${GENOME} ${IDX}
fi

# Collect trimmed paired-end read files
FWD_FILES=$(ls SRR*_fpaired.fq.gz | paste -sd, -)
REV_FILES=$(ls SRR*_rpaired.fq.gz | paste -sd, -)

# Run HISAT2 and save directly as BAM file
hisat2 -p 8 -q -x ${IDX} \
  -1 ${FWD_FILES} -2 ${REV_FILES} \
  2> ${NAME}.err | samtools view -@ 8 -bS - > ${NAME}.bam


squeue -u $USER # to chcek the status, R = running
less rnaseq_alignment.err - # gives the log of data </pre>

## Gene Prediction Using BRAKER3
After generating a sorted BAM file from the RNA-seq alignments and softmasked genome, BRAKER3 is used to predict gene models by integrating ab initio predictions with RNA-seq evidence and optional protein hints.

<pre> #!/bin/bash -l
#SBATCH -J braker3
#SBATCH -o slurm_%j_braker3.out
#SBATCH -e slurm_%j_braker3.err
#SBATCH -t 48:00:00
#SBATCH -c 32

T=32
SORTED_BAM=rnaseq_alignment_sorted.bam
GENOME=genome.cleaned.fa
PROT_DB=Fungi.fa
WD=$(basename -s .bam ${SORTED_BAM})_$(basename -s .fa ${PROT_DB})

singularity exec -B ${PWD}:${PWD},${HOME} ${HOME}/braker3.sif braker.pl \
  --genome=${GENOME} \
  --prot_seq=${PROT_DB} \
  --bam=${SORTED_BAM} \
  --workingdir=${WD} \
  --threads=${T} \
  --gff3
</pre>

## Gene Model Completeness Assessment with BUSCO
To assess the completeness of predicted proteins, BUSCO is used with the Basidiomycota lineage dataset.

<pre>busco -i braker.aa \
  -l basidiomycota_odb10 \
  -o busco_output \
  -m protein \
  --cpu 8</pre>
Output: Summary of complete, fragmented, and missing orthologs in busco_output/

## Annotation Comparison with GFFCompare
To compare the structural similarity between BRAKER3 and the reference (FungiDB) annotation:

<pre># Convert GFF3 to GTF (if needed)
gffread braker.gff3 -T -o braker.gtf
gffread fungidb_reference.gff3 -T -o fungidb_reference.gtf

# Run GFFCompare
gffcompare -r fungidb_reference.gtf -o compare_out braker.gtf </pre>





