#!/bin/bash
#SBATCH --job-name=tet3_rna_seq_processing
#SBATCH --output=tet3_rna_seq_processing_%j.out
#SBATCH --error=tet3_rna_seq_processing_%j.err
#SBATCH --time=48:00:00       # Adjust time as needed
#SBATCH --mem-per-cpu=3988    # Memory per CPU core
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24    # Number of CPU cores per task

# Load necessary modules
module load Mamba/4.14.0-0

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bioinf

# Ensure project bin is in PATH (/path/to/project/bin)
export PATH="/path/to/project/bin:$PATH"

# Set up environment
set -e
set -o pipefail

echo "Job started on $(date)"
echo "Running on node $(hostname)"

# Define variables
dir="/path/to/project/rnaseq"
genome="/path/to/genome/v4_Mt_Pt_TET3/S_lycopersicum_chromosomes.4.00_Mt_Pt_TET3.fa"
annotation="/path/to/genome/v4_Mt_Pt_TET3/ITAG4.1_gene_models_Mt_Pt_TET3.gtf"

# Create necessary directories
mkdir -p "$dir/qa/1_raw"
mkdir -p "$dir/qa/2_trimmed"
mkdir -p "$dir/data/2_trimmed"
mkdir -p "$dir/data/3_mapped"

# Function for FastQC on raw data
run_fastqc_raw() {
    echo "Running FastQC on raw data..."
    cd "$dir"
    fastqc -t 16 --quiet -o qa/1_raw/ data/1_raw/*fq.gz
    multiqc --filename qa/1_raw/01_raw_FastQC_multiqc qa/1_raw
    echo "FastQC on raw data completed."
}

# Function for trimming with Fastp
run_fastp() {
    echo "Running Fastp trimming..."
    cd "$dir/data/1_raw"
    for i in *.R1.fq.gz; do
        fastp --thread 16 --cut_right --cut_window_size 4 --cut_mean_quality 20 --length_required 40 \
            --overrepresentation_analysis \
            -i "$i" -o "../2_trimmed/${i/.R1./.R1.trimmed.}" \
            -h "../../qa/2_trimmed/${i%.R1.*}_fastp.html" \
            -j "../../qa/2_trimmed/${i%.R1.*}_fastp.json"
    done
    echo "Fastp trimming completed."
    # MultiQC report for Fastp
    cd "$dir"
    multiqc --filename qa/2_trimmed/02_trimmed_fastp_multiqc qa/2_trimmed -f
}

# Function for FastQC on trimmed data
run_fastqc_trimmed() {
    echo "Running FastQC on trimmed data..."
    cd "$dir"
    fastqc -t 16 --quiet -o qa/2_trimmed/ data/2_trimmed/*.fq.gz
    multiqc --filename qa/2_trimmed/02_trimmed_FastQC_multiqc qa/2_trimmed
    echo "FastQC on trimmed data completed."
}

# Function for STAR genome index generation
generate_star_index() {
    echo "Generating STAR genome index..."
    STAR --runMode genomeGenerate \
         --runThreadN 24 \
         --genomeDir "$(dirname "$genome")" \
         --genomeFastaFiles "$genome" \
         --sjdbGTFfile "$annotation" \
         --sjdbOverhang 74 \
         --genomeSAindexNbases 13
    echo "STAR genome index generation completed."
}

# Function for mapping with STAR
run_star_mapping() {
    echo "Running STAR mapping..."
    cd "$dir/data/2_trimmed"
    time_file="../timing_mapping_STAR.tsv"
    printf "sample\tstart\tend\n" > "$time_file"
    for i in *R1.trimmed.fq.gz; do
        sample_name="${i%.R1.trimmed.fq.gz}"
        printf "$sample_name\t$(date '+%Y/%m/%d %H:%M:%S')\t" >> "$time_file"
        STAR --runThreadN 12 \
             --genomeDir "$(dirname "$genome")" \
             --readFilesIn "$i" \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "../3_mapped/${sample_name}." \
             --outWigType bedGraph \
             --outWigStrand Unstranded \
             --quantMode GeneCounts
        printf "$(date '+%Y/%m/%d %H:%M:%S')\n" >> "$time_file"
    done
    echo "STAR mapping completed."
}

# Execute functions
# Uncomment the functions you want to run
#run_fastqc_raw
#run_fastp
#run_fastqc_trimmed
generate_star_index
run_star_mapping

echo "All steps completed."
echo "Job finished on $(date)"
