#!/bin/bash
#SBATCH --job-name=tet3_rna_seq_for_TEtranscripts
#SBATCH --output=tet3_rna_seq_for_TEtranscripts_%j.out
#SBATCH --error=tet3_rna_seq_for_TEtranscripts_%j.err
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

# Function for mapping with STAR (modified for TEtranscripts)
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
             --outSAMtype BAM SortedByCoordinate Unsorted \
             --outFileNamePrefix "../3_mapped/${sample_name}_" \
             --outFilterMultimapNmax 100 \
             --winAnchorMultimapNmax 200 \
             --outSAMunmapped Within \
             --twopassMode Basic \
             --outWigType bedGraph \
             --outWigStrand Unstranded \
             --quantMode GeneCounts
        printf "$(date '+%Y/%m/%d %H:%M:%S')\n" >> "$time_file"
    done
    echo "STAR mapping completed."
}

# Execute functions
run_star_mapping

echo "All steps completed."
echo "Job finished on $(date)"