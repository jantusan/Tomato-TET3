#!/bin/bash
#SBATCH --job-name=TEcount_analysis
#SBATCH --output=TEcount_analysis_%j.out
#SBATCH --error=TEcount_analysis_%j.err
#SBATCH --time=12:00:00       # Adjust time as needed
#SBATCH --mem-per-cpu=8000    # Memory per CPU core
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4     # Number of CPU cores per task

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
bam_dir="$dir/data/3_mapped"
output_dir="$dir/TEcount_output"
annotation="/path/to/genome/v4_Mt_Pt_TET3/ITAG4.1_gene_models_Mt_Pt_TET3.gtf"
te_annotation="/path/to/genome/ITAG4.0/TEtranscripts/solLyc_ITAG4.0_RepeatModeler_repeats_light_TE.gtf"

# Create output directory
mkdir -p "$output_dir"

# Function to run TEcount on each BAM file
run_TEcount() {
    echo "Running TEcount..."
    cd "$bam_dir"

    # List of BAM files
    bam_files=(
        "WT_1_Aligned.out.bam"
        "WT_2_Aligned.out.bam"
        "TET3_1_Aligned.out.bam"
        "TET3_2_Aligned.out.bam"
        "NTS_1_Aligned.out.bam"
        "NTS_2_Aligned.out.bam"
    )

    # Loop over BAM files and run TEcount
    for bam_file in "${bam_files[@]}"; do
        sample_name="${bam_file%_Aligned.out.bam}"
        output_prefix="$output_dir/${sample_name}"

        # Check if BAM file exists
        if [ ! -f "$bam_file" ]; then
            echo "Error: BAM file $bam_file does not exist."
            exit 1
        fi

        # Run TEcount
        TEcount --format BAM \
                --stranded no \
                --mode multi \
                -b "$bam_file" \
                --GTF "$annotation" \
                --TE "$te_annotation" \
                --project "$output_prefix" \
                --sortByPos

        echo "TEcount completed for $sample_name."
    done
    echo "All TEcount analyses completed."
}

# Execute function
run_TEcount

echo "Job finished on $(date)"