#!/bin/bash


####################################  
# Setting up the pipeline script
#################################### 
# Defining the working directory
export ROOT="/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project"
# export ROOT="$(dirname "$(dirname "$(realpath "$0")")")"

export script_dir="$ROOT/code"
# export pipeline_scripts_dir="$script_dir/pipeline_scripts"
export pipeline_scripts_dir="$script_dir" # Ändrar senare

remove_files_scripts_dir="$script_dir/remove_files_scripts"

export data_dir="$ROOT/data"
export raw_data_dir="$ROOT/data/raw"
export results_dir="$ROOT/results"

######################################  
####### Defining parameter values #######
######################################
# export create_gene_catalog=TRUE
export create_gene_catalog=FALSE


# Get the number of logical cores available
cores=$(nproc)
# Set the maximum number of neutral model simulations to be run in parallel
# The value is set by default as the amount of logical cores available, but can be manually altered.
export max_parallel_jobs_neutral_model_simulations=$((cores / 1))
# export max_parallel_jobs_neutral_model_simulations=20


# -------------[ Load Configuration ]-------------
CONFIG_FILE="$script_dir/config.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    echo -e "${GREEN}[INFO]${NC} Loading configuration from ${CONFIG_FILE}"
    source "$CONFIG_FILE"
else
    echo -e "${RED}[ERROR]${NC} Could not find config file at ${CONFIG_FILE}"
    exit 1
fi

# # If chr_specific_recombination_rate=FALSE, the mean recombination rate across all autosomal chromosomes will be used to determine the genetic length of the modeled chromosome:
# export average_recombination_rate=$(echo "scale=4; ($(IFS=+; echo "${chromosome_recombination_rates_cM_per_Mb[*]}")) / ${#chromosome_recombination_rates_cM_per_Mb[@]}" | bc)
# # If chr_specific_recombination_rate=TRUE, the chromosome specific recombination rate of the simulated chromosome will be used: 
# export model_chromosome_recombination_rate="${chromosome_recombination_rates_cM_per_Mb[${chr_simulated}]}"

# Function to handle user interruption
handle_interrupt() {
    echo "Pipeline interrupted. Exiting."
    # Could potentially clean up the files created up until the script termination here
    exit 1
}
# Trap the SIGINT signal (Ctrl+C) and call the handle_interrupt function
trap 'handle_interrupt' SIGINT

# Create a runtime log file
runtime_log="$script_dir/pipeline_runtime.txt"
# Remove the existing runtime logfile if it exists
if [ -e "$runtime_log" ]; then
    rm "$runtime_log"
fi

# Start the timer 
pipeline_start=$(date +%s)
echo "Pipeline Runtimes:" > $runtime_log
# Changing the working directory
cd $ROOT


#���������������������
#� Pipeline Run �
#���������������������

# # Step 0 - Merging multiple coordinate files for wingen format (Optional)
# # Note: 
# run_coordinate_mapping=TRUE
run_coordinate_mapping=FALSE
if [[ "$run_coordinate_mapping" == "TRUE" ]]; then
    step=0
    script_name="1_Coordinate_mapping.Rmd"
    echo "Step $step: Running $script_name"
    source "$pipeline_scripts_dir/$script_name"
    end_step0=$(date +%s)
    runtime_step0=$((end_step0-pipeline_start))
    echo "Step $step: $script_name Runtime: $runtime_step0 seconds" >> $runtime_log
    ((step++))
else
    echo "Coordinate mapping for merging multiple coordinate files set to FALSE (run_coordinate_mapping=${run_coordinate_mapping})."
fi

# # Step 3 - Running the wingen pipeline
# # Note: 
# # /filter_vcf_by_chr_or_gene.sh  chrö  chr8:1000000-5000000  
# input_vcf_file="herring.vcf.gz" # Input path for the vcf file
export wingen_output_dir=$results_dir
mkdir -p $wingen_output_dir
export sample_coordinates_file="${raw_data_dir}/genetic_variation_herring/125_individual_info_list_shared_with_media_coordinates.xlsx" # Output path for the filtered vcf file
export wingen_script_name="2_Herring_Subpopulation_Wingen.Rmd"

export input_vcf_file="${raw_data_dir}/genetic_variation_herring/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz" # Input path for the vcf file

if [[ "$create_gene_catalog" == "TRUE" ]]; then
    export script_name="2_wingen_create_gene_catalog.sh"
    echo "Step $step: Gene Catalog Creation set to TRUE (create_gene_catalog=${create_gene_catalog}). Running $script_name"
    start_step1=$(date +%s)
    source "$pipeline_scripts_dir/$script_name"
    end_step1=$(date +%s)
    runtime_step1=$((end_step1-start_step1))
    echo "Step $step: $script_name Runtime: $runtime_step1 seconds" >> $runtime_log
    ((step++))
else
    echo "Gene Catalog Creation set to False (create_gene_catalog=${create_gene_catalog})."
    # Step 1 - Filtering vcf file for wingen pipeline
    # Note: 
    # /filter_vcf_by_chr_or_gene.sh  chr  chr8:1000000-5000000
    export vcf_filter_mode="chr" # "gene"/"chr" 
    vcf_filter_value="chr1:1-1000000"

    export vcf_filter_mode="gene" # "gene"/"chr" 
    vcf_filter_value="ENSCHAG00000012345"

    export filtered_vcf_file="${raw_data_dir}/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf" # Output path for the filtered vcf file
    compress_vcf_file="true"   # optional 5th arg; not given -> "false"
    step=1
    script_name="filter_vcf_by_chr_or_gene.sh"
    echo "Step $step: Running $script_name"
    start_step1=$(date +%s)
    source "$pipeline_scripts_dir/$script_name" $input_vcf_file $vcf_filter_mode $vcf_filter_value $filtered_vcf_file $compress_vcf_file
    end_step1=$(date +%s)
    runtime_step1=$((end_step1-start_step1))
    echo "Step $step: $script_name Runtime: $runtime_step1 seconds" >> $runtime_log
    ((step++))
    export model_chromosome_physical_length_bp=$model_chromosome_physical_length_bp # Value imported from 1_1_plink_preprocessing_empirical_data.sh
    # export selected_chr_snp_density_mb=$(echo "$selected_chr_preprocessed_snp_density_mb * $num_markers_raw_empirical_dataset_scaling_factor" | bc) # Exempel på numeriska värden att exportera
    # echo "selected_chr_snp_density_mb: $selected_chr_snp_density_mb"

    echo "Step $step: Gene Catalog Creation set to False (create_gene_catalog=${create_gene_catalog}). Running $wingen_script_name"
    source "$pipeline_scripts_dir/$wingen_script_name"
    end_step3=$(date +%s)
    runtime_step3=$((end_step3-end_step2))
    echo "Step $step: $wingen_script_name Runtime: $runtime_step3 seconds" >> $runtime_log
    ((step++))
    # export model_chromosome_physical_length_bp=$model_chromosome_physical_length_bp # Value imported from 1_1_plink_preprocessing_empirical_data.sh
    # export selected_chr_snp_density_mb=$(echo "$selected_chr_preprocessed_snp_density_mb * $num_markers_raw_empirical_dataset_scaling_factor" | bc) # Exempel på numeriska värden att exportera
    # echo "selected_chr_snp_density_mb: $selected_chr_snp_density_mb"
fi


final_step=$(date +%s)
# Calculate total runtime
total_runtime=$((final_step-pipeline_start))
echo "Total Pipeline Runtime: $total_runtime seconds" >> $runtime_log
# Print runtimes
for ((i=1; i<step; i++)); do
    script_name=$(eval echo \$script_name$i)
    runtime=$(eval echo \$runtime_step$i)
    echo "Step $i: $script_name Runtime: $runtime seconds"
done
echo "Total Pipeline Runtime: $total_runtime seconds"

# cd $ROOT
# source "$remove_files_scripts_dir/remove_individual_coverage_and_ROH_pipeline_files.sh"
