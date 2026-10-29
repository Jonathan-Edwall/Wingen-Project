#!/bin/bash -l

# Start the timer 
script_start=$(date +%s)

####################################  
# Defining the working directory
#################################### 

# ROOT=/ROOT/jonathan
conda_setup_script_path="/proj/uppmax2025-2-125/nobackup/jonathan/anaconda3/etc/profile.d/conda.sh"  # Variable defined in config.sh
source $conda_setup_script_path  # Source Conda initialization script
# Activate the conda environment
conda activate wingen_env


# export ROOT="/proj/uppmax2025-2-125/nobackup/jonathan/Wingen-Project" # Variable defined in the main pipeline script!

# cd $ROOT
# pipeline_script_dir=$ROOT/code/pipeline_scripts

# export script_dir="${ROOT}/code" # Variable defined in the main pipeline script!
# export pipeline_scripts_dir="${script_dir}" # Variable defined in the main pipeline script!

# export results_dir="${ROOT}/results" # Variable defined in the main pipeline script!
# export raw_data_dir="${ROOT}/data/raw" # Variable defined in the main pipeline script!
export wingen_script_name="2_Herring_Wingen.Rmd"

#################################### 
# Defining Simulation parameters
#################################### 
# data_dir=$ROOT/data # Variable defined in the main pipeline script!
# vcf_file_dir_gene_catalog="${raw_data_dir}/genetic_variation_herring/temp"
# mkdir -p $vcf_file_dir_gene_catalog

output_dir_gene_catalog="${ROOT}/results_gene_catalog"
# output_dir_gene_catalog="${results_dir}/gene_catalog"
mkdir -p $output_dir_gene_catalog
export wingen_output_dir=$output_dir_gene_catalog # Variable defined in the main pipeline script!

# export chr_simulated="chr3" # Defined in the main script
n_genes=1 # Defined in the main script


############
# Temp
############
# Input path for the vcf file
input_vcf_file="${raw_data_dir}/genetic_variation_herring/head_chr1_1-1e6_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz" # Variable defined in the main pipeline script!
# input_vcf_file="${raw_data_dir}/genetic_variation_herring/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz" # Variable defined in the main pipeline script!

vcf_file_dir_gene_catalog="${raw_data_dir}/genetic_variation_herring/temp"

mkdir -p $vcf_file_dir_gene_catalog

export calc_gen_diversity_stats="TRUE"
# Latitude Coordinates
export LAT_MIN="30"
export LAT_MAX="72"
# Longitude Coordinates
export LON_MIN="-25"
export LON_MAX="60" # Ural mountains in Russia roughly

# export GRD_DISAGG_FACTOR="1"
export GRD_DISAGG_FACTOR="7"
export RESOLUTION_KM="5"
# export MAX_INTERP_DIST_MODE
# export MAX_INTERP_DIST_KM
# export MAX_INTERP_DIST_SCALE

# export NMAX_NEIGHBOURS
# export USE_WEIGHTED_KRIGING

# export MIN_SAMPLE_COUNT_MODE
# export MIN_SAMPLE_COUNT
# export HEATMAP_COLOR_PALETTE

# gene_to_filter="ENSCHAG00000012345"
# filtered_vcf_file="${vcf_file_dir_gene_catalog}/ENSCHAG00000012345_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"

export sample_coordinates_file="${raw_data_dir}/genetic_variation_herring/125_individual_info_list_shared_with_media_coordinates.xlsx"

export study_species="herring_sentieon_125ind"
# 'gene_list_source' = "array" (default) | "txt" | "bed". Variable defined in the main pipeline script!
#     "array" - use the hardcoded 'genes' array immediately below, as-is (unchanged
#               default behaviour).
#     "txt"   - read gene IDs from gene_list_file, one per line, into 'genes'. This is
#               exactly the format snpeff_vcf_detect_unique_genes.sh writes with
#               format=txt, so its output can be pointed at directly as gene_list_file.
#     "bed"   - read gene IDs from gene_list_file's 4th (gene_id) column into 'genes'.
#               This is exactly the format snpeff_vcf_detect_unique_genes.sh writes with
#               format=bed (chrom, start, end, gene_id), or any hand-built BED of genes.
#               A header row (e.g. "chrom start end gene_id" from a hand-made or
#               downloaded BED) is skipped automatically: real BED data rows always have
#               a numeric start in column 2, so a row that doesn't is dropped rather than
#               turned into a bogus "gene" - snpeff_vcf_detect_unique_genes.sh never
#               writes a header itself, so this only matters for BED files from elsewhere.
#
# 'gene_list_file' = required only if gene_list_source is "txt" or "bed". Path to that file.
#     Note this is independent of gene_bed_file below (used for coordinate lookup when
#     vcf_annotation_format="bed") - they're read separately, so nothing forces them to
#     match, but you can point both at the same file if that's the one BED you have.

gene_list_source="${gene_list_source:-array}" # "array" / "txt" / "bed"

# gene_list_source="txt" # temp 2026-09-14

gene_list_file="${gene_list_file:-}" # Variable defined in the main pipeline script!
# gene_list_file="${raw_data_dir}/genetic_variation_herring/unique_genes.txt" # temp 2026-09-14

genes=("ENSCHAG00000000003" "ENSCHAG00000000004")   # only used when gene_list_source="array"

case "$gene_list_source" in
    array)
        : # 'genes' array above is used as-is
        ;;
    txt)
        if [ ! -f "$gene_list_file" ]; then
            echo "ERROR: gene_list_source='txt' requires gene_list_file to point to an existing file (got '$gene_list_file')" >&2
            exit 1
        fi
        # Extracting the gene catalog genes defined in the txt-file. grep -v drops any blank lines rather than turning them into empty "genes"
        mapfile -t genes < <(grep -v '^[[:space:]]*$' "$gene_list_file")
        ;;
    bed)
        if [ ! -f "$gene_list_file" ]; then
            echo "ERROR: gene_list_source='bed' requires gene_list_file to point to an existing BED file (got '$gene_list_file')" >&2
            exit 1
        fi
        # Extracting the bed-file genes. $2 ~ /^[0-9]+$/ requires a numeric start column,
        # which drops a header row (its column 2 would be text like "start", not a
        # number) instead of reading it in as a gene; grep -v drops any blank lines too.
        mapfile -t genes < <(awk -F'\t' '$2 ~ /^[0-9]+$/ {print $4}' "$gene_list_file" | grep -v '^[[:space:]]*$')
        ;;
    *)
        echo "ERROR: gene_list_source must be 'array', 'txt', or 'bed' (got '$gene_list_source')" >&2
        exit 1
        ;;
esac

# 'gene_list_limit' = optional. If set, caps the gene list (from whichever source above)
# to its first N genes - handy for a quick toy-example test run without needing a
# separate small file. Leave unset/empty to process the full list.
gene_list_limit="${gene_list_limit:-}"
# gene_list_limit="2" # temp 2026-09-14
if [ -n "$gene_list_limit" ]; then
    if ! [[ "$gene_list_limit" =~ ^[0-9]+$ ]]; then
        echo "ERROR: gene_list_limit must be a positive integer if set (got '$gene_list_limit')" >&2
        exit 1
    fi
    genes=("${genes[@]:0:$gene_list_limit}")
    echo "gene_list_limit=$gene_list_limit set - truncated gene list to the first $gene_list_limit gene(s) for this run"
fi

# 'vcf_annotation_format' = "snpeff" (default) or "bed". Variable defined in the main pipeline script!
#     "snpeff" - filter_vcf_by_chr_or_gene.sh is called in "gene" mode, gene ID as-is.
#                Despite the name, this no longer strictly requires input_vcf_file to be
#                SnpEff-annotated: gene mode itself now checks the header and either
#                matches SnpEff's ANN (anchored) or falls back to a plain INFO-column
#                substring search for some other tool's gene tagging (with a warning,
#                since that fallback has no known field structure to anchor to). So read
#                this value as "resolve the gene via gene mode" rather than "this file is
#                definitely SnpEff".
#     "bed"    - input_vcf_file has NO usable in-file gene tagging at all. gene_bed_file
#                (below) is used to look up each gene's coordinates, and
#                filter_vcf_by_chr_or_gene.sh is called in "chr" mode with that region
#                instead - bypassing INFO entirely regardless of what's in it.
#
# 'gene_bed_file' = required only if vcf_annotation_format="bed". A BED file
#     (chrom, start, end, gene_id, ...) mapping each gene to its genomic span, in
#     input_vcf_file's coordinate system. One row per gene is assumed.
vcf_annotation_format="${vcf_annotation_format:-snpeff}"

# ---- Validate configuration up front, before any filtering work starts ----
if [[ "$vcf_annotation_format" != "snpeff" && "$vcf_annotation_format" != "bed" ]]; then
    echo "ERROR: vcf_annotation_format must be 'snpeff' or 'bed' (got '$vcf_annotation_format')" >&2
    exit 1
fi

if [ "$vcf_annotation_format" == "bed" ] && [ ! -f "${gene_bed_file:-}" ]; then
    echo "ERROR: vcf_annotation_format='bed' requires gene_bed_file to point to an existing BED file (got '${gene_bed_file:-<unset>}')" >&2
    exit 1
fi

# 'gene_scope' = "all" (default) | "body" | "exonic" | "cds" (variable is defined in themain pipeline script!)
# This variable determines if the vcf filtering in the gene catalog should be narrowed down to only include gene body ("body"), "exonic" or "cds".
# for narrowed down gene filtering, the vcf file needs to have ANN annotation in the INFO column for filtering and it is assumed that the input vcf file has been processed through snpEff.
# see filter_vcf_by_scope.sh for why a coordinate-based narrowing can't do this correctly for exonic/cds- exons and CDS
# (disjoint intervals with introns as gaps between them, and a single [start,end) interval can't represent that.
# Note: This 'gene_scope' variable is only used when vcf_annotation_format="snpeff" (gene mode/ANN-based).

gene_scope="${gene_scope:-all}"
# gene_scope="exonic" # temp 2026-09-14

if [[ "$gene_scope" != "all" && "$gene_scope" != "body" && "$gene_scope" != "exonic" && "$gene_scope" != "cds" ]]; then
    echo "ERROR: gene_scope must be 'all', 'body', 'exonic', or 'cds' (got '$gene_scope')" >&2
    exit 1
fi

if [ "$gene_scope" != "all" ] && [ "$vcf_annotation_format" == "bed" ]; then
    echo "ERROR: gene_scope='$gene_scope' requires vcf_annotation_format='snpeff' - the 'bed'" >&2
    echo "       pathway has no ANN in the vcf to check scope against. Use gene_scope='all'" >&2
    echo "       (the default) with vcf_annotation_format='bed', or switch to 'snpeff'." >&2
    exit 1
fi

if [ ${#genes[@]} -eq 0 ]; then
    echo "ERROR: No genes to process. Define the 'genes' array in the main pipeline script, e.g.:" >&2
    echo "         genes=(\"ENSCHAG00000012345\" \"ENSCHAG00000067890\")" >&2
    echo "       ...or set gene_list_source=\"txt\"/\"bed\" with gene_list_file pointing at a list." >&2
    exit 1
fi



n_genes=${#genes[@]}

#����������������������������������������������������������������������������
# Function: 
# Gene Catalog Creation in WINGEN
#
###Input: 
# 
###Output: 
#����������������������������������������������������������������������������
# cd $output_dir_gene_catalog

# Number of parallel jobs to run at a time
max_parallel_jobs_gene_catalog_creation=2 # Variable defined in the main script!

rmd_script_full_path="${pipeline_scripts_dir}/${wingen_script_name}"

# Look up the genomic region (chr:start-end) for a gene from a BED file.
# BED is 0-based half-open; bcftools regions are 1-based closed, so start gets +1.
# Assumes one row per gene (a gene-level BED, not one row per exon/transcript).
get_gene_region_from_bed() {
    local gene_id=$1
    local bed_file=$2
    local bed_line chrom start end

    bed_line=$(awk -F'\t' -v gene="$gene_id" '$4 == gene {print; exit}' "$bed_file")

    if [ -z "$bed_line" ]; then
        echo "ERROR: Gene '$gene_id' not found in BED file: $bed_file" >&2
        exit 1
    fi

    chrom=$(cut -f1 <<< "$bed_line")
    start=$(cut -f2 <<< "$bed_line")
    end=$(cut -f3 <<< "$bed_line")

    echo "${chrom}:$((start + 1))-${end}"
}

# Function to run a single simulation
create_gene_vcf_file() {
    local gene_to_filter=$1
    local filtered_vcf_file=$2 # Path to the new vcf file

    # export vcf_filter_mode="gene" # "gene"/"chr"
    compress_vcf_file="true"   # optional 5th arg; not given -> "false"
    export input_vcf_file="$input_vcf_file" # Variable defined in the main pipeline script!

    if [ "$vcf_annotation_format" == "bed" ]; then
        export vcf_filter_mode="chr"
        filter_value=$(get_gene_region_from_bed "$gene_to_filter" "$gene_bed_file") # chromosome/region, e.g. chr8:1000000-5000000
        echo "Gene $gene_to_filter -> region $filter_value (from $gene_bed_file)"
    else
        export vcf_filter_mode="gene"
        filter_value=$gene_to_filter # gene, e.g. "ENSCHAG00000000001"
    fi
    script_name="filter_vcf_by_chr_or_gene.sh"
    echo "Running $script_name for gene=${gene_to_filter}"
    start_filter=$(date +%s)
    # Arg 1 (INPUT): $input_vcf_file = Input path of the vcf file
    # Arg 2 (MODE): $vcf_filter_mode = Filter mode. I.e: "gene"/"chr"
    # Arg 3 (VALUE): $filter_value = resolved value to filter on.
    #                * For "gene" mode this is the gene ID itself
    #                * For "chr" mode this is the region string looked up from the gene_bed_file 
    # Arg 4 (OUTPUT): $filtered_vcf_file = Output path for the filtered vcf file
    # Arg 5 (COMPRESS), Optional: $compress_vcf_file = Boolean value determining if the output vcf file should be binary or not (bcf or vcf)

    source "$pipeline_scripts_dir/$script_name" $input_vcf_file $vcf_filter_mode $filter_value $filtered_vcf_file $compress_vcf_file
    end_filter=$(date +%s)
    echo "Filtering runtime for gene $gene_to_filter: $((end_filter - start_filter)) seconds"


}


# Function to run a single simulation
create_gene_catalog() {
    local counter=$1
    local gene=$2
    local knit_document_check=$3 # Variable that controls the knitting of the .rmd file
    export sample_coordinates_file=$sample_coordinates_file # Variable defined in the main pipeline script! 

    # export chr_specific_recombination_rate="$chr_specific_recombination_rate" # Variable defined in the main pipeline script! 
    gene_to_filter=$gene
    export filtered_vcf_file="${vcf_file_dir_gene_catalog}/${gene_to_filter}_filtered.vcf.gz" # Needed for the wingen script

    # Ideally output would be this, where the input filename is used and the filtered gene is used as a prefix in the filtered_vcf_file:
    # gene_to_filter="ENSCHAG00000012345"
    # input_vcf_file="${raw_data_dir}/genetic_variation_herring/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz" # Input path for the vcf file
    # filtered_vcf_file="${vcf_file_dir_gene_catalog}/ENSCHAG00000012345_herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.biallelic.vcf.gz"

    create_gene_vcf_file $gene_to_filter $filtered_vcf_file
    
    # Additional filtering step (filter_vcf_by_scope.sh) used if 'gene_scope' = body/exonic/cds, to ensure that only the request scope of the gene is kept.
    # Note: This part is deliberately not wrapped in "||". By letting a failure here propagate via set -e (inherited from filter_vcf_by_chr_or_gene.sh
    # having already been sourced above), lets the existing wait -n/ failed_job_count tracking in the main loop correctly see if this gene as failed;
    # catching it locally here would "silently swallow" that instead.
    export gene_scope="$gene_scope"
    if [ "$gene_scope" != "all" ]; then
        scope_vcf_file="${vcf_file_dir_gene_catalog}/${gene_to_filter}_${gene_scope}.vcf.gz"
        echo "Narrowing gene=${gene_to_filter} to scope=${gene_scope}"
        bash "$pipeline_scripts_dir/filter_vcf_by_scope.sh" "$filtered_vcf_file" "$gene_to_filter" "$gene_scope" "$scope_vcf_file" true
        export filtered_vcf_file="$scope_vcf_file"
    fi

    
    export vcf_filter_value="$gene_to_filter"
    # export input_vcf_file="${filtered_vcf_file}" # Needed for the wingen script

    # Gene-specific intermediates dir: without this, concurrent renders of the same
    # .Rmd clobber each other's .knit.md - required now that max_parallel_jobs_gene_catalog_creation > 1.
    gene_render_dir="${output_dir_gene_catalog}/.render_${gene_to_filter}"
    mkdir -p "$gene_render_dir"
    # Rscript -e "rmarkdown::render('$rmd_script_full_path')"
    # /usr/bin/time -v gives per-gene wall time / peak memory / CPU time - sacct/jobstats
    # only see the whole job, not individual genes within it. Written to
    # output_dir_gene_catalog (not gene_render_dir) since gene_render_dir gets deleted
    # by the cleanup step below.
    /usr/bin/time -v Rscript -e "rmarkdown::render('$rmd_script_full_path', run_pandoc=FALSE, intermediates_dir='$gene_render_dir')" 2> "${output_dir_gene_catalog}/resource_usage_${gene_to_filter}.txt" # Run the .rmd script without knitting!

    # Rscript -e "rmarkdown::render('$rmd_script_full_path')"
    # Rscript -e "rmarkdown::render('$rmd_script_full_path', run_pandoc=FALSE, intermediates_dir='$gene_render_dir')" # Run the .rmd script without knitting!

# e    if [ "$knit_document_check" -eq 1 ]; then
#         Rscript -e "rmarkdown::render('$rmd_script_full_path')"
#     else
#         Rscript -e "rmarkdown::render('$rmd_script_full_path', run_pandoc=FALSE)" # Run the .rmd script without knitting!
#     fi         

    echo "Gene $counter of $n_genes genes completed"
}

# Running through every gene in $genes, with job control
failed_job_count=0

for ((i=0; i<n_genes; i++))
do
    counter=$((i + 1))
    gene="${genes[$i]}"

    if [ "$counter" -eq 1 ]; then
        knit_document_check=1  # Only knit for the first gene
    else
        knit_document_check=0  # Just run the script for all other cases
    fi

    create_gene_catalog $counter "$gene" $knit_document_check &  # Run the job in the background
    
    # If the number of background jobs reaches the limit, wait for one to finish
    while [ $(jobs -r | wc -l) -ge $max_parallel_jobs_gene_catalog_creation ]; do
        wait -n || failed_job_count=$((failed_job_count + 1))
    done
done

# Wait for all remaining jobs to finish, still counting any failures
while [ "$(jobs -r | wc -l)" -gt 0 ]; do
    wait -n || failed_job_count=$((failed_job_count + 1))
done

if [ "$failed_job_count" -gt 0 ]; then
    echo "ERROR: $failed_job_count of $n_genes gene job(s) failed - see FAILED lines above" >&2
    exit 1
fi

# Removing the generated per-gene render directories (.knit.md and other rmarkdown intermediates)
rm -rf "${output_dir_gene_catalog}"/.render_*

# Ending the timer 
script_end=$(date +%s)
# Calculating the script_runtime of the script
script_runtime=$((script_end-script_start))

echo "Gene Catalog Created"
echo "The outputfiles are stored in: $output_dir_gene_catalog"
echo "Runtime: $script_runtime seconds"
