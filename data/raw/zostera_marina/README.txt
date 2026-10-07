### General information
 
Authors: Stefanie R. Ries, Ellika Faust, Kerstin Johannesson, Per R. Jonsson, Per-Olav Moksnes, Ricardo T. Pereyra, Marlene Jahnke
Contact e-mail:stefanie.ries@gu.se
DOI: 10.17044/scilifelab.24610356
Licence: CC BY 4.0
This readme file was last updates: 21-12-2023

Please cite as: Ries, S. R., Faust, E., Johannesson, K., Jonsson, P. R., Moksnes, P.-O., Pereyra, R. T., et al. (2023). Genetic structure and diversity of the seagrass Zostera marina along a steep environmental gradient, with implications for genetic monitoring. Frontiers in Climate 5. doi: 10.3389/fclim.2023.1303337. 
SciLifeLab. Dataset. https://doi.org/10.17044/scilifelab.24610356

### Dataset description

This repository includes data files and scripts used in the publication "Genetic structure and diversity of the seagrass Zostera marina along a steep environmental gradient, with implications for genetic monitoring" in the special Issue "Seagrasses and Climate Change: Challenges and Opportunities" in the journal Frontiers Climate 2023.

This publication assessed 15 sampling sites along the Swedish coast and examined clonality, genetic variation and population structure using 2,138 Single Nucleotide Polymorphisms (SNPs) derived from 2b-RAD sequencing. Data in this repository includes final filtered vcf files, the bioinformatic pipeline, the laboratory protocol, sampling site information, and scripts for the data analysis.

Data in this repository includes: 
final filtered vcf files for all, neutral and outlier loci.
Raw sequence data is archived in GenBank's SRA: BioProject ID PRJNA1043091 and BioSample 
IDs SAMN38331216 - SAMN38331474. 

### Available variables 

File list and usage notes:

'01.Read_convert.rmd': R script for converting vcf file into genind and genlight files.

'02.poppr_clones.rmd': R script for clonal analysis.

'03.Basic_stat.rmd': R script for basic statistics of the data.

'04.mlg_diversity.rmd': R script for genetic diversity statistics for the data containung no clones.

'05.maps.rmd': R script for making maps.

'06.PCA.rmd': R script for making principal component analysis.

'07.new_LEA.rmd': R script for doing structure plots.

'08.Sea_Distance.rmd': R script for calculating sea distance between sampling sites.

'09.genetic_div_map.rmd': R script for making interpolated heatmaps with the genetic diversity values.

'pop_info.csv': sampling site information with coordinates, site names and abbreviations

'zostera_monitoring_230504_5771_SNPs_all_with_replicates.vcf.zip': vcf file with 259 samples including clones and replicates

'zostera_monitoring_230504_2138_SNPs_all_no_replicates.vcf.zip': vcf file with 259 samples including clones and no replicates

'zostera_monitoring_230504_2138_SNPs_81_MLLs.vcf.zip': vcf file with 81 clonal lineages no replicates

'Protocol_for_Illumina_2b-RAD_Ries_2023.pdf': adapted protocol for the 2bRAD library preperation

'2bRAD_Zostera_marina_pipeline_2023-09-25.sh': bioinformatic pipeline