library(rtracklayer)
library(GenomicRanges)

#### Settings ####
snp_file <- "~/Documents/Herring_Array/Herring_VCF_SNP_DESIGN_ORIGINAL/List_CHROM_POS_for_n6517_SNP.txt"
gtf_file <- "/Users/jakgo991/Documents/Herring_Array/Genomic_raw_data/Clupea_harengus.Ch_v2.0.2.108.gtf"
gff3_file <- "~/Documents/Herring_Array/Herring_VCF_SNP_DESIGN_ORIGINAL/Clupea_harengus.Ch_v2.0.2.108.gff3"
out_file <- "~/Documents/Herring_Array/Herring_VCF_SNP_DESIGN_ORIGINAL/ALL_SNP_Gene_Annotations.txt"

#### Read SNPs (all chromosomes) ####
snps <- read.delim(snp_file)            # columns: CHROM, POS
snps$CHROM  <- as.character(snps$CHROM)
snps$SNP_ID <- paste(snps$CHROM, snps$POS, sep = ":")

snp_gr <- GRanges(seqnames = snps$CHROM,
                  ranges   = IRanges(start = snps$POS, end = snps$POS),
                  strand   = "*")

#### Read GTF, keep genes only ####
gtf   <- import(gtf_file)
genes <- gtf[gtf$type == "gene"]

# Sanity check: chromosome names must match between SNPs and GTF
shared <- intersect(seqlevels(snp_gr), seqlevels(genes))
message(length(shared), " of ", length(seqlevels(snp_gr)), " SNP chromosomes found in GTF")
if (length(shared) == 0) stop("No matching seqnames: check chromosome naming (e.g. '15' vs 'chr15')")

#### SNPs that fall inside a gene (no distance window) ####
hits <- findOverlaps(snp_gr, genes, ignore.strand = TRUE)

genic <- data.frame(
  SNP_ID       = snps$SNP_ID[queryHits(hits)],
  CHROM        = snps$CHROM[queryHits(hits)],
  POS          = snps$POS[queryHits(hits)],
  gene_id      = genes$gene_id[subjectHits(hits)],
  gene_start   = start(genes)[subjectHits(hits)],
  gene_end     = end(genes)[subjectHits(hits)],
  gene_strand  = as.character(strand(genes))[subjectHits(hits)],
  gene_biotype = if ("gene_biotype" %in% names(mcols(genes))) genes$gene_biotype[subjectHits(hits)] else NA,
  stringsAsFactors = FALSE
)

#### Gene names and descriptions from the Ensembl release 108 GFF3 (no BioMart needed) ####
gff <- import(gff3_file)
gff_genes <- gff[gff$type %in% c("gene", "ncRNA_gene", "pseudogene")]
gff_name <- if ("Name" %in% names(mcols(gff_genes))) as.character(gff_genes$Name) else rep(NA_character_, length(gff_genes))
gff_desc <- if ("description" %in% names(mcols(gff_genes))) as.character(gff_genes$description) else rep(NA_character_, length(gff_genes))
m <- match(genic$gene_id, gff_genes$gene_id)
genic$gene_name   <- gff_name[m]
genic$description <- gff_desc[m]
if ("gene_name" %in% names(mcols(genes))) {   # fall back to GTF gene_name
  gtf_names <- setNames(genes$gene_name, genes$gene_id)
  miss <- is.na(genic$gene_name) | genic$gene_name == ""
  genic$gene_name[miss] <- gtf_names[genic$gene_id[miss]]
}
genic$annotation <- "genic"

#### SNPs in no gene are flagged intergenic ####
inter_idx <- setdiff(seq_along(snp_gr), queryHits(hits))
intergenic <- data.frame(
  SNP_ID = snps$SNP_ID[inter_idx],
  CHROM  = snps$CHROM[inter_idx],
  POS    = snps$POS[inter_idx],
  annotation = ifelse(snps$CHROM[inter_idx] %in% seqlevels(genes), "intergenic", "chrom_not_in_GTF"),
  stringsAsFactors = FALSE
)

#### Combine (one row per SNP-gene pair; one row per intergenic SNP) ####
intergenic[setdiff(names(genic), names(intergenic))] <- NA
ann <- rbind(genic, intergenic[names(genic)])
ann <- ann[order(match(ann$SNP_ID, snps$SNP_ID)), ]

write.table(ann, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

message(nrow(snps), " SNPs: ", length(unique(genic$SNP_ID)), " genic, ", sum(intergenic$annotation == "intergenic"), " intergenic, ",
        sum(intergenic$annotation == "chrom_not_in_GTF"), " on chromosomes/scaffolds absent from the GTF")
