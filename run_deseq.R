library(DESeq2)
library(data.table)
library(ggplot2)

# dds_file <- snakemake@input[["dds"]]
# inj_res_file <- snakemake@output[["injury_results"]]
# group1_res_file <- snakemake@output[["group1_results"]]
# group2_res_file <- snakemake@output[["group2_results"]]
# threads <- snakemake@threads

# dev
dds_file <- "tmp/dds.Rds"
threads <- 8

BiocParallel::register(
  BiocParallel::MulticoreParam(workers = threads))


dds <- readRDS(dds_file)


#GENERAL FILTERING
#Getting normalized read count / library size
dds <- estimateSizeFactors(dds)

# remove genes with low expression
keep_genes <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 3


#FILTERING FOR MICROGLIA

# filter out Bulk samples here 
cd45_keep_samples <- colnames(dds)[grepl("CD45", colnames(dds))]
cd45_dds_filtered <- dds[keep_genes, cd45_keep_samples]
design(cd45_dds_filtered) <- ~ sample_name
cd45_dds_filtered <- DESeq(cd45_dds_filtered, test = "LRT", reduced = ~ 1)


#KEEPING BULK SAMPLES
all_dds_filtered <- dds[keep_genes]
colnames(all_dds_filtered)
design(all_dds_filtered) <- ~ group


all_dds_filtered <- DESeq(all_dds_filtered)
colnames(all_dds_filtered)
bulk_vs_cd45_results <- data.table(results(all_dds_filtered,
                                   name = "group_CD45_vs_Bulk",
                                   tidy = TRUE,
                                   lfcThreshold = log(1.5, 2),
                                   alpha = 0.1))
all_stabilised_dds <- varianceStabilizingTransformation(all_dds_filtered, blind=TRUE)
colnames(all_stabilised_dds)
pca <- prcomp(t(assay(all_stabilised_dds)), center = TRUE)

pca$x[,1]

plotPCA(all_stabilised_dds, intgroup = c("group"))


#DBSCAN 
#fuzzy c-means
