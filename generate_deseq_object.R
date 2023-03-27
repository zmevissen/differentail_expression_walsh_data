library(data.table)
library(DESeq2)

dds_file <- "tmp/dds.Rds"


file_names <- list.files("output/star/pass2", full.names = TRUE)



# name the paths to each file with the sample name 
names(file_names) <- sapply(strsplit(basename(file_names), ".", fixed = TRUE), function(x) unlist(x)[[1]])



star_list <- lapply(file_names, fread)
star_counts_raw <- rbindlist(star_list, idcol = "sample")

# remove unmapped reads
star_counts <- star_counts_raw[!startsWith(V1, "N_")]

# column 1:  gene ID
# column 2:  counts for unstranded RNA-seq
# column 3:  counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4:  counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

# TODO CHECK IF STRANDED SEQUENCING

x <- star_counts[, lapply(.SD, sum), .SDcols = c("V2", "V3", "V4"),
                 by = sample]


sample_data <- unique(star_counts[,.(sample, group=gsub("^.*-", "", sample))])
sample_data[, sample_name:= sample]

#What does this do?
# dcast converts long to wide
count_dt <- dcast(star_counts, V1 ~ sample, value.var = "V4")


#What does this do?
#DESeq wants the count_dt column names to be the same as the sample_data row names 
# Here we are ordering count_dt columns by sample_data row order
setcolorder(count_dt, sample_data$sample)



countData <- as.matrix(data.frame(count_dt, row.names = "V1"))


colData <- data.frame(sample_data, row.names = "sample")

#R doesnt like "-" in column names of tables so need to sub "-" for "." in the name
row.names(colData) <- sub("-", ".", row.names(colData))

colData$group <- as.factor(colData$group)
colData$sample_name <- as.factor(colData$sample_name)
# Generate DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ 1)

saveRDS(dds, dds_file)


