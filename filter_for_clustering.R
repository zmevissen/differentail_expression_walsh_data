library(data.table)
library("dbscan")

install.packages("ClusterR")
install.packages("cluster")

# Loading package
library(ClusterR)
library(cluster)



file_names <- list.files("output/star/pass2", full.names = TRUE)



# name the paths to each file with the sample name 
names(file_names) <- sapply(strsplit(basename(file_names), ".", fixed = TRUE), function(x) unlist(x)[[1]])



star_list <- lapply(file_names, fread)
star_counts_raw <- rbindlist(star_list, idcol = "sample")

# remove unmapped reads
star_counts <- star_counts_raw[!startsWith(V1, "N_")]

count_dt <- dcast(star_counts, V1 ~ sample, value.var = "V4")
countData <- as.matrix(data.frame(count_dt, row.names = "V1"))

countData <- countData[,colnames(countData)[grepl("CD45", colnames(dds))]]
countData <- countData[0:1000,]

kmeansData <- kmeans(countData, centers = 50, nstart = 20)

kmeansData$cluster

#VISUALISATION
y_kmeans <- kmeansData$cluster
clusplot(countData,
         y_kmeans,
         lines = 0,
         shade = TRUE,
         color = TRUE,
         labels = 2,
         plotchar = FALSE,
         span = TRUE)




#db <- dbscan(countData[0:1000,], eps = 0.4, minPts = 4)
#pairs(countData[0:1000,], col = db$cluster + 1L)

