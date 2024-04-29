#install.packages("BiocManager")

#BiocManager::install("SC3")
#BiocManager::install("scater")


###Prep count matrix and metadata
plate<-gsub('^(.*?_.*?_.*?)_.*', "\\1",  colnames(counts))

plate.md<-cbind(colnames(counts),plate)
plate.md<-data.frame(plate.md)
rownames(plate.md)<-unique(plate.md$V1)

celltype<-ifelse(grepl("d0",plate), "d0",
                 ifelse(grepl("d4",plate), 'd4', 
                        ifelse(grepl("d8",plate), 'd8',
                               'd12')))

celltype.md<-cbind(colnames(counts),celltype)
celltype.md<-data.frame(celltype.md)
rownames(celltype.md)<-celltype.md$V1
counts<-as.matrix(counts)



browseVignettes("SC3")

library(SingleCellExperiment)
library(SC3)
library(scater)
install.packages("mclust")
library(mclust)
library(dplyr)


d0_matches <- grep("d0", colnames(mouseTags), value = TRUE)
d0_matches_clean <- sub(".*d0.*", "d0", d0_matches)

d4_matches <- grep("d4", colnames(mouseTags), value = TRUE)
d4_matches_clean <- sub(".*d4.*", "d4", d4_matches)

d8_matches <- grep("d8", colnames(mouseTags), value = TRUE)
d8_matches_clean <- sub(".*d8.*", "d8", d8_matches)

d12_matches <- grep("d12", colnames(mouseTags), value = TRUE)
d12_matches_clean <- sub(".*d12.*", "d12", d12_matches)

# Combine all matches into a single vector
cellType <- c(d0_matches_clean, d4_matches_clean, d8_matches_clean, d12_matches_clean)
c(d4_matches_clean, d8_matches_clean, d12_matches_clean, d0_matches_clean)
cellType <- factor(cellType)

# Manually assign values to types
types <- c("d0", "d4", "d8", "d12")

## Assign each cell type a color
scols <- c("red","blue","green","brown","pink","purple","darkgreen","grey")
cols <- rep(NA,length(cellType))
for (i in 1:length(cols)){
  cols[i] <- scols[which(types==cellType[i])]
}
cols

# Check for constant or zero variance columns
constant_zero_variance_cols <- sapply(mouseTags, function(x) length(unique(x)) == 1)

# Remove constant or zero variance columns
data_clean <- mouseTags[, !constant_zero_variance_cols]

# Perform PCA on the preprocessed data
#pca_result <- prcomp(data_clean, center = TRUE, scale. = TRUE)


# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data_clean),
    logcounts = log2(as.matrix(data_clean) + 1)
  ),
  #colData = cols
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- runPCA(sce)
plotPCA(
  sce
)

#plotPCA(sce, colour_by = "cell_type1")
sce <- sc3(sce, ks = 2:5, biology = TRUE)
sc3_plot_consensus(sce, k = 5)
sc3_plot_markers(sce, k = 3)


sc3_export_results_xls(sce)

########################################################################
humanTags <- read.delim("counts.human.tsv")

# Get the number ofhdata)
num_rows <- nrow(humanTags)

# Calculate the midpoint
midpoint <- floor(num_rows / 2)

# Split the data into two halves
first_half <- humanTags[1:midpoint, ]
second_half <- humanTags[(midpoint + 1):num_rows, ]
# Shuffle the data
shuffled_data <- humanTags[sample(nrow(humanTags)), ]

# Split the shuffled data into two halves
first_half <- shuffled_data[1:midpoint, ]
second_half <- shuffled_data[(midpoint + 1):num_rows, ]

filtered_data <- humanTags %>%
  select(-starts_with("M"))

# Assuming 'data' is your original data frame and 'sd_noise' is the standard deviation of the Gaussian noise
# 'mean_noise' is the mean, usually 0 for Gaussian noise

# Generate Gaussian noise
noise <- rnorm(nrow(humanTags))

# Add noise to your data
noisy_data <- humanTags + noise

info <- colnames(humanTags)
# Filter on the celltype prefixs
# Match Fibroblasts
# Match Fibroblasts
fibroblasts_matches <- grep("fibroblasts", tolower(colnames(humanTags)), value = TRUE)
fibroblasts_matches_clean <- sub("_.*", "", fibroblasts_matches)

# Match HEK293T
hek293t_matches <- grep("hek293t", tolower(colnames(humanTags)), value = TRUE)
hek293t_matches_clean <- sub("_.*", "", hek293t_matches)

# Match MCF7
mcf7_matches <- grep("mcf7", tolower(colnames(humanTags)), value = TRUE)
mcf7_matches_clean <- sub("_.*", "", mcf7_matches)

# Combine all matches into a single vector
cellType <- c(fibroblasts_matches_clean, hek293t_matches_clean, mcf7_matches_clean)
cellType <- factor(cellType)

# Manually assign values to types
types <- c("fibroblasts", "hek293t", "mcf7")

## Assign each cell type a color
scols <- c("red","blue","green","brown","pink","purple","darkgreen","grey")
cols <- rep(NA,length(cellType))
for (i in 1:length(cols)){
  cols[i] <- scols[which(types==cellType[i])]
}

# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(humanTags),
    logcounts = log2(as.matrix(humanTags) + 1)
  ),
  colData = cols
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- runPCA(sce)

# Think about changing parameters for each step2
sce <- sc3(sce, ks = 2:4, biology = TRUE)
adjustedRandIndex(colData(sce)$sc3_3_clusters, cols) # 1

twotypessce <- sce
no_svm_labels_twotypes <- colData(twotypessce)$sc3_3_clusters # 1
adjustedRandIndex(no_svm_labels_twotypes, cols) # [1] 0.996414
sc3_plot_consensus(
  twotypessce, k = 3, 
  show_pdata = c(
    "X", 
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)

first_half_sce = sce
sc3_plot_consensus(
  twotypessce, k = 3, 
  show_pdata = c(
    "X", 
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)
no_svm_labels_firsthalf <- colData(first_half_sce)$sc3_3_clusters # 1
adjustedRandIndex(no_svm_labels_firsthalf, cols)


new_sce = sce
no_svm_labels <- colData(sce)$sc3_3_clusters
adjustedRandIndex(no_svm_labels, cols) # 1
# Assuming 'sce' is your SingleCellExperiment object

# Step 1: Extract the count data from the SingleCellExperiment object
count_data <- assays(sce)$sc3_3_clusters
sd(count_data)

sc3_plot_markers(sce, k = 3)

sc3_plot_consensus(
  sce, k = 3, 
  show_pdata = c(
    "X", 
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)

sc3_plot_de_genes(sce, k = 3) # needs to be 3 for human
sc3_plot_consensus(sce, k = 3)
sc3_export_results_xls(sce)

k <- sc3_estimate_k(sce)
sc3_estimate_k(sce)
str(metadata(sce)$sc3)
sc3_plot_cluster_stability(sce, k = 3)
plotPCA(
  sce,
  colour_by = "X"
)
