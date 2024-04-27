# Name of the package you want to install
package_name <- "VCCRI/CIDR"

# Check if the package is already installed
if (!require(package_name, character.only = TRUE)) {
  # If not installed, install the package
  devtools::install_github("VCCRI/CIDR")
}

library(cidr)

## Read in data
humanTags <- read.delim("counts.human.tsv")
# -------------------------------------------------
## Read in annotation
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

##  Standard principal component analysis using prcomp
priorTPM <- 1
humanTags10 <- humanTags[rowSums(humanTags)>10,]
humanTags10_lcpm <- log2(t(t(humanTags10)/colSums(humanTags10))*1000000+priorTPM)
pca <- prcomp(t(humanTags10_lcpm))
plot.new()
legend("center", legend=types, col=scols, pch=1)
plot(pca$x[,c(1,2)], col=cols, pch=1,
     xlab="PC1", ylab="PC2",
     main="Principal Component Analysis (prcomp)")

#### CIDR ######
################
scHuman <- scDataConstructor(as.matrix(humanTags))
scHuman <- determineDropoutCandidates(scHuman)
scHuman <- wThreshold(scHuman)
scHuman <- scDissim(scHuman)
scHuman <- scPCA(scHuman)
scHuman <- nPC(scHuman)
nCluster(scHuman)
scHuman <- scCluster(scHuman)

## Two dimensional visualization plot output by CIDR
## Different colors denote the cell types annotated by the human brain single-cell RNA-Seq study
## Different plotting symbols denote the clusters output by CIDR
plot.new()
legend("center", legend=types, col=scols, pch=15)
plot(scHuman@PC[,c(1,2)], col=cols, pch=scHuman@clusters,
     main="CIDR", xlab="PC1", ylab="PC2")
scHuman@dscThreshold
## Use Adjusted Rand Index to measure the accuracy of CIDR clustering
ARI_CIDR <- adjustedRandIndex(scHuman@clusters,cols)
ARI_CIDR
## 0.8977449

## This section shows how to alter the parameters of CIDR
scHuman@nPC
## The default number of principal coordinates used in clustering is 4
## Use 6 instead of 4 principal coordinates in clustering
scHuman <- scCluster(scHuman, nPC=6)

plot(scHuman@PC[,c(1,2)],
     col=cols, pch=scHuman@clusters, main="CIDR (nPC=6)",
     xlab="PC1", ylab="PC2")
ARI_CIDR <- adjustedRandIndex(scHuman@clusters,cols)
ARI_CIDR
## 0.8977449

## This section shows how to alter the parameters of CIDR
scHuman@nPC
## The default number of principal coordinates used in clustering is 4
## Use 6 instead of 4 principal coordinates in clustering
scHuman <- scCluster(scHuman, nPC=6)

plot(scHuman@PC[,c(1,2)],
     col=cols, pch=scHuman@clusters, main="CIDR (nPC=6)",
     xlab="PC1", ylab="PC2")

scHuman@nCluster
## The default number of clusters is 6
## Examine the Calinski-Harabasz Index versus Number of Clusters plot
nCluster(scHuman)
## Try 10 clusters
scHuman <- scCluster(scHuman, nCluster=10)

plot(scHuman@PC[,c(1,2)], col=cols, pch=scHuman@clusters,
     main="CIDR (nPC=6, nCluster=10)", xlab="PC1", ylab="PC2")

ARI_CIDR <- adjustedRandIndex(scHuman@clusters,cols)
ARI_CIDR
## 0.7294235

## Use a different imputation weighting threshold
## Default 6.55
scHuman@wThreshold <- 8
scHuman <- scDissim(scHuman)
scHuman <- scPCA(scHuman)
scHuman <- nPC(scHuman)
scHuman <- scCluster(scHuman)

plot(scHuman@PC[,c(1,2)], col=cols, pch=scHuman@clusters,
     main="CIDR (wThreshold=8)", xlab="PC1", ylab="PC2")

ARI_CIDR <- adjustedRandIndex(scHuman@clusters,cols)
ARI_CIDR
## 0.8101352
