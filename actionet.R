require(ACTIONet)
data <- read.delim("GSE151334_counts.mouse.tsv")
cz <- sapply(data, function(x) length(unique(x)) ==1)

datac <- data[, !cz]
ace = import.ace.from.counts(as.matrix(datac))
ace = normalize.ace(ace)
# Params: dim
ace = reduce.ace(ace)
#Params: depth = 30
#ACTION.out = run_ACTION(S_r, k_max = depth)
ace = run_ACTION(ace)

# Annotate cell-types
data("curatedMarkers_human")
markers = curatedMarkers_human$Blood$PBMC$Ding2019$marker.genes
annot.out = annotate.cells.using.markers(ace, markers)
ace$celltypes = annot.out$Labels

# Visualize output
plot.ACTIONet(ace)


################################################################
hdata <- read.delim("GSE151334_counts.human.tsv")
cz <- sapply(hdata, function(x) length(unique(x)) ==1)

hdatac <- hdata[, !cz]
ace = import.ace.from.counts(as.matrix(hdatac))
ace = normalize.ace(ace)
# Params: dim
ace = reduce.ace(ace)
#Params: depth = 30
#ACTION.out = run_ACTION(S_r, k_max = depth)
ace = run_ACTION(ace)

# Annotate cell-types
data("curatedMarkers_human")
markers = curatedMarkers_human$Blood$PBMC$Ding2019$marker.genes
annot.out = annotate.cells.using.markers(ace, markers)
ace$celltypes = annot.out$Labels

# Visualize output
plot.ACTIONet(ace, ace$louvain, coordinate_attr = "ACTIONet2D")


