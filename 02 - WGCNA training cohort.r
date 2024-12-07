
load("phenowgcna.rda")

## WGCNA
library(WGCNA)
matrix<-as.matrix(t(data))

# filtration of the matrix
gsg <- goodSamplesGenes(matrix, verbose = 3)
gsg$allOK 
if (!gsg$allOK) {
  matrix <- matrix[gsg$goodSamples, gsg$goodGenes]
}
matrix<-scale(matrix)
# Cluster and detection of outliers
sampleTree <- hclust(dist(matrix), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")



# Choose a soft threshold power- USE A SUPERCOMPUTER IRL ------------------------------------

powers = c(c(1:10), seq(from =4, to=30, by=2)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(matrix, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")


#build a adjacency "correlation" matrix
enableWGCNAThreads()
softPower =16
adjacency = adjacency(matrix, power = softPower, type = "signed") #specify network type
head(adjacency)

# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM



# Generate Modules --------------------------------------------------------
# Generate Modules --------------------------------------------------------

library(flashClust)
# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=4, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(matrix, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")
table(dynamicMods)

TOMplot(dissTOM,geneTree,dynamicColors)


## mds plot
cmd1=cmdscale(as.dist(dissTOM),2)

plot(cmd1,col=dynamicColors,main="MDS plot",xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")


#plots tree showing how the eigengenes cluster together
#INCLUE THE NEXT LINE TO SAVE TO FILE
MEDissThres = 0.0
merge = mergeCloseModules(matrix, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs


#plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="cluster.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

pheno%>%select(!ends_with("_NA"))->pheno
#Define number of genes and samples
nGenes = ncol(matrix)
nSamples = nrow(matrix)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(matrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, pheno, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)





sizeGrWindow(12, 12)
par(mar= c(20, 10, 3, 3))
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(pheno),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.7,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))





# compute   intramodulle connectivity
moduleColors <- dynamicColors
geneModuleMembership <- as.data.frame(cor(matrix, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
write.csv(geneModuleMembership,file="genemodulecor.csv",row.names=T)
write.csv(MMPvalue,file="genemodulepvalue.csv",row.names=T)

# correlation between phenotype and genes
geneTraitSignificance <- as.data.frame(cor(matrix, pheno, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
write.csv(geneTraitSignificance,file="genetraitcor.csv",row.names=T)
write.csv(GSPvalue,file="genetraitpvalue.csv",row.names=T)



# identify hub genes
for (module in unique(moduleColors)) {
	df<-geneModuleMembership
	colnames(df)<-sub("^[^ME]*ME","",colnames(df))
	moduleGenes <- (moduleColors == module)
	cat("\nModule: ", module)
	cat("\nTop hub genes based on module membership and trait significance:\n")
 	hubGenes <- rownames(t(matrix))[moduleGenes]
	df<-df[row.names(df)%in%hubGenes,] 
	print(df%>%arrange(desc(!!sym(module))))
}




## connectivity for grey modulle
library(WGCNA)
library(igraph)


# modulle choice
module <- "grey"
moduleGenes <- (moduleColors == module)

moduleGenes <- rownames(t(matrix))[moduleGenes]
moduleExprs <- matrix[,colnames(matrix)%in%moduleGenes]

# compute (TOM)
TOM <- TOMsimilarityFromExpr(moduleExprs, power = 4)
names<-as.list(moduleGenes)
row.names(TOM) <- names
colnames(TOM)<- names

moduleTOM <- TOM[moduleGenes, moduleGenes]


# convert in igraph object
moduleGraph <- graph_from_adjacency_matrix(moduleTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
moduleGraph <- simplify(moduleGraph, remove.multiple = TRUE, remove.loops = TRUE)

# connectivity of genes
geneConnectivity <- colSums(moduleTOM) - 1
topGenes <- sort(geneConnectivity, decreasing = TRUE)
print(head(topGenes))

# nodes
V(moduleGraph)$name <- colnames(moduleExprs)
V(moduleGraph)$size <- geneConnectivity * 5 / max(geneConnectivity)

# network
plot(moduleGraph, vertex.label = V(moduleGraph)$name, vertex.label.cex = 1,vertex.label.dist=1, 
vertex.size = V(moduleGraph)$size, edge.width = E(moduleGraph)$weight * 2, 
layout = layout_nicely)

## gene names from grey modulle
grey<-colnames(matrix)[moduleColors=="grey"]
grey<-data.frame(colnames(matrix)[moduleColors=="grey"])
write.table(grey,file="module_grey.txt")

## gene names from turquoise modulle
turquoise<-colnames(matrix)[moduleColors=="turquoise"]
turquoise<-data.frame(colnames(matrix)[moduleColors=="turquoise"])
write.table(turquoise,file="module_turquoise.txt")

## gene names from green modulle
green<-colnames(matrix)[moduleColors=="green"]
green<-data.frame(colnames(matrix)[moduleColors=="green"])
write.table(green,file="module_green.txt")

## gene names from red modulle
red<-colnames(matrix)[moduleColors=="red"]
red<-data.frame(colnames(matrix)[moduleColors=="red"])
write.table(red,file="module_red.txt")

## gene names from blue modulle
blue<-colnames(matrix)[moduleColors=="blue"]
blue<-data.frame(colnames(matrix)[moduleColors=="blue"])
write.table(blue,file="module_blue.txt")


## gene names from black modulle
black<-colnames(matrix)[moduleColors=="black"]
black<-data.frame(colnames(matrix)[moduleColors=="black"])
write.table(black,file="module_black.txt")


## gene names from pink modulle
pink<-colnames(matrix)[moduleColors=="pink"]
pink<-data.frame(colnames(matrix)[moduleColors=="pink"])
write.table(pink,file="module_pink.txt")


## gene names from brown modulle
brown<-colnames(matrix)[moduleColors=="brown"]
brown<-data.frame(colnames(matrix)[moduleColors=="brown"])
write.table(brown,file="module_brown.txt")


## gene names from yellow modulle
yellow<-colnames(matrix)[moduleColors=="yellow"]
yellow<-data.frame(colnames(matrix)[moduleColors=="yellow"])
write.table(yellow,file="module_yellow.txt")







library(clusterProfiler)
library(org.Hs.eg.db)  

gene_lists <- list(
    turquoise = turquoise[[1]],
    green = green[[1]],
	red = red[[1]],
	pink = pink[[1]],
	brown = brown[[1]],
	black = black[[1]],
	yellow = yellow[[1]]
	)



go_enrichment <- compareCluster(
    geneCluster = gene_lists, 
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,  
    keyType = "SYMBOL",  
    ont = "BP",            
    pvalueCutoff = 0.05,  
    qvalueCutoff = 0.05    
)


library(ggplot2)
dotplot(go_enrichment, showCategory = 4) +
    ggtitle("GO::BP Enrichment WGCNA ferrptosis modules in Medulloblastoma")

  # connexion between terms
cnetplot(go_enrichment, showCategory = 3)

library(enrichplot)
go_enrichment <- pairwise_termsim(go_enrichment)
emapplot(go_enrichment,showCategory = 8)





