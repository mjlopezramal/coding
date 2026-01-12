In Linux  /home/mjlopez/… creamos 
Create the folder and extend permissions:
```
sudo mkdir -p /home/mjlopez/transcriptome_class_linux
sudo chmod 777 /home/mjlopez/transcriptome_class_linux
```

Run docker in Ubuntu and bid the local “transcriptome_class_linux” to a container folder “/data”
```
docker run -it -v /home/mjlopez/transcriptome_class_linux:/data \
  -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY \
  csicunam/bioinformatics_iamz:latest
cd /data
```
Create 2 folders for the report: report and report_filtered
```
mkdir report report_filtered
```
Create a link between the folders:
```
ln -s /home/vep/test_data/transcripts/brachy brachy
```
Run FastQC in raw reads (in Docker) to créate the report, this document uses as input the 2 files *.gz compressed and we analyzed at the same time. The qualiity reports are saved in the report folder:
```
fastqc brachy/01_RNAseq_raw/*.gz -o report
```
we find te results in htm format in the folder report and zipped. We can open it in Windows. 

```
vep@81b019310d08:/data$ fastqc brachy/01_RNAseq_raw/*.gz -o report
```

##TASK 02: Filter the sequences of ABR3_W_03 sample. Select the most appropiate parameter base don the QC report to improve their quality. We need to do this for the two kind of data, Water content and drought.
For this we can use the software Trimmomatic v.0.39.  For this we créate a new folder.
```
vep@81b019310d08:/data$ mkdir 02_RNAseq_filtered
vep@81b019310d08:/data$ ll
```
Trimming occurs in the order which the steps are specified on the command line.
It is recommended in most cases that adapter clipping, if required, is done as early as possible.

```
vep@81b019310d08:/data$ TrimmomaticSE -threads 4 -phred33 \
brachy/01_RNAseq_raw/ABR3_W_03.fastq.gz \
02_RNAseq_filtered/ABR3_W_03_trimmed.fastq.gz \
ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Con esto ya tenemos un archivo para comprobar la calidad (creo) sin los adaptadores ni las lecutras de mala calidad.
```
vep@81b019310d08:/data$ TrimmomaticSE -threads 4 -phred33 \
brachy/01_RNAseq_raw/ABR3_D_03.fastq.gz \
02_RNAseq_filtered/ABR3_D_03_trimmed.fastq.gz \
ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Run FastQC in filtered/trimmed reads
```
vep@81b019310d08:/data$ fastqc 02_RNAseq_filtered/*.gz -o report_filtered
```
Esto genera unos archivos comprimidos que contienen las secuencias y la calidad. 
Hasta aquí tenemos listo el filtrado, la lectura y podemos hacer el mapeo con kalisto:
Para hacer las cosas ideales, lo suyo es hacer este otro quality check y compararlo con el anterior.

#Task 04: Perform pseudo-alignments of the reads from sample ABR_W_03 against the reference transcriptome (Bdistachyon_314_v3.1.transcript.TagSeq500b.fa.gz9 and indicate the expresión of the Bradi3g51200 gene in est_counts and TPM.
# Index reference the reference transcriptome (kallisto index):
```
vep@81b019310d08:/data$ kallisto index -i Bdistachyon_314_v3.1.transcript.TagSeq500b.idx \
brachy/00_Reference_transcriptome/Bdistachyon_314_v3.1.transcript.TagSeq500b.fa.gz
```
# Quantification for one sample (kallisto quant). Vamos a crear los mapas.
We crearte the folder is wich we are going to sabe the results of Kallisto.
```
vep@81b019310d08:/data$ mkdir 03_quant && mkdir 03_quant/ABR3_W_03
vep@81b019310d08:/data$ kallisto quant -i Bdistachyon_314_v3.1.transcript.TagSeq500b.idx \
-o 03_quant/ABR3_W_03 \
-b 100 --single -l 100 -s 20 <(zcat 02_RNAseq_filtered/ABR3_W_03_trimmed.fastq.gz)
```
Check abundances: more 03_quant/ABR3_W_03/abundance.tsv
Here we can find the last columna, tpm, with the information we need for next steps.

#For Sleuth and WGCNA:
Into /data (Linux WSL) we are going to copy the folders
```
vep@81b019310d08:/data$ cp -r /home/vep/test_data/transcripts/brachy/05_WGCNA/ .
vep@81b019310d08:/data$ cp -r /home/vep/test_data/transcripts/brachy/04_Sleuth/ .
vep@81b019310d08:/data$ chmod -R 777 04_Sleuth 05_WGCNA
vep@81b019310d08:/data$
```
Copy Rstudio scripts:
Copy Rstudio script "Practical_Sleuth.rmd" in \\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\04_Sleuth
Copy Rstudio script "Practical_WGCNA_W_dataset" in \\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\05_WGCNA
Copy Rstudio script "Practical_WGCNA_D_dataset" in \\wsl.localhost\Ubuntu-22.04\home\rsancho\transcriptome_class\05_WGCNA

#########################
######### WGCNA ########
#########################

I've been using the practical file as template for work. The s are answered here but the complete code is bellow.
**QUESTIONS**
### 3.1) How many samples & isoforms are included in TPM_counts_Drought_W_dataset.csv?
The TPM_counts_Drought_W_dataset.csv file includes expression values for 9,940 isoforms across 39 watered (W) RNA-seq samples + target

### 3.2) How many samples are discarded after outlier analysis?
Using:
```
cutreeStatic(sampleTreeW, cutHeight=outliersW, minSize=20)
table(clustW)
```
13 samples were discarded and 26 were kept.

### 3.3) What power value have you set as appropriate for calculating adjacency?

```
sftW$powerEstimate 
[1] 8 
```
the appropiate value is 8

### 3.4) How many co-expression modules are established before and how many after the module merging process?
```
#After you cut the gene dendrogram dynamically (with cutreeDynamic), you get the initial modules: 
dynamicMods <- cutreeDynamic(dendro = geneTree_W, distM = TOM_diss_W,                              deepSplit = 2, pamRespectsDendro = FALSE,                              minClusterSize = 30) table(dynamicMods) 
#The number of unique values in dynamicMods = number of modules before merging: 
length(unique(dynamicMods)) 
#After merging similar modules (using mergeCloseModules): 
mergedModules <- mergeCloseModules(datExprW, dynamicMods, cutHeight = 0.25, verbose = 3) 
mergedColors <- mergedModules$colors 
length(unique(mergedColors))
```
Before merging, 41 modules were detected; after merging, 36 modules remained.

### 3.5) What is the hub isoform (or hub gene) of the cyan module?
```
#identify genes in the cyan module: 
colorOfInterest <- "cyan" 
genesInModule <- names(datExprW)[modulecolors_W == colorOfInterest]

length(genesInModule) #Deberian salirme 106 genes

hubGene <- genesInModule[which.max(kME[genesInModule, paste0("kME", colorOfInterest)])] 
hubGene 
```
The hub isoform of the cyan module is Bradi4g24367.1

### 3.6) According to the module-trait association heat map, which module has the highest positive correlation with the "blwgrd (below ground biomass)" trait?
```
# Load expression data W
load("datExpr_W.RData")

# Load network data W
load("net_W.RData")

# Prepare trait W data set
traitdata <- read.csv("TRAITS_W.csv", header = TRUE, stringsAsFactors = FALSE)
head(traitdata)

# Keep samples and traits
allTraits <- traitdata[, c(4, 6:17) ] # keep samples and traits 
head(allTraits)
# variables/traits numeric
allTraits[, 2:ncol(allTraits)] <- lapply(allTraits[, 2:ncol(allTraits)], 
                                         function(x) as.numeric(as.character(x)))

# check classes
sapply(allTraits, class)
head(allTraits)

# Make sure the rownames of datTraits match the rownames of datExprW
# allTraits must have a unique 'sample' column
allTraits_unique <- allTraits[!duplicated(allTraits$sample), ]  

# Keep only the samples present in datExprW
datTraits <- allTraits_unique[match(rownames(datExprW), allTraits_unique$sample), -1, drop=FALSE]

# Set rownames to sample names
rownames(datTraits) <- rownames(datExprW)

# Check dimensions
dim(datTraits)       # should be nSamples x nTraits
dim(MEs_W)           # should be nSamples x nModules


#Para calcular module-trait correlation: 
moduleTraitCor <- cor(MEs_W, datTraits, use="p")  

# Correlation with 'blwgrd'
corBlwgrd <- moduleTraitCor[, "blwgrd"]

# Module with highest positive correlation
maxModule <- names(corBlwgrd)[which.max(corBlwgrd)]
maxModule
```

TThe hightest positive correlation value is for MEgreenyellow = 0.39212084
The hightest negative correlation value is forMElightcyan = -0.52028579

**COMPLETE TEST**
```
{r setup, include=FALSE}

knitr::opts_chunk$set(echo = TRUE, warning = TRUE, message = TRUE)
```

```{r install_WGCNA_v.1.73}

# https://cran.r-project.org/web/packages/WGCNA/index.html
# https://cran.r-project.org/web/packages/WGCNA/WGCNA.pdf

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GO.db")
BiocManager::install("impute")
BiocManager::install("preprocessCore")

install.packages("WGCNA")

```

```
BiocManager::install("WGCNA", dependencies = TRUE, force = TRUE)
library(WGCNA)
#help(package = 'WGCNA')
```

```
setwd(".")
getwd()

# Water (W)

W_dataset = read.csv("TPM_counts_Drought_W_dataset.csv")
dim(W_dataset) #casi 10000 genes y 40 muestras
names(W_dataset)
#Here we can see that there are 9940 genes + 39 isoforms

# Transpose rows/columns in W dataset

datExprW = as.data.frame(t(W_dataset[, -c(1)]));
names(datExprW) = W_dataset$target_id;
rownames(datExprW) = names(W_dataset)[-c(1)];
```

# Iterative filtering of samples and genes with too many missing entries

```{r goodgenes_D}
# remove genes with too many missing samples in W
#Es una forma de eliminar outlayers

gsgW = goodSamplesGenes(datExprW, verbose = 3);

summary(gsgW)

gsgW$allOK # if TRUE, all genes are ok

if (!gsgW$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsgD$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExprW)[!gsgD$goodGenes], collapse = ", ")));
  if (sum(!gsgW$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExprW)[!gsgW$goodSamples], collapse = ", ")));
  datExprW = datExprW[gsgW$goodSamples, gsgW$goodGenes]
  dim(datExprW)
}
```
Como aparece TRUE quiere decir que no hay outlyers. Todas las muestras están bien. Mirando los genes con 0 expresion.

## Identifying outlier samples

## You can identify outlier samples by using hierarchical clustering.

```{r D_outliers}

# Run all chunk
# cluster and plot the samples (genes that will come later) to see if there are any obvious outliers

#Cosas que me voy encontrando: 
#any(is.na(datExprW))
#any(is.nan(as.matrix(datExprW)))

#zeroVarRows <- which(apply(datExprW, 1, var, na.rm = TRUE) == 0)
#zeroVarCols <- which(apply(datExprW, 2, var, na.rm = TRUE) == 0)

#length(zeroVarRows)
#length(zeroVarCols)

#datExprW <- datExprW[complete.cases(datExprW), ]

#datExprW <- datExprW[apply(datExprW, 1, var, na.rm = TRUE) != 0, ]
#datExprW <- datExprW[, apply(datExprW, 2, var, na.rm = TRUE) != 0]

# W dataset outliers detection
sampleTreeW = hclust(dist(datExprW), method = "average")


png("dendrogram_W.png", width = 1600, height = 1200, units = "px", pointsize = 32, bg = "white")

par(cex = 0.6)
par(mar = c(0,4,2,0))

plot(sampleTreeW,
     main = "W Sample clustering to detect outliers",
     sub = "",
     xlab = "",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2)

outliersW <- 300000
abline(h = outliersW, col = "red")

#dev.off()
```

# Save filtered dataset
Con esto vamos a cargarnos todos los outliers detectados antes
```{r filter_dataset}

## keep only desired samples in W dataset

clustW = cutreeStatic(sampleTreeW, cutHeight=outliersW, minSize = 20)
table(clustW)
keepSamplesW = (clustW==1)
datExprW = datExprW[keepSamplesW, ]
nGenesW = ncol(datExprW)
nSamplesW = nrow(datExprW)

## keep filtering dataset

save(datExprW, file = "datExpr_W.RData")
```
Sale que nos hemos cargado 13 muestras y nos quedan 26 para cosntruir la coexpresion network.

# Soft-thresholding powers selection for network construction
Ahora tenemos que decinir el power, el valor de similaridad beta = Hemos elegido usar unsigned.

```{r D_power}

# Choose a set of soft-thresholding powers

powers = c(c(1:13), seq(from = 14, to=20, by=2))

# Analysis of scale free topology for soft-thresholding

# Unsigned: Both strong positive and strong negative correlations 
#           are considered connections
# Signed: Only positive correlations contribute to network connectivity.
#         Negative correlations reduce similarity toward zero, but correlations 
#         near zero still retain some intermediate similarity.
# Signed hybrid: Positive correlations are treated the same as 
#                in signed networks, but negative correlations are 
#               truncated to zero (completely ignored).

sftW = pickSoftThreshold(datExprW,
                         powerVector = powers,
                         networkType = "unsigned", # default
                         verbose = 5)
```

```{r Plot_D_power}

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power D

#Esto es mio:
# Plot the scale-free topology fit index
plot(sftW$fitIndices[,1], -sign(sftW$fitIndices[,3])*sftW$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     main="Scale independence (W)",
     pch=19,       # plot points
     col="blue")   # color points

# Add the power labels
text(sftW$fitIndices[,1], -sign(sftW$fitIndices[,3])*sftW$fitIndices[,2],
     labels=powers, cex=cex1, col="red")

# Add a horizontal line for the R^2 cutoff
abline(h=0.85, col="red", lty=2)

head(sftW$fitIndices)
```
De acuerdo con R^2 el mejor valor es 5 porque ya tenemos un valor de R2 = 0.8 o más

## Optimal D power = 5
Usamos ese 5 para construir la matriz de adyacencias. \# Unsigned net for D dataset En este paso construimos las redes y definimos los modelos para esas redes (no se si modelos o modulos) teniendo en cuenta el perfil de expresion de las distintas isoformas, si es más parecido o menos.

```{r net_D}

# load dataset
load(file = "datExpr_W.RData");
genes = colnames(datExprW)

# Calling the Adjacency Function
adjacency_W <- adjacency(datExprW, power = 5)

# Topological Overlap Matrix
TOM_W <- TOMsimilarity(adjacency_W)   ## This step takes a long time (5 min) ##
sftW$powerEstimate

# convert this matrix into a dissimilarity matrix 
TOM_diss_W <- 1-TOM_W

# Hierarchical Clustering Analysis

# Creating the dendrogram 
geneTree_W <- hclust(as.dist(TOM_diss_W), method = "average")

# Plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree_W, xlab="", sub="", 
     main = "WATER: Gene clustering on TOM-based dissimilarity", 
labels = FALSE, hang = 0.04)

# Identification of modules D
modules_W <- cutreeDynamic(dendro = geneTree_W, 
                           deepSplit = 2, 
                           pamRespectsDendro = FALSE, 
                           minClusterSize = 30)

table(modules_W) #returns a table of the counts of factor levels in an object. 
            # In this case how many genes are assigned to each created module.

modulecolors_W <- labels2colors(modules_W) # assigns each module number a color

table(modulecolors_W) # returns the counts for each color
                      # (aka the number of genes within each module)


# Plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree_W, modulecolors_W,"Module W",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors for D network")


# Module Eigengene Identification
## A ME (Module Eigengene) is the standardized gene expression profile for 
## a given module.

## An eigengene is the gene whose expression is representative of the the 
## majority of genes expressed within a module.
MElist_W <- moduleEigengenes(datExprW, colors = modulecolors_W) 
MEs_W <- MElist_W$eigengenes 

head(MEs_W)
```
Los colores indican el modelo de coexpresion al que pertenecen. al clicar en cada color podemos ver los diferentes identificadores de los genes coexpresados. Basicamente es como ver un dendograma.

```{r merging}
# Module Merging
ME_diss_W = 1-cor(MElist_W$eigengenes, 
                         use="complete") # Calculate eigengene dissimilarity

METree_W = hclust(as.dist(ME_diss_W), 
                  method = "average") # Clustering eigengenes 
par(mar = c(0,4,2,0)) # setting margin sizes
par(cex = 0.6); # scaling the graphic
plot(METree_W)
abline(h=.25, col = "red") # a height of .25 corresponds to correlation of .75

merge_W <- mergeCloseModules(datExprW, modulecolors_W, cutHeight = .25)


# The merged module colors, assigning one color to each module
mergedcolors_W = merge_W$colors

# Eigengenes of the new merged modules
mergedMEs_W = merge_W$newMEs

# Plot the merged dendrogram

plotDendroAndColors(geneTree_W, cbind(modulecolors_W, mergedcolors_W), 
c("Original Module W", "Merged Module W"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Dendrogram and modules for original and merged modules in W dataset")

# Rename to modulecolors
modulecolors_W = mergedcolors_W

# Construct numerical labels corresponding to the colors
colororder_W = c("grey", standardColors(50));
modulelabels_W = match(modulecolors_W, colororder_W)-1;
MEs_W = mergedMEs_W;
table(modulelabels_W) # Merged modules = 40; original modules = 60
table(modulecolors_W)

# Save module colors and labels for use in subsequent parts
save(MEs_W, modulelabels_W, modulecolors_W, geneTree_W, file = "net_W.RData")

#After you cut the gene dendrogram dynamically (with cutreeDynamic), you get the initial modules:
dynamicMods <- cutreeDynamic(dendro = geneTree_W, distM = TOM_diss_W,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)
table(dynamicMods)
#The number of unique values in dynamicMods = number of modules before merging:
length(unique(dynamicMods))
#After merging similar modules (using mergeCloseModules):
mergedModules <- mergeCloseModules(datExprW, dynamicMods, cutHeight = 0.25, verbose = 3)
mergedColors <- mergedModules$colors
length(unique(mergedColors))
```

# Export genes-modules
Para saber cuales isoformas están en cada modulo usamos esto:
```{r export_modules_genes}

# Export lists of genes within each module of D
load("datExpr_D.RData")
load("net_D.RData")

# Extract gene names
geneNames <- colnames(datExprD)

# Create gene-module table
geneModuleMembership_D <- data.frame(
  Gene = geneNames,
  ModuleColor = modulecolors_D
)
head(geneModuleMembership_D)

# Prepare data frame and export to tsv file
genes_by_module_D <- split(geneModuleMembership_D$Gene, geneModuleMembership_D$ModuleColor)


genes_by_module_df_D <- data.frame(
  ModuleColor = rep(names(genes_by_module_D), times = sapply(genes_by_module_D, length)),
  Gene = unlist(genes_by_module_D, use.names = FALSE)
)

head(genes_by_module_df_D)

# Export to tsv
write.table(genes_by_module_df_D,
            file = "Genes_per_module_D.tsv",
            sep = "\t",        # sep tab
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

```
Que genera un archivo que tiene dos columnas, la primera es el color que tiene el modlu al que pertenece y luego la isoforma.

# Top hub gene per module
Hub genes son los genes más interconectados que normalmente tienen una funcion biologica super importante.

```{r Top_hub_gene}
# The single top hub gene per module:
chooseTopHubInEachModule(
   datExprD, 
   modulecolors_D, 
   omitColors = "grey", 
   power = 8, 
   type = "unsigned")




#identify genes in the cyan module: 
colorOfInterest <- "cyan"  # or whatever exists in your data
genesInModule <- names(datExprW)[modulecolors_W == colorOfInterest]

# Make sure there are genes
length(genesInModule)
#Deberian salirme 106 genes

hubGene <- genesInModule[which.max(kME[genesInModule, paste0("kME", colorOfInterest)])]
hubGene
```

# Module-trait associations
```{r trait_dataset}
# Load expression data W
load("datExpr_W.RData")

# Load network data W
load("net_W.RData")

# Prepare trait W data set
traitdata <- read.csv("TRAITS_W.csv", header = TRUE, stringsAsFactors = FALSE)
head(traitdata)

# Keep samples and traits
allTraits <- traitdata[, c(4, 6:17) ] # keep samples and traits 
head(allTraits)
# variables/traits numeric

allTraits[, 2:ncol(allTraits)] <- lapply(allTraits[, 2:ncol(allTraits)], 
                                         function(x) as.numeric(as.character(x)))

# check classes

sapply(allTraits, class)

head(allTraits)

# Make sure the rownames of datTraits match the rownames of datExprW
# allTraits must have a unique 'sample' column
allTraits_unique <- allTraits[!duplicated(allTraits$sample), ]  

# Keep only the samples present in datExprW
datTraits <- allTraits_unique[match(rownames(datExprW), allTraits_unique$sample), -1, drop=FALSE]

# Set rownames to sample names
rownames(datTraits) <- rownames(datExprW)

# Check dimensions
dim(datTraits)       # should be nSamples x nTraits
dim(MEs_W)           # should be nSamples x nModules

#Para calcular module-trait correlation: 
moduleTraitCor <- cor(MEs_W, datTraits, use="p")  
moduleTraitCor$blwgrd

# Correlation with 'blwgrd'
corBlwgrd <- moduleTraitCor[, "blwgrd"]

# Module with highest positive correlation
maxModule <- names(corBlwgrd)[which.max(corBlwgrd)]
maxModule

corBlwgrd <- moduleTraitCor[, "blwgrd"]  # todas las correlaciones con blwgrd
corBlwgrd
#The highest positive value is greenyellow

# Match the trait data to the expression data by the sample name
Samples <- rownames(datExprD)
traitRows <- match(Samples, allTraits$sample)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]
# Eigengene gene significance: correlation of the trait with the identified 
# module eigengenes

# Define numbers of genes and samples
nGenes = ncol(datExprW)
nSamples = nrow(datExprW)
module.trait.correlation = cor(MEs_W, 
                               datTraits, 
                               use = "p") # Pearson correlation coefficient

# p-value associated with the correlation
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) 

# Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
pdf("heatmap_traits_modules_W.pdf", width = 8, height = 12) # save in pdf
```
