
#############################################################
####June 2nd, 2021###########################################
####Now it's time to run RUVSeq on my own CLUSTERS data! ####
#############################################################

#Mixing the RUVSeq_TvagDEx_Annotated.R and DEx+EVs.R scripts :)



# Loading packages --------------------------------------------------------

library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
#install.packages("devtools")
#devtools::install_github("G-Thomson/Manu")
library(Manu)
library(RUVSeq)
library(tibble)
library(DESeq2)
library(ggplot2)
library(data.table)

# Reading countdata and sampleinfo in -------------------------------------

results <- read.delim(
  "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\ShortStack2\\ShortStack_1600395761\\Results.txt")
#View(results_EVs)
dim(results)

counts <- read.delim(
  "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\ShortStack2\\ShortStack_1600395761\\Counts.txt")
#View(counts_EVs)
dim(counts)

head(counts)
head(results)

#Selecting only the columns we'll need:
colnames(counts)
countssubsetted <- subset(counts, select=-c(Locus, main))
colnames(countssubsetted)
rownames(countssubsetted)
head(countssubsetted)

#Renaming them
names(countssubsetted) <- c("Clusters", "L1_1", "L1_6", "L1_7", "L1_8", "L2_1", "L2_6", "L2_7", "L2_8")


#Filtering the counts to remove the fragments that haven't been assigned to any cluster
library(tidyr)
countssubsetted$Clusters[countssubsetted$Clusters == ""] <- NA
countssubsetted <- drop_na(countssubsetted, Clusters)
dim (countssubsetted)


#Load the sampleinfo file
sampleinfo <- read.table("H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/featureCounts/WithBoeysHelp/Sampleinfo.txt", 
                         header=TRUE)
sampleinfo$Lane <- as.factor (sampleinfo$Lane)
colnames(sampleinfo)

#storing the cluster ids as rownames
rownames(countssubsetted) <- countssubsetted[,1] 
countssubsetted <- subset(countssubsetted, select= -c(Clusters))
dim(countssubsetted)
rownames(countssubsetted)
head(countssubsetted)

#Making sure that the names of the samples are the same both in counts and sample info data sets
table(colnames(countssubsetted) == sampleinfo$Sample)


# Filtering low counts ----------------------------------------------------
#From RUVSeq_TvagDEx_Annotated.R
#June 2nd, 2021

#Following Robinson 2021 parameters, where they filter samples with more than
#20 counts in at least 2 samples


filter <- apply(countssubsetted, 1, function(x) length(x[x>20])>=2)
filtered <- countssubsetted[filter,]
dim(filtered)
dim(countssubsetted)


#After the filtering, we went from 4640  to 4629 clusters
4629/4640 #only 1% filtered out!

#We store the data in a SeqExpressionSet object from the EDASeq package.

x <- as.factor(rep(c("Cell", "EV", "EV", "EV", "Cell", "EV", "EV", "EV"), each=1))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x,row.names=colnames(filtered)))
set


# Exploratory analysis ---------------------------------------------------
#RLE=log-ratio of read count to median read count across sample

par(mfrow = c(1,1))
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4,4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#The betweenLaneNormalization is a function of EDASeq to normalize the data
#using upper-quartile (UQ) normalization
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#par(mfrow = c(1,1))


# !!! Empirical control genes ---------------------------------------------

#If no genes are known prior to be influenced by the covariates of interest,
#one can obtain a set of "in-silico empirical" negative controls, e.g., least 
#significantly DE genes based on a first-pass DE analysis performed prior to 
#RUVg normalization.

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:100000]))]

#Here, we consider all but the top 500 genes as ranked by edgeR p-values

set2 <- RUVg(set, empirical, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)


# Differential expression analysis with DESeq2 -----------------------------

dds <- DESeqDataSetFromMatrix(countData=counts(set2),
                              colData=pData(set2),
                              design = ~ W_1 + x)
dds <- DESeq(dds)
res <- results(dds)
res


############################################################################
##Using the DEX_loops.R script to go further with my analysis###############
############################################################################
#June 2n, 2021


###########################################################
#Based on the summaries that I got from the analysis above, I got these charts:
names (y)

y$samples$lib.size
barplot(y$samples$lib.size, names=colnames(y), las=2)

title ("Barplot of library sizes for cells and EVs clusters")

logcounts <- cpm(y, log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2)
abline (h=median(logcounts), col="blue")
title("Boxplots of logCPMs for all lanes (normalised)")

#Coloring the MDS plot based on the origin of the RNA
par(mfrow=c(1,2))
sampleinfo$Origin <- as.factor(sampleinfo$Origin)
levels(sampleinfo$Origin)
col.origin <- c("purple", "orange")[sampleinfo$Origin]
data.frame(sampleinfo$Origin, col.origin)

plotMDS(y, col= alpha(col.origin, 0.7), cex=3, pch=17)
legend("topright", fill=c("purple", "orange"), legend=levels(sampleinfo$Origin))
title("Origin")

#Coloring the MDS plot based on the lane of the RNA
levels(sampleinfo$Lane)
col.lane <- c("blue", "pink", "gray", "red")[sampleinfo$Lane]
data.frame(sampleinfo$Lane, col.lane)
plotMDS(y, col= alpha(col.lane, 0.7), cex=3, pch=17)
legend("topright", fill=c("blue", "pink", "gray", "red"), legend=levels(sampleinfo$Lane))
title("Lanes")


#Hierachical clustering with heatmaps
#Estimating the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))#[1:500]
head(select_var)

highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

mypalette <- brewer.pal(8, "RdYlBu")
morecols <- colorRampPalette(mypalette)

col.cell <- c("purple", "orange")[sampleinfo$Origin]
library(gplots)
heatmap.2(highly_variable_lcpm, col=rev(morecols(50)), trace="none", 
          main="Top 50 most variable clusters across samples", 
          colSideColors=col.cell, scale="row",
          cexRow = 1,
          margins = c(4, 8))


#Decided to add a better heatmap, FEb 18th 2022
library(pheatmap)

pheatmap(highly_variable_lcpm[1:50,] , cutree_cols = 2, main = "RUVSeq top 50 clusters")
pheatmap(highly_variable_lcpm[450:500,] , cutree_cols = 2, main = "RUVSeq bottom 50 clusters")


pheatmap(highly_variable_lcpm[450:500,] , cutree_cols = 2, main = "RUVSeq bottom 50 clusters", cellwidth=10)
pheatmap(highly_variable_lcpm[1:50,] , cutree_cols = 2, main = "RUVSeq top 50 clusters", cellwidth = 10)



pheatmap(highly_variable_lcpm[400:500,] , cutree_cols = 2, main = "RUVSeq bottom 100 clusters", cellwidth=10)
pheatmap(highly_variable_lcpm[1:100,] , cutree_cols = 2, main = "RUVSeq top 100 clusters", cellwidth = 10)

#png(file="High_var_genes_heatmap.png")



# SKIP THIS ---------------------------------------------------------------

#DON'T DO THIS - IT'S ALREADY NORMALIZED Apply normalisation to DGEList object
#y <- calcNormFactors(y)
par(mfrow=c(1,2))
plotMD(logcounts,column=1)
abline(h=0,col="blue")
plotMD(logcounts, column=4)
abline(h=0, col="blue")

par(mfrow=c(1,2))
plotMD(y,column=1)
abline(h=0,col="blue")
plotMD(y, column=4)
abline(h=0, col="blue")

group <- paste(sampleinfo$Origin, sampleinfo$Lane, sep=".")
group <- factor(group)
group

design <- model.matrix(~ 0 + group)
design
colnames(design) <- levels (group)
design

par(mfrow=c(1,1))
v <- voom(y, design, plot=TRUE)
v

par (mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2, main="Unnormalised logCPM")
abline(h=median(logcounts), col="red")
boxplot(v$E, xlab="", ylab= "Log2 counts per million", las=2, main="Voom transformed logCPM")
abline(h=median(v$E), col="red")

#Testing for differential expression
fit <- lmFit(v)
names(fit)

#Comparing differential expression among all groups
#Cells vs EVs (Lane6)
cont.matrix6 <- makeContrasts(Cell.1VsEV.6 = Cell.1 - EV.6, levels=design)
cont.matrix6
fit.cont6 <- contrasts.fit(fit, cont.matrix6)
fit.cont6 <-eBayes(fit.cont6)
dim(fit.cont6)
summa.fit6 <- decideTests(fit.cont6)
summary(summa.fit6)
topTable(fit.cont6, coef="Cell.1VsEV.6", sort.by="p")
limma.res6 <- topTable(fit.cont6, coef="Cell.1VsEV.6", sort.by="p", n="Inf")
#write.csv(limma.res6, file=
          #  "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\CellvsLane6_Clusters.csv",
          #row.names=TRUE)
#write.table(limma.res6, file=
           #   "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\CellvsLane6_Clusters.txt", 
          #  row.names=TRUE)
View(limma.res6)


# Making some more analysis with Limma after  DEx -------------------------

par(mfrow=c(1,1))
plotMD(fit.cont6,coef=1,status=summa.fit6[,"Cell.1VsEV.6"])
volcanoplot(fit.cont6,coef=1,highlight=100,names=fit.cont6$genes$SYMBOL)



# Interactive version of the volcano plot using Glimma package: -----------

#install.packages("Glimma")

library(Glimma)
group2 <- group
levels(group2) <- c("Cells", "EV_6", "EV_7", "EV_8")
glXYPlot(x=fit.cont6$coefficients[,1], y=fit.cont6$lods[,1],
         xlab="logFC", ylab="-log of the p value", main="Clusters in all lanes (RUVSeq threshold)", 
         counts=y$counts, groups=group2, status=summa.fit6[,1],
         anno=fit.cont6$genes, side.main="SYMBOL", folder="volcano")


#Testing relative to a threshold(TREAT)
#This allows the cut-off to be done taking into account both the p-value and the fold change (FC),
#so that the results are even more reliable.
#I`m deciding on the FC of 1, which means it`ll only select genes that are 2x or 0.5x expressed.

fit.treat6 <- treat(fit.cont6, lfc=1)
res.treat6 <- decideTests(fit.treat6)
summary(res.treat6)
toptable6 <- topTable(fit.treat6, coef=1, sort.by="p")
#I can now compare to the previous MAplot to see the difference in the number of genes that made it through the cutoff.

par(mfrow=c(1,2))
plotMD(fit.cont6,coef=1,status=summa.fit6[,"Cell.1VsEV.6"])
plotMD(fit.treat6,coef=1,status=res.treat6[,"Cell.1VsEV.6"])

#There`s also an interactive version of it in the Glimma package:
glMDPlot(fit.treat6, coef=1, counts=y$counts, groups=group2, 
         status=res.treat6, side.main="SYMBOL", main="Clusters in all lanes (RUVSeq threshold)",
         folder="md")


# Creating significant topTables -----------------------------------------
#June 1st, 2021


#What I could do now to come up with the list is to go through the limma.res6
#dataset and select only the ones with significant adjusted p-value. Then,
#use the rownames to filter out the counts table so I can have a list of the counts
#for each lane and how significant they could be.

significanttable6 <- filter(limma.res6, "adj.P.Val" <= 0.05)
names <- rownames(significanttable6)
significantcounts <- filtered
significantcounts <- subset(significantcounts, rownames(significantcounts) %in% names)


#And to know what are the sequences of the most abundant clusters, I can filter
#the results dataset:


#storing the cluster ids as rownames

rownames(results) <- results[,2] 
significantresults <- subset(results, rownames(results) %in% names)
significantresults <- subset(significantresults, select= -c(Name))
dim(significantresults)

significantresultscounts <- merge(significantcounts, significantresults,  by="row.names", all=TRUE)  # merge by row names (by=0 or by="row.names")

#write.csv(significantresultscounts, file=
      #      "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\ALLSignificantclusters.csv",
       #   row.names=TRUE)
#write.table(significantresultscounts, file=
 #             "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\ALLSignificantclusters.txt", 
  #          row.names=TRUE)

#But what if I know select only the clusters where the sum of Lanes 6, 7 and 8 are higher than from Lane 1?
significantresultscounts$SumEVs <- apply(significantresultscounts[3:5], 1, sum)
significantresultscounts$SumEVs2 <- apply(significantresultscounts[7:9], 1, sum)
significantresultscounts$EVsTOTAL <- apply(significantresultscounts[30:31], 1, sum)


significantresultscounts$CellsTOTAL <- significantresultscounts$L1_1 + significantresultscounts$L2_1
significantresultscounts <- subset(significantresultscounts, select= -c(SumEVs, SumEVs2))


#AND, FINALLY:

enrichedEVclusters <- filter(significantresultscounts, 
                         EVsTOTAL > CellsTOTAL)


#write.csv(enrichedEVclusters, file=
 #           "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\enrichedEVclustersSUM.csv",
  #        row.names=TRUE)
#write.table(enrichedEVclusters, file=
 #             "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\enrichedEVclustersSUM.txt", 
  #          row.names=TRUE)



##################################################################################
##################################################################################
##############################MY JOB HERE IS DONE!################################
##################################################################################
##################################################################################

#June 2nd, 2021
#Ops. I actually didn't think this through.
#I want the average of columns 3-5, the average of columns 7-9 and then the average
#between the averages. Same with cells. And then I can actually see what's more 
#present in EVs than in cells.

setDT(significantresultscounts)
significantresultscounts$AverageEVs <- significantresultscounts[,.(rMean=rowMeans(.SD,na.rm = T)),.SDcols = c('L1_6','L1_7', 'L1_8', 'L2_6', 'L2_7', 'L2_8')]
significantresultscounts$AverageCells <- significantresultscounts[,.(rMean=rowMeans(.SD,na.rm = T)),.SDcols = c('L1_1','L2_1')]
#Removing the decimals
significantresultscounts$AverageEVs <- floor(significantresultscounts$AverageEVs)
significantresultscounts$AverageCells <- floor(significantresultscounts$AverageCells)
  

#Now I got it:

enrichedEVclusters <- filter(significantresultscounts, 
                             AverageEVs > AverageCells)


#write.csv(enrichedEVclusters, file=
 #           "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\enrichedEVclustersAVERAGE.csv",
  #        row.names=TRUE)
#write.table(enrichedEVclusters, file=
 #             "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\enrichedEVclustersAVERAGE.txt", 
  #          row.names=TRUE)

#Adding a foldchange column

enrichedEVclusters$foldchange <- enrichedEVclusters$AverageEVs / enrichedEVclusters$AverageCells
TOGsenrichedEVclusters <- enrichedEVclusters[enrichedEVclusters$MajorRNA %like% "^GGG",] 
TOGsenrichedEVclusters2 <- enrichedEVclusters[enrichedEVclusters$MajorRNA %like% "^.GGG",] 

samecols <- intersect(colnames(TOGsenrichedEVclusters),colnames(TOGsenrichedEVclusters2))
TOGsandxTOGs <- merge(TOGsenrichedEVclusters, TOGsenrichedEVclusters2,  by=samecols, all=TRUE)

#write.csv(TOGsandxTOGs, file=
 #           "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\TOGsandxTOGs.csv",
  #        row.names=TRUE)
#write.table(TOGsandxTOGs, file=
 #             "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts\\WithBoeysHelp\\Output toptables\\TOGsandxTOGs.txt", 
  #          row.names=TRUE)



##############################################
#TRying to plot a heatmap with the results from RUVSeq!!!!

significantcounts2 <- apply(significantcounts, 1, var)
head(significantcounts2)
select_significant <- names(sort (significantcounts2, decreasing = TRUE))[1:500]
head(select_significant)

highly_variable_selected <- logcounts[select_significant,]
dim(highly_variable_selected)
head(highly_variable_selected)

pheatmap(highly_variable_lcpm[450:500,] , cutree_cols = 2, main = "RUVSeq bottom 50 clusters", cellwidth=10)
pheatmap(highly_variable_lcpm[1:50,] , cutree_cols = 2, main = "RUVSeq top 50 clusters", cellwidth = 10)


#Hierachical clustering with heatmaps
#Estimating the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))#[1:500]
head(select_var)

highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

mypalette <- brewer.pal(8, "RdYlBu")
morecols <- colorRampPalette(mypalette)

col.cell <- c("purple", "orange")[sampleinfo$Origin]
library(gplots)
heatmap.2(highly_variable_lcpm, col=rev(morecols(50)), trace="none", 
          main="Top 50 most variable clusters across samples", 
          colSideColors=col.cell, scale="row",
          cexRow = 1,
          margins = c(4, 8))


# Better heatmap - Feb 21st, 2022 -----------------------------------------

#Decided to add a better heatmap, FEb 18th 2022
library(pheatmap)

pheatmap(highly_variable_lcpm[1:50,] , cutree_cols = 2, main = "RUVSeq top 50 clusters")
pheatmap(highly_variable_lcpm[450:500,] , cutree_cols = 2, main = "RUVSeq bottom 50 clusters")


pheatmap(highly_variable_lcpm[450:500,] , cutree_cols = 2, main = "RUVSeq bottom 50 clusters", cellwidth=10)
pheatmap(highly_variable_lcpm[1:50,] , cutree_cols = 2, main = "RUVSeq top 50 clusters", cellwidth = 10)



pheatmap(highly_variable_lcpm[400:500,] , cutree_cols = 2, main = "RUVSeq bottom 100 clusters", cellwidth=10)
pheatmap(highly_variable_lcpm[1:100,] , cutree_cols = 2, main = "RUVSeq top 100 clusters", cellwidth = 10)



#I still want to get to the toptab;e I generated from RUVSeq but it isn't working,
#so I'll load the output and transform it to make a heatmap from it, rather than
#doing it the other way around. 



RUVSeqOutput <- read.delim(
  "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/featureCounts/WithBoeysHelp/Output toptables/50ComparisonRUVSEqALLAugust24th-Feb21st22.txt",
  header = FALSE)
rownames(RUVSeqOutput) <- RUVSeqOutput$V1
names <- rownames(RUVSeqOutput)


RUVSeqhighly_variable_lcpm <- subset(highly_variable_lcpm, rownames(highly_variable_lcpm) %in% names)

pheatmap(RUVSeqhighly_variable_lcpm, cutree_cols =2, main = "RUVSeq clusters", cellwidth = 10)


#Now, for the ones that are more expressed in the EVs than in the cells (all rRNA)

MoreEVs <- read.delim(
  "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/featureCounts/WithBoeysHelp/Output toptables/MoreEVs.txt",
  header = FALSE)
rownames(MoreEVs) <- MoreEVs$V1
MoreEVsnames <- rownames(MoreEVs)
MoreEVshighly_variable_lcpm <- subset(highly_variable_lcpm, rownames(highly_variable_lcpm) %in% MoreEVsnames)

pheatmap(MoreEVshighly_variable_lcpm, cutree_cols = 2, main = "RUVSeq clusters, more in the EVs", cellwidth = 10)

#Save as 694 x 1149