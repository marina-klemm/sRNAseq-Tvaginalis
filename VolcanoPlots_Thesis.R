##Volcano plots of tRNA in cells and EVs
#For that, I need to mix RUVSeq_TvagDEx_Clusters.r
#but using the files from DEx_tRNA_ALL.R

# Load packages -----------------------------------------------------------

library(magrittr)
library(airway)
library(EnhancedVolcano)
library(DESeq2)
library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
#install.packages("devtools")
#devtools::install_github("G-Thomson/Manu")
library(Manu)
library(RUVSeq)
library(ggplot2)
library(VennDiagram)
library(compare)
library(dplyr)
library(gplots)
library(Glimma)
library(tidyverse)
library(ggrepel)
library(hms)
library(kableExtra)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)

# Importing the sampleinfo and counts files, from AScriptToRuleThemAll.R -------------------------------

filelist <- list.files(path=
                         "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts2",
                       pattern="*.txt", full.names=TRUE)
datalist <- lapply(filelist, function(x)read.table(x, header=T))

datasetnames = c("CDS_all", "CDS_cells", "CDS_EVs",  "rRNA_all", "rRNA_cells", "rRNA_EVs", "tRNA_all", "tRNA_cells", "tRNA_EVs")
names(datalist) <- datasetnames

#Checking if that worked
#View(datalist$CDS_all)
#colnames(datalist$CDS_all)



#Subsetting the count matrices so they only have the essential information for this analysis
subsetting <- function(datalist) {
  subset(datalist, select=-c(Chr, Start, End, Strand, Length))
}
subsetted <- lapply(datalist, subsetting) 

#checking if it worked
#colnames(subsetted$CDS_EVs)
#colnames(subsetted$CDS_all)
#colnames(subsetted$CDS_cells)

#Renaming the columns of the count matrices to make it match the sampleinfo file. The loop won't work on this case because
#each of the datasets has a different size, according to the number of samples on them!
names(subsetted$CDS_all) <- c("Genelist", "L1_1", "L1_6", "L1_7", "L1_8", "L2_1", "L2_6", "L2_7", "L2_8")
names(subsetted$tRNA_all) <- c("Genelist", "L1_1", "L1_6", "L1_7", "L1_8", "L2_1", "L2_6", "L2_7", "L2_8")
names(subsetted$rRNA_all) <- c("Genelist", "L1_1", "L1_6", "L1_7", "L1_8", "L2_1", "L2_6", "L2_7", "L2_8")

names(subsetted$CDS_cells) <- c("Genelist", "L1_1", "L2_1")
names(subsetted$tRNA_cells) <- c("Genelist", "L1_1", "L2_1")
names(subsetted$rRNA_cells) <- c("Genelist", "L1_1", "L2_1")

names(subsetted$CDS_EVs) <- c("Genelist", "L1_6", "L1_7", "L1_8", "L2_6", "L2_7", "L2_8")
names(subsetted$tRNA_EVs) <- c("Genelist", "L1_6", "L1_7", "L1_8", "L2_6", "L2_7", "L2_8")
names(subsetted$rRNA_EVs) <- c("Genelist", "L1_6", "L1_7", "L1_8", "L2_6", "L2_7", "L2_8")

CDS_all <- as.data.frame(subsetted$CDS_all)
CDS_cells <- as.data.frame(subsetted$CDS_cells)
CDS_EVs <- as.data.frame(subsetted$CDS_EVs)

rRNA_all <- as.data.frame(subsetted$rRNA_all)
rRNA_cells <- as.data.frame(subsetted$rRNA_cells)
rRNA_EVs <- as.data.frame(subsetted$rRNA_EVs)

tRNA_all <- as.data.frame(subsetted$tRNA_all)
tRNA_cells <- as.data.frame(subsetted$tRNA_cells)
tRNA_EVs <- as.data.frame(subsetted$tRNA_EVs)

#Importing the sampleinfo and editing it to keep only the EV data
sampleinfo <- read.table(
  "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/featureCounts/WithBoeysHelp/sampleinfo.txt", header=TRUE)
sampleinfo_EVs <- sampleinfo[which(sampleinfo$Origin == "EV"),]
#View(sampleinfo_EVs)
sampleinfo$Lane <- as.factor (sampleinfo$Lane)

rownames(tRNA_all) <- tRNA_all[,1] #storing the gene ids as rownames
tRNA_all <- subset(tRNA_all, select = - c(Genelist)) #removing the first column
table(colnames(tRNA_all) == sampleinfo$Sample) #checking if the names correspond
#View(tRNA_all)



# Calculating the threshold -----------------------------------------------

#I decided to start looking at DEx of tRNA using only the EVs, then comparing with only cells and realized I should be doing all of the tRNAs
#together, otherwise it doesn't make much sense, really. 

dim(tRNA_all)
countdata <- tRNA_all

#Filtering to remove lowly expressed genes
par(mfrow=c(1,1))
myCPM <- cpm(countdata)
head(myCPM)
thresh <- myCPM > 0.8
head(thresh)
table (rowSums(thresh))
keep <- rowSums(thresh) >= 2
counts.keep <- countdata[keep,]
summary (keep)
plot(myCPM[,1],countdata[,1],ylim=c(0,15),xlim=c(0,10), main="Threshold for tRNA genes")
abline (v=0.8)
abline (h=10)


# Making the DGEList ------------------------------------------------------

y <- DGEList(counts.keep)
y

names (y)

y$samples$lib.size
barplot(y$samples$lib.size, names=colnames(y), las=2)
title ("Barplot of library sizes")

logcounts <- cpm(y, log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2)
abline (h=median(logcounts), col="blue")
title("Boxplots of logCPMs (unnormalised)")

#kereru <- get_pal("Kereru")
#hoiho <- get_pal("Hoiho")
#kaka <- get_pal("Kaka")
#takahe <- get_pal("Takahe")


# MDS plot ----------------------------------------------------------------

#Coloring the MDS plot based on the origin of the RNA
par(mfrow=c(1,2), oma=c(0,0,2,0))
sampleinfo$Origin <- as.factor(sampleinfo$Origin)
levels(sampleinfo$Origin)
col.origin <- get_pal("Takahe")[sampleinfo$Origin]
data.frame(sampleinfo$Origin, col.origin)

plotMDS(y, col= alpha(col.origin, 0.7), cex=3, pch=17)
legend("topright", fill= alpha(get_pal("Takahe"), 0.7), legend=levels(sampleinfo$Origin), cex=1.5)
title("Origin", cex.main=1.5)

#Coloring the MDS plot based on the lane of the RNA
levels(sampleinfo$Lane)
col.lane <- get_pal("Kereru")[sampleinfo$Lane]
data.frame(sampleinfo$Lane, col.lane)
plotMDS(y, col= alpha(col.lane, 0.7), cex=3, pch=17)
legend("topright", fill= alpha(get_pal("Kereru"), 0.7), legend=levels(sampleinfo$Lane), cex=1.5)
title("Lane", cex.main=1.5)
mtext("MDS plots for tRNA", outer=TRUE, cex=2)




# Now, adding stuff from the RUVSeq_TvagDEx_Clusters.R --------------------

#Load the sampleinfo file
sampleinfo <- read.table("H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/featureCounts/WithBoeysHelp/Sampleinfo.txt", 
                         header=TRUE)
sampleinfo$Lane <- as.factor (sampleinfo$Lane)
colnames(sampleinfo)




# tRNAs -------------------------------------------------------------------


#Whenever I used countssubsetted in RUVSeq dataset, use tRNA_all instead.

rownames(tRNA_all)
head(tRNA_all)

#Making sure that the names of the samples are the same both in counts and sample info data sets
table(colnames(tRNA_all) == sampleinfo$Sample)


# Filtering low counts ----------------------------------------------------
#From RUVSeq_TvagDEx_Annotated.R
#June 2nd, 2021

#Following Robinson 2021 parameters, where they filter samples with more than
#20 counts in at least 2 samples


filter <- apply(tRNA_all, 1, function(x) length(x[x>20])>=2)
filtered <- tRNA_all[filter,]
dim(filtered)
dim(tRNA_all)


#After the filtering, we went from 4640  to 4629 clusters
251/468 #47% filtered out!

#We store the data in a SeqExpressionSet object from the EDASeq package.

x <- as.factor(rep(c("Cell", "EV", "EV", "EV", "Cell", "EV", "EV", "EV"), each=1))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x,row.names=colnames(filtered)))
set


# Exploratory analysis ---------------------------------------------------
#RLE=log-ratio of read count to median read count across sample

par(mfrow = c(1,2))
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
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:250]))]

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


# Creating a table with the top ones --------------------------------------

# A short function for outputting the tables
#From https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

head(res) %>% 
  knitr_table()

resdf <- as.data.frame(res)
colnames(resdf)

p1 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion  
  geom_point(size =1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1

resdf <- resdf %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & pvalue <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & pvalue <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

head(resdf) %>% 
  knitr_table()


p2 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Expression), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

resdf %>% 
  count(Expression) %>% 
  knitr_table()




resdf <- resdf %>% 
  mutate(
    Significance = case_when(
      abs(log2FoldChange) >= log(2) & pvalue <= 0.05 & pvalue > 0.01 ~ "pvalue 0.05", 
      abs(log2FoldChange) >= log(2) & pvalue <= 0.01 & pvalue > 0.001 ~ "pvalue 0.01",
      abs(log2FoldChange) >= log(2) & pvalue <= 0.001 ~ "pvalue 0.001", 
      TRUE ~ "Unchanged")
  )
head(resdf) %>% 
  knitr_table()

p3 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Significance), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3

resdf %>% 
  count(Expression, Significance) %>% 
  knitr_table()

top <- 10
top_genes <- bind_rows(
  resdf %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>% 
    head(top),
  resdf %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>% 
    head(top)
)

rownames(top_genes) <- gsub("-1", "", rownames(top_genes))

top_genes %>% 
  knitr_table()

#resdf$Genes <- rownames(resdf)
update_geom_defaults("text", list(size = 15))
global_size = 16

p3 <-  p3  + 
  theme(text=element_text(size=global_size),
        plot.title = element_text(hjust=0.5)) +
  labs(title = "tRNAs") + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FoldChange, -log(pvalue,10), label = rownames(top_genes)),
                   size = 4)
p3
  
p3

SAVE
VolcanotRNAsline
909
703

# Save plot to a specified folder
#ggsave(file = "H:/Documents/Documents/UniversityofAuckland/PhDthesis/LaTeXfiles/Scripts/VolcanoPlots_Thesis.R/VolcanotRNA.jpeg", plot = p3, width = 883, height = 552, dpi = 300, units = "px")


top_genes <- as.data.frame(top_genes)
#write.csv(top_genes, file=
 #           "H:/Documents/Documents/UniversityofAuckland/PhDthesis/LaTeXfiles/Scripts/VolcanoPlots_Thesis.R/TopTabletRNA.csv")







# CDS ---------------------------------------------------------------------

#Whenever I used countssubsetted in RUVSeq dataset, use CDS_all instead.

rownames(CDS_all)
head(CDS_all)

#Making sure that the names of the samples are the same both in counts and sample info data sets
rownames(CDS_all) <- CDS_all$Genelist
CDS_all <- CDS_all[,-1]
table(colnames(CDS_all) == sampleinfo$Sample)


# Filtering low counts ----------------------------------------------------
#From RUVSeq_TvagDEx_Annotated.R
#June 2nd, 2021

#Following Robinson 2021 parameters, where they filter samples with more than
#20 counts in at least 2 samples


filter <- apply(CDS_all, 1, function(x) length(x[x>20])>=2)
filtered <- CDS_all[filter,]
dim(filtered)
dim(CDS_all)


#After the filtering, we went from 4640  to 4629 clusters
5264/97530 #95% filtered out!

#We store the data in a SeqExpressionSet object from the EDASeq package.

x <- as.factor(rep(c("Cell", "EV", "EV", "EV", "Cell", "EV", "EV", "EV"), each=1))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x,row.names=colnames(filtered)))
set


# Exploratory analysis ---------------------------------------------------
#RLE=log-ratio of read count to median read count across sample

par(mfrow = c(1,2))
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
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

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


# Creating a table with the top ones --------------------------------------

# A short function for outputting the tables
#From https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

head(res) %>% 
  knitr_table()

resdf <- as.data.frame(res)
colnames(resdf)

p1 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion  
  geom_point(size =1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1

resdf <- resdf %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & pvalue <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & pvalue <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

head(resdf) %>% 
  knitr_table()


p2 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Expression), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

resdf %>% 
  count(Expression) %>% 
  knitr_table()




resdf <- resdf %>% 
  mutate(
    Significance = case_when(
      abs(log2FoldChange) >= log(2) & pvalue <= 0.05 & pvalue > 0.01 ~ "pvalue 0.05", 
      abs(log2FoldChange) >= log(2) & pvalue <= 0.01 & pvalue > 0.001 ~ "pvalue 0.01",
      abs(log2FoldChange) >= log(2) & pvalue <= 0.001 ~ "pvalue 0.001", 
      TRUE ~ "Unchanged")
  )
head(resdf) %>% 
  knitr_table()

p3 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Significance), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3

resdf %>% 
  count(Expression, Significance) %>% 
  knitr_table()

top <- 10
top_genes <- bind_rows(
  resdf %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>% 
    head(top),
  resdf %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>% 
    head(top)
)

rownames(top_genes) <- gsub("-1.*", "", rownames(top_genes))

top_genes %>% 
  knitr_table()

#resdf$Genes <- rownames(resdf)
update_geom_defaults("text", list(size = 15))
global_size = 16

p3 <-  p3  + 
  theme(text=element_text(size=global_size),
        plot.title = element_text(hjust=0.5)) +
  labs(title = "CDS") + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FoldChange, -log(pvalue,10), label = rownames(top_genes)),
                   size = 4)
p3

p3

SAVE
VolcanoCDSline
909
703
p3


top_genes <- as.data.frame(top_genes)
#write.csv(top_genes, file=
         #  "H:/Documents/Documents/UniversityofAuckland/PhDthesis/LaTeXfiles/Scripts/VolcanoPlots_Thesis.R/TopTableCDS.csv")









# rRNA ---------------------------------------------------------------------

#Whenever I used countssubsetted in RUVSeq dataset, use rRNA_all instead.

rownames(rRNA_all)
head(rRNA_all)

#Making sure that the names of the samples are the same both in counts and sample info data sets
rownames(rRNA_all) <- rRNA_all$Genelist
rRNA_all <- rRNA_all[,-1]
table(colnames(rRNA_all) == sampleinfo$Sample)


# Filtering low counts ----------------------------------------------------
#From RUVSeq_TvagDEx_Annotated.R
#June 2nd, 2021

#Following Robinson 2021 parameters, where they filter samples with more than
#20 counts in at least 2 samples


filter <- apply(rRNA_all, 1, function(x) length(x[x>20])>=2)
filtered <- rRNA_all[filter,]
dim(filtered)
dim(rRNA_all)


#After the filtering, we went from 4640  to 4629 clusters
350/668 #48% filtered out!

#We store the data in a SeqExpressionSet object from the EDASeq package.

x <- as.factor(rep(c("Cell", "EV", "EV", "EV", "Cell", "EV", "EV", "EV"), each=1))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x,row.names=colnames(filtered)))
set


# Exploratory analysis ---------------------------------------------------
#RLE=log-ratio of read count to median read count across sample

par(mfrow = c(1,2))
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
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:300]))]

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


# Creating a table with the top ones --------------------------------------

# A short function for outputting the tables
#From https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

head(res) %>% 
  knitr_table()

resdf <- as.data.frame(res)
colnames(resdf)

p1 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion  
  geom_point(size =1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1

resdf <- resdf %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & pvalue <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & pvalue <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

head(resdf) %>% 
  knitr_table()


p2 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Expression), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

resdf %>% 
  count(Expression) %>% 
  knitr_table()




resdf <- resdf %>% 
  mutate(
    Significance = case_when(
      abs(log2FoldChange) >= log(2) & pvalue <= 0.05 & pvalue > 0.01 ~ "pvalue 0.05", 
      abs(log2FoldChange) >= log(2) & pvalue <= 0.01 & pvalue > 0.001 ~ "pvalue 0.01",
      abs(log2FoldChange) >= log(2) & pvalue <= 0.001 ~ "pvalue 0.001", 
      TRUE ~ "Unchanged")
  )
head(resdf) %>% 
  knitr_table()

p3 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Significance), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3

resdf %>% 
  count(Expression, Significance) %>% 
  knitr_table()

top <- 10
top_genes <- bind_rows(
  resdf %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>% 
    head(top),
  resdf %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>% 
    head(top)
)

rownames(top_genes) <- gsub("-1.*", "", rownames(top_genes))

top_genes %>% 
  knitr_table()

#resdf$Genes <- rownames(resdf)
update_geom_defaults("text", list(size = 15))
global_size = 16

p3 <-  p3  + 
  theme(text=element_text(size=global_size),
        plot.title = element_text(hjust=0.5)) +
  labs(title = "rRNAs") + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FoldChange, -log(pvalue,10), label = rownames(top_genes)),
                   size = 4)
p3
p3

SAVE
VolcanorRNA
909
703
p3


top_genes <- as.data.frame(top_genes)
#write.csv(top_genes, file=
 #           "H:/Documents/Documents/UniversityofAuckland/PhDthesis/LaTeXfiles/Scripts/VolcanoPlots_Thesis.R/TopTablerRNA.csv")














# What if I merge all the annotated datasets? -----------------------------

genes_all <- rbind(CDS_all, tRNA_all, rRNA_all)
#rownames(genes_all) <- gsub("-1.*", "", rownames(genes_all))


# genes -------------------------------------------------------------------


#Whenever I used countssubsetted in RUVSeq dataset, use genes_all instead.

rownames(genes_all)
head(genes_all)

#Making sure that the names of the samples are the same both in counts and sample info data sets
table(colnames(genes_all) == sampleinfo$Sample)


# Filtering low counts ----------------------------------------------------
#From RUVSeq_TvagDEx_Annotated.R
#June 2nd, 2021

#Following Robinson 2021 parameters, where they filter samples with more than
#20 counts in at least 2 samples


filter <- apply(genes_all, 1, function(x) length(x[x>20])>=2)
filtered <- genes_all[filter,]
dim(filtered)
dim(genes_all)


#After the filtering, we went from 4640  to 4629 clusters
5865/98666 #95% filtered out!

#We store the data in a SeqExpressionSet object from the EDASeq package.

x <- as.factor(rep(c("Cell", "EV", "EV", "EV", "Cell", "EV", "EV", "EV"), each=1))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x,row.names=colnames(filtered)))
set


# Exploratory analysis ---------------------------------------------------
#RLE=log-ratio of read count to median read count across sample

par(mfrow = c(1,2))
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
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

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


# Creating a table with the top ones --------------------------------------

# A short function for outputting the tables
#From https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

head(res) %>% 
  knitr_table()

resdf <- as.data.frame(res)
colnames(resdf)

p1 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion  
  geom_point(size =1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1

resdf <- resdf %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & pvalue <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & pvalue <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

head(resdf) %>% 
  knitr_table()


p2 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Expression), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

resdf %>% 
  count(Expression) %>% 
  knitr_table()




resdf <- resdf %>% 
  mutate(
    Significance = case_when(
      abs(log2FoldChange) >= log(2) & pvalue <= 0.05 & pvalue > 0.01 ~ "pvalue 0.05", 
      abs(log2FoldChange) >= log(2) & pvalue <= 0.01 & pvalue > 0.001 ~ "pvalue 0.01",
      abs(log2FoldChange) >= log(2) & pvalue <= 0.001 ~ "pvalue 0.001", 
      TRUE ~ "Unchanged")
  )
head(resdf) %>% 
  knitr_table()

p3 <- ggplot(resdf, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Significance), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3

resdf %>% 
  count(Expression, Significance) %>% 
  knitr_table()

top <- 10
top_genes <- bind_rows(
  resdf %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>% 
    head(top),
  resdf %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(pvalue, desc(abs(log2FoldChange))) %>% 
    head(top)
)

rownames(top_genes) <- gsub("-1", "", rownames(top_genes))

top_genes %>% 
  knitr_table()


#resdf$Genes <- rownames(resdf)
update_geom_defaults("text", list(size = 15))
global_size = 16

p3 <-  p3  + 
  theme(text=element_text(size=global_size),
        plot.title = element_text(hjust=0.5)) +
  labs(title = "All RNAs") + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FoldChange, -log(pvalue,10), label = rownames(top_genes)),
                   size = 4)
p3


SAVE
VolcanoAllRNAs
909
703

top_genes <- as.data.frame(top_genes)


#write.csv(top_genes, file=
 #           "H:/Documents/Documents/UniversityofAuckland/PhDthesis/LaTeXfiles/Scripts/VolcanoPlots_Thesis.R/TopTableGenes.csv")



# What are the annotations of the merged top tables? ----------------------


mergedTopTable <- read.csv(
  "H:/Documents/Documents/UniversityofAuckland/PhDthesis/LaTeXfiles/Scripts/VolcanoPlots_Thesis.R/MergedTopTables.csv")

rownames(mergedTopTable) <- mergedTopTable$X
colnames(mergedTopTable)[colnames(mergedTopTable) == "X"] <- "Gene"

head(mergedTopTable)

#Then, I need to load the fasta file, edit it, and add the geneproduct to it.
library(Biostrings)
library(tidyr)

trichFasta <- readDNAStringSet("H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/TrichDB-46_TvaginalisG3_AnnotatedTranscripts.fasta")
name <- names(trichFasta)
seq <- paste(trichFasta)
df <- data.frame(name, seq)
trichgenome <- df

#write.table(df, file="H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/ShortStack2/ShortStack_1600395761/TOGs_FASTA.txt", row.names=TRUE)

trichgenome$name <- gsub('\\|', '-', trichgenome$name)

trichgenome <- trichgenome %>%
  separate(name, c("rna", "gene", "organism", "geneproduct", "transcript", "location", "length", "SO", "SO2", "is"), " - ")

rownames(trichgenome) <- trichgenome$rna

head(trichgenome)

trichgenome$geneproduct <- gsub("gene_product=", "", trichgenome$geneproduct)
rownames(trichgenome) <- gsub("-1", "", rownames(trichgenome))

annotatedTopTable <- merge(mergedTopTable, trichgenome, by="row.names")

# Count number of elements of each type in "geneproduct" column
count_geneproduct <- table(annotatedTopTable$geneproduct)
count_geneproduct <- as.data.frame(count_geneproduct)
#write.csv(count_geneproduct,
 #         "H:/Documents/Documents/UniversityofAuckland/PhDthesis/LaTeXfiles/Scripts/VolcanoPlots_Thesis.R/count_geneproduct.csv")

head(annotatedTopTable)

# Count number of elements of each type in "geneproduct" column based on "Expression"
count_geneproduct <- aggregate(Expression ~ geneproduct, data = annotatedTopTable, function(x) table(x))

head(count_geneproduct)



# Ok, can I see if cells or EVs are more expressed? -----------------------

# Create MA plot
res <- results(dds)
res$log2FoldChange <- ifelse(res$log2FoldChange == -Inf, -10, res$log2FoldChange)
res$log2FoldChange <- ifelse(res$log2FoldChange == Inf, 10, res$log2FoldChange)
plotMA(res, ylim = c(-5, 5))

# Compare test groups
test_group1 <- "Cell"
test_group2 <- "EV"
res_test_group1 <- results(dds, contrast = c("x", test_group1, "EV"))
res_test_group2 <- results(dds, contrast = c("x", test_group2, "Cell"))
plotMA(res, ylim = c(-5, 5))
points(res_test_group1$log2FoldChange, res_test_group1$log10pvalue, col = "blue", cex = 0.5)
points(res_test_group2$log2FoldChange, res_test_group2$log10pvalue, col = "red", cex = 0.5)
legend("topright", legend = c(test_group1, test_group2), col = c("blue", "red"), pch = 1)




# Adding heatmaps ---------------------------------------------------------

library(pheatmap)



