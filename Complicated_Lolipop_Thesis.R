#March 8th, 2023

#In this script, I want to make lolipop charts with the trimers that start
#the fragment sequences, and also the percentage of each that is of tRNA,
#rRNA and CDS.
#I will need the dataset from RNAtypesAndSizesPieChart, since that scrip
#is annotated, and also the code to make the loliplot based on the MajorRNA
#from SplittingClustersBasedOnFirst3Nt



# Loading the packages ----------------------------------------------------
library(KernSmooth)
library(tidyverse)
library(Biostrings)
library(GenomicRanges)
library(DESeq2)
library(data.table)
library(dplyr)
library(scales)
library(treemap)
library(ggpubr)

# Importing the dataset ---------------------------------------------------


fastaFile <- readDNAStringSet("H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/TrichDB-46_TvaginalisG3_AnnotatedTranscripts.fasta")
name <- names(fastaFile)
seq <- paste(fastaFile)
df <- data.frame(name, seq) #making it look nicer

##ShortStack Results

filelist <- list.files(path=
                         "H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\featureCounts2",
                       pattern="*.txt", full.names=TRUE)
datalist <- lapply(filelist, function(x)read.table(x, header=T))

datasetnames = c("CDS_all", "CDS_cells", "CDS_EVs",  "rRNA_all", "rRNA_cells", "rRNA_EVs", "tRNA_all", "tRNA_cells", "tRNA_EVs")
names(datalist) <- datasetnames

Reads_all <- read.delim("H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\ShortStack3\\ResultsFiles\\Results_all.txt")
Reads_EVs <- read.delim("H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\ShortStack3\\ResultsFiles\\Results_EVs.txt")
Reads_cells <- read.delim("H:\\Documents\\Documents\\UniversityofAuckland\\Bioinformatics\\TVseqNESI\\ShortStack3\\ResultsFiles\\Results_cells.txt")



#Subsetting the count matrices so they only have the essential information for this analysis
subsetting <- function(datalist) {
  subset(datalist, select=-c(Chr, Start, End, Strand, Length))
}
subsetted <- lapply(datalist, subsetting) 

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


# Editing the df file for all the following uses of it ------------------------------------


df2 <- df
df2$name <- gsub('\\|', '-', df2$name) #changing the | for - to make it easier to split the columns
df2 <- df2 %>%
  separate(name, c("rna", "gene", "organism", "geneproduct", "transcript", "location", "length", "SO", "SO2", "is"), " - ") #splitting the columns
df2types <- table(df2$geneproduct)#counting how many occurences of each type of gene product are there
df2types <- as.data.frame(df2types) #making it into a data frame so I can look at it
df2types$Var1 <- gsub('gene_product=', '', df2types$Var1) #removing the prefix of each of them
#df2types$Var1 <- gsub('.{16}$', '', df2types$Var1) #removing the last 16 letters
#df2types$Var1 <- gsub('^.{6}', '',df2types$Var1) #removing the first 6 letters

df2$rna <- gsub('.{2}$', '', df2$rna)
rownames(df2) <- df2[,1]
colnames(df2)
df2 <- subset(df2, select=-c(rna, gene, organism, transcript, location, SO, SO2, is))
df2$geneproduct <- gsub('^.{13}', '', df2$geneproduct) #removing the first 13 letters
df2$length <- gsub('^.{7}', '', df2$length) #removing the first 7 letters
str(df2)
df2$length <- as.numeric(df2$length)


# Filtering out the lowly expressed counts --------------------------------
rownames(CDS_EVs) <- CDS_EVs$Genelist
rownames(tRNA_EVs) <- tRNA_EVs$Genelist
rownames(rRNA_EVs) <- rRNA_EVs$Genelist

CDS_EVs <- subset(CDS_EVs, select= -c(Genelist))
tRNA_EVs <- subset(tRNA_EVs, select = -c(Genelist))
rRNA_EVs <- subset(rRNA_EVs, select = -c(Genelist))

rownames(CDS_cells) <- CDS_cells$Genelist
rownames(tRNA_cells) <- tRNA_cells$Genelist
rownames(rRNA_cells) <- rRNA_cells$Genelist

CDS_cells <- subset(CDS_cells, select= -c(Genelist))
tRNA_cells <- subset(tRNA_cells, select = -c(Genelist))
rRNA_cells <- subset(rRNA_cells, select = -c(Genelist))


##For CDs
filter <- apply(CDS_EVs, 1, function(x) length(x[x>20])>=6)
filteredCDS_EVs <- CDS_EVs[filter,]
dim(CDS_EVs)
dim(filteredCDS_EVs)

filter <- apply(CDS_cells, 1, function(x) length(x[x>20])>=2)
filteredCDS_cells <- CDS_cells[filter,]
dim(CDS_cells)
dim(filteredCDS_cells)

##For tRNA
filter <- apply(tRNA_EVs, 1, function(x) length(x[x>20])>=6)
filteredtRNA_EVs <- tRNA_EVs[filter,]
dim(tRNA_EVs)
dim(filteredtRNA_EVs)

filter <- apply(tRNA_cells, 1, function(x) length(x[x>20])>=2)
filteredtRNA_cells <- tRNA_cells[filter,]
dim(tRNA_cells)
dim(filteredtRNA_cells)


##For rRNA
filter <- apply(rRNA_EVs, 1, function(x) length(x[x>20])>=6)
filteredrRNA_EVs <- rRNA_EVs[filter,]
dim(rRNA_EVs)
dim(filteredrRNA_EVs)

filter <- apply(rRNA_cells, 1, function(x) length(x[x>20])>=2)
filteredrRNA_cells <- rRNA_cells[filter,]
dim(rRNA_cells)
dim(filteredrRNA_cells)


# Adding a type column to the datasets -----------

filteredCDS_cells$type <- rep("CDs", nrow(filteredCDS_cells))
filteredCDS_EVs$type <- rep("CDs", nrow(filteredCDS_EVs))
filteredrRNA_cells$type <- rep("rRNA", nrow(filteredrRNA_cells))
filteredrRNA_EVs$type <- rep("rRNA", nrow(filteredrRNA_EVs))
filteredtRNA_cells$type <- rep("tRNA", nrow(filteredtRNA_cells))
filteredtRNA_EVs$type <- rep("tRNA", nrow(filteredtRNA_EVs))

# Merging the edited datasets together ------------------------------------

#Merging CELL
mergedfilteredCELL <- rbind(filteredCDS_cells, filteredrRNA_cells)
mergedfilteredCELL <- rbind(mergedfilteredCELL, filteredtRNA_cells)

#Merging EVs
mergedfilteredEVs <- rbind(filteredCDS_EVs, filteredrRNA_EVs)
mergedfilteredEVs <- rbind(mergedfilteredEVs, filteredtRNA_EVs)


# New dataframe just with the amount of each type for the pie chart -------
typesEVs <- data.frame(
  type = c("CDs", "rRNA", "tRNA"),
  value = c(nrow(filteredCDS_EVs), nrow(filteredrRNA_EVs), nrow(filteredtRNA_EVs))
)
head(typesEVs)

typescells <- data.frame(
  type = c("CDs", "rRNA", "tRNA"),
  value = c(nrow(filteredCDS_cells), nrow(filteredrRNA_cells), nrow(filteredtRNA_cells))
)
head(typescells)


# Adding the fragment sequence to the results file ------------------------

####CELLS#####
##CDS
rownames(df2)
rownames(filteredCDS_cells)

rownames(filteredCDS_cells) <- make.unique(gsub("-.*", "", rownames(filteredCDS_cells)))

# Find the matching row indices and add the "seq" column to filteredCDS_cells
match_rows <- match(rownames(filteredCDS_cells), rownames(df2))
filteredCDS_cells$seq <- df2$seq[match_rows]

# Print the updated dataset
filteredCDS_cells

##tRNA
rownames(df2)
rownames(filteredtRNA_cells)

rownames(filteredtRNA_cells) <- make.unique(gsub("-.*", "", rownames(filteredtRNA_cells)))

# Find the matching row indices and add the "seq" column to filteredtRNA_cells
match_rows <- match(rownames(filteredtRNA_cells), rownames(df2))
filteredtRNA_cells$seq <- df2$seq[match_rows]

# Print the updated dataset
filteredtRNA_cells

##rRNA
rownames(df2)
rownames(filteredrRNA_cells)

rownames(filteredrRNA_cells) <- make.unique(gsub("-.*", "", rownames(filteredrRNA_cells)))

# Find the matching row indices and add the "seq" column to filteredrRNA_cells
match_rows <- match(rownames(filteredrRNA_cells), rownames(df2))
filteredrRNA_cells$seq <- df2$seq[match_rows]

# Print the updated dataset
filteredrRNA_cells


####EVs#####
##CDS
rownames(df2)
rownames(filteredCDS_EVs)

rownames(filteredCDS_EVs) <- make.unique(gsub("-.*", "", rownames(filteredCDS_EVs)))

# Find the matching row indices and add the "seq" column to filteredCDS_EVs
match_rows <- match(rownames(filteredCDS_EVs), rownames(df2))
filteredCDS_EVs$seq <- df2$seq[match_rows]

# Print the updated dataset
filteredCDS_EVs

##tRNA
rownames(df2)
rownames(filteredtRNA_EVs)

rownames(filteredtRNA_EVs) <- make.unique(gsub("-.*", "", rownames(filteredtRNA_EVs)))

# Find the matching row indices and add the "seq" column to filteredtRNA_EVs
match_rows <- match(rownames(filteredtRNA_EVs), rownames(df2))
filteredtRNA_EVs$seq <- df2$seq[match_rows]

# Print the updated dataset
filteredtRNA_EVs

##rRNA
rownames(df2)
rownames(filteredrRNA_EVs)

rownames(filteredrRNA_EVs) <- make.unique(gsub("-.*", "", rownames(filteredrRNA_EVs)))

# Find the matching row indices and add the "seq" column to filteredrRNA_EVs
match_rows <- match(rownames(filteredrRNA_EVs), rownames(df2))
filteredrRNA_EVs$seq <- df2$seq[match_rows]

# Print the updated dataset
filteredrRNA_EVs


## Great! Ok, now I will need to plot the lollipop of the annotated fragments 
#and merge it with the percentage of each type of RNA per trimer.



# Then, I need to merge the filtered cells and filtered EVs ---------------

merged_filteredCells <- rbind(filteredCDS_cells, filteredrRNA_cells, filteredtRNA_cells)
merged_filteredEVs <- rbind(filteredCDS_EVs, filteredrRNA_EVs, filteredtRNA_EVs)

# Calculate the column means for numeric rows only
numeric_row_means <- rowMeans(merged_filteredCells[, sapply(merged_filteredCells, is.numeric)], na.rm = TRUE)
merged_filteredCells$average <- numeric_row_means

# Calculate the column means for numeric rows only
numeric_row_means <- rowMeans(merged_filteredEVs[, sapply(merged_filteredEVs, is.numeric)], na.rm = TRUE)
merged_filteredEVs$average <- numeric_row_means



# Looking at the proportion of each possible 3 letters that start  --------

rankedCells <- merged_filteredCells[order(merged_filteredCells$average, decreasing=TRUE),]
rankedCells <- merged_filteredCells[order(merged_filteredCells$average, decreasing=TRUE),]

rankedEVs <- merged_filteredEVs[order(merged_filteredEVs$average, decreasing=TRUE),]
rankedEVs <- merged_filteredEVs[order(merged_filteredEVs$average, decreasing=TRUE),]


rankedCells$three <- substr(rankedCells$seq, 0, 3)
rankedCells$four <- substr(rankedCells$seq, 0, 4)

rankedEVs$three <- substr(rankedEVs$seq, 0, 3)
rankedEVs$four <- substr(rankedEVs$seq, 0, 4)


# Counting the first letters frequency for cells only: --------------------

unique(rankedCells$three)

ids <- unique(rankedCells$three)
ids <- ids[order(ids, decreasing=FALSE)]
dfList <- list()
for (i in ids) {
  dfname <- paste0(i)
  dfList[[dfname]] <- rankedCells[rankedCells$three == i,]
}
dfList



df.lst <- dfList
newdf <- NULL
for (df in df.lst) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  newdf <- rbind(newdf, colSums(df, na.rm=TRUE))  
}
rownames(newdf) <- ids
newdf
newdf <- as.data.frame(newdf)
sum(newdf$average)

newdfCells3 <- newdf[order(newdf$average, decreasing=TRUE),]


unique(rankedCells$four)

ids <- unique(rankedCells$four)
ids <- ids[order(ids, decreasing=FALSE)]
dfList <- list()
for (i in ids) {
  dfname <- paste0(i)
  dfList[[dfname]] <- rankedCells[rankedCells$four == i,]
}
dfList



df.lst <- dfList
newdf <- NULL
for (df in df.lst) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  newdf <- rbind(newdf, colSums(df, na.rm=TRUE))  
}
rownames(newdf) <- ids
newdf
newdf <- as.data.frame(newdf)
sum(newdf$average)

newdfCells4 <- newdf[order(newdf$average, decreasing=TRUE),]


# Counting the first letters frequency for EVs only: --------------------

unique(rankedEVs$three)

ids <- unique(rankedEVs$three)
ids <- ids[order(ids, decreasing=FALSE)]
dfList <- list()
for (i in ids) {
  dfname <- paste0(i)
  dfList[[dfname]] <- rankedEVs[rankedEVs$three == i,]
}
dfList



df.lst <- dfList
newdf <- NULL
for (df in df.lst) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  newdf <- rbind(newdf, colSums(df, na.rm=TRUE))  
}
rownames(newdf) <- ids
newdf
newdf <- as.data.frame(newdf)
sum(newdf$average)

newdfEVs3 <- newdf[order(newdf$average, decreasing=TRUE),]


unique(rankedEVs$four)

ids <- unique(rankedEVs$four)
ids <- ids[order(ids, decreasing=FALSE)]
dfList <- list()
for (i in ids) {
  dfname <- paste0(i)
  dfList[[dfname]] <- rankedEVs[rankedEVs$four == i,]
}
dfList



df.lst <- dfList
newdf <- NULL
for (df in df.lst) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  newdf <- rbind(newdf, colSums(df, na.rm=TRUE))  
}
rownames(newdf) <- ids
newdf
newdf <- as.data.frame(newdf)
sum(newdf$average)

newdfEVs4 <- newdf[order(newdf$average, decreasing=TRUE),]


# Now, plotting it --------------------------------------------------------

newdfEVs3
newdfCells3
newdfEVs4
newdfCells4

dim(newdfCells3)
dim(newdfEVs3)

#I need to change the datasets slightly to be able to plot them:

newdfCells3$type <- rownames(newdfCells3)

newdfCells4$type <- rownames(newdfCells4)

newdfEVs3$type <- rownames(newdfEVs3)

newdfEVs4$type <- rownames(newdfEVs4)



##Adding percentage
percentagecells3 <- newdfCells3 %>% mutate (Percentage = round(average/sum(average)*100,2))
percentagecells4 <- newdfCells4 %>% mutate (Percentage = round(average/sum(average)*100,2))

percentageEVs3 <- newdfEVs3 %>% mutate (Percentage = round(average/sum(average)*100,2))
percentageEVs4 <- newdfEVs4 %>% mutate (Percentage = round(average/sum(average)*100,2))

newdfCells3$percentage <- percentagecells3$Percentage
newdfCells4$percentage <- percentagecells4$Percentage

newdfEVs3$percentage <- percentageEVs3$Percentage
newdfEVs4$percentage <- percentageEVs4$Percentage


ggplot( data = newdfCells3[1:5,], aes ( x = type, y = average)) +
  geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e")

ggplot( data = newdfCells4[1:5,], aes ( x = type, y = average)) +
  geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e")

ggplot( data = newdfEVs3[1:5,], aes ( x = type, y = average)) +
  geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e")

ggplot( data = newdfEVs4[1:5,], aes ( x = type, y = average)) +
  geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e")


#Maybe I can do a lollipop plot!

Cells3 <- ggplot( data = newdfCells3[1:5,], aes ( x = type, y = average)) +
  geom_segment( aes(x=type ,xend=type, y=0, yend=average), color="grey",
                size = 3) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_point(size=6, color="aquamarine3") +
  coord_flip() +
  theme_minimal(base_size=16) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Cells") +
  geom_text(label = with(newdfCells3[1:5,], paste(paste0(percentage, '%'))), vjust=2) 

#Save as: 
795
449


EVs3 <- ggplot( data = newdfEVs3[1:5,], aes ( x = type, y = average)) +
  geom_segment( aes(x=type ,xend=type, y=0, yend=average), color="grey",
                size = 3) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_point(size=6, color="darkorange") +
  coord_flip() +
  theme_minimal(base_size=16) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("EVs") +
  geom_text(label = with(newdfEVs3[1:5,], paste(paste0(percentage, '%'))), vjust=2) 


ggarrange(Cells3, EVs3,
          ncol = 1, nrow = 2)





#Save as: 
795
449


#Maybe I can do a lollipop plot!

Cells4 <- ggplot( data = newdfCells4[1:5,], aes ( x = type, y = average)) +
  geom_segment( aes(x=type ,xend=type, y=0, yend=average), color="grey",
                size = 3) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_point(size=6, color="aquamarine3") +
  coord_flip() +
  theme_minimal(base_size=16) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Cells") +
  geom_text(label = with(newdfCells4[1:5,], paste(paste0(percentage, '%'))), vjust=2) 

#Save as: 
795
449


EVs4 <- ggplot( data = newdfEVs4[1:5,], aes ( x = type, y = average)) +
  geom_segment( aes(x=type ,xend=type, y=0, yend=average), color="grey",
                size = 3) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_point(size=6, color="darkorange") +
  coord_flip() +
  theme_minimal(base_size=16) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("EVs") +
  geom_text(label = with(newdfEVs4[1:5,], paste(paste0(percentage, '%'))), vjust=2) 


ggarrange(Cells4, EVs4,
          ncol = 1, nrow = 2)





# Reordering the plot: ----------------------------------------------------
Cells3 <- newdfCells3[1:5,] %>%
  mutate(type = fct_reorder(type, average)) %>%
  ggplot(aes ( x = type, y = average)) +
  geom_segment( aes(x=type ,xend=type, y=0, yend=average), color="grey",
                size = 3) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_point(size=6, color="aquamarine3") +
  coord_flip() +
  theme_minimal(base_size=16) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Cells") +
  geom_text(label = with(newdfCells3[1:5,], paste(paste0(percentage, '%'))), vjust=2) 



EVs3 <- newdfEVs3[1:5,] %>%
  mutate(type = fct_reorder(type, average)) %>%
  ggplot( aes ( x = type, y = average)) +
  geom_segment( aes(x=type ,xend=type, y=0, yend=average), color="grey",
                size = 3) +
  scale_y_continuous(#limits = c(0,13000000),
    labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_point(size=6, color="darkorange") +
  coord_flip() +
  theme_minimal(base_size=16) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("EVs") +
  geom_text(label = with(newdfEVs3[1:5,], paste(paste0(percentage, '%'))), vjust=2) 


ggarrange(Cells3, EVs3,
          ncol = 1, nrow = 2)

#Save as: 
882
639


Cells4 <- newdfCells4[1:5,] %>%
  mutate(type = fct_reorder(type, average)) %>%
  ggplot(aes ( x = type, y = average)) +
  geom_segment( aes(x=type ,xend=type, y=0, yend=average), color="grey",
                size = 3) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_point(size=6, color="aquamarine3") +
  coord_flip() +
  theme_minimal(base_size=16) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Cells") +
  geom_text(label = with(newdfCells4[1:5,], paste(paste0(percentage, '%'))), vjust=2) 



EVs4 <- newdfEVs4[1:5,] %>%
  mutate(type = fct_reorder(type, average)) %>%
  ggplot( aes ( x = type, y = average)) +
  geom_segment( aes(x=type ,xend=type, y=0, yend=average), color="grey",
                size = 3) +
  scale_y_continuous(#limits = c(0,13000000),
    labels = unit_format(unit = "M", scale = 1e-6)) +
  geom_point(size=6, color="darkorange") +
  coord_flip() +
  theme_minimal(base_size=16) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("EVs") +
  geom_text(label = with(newdfEVs4[1:5,], paste(paste0(percentage, '%'))), vjust=2) 


ggarrange(Cells4, EVs4,
          ncol = 1, nrow = 2)



# Ok, not quite what I wanted because this is with ------------------------
#the annotated genes and not with the fragments.
#Regardless, how many from each trimer are of each annotation?
#The datasets I want are:

head(rankedCells)
rankedCellsU <- rankedCells
rankedCellsU$seq <- gsub("T", "U", rankedCells$seq)
rankedCellsU$three <- gsub("T", "U", rankedCells$three)
rankedCellsU$four <- gsub("T", "U", rankedCells$four)


rankedEVsU <- rankedEVs
rankedEVsU$seq <- gsub("T", "U", rankedEVs$seq)
rankedEVsU$three <- gsub("T", "U", rankedEVs$three)
rankedEVsU$four <- gsub("T", "U", rankedEVs$four)

Cells3Type <- as.data.frame(table(rankedCellsU$three, rankedCellsU$type))
Cells4Type <- as.data.frame(table(rankedCellsU$four, rankedCellsU$type))

EVs3Type <- as.data.frame(table(rankedEVsU$three, rankedEVsU$type))
EVs4Type <- as.data.frame(table(rankedEVsU$four, rankedEVsU$type))


# Selecting only the top fragments and calculating % ----------------------
#The data from the top fragments are in "H:\Documents\Documents\UniversityofAuckland\PhDthesis\LaTeXfiles\Scripts\SplittingClustersBasedOnFirst3Nt.R\"

subCells3 <- subset(Cells3Type, Var1 %in% c("UCC", "GCA", "GGG", "GCC", "CGA"))
subCells4 <- subset(Cells4Type, Var1 %in% c("UCCG", "GCAC", "GGGA", "GCCU", "CGAU"))


subEVs3 <- subset(EVs3Type, Var1 %in% c("UCC", "GCA", "GGG", "GCC", "AGG"))
subEVs4 <- subset(EVs4Type, Var1 %in% c("GCAC", "UCCG", "GGGA", "GCCU", "AGGG"))




Cells3Type <- as.data.frame(table(rankedCells$three, rankedCells$type))
Cells4Type <- as.data.frame(table(rankedCells$four, rankedCells$type))

EVs3Type <- as.data.frame(table(rankedEVs$three, rankedEVs$type))
EVs4Type <- as.data.frame(table(rankedEVs$four, rankedEVs$type))


# Selecting only the top fragments and calculating % ----------------------
#The data from the top fragments are in "H:\Documents\Documents\UniversityofAuckland\PhDthesis\LaTeXfiles\Scripts\SplittingClustersBasedOnFirst3Nt.R\"

subCells3 <- subset(Cells3Type, Var1 %in% c("TCC", "GCA", "GGG", "GCC", "ATG"))
subCells4 <- subset(Cells4Type, Var1 %in% c("TCCG", "GCAC", "GGGA", "GCCT", "GGAC"))


subEVs3 <- subset(EVs3Type, Var1 %in% c("GCA", "TCC", "ATG", "GGG", "GCC"))
subEVs4 <- subset(EVs4Type, Var1 %in% c("GCAC", "TCCG", "GGGA", "GCCT", "ATGA"))


# subEVs3 ----------------------------------------------------------------



# calculate the percentage of each Var2 for each unique Var1
subEVs3_perc <- subEVs3 %>%
  group_by(Var1) %>%
  mutate(perc = Freq/sum(Freq) * 100)



# create separate bar plots for each unique Var1
ggplot(subEVs3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))


subEVs3_perc$Var1 <- gsub("ATG", "AGG", subEVs3_perc$Var1)
subEVs3_perc$Var1 <- gsub("TCC", "UCC", subEVs3_perc$Var1)

subEVs3_perc$Var1 <- factor(subEVs3_perc$Var1, levels=c("GCA", "UCC", "GGG", "GCC", "AGG"))


ggplot(subEVs3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + theme(legend.position="none")

ggplot(subEVs3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + theme(legend.position="none") +
  ggtitle("EVs")



ggplot(subEVs3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_y", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  #theme(axis.text.x = element_text(vjust = 0.5)) +
  # coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + 
  theme(legend.position="none",
        text= element_text(size=16)) +
  ggtitle("EVs") 

# subCells3 ----------------------------------------------------------------



# calculate the percentage of each Var2 for each unique Var1
subCells3_perc <- subCells3 %>%
  group_by(Var1) %>%
  mutate(perc = Freq/sum(Freq) * 100)



# create separate bar plots for each unique Var1
ggplot(subCells3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))


subCells3_perc$Var1 <- gsub("ATG", "CGA", subCells3_perc$Var1)
subCells3_perc$Var1 <- gsub("TCC", "UCC", subCells3_perc$Var1)

subCells3_perc$Var1 <- factor(subCells3_perc$Var1, levels=c("UCC", "GCA", "GGG", "GCC", "CGA"))


ggplot(subCells3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + theme(legend.position="none")

ggplot(subCells3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + theme(legend.position="none") +
  ggtitle("Cells")



# For saving --------------------------------------------------------------


PlotSubCells3 <- ggplot(subCells3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  #theme(axis.text.x = element_text(vjust = 0.5)) +
 # coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + 
  theme(legend.position="none",
        text= element_text(size=20),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), 
        axis.title.y = element_blank()) 

PlotSubEVs3 <- ggplot(subEVs3_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  #theme(axis.text.x = element_text(vjust = 0.5)) +
  # coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + 
  theme(legend.position="none",
        text= element_text(size=20),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), 
        axis.title.y = element_blank()) 


# This is done after running SplittingClustersBased until 3 ----------


ggarrange(Cells, PlotSubCells3,
          ncol = 2, nrow = 1)

ggarrange(EVs, PlotSubEVs3,
          ncol=2, nrow = 1)

SAVE
1185
705


# subCells4 ----------------------------------------------------------------



# calculate the percentage of each Var2 for each unique Var1
subCells4_perc <- subCells4 %>%
  group_by(Var1) %>%
  mutate(perc = Freq/sum(Freq) * 100)



# create separate bar plots for each unique Var1
ggplot(subCells4_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))


subCells4_perc$Var1 <- gsub("TCCG", "UCCG", subCells4_perc$Var1)
subCells4_perc$Var1 <- gsub("GCCT", "GCCU", subCells4_perc$Var1)
subCells4_perc$Var1 <- gsub("GGAC", "CGAU", subCells4_perc$Var1)

subCells4_perc$Var1 <- factor(subCells4_perc$Var1, levels=c("UCCG", "GCAC", "GGGA", "GCCU", "CGAU"))


ggplot(subCells4_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + theme(legend.position="none")

ggplot(subCells4_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + theme(legend.position="none") +
  ggtitle("Cells")



PlotSubCells4<- ggplot(subCells4_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_y", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  #theme(axis.text.x = element_text(vjust = 0.5)) +
  # coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + 
  theme(legend.position="none",
        text= element_text(size=20),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), 
        axis.title.y = element_blank()) 


# This is done after running SplittingClustersBased until 4 ----------
Cells
PlotSubCells4

ggarrange(Cells, PlotSubCells4,
          ncol = 2, nrow = 1)

SAVE
1185
705



# subEVs4 ----------------------------------------------------------------

# calculate the percentage of each Var2 for each unique Var1
subEVs4_perc <- subEVs4 %>%
  group_by(Var1) %>%
  mutate(perc = Freq/sum(Freq) * 100)



# create separate bar plots for each unique Var1
ggplot(subEVs4_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))



subEVs4_perc$Var1 <- gsub("ATGA", "UCCG", subEVs4_perc$Var1)
subEVs4_perc$Var1 <- gsub("GCCT", "GCCU", subEVs4_perc$Var1)
subEVs4_perc$Var1 <- gsub("GGGA", "AGGG", subEVs4_perc$Var1)
subEVs4_perc$Var1 <- gsub("TCCG", "GGGA", subEVs4_perc$Var1)


subEVs4_perc$Var1 <- factor(subEVs4_perc$Var1, levels=c("GCAC", "UCCG", "GGGA", "GCCU", "AGGG"))


ggplot(subEVs4_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + theme(legend.position="none")

ggplot(subEVs4_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_x", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + theme(legend.position="none") +
  ggtitle("EVs")



PlotSubEVs4<- ggplot(subEVs4_perc, aes(x = Var2, y = perc, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_grid(Var1 ~ ., scales = "free_y", 
             switch = "y") +
  labs(title = "", 
       x = "", y = "Percentage") +
  #theme(axis.text.x = element_text(vjust = 0.5)) +
  # coord_flip() +
  geom_text(aes(label = paste0(round(perc, 1), "%")), 
            position = position_stack(vjust = 0.5))  + 
  theme(legend.position="none",
        text= element_text(size=20),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), 
        axis.title.y = element_blank()) 


# This is done after running SplittingClustersBased until 4 ----------
EVs
PlotSubEVs4

ggarrange(EVs, PlotSubEVs4,
          ncol = 2, nrow = 1)



SAVE
1185
705