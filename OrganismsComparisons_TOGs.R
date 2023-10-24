#May 24th, 2021
#Comparing how many TOGs (4+) other organisms have on their genome in comparison to Tvag.


# Packages ----------------------------------------------------------------

library(data.table)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(dplyr)
library(tidyverse)
#remotes::install_github("asteves/tayloRswift")
library(tayloRswift)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Homo.sapiens")
#library(Homo.sapiens)
library(Biostrings)
library(ggtext)

##HUMAN:

# Turn the file into a data frame -----------------------------------------





#It was really difficult to try do load the whole of the human genome to R but I found the tRNA database and got data from
#there!! :) http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/
humanFasta <- readDNAStringSet("H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/hg19-tRNAs.fa")
name <- names(humanFasta)
seq <- paste(humanFasta)
df <- data.frame(name, seq)
humangenome <- df


# Split the first column based on commas ----------------------------------




#write.table(df, file="H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/ShortStack2/ShortStack_1600395761/HumanEditedFASTA.txt", row.names=TRUE)


# How many tRNA transcripts are there? From fasta do data frame ------------------------------------


#Making the human tRNAs dataset look nice
dim(humangenome)
humangenome <- humangenome %>%
  separate(name, c("name", "type", "etc"), "\\)")
humangenome$type <- paste0(humangenome$type, ")")
dim(humangenome)
humantRNAtypes <- table(humangenome$type)
humantRNAtypes <- as.data.frame(humantRNAtypes)

#write.table(humangenome, file = "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/humantRNAtypes.txt",
           # append=FALSE, quote=FALSE, sep = "\t")

#Making the human TOGs dataset look nice
humanTOGs <- humangenome[humangenome$seq %like% "^GGG",]

dim(humanTOGs)
humanTOGstypes <- table(humanTOGs$type)
humanTOGstypes <- as.data.frame(humanTOGstypes)

#Augusto thinks I don't need to actually separate the different codons of tRNA so I should manipulate this
#table a bit differently now. Maybe separate the name from the parenthesis and then table it again?

humantRNAtypes <- humantRNAtypes %>%
  separate(Var1, c("type", "codon"), " ")

humantRNAtypes <- as.data.frame(humantRNAtypes %>%
    group_by(codon) %>%
    summarise(Freq = (sum(Freq))) )

#Perfect!! Now, for the TOGs:

humanTOGstypes <- humanTOGstypes %>%
  separate(Var1, c("type", "codon"), " ")

humanTOGstypes <- as.data.frame(humanTOGstypes %>%
                                  group_by(codon) %>%
                                  summarise(Freq = (sum(Freq))) )

#Now, creating a new dataset so I can merge both:
humanTOGstypes2 <- humanTOGstypes
humanTOGstypes2$codon <- paste0(humanTOGstypes2$codon, " (TOGs)")
comparisonhumanTOGs <- rbind(humantRNAtypes, humanTOGstypes2)




# Beautiful ggplots of it :) ----------------------------------------------



ggplot(data = comparisonhumanTOGs %>% gather(codon, Freq), 
       aes(x = codon, y = Freq, fill = Freq)) + 
  theme_minimal() +
  geom_bar(stat = 'identity', position = 'dodge', width=0.7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_taylor(palette = "lover", guide = "none", discrete = FALSE) +
  xlab("") + ylab("Count") + ggtitle("Types of Homo sapiens RNAs") +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.text.x = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.position = "none") 


#Maybe I should select only the RNAs that have TOGs rather than all types of RNAs...
comparisonhumanTOGs2 <- comparisonhumanTOGs %>% filter(
  codon == "Ala" | codon == "Arg" | codon == "Cys" | codon == "Val" |
    codon == "Ala (TOGs)" | codon == "Arg (TOGs)" | codon == "Cys (TOGs)" | codon == "Val (TOGs)"
)

comparisonhumanTOGs2$codon <- as.factor(comparisonhumanTOGs2$codon)

#lover = c("#b8396b", "#ffd1d7", "#fff5cc", "#76bae0", "#b28f81", "#54483e")
#lover

comparisonhumansplot <- ggplot(data = comparisonhumanTOGs2 %>% gather(codon, Freq), aes(x = codon, y = Freq, fill = Freq))+ 
  theme_minimal() +
  geom_bar(stat = 'identity', position = 'dodge', width=0.7, fill=c(rep("#b2d8d8", 4), rep("#76bae0", 4))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("") + ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.text.x = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.position = "none")


#What if I want to add the number of original tRNAs and original TOGs to compare?

codon <- c("tRNA genes", "tRNA genes \n(with TOGs)")
Freq <- c(419, 71)
moreinfohuman <- data.frame(codon, Freq)

#comparisonhumanTOGs3 <- rbind(comparisonhumanTOGs2, moreinfohuman)
#The above group makes the graph too squishy, so I think it's best to separate it!



totalhumansplot <- ggplot(data = moreinfohuman %>% gather(codon, Freq), aes(x = codon, y = Freq, fill = Freq))+ 
  theme_minimal() +
  geom_bar(stat = 'identity', position = 'dodge', width=0.7, fill=c(rep("#b2d8d8", 1), rep("#76bae0", 1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("") + ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.text.x = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.position = "none")


figure <- ggarrange(totalhumansplot, comparisonhumansplot,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)

annotate_figure(figure, 
                top=text_grob("Presence of TOGs in Homo sapiens RNAs",
                              face = "bold", size=16))








##Giardia:
#From https://giardiadb.org/giardiadb/app/downloads/Current_Release/GintestinalisAssemblageAWB/fasta/data/

giardiaFasta <- readDNAStringSet("H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/GiardiaDB-52_GintestinalisAssemblageAWB_AnnotatedTranscripts.fasta")
name <- names(giardiaFasta)
seq <- paste(giardiaFasta)
df <- data.frame(name, seq)
giardiagenome <- df

# How many tRNA and TOGs are there in Giardia? ------------------------------------


dim(giardiagenome)
giardiagenome$name <- gsub('\\|', '-', giardiagenome$name)
giardiagenome <- giardiagenome %>%
  separate(name, c("rna", "gene", "organism", "geneproduct", "transcript", "location", "length", "SO", "SO2", "is"), " - ")


giardiatRNA <- giardiagenome[giardiagenome$geneproduct %like% "gene_product=tRNA",]
dim(giardiatRNA)

giardiaTOGs <- giardiatRNA[giardiatRNA$seq %like% "^GGG",]
dim(giardiaTOGs)
9/73

giardiatRNAtypes <- table(giardiatRNA$geneproduct)
giardiatRNAtypes <- as.data.frame(giardiatRNAtypes)
giardiatRNAtypes$Var1 <- gsub('gene_product=tRNA-', '', giardiatRNAtypes$Var1)

giardiaTOGstypes <- table(giardiaTOGs$geneproduct)
giardiaTOGstypes <- as.data.frame(giardiaTOGstypes)
giardiaTOGstypes$Var1 <- gsub('gene_product=tRNA-', '', giardiaTOGstypes$Var1)


giardiaTOGstypes$Var1 <- paste0(giardiaTOGstypes$Var1, " (TOGs)")
comparisongiardiaTOGs <- rbind(giardiatRNAtypes,giardiaTOGstypes)


comparisongiardiaTOGs2 <- comparisongiardiaTOGs
comparisongiardiaTOGs2 <- comparisongiardiaTOGs %>% filter(
  Var1 == "Ala" | Var1 == "Pro" | Var1 == "Cys" | Var1 == "Trp" |
    Var1 == "Ala (TOGs)" | Var1 == "Pro (TOGs)" | Var1 == "Cys (TOGs)" | Var1 == "Trp (TOGs)"
)

comparisongiardiaplot <- ggplot(data = comparisongiardiaTOGs2 %>% gather(Var1, Freq), aes(x = Var1, y = Freq, fill = Freq))+ 
  theme_minimal() +
  geom_bar(stat = 'identity', position = 'dodge', width=0.7, fill=c(rep("#b2d8d8", 4), rep("#76bae0", 4))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("") + ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.position = "none")


Var1 <- c("tRNA genes", "tRNA genes \n(with TOGs)")
Freq <- c(73, 9)
moreinfogiardia <- data.frame(Var1, Freq)

totalgiardiaplot <- ggplot(data = moreinfogiardia %>% gather(Var1, Freq), aes(x = Var1, y = Freq, fill = Freq))+ 
  theme_minimal() +
  geom_bar(stat = 'identity', position = 'dodge', width=0.7, fill=c(rep("#b2d8d8", 1), rep("#76bae0", 1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("") + ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        axis.text.x = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        legend.position = "none")



figure <- ggarrange(totalgiardiaplot, comparisongiardiaplot,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)

annotate_figure(figure, 
top=text_grob("Presence of TOGs in Giardia intestinalis RNAs",
              face="bold", size=16))




# T.cruzi -----------------------------------------------------------------


##Tcruzi:
#From https://tritrypdb.org/tritrypdb/app/downloads/Current_Release/TcruziBrazilA4/fasta/data/

trypaFasta <- readDNAStringSet("H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/TriTrypDB-52_TcruziBrazilA4_AnnotatedTranscripts.fasta")
name <- names(trypaFasta)
seq <- paste(trypaFasta)
trypagenome <- data.frame(name, seq)

#It does not look like the Tcruzi genome is very well annotated, couldn't find any hit for tRNA-Ala, for example.


# Comparing the TOGs with each other! -------------------------------------


install.packages("devtools")
devtools::install_github("TS404/AlignStat")
library("AlignStat")


# How to create output files of human, giardia and trich TOGs in .fasta -------------------

#26th May, 2021
head(humanFasta)
glimpse(humanFasta)
unlist(humanFasta)
humanFastaTOGs <- humanFasta[humanFasta$seq %like% "^GGG",]
writeXStringSet(humanFastaTOGs,
                "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/humanTOGs.fasta",
                append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
head(humanFastaTOGs)

head(giardiaFasta)
glimpse(giardiaFasta)
giardiaFastaTOGs <- giardiaFasta[giardiaFasta$seq %like% "^GGG",]
writeXStringSet(giardiaFastaTOGs,
                "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/giardiaTOGs.fasta",
                append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
head(giardiaFastaTOGs)

#The above did not work. Maybe just go back to blasting it? 
#Selecting only tRNA-Ala and tRNA-Cys to start. Blasting on TrichDB
humanTOGsAla <- humanTOGs %>% filter(
   type == " Ala (AGC)" | type == " Ala (CGC)" | type == " Ala (TGC)"
)

#write.table(humanTOGs, 
          #file = "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/humanTOGs.txt",
          #append=FALSE,
          #quote=FALSE,
          #sep = "\t")


#write.table(giardiaTOGs, 
           # file = "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/giardiaTOGs.txt",
          #  append=FALSE,
           # quote=FALSE,
          #  sep = "\t")

#write.table(trichTOGs, 
           # file = "H:/Documents/Documents/UniversityofAuckland/Bioinformatics/TVseqNESI/AlignmentOnR/trichTOGs.txt",
          #  append=FALSE,
          #  quote=FALSE,
          #  sep = "\t")
