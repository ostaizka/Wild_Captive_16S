#Load the data
#Genus level

library(reshape2)
library(gplots)
library(purrr)
library(hilldiv)
library(ape)
library("phyloseq")
library("ggplot2")
library(vegan)
library(DESeq2)
library(microbiome)
library("RColorBrewer")
library(tidyverse)
library(mctoolsr)
library(metacoder)
library(metagenomeSeq)
library(metagMisc)
library(phytools)
library(dplyr)
library("ggpubr")
library(dmetar)
library(meta)

#Metadata from airtable
metadata <- read.table("/Users/ostaizka/reanalysis/Datafiles-Grid.csv", header=T, row.names=1,check.names=F, sep="\t")
names(metadata)[names(metadata) == "Species (from Individuals)"] <- "Species"
names(metadata)[names(metadata) == "Origin (from Individual)"] <- "Origin"
names(metadata)[names(metadata) == "Study (from Individuals)"] <- "Study"
#colnames(metadata) <- c("Datafiles", "Selected", "Origin", "Species", "Order", "Class", "Latitude", "Longitude", "Study", "16S_region", "F_Primer", 
#                        "R_Primer", "Read_length", "Archived_F", "Processing_status", "Status")
colnames(metadata.filtered) <- gsub(" ","_",colnames(metadata.filtered))
metadata$Origin <- gsub("Zoo","Captivity",metadata$Origin)
metadata$Origin <- gsub("Laboratory","Captivity",metadata$Origin)
metadata$Origin <- gsub("Domestic","Captivity",metadata$Origin)
hierarchy_all <- tibble::rownames_to_column(metadata, "Samples")

#code names
setwd("~/reanalysis/Species")
list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}
code.list <- list.dirs()

####GET TABLES####
setwd("~/reanalysis/Species")
for (code in code.list){
  print(code)
#Load files at genus level
files.g = list.files(path=code, pattern="*.tax_genus.txt")
names.g = gsub(".tax_genus.txt","",files.g)
filelist.g = lapply(paste(code,files.g,sep="/"), function(x)read.table(x, header=T,sep="\t"))
#Note some erroneous files can prevent all files from loading (remove erroneous ones)
#Transform list to matrix
filelist.g %>% purrr::reduce(full_join, by = "Taxa") -> count.table.g
rownames(count.table.g) <- count.table.g[,1]
count.table.g <- count.table.g[,-1]
#Name columns
colnames(count.table.g) <- names.g
#Replace NAs by 0s
count.table.g[is.na(count.table.g)] <- 0
write.table(count.table.g, paste("/Users/ostaizka/reanalysis/table/count_genus_",code,".tsv",sep=""))
}
for (code in code.list){
  print(code)
#Load files at family level
files.f = list.files(path=code, pattern="*.tax_family.txt")
names.f = gsub(".tax_family.txt","",files.f)
filelist.f = lapply(paste(code,files.f,sep="/"), function(x)read.table(x, header=T,sep="\t"))
#Note some erroneous files can prevent all files from loading (remove erroneous ones)
#Transform list to matrix
filelist.f %>% purrr::reduce(full_join, by = "Taxa") -> count.table.f
rownames(count.table.f) <- count.table.f[,1]
count.table.f <- count.table.f[,-1]
#Name columns
colnames(count.table.f) <- names.f
#Replace NAs by 0s
count.table.f[is.na(count.table.f)] <- 0
write.table(count.table.f, paste("/Users/ostaizka/reanalysis/table/count_family_",code,".tsv",sep=""))
}
for (code in code.list){
  print(code)
#Load files at phylum
files.p = list.files(path=code, pattern="*.tax_phylum.txt")
names.p = gsub(".tax_phylum.txt","",files.p)
filelist.p = lapply(paste(code,files.p,sep="/"), function(x)read.table(x, header=T,sep="\t"))
#Note some erroneous files can prevent all files from loading (remove erroneous ones)
#Transform list to matrix
filelist.p %>% purrr::reduce(full_join, by = "Taxa") -> count.table.p
rownames(count.table.p) <- count.table.p[,1]
count.table.p <- count.table.p[,-1]
#Name columns
colnames(count.table.p) <- names.p
#Replace NAs by 0s
count.table.p[is.na(count.table.p)] <- 0
write.table(count.table.p, paste("/Users/ostaizka/reanalysis/table/count_phylum_",code,".tsv",sep=""))
}

####Filtering tables####
for (code in code.list){
  print(code)
  table.g <- read.table(paste("/Users/ostaizka/reanalysis/table/count_genus_",code,".tsv",sep=""))  
  assign(paste("table.g.",code,sep=""),get("table.g"))
}
for (code in code.list){
  print(code)
  final.table.p <- read.table(paste("/Users/ostaizka/reanalysis/table/count_phylum_",code,".tsv",sep=""))  
  assign(paste("final.table.p.",code,sep=""),get("final.table.p"))
}
for (code in code.list){
  print(code)
  final.table.f <- read.table(paste("/Users/ostaizka/reanalysis/table/count_family_",code,".tsv",sep=""))  
  assign(paste("final.table.f.",code,sep=""),get("final.table.f"))
}

for (code in code.list){
counts.filtdepth <- depth_filt(get(paste("table.g.",code,sep="")),1000)
samples.kept <- colnames(counts.filtdepth)
filtered.g <- counts.filtdepth[rowSums(counts.filtdepth)>0, ]
assign(paste("filtered.g.",code,sep=""),get("filtered.g"))
write.table(get(paste("filtered.g.",code,sep="")),paste("/Users/ostaizka/reanalysis/table/countfiltdepth_Genus_",code,".tsv",sep=""))
final.table.p <- get(paste("final.table.p.",code,sep=""))
filtered.p <- final.table.p[,colnames(final.table.p) %in% samples.kept]
assign(paste("filtered.p.",code,sep=""),get("filtered.p"))
write.table(get(paste("filtered.p.",code,sep="")),paste("/Users/ostaizka/reanalysis/table/countfiltered_Phylum_",code,".tsv",sep=""))
final.table.f <- get(paste("final.table.f.",code,sep=""))
filtered.f <- final.table.f[,colnames(final.table.f) %in% samples.kept]
assign(paste("filtered.f.",code,sep=""),get("filtered.f"))
write.table(get(paste("filtered.f.",code,sep="")),paste("/Users/ostaizka/reanalysis/table/countfiltered_Family_",code,".tsv",sep=""))
}

## Remove NA taxa
taxonlist_NA <- read.table("NA_taxa.txt")
taxonlist_NA <-  as.character(taxonlist_NA[,1])

for (code in code.list){
  countfiltdepth <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltdepth_Genus_",code,".tsv",sep=""))  
  countfiltered_Genus <- countfiltdepth[! rownames(countfiltdepth) %in% taxonlist_NA,]
  assign(paste("countfiltered_Genus_",code,sep=""),get("countfiltered_Genus"))
  write.table(get(paste("countfiltered_Genus_",code,sep="")),paste("/Users/ostaizka/reanalysis/table/countfiltered_Genus_",code,".tsv",sep=""))
}

####Diversity analysis based on abundace####
taxa = "Genus"
summary_q0 <- c()
for (code in code.list){
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divq0 <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  divq0_table <- divq0$data
  divq0_wild <- divq0_table[divq0_table$Group == "Wild",]
  n_wild <- nrow(divq0_wild)
  mean_wild <- mean(divq0_wild[,c("Value")])
  sd_wild <- sd(divq0_wild[,c("Value")])
  divq0_captive <- divq0_table[divq0_table$Group == "Captivity",]
  n_captive <- nrow(divq0_captive)
  mean_captive <- mean(divq0_captive[,c("Value")])
  sd_captive <- sd(divq0_captive[,c("Value")])
  summary_q0 <- rbind(summary_q0,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
  colnames(summary_q0) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
  rownames(summary_q0) <- code.list
  write.table(summary_q0, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq0_",taxa,".tsv",sep=""))
  
  #q=1
  summary_q1 <- c()
  for (code in code.list){
    final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
    hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
    divq1 <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)])
    divq1_table <- divq1$data
    divq1_wild <- divq1_table[divq1_table$Group == "Wild",]
    n_wild <- nrow(divq1_wild)
    mean_wild <- mean(divq1_wild[,c("Value")])
    sd_wild <- sd(divq1_wild[,c("Value")])
    divq1_captive <- divq1_table[divq1_table$Group == "Captivity",]
    n_captive <- nrow(divq1_captive)
    mean_captive <- mean(divq1_captive[,c("Value")])
    sd_captive <- sd(divq1_captive[,c("Value")])
    summary_q1 <- rbind(summary_q1,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
  }
  colnames(summary_q1) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
  rownames(summary_q1) <- code.list
  write.table(summary_q1, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq1_",taxa,".tsv",sep=""))
#q=2
  summary_q2 <- c()
  for (code in code.list){
    final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
    hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
    divq2 <- div_test(final.table,qvalue=2,hierarchy=hierarchy[,c(1,5)])
    divq2_table <- divq2$data
    divq2_wild <- divq2_table[divq2_table$Group == "Wild",]
    n_wild <- nrow(divq2_wild)
    mean_wild <- mean(divq2_wild[,c("Value")])
    sd_wild <- sd(divq2_wild[,c("Value")])
    divq2_captive <- divq2_table[divq2_table$Group == "Captivity",]
    n_captive <- nrow(divq2_captive)
    mean_captive <- mean(divq2_captive[,c("Value")])
    sd_captive <- sd(divq2_captive[,c("Value")])
    summary_q2 <- rbind(summary_q2,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
  }
  colnames(summary_q2) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
  rownames(summary_q2) <- code.list
write.table(summary_q2, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq2_",taxa,".tsv",sep=""))

#Beta diversity
s <- div_part(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])

####Meta-analysis####
#q=0
summary_q0
summary_q0 <- as.data.frame(summary_q0)
meta_q0_ready <- tibble::rownames_to_column(summary_q0,"Author")
rownames(meta_q0_ready) <- meta_q0_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_q0_ready <- meta_q0_ready[sp_sorted,]
meta_q0.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                       data = meta_q0_ready,
                       studlab = paste(Author),
                       comb.fixed = FALSE,
                       comb.random = TRUE,
                       method.tau = "SJ",
                       hakn = TRUE,
                       prediction = TRUE,
                       sm = "SMD",
                       )
pdf(paste("~/reanalysis/Results/Diversity/",taxa,"/forest_q0",taxa,".pdf",sep=""),width=13, height=8)
forest(meta_q0.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()
saveRDS(meta_q0.raw, paste("~/reanalysis/Results/Diversity/",taxa,"/meta_q0_",taxa,".RData",sep=""))


#The box represents the effect size (in this case, OR) of each study. 
#The bigger the box means the study weighted more (i.e., bigger sample size).
#The diamond shape represents the pooled OR of all studies. We can see the diamond cross the vertical line OR = 1, which indicates no significance for the association as the diamond almost equalized in both sides. We can confirm this also from the 95% confidence interval that includes one and the p value > 0.05.
#For heterogeneity, we see that I2 = 0%, which means no heterogeneity is detected; the study is relatively homogenous (it is rare in the real study). To evaluate publication bias related to the meta-analysis of adverse events of arthralgia, we can use the metabias function from the R meta package

#q=1
summary_q1 <- as.data.frame(summary_q1)
meta_q1_ready <- tibble::rownames_to_column(summary_q1,"Author")
rownames(meta_q1_ready) <- meta_q1_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_q1_ready <- meta_q1_ready[sp_sorted,]
meta_q1.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                        data = meta_q1_ready,
                        studlab = paste(Author),
                        comb.fixed = FALSE,
                        comb.random = TRUE,
                        method.tau = "SJ",
                        hakn = TRUE,
                        prediction = TRUE,
                        sm = "SMD")
pdf(paste("~/reanalysis/Results/Diversity/",taxa,"/forest_q1",taxa,".pdf",sep=""),width=13, height=8)
forest(meta_q1.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()
saveRDS(meta_q1.raw,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/meta_q1_",taxa,".RData",sep=""))

#q=2
summary_q2
summary_q2 <- as.data.frame(summary_q2)
meta_q2_ready <- tibble::rownames_to_column(summary_q2,"Author")
rownames(meta_q2_ready) <- meta_q2_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_q2_ready <- meta_q2_ready[sp_sorted,]
meta_q2.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                        data = meta_q2_ready,
                        studlab = paste(Author),
                        comb.fixed = FALSE,
                        comb.random = TRUE,
                        method.tau = "SJ",
                        hakn = TRUE,
                        prediction = TRUE,
                        sm = "SMD")
pdf(paste("~/reanalysis/Results/Diversity/",taxa,"/forest_q2",taxa,".pdf",sep=""),width=13, height=8)
forest(meta_q2.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()
saveRDS(meta_q2.raw,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/meta_q2_",taxa,".RData",sep=""))


#q=0 subsample (Primates)
taxa = "Genus"
summary_q0.subset <- c()
sublist <- c("RHBR","PYNE","PAAN","PATR","GOGO")
for (code in sublist){
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divq0 <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  divq0_table <- divq0$data
  divq0_wild <- divq0_table[divq0_table$Group == "Wild",]
  n_wild <- nrow(divq0_wild)
  mean_wild <- mean(divq0_wild[,c("Value")])
  sd_wild <- sd(divq0_wild[,c("Value")])
  divq0_captive <- divq0_table[divq0_table$Group == "Captivity",]
  n_captive <- nrow(divq0_captive)
  mean_captive <- mean(divq0_captive[,c("Value")])
  sd_captive <- sd(divq0_captive[,c("Value")])
  summary_q0.subset <- rbind(summary_q0.subset,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_q0.subset) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_q0.subset) <- sublist
write.table(summary_q0.subset, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq0_subset",taxa,".tsv",sep=""))

summary_q0.subset <- as.data.frame(summary_q0.subset)
meta_q0sub_ready <- tibble::rownames_to_column(summary_q0.subset,"Author")
rownames(meta_q0sub_ready) <- meta_q0sub_ready$Author
sp_sorted <- c("RHBR","PYNE","PAAN","PATR","GOGO")
meta_q0sub_ready <- meta_q0sub_ready[sp_sorted,]
meta_q0_sub.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                            data = meta_q0sub_ready,
                            studlab = paste(Author),
                            comb.fixed = FALSE,
                            comb.random = TRUE,
                            method.tau = "SJ",
                            hakn = TRUE,
                            prediction = TRUE,
                            sm = "SMD",
)
pdf(paste("~/reanalysis/Results/Diversity/",taxa,"/forest_q0_sub_",taxa,".pdf",sep=""),width=13, height=8)
forest(meta_q0_sub.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()
saveRDS(meta_q0_sub.raw, paste("~/reanalysis/Results/Diversity/",taxa,"/meta_q0_sub_",taxa,".RData",sep=""))

#q1tree subset (Primates)
taxa = "Genus"
summary_q1tree_sub <- c()
sublist <- c("RHBR","PYNE","PAAN","PATR","GOGO")

for (code in sublist){
  print(code)
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  rownames(final.table) <- gsub(" ","_",rownames(final.table))
  rownames(final.table) <- gsub("[","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("]","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("(","-",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub(")","-",rownames(final.table),fixed = TRUE)
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divq1 <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)], tree=tree_filtered)
  divq1_table <- divq1$data
  divq1_wild <- divq1_table[divq1_table$Group == "Wild",]
  n_wild <- nrow(divq1_wild)
  mean_wild <- mean(divq1_wild[,c("Value")])
  sd_wild <- sd(divq1_wild[,c("Value")])
  divq1_captive <- divq1_table[divq1_table$Group == "Captivity",]
  n_captive <- nrow(divq1_captive)
  mean_captive <- mean(divq1_captive[,c("Value")])
  sd_captive <- sd(divq1_captive[,c("Value")])
  summary_q1tree_sub <- rbind(summary_q1tree_sub,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_q1tree_sub) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_q1tree_sub) <- sublist
write.table(summary_q1tree_sub, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq1tree_sub",taxa,".tsv",sep=""))

#q=1
summary_q1tree_sub <- as.data.frame(summary_q1tree_sub)
meta_q1tree_sub_ready <- tibble::rownames_to_column(summary_q1tree_sub,"Author")
rownames(meta_q1tree_sub_ready) <- meta_q1tree_sub_ready$Author
sp_sorted <- c("RHBR","PYNE","PAAN","PATR","GOGO")
meta_q1tree_sub_ready <- meta_q1tree_sub_ready[sp_sorted,]
meta_q1_sub_tree.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                             data = meta_q1tree_sub_ready,
                             studlab = paste(Author),
                             comb.fixed = FALSE,
                             comb.random = TRUE,
                             method.tau = "SJ",
                             hakn = TRUE,
                             prediction = TRUE,
                             sm = "SMD")
pdf(paste("~/reanalysis/Results/Diversity/",taxa,"/forest_q1_sub_tree",taxa,".pdf",sep=""),width=13, height=8)
forest(meta_q1_sub_tree.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()
saveRDS(meta_q1_sub_tree.raw,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/meta_q1_sub_tree",taxa,".RData",sep=""))



####Correlation####
#q=0
#open the meta_qX.raw
species_list <- meta_q0.raw$studlab
smd_value <- meta_q0.raw$TE
smd_matrix <- cbind(species_list,smd_value)
smd_table <- as.data.frame(smd_matrix)
smd_table <- tibble::column_to_rownames(smd_table, "species_list")
merged <- merge(summary_q0, smd_table, by = "row.names")
merged$n_sum <- merged$n_wild+merged$n_captive
n <- merged$n_sum
smd <- as.numeric(as.character(merged$smd_value))
cor.test(n,smd)
plot(n,smd)
abline(lm(smd ~ n, data = merged), col = "red")

#q=1
species_list <- meta_q1.raw$studlab
smd_value <- meta_q1.raw$TE
smd_matrix <- cbind(species_list,smd_value)
smd_table <- as.data.frame(smd_matrix)
smd_table <- tibble::column_to_rownames(smd_table, "species_list")
merged <- merge(summary_q1, smd_table, by = "row.names")
merged$n_sum <- merged$n_wild+merged$n_captive
n <- merged$n_sum
smd <- as.numeric(as.character(merged$smd_value))
cor.test(n,smd)
plot(n,smd)
abline(lm(smd ~ n, data = merged), col = "red")

#########PERMANOVA#############
#Permutest results
#Pairwise (dis)similarity computation based on beta diversities
taxa = "Genus"
#q=0
#code
permanovaq0_results <- c()
permutesq0_results <- c()
for (code in code.list){
filtered.g <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
##Filter the hierarchy
hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(filtered.g)),]
#Check whether both lists are identical
identical(sort(colnames(filtered.g)),sort(as.character(hierarchy[,1])))
##Filter the metadata
samples.kept <- colnames(filtered.g)
metadata.filtered <- metadata[rownames(metadata) %in% samples.kept,]
metadata.filtered1 <- tibble::rownames_to_column(metadata.filtered, "Code")
row.names(metadata.filtered1) <- metadata.filtered1$Datafiles
pairdis.q0 <- pair_dis(filtered.g,qvalue=0,hierarchy=hierarchy[,c(1,5)])
u0n <- pairdis.q0$L1_UqN
u0n.dist <- as.dist(u0n)
ps.disper.u0n.origin <- betadisper(u0n.dist, metadata.filtered1$Origin)
permutestq0 <- permutest(ps.disper.u0n.origin, pairwise = TRUE)
permutesq0_results <- append(permutesq0_results,permutestq0[1])
permanovaq0 <- adonis(u0n.dist ~ Origin, data =metadata.filtered1, permutations = 999)
print(code)
permanovaq0_results <- append(permanovaq0_results,permanovaq0[1])
}
names(permanovaq0_results) <- code.list
names(permutesq0_results) <- code.list
#write.table(permanovaq0_results,"permanova_q0.tsv",quote = FALSE, sep=";")
#write.table(permutesq0_results,"permutest_q0.tsv",quote = FALSE, sep=";")
saveRDS(permanovaq0_results,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/permanova_q0_",taxa,".RData",sep=""))
saveRDS(permutesq0_results,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/permutest_q0_",taxa,".RData",sep=""))
permanovaq0.genus <- readRDS("/Users/ostaizka/reanalysis/Results/Diversity/Genus/permanova_q0_Genus.RData")

#q=1
permanovaq1_results <- c()
permutesq1_results <- c()
for (code in code.list){
  filtered.g <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep="")) 
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(filtered.g)),]
  samples.kept <- colnames(filtered.g)
  metadata.filtered <- metadata[rownames(metadata) %in% samples.kept,]
  metadata.filtered1 <- tibble::rownames_to_column(metadata.filtered, "Code")
  row.names(metadata.filtered1) <- metadata.filtered1$Datafiles
  pairdis.q1 <- pair_dis(filtered.g,qvalue=1,hierarchy=hierarchy[,c(1,5)])
  u1n <- pairdis.q1$L1_UqN
  u1n.dist <- as.dist(u1n)
  ps.disper.u1n.origin <- betadisper(u1n.dist, metadata.filtered1$Origin)
  permutestq1 <- permutest(ps.disper.u1n.origin, pairwise = TRUE)
  permutesq1_results <- append(permutesq1_results,permutestq1[1])
  permanovaq1 <- adonis(u1n.dist ~ Origin, data =metadata.filtered1, permutations = 999)
  permanovaq1_results <- append(permanovaq1_results,permanovaq1[1])
}
names(permanovaq1_results) <- code.list
names(permutesq1_results) <- code.list
#write.table(permanovaq1_results,"permanova_q1.tsv",quote = FALSE, sep=";")
#write.table(permutesq1_results,"permutest_q1.tsv",quote = FALSE, sep=";")
saveRDS(permanovaq1_results,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/permanova_q1_",taxa,".RData",sep=""))
saveRDS(permutesq1_results,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/permutest_q1_",taxa,".RData",sep=""))


#####Shapiro-Wilk’s method for normality test#####
#Preliminary test to check independent t-test assumptions: normality(Shapiro) & homogeneity in variances(F-test)
taxa ="Phylum"
#wild
for (code in code.list){
filtered.g <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep="")) 
hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(filtered.g)),]
hierarchy_species <- hierarchy[,c(1,5)]
filtered.pecent <- tss(filtered.g)*100
filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
wild.capt <- merge(table.W,hierarchy_species,by="Samples")
wild.capt <- wild.capt[,-c(1)]
#all.columns <- colnames(wild.capt)
#taxa.all <- all.columns[-length(all.columns)]
wild <- wild.capt[wild.capt$Origin == "Wild",]
wild <- wild[,-ncol(wild)]
wild <- wild[,colSums(wild) > 0]
wild.taxa <- colnames(wild)
shapiro_wild <- c()
for (y in wild.taxa){
  res.shapiro.wild <- shapiro.test(wild[,y])
  shapiro_wild <- rbind(shapiro_wild,c(res.shapiro.wild[[1]],Pvalue=res.shapiro.wild[[2]]))
}
rownames(shapiro_wild) <- wild.taxa
assign(paste("shapiro_wild_",code,sep=""),get("shapiro_wild"))
write.table(get(paste("shapiro_wild_",code,sep="")),paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/shapiro_wild_",taxa,"_",code,".tsv",sep=""))
}
#captive
for (code in code.list){
  filtered.g <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep="")) 
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(filtered.g)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(filtered.g)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  captive <- wild.capt[wild.capt$Origin == "Captivity",]
  captive <- captive[,-ncol(captive)]
  captive <- captive[,colSums(captive) > 0]
  captive.taxa <- colnames(captive)
  shapiro_captive <- c()
  for (y in captive.taxa){
    res.shapiro.captive <- shapiro.test(captive[,y])
    shapiro_captive <- rbind(shapiro_captive,c(res.shapiro.captive[[1]],Pvalue=res.shapiro.captive[[2]]))
  }
  rownames(shapiro_captive) <- captive.taxa
  assign(paste("shapiro_captive_",code,sep=""),get("shapiro_captive"))
  write.table(get(paste("shapiro_captive_",code,sep="")),paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/shapiro_captive_",taxa,"_",code,".tsv",sep=""))
}

####Unpaired Two-Samples Wilcoxon Test at different level####
taxa = "Phylum"
for (code in code.list){
  filtered.g <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep="")) 
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(filtered.g)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(filtered.g)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  wild.capt.no <- wild.capt[,-ncol(wild.capt)]
  all.taxa <- colnames(wild.capt.no)
  Wilcox_result <- c()
  for (y in all.taxa){
  res.wilcox <- wilcox.test(wild.capt[,y] ~ Origin, data = wild.capt,
                            exact = FALSE, alternative = "less")
  Wilcox_result <- rbind(Wilcox_result,c(res.wilcox[[1]],Pvalue=res.wilcox[[3]]))
}
rownames(Wilcox_result) <- all.taxa
assign(paste("Wilcox_result_",code,sep=""),get("Wilcox_result"))
write.table(get(paste("Wilcox_result_",code,sep="")),paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/Wilcox_result_",taxa,"_",code,".tsv",sep=""))
}

#####Shared taxa#####
taxa="Genus"
Venn_result <- c()
for (code in code.list){
  filtered.g <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep="")) 
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(filtered.g)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(filtered.g)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  #all.columns <- colnames(wild.capt)
  #taxa.all <- all.columns[-length(all.columns)]
  wild <- wild.capt[wild.capt$Origin == "Wild",]
  wild <- wild[,-ncol(wild)]
  wild <- wild[,colSums(wild) > 0]
  wild.taxa <- colnames(wild)
  captive <- wild.capt[wild.capt$Origin == "Captivity",]
  captive <- captive[,-ncol(captive)]
  captive <- captive[,colSums(captive) > 0]
  captive.taxa <- colnames(captive)
  both.taxa <- intersect(wild.taxa,captive.taxa)
  wildonly.taxa <- wild.taxa[!wild.taxa %in% captive.taxa]
  captiveonly.taxa <- captive.taxa[! captive.taxa %in% wild.taxa]
  length(wildonly.taxa)
  length(both.taxa)
  length(captiveonly.taxa)
  shared <- cbind(length(wildonly.taxa),length(both.taxa),length(captiveonly.taxa))
  Venn_result <- rbind(Venn_result,shared)
}
colnames(Venn_result) <- c("wild_only", "both", "captive_only")
rownames(Venn_result) <- code.list
write.table(Venn_result,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/shared_",taxa,".tsv",sep=""))

### Beta diversity - incidence-based (Antton)

taxa = "Genus"
summary_beta_all <- c()
for (code in code.list){
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  final.table.inc <- to.incidence(final.table,hierarchy=hierarchy[,c(1,5)])
  divpart_q0 <- div_part(final.table.inc,qvalue=0)
  betadis_q0 <- beta_dis(beta=divpart_q0$Beta,qvalue=0,N=2)
  divpart_q1 <- div_part(final.table.inc,qvalue=1)
  betadis_q1 <- beta_dis(beta=divpart_q1$Beta,qvalue=1,N=2)
  summary_beta_all <- rbind(summary_beta_all,c(divpart_q0$Beta,divpart_q1$Beta,betadis_q0$CqN,betadis_q0$UqN,betadis_q1$CqN,betadis_q1$VqN,betadis_q1$SqN))
 }
colnames(summary_beta_all) <- c("beta_q0", "beta_q1","CqN_q0","UqN_q0","CqN_q1","VqN_q1","SqN_q1")
rownames(summary_beta_all) <- code.list
write.table(summary_beta_all,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/summary_beta_all",taxa,".tsv",sep=""))

####Species diversity analysis based on abundace####
taxa = "Genus"
summary_q0_all <- c()
for (code in code.list){
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divq0 <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  divq0_table <- divq0$data
  n_all <- nrow(divq0_table)
  mean_all <- mean(divq0_table[,c("Value")])
  sd_all <- sd(divq0_table[,c("Value")])
  summary_q0_all <- rbind(summary_q0_all,c(n_all,mean_all,sd_all))
}
colnames(summary_q0_all) <- c("n_q0", "mean_q0", "sd_q0")
rownames(summary_q0_all) <- code.list
#write.table(summary_q0_all, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq0_all_",taxa,".tsv",sep=""))

taxa = "Genus"
summary_q1_all <- c()
for (code in code.list){
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divq1 <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)])
  divq1_table <- divq1$data
  n_all <- nrow(divq1_table)
  mean_all <- mean(divq1_table[,c("Value")])
  sd_all <- sd(divq1_table[,c("Value")])
  summary_q1_all <- rbind(summary_q1_all,c(n_all,mean_all,sd_all))
}
colnames(summary_q1_all) <- c("n_q1", "mean_q1", "sd_q1")
rownames(summary_q1_all) <- code.list
#write.table(summary_q1_all, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq1_all_",taxa,".tsv",sep=""))
summary_all <- merge(summary_q0_all,summary_q1_all, by= "row.names")
summary_all <- tibble::column_to_rownames(summary_all, "Row.names")
write.table(summary_all, paste("~/reanalysis/Results/Diversity/",taxa,"/summary_all_",taxa,".tsv",sep=""))


#### Merging all tables in one to analyse the differences between species ####
#all individuals of different species in one table
#Genus level
setwd("~/reanalysis")
  #Load files at genus level
  files.g = list.files(path="table",pattern="countfiltered_Genus_*")
  #names.g = gsub("countfiltered_Genus_","",files.g)
  #names.g = gsub(".tsv","",names.g)
  filelist.g = lapply(paste("table",files.g,sep="/"), function(x)read.table(x, header=T))
  rownameslist = lapply(filelist.g, function(x) rownames(x))
  for (l in c(1:length(filelist.g))){
    filelist.g[[l]]$taxa <- rownames(filelist.g[[l]])
  }
  #Note some erroneous files can prevent all files from loading (remove erroneous ones)
  #Transform list to matrix
  filelist.g %>% purrr::reduce(full_join, by = "taxa") -> count.table.all.g
  rownames(count.table.all.g) <- count.table.all.g[,"taxa"]
  count.table.all.g <- count.table.all.g[,colnames(count.table.all.g) != "taxa"]
  #Name columns
  #colnames(count.table.all.g) <- names.g
  #Replace NAs by 0s
  count.table.all.g[is.na(count.table.all.g)] <- 0
  write.table(count.table.all.g, "/Users/ostaizka/reanalysis/table/count_Genus_all.tsv")

#Family level
setwd("~/reanalysis")
  #Load files at genus level
  files.f = list.files(path="table",pattern="countfiltered_Family_*")
  filelist.f = lapply(paste("table",files.f,sep="/"), function(x)read.table(x, header=T))
  rownameslist = lapply(filelist.f, function(x) rownames(x))
  for (l in c(1:length(filelist.f))){
    filelist.f[[l]]$taxa <- rownames(filelist.f[[l]])
  }
  #Note some erroneous files can prevent all files from loading (remove erroneous ones)
  #Transform list to matrix
  filelist.f %>% purrr::reduce(full_join, by = "taxa") -> count.table.all.f
  rownames(count.table.all.f) <- count.table.all.f[,"taxa"]
  count.table.all.f <- count.table.all.f[,colnames(count.table.all.f) != "taxa"]
  #Name columns
  #colnames(count.table.all.f) <- names.f
  #Replace NAs by 0s
  count.table.all.f[is.na(count.table.all.f)] <- 0
  write.table(count.table.all.f,"/Users/ostaizka/reanalysis/table/count_Family_all.tsv")

#Phylum level
setwd("~/reanalysis")
  #Load files at genus level
  files.p = list.files(path="table",pattern="countfiltered_Phylum_*")
  filelist.p = lapply(paste("table",files.p,sep="/"), function(x)read.table(x, header=T))
  rownameslist = lapply(filelist.p, function(x) rownames(x))
  for (l in c(1:length(filelist.p))){
    filelist.p[[l]]$taxa <- rownames(filelist.p[[l]])
  }
  #Note some erroneous files can prevent all files from loading (remove erroneous ones)
  #Transform list to matrix
  filelist.p %>% purrr::reduce(full_join, by = "taxa") -> count.table.all.p
  rownames(count.table.all.p) <- count.table.all.p[,"taxa"]
  count.table.all.p <- count.table.all.p[,colnames(count.table.all.p) != "taxa"]
  #Name columns
  #colnames(count.table.all.p) <- names.p
  #Replace NAs by 0s
  count.table.all.p[is.na(count.table.all.p)] <- 0
  write.table(count.table.all.p,"/Users/ostaizka/reanalysis/table/count_Phylum_all.tsv")

#Divtest
hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(count.table.all.g)),]
divq0.all <- div_test(count.table.all.g,qvalue=0,hierarchy=hierarchy[,c(1,3)])

#Neutral pairdis
pairdis.q0.all <- pair_dis(count.table.all.g,qvalue=0,hierarchy=hierarchy[,c(1,3)])
saveRDS(pairdis.q0.all,"/Users/ostaizka/reanalysis/Results/Diversity/pairdis_q0_all.RData")
pairdis.q1.all <- pair_dis(count.table.all.g,qvalue=1,hierarchy=hierarchy[,c(1,3)])
saveRDS(pairdis.q1.all,"/Users/ostaizka/reanalysis/Results/Diversity/pairdis_q1_all.RData")

#pair_dis_plot(pairdis.q0.all$L1_UqN,hierarchy=hierarchy[,c(1,5)])
dis_nmds(pairdis.q1.all_L1_UqN,hierarchy=hierarchy[hierarchy$Samples != "D12719",c(1,19)], centroids=TRUE, runs=200)
saveRDS(pairdis.q0.all,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/permanova_q1_",taxa,".RData",sep=""))

## permanova
##Filter the metadata
samples.kept <- colnames(count.table.all.g)
metadata.filtered <- metadata[rownames(metadata) %in% samples.kept,]
#metadata.filtered1 <- tibble::rownames_to_column(metadata.filtered, "Code")
#row.names(metadata.filtered1) <- metadata.filtered1$Datafiles
#q0
u0n <- pairdis.q0.all$L1_UqN
u0n.dist <- as.dist(u0n)
ps.disper.u0n.origin <- betadisper(u0n.dist, metadata.filtered$Species)
permutestq0 <- permutest(ps.disper.u0n.origin, pairwise = TRUE)
permanovaq0 <- adonis(u0n.dist ~ Species, data =metadata.filtered, permutations = 999)
#not working u0n.pairwise <- calc_pairwise_permanovas(u0n.dist, metadata.filtered, "Species")
adonis(u0n.dist ~ Origin + Species + Sample_type + Type + S_region + F_Primer, data =metadata.filtered, permutations = 999)
colnames(metadata.filtered) <- gsub("16","",colnames(metadata.filtered))
colnames(metadata.filtered) <- gsub(" ","",colnames(metadata.filtered))
colnames(metadata.filtered) <- gsub("-","_",colnames(metadata.filtered))
#q1
u1n <- pairdis.q1.all$L1_UqN
u1n.dist <- as.dist(u1n)
ps.disper.u1n.origin <- betadisper(u1n.dist, metadata.filtered$Origin)
permutestq1 <- permutest(ps.disper.u1n.origin, pairwise = TRUE)
permanovaq1 <- adonis(u1n.dist ~ Origin, data =metadata.filtered, permutations = 999)
#not working u1n.pairwise <- calc_pairwise_permanovas(u1n.dist, metadata.filtered, "Species")

###### Phylogenetic ######
capwild.tree <- read.tree("~/reanalysis/genustree.tre")
#fishing.bats.tree <- force.ultrametric(fishing.bats.tree, method = "extend")
#tree_filtered <- match_data(final.table,capwild.tree,output="tree")

taxa = "Genus"
summary_q0tree <- c()
for (code in code.list){
  print(code)
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  rownames(final.table) <- gsub(" ","_",rownames(final.table))
  rownames(final.table) <- gsub("[","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("]","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("(","-",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub(")","-",rownames(final.table),fixed = TRUE)
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divq0 <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)], tree=tree_filtered)
  divq0_table <- divq0$data
  divq0_wild <- divq0_table[divq0_table$Group == "Wild",]
  n_wild <- nrow(divq0_wild)
  mean_wild <- mean(divq0_wild[,c("Value")])
  sd_wild <- sd(divq0_wild[,c("Value")])
  divq0_captive <- divq0_table[divq0_table$Group == "Captivity",]
  n_captive <- nrow(divq0_captive)
  mean_captive <- mean(divq0_captive[,c("Value")])
  sd_captive <- sd(divq0_captive[,c("Value")])
  summary_q0tree <- rbind(summary_q0tree,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_q0tree) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_q0tree) <- code.list
write.table(summary_q0tree, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq0tree_",taxa,".tsv",sep=""))

taxa = "Genus"
summary_q1tree <- c()
for (code in code.list){
  print(code)
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  rownames(final.table) <- gsub(" ","_",rownames(final.table))
  rownames(final.table) <- gsub("[","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("]","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("(","-",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub(")","-",rownames(final.table),fixed = TRUE)
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divq1 <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)], tree=tree_filtered)
  divq1_table <- divq1$data
  divq1_wild <- divq1_table[divq1_table$Group == "Wild",]
  n_wild <- nrow(divq1_wild)
  mean_wild <- mean(divq1_wild[,c("Value")])
  sd_wild <- sd(divq1_wild[,c("Value")])
  divq1_captive <- divq1_table[divq1_table$Group == "Captivity",]
  n_captive <- nrow(divq1_captive)
  mean_captive <- mean(divq1_captive[,c("Value")])
  sd_captive <- sd(divq1_captive[,c("Value")])
  summary_q1tree <- rbind(summary_q1tree,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_q1tree) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_q1tree) <- code.list
write.table(summary_q1tree, paste("~/reanalysis/Results/Diversity/",taxa,"/summaryq1tree_",taxa,".tsv",sep=""))

####Meta-analysis####
#q=0
summary_q0tree
summary_q0tree <- as.data.frame(summary_q0tree)
meta_q0tree_ready <- tibble::rownames_to_column(summary_q0tree,"Author")
rownames(meta_q0tree_ready) <- meta_q0tree_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_q0tree_ready <- meta_q0tree_ready[sp_sorted,]
meta_q0_tree.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                             data = meta_q0tree_ready,
                             studlab = paste(Author),
                             comb.fixed = FALSE,
                             comb.random = TRUE,
                             method.tau = "SJ",
                             hakn = TRUE,
                             prediction = TRUE,
                             sm = "SMD",
)
pdf(paste("~/reanalysis/Results/Diversity/",taxa,"/forest_q0_tree",taxa,".pdf",sep=""),width=13, height=8)
forest(meta_q0_tree.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()
saveRDS(meta_q0_tree.raw, paste("~/reanalysis/Results/Diversity/",taxa,"/meta_q0_tree",taxa,".RData",sep=""))
meta_q1_tree <- readRDS("~/reanalysis/Results/Diversity/Genus/meta_q1_treeGenus.RData")
meta_q0 <- readRDS("~/reanalysis/Results/Diversity/Genus/meta_q0_Genus.RData")

#The box represents the effect size (in this case, OR) of each study. 
#The bigger the box means the study weighted more (i.e., bigger sample size).
#The diamond shape represents the pooled OR of all studies. We can see the diamond cross the vertical line OR = 1, which indicates no significance for the association as the diamond almost equalized in both sides. We can confirm this also from the 95% confidence interval that includes one and the p value > 0.05.
#For heterogeneity, we see that I2 = 0%, which means no heterogeneity is detected; the study is relatively homogenous (it is rare in the real study). To evaluate publication bias related to the meta-analysis of adverse events of arthralgia, we can use the metabias function from the R meta package

#Check metaanalysis_example.R for editing
## Detecting outliers & influential cases
#library(dmetar)
#meta_q0.raw$lower.random
#meta_q0.raw$upper.random
#find.outliers(meta_q0.raw)

#q=1
summary_q1tree <- as.data.frame(summary_q1tree)
meta_q1tree_ready <- tibble::rownames_to_column(summary_q1tree,"Author")
rownames(meta_q1tree_ready) <- meta_q1tree_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_q1tree_ready <- meta_q1tree_ready[sp_sorted,]
meta_q1_tree.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                             data = meta_q1tree_ready,
                             studlab = paste(Author),
                             comb.fixed = FALSE,
                             comb.random = TRUE,
                             method.tau = "SJ",
                             hakn = TRUE,
                             prediction = TRUE,
                             sm = "SMD")
pdf(paste("~/reanalysis/Results/Diversity/",taxa,"/forest_q1_tree",taxa,".pdf",sep=""),width=13, height=8)
forest(meta_q1_tree.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()
saveRDS(meta_q1_tree.raw,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/meta_q1_tree",taxa,".RData",sep=""))

#Tree
count.table.all.g <- read.table("/Users/ostaizka/reanalysis/table/count_Genus_all.tsv")
rownames(count.table.all.g) <- gsub(" ","_",rownames(count.table.all.g))
rownames(count.table.all.g) <- gsub("[","",rownames(count.table.all.g),fixed = TRUE)
rownames(count.table.all.g) <- gsub("]","",rownames(count.table.all.g),fixed = TRUE)
rownames(count.table.all.g) <- gsub("(","-",rownames(count.table.all.g),fixed = TRUE)
rownames(count.table.all.g) <- gsub(")","-",rownames(count.table.all.g),fixed = TRUE)
hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(count.table.all.g)),]
capwild.tree <- read.tree("~/reanalysis/genustree.tre")
tree_filtered <- match_data(count.table.all.g,capwild.tree,output="tree")

#Phylogenetic divtest species
divq0.all.tree <- div_test(count.table.all.g,qvalue=0,hierarchy=hierarchy[,c(1,3)],tree = tree_filtered)
saveRDS(divq0.all.tree,"/Users/ostaizka/reanalysis/Results/Diversity/divq0.all.tree.RData")
divq1.all.tree <- div_test(count.table.all.g,qvalue=1,hierarchy=hierarchy[,c(1,3)],tree = tree_filtered)
saveRDS(divq1.all.tree,"/Users/ostaizka/reanalysis/Results/Diversity/divq1.all.tree.RData")

#Correlation: diversity vs sequencing depth
diversity1.tree <- divq1.all.tree$data
seq.depth <- colSums(count.table.all.g)
seq.depth1 <- as.data.frame(seq.depth)
seq.depth1 <- tibble::rownames_to_column(seq.depth1, "Sample")
table <- merge(diversity1.tree,seq.depth1, by="Sample")
div <- table$Value
seq <- table$seq.depth
smd <- as.numeric(as.character(merged$smd_value))
cor.test(div,seq)
plot(div,seq)
abline(lm(seq ~ div, data = table), col = "red")





#Phylogenetic pairdis
pairdis.q0.all.tree <- pair_dis(count.table.all.g,qvalue=0,hierarchy=hierarchy[,c(1,3)],tree = tree_filtered)
saveRDS(pairdis.q0.all.tree,"/Users/ostaizka/reanalysis/Results/Diversity/pairdisq0_all_tree1.RData")
pairdis.q1.all.tree <- pair_dis(count.table.all.g,qvalue=1,hierarchy=hierarchy[,c(1,3)],tree = tree_filtered)
saveRDS(pairdis.q1.all.tree,"/Users/ostaizka/reanalysis/Results/Diversity/pairdisq1_all_tree1.RData")
dis_nmds(pairdis.q1.all.tree$L1_UqN,hierarchy=hierarchy[,c(1,3)], centroids=TRUE)
dis_nmds(pairdis.q1.all.tree$L1_UqN,hierarchy=hierarchy[hierarchy$Samples != "D12719",c(1,19)], centroids=TRUE, runs=200)

##Phylogenetic permanova
##Filter the metadata
samples.kept <- colnames(count.table.all.g)
metadata.filtered <- metadata[rownames(metadata) %in% samples.kept,]
metadata.filtered1 <- tibble::rownames_to_column(metadata.filtered, "Code")
row.names(metadata.filtered1) <- metadata.filtered1$Datafiles
#q0
u0n.tree <- pairdis.q0.all.tree$L1_UqN
u0n.dist.tree <- as.dist(u0n.tree)
ps.disper.u0n.tree <- betadisper(u0n.dist.tree, metadata.filtered$Species)
permutestq0.tree <- permutest(ps.disper.u0n.tree, pairwise = TRUE)
permanovaq0.tree <- adonis(u0n.dist.tree ~ Species, data =metadata.filtered, permutations = 999)
#not working u0n.pairwise <- calc_pairwise_permanovas(u0n.dist.tree, metadata.filtered, "Study")

#q1
u1n.tree <- pairdis.q1.all.tree$L1_UqN
u1n.dist.tree <- as.dist(u1n.tree)
ps.disper.u1n.tree <- betadisper(u1n.dist.tree, metadata.filtered$Species)
permutestq1.tree <- permutest(ps.disper.u1n.tree, pairwise = TRUE)
permanovaq1.tree <- adonis(u1n.dist.tree ~ Origin * Species, data =metadata.filtered, permutations = 999)
permanovaq1.tree <- adonis(u1n.dist.tree ~ Origin + Species, strata=metadata.filtered$Species, data =metadata.filtered, permutations = 999)
permanovaq1.tree <- adonis(u1n.dist.tree ~ Origin * Species, strata=metadata.filtered$Species, data =metadata.filtered, permutations = 999)

#not working u1n.pairwise <- calc_pairwise_permanovas(u1n.dist, metadata.filtered1, "Species")


##### Did captivity enrich some taxa?####
taxa="Genus"
incid_result_table <- c()
incid_ratio_table <- c()
count.table.all <- read.table(paste("/Users/ostaizka/reanalysis/table/count_",taxa,"_all.tsv",sep = ""))
hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(count.table.all)),]

for (code in code.list){
  filtered.table.w <- count.table.all[,hierarchy[(hierarchy$Species == code) & (hierarchy$Origin == "Wild"),"Samples"]]
  wild <- to.incidence(filtered.table.w)
  filtered.table.c <- count.table.all[,hierarchy[(hierarchy$Species == code) & (hierarchy$Origin == "Captivity"),"Samples"]]
  captive <- to.incidence(filtered.table.c)
  incid_result <- cbind(wild,captive)
  incid_result1 <- as.data.frame(cbind(wild=incid_result[,1]*100/ncol(filtered.table.w),captive=incid_result[,2]*100/ncol(filtered.table.c)))
  #Ratio
  incid_ratio_vector <- incid_result1$captive / incid_result1$wild
  incid_ratio_table <- cbind(incid_ratio_table,incid_ratio_vector)
  #Data
  colnames(incid_result1) <- c(paste(code,"wild",sep="_"),paste(code,"captive",sep="_"))
  if (is.null(nrow(incid_result_table))){
  incid_result_table <- incid_result1
  }else{
  incid_result_table <- cbind(incid_result_table,incid_result1)
  }
}
rownames(incid_ratio_table) <- rownames(incid_result_table)
colnames(incid_ratio_table) <- code.list
#write.table(incid_ratio_table,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/incid_ratio_",taxa,".tsv",sep=""))

incid_ratio_table2 <- incid_ratio_table
incid_ratio_table2[incid_ratio_table2 > 1] <- 1
incid_ratio_table2[incid_ratio_table2 < 1] <- 0
captive.enriched <- rowSums(incid_ratio_table2,na.rm=TRUE) / rowSums(!is.na(incid_ratio_table2)) * 100
captive.enriched[captive.enriched == 100]








##########################################
#####Shared taxa#####
taxa="Phylum"
both_shared <- c()
wild_only <- c()
captive_only <- c()

for (code in code.list){
  filtered.g <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep="")) 
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(filtered.g)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(filtered.g)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  #all.columns <- colnames(wild.capt)
  #taxa.all <- all.columns[-length(all.columns)]
  wild <- wild.capt[wild.capt$Origin == "Wild",]
  wild <- wild[,-ncol(wild)]
  wild <- wild[,colSums(wild) > 0]
  wild.taxa <- colnames(wild)
  captive <- wild.capt[wild.capt$Origin == "Captivity",]
  captive <- captive[,-ncol(captive)]
  captive <- captive[,colSums(captive) > 0]
  captive.taxa <- colnames(captive)
  both.taxa <- intersect(wild.taxa,captive.taxa)
  both_shared <- rbind(both_shared,c(both.taxa, code))
  wildonly.taxa <- wild.taxa[!wild.taxa %in% captive.taxa]
  wild_only <- rbind(wild_only,c(wildonly.taxa,code))
  captiveonly.taxa <- captive.taxa[! captive.taxa %in% wild.taxa]
  captive_only <- rbind(captive_only,c(captiveonly.taxa,code))
}
  
  
  
  rownames(q) <- code
  filtered.f <- final.table.f[,colnames(final.table.f) %in% samples.kept]
  
  
  
  assign(paste("filtered.f.",code,sep=""),get("filtered.f"))
#Shared taxa list
for (code in code.list){
  filtered.g <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep="")) 
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(filtered.g)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(filtered.g)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  #all.columns <- colnames(wild.capt)
  #taxa.all <- all.columns[-length(all.columns)]
  wild <- wild.capt[wild.capt$Origin == "Wild",]
  wild <- wild[,-ncol(wild)]
  wild <- wild[,colSums(wild) > 0]
  wild.taxa <- colnames(wild)
  captive <- wild.capt[wild.capt$Origin == "Captivity",]
  captive <- captive[,-ncol(captive)]
  captive <- captive[,colSums(captive) > 0]
  captive.taxa <- colnames(captive)
  ItemsList <- venn(list(wild.taxa, captive.taxa), show.plot = FALSE)
  saveRDS(taxa_shared,paste("/Users/ostaizka/reanalysis/Results/Diversity/",taxa,"/shared_list",taxa,"_",code,".RData",sep=""))
  }          

#attr(x = ItemsList, "intersections")
# get specific intersections
#attr(x = ItemsList, "intersections")$`A:B`










####Diversity analysis based on abundace of each species####
#q=0
taxa = "Genus"
summary_all_q0 <- c()
for (code in code.list){
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divq0 <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  divq0_all <- divq0$data
  n_all <- nrow(divq0_all)
  mean_all <- mean(divq0_all[,c("Value")])
  sd_all <- sd(divq0_all[,c("Value")])
  summary_all_q0 <- rbind(summary_all_q0,c(n_all,mean_all,sd_all))
}
colnames(summary_all_q0) <- c("N", "Mean", "SD")
rownames(summary_all_q0) <- code.list
write.table(summary_all_q0, paste("~/reanalysis/Results/Diversity/",taxa,"/all_summaryq0_",taxa,".tsv",sep=""))
#q=1
taxa = "Genus"
summary_all_q1 <- c()
for (code in code.list){
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))  
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divq1 <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)])
  divq1_all <- divq1$data
  n_all <- nrow(divq1_all)
  mean_all <- mean(divq1_all[,c("Value")])
  sd_all <- sd(divq1_all[,c("Value")])
  summary_all_q1 <- rbind(summary_all_q1,c(n_all,mean_all,sd_all))
}
colnames(summary_all_q1) <- c("N", "Mean", "SD")
rownames(summary_all_q1) <- code.list
write.table(summary_all_q1, paste("~/reanalysis/Results/Diversity/",taxa,"/all_summaryq1_",taxa,".tsv",sep=""))








###not finishhed, pairdis by alltogether
code.list[1]
table1
code.list[2]
table2
all_together <- merge(table1, table1, by = "row.names")
code.list[3]
all_together <- merge(all_together, table3, by = "row.names")
code.list[4]
all_together <- merge(all_together, table4, by = "row.names")



for (code in code.list){
  final.table <- read.table(paste("/Users/ostaizka/reanalysis/table/countfiltered_",taxa,"_",code,".tsv",sep=""))


  
  
  
