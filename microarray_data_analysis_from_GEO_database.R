#                      
#             ****Microarray data analysis pipeline--Affymetrix****
#                       ****from CEL files to DEGs****

# Elif Duz
# 28.03.25

# libraries
library(affy)
library(GEOquery)
library(tidyverse)
library(readxl)
library(limma)
library(readxl)
library(gprofiler2)
library(ggplot2)
library(ggpubr)

gseid= "GSE42669"  #----> define GSEID in here
gse <- getGEO(gseid, GSEMatrix = TRUE)

#for multiple GPL file: --------> check the gse file
gse=gse[[paste0(gseid, "-GPL96_series_matrix.txt.gz")]]

#option 1:
# RMA normalization
# download data from CEL files 
# raw data is used for rma normalization

getGEOSuppFiles(gseid)
folder_name <- paste0(gseid, "data")# ---> create a folder same name with data
untar(paste0(gseid, "/", gseid, "_RAW.tar"), exdir = folder_name)
raw.data <- ReadAffy(celfile.path = folder_name)
normalized.data <- rma(raw.data)
normalized.expr <- as.data.frame(exprs(normalized.data))
#------------------------------------------------------
#option 2:
#if you use the data directly from GEO database wo normalization:
#you can use this option for illumina datasets
normalized.expr <- data.frame(gse[[paste0(gseid, "_series_matrix.txt.gz")]]@assayData[["exprs"]])
#if you want to apply log2 transformation: --------> check the data for log transformation
normalized.expr= log2(normalized.expr)
#------------------------------------------------------

feature.data <- gse[[paste0(gseid, "_series_matrix.txt.gz")]]@featureData@data
feature.data <- feature.data[,c(1,10)] #----> check the column names (AccessionID and Gene Symbol)
feature.data$ID= as.character(feature.data$ID) #convert the numerical values into character vector

#option 3
# if you can not connect GEOdatabases you can download CEL file from GEO webpage store zip file in a folder with same name the GSE data.
# get probeID from geo database GPL files

folder_name <- paste0(gseid, "data")# ---> create a folder same name with data
untar(paste0(gseid, "/", gseid, "_RAW.tar"), exdir = folder_name)
raw.data <- ReadAffy(celfile.path = folder_name)
#go back to RMA normalization section and create normalized.expr data

feature.data= read.delim2("GPL570-55999.txt", comment.char="#")# ---> get the feature data from GEO database manually
feature.data= gse@featureData@data # get the feature data from gse data
feature.data <- feature.data[,c(1,10)] #----> check the column names (AccessionID and Gene Symbol)
feature.data$ID= as.character(feature.data$ID) #convert the numerical values into character vector
#------------------------------------------------------

# combine probe ID and gene symbols
combined <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')

data= data.frame(genesym= combined[, length(combined)], combined[,2:(length(combined)-1)]) #---> first and last columns are removed

# option 1:
# if data have multiple gene symbols for one gene;
data2= separate_rows(data, c(genesym), sep="///", convert=TRUE) # either this way
data2= separate_rows(data, genesym, sep="\\|", convert=TRUE) # or that one
#------------------------------------------------------
# option 2:
# if gene symbols should seperated from multiple identifiers;
genesymb= data.frame(genesym=feature.data$gene_assignment)
genesymb=genesymb %>% separate(genesym, into = c("Column1", "genesym", "Column3"), sep = "//")
data2= data.frame(genesym= genesymb$genesym, data[,2:length(data)])
#------------------------------------------------------

data2$genesym=trimws(toupper(data2$genesym)) #trim the gene symbols and use upper characters

# if column names include unnecessary characters ----> check the column names for consistency
colnames(data2)= sub("_.*", "", colnames(data2))

# unique the data based on gene symbols, for microarray select the row with max expression
df_unique <- data2 %>%
rowwise() %>%
  mutate(Total = sum(c_across(starts_with("GSM")))) %>% 
  group_by(genesym) %>%
  filter(Total == max(Total)) %>%
  ungroup() %>%
  dplyr::select(-Total) 
# if this code do not work try this:
# df_unique=data2 %>%group_by(genesym) %>%summarise_all(max)

df_unique <- df_unique[!(is.na(df_unique$genesym) | df_unique$genesym == ""), ] #remove the nameless row
df_unique= unique(df_unique)

# prepare metadata in excel file name the sheet wrt gseid and upload it
# order the data, disease vs control
metadata <- read_excel("celllines.xlsx", sheet = gseid)

# if you want to remove some information from metadata:
# metadata= na.omit(metadata)
# metadata= metadata[metadata$resistance!= "-",]

# type column represents the disease status, arrange it based on status----> check the disease status column for comparison
df_unique= data.frame(genesym= toupper(df_unique$genesym),df_unique[, metadata[metadata$type=="TMZ2",]$Accession],df_unique[, metadata[metadata$type=="TMZ1",]$Accession],df_unique[, metadata[metadata$type=="control",]$Accession])
df_unique <- df_unique[rowSums(df_unique != 0, na.rm = TRUE) > 0, ] #at least one row (>0) that total expression is different than 0

rnames= df_unique$genesym
df_unique$genesym <- NULL
rownames(df_unique) <- rnames

# perform PCA for outlier detection
sigmatrix= data.frame(status=factor(metadata$type),rownames=colnames(df_unique))

pcaPlot <- function(data, Key, title){
  require(ggplot2)
  df <- as.matrix(data)
  # remove the NA values
  df[is.infinite(df)] <- NA
  df <- na.omit(df)  
  # remove the rows with 0 in all samples
  df <- df[rowSums(df != 0, na.rm = TRUE) > 0, ]
  # remove the 0 varience
  df <- df[apply(df, 1, var, na.rm = TRUE) > 0, ]
  
  pcaanal <- prcomp(t(df))
  
  pcaplot <- data.frame(PC1 = pcaanal$x[,1], 
                        PC2 = pcaanal$x[,2], 
                        Sample = rownames(pcaanal$x),
                        Key = Key)  

  ggplot(data = pcaplot, aes(x = PC1, y = PC2, color = Key)) +
    geom_point(size = 3) +
    geom_text(aes(label = Sample), hjust = 0.6, vjust = 0, size = 4) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_bw() +
    theme(legend.title = element_blank())  
}
pcaPlot(df_unique,sigmatrix$status,paste0(gseid,"_PCA"))
# ----> to check the cumulative proportions:
df=as.matrix(df_unique)
pcaanal=prcomp(t(df))
summary(pcaanal)

# save the PCA plot
png(paste0(gseid, "_PCA.png"), width = 3000, height = 1500, res = 300)
pcaPlot(df_unique,sigmatrix$status,paste0(gseid,"_PCA"))
dev.off()

# save the normalized data
save(metadata, df_unique, file= paste0(gseid, "_normalized_unique.RData"))

# option 
# if any of the samples removed from the dataset, rearrange both data and metadata:
metadata= metadata[metadata$type!= "TMZ2",]
rnames= rownames(df_unique)
df_unique= df_unique[, metadata$Accession]
rownames(df_unique)= rnames
#------------------------------------------------------

# DEG analysis with limma 
# option 1:unpaired samples
sigmatrix= data.frame(status=factor(metadata$type),rownames=colnames(df_unique))
disease_control_list <- sigmatrix

group <- factor(metadata$type, levels = c("TMZ1", "control")) #-----> check the levels and rearrange them (disease vs control)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(df_unique, design)
contrast_matrix <- makeContrasts(TMZ1_vs_control = TMZ1 - control, levels = design) #----> rearrange this part, disease vs control
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
result <- topTable(fit2, adjust = "BH", number = Inf)
#------------------------------------------------------

# option 2: paired samples

samples <- data.frame(
  patient <- factor(c(1:3 ,1:3)),  # match the samples first three samples are recurrent, last three samples from primary 
  group <- factor(metadat$type2))

design= model.matrix(~ patient + group, data=samples)

fit <- lmFit(mesenchymaldata, design)
fit <- eBayes(fit)
results <- topTable(fit, coef="grouprecurrent", adjust="BH", number=Inf) #----> rearrange this part (grouprecurrent) recurrent is diseased group 
#------------------------------------------------------

DEGs= result[result$P.Value <0.05,]#----> define p-value for significance, you can choose adjusted or non-adjusted p-values based on number of DEGs
result$FC= 2^(result$logFC)

write.table(result, file= paste0(gseid, "_DEGs.xls"))

results2 <- data.frame(
  gene = rownames(result),      # Gene names
  logFC = result$logFC,         # Log2 fold change
  pvalue = result$P.Value       # P-values
)

results2$negLogP <- -log10(results2$pvalue)

results2$color <- ifelse(results2$logFC > 0.5849625 & results2$pvalue< 0.05, "upregulated", 
                         ifelse(results2$logFC < -0.5849625 & results2$pvalue < 0.05, "downregulated", "not significant"))
top_genes <- results2[order(results2$pvalue), ][1:1000, ] #---> you can select the number of genes shown in graph

# Volcano plot
png(paste0(gseid, "_volcano.png"), width = 2000, height = 1300, res = 300) #-------> check the graph and arrange the size

ggplot(results2, aes(x = logFC, y = negLogP, color = color)) +
  geom_point(alpha = 0.6, size = 3) +  
  scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "not significant" = "grey")) + 
  labs(title = paste0(gseid," Volcano Plot"),
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = c(-0.5849625, 0.5849625), linetype = "dashed", color = "blue") +  
  geom_text(data = top_genes, aes(label = gene), vjust = 1.5, size = 3, check_overlap = TRUE)  

dev.off()

# enrichment analysis
upgenes= na.omit(rownames(result[result$P.Val<0.05 & result$FC>1.5,]))
downgenes= na.omit(rownames(result[result$P.Val<0.05 & result$FC<0.6666,]))

# convert gene symbols to ENSID for enrichment analysis
ENSID= gconvert(as.vector(upgenes))$target
module= gost(as.vector(ENSID), evcodes=TRUE)$result
module2up= module[(module$source != "TF") & (module$source != "HPA") & (module$source != "MIRNA") ,]
GOtermsup= module2up[grep("GO:", module2up$source),]
othersup= module2up[grep("GO:", module2up$source, invert=TRUE),]
moduleup= module2up[,c(9,10,11,3,6,16)]

ENSID= gconvert(as.vector(downgenes))$target
module= gost(as.vector(ENSID), evcodes=TRUE)$result
module2down= module[(module$source != "TF") & (module$source != "HPA") & (module$source != "MIRNA") ,]
GOtermsdown= module2down[grep("GO:", module2down$source),]
othersdown= module2down[grep("GO:", module2down$source, invert=TRUE),]
moduledown= module2down[,c(9,10,11,3,6,16)]

moduleup$sign= "upregulated"
moduledown$sign= "downregulated"
resultenrichment= rbind(moduleup, moduledown)
write.table(resultenrichment, paste0( gseid,"_enrichment", ".xls"))

# if you want to check the expression value of selected genes in disease and control grupes:
# for two conditions:

genes = data.frame(genes=c("SYTL2", "SPINK6", "MMP11", "ESRRA","TPRG1", "NSMF")) #---> get the genes from excel table or just write their name
genes= unique(genes$genes)
df_unique= data.frame(genesym= rownames(df_unique), df_unique)
selected= df_unique[df_unique$genesym %in% genes, ]
metadata_selected= metadata[, c("Accession","type")] #------> Accession and type

data_long <- selected %>%
  pivot_longer(cols = -genesym, names_to = "Sample", values_to = "Expression")

data_long <- data_long %>%
  left_join(metadata_selected, by = c("Sample" = "Accession"))

colnames(data_long)[4]= "type2" #--------> one of the column represents the disease status

# if you want to rearrange the disease status names
data_long$type2[which(data_long$type2=="TMZ1")]="TMZ2"

ggplot(data_long, aes(x = type2, y = Expression, fill = type2)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  stat_compare_means(comparisons = list(c("control", "TMZ2")), #-------> change the status names
                     method = "t.test", 
                     label = "p.signif", 
                     tip.length = 0.05) +
  facet_wrap(~ genesym, scales = "free_y") +
  theme_minimal() +
  theme(strip.text = element_text(size = 15, color = "black")) + # the size and the color of the gene symbols
  labs(title = gseid,
       x = "Group",
       y = "Expression Level")
