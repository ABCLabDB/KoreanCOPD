####################################
## Korean COPDGene project        ##
## Replication code               ##
## 2023-10-04                     ##
## Minseok Seo                    ##
####################################

## Load packages
setRepositories(ind = 1:7)

library(data.table)
library(edgeR)
library(ggplot2)
library(ggsci)
library(limma)
library(ggpubr)
library(dplyr)
library(stringr)
library(corrplot)

## Set Working Dir.
WORK_DIR <- "C:/Project/2.KorCOPDGene/5.Github"
setwd(WORK_DIR)


## Load misc. file
masterFile <- data.frame(fread("./0.Data/masterFile.txt", head = T, sep = "\t"))
dim(masterFile)


## Load raw count data
count <- data.frame(fread("./0.Data/RawCountMatrix.txt", head = T, stringsAsFactors = F))
geneSymbol <- count[,1]
count <- count[,-1] # remove gene id from the matrix
dim(count)


## Removing non-expressed genes 
indexRemoval <- which(rowSums(count) == 0)

count <- count[-indexRemoval,]
geneSymbol <- geneSymbol[-indexRemoval]

dim(count)
length(geneSymbol)

rm(indexRemoval)


## Get and mapping gene annotation
geneAnno <- data.frame(fread("./0.Data/mart_export_GRCh38.p13.txt", head = T, stringsAsFactors = F))

temp <- left_join(data.frame(geneSymbol), geneAnno, by = "geneSymbol")
all.equal(geneSymbol, temp$geneSymbol)

geneAnno <- temp
rm(temp)



## Drawing Figure 2A
plotData <- data.frame(sampleID = masterFile$sampleID,
                       totalReads = masterFile$QC_Total_sequences,
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggbarplot(plotData, x="sampleID", y="totalReads", fill = "Phase", palette = "ucscgb") +
  geom_hline(yintercept = 20000000, col = "red", linetype = 3, size = 1.2) + xlab("Samples") + ylab("# of reads generated")

ggsave("Fig2A.pdf", height = 5, width = 14)


## Drawing Figure 2B
plotData <- data.frame(sampleID = masterFile$sampleID,
                       GCRatio = masterFile$QC_Ratio_of_GC,
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="GCRatio", fill = "Phase", palette = "ucscgb",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3),
          ylim = c(30, 70)) +
  geom_hline(yintercept = 60, col = "grey", linetype = 2, size = 1.2) +
  geom_hline(yintercept = 40, col = "grey", linetype = 2, size = 1.2) +
  theme(legend.position="none") +
  ylab("GC Ratio") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig2B.pdf", height = 6, width = 5)


## Drawing Figure 2C
plotData <- data.frame(sampleID = masterFile$sampleID,
                       PhreadScore = masterFile$QC_Mean_of_mean_per_base_sequence,
                       Phase = masterFile$Final.Phase)


ggboxplot(plotData, x="Phase", y="PhreadScore", fill="Phase", palette = "ucscgb",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3),
          ylim = c(25, 40)) +
  geom_hline(yintercept = 30, col = "red", linetype = 2, size = 1.2) +
  theme(legend.position="none") +
  ylab("Phred Score") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig2C.pdf", height = 6, width = 5)



## Drawing Figure 2D
plotData <- data.frame(sampleID = masterFile$sampleID,
                       RatioN = masterFile$QC_result_Mean_base_N,
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="RatioN", fill = "Phase", palette = "ucscgb",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3),
          ylim = c(0, 0.5)) +
  theme(legend.position="none") +
  ylab("Proportion of 'N' (%)") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig2D.pdf", height = 6, width = 5)


## Drawing Figure 2E
plotData <- data.frame(sampleID = masterFile$sampleID,
                       IlluminaAdapter = masterFile$Illumina.Universal.Adapter,
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="IlluminaAdapter", fill = "Phase", palette = "ucscgb",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3)) +
  theme(legend.position="none") +
  ylab("Illumina Adapter (%)") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig2E.pdf", height = 6, width = 5)


## Drawing Figure 3A
plotData <- data.frame(sampleID = masterFile$sampleID,
                       survivingReads = as.numeric(gsub("%", "", masterFile$surviving_percent)),
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="survivingReads", fill = "Phase", palette = "simpsons",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3),
          ylim = c(0.95, 1.0)) +
  theme(legend.position="none") +
  ylab("Clean reads (%)") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig3A.pdf", height = 4, width = 6)



## Drawing Figure 3B
plotData <- data.frame(sampleID = masterFile$sampleID,
                       Forward_only_surviving_percent = as.numeric(gsub("%", "", masterFile$Forward_only_surviving_percent)),
                       Phase = masterFile$Final.Phase)


plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="Forward_only_surviving_percent", fill = "Phase", palette = "simpsons",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3)) +
  theme(legend.position="none") +
  ylab("Only forward reads surviving (%)") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig3B.pdf", height = 4, width = 6)



## Drawing Figure 3C
plotData <- data.frame(sampleID = masterFile$sampleID,
                       Reverse_only_surviving_percent = as.numeric(gsub("%", "", masterFile$Reverse_only_surviving_percent)),
                       Phase = masterFile$Final.Phase)


plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="Reverse_only_surviving_percent", fill = "Phase", palette = "simpsons",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3)) +
  theme(legend.position="none") +
  ylab("Only reverse reads surviving (%)") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig3C.pdf", height = 4, width = 6)



## Drawing Figure 3D
plotData <- data.frame(sampleID = masterFile$sampleID,
                       Trimmed_input_read = masterFile$Trimmed_input_read,
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggbarplot(plotData, x="sampleID", y="Trimmed_input_read", fill = "Phase", palette = "simpsons") +
  geom_hline(yintercept = 20000000, col = "red", linetype = 3, size = 1.2) + 
  xlab("Samples") + ylab("# of trimmed reads")

ggsave("Fig3D.pdf", height = 5, width = 18)


## Drawing Figure 4A
plotData <- data.frame(sampleID = masterFile$sampleID,
                       Overall_Mapping_Rate = as.numeric(gsub("%", "", masterFile$Overall_Mapping_Rate)),
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="Overall_Mapping_Rate", fill = "Phase", palette = "rickandmorty",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3),
          ylim = c(0.9, 1.0)) +
  geom_hline(yintercept = 0.95, col = "red", linetype = 2, size = 1.2) +
  theme(legend.position="none") +
  ylab("Overall Mapping Rate (%)") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig4A.pdf", height = 4, width = 6)



## Drawing Figure 4B
plotData <- data.frame(sampleID = masterFile$sampleID,
                       Aligned_concordantly_more_1_times_Rate = as.numeric(gsub("%", "", masterFile$Aligned_concordantly_more_1_times_Rate)),
                       Phase = masterFile$Final.Phase)


plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="Aligned_concordantly_more_1_times_Rate", fill = "Phase", palette = "rickandmorty",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3)) +
  theme(legend.position="none") +
  ylab("Reads mapped to multiple region (%)") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig4B.pdf", height = 4, width = 6)



## Drawing Figure 4C
plotData <- data.frame(sampleID = masterFile$sampleID,
                       Assignment.Rate = as.numeric(gsub("%", "", masterFile$Assignment.Rate)),
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggboxplot(plotData, x="Phase", y="Assignment.Rate", fill = "Phase", palette = "rickandmorty",
          add="jitter", add.params = list(color = "Phase", size = 3, alpha = 0.3)) +
  geom_hline(yintercept = 0.3, col = "grey", linetype = 3, size = 1.2) +
  theme(legend.position="none") +
  ylab("Proportion of reads quantified (%)") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Fig4C.pdf", height = 4, width = 6)


## Drawing Figure 4D
plotData <- data.frame(sampleID = masterFile$sampleID,
                       Assigned = masterFile$Assigned,
                       Phase = masterFile$Final.Phase)

plotData <- plotData[order(plotData$Phase),]

ggbarplot(plotData, x="sampleID", y="Assigned", fill = "Phase", palette = "rickandmorty") +
  geom_hline(yintercept = 20000000, col = "red", linetype = 3, size = 1.2) +
  xlab("Samples") + ylab("# of quantified reads")

ggsave("Fig4D.pdf", height = 5, width = 18)



## Defining experimental variables for statistical tests
Phase <- factor(masterFile$Final.Phase)
SeqPlatform <- masterFile$Final.Phase

SeqPlatform[which(SeqPlatform == "Phase1")] <- "Nova-seq"
SeqPlatform[which(SeqPlatform == "Phase2")] <- "Nova-seq"
SeqPlatform[which(SeqPlatform == "Phase3")] <- "Nova-seq"
SeqPlatform[which(SeqPlatform == "Phase4")] <- "Nova-seq"

SeqPlatform[which(SeqPlatform == "Phase5")] <- "MGI-seq"
SeqPlatform[which(SeqPlatform == "Phase6")] <- "MGI-seq"
SeqPlatform[which(SeqPlatform == "Phase7")] <- "MGI-seq"

SeqPlatform <- factor(SeqPlatform)


## DEG identification (NovaSeq vs MGISeq)
design <- model.matrix(~SeqPlatform + COPD_Grade + current_smoke + gender + age + Ht + Wt, data = masterFile)
data.frame(colnames(design))
colnames(design)

Y <- DGEList(counts=count, gene=geneSymbol)
Y <- calcNormFactors(Y, method="TMM")

logCPM <- cpm(Y, normalized.lib.sizes=TRUE, log=T)


## Voom transformation and eBayes
V <- voom(Y, design, plot = T)

fit <- lmFit(V, design)
fit <- eBayes(fit)


## Get result H0: two sequencing platforms are equal
fit_Platdiff <- contrasts.fit(fit, coef = 2)
fit_Platdiff <- eBayes(fit_Platdiff)

Result_Platdiff <- topTable(fit_Platdiff, number = nrow(logCPM), sort.by = "none")

sum(Result_Platdiff$adj.P.Val <= 0.05)


## Exporting Stat results
output <- data.frame(geneSymbol,
                     geneAnno,
                     logFC = Result_Platdiff$logFC,
                     P = Result_Platdiff$P,
                     adjP = Result_Platdiff$adj.P.Val)

write.table(output, "Result_DEGs_SeqPlatform.txt", sep = "\t", quote = F, row.names = F)
rm(output)


## Drawing MDS plot
MDS_data <- plotMDS(Y)

plot_data <- data.frame(Phase,
                        Gender = Sex,
                        currentSmoke = CurrentSmoke,
                        COPDGrade,
                        X=MDS_data$x, Y=MDS_data$y)


p1 <- ggplot(plot_data, aes(x=X, y=Y, color=Phase)) +
  geom_point(size = 5, alpha = 0.8) +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  theme_classic()+
  scale_color_nejm()

p1

ggsave("Fig5A_withLegend.pdf", width = 6, height = 6)



## DEG identification (Associated with the currentSmoke)
design <- model.matrix(~current_smoke + COPD_Grade + Phase + gender + age + Ht + Wt, data = masterFile)
data.frame(colnames(design))
colnames(design)

Y <- DGEList(counts=count, gene=geneSymbol)
Y <- calcNormFactors(Y, method="TMM")

logCPM <- cpm(Y, normalized.lib.sizes=TRUE, log=T)

## Limma modeling
V <- voom(Y, design, plot = T)

fit <- lmFit(V, design)
fit <- eBayes(fit)


## Get result for currentSmoke
fit_currentSmoke <- contrasts.fit(fit, coef = 2)
fit_currentSmoke <- eBayes(fit_currentSmoke)

Result_currentSmoke <- topTable(fit_currentSmoke, number = nrow(logCPM), sort.by = "none")

sum(Result_currentSmoke$adj.P.Val <= 0.05)


## Exporting results
output <- data.frame(geneSymbol,
                     geneAnno,
                     logFC_currentSmoke = Result_currentSmoke$logFC,
                     P_currentSmoke = Result_currentSmoke$P,
                     adjP_currentSmoke = Result_currentSmoke$adj.P.Val)

write.table(output, "Result_DEGs_CurrentSmoke.txt", sep = "\t", quote = F, row.names = F)


## Drawing boxplots for top 8 DEGs

targetIdx <- order(Result_currentSmoke$P.Value)[1:8]


if(geneAnno$GeneName[targetIdx[1]] != ""){
  geneName <- geneAnno$GeneName[targetIdx[1]]
}else{
  geneName <- geneAnno$geneSymbol[targetIdx[1]]
}

plotData <- data.frame(geneExp = as.numeric(logCPM[targetIdx[1],]),
                       CurrentSmoke,
                       geneSymbol = rep(geneName, ncol(logCPM)))


for(i in 2:8){
  if(geneAnno$GeneName[targetIdx[i]] != ""){
    geneName <- geneAnno$GeneName[targetIdx[i]]
  }else{
    geneName <- geneAnno$geneSymbol[targetIdx[i]]
  }
  
  temp <- data.frame(geneExp = as.numeric(logCPM[targetIdx[i],]),
                     CurrentSmoke,
                     geneSymbol = rep(geneName, ncol(logCPM)))
  
  plotData <- rbind(plotData, temp)
}

plotData$geneSymbol <- factor(plotData$geneSymbol,
                              levels = c("GPR15", "LRRN3", "PID1", "SASH1", "ENSG00000227240", "ENSG00000286285", "SEMA6B", "CLDND1"))


colnames(plotData)

ggboxplot(plotData, x = "CurrentSmoke", y = "geneExp", fill = "CurrentSmoke", palette = "simpsons", add = "jitter",
          ylab="TMM noramlized logCPM", xlab = "Current smoking Status", facet.by = "geneSymbol", scales = "free_y", nrow=2,
          add.params = list(alpha=0.1, size=4, color = "CurrentSmoke")) +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Never/Ever Smoker",
                            "Current Smoker")) +
  geom_violin(alpha = 0.2)

ggsave("Fig5C.pdf", width = 14, height = 4)

