###########################################################################################
###########################################################################################
###########################################################################################
################################## Locatelli   ############################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

#Commands copied from prep.R and merge here to better understand the process and steps
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli")
#cp BetaValues_CTRL_MEDIAN_UCSC.txt BetaValues_CTRL_MEDIAN.gtf  #These are Avg of Control of 20 Individulas
#cp BetaValues_Locatelli_Lorenzo_UCSC.txt BetaValues_Locatelli_Lorenzo.gtf
#cp BetaValues_Locatelli_Mario_UCSC.txt BetaValues_Locatelli_Mario.gtf
#cp BetaValues_Locatelli_Madre_UCSC.txt BetaValues_Locatelli_Madre.gtf
#cp BetaValues_Locatelli_Padre_UCSC.txt BetaValues_Locatelli_Padre.gtf
#sed -i '1d'  *.gtf

control=read.table("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/BetaValues_CTRL_MEDIAN.gtf", header = FALSE, stringsAsFactors = FALSE)
lorenzo=read.table("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/BetaValues_Locatelli_Lorenzo.gtf", header = FALSE, stringsAsFactors = FALSE)
mario=read.table("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/BetaValues_Locatelli_Mario.gtf", header = FALSE, stringsAsFactors = FALSE)
madre=read.table("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/BetaValues_Locatelli_Madre.gtf", header = FALSE, stringsAsFactors = FALSE)
padre=read.table("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/BetaValues_Locatelli_Padre.gtf", header = FALSE, stringsAsFactors = FALSE)
head(control)
clmmp = cbind(control, lorenzo, mario, madre, padre)
head(clmmp)
clmmp = data.frame(clmmp)
head(clmmp)
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli")
#write.table(clmmp, "CLMMP.txt" , append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#Prep file for PCA: awk '{print $1"%"$2"%"$3"\t"$4"\t"$8"\t"$12"\t"$16"\t"$20}' CLMMP.txt > CLMMP_pca.txt
#More easy way: 
CLMMP <- clmmp[, c(1:4, 8, 12, 16, 20)]
head(CLMMP)
CLMMP_PCA <- cbind(paste(CLMMP[,1],CLMMP[,2],CLMMP[,3], sep = "%"), CLMMP[,c(4:8)])
colnames(CLMMP_PCA) <- c("ID","Control", "Lorenzo", "Mario", "Madre", "Padre")
head(CLMMP_PCA)
#write.table(CLMMP_PCA, "CLMMP_pca.txt", quote = FALSE, append = FALSE ,row.names = T)
#Plot PCA by group color and labelling
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli")
mydata1=read.table("CLMMP_pca.txt", header = TRUE, stringsAsFactors = FALSE)
rownames(mydata1)
mydata1[,1]
rownames(mydata1)=mydata1[,1]
rownames(mydata1)
colnames(mydata1)
mydata1 = mydata1[,-1]
mydata1 <- data.frame(mydata1)
head(mydata1)
str(mydata1)
head(mydata1)
df <- as.data.frame(mydata1)
head(df)
df_pca <- prcomp(df)
df_out <- as.data.frame(df_pca$rotation)
head(df_out)
df_out["Color"] <- c("Control", "Lorenzo","Mario","Madre","Padre")
percentage <- round(df_pca$sdev^2 / sum(df_pca$sdev^2) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=df_out$Color, label=df_out$Color))
p<-p+geom_point()+theme + xlab(percentage[1]) + ylab(percentage[2])
p <- p+theme_bw() 
p
ggsave("CLMMP_PCA.svg", width=10*1.25, height=10*1.25, units="cm", dpi=96)

pca <- prcomp(mydata1, scale=T)
plot(df_pca$rotation[,1],df_pca$rotation[,2])
text(df_pca$rotation[,1],df_pca$rotation[,2], colnames(mydata1) , pos = 2, col = "black")
#ankitv@ankitv-lab-riccio:/media/ankitv/Archivio2/ankit/Locatelli$ sort -k1,1 -k2,2n CLMMP_readysep.txt | awk '{print $1"\t"$2"\t"$3"\t"$4}' > Control_sort.txt
#ankitv@ankitv-lab-riccio:/media/ankitv/Archivio2/ankit/Locatelli$ sort -k1,1 -k2,2n CLMMP_readysep.txt | awk '{print $1"\t"$2"\t"$3"\t"$5}' > Lorenzo_sort.txt
#ankitv@ankitv-lab-riccio:/media/ankitv/Archivio2/ankit/Locatelli$ sort -k1,1 -k2,2n CLMMP_readysep.txt | awk '{print $1"\t"$2"\t"$3"\t"$6}' > Mario_sort.txt
#ankitv@ankitv-lab-riccio:/media/ankitv/Archivio2/ankit/Locatelli$ sort -k1,1 -k2,2n CLMMP_readysep.txt | awk '{print $1"\t"$2"\t"$3"\t"$7}' > Madre_sort.txt
#ankitv@ankitv-lab-riccio:/media/ankitv/Archivio2/ankit/Locatelli$ sort -k1,1 -k2,2n CLMMP_readysep.txt | awk '{print $1"\t"$2"\t"$3"\t"$8}' > Padre_sort.txt
#bedtools intersect -wa -wb -a CLMMP_readysep.txt -b ~/ref_av/human_ICR.bed > CLMMP_readysep_humanICR.txt
#Take average using excel sheet: see excel sheet and save as CLMMP_heatmap.txt
#Heatmap using control of 20 Individuals
library(gplots)
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/")
data <- read.table("CLMMP_heatmap.txt",header=TRUE)
rownames(data)
rownames(data)=data[,1]
rownames(data)
colnames(data)
data = data[,-1]
data = as.matrix(data)
head(data)
dim(data)
colfunc <- colorRampPalette(c("#291BEB","white", "red"))
heatmap.2(data, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))

setwd("/media/ankitv/Archivio2/ankit/Locatelli/scatter/2kb+1")
mydata=read.table("CM_merge.scatter.txt", header = TRUE, stringsAsFactors = FALSE)
rownames(mydata)
mydata[,1]
rownames(mydata)=mydata[,1]
rownames(mydata)
colnames(mydata)
mydata = mydata[,-1]
mydata <- data.frame(mydata)
head(mydata)
x = mydata$Control
y = mydata$Mario
df <- data.frame(x = mydata$Control, y = mydata$Mario,
                 d = densCols(x, y, colramp = colorRampPalette(c("blue", "green", "red"))))
p <- ggplot(df) + geom_hex(aes(x, y), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
print(p)
ggsave("CM_merge.scatter.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)
ggsave("CM_merge.scatter.png", width=15*1.25, height=15*1.25, units="cm", dpi=96)


setwd("/media/ankitv/Archivio2/ankit/Locatelli/scatter/2kb+1")
mydata=read.table("LM_merge.scatter.txt", header = TRUE, stringsAsFactors = FALSE)
rownames(mydata)
mydata[,1]
rownames(mydata)=mydata[,1]
rownames(mydata)
colnames(mydata)
mydata = mydata[,-1]
mydata <- data.frame(mydata)
head(mydata)
x = mydata$Lorenzo
y = mydata$Mario
df <- data.frame(x = mydata$Lorenzo, y = mydata$Mario,
                 d = densCols(x, y, colramp = colorRampPalette(c("blue", "green", "red"))))
p <- ggplot(df) + geom_hex(aes(x, y), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
print(p)
ggsave("LM_merge.scatter.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)
ggsave("LM_merge.scatter.png", width=15*1.25, height=15*1.25, units="cm", dpi=96)

##################################     Data analysed with control from Dec 2018 6 Individuals  ################################
LMMP <- CLMMP_PCA[,c(1,3:6)]
head(LMMP)
#write.table(LMMP, "LMMP.txt" , append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

Control6_dec2018 =read.table("CONTROLS_NEW_TAB_dec2018.txt", header = TRUE, stringsAsFactors = FALSE)
head(Control6_dec2018)
LMMPC6 <- cbind(LMMP, Control6_dec2018)
head(LMMPC6)
LMMPC6b <- LMMPC6[,c(1:5, 9:14)]
head(LMMPC6b)
colnames(LMMPC6b) <- c("Chr","Lorenzo", "Mario", "Madre", "Padre", "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5", "CTRL6")
#write.table(LMMPC6b, "PCA_Locatelli_family_Controls-dec2018.txt" , append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
LMMPC6a <- LMMPC6[,c(6, 9:14, 2:5)]
head(LMMPC6a)
#write.table(LMMPC6a, "PCA_Locatelli_family_Controls-dec2018_probechr.txt" , append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
roundLMMPC6a <- data.frame(cbind(LMMPC6a[,1], round(LMMPC6a[,c(2:11)], 3)))
roundLMMPC6a <- roundLMMPC6a[,c(1:8, 10,9,11)]
colnames(roundLMMPC6a) <- c("ID_REF","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","L1","L2","L3","L4")
head(roundLMMPC6a)
#Extract rounded off file for uploading to GEO
#write.table(roundLMMPC6a, "/media/ankitv/Archivio2/ankit/BWS_array/GEOSubmission/PCA_Locatelli_family_Controls-dec2018_probechr.txt" , append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
##awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($6+$7+$8+$9+$10+$11)/6}' PCA_Locatelli_family_Controls-dec2018.txt  > PCA_Locatelli_family_Controls-dec2018.avg%.txt
#Column name as Chr Lorenzo Mario Madre Padre Control
#awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($6+$7+$8+$9+$10+$11)/6}' PCA_Locatelli_family_Controls-dec2018.txt | awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > PCA_Locatelli_family_Controls-dec2018.avg.txt
#bedtools intersect -wa -wb -a PCA_Locatelli_family_Controls-dec2018.avg.txt -b ~/ref_av/human_ICR.bed > PCA_Locatelli_family_Controls-dec2018.avg_human_ICR.txt
#awk '{print $12"\t"$13"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' PCA_Locatelli_family_Controls-dec2018.avg_human_ICR.txt > PCA_Locatelli_family_Controls-dec2018.avg_human_ICR.rearranged.txt
probename_chrpos_EPIC <- read.table("/media/ankitv/Archivio2/ankit/BWS_array/probename_chrpos_EPIC.txt",header=TRUE)
head(probename_chrpos_EPIC)

library(ggpubr)
library(ggplot2)
library(plyr)
#Add column names as follows: DMR DMRType Lorenzo Mario	Madre	Padre	Control
#PCA with 6 Controls Avg Human ICR
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli")
data1 <- read.table("PCA_Locatelli_family_Controls-dec2018.avg_human_ICR.rearranged.txt", header = TRUE)
rownames(data1)
colnames(data1)
head(data1)
Count_dataICR1 <- count(data1, "DMR")
head(Count_dataICR1)
aggregate1 = aggregate(data1[,3:7],by=list(data1$DMR), median)
head(aggregate1, 2)
#Plot PCA by group color and labelling
dataICR1=aggregate1
rownames(dataICR1)
dataICR1[,1]
rownames(dataICR1)=dataICR1[,1]
rownames(dataICR1)
colnames(dataICR1)
dataICR1 = dataICR1[,-1]
head(dataICR1)
dim(dataICR1)
dataICR1_counted <- cbind(dataICR1, Count_dataICR1)
head(dataICR1_counted)
dim(dataICR1_counted)
dataICR1_counted1 <- dataICR1_counted[which(dataICR1_counted$freq >= 3),]
head(dataICR1_counted1)
dim(dataICR1_counted1)
##write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
dataICR1_counted2 <- dataICR1_counted1[,1:5]
head(dataICR1_counted2)
summary(dataICR1_counted2)
#df <- as.data.frame(dataICR1)
dfICR <- dataICR1_counted2
head(dfICR)
dfICR_pca <- prcomp(dfICR)
dfICR_out <- as.data.frame(dfICR_pca$rotation)
head(dfICR_out)
dfICR_out["Color"] <- c("Lorenzo" ,  "Mario"  , "Madre" ,  "Padre", "Control")
percentageICR <- round(dfICR_pca$sdev^2 / sum(dfICR_pca$sdev^2) * 100, 2)
percentageICR <- paste( colnames(dfICR_out), "(", paste( as.character(percentageICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
head(dfICR_out)
#Color By hroup as mentioned in file color
p<-ggplot(dfICR_out,aes(x=PC1,y=PC2,color=dfICR_out$Color, label=dfICR_out$Color))
p<-p+geom_point()+theme + xlab(percentageICR[1]) + ylab(percentageICR[2])
p <- p+theme_bw()
p
ggsave("/media/ankitv/Archivio2/ankit/BWS_array/Median/PCA_Loca.humanICR.averaged.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

################################### Heatmap human ICR ######################################
#Heatmap with 6 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
dataICR2 = as.matrix(dataICR1_counted2)
head(dataICR2)
dim(dataICR2)
#colnames(dataICR2) <- c("Lorenzo",	"Mario",	"Madre",	"Padre",	"Control")
colfunc <- colorRampPalette(c("navy","white", "darkred"))
dataICR3 <- dataICR2[,c(5,4,3,1,2)]
head(dataICR3)
write.table(dataICR3, "/media/ankitv/Archivio2/ankit/BWS_array/Median/Heatmap6_Locatelli_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
#Control Samples
#heatmap.2(dataICR2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#Let it cluster
svg(filename="/media/ankitv/Archivio2/ankit/BWS_array/Median/Heatmap_Locatelli_family_Controls-dec2018.avg_human_ICR.rearranged_SVG.svg", width=5, height=10, pointsize=12)
#heatmap.2(dataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(dataICR3, col = colfunc, Colv = "NA",dendrogram = c("none"), trace = "none", keysize=0.8, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

svg(filename="/media/ankitv/Archivio2/ankit/BWS_array/Median/Heatmap_Locatelli_family_Controls-dec2018.avg_human_ICR.rearranged_SVG_clusterboth.svg", width=5, height=10, pointsize=12)
#heatmap.2(dataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(dataICR3, col = colfunc, dendrogram = c("both"), trace = "none", key=TRUE, density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

#Box and Violin Plot human ICR using 6 Controls from Dec 2018
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli")
data1 <- read.table("PCA_Locatelli_family_Controls-dec2018.avg_human_ICR.rearranged.txt", header = TRUE)
rownames(data1)
colnames(data1)
head(data1)
Count_dataICR1 <- count(data1, "DMR")
head(Count_dataICR1)
aggregate1 = aggregate(data1[,3:7],by=list(data1$DMR), median)
head(aggregate1, 2)
dataICR1=aggregate1
rownames(dataICR1)
dataICR1[,1]
rownames(dataICR1)=dataICR1[,1]
rownames(dataICR1)
colnames(dataICR1)
dataICR1 = dataICR1[,-1]
head(dataICR1)
dim(dataICR1)
#write.table(dataICR1, "aggregated_dataICR1.txt", sep="\t", quote = FALSE, append = FALSE)
dataICR1_counted <- cbind(dataICR1, Count_dataICR1)
head(dataICR1_counted)
dim(dataICR1_counted)
dataICR1_counted1 <- dataICR1_counted[which(dataICR1_counted$freq >= 3),]
head(dataICR1_counted1)
dim(dataICR1_counted1)
##write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
dataICR1_counted2 <- dataICR1_counted1[,1:5]
head(dataICR1_counted2)
VioICR <- data.frame(dataICR1_counted2[,c(5,4,3,1,2)])
head(VioICR)
write.table(VioICR, "/media/ankitv/Archivio2/ankit/BWS_array/Median/VioICR.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
boxplot(VioICR, main="Average methylation at human ICRs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkred", "darkorange", "darkgreen", "blue"))
VioICR2 <- stack(VioICR)
head(VioICR2)
colnames(VioICR2) <- c("Methylation", "Datasets")
head(VioICR2)
ggplot(VioICR2, aes(x=Datasets, y=Methylation, color=Datasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("black","blue", "magenta", "red", "orange"))
ggsave("/media/ankitv/Archivio2/ankit/BWS_array/Median/Violin_plot_Locatelli_humanICR_6Control.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#Box and Violin Plot Global methylation using 6 Controls avg from Dec 2018
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/")
dataGlob=read.table("PCA_Locatelli_family_Controls-dec2018.avg%.txt", header = TRUE, stringsAsFactors = FALSE)
rownames(dataGlob)
dataGlob[,1]
rownames(dataGlob)=dataGlob[,1]
rownames(dataGlob)
colnames(dataGlob)
dataGlob1 = dataGlob[,-1]
dataGlob1 <- data.frame(dataGlob1)
head(dataGlob1)
str(dataGlob1)
VioGlob <- dataGlob1[,c(5,4,3,1,2)]
head(VioGlob)
boxplot(VioGlob, main="Global average methylation", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkred", "darkorange", "darkgreen", "blue"))
VioGlob2 <- stack(VioGlob)
colnames(VioGlob2) <- c("Methylation", "Datasets")
head(VioGlob2)
ggplot(VioGlob2, aes(x=Datasets, y=Methylation, color=Datasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("black","blue", "magenta", "red", "orange"))
ggsave("/media/ankitv/Archivio2/ankit/BWS_array/Median/Violin_plot_Locatelli_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#Scatter
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/scatter/2kb+1")
mydata=read.table("CL_merge.scatter.txt", header = TRUE, stringsAsFactors = FALSE)
rownames(mydata)
mydata[,1]
rownames(mydata)=mydata[,1]
rownames(mydata)
colnames(mydata)
mydata = mydata[,-1]
mydata <- data.frame(mydata)
head(mydata)
x = mydata$Control
y = mydata$Lorenzo
df <- data.frame(x = mydata$Control, y = mydata$Lorenzo,
                 d = densCols(x, y, colramp = colorRampPalette(c("blue", "green", "red"))))
p <- ggplot(df) + geom_hex(aes(x, y), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
print(p)
ggsave("CL_merge.scatter.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)
ggsave("CL_merge.scatter.png", width=15*1.25, height=15*1.25, units="cm", dpi=96)


#Scatter
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/scatter/2kb+1/rep_CpG_removed")
Loca_repCpGminus.2kb_1 =read.table("PCA_Locatelli_family_Controls-dec2018.sep.repCpGminus.2kbp.txt", header = FALSE, stringsAsFactors = FALSE)
head(Loca_repCpGminus.2kb_1)
aggregate_Loca_repCpGminus.2kb_1 = aggregate(Loca_repCpGminus.2kb_1[,5:14],by=list(Loca_repCpGminus.2kb_1$V1), median)
head(aggregate_Loca_repCpGminus.2kb_1)
aggregate_Loca_repCpGminus.2kb_1_sorted = aggregate_Loca_repCpGminus.2kb_1[order(aggregate_Loca_repCpGminus.2kb_1$Group.1),]
head(aggregate_Loca_repCpGminus.2kb_1_sorted)
rownames(aggregate_Loca_repCpGminus.2kb_1_sorted)
rownames(aggregate_Loca_repCpGminus.2kb_1_sorted)=aggregate_Loca_repCpGminus.2kb_1_sorted[,1]
rownames(aggregate_Loca_repCpGminus.2kb_1_sorted)
colnames(aggregate_Loca_repCpGminus.2kb_1_sorted)
aggregate_Loca_repCpGminus.2kb_1_sorted = aggregate_Loca_repCpGminus.2kb_1_sorted[,-1]
aggregate_Loca_repCpGminus.2kb_1_sorted = as.matrix(aggregate_Loca_repCpGminus.2kb_1_sorted)
head(aggregate_Loca_repCpGminus.2kb_1_sorted)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted)
Count_Loca_repCpGminus.2kb_1 <- count(Loca_repCpGminus.2kb_1, "V1")
aggregate_Loca_repCpGminus.2kb_1_sorted_counted <- cbind(aggregate_Loca_repCpGminus.2kb_1_sorted, Count_Loca_repCpGminus.2kb_1)
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted)

#CpG Nocutoff
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_Nocutoff <- aggregate_Loca_repCpGminus.2kb_1_sorted_counted
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_Nocutoff)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_Nocutoff)
#write.table(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_Nocutoff, "aggregate_Loca_repCpGminus.2kb_1_sorted_counted_Nocutoff.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff <- aggregate_Loca_repCpGminus.2kb_1_sorted_counted_Nocutoff[,1:10]
colnames(aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff) <- c("Lorenzo","Mario","Madre","Padre","CTRL1","CTRL2","CTRL3","CTRL4","CTRL5","CTRL6")
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff)
colfunc <- colorRampPalette(c("#291BEB", "red"))
ggplot(aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff) + geom_hex(aes(x=Lorenzo, y=Mario), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
ggsave("scatter_Lorenzo_Mario.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)
#write.table(aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff, "aggregate_Loca_repCpGminus.2kb_1_sorted_counted2_Nocutoff.txt", quote = FALSE, append = FALSE, sep = "\t")


#Lorenzo_Mario
ggplot(aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff) + geom_hex(aes(x=Lorenzo, y=Mario), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
ggsave("scatter_Lorenzo_Mario.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)

#--------------------------Average control vs patients -------------------------------------

head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff)
myP <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted1_Nocutoff)
head(myP)
avg_CL <- cbind(myP[1],rowMeans(myP[5:10]))
head(avg_CL)
colnames(avg_CL) <- c("Lorenzo","CTRLN")
head(avg_CL)
ggplot(avg_CL) + geom_hex(aes(x=CTRLN, y=Lorenzo), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
ggsave("Nscatter_Lorenzo_CTRLN.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)

avg_CM <- cbind(myP[2],rowMeans(myP[5:10]))
head(avg_CM)
colnames(avg_CM) <- c("Mario","CTRLN")
head(avg_CM)
ggplot(avg_CM) + geom_hex(aes(x=CTRLN, y=Mario), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
ggsave("Nscatter_Mario_CTRLN.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)


###################### NEW DMR/ Hypomethylation regions in Locatelli data: CH t-test #####################à
#Hypomethylation CH-t-test:
#Remove repeated CpGs: fgrep -f CpGmayrepeated_cordsmerge.txt PCA_Locatelli_family_Controls.txt -v > PCA_Locatelli_family_Controls.repCpGminus.txt
#awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27}' PCA_Locatelli_family_Controls.repCpGminus.txt | grep chr > PCA_Locatelli_family_Controls.repCpGminus.sep.txt
#bedtools intersect -wa -wb -a hg19.2kb+1_corrected.txt -b PCA_Locatelli_family_Controls.repCpGminus.sep.txt | awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30}' > PCA_Locatelli_family_Controls.repCpGminus.sep.2kbp.txt

setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/scatter/2kb+1/rep_CpG_removed")
Loca_repCpGminus.2kb_1 =read.table("PCA_Locatelli_family_Controls.repCpGminus.sep.2kbp.txt", header = FALSE, stringsAsFactors = FALSE)
head(Loca_repCpGminus.2kb_1)

#Calculate median methylation in each 2kb bins (median of beta-values)
aggregate_Loca_repCpGminus.2kb_1 = aggregate(Loca_repCpGminus.2kb_1[,5:28],by=list(Loca_repCpGminus.2kb_1$V1), median)
head(aggregate_Loca_repCpGminus.2kb_1)
dim(aggregate_Loca_repCpGminus.2kb_1)
aggregate_Loca_repCpGminus.2kb_1_sorted = aggregate_Loca_repCpGminus.2kb_1[order(aggregate_Loca_repCpGminus.2kb_1$Group.1),]
head(aggregate_Loca_repCpGminus.2kb_1_sorted)
rownames(aggregate_Loca_repCpGminus.2kb_1_sorted)
rownames(aggregate_Loca_repCpGminus.2kb_1_sorted)=aggregate_Loca_repCpGminus.2kb_1_sorted[,1]
aggregate_Loca_repCpGminus.2kb_1_sorted = aggregate_Loca_repCpGminus.2kb_1_sorted[,-1]
aggregate_Loca_repCpGminus.2kb_1_sorted = as.matrix(aggregate_Loca_repCpGminus.2kb_1_sorted)
head(aggregate_Loca_repCpGminus.2kb_1_sorted)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted)

#Count number of CpGs in each regions and combine information with mean methylation
Count_Loca_repCpGminus.2kb_1 <- count(Loca_repCpGminus.2kb_1, "V1") 
aggregate_Loca_repCpGminus.2kb_1_sorted_counted <- cbind(aggregate_Loca_repCpGminus.2kb_1_sorted, Count_Loca_repCpGminus.2kb_1)
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted)
tail(aggregate_Loca_repCpGminus.2kb_1_sorted_counted)

#Extract 2kb regions with minimum 3CpGs
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff <- aggregate_Loca_repCpGminus.2kb_1_sorted_counted[which(aggregate_Loca_repCpGminus.2kb_1_sorted_counted$freq >= 3),]
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff)
#write.table(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff, "aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff.txt", quote = FALSE, append = FALSE, sep = "\t", row.names = T)

#Reassign column names
colnames(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff) <- c("Lorenzo","Mario","Madre","Padre","CTRL1","CTRL2","CTRL3","CTRL4","CTRL5","CTRL6","CTRL7","CTRL8","CTRL10","CTRL11","CTRL12","CTRL13","CTRL14","CTRL15","CTRL16","CTRL17","CTRL18","CTRL19","CTRL20","CTRL21","Regions","freq")
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff)
#write.table(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff, "aggregate_Loca20_repCpGminus.2kb_1_sorted_counted_3cutoff.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff)
# Mean = X, 
#In apply function, 1 is row and 2 is column, This command is validated by simulation data, matrix <- matrix(1:12, nrow = 3)
#Rowwise effect.apply(matrix,1,mean), Columnwise.effect.apply(matrix,2,mean)
Control20mean.2kb <- data.frame(apply(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff[,5:24],1,mean)) #This command is validated by excel sheet manually
colnames(Control20mean.2kb) <- "Control20mean.2kb"
head(Control20mean.2kb)
#Rowwise standard deviation, apply(matrix,1,sd)
Control20sd.2kb <- data.frame(apply(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff[,5:24],1,sd)) # Standard Deviation = s, #This command is validated by excel sheet manually see copy. 1.7.19
colnames(Control20sd.2kb) <- "Control20sd.2kb"
head(Control20sd.2kb)
#Alignment checked
#Alignment_check <- cbind(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff[,1:4],
#                                                                        rownames(Control20mean.2kb),
#                                                                        Control20mean.2kb,
#                                                                        rownames(Control20sd.2kb),
#                                                                        Control20sd.2kb,
#                                                                        aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff[,25:26])
##write.table(Alignment_check, "Alignment_check.txt", quote = FALSE, append = FALSE, sep = "\t", row.names = T)
#awk '{print $1"\t"$6"\t"$8"\t"$10}' Alignment_check.txt | head/tail

aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd <- cbind(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff[,1:4],
                                                                        Control20mean.2kb,
                                                                        Control20sd.2kb,
                                                                        aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff[,25:26])

head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
n <- 20
v <- sqrt((n + 1)/n)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["DenoBase"] <- Control20sd.2kb * v #Calculate Denominator
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
#Apply CH t-test
# tch = (x* -X) / s * sqrt((n+1)/n) 


#Lorenzo
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["xCL"] <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$Lorenzo - aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$Control20mean.2kb)  #x* - X, where x* =Mean methylation of each 2kb regions in the case sample
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["txCL"] <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$xCL/aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$DenoBase)				#Calculate t-values, txCL/DenoBase
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
#pt(tval), df=degfree)) #one-tailed  #Dr Claudia Suggested and explained to me
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["ptxCL"] <- data.frame(pt((aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$txCL), df=n-1))	#Calculate p-values
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["padjtxCL"] <- data.frame(p.adjust(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$ptxCL, method = "BH"))		#Calculate p-adj values, BH multiplicity correction
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)

#Mario
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["xCM"] <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$Mario - aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$Control20mean.2kb)  #x* - X, where x* =Mean methylation of each 2kb regions in the case sample
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["txCM"] <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$xCM/aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$DenoBase)					#Calculate t-values, txCM/DenoBase
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["ptxCM"] <- data.frame(pt((aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$txCM),df=n-1))		#Calculate p-values
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["padjtxCM"] <- data.frame(p.adjust(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$ptxCM, method = "BH"))		#Calculate p-adj values, BH multiplicity correction
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)

#Madre
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["xCMd"] <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$Madre - aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$Control20mean.2kb)  #x* - X, where x* =Mean methylation of each 2kb regions in the case sample
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["txCMd"] <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$xCMd/aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$DenoBase)				#Calculate t-values, txCMd/DenoBase
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["ptxCMd"] <- data.frame(pt((aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$txCMd),df=n-1))		#Calculate p-values
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["padjtxCMd"] <- data.frame(p.adjust(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$ptxCMd, method = "BH"))		#Calculate p-adj values, BH multiplicity correction
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)

#Padre
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["xCPd"] <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$Padre - aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$Control20mean.2kb)  #x* - X, where x* =Mean methylation of each 2kb regions in the case sample
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["txCPd"] <- data.frame(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$xCPd/aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$DenoBase)				#Calculate t-values, txCPd/DenoBase
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["ptxCPd"] <- data.frame(pt((aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$txCPd),df=n-1))		#Calculate p-values
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd["padjtxCPd"] <- data.frame(p.adjust(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$ptxCPd, method = "BH"))		#Calculate p-adj values, BH multiplicity correction
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd)

#Export Lorenzo p<0.05
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pLorenzoless0.05 <- aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd[which(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$padjtxCL < 0.05),]
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pLorenzoless0.05)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pLorenzoless0.05)
write.table(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pLorenzoless0.05, "/media/ankitv/Archivio2/ankit/BWS_array/Median/aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pLorenzoless0.05.txt", quote = FALSE, append = FALSE, sep = "\t", row.names =T)

#Export Mario p<0.05
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMarioless0.05 <- aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd[which(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$padjtxCM < 0.05),]
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMarioless0.05)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMarioless0.05)
write.table(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMarioless0.05, "/media/ankitv/Archivio2/ankit/BWS_array/Median/aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMarioless0.05.txt", quote = FALSE, append = FALSE, sep = "\t", row.names =T)

#Export Madre p<0.05
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMadreless0.05 <- aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd[which(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$padjtxCMd < 0.05),]
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMadreless0.05)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMadreless0.05)
write.table(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMadreless0.05, "/media/ankitv/Archivio2/ankit/BWS_array/Median/aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pMadreless0.05.txt", quote = FALSE, append = FALSE, sep = "\t", row.names =T)

#Export Padre p<0.05
aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pPadreless0.05 <- aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd[which(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd$padjtxCPd < 0.05),]
head(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pPadreless0.05)
dim(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pPadreless0.05)
write.table(aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pPadreless0.05, "/media/ankitv/Archivio2/ankit/BWS_array/Median/aggregate_Loca_repCpGminus.2kb_1_sorted_counted_3cutoff_meansd_pPadreless0.05.txt", quote = FALSE, append = FALSE, sep = "\t", row.names =T)



#Statistical Analysis for Locatelli family
dim(VioICR)
summary(VioICR)
#Data is paired
wilcox.test(VioICR$Control, VioICR$Lorenzo, paired = TRUE, alternative = "two.sided")
#wilcox.test(VioICR$Control, VioICR$Lorenzo, paired = FALSE)
wilcox.test(VioICR$Control, VioICR$Mario, paired = TRUE, alternative = "two.sided")
#wilcox.test(VioICR$Control, VioICR$Mario, paired = FALSE)
wilcox.test(VioICR$Control, VioICR$Madre, paired = TRUE, alternative = "two.sided")
wilcox.test(VioICR$Control, VioICR$Padre, paired = TRUE, alternative = "two.sided")


dim(VioGlob)
summary(VioGlob)
#Data is paired
wilcox.test(VioGlob$Control, VioGlob$Lorenzo, paired = TRUE, alternative = "two.sided")
#wilcox.test(VioGlob$Control, VioGlob$Lorenzo, paired = FALSE)
wilcox.test(VioGlob$Control, VioGlob$Mario, paired = TRUE, alternative = "two.sided")
#wilcox.test(VioGlob$Control, VioGlob$Mario, paired = FALSE)
wilcox.test(VioGlob$Control, VioGlob$Madre, paired = TRUE, alternative = "two.sided")
wilcox.test(VioGlob$Control, VioGlob$Padre, paired = TRUE, alternative = "two.sided")


#Individula p.values for 6 Controls
#bedtools intersect -wa -wb -a PCA_Locatelli_family_Controls-dec2018.sep.txt -b ~/ref_av/human_ICR.bed > PCA_Locatelli_family_Controls-dec2018.human_ICR.txt
#awk '{print $17"\t"$18"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' PCA_Locatelli_family_Controls-dec2018.human_ICR.txt > PCA_Locatelli_family_Controls-dec2018.human_ICR.rearranged.txt

#Box and Violin Plot human ICR using 6 Controls from Dec 2018
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli")
Indiv <- read.table("PCA_Locatelli_family_Controls-dec2018.human_ICR.rearranged.txt", header = TRUE)
rownames(Indiv)
colnames(Indiv)
head(Indiv)
Count_IndivICR1 <- count(Indiv, "DMR")
head(Count_IndivICR1)
aggregate1 = aggregate(Indiv[,3:12],by=list(Indiv$DMR), median)
head(aggregate1, 2)
IndivICR1=aggregate1
rownames(IndivICR1)
IndivICR1[,1]
rownames(IndivICR1)=IndivICR1[,1]
rownames(IndivICR1)
colnames(IndivICR1)
IndivICR1 = IndivICR1[,-1]
head(IndivICR1)
dim(IndivICR1)
#write.table(IndivICR1, "aggregated_IndivICR1.txt", sep="\t", quote = FALSE, append = FALSE)
IndivICR1_counted <- cbind(IndivICR1, Count_IndivICR1)
head(IndivICR1_counted)
dim(IndivICR1_counted)
IndivICR1_counted1 <- IndivICR1_counted[which(IndivICR1_counted$freq >= 3),]
head(IndivICR1_counted1)
dim(IndivICR1_counted1)
##write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
IndivICR1_counted2 <- IndivICR1_counted1[,1:10]
head(IndivICR1_counted2)
VioIndivICR <- data.frame(IndivICR1_counted2[,c(5:10,4,3,1,2)])
head(VioIndivICR)
#write.table(VioIndivICR, "VioIndivICR.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
boxplot(VioIndivICR, main="Average methylation at human ICRs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkred", "darkorange", "darkgreen", "blue"))
VioIndivICR2 <- stack(VioIndivICR)
head(VioIndivICR2)
colnames(VioIndivICR2) <- c("Methylation", "Datasets")
head(VioIndivICR2)
ggplot(VioIndivICR2, aes(x=Datasets, y=Methylation, color=Datasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("black","darkblue", "violet", "darkred", "green","pink","blue", "magenta", "red", "orange"))
ggsave("/media/ankitv/Archivio2/ankit/BWS_array/Median/Violin_plot_Locatelli_humanICR_6Indivi.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)


###########################################################################################
###########################################################################################
################################ Mackay   ################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################


#https://rdrr.io/bioc/lumi/man/m2beta.html
library(lumi)
#MacKay
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/Mackay")
Mack <- read.table("NLRP5_mut_mackay_Mvalues.txt", header = TRUE)
head(Mack)
dim(Mack)
rownames(Mack)
Mack[,1]
rownames(Mack)=Mack[,1]
rownames(Mack)
colnames(Mack)
Mack = Mack[,-1]
head(Mack)
Mackm2beta = m2beta(Mack)
head(Mackm2beta)
#write.table(Mackm2beta, "NLRP5_mut_mackay.Betavalues.txt", sep="\t", quote = FALSE, append = FALSE)
#Use excel to unmerged cells and fill the coordinates of hypo and hyper regions 
#bedtools intersect -wa -wb -a NLRP5_mut_mackay.Betavalues.cords.sep.txt -b ~/ref_av/human_ICR.bed > NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.txt
#Rearrange data and save as NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.allFamilyin.txt
#Separate families
#grep F1P1 NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.allFamilyin.txt > NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F1P1.txt
#grep F1P2 NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.allFamilyin.txt > NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F1P2.txt
#grep F2P1 NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.allFamilyin.txt > NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F2P1.txt
#grep F3P NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.allFamilyin.txt > NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F3P.txt
#grep F5P NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.allFamilyin.txt > NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F5P.txt

humanICRorder <- read.table("humanICRorder.txt", header = FALSE)
colnames(humanICRorder) <- "id"

F1P1 <- read.table("NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F1P1.txt", header = FALSE)
head(F1P1)
aggregate_F1P1 <- aggregate(F1P1[,2:3], by = list(F1P1$V4), median)
head(aggregate_F1P1)
colnames(aggregate_F1P1) <- c("id", "ControlF1P1", "CaseF1P1")
#write.table(aggregate_F1P1, "aggregate_F1P1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
aggregate_F1P1_merge <- merge(humanICRorder, aggregate_F1P1, by="id", all.x=TRUE)
#write.table(aggregate_F1P1_merge, "aggregate_F1P1_merge.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

F1P2 <- read.table("NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F1P2.txt", header = FALSE)
head(F1P2)
aggregate_F1P2 <- aggregate(F1P2[,2:3], by = list(F1P2$V4), median)
head(aggregate_F1P2)
colnames(aggregate_F1P2) <- c("id", "ControlF1P2", "CaseF1P2")
#write.table(aggregate_F1P2, "aggregate_F1P2.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
aggregate_F1P2_merge <- merge(humanICRorder, aggregate_F1P2, by="id", all.x=TRUE)
#write.table(aggregate_F1P2_merge, "aggregate_F1P2_merge.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

F2P1 <- read.table("NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F2P1.txt", header = FALSE)
head(F2P1)
aggregate_F2P1 <- aggregate(F2P1[,2:3], by = list(F2P1$V4), median)
head(aggregate_F2P1)
colnames(aggregate_F2P1) <- c("id", "ControlF2P1", "CaseF2P1")
#write.table(aggregate_F2P1, "aggregate_F2P1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
aggregate_F2P1_merge <- merge(humanICRorder, aggregate_F2P1, by="id", all.x=TRUE)
#write.table(aggregate_F2P1_merge, "aggregate_F2P1_merge.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

F3P <- read.table("NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F3P.txt", header = FALSE)
head(F3P)
aggregate_F3P <- aggregate(F3P[,2:3], by = list(F3P$V4), median)
head(aggregate_F3P)
colnames(aggregate_F3P) <- c("id", "ControlF3P", "CaseF3P")
#write.table(aggregate_F3P, "aggregate_F3P.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
aggregate_F3P_merge <- merge(humanICRorder, aggregate_F3P, by="id", all.x=TRUE)
#write.table(aggregate_F3P_merge, "aggregate_F3P_merge.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

F5P <- read.table("NLRP5_mut_mackay.Betavalues.cords.sep.humanICR.F5P.txt", header = FALSE)
head(F5P)
aggregate_F5P <- aggregate(F5P[,2:3], by = list(F5P$V4), median)
head(aggregate_F5P)
#write.table(aggregate_F5P_merge, "aggregate_F5P_merge.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
colnames(aggregate_F5P) <- c("id", "ControlF5P", "CaseF5P")
#write.table(aggregate_F5P, "aggregate_F5P.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
aggregate_F5P_merge <- merge(humanICRorder, aggregate_F5P, by="id", all.x=TRUE)
#write.table(aggregate_F5P_merge, "aggregate_F5P_merge.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

TableMack <- data.frame(cbind(aggregate_F1P1_merge, aggregate_F1P2_merge, aggregate_F2P1_merge, aggregate_F3P_merge, aggregate_F5P_merge))
head(TableMack)
TableMack1 <- TableMack[,c(1:3, 5:6,8:9,11:12,14:15)]
#write.table(TableMack1, "TableMack1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)


#Load data of new sample sent by Mackay: Luciano Processed th data. He lablled the header wrongly so I removed prefix "Mother_"from each samples
#Data was rearanged using excel
MACKAY_BEGEM <- read.csv("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/Mackay/Data_for_DJGM/luciano_processed/wetransfer-8eb990/BETA_DATASET_MACKAY_BEGEM.csv", header = TRUE)
head(MACKAY_BEGEM)
MACKAY_BEGEM_re <- MACKAY_BEGEM[,c(20:21, 21, 1:9)]
head(MACKAY_BEGEM_re)
colnames(MACKAY_BEGEM_re) <- c("Chr","Start","End","ProbeID","AO","JL","MW","NS","DB","WA","SF","RM")
head(MACKAY_BEGEM_re)
write.table(MACKAY_BEGEM_re, "/media/ankitv/Archivio2/ankit/BWS_array/Median/MACKAY_BEGEM_re.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MACKAY_BEGEM_re.txt > MACKAY_BEGEM_rechr.txt
#Remove header manually MACKAY_BEGEM_rechr.txt
##Extract Mackay_begem Probes: awk '{print $4}' MACKAY_BEGEM_re.txt | grep cg > MACKAY_BEGEM_re.ProbeID
#bedtools intersect -wa -wb -a MACKAY_BEGEM_rechr.txt -b ~/ref_av/human_ICR.bed > MACKAY_BEGEM_rechr_humanICR.txt

############################################ Data #################################################
#Prep Mackay data (Luciano Analysis)
#sort -k1,1 -u Bens_EPIC.humanICR.ProbeID | wc -l 690
#Overlap Bens_EPIC on Mackay_Begem probeID
#fgrep -f Bens_EPIC.humanICR.ProbeID MACKAY_BEGEM_re.ProbeID >  Bens_EPIC_MACKAY_BEGEM_re.ProbeID
#fgrep -f Bens_EPIC.humanICR.ProbeID MACKAY_BEGEM_re.ProbeID | wc -l 621
#Extract shared probes
#Mackay_Begem
#fgrep -f Bens_EPIC_MACKAY_BEGEM_re.ProbeID MACKAY_BEGEM_rechr_humanICR.txt  -w  | sort -k1,1 -k2,2n > MACKAY_BEGEM_rechr_humanICR_shared_Probed_overlappedMACKAY.txt
#710 MACKAY_BEGEM_rechr_humanICR.txt
#625 MACKAY_BEGEM_rechr_humanICR_shared_Probed_overlappedMACKAY.txt
#Locatelli
#Extract Locatelli + Control Data with Shared Probes
#fgrep -f Bens_EPIC_MACKAY_BEGEM_re.ProbeID PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR.txt  -w  | sort -k1,1 -k2,2n > PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR_shared_Probed_overlappedMACKAY.txt
#766 PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR.txt
#625 PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR_shared_Probed_overlappedMACKAY.txt
#Bens
#Extract Bens Data with shared Probes
#fgrep -f Bens_EPIC_MACKAY_BEGEM_re.ProbeID Bens/Samples/GSM_Bens_Betavalues_PCA.sort.chr_humanICR.txt   -w  | sort -k1,1 -k2,2n > GSM_Bens_Betavalues_PCA.sort.chr_humanICR_shared_Probe_overlappedMACKAY.txt
#789 Bens/Samples/GSM_Bens_Betavalues_PCA.sort.chr_humanICR.txt
#625 GSM_Bens_Betavalues_PCA.sort.chr_humanICR_shared_Probe_overlappedMACKAY.txt

setwd("/media/ankitv/Archivio2/ankit/BWS_array")
rearranged_Locatelli <- read.table("PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR_shared_Probed_overlappedMACKAY.txt", header = FALSE)
head(rearranged_Locatelli)
head(rearranged_Locatelli[,c(20,6:15)])
rearranged_Locatelli <- rearranged_Locatelli[,c(20,6:15)]
head(rearranged_Locatelli, 10)
Count_rearranged_Locatelli <- count(rearranged_Locatelli, "V20")
head(Count_rearranged_Locatelli)
aggregate_Locatelli = aggregate(rearranged_Locatelli[,2:11],by=list(rearranged_Locatelli$V20), median)
head(aggregate_Locatelli, 3)
aggregate_Locatelli_counted <- cbind(aggregate_Locatelli, Count_rearranged_Locatelli)
head(aggregate_Locatelli_counted)
dim(aggregate_Locatelli_counted)
aggregate_Locatelli_counted1 <- aggregate_Locatelli_counted[which(aggregate_Locatelli_counted$freq >= 3),]
head(aggregate_Locatelli_counted1)
dim(aggregate_Locatelli_counted1)
colnames(aggregate_Locatelli_counted1) <- c("DMR", "Lorenzo", "Mario", "Madre", "Padre", "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5", "CTRL6", "DMR", "Freq")
head(aggregate_Locatelli_counted1)

rearranged_bens <- read.table("GSM_Bens_Betavalues_PCA.sort.chr_humanICR_shared_Probe_overlappedMACKAY.txt", header = FALSE)
head(rearranged_bens)
head(rearranged_bens[,c(104,5:99)])
rearranged_bens <- rearranged_bens[,c(104,5:99)]
head(rearranged_bens, 1)
Count_rearranged_bens <- count(rearranged_bens, "V104")
head(Count_rearranged_bens)
aggregate_bens = aggregate(rearranged_bens[,2:96],by=list(rearranged_bens$V104), median)
head(aggregate_bens, 1)
aggregate_bens_counted <- cbind(aggregate_bens, Count_rearranged_bens)
head(aggregate_bens_counted)
tail(aggregate_bens_counted)
dim(aggregate_bens_counted)
aggregate_bens_counted1 <- aggregate_bens_counted[which(aggregate_bens_counted$freq >= 3),]
head(aggregate_bens_counted1)
dim(aggregate_bens_counted1)
colnames(aggregate_bens_counted1) <- c("DMR", "ID1_TNDM","ID2_1_SRS","ID2_2_SRS","ID3_BWS.H19.Hyper","ID4_SRS","ID5_1_BWS.H19.Hyper","ID5_2_BWS.H19.Hyper","ID6_BWS.Kcnq1ot1.Hypo","ID7_MLID.BWS","ID8_MLID.SRS","ID9_AS","ID10_PWS","ID11_PWS","ID12_PWS","ID13_AS.mosaic","ID14_AS.mosaic","ID15_AS","ID16_AS","ID17_SRS","ID18_BWS.Kcnq1ot1.Hypo","ID19_1_BWS.Kcnq1ot1.Hypo","ID19_2_BWS.Kcnq1ot1.Hypo","ID20_1_BWS.UPD","ID20_2_BWS.UPD","ID21_MLID.SRS","ID22_1_TS.MILD","ID22_2_TS.MILD","ID22_3_TS.MILD","ID23_MLID.BWS","ID24_MLID.BWS","ID25_1_MLID.SRS","ID25_2_MLID.SRS","ID26_MLID.BWS","ID27_MLID.BWS","ID28_AS.mosaic.MILD","ID29_MLID.BWS","ID30_MLID.BWS","ID31_MLID.SRS","ID32_MLID.BWS","ID33_MLID.BWS","ID34_1_MLID.SRS","ID34_2_MLID.SRS","ID35_1_MLID.SRS","ID35_2_MLID.SRS","ID36_1_MLID.PHPIb","ID36_2_MLID.PHPIb","ID37_1_MLID.SRS.&.BWS","ID37_2_MLID.SRS.&.BWS","CB1","CB2","CB3","CB4","CB5","CB6","CB7","CB8","CB9","CB10","CB11","CB12","CB13","CB14","CB15","CB16","CB17","CB18","CB19","CB20","CB21","CB22","CB23","CB24","CB25","CB26","CB27","CB28","CB29","CB30","CB31","CB32_1","CB32_2","CB33_1","CB33_2","CB33_3","CB34_1","CB34_2","CB35_1","CB35_2","CB36_1","CB36_2","CB37_1","CB37_2","CB38_2","CB38_1","CB39","DMR","Freq")
head(aggregate_bens_counted1)

rearranged_mackay_begem <- read.table("MACKAY_BEGEM_rechr_humanICR_shared_Probed_overlappedMACKAY.txt", header = FALSE)
head(rearranged_mackay_begem)
head(rearranged_mackay_begem[,c(17,5:12)])
rearranged_mackay_begem <- rearranged_mackay_begem[,c(17,5:12)]
head(rearranged_mackay_begem, 10)
Count_rearranged_mackay_begem <- count(rearranged_mackay_begem, "V17")
head(Count_rearranged_mackay_begem)
aggregate_mackay_begem = aggregate(rearranged_mackay_begem[,2:9],by=list(rearranged_mackay_begem$V17), median)
head(aggregate_mackay_begem, 3)
aggregate_mackay_begem_counted <- cbind(aggregate_mackay_begem, Count_rearranged_mackay_begem)
head(aggregate_mackay_begem_counted)
dim(aggregate_mackay_begem_counted)
aggregate_mackay_begem_counted1 <- aggregate_mackay_begem_counted[which(aggregate_mackay_begem_counted$freq >= 3),]
head(aggregate_mackay_begem_counted1)
dim(aggregate_mackay_begem_counted1)
colnames(aggregate_mackay_begem_counted1) <- c("DMR","AO","JL","MW","NS","DB","WA","SF","RM","DMR","Freq")
head(aggregate_mackay_begem_counted1)

Merged_LBK <- cbind(aggregate_Locatelli_counted1, aggregate_bens_counted1, aggregate_mackay_begem_counted1)
head(Merged_LBK)
dim(Merged_LBK)
#write.table(Merged_LBK , "Merged_LBK.txt", sep="\t", quote = FALSE, append = FALSE)
Merged_LBK1 <- Merged_LBK[,c(1:11, 15:109, 113:120)]
head(Merged_LBK1)
#write.table(Merged_LBK1, "Merged_LBK1.txt", sep="\t", quote = FALSE, append = FALSE)
#Copy in excel and take the average of duplicated samples and save as Merged_LBK1_dedup.txt

################################### Heatmap human ICR ######################################àà
library(ggplot2)
library(gplots)
head(Merged_LBK1,1)
rownames(Merged_LBK1)
Merged_LBK1[,1]
rownames(Merged_LBK1)=Merged_LBK1[,1]
rownames(Merged_LBK1)
colnames(Merged_LBK1)
Merged_LBK1 = Merged_LBK1[,-1]
head(Merged_LBK1)
colnames(Merged_LBK1)
Merged_LBK1 = as.matrix(Merged_LBK1)
dim(Merged_LBK1)
head(Merged_LBK1,1)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
#Col cluster
svg(filename="Heatmap_Merged_LBK1dedup.reorderICR_heamap.svg", width=10, height=10, pointsize=12)
heatmap.2(Merged_LBK1,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="column", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(0,1, length.out = 100))
dev.off()
Merged_LBK1 <- data.frame(Merged_LBK1)
#Take average here instead of excel sheet
Merged_LBK2 <- data.frame((cbind(Merged_LBK1[,c(1:11)], 
                                 rowMeans(Merged_LBK1[,12:13]), 
                                 Merged_LBK1[,c(14:15)], 
                                 rowMeans(Merged_LBK1[,16:17]), 
                                 Merged_LBK1[,c(18:30)], 
                                 rowMeans(Merged_LBK1[,31:32]), 
                                 rowMeans(Merged_LBK1[,33:34]), 
                                 Merged_LBK1[,c(35)], 
                                 rowMeans(Merged_LBK1[,36:38]), 
                                 Merged_LBK1[,c(39:40)], 
                                 rowMeans(Merged_LBK1[,41:42]), 
                                 Merged_LBK1[,c(43:50)], 
                                 rowMeans(Merged_LBK1[,51:52]), 
                                 rowMeans(Merged_LBK1[,53:54]), 
                                 rowMeans(Merged_LBK1[,55:56]), 
                                 rowMeans(Merged_LBK1[,57:58]), 
                                 Merged_LBK1[,c(59:89)], 
                                 rowMeans(Merged_LBK1[,90:91]), 
                                 rowMeans(Merged_LBK1[,92:94]), 
                                 rowMeans(Merged_LBK1[,95:96]), 
                                 rowMeans(Merged_LBK1[,97:98]), 
                                 rowMeans(Merged_LBK1[,99:100]), 
                                 rowMeans(Merged_LBK1[,101:102]), 
                                 rowMeans(Merged_LBK1[,103:104]), 
                                 Merged_LBK1[,c(105:113)])))
head(Merged_LBK2)
#write.table(Merged_LBK2, "Merged_LBK2.txt", sep="\t", quote = FALSE, append = FALSE)
colnames(Merged_LBK2) <- c("Lorenzo", "Mario", "Madre", "Padre", "CL1", "CL2", "CL3", "CL4", "CL5", "CL6", "ID1_TNDM","ID2_SRS","ID3_BWS.H19.Hyper","ID4_SRS","ID5_BWS.H19.Hyper","ID6_BWS.Kcnq1ot1.Hypo","ID7_MLID.BWS","ID8_MLID.SRS","ID9_AS","ID10_PWS","ID11_PWS","ID12_PWS","ID13_AS.mosaic","ID14_AS.mosaic","ID15_AS","ID16_AS","ID17_SRS","ID18_BWS.Kcnq1ot1.Hypo","ID19_BWS.Kcnq1ot1.Hypo","ID20_BWS.UPD","ID21_MLID.SRS","ID22_TS.MILD","ID23_MLID.BWS","ID24_MLID.BWS","ID25_MLID.SRS","ID26_MLID.BWS","ID27_MLID.BWS","ID28_AS.mosaic.MILD","ID29_MLID.BWS","ID30_MLID.BWS","ID31_MLID.SRS","ID32_MLID.BWS","ID33_MLID.BWS","ID34_MLID.SRS","ID35_MLID.SRS","ID36_MLID.PHPIb","ID37_MLID.SRS.BWS","CB1","CB2","CB3","CB4","CB5","CB6","CB7","CB8","CB9","CB10","CB11","CB12","CB13","CB14","CB15","CB16","CB17","CB18","CB19","CB20","CB21","CB22","CB23","CB24","CB25","CB26","CB27","CB28","CB29","CB30","CB31","CB32","CB33","CB34","CB35","CB36","CB37","CB38","CB39","AO","JL","MW","NS","DB","WA","SF","RM")
head(Merged_LBK2)
#write.table(Merged_LBK2, "Merged_LBK2dedup.txt", sep="\t", quote = FALSE, append = FALSE)

#Extract MLID cases from Bens
head(Merged_LBK2)
Merged_LBK2mlid <- as.matrix(Merged_LBK2[,c(1:2, 5:10, 31:94)])
head(Merged_LBK2mlid)
#Heatmap human ICR merged, dedup data
svg(filename="/media/ankitv/Archivio2/ankit/BWS_array/Median/Heatmap_Merged_LBK1mlid_dedup.svg", width=15, height=10, pointsize=12)
heatmap.2(Merged_LBK2mlid,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="column", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(0,1, length.out = 100))
dev.off()
library(pheatmap)
my_sample_col <- data.frame(Annotations= c("A_Lorenzo","A_Mario","yControl","yControl","yControl","yControl","yControl","yControl","B_MLID.SRS","B_MLID.TS","B_MLID.BWS","B_MLID.BWS","B_MLID.SRS","B_MLID.BWS","B_MLID.BWS","B_MLID.AS.mosaic","B_MLID.BWS","B_MLID.BWS","B_MLID.SRS","B_MLID.BWS","B_MLID.BWS","B_MLID.SRS","B_MLID.SRS","B_MLID.PHPIb","B_MLID.SRS.BWS","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","C_AO","C_JL","C_MW","C_NS","C_DB","C_WA","C_SF","C_RM"))
row.names(my_sample_col) <- colnames(Merged_LBK2mlid)
#my_sample_col$Annotations <- factor(my_sample_col$Annotations, levels = c("A_Lorenzo","A_Mario","zControl","zControl","zControl","zControl","zControl","zControl","B_MLID/SRS","B_MLID/TS","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/BWS","B_MLID/BWS","B_MLID/AS.mosaic","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/SRS","B_MLID/PHPIb","B_MLID/SRS.BWS","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","C_AO","C_JL","C_MW","C_NS","C_DB","C_WA","C_SF","C_RM"))
my_colour = list(Annotations = c(A_Lorenzo = "pink", A_Mario = "#FA8072", B_MLID.AS.mosaic = "Red", B_MLID.BWS = "Orange", B_MLID.PHPIb = "#00FFFF", B_MLID.SRS = "#8B0000", B_MLID.SRS.BWS = "#CC99FF", B_MLID.TS = "#666600", C_AO = "#99F38F", C_DB = "#2CC71B", C_JL = "#cccc00", C_MW = "#993399", C_NS = "#cc6600", C_RM = "#CC0000", C_SF = "#006600", C_WA = "#FFCC00", yControl = "#808080", zControl = "#D3D3D3"))
breaksList = seq(0, 1, by =0.01)
pheatmap(Merged_LBK2mlid,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)),
         annotation_colors = my_colour,
         fontsize = 8,
         annotation_col = my_sample_col,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)
#################################################### Transform data #####################################à

Merged_LBK3mlid <- as.matrix(Merged_LBK2[,c(1:10, 31:94)])
tMerged_LBK3mlid = t(Merged_LBK3mlid)
head(tMerged_LBK3mlid)
dim(tMerged_LBK3mlid)
rownames(tMerged_LBK3mlid)
tMerged_LBK3mlid = data.frame(tMerged_LBK3mlid)
head(tMerged_LBK3mlid)
#write.table(tMerged_LBK3mlid , "tMerged_LBK3mliddedup.txt", sep="\t", quote = FALSE, append = FALSE)
tMerged_LBK3mlid["Color"] <- c("A_Lorenzo","A_Mario","A_Mradre", "A_Padre","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","B_MLID/SRS","B_MLID/TS","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/BWS","B_MLID/BWS","B_MLID/AS.mosaic","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/SRS","B_MLID/PHPIb","B_MLID/SRS.BWS","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","A_Control_Bens","C_AO","C_JL","C_MW","C_NS","C_DB","C_WA","C_SF","C_RM")
head(tMerged_LBK3mlid)
dim(tMerged_LBK3mlid)
dfx <-tMerged_LBK3mlid[c(1:41)]
head(dfx)
PC<-prcomp(dfx)
head(PC)
PCi<-data.frame(PC$x,Color=tMerged_LBK3mlid$Color)
percentageICR <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentageICR <- paste( colnames(PCi), "(", paste( as.character(percentageICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

#df <- as.data.frame(dataMergedLBJM1)
p<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentageICR[1]) + ylab(percentageICR[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("#D3D3D3","#808080","pink","#FA8072","Magenta","Blue","Red","Orange","#00FFFF","#8B0000","#CC99FF","#666600","#99F38F","#2CC71B", "#cccc00", "#993399", "#cc6600", "#CC0000", "#006600", "#FFCC00"))+
  scale_shape_manual(values=c(1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
p <- p+theme_bw()
p
ggsave("/media/ankitv/Archivio2/ankit/BWS_array/Median/PCA-TMergedLBK.humanICR.averageddedup.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96)


#F.Rezwan Analysis
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Locatelli/Mackay/Data_for_DJGM/rezwan_processed/Data")
load("beta_QN_autosome.studyID.combat.RData")
head(beta)
dim(beta)
colnames(beta) <- c("NS", "AO", "JL", "CTRL_1_1", "CTRL001b", "CTRL003", "DB", "RM", "SF", "WA", "MW")
head(beta)
beta1 <- data.frame(cbind(rownames(beta), beta[,c(4:6, 1:3, 7:11)]))
beta1 <- beta1[order(beta1$V1),]
colnames(beta1) <- c("Probe", "CTRL_1_1", "CTRL001b", "CTRL003", "NS", "AO", "JL", "DB", "RM", "SF", "WA", "MW")
head(beta1)
dim(beta1)
#write.table(beta1[,1], "beta1.id", sep="\t", quote = FALSE, append = FALSE, row.names = F)
#write.table(beta1, "beta1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)


#LOAD MANIFEST "cg" probes only
#grep cg HumanMethylation450_15017482_v1-2.csv > HumanMethylation450_15017482_v1-2_cg.csv
manifest <- read.csv("HumanMethylation450_15017482_v1-2_cg.csv", header = FALSE)
head(manifest, 30)
dim(manifest)
MANIFEST <- manifest[,c(1, 12:13, 13)]
head(MANIFEST)
dim(MANIFEST)
MANIFEST <- MANIFEST[order(MANIFEST$V1),]
head(MANIFEST, 50)
colnames(MANIFEST) <- c("Probe", "Chr", "Start", "End")
#write.table(MANIFEST, "MANIFEST.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)
#write.table(MANIFEST[,1], "MANIFEST.id", sep="\t", quote = FALSE, append = FALSE, row.names = F)
#awk '{print $1"\t""chr"$2"\t"$3"\t"$4}' MANIFEST.txt > MANIFEST_chr.txt
#Extract shared ids between beta-data and manifest file, fgrep -f beta1.id MANIFEST.id > shared_beat1_manifest.id
#Extract manifest info using shared probes: fgrep -f shared_beat1_manifest.id MANIFEST_chr.txt > MANIFEST_chr_sharedprobesb.txt
#Extract beta-data info using shared probes: fgrep -f shared_beat1_manifest.id beta1.txt > beta1_sharedprobesb.txt
beta1_sharedprobesb =read.table("beta1_sharedprobesb.txt", header = FALSE, stringsAsFactors = FALSE)
MANIFEST_chr_sharedprobesb =read.table("MANIFEST_chr_sharedprobesb.txt", header = FALSE, stringsAsFactors = FALSE)
rezwan_beta1_manifest <- cbind(MANIFEST_chr_sharedprobesb, beta1_sharedprobesb)
head(rezwan_beta1_manifest)
tail(rezwan_beta1_manifest)
dim(rezwan_beta1_manifest)
count(data.frame(intersect(rezwan_beta1_manifest[,1], rezwan_beta1_manifest[,5]))) #Check counts oerlap of probe, All overlapped and aligned
rezwan_beta1_manifest <- rezwan_beta1_manifest[,c(2:16)]
colnames(rezwan_beta1_manifest) <- c("Chr","Star","End","Probe", "CTRL_1_1", "CTRL001b", "CTRL003", "NS", "AO", "JL", "DB", "RM", "SF", "WA", "MW")
head(rezwan_beta1_manifest)
#write.table(rezwan_beta1_manifest, "rezwan_beta1_manifest.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F, col.names = F)
##Extract rezwan_beta1 Probes: awk '{print $4}' rezwan_beta1_manifest.txt > rezwan_beta1_manifest.ProbeID.
#Note rezwan_beta1_manifest.ProbeID is same as shared_beat1_manifest.id
#bedtools intersect -wa -wb -a rezwan_beta1_manifest.txt -b ~/ref_av/human_ICR.bed > rezwan_beta1_manifest_humanICR.txt
#Prep Mackay data (Rezwan Analysis)
#See R script for derivation of file  rezwan_beta1_manifest.ProbeID
#cp Locatelli/Mackay/Data_for_DJGM/rezwan_processed/Data/rezwan_beta1_manifest.ProbeID ./
#Overlap Bens_EPIC on rezwan_beta1 probeID
#fgrep -f Bens_EPIC.humanICR.ProbeID rezwan_beta1_manifest.ProbeID >  Bens_EPIC_rezwan_beta1_manifest.ProbeID #688
#fgrep -f Bens_EPIC.humanICR.ProbeID rezwan_beta1_manifest.ProbeID | #wc -l
#Extract shared probes
#rezwan_beta1
#fgrep -f Bens_EPIC_rezwan_beta1_manifest.ProbeID rezwan_beta1_manifest_humanICR.txt -w | sort -k1,1 -k2,2n > rezwan_beta1_rechr_humanICR_shared_Probed_overlappedrezwan.txt
#wc -l rezwan_beta1_manifest_humanICR.txt  , 785 rezwan_beta1_manifest_humanICR.txt
#wc -l rezwan_beta1_rechr_humanICR_shared_Probed_overlappedrezwan.txt, 693
#Locatelli
#Extract Locatelli + Control Data with Shared Probes
#fgrep -f Bens_EPIC_rezwan_beta1_manifest.ProbeID PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR.txt -w | sort -k1,1 -k2,2n > PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR_shared_Probed_overlappedrezwan.txt
#wc -l PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR.txt, 766
#wc -l PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR_shared_Probed_overlappedrezwan.txt, 693
#Bens
#Extract Bens Data with shared Probes
#fgrep -f Bens_EPIC_rezwan_beta1_manifest.ProbeID Bens/Samples/GSM_Bens_Betavalues_PCA.sort.chr_humanICR.txt -w | sort -k1,1 -k2,2n > GSM_Bens_Betavalues_PCA.sort.chr_humanICR_shared_Probe_overlappedrezwan.txt
#wc -l Bens/Samples/GSM_Bens_Betavalues_PCA.sort.chr_humanICR.txt, 789
#wc -l GSM_Bens_Betavalues_PCA.sort.chr_humanICR_shared_Probe_overlappedrezwan.txt, 693

############################################ Data #################################################
#Import ICR intersected probes for each of the samples
setwd("/media/ankitv/Archivio2/ankit/BWS_array")
rearranged_Locatelli_rezwan <- read.table("PCA_Locatelli_family_Controls-dec2018.chrPos_humanICR_shared_Probed_overlappedrezwan.txt", 
                                          header = FALSE)
head(rearranged_Locatelli_rezwan)
head(rearranged_Locatelli_rezwan[,c(20,6:15)])
rearranged_Locatelli_rezwan <- rearranged_Locatelli_rezwan[,c(20,6:15)]
head(rearranged_Locatelli_rezwan, 10)
Count_rearranged_Locatelli_rezwan <- count(rearranged_Locatelli_rezwan, "V20")
head(Count_rearranged_Locatelli_rezwan)
aggregate_Locatelli_rezwan = aggregate(rearranged_Locatelli_rezwan[,2:11],by=list(rearranged_Locatelli_rezwan$V20), median)
head(aggregate_Locatelli_rezwan, 3)
aggregate_Locatelli_rezwan_counted <- cbind(aggregate_Locatelli_rezwan, Count_rearranged_Locatelli_rezwan)
head(aggregate_Locatelli_rezwan_counted)
dim(aggregate_Locatelli_rezwan_counted)
aggregate_Locatelli_rezwan_counted1 <- aggregate_Locatelli_rezwan_counted[which(aggregate_Locatelli_rezwan_counted$freq >= 3),]
head(aggregate_Locatelli_rezwan_counted1)
dim(aggregate_Locatelli_rezwan_counted1)
colnames(aggregate_Locatelli_rezwan_counted1) <- c("DMR", "Lorenzo", "Mario", "Madre", "Padre", "CTRL1", "CTRL2", "CTRL3", "CTRL4", "CTRL5", "CTRL6", "DMR", "Freq")
head(aggregate_Locatelli_rezwan_counted1)

rearranged_bens_rezwan <- read.table("GSM_Bens_Betavalues_PCA.sort.chr_humanICR_shared_Probe_overlappedrezwan.txt", header = FALSE)
head(rearranged_bens_rezwan)
head(rearranged_bens_rezwan[,c(104,5:99)])
rearranged_bens_rezwan <- rearranged_bens_rezwan[,c(104,5:99)]
head(rearranged_bens_rezwan, 1)
Count_rearranged_bens_rezwan <- count(rearranged_bens_rezwan, "V104")
head(Count_rearranged_bens_rezwan)
aggregate_bens_rezwan = aggregate(rearranged_bens_rezwan[,2:96],by=list(rearranged_bens_rezwan$V104), median)
head(aggregate_bens_rezwan, 1)
aggregate_bens_rezwan_counted <- cbind(aggregate_bens_rezwan, Count_rearranged_bens_rezwan)
head(aggregate_bens_rezwan_counted)
tail(aggregate_bens_rezwan_counted)
dim(aggregate_bens_rezwan_counted)
aggregate_bens_rezwan_counted1 <- aggregate_bens_rezwan_counted[which(aggregate_bens_rezwan_counted$freq >= 3),]
head(aggregate_bens_rezwan_counted1)
dim(aggregate_bens_rezwan_counted1)
colnames(aggregate_bens_rezwan_counted1) <- c("DMR", "ID1_TNDM","ID2_1_SRS","ID2_2_SRS","ID3_BWS.H19.Hyper","ID4_SRS","ID5_1_BWS.H19.Hyper","ID5_2_BWS.H19.Hyper","ID6_BWS.Kcnq1ot1.Hypo","ID7_MLID.BWS","ID8_MLID.SRS","ID9_AS","ID10_PWS","ID11_PWS","ID12_PWS","ID13_AS.mosaic","ID14_AS.mosaic","ID15_AS","ID16_AS","ID17_SRS","ID18_BWS.Kcnq1ot1.Hypo","ID19_1_BWS.Kcnq1ot1.Hypo","ID19_2_BWS.Kcnq1ot1.Hypo","ID20_1_BWS.UPD","ID20_2_BWS.UPD","ID21_MLID.SRS","ID22_1_TS.MILD","ID22_2_TS.MILD","ID22_3_TS.MILD","ID23_MLID.BWS","ID24_MLID.BWS","ID25_1_MLID.SRS","ID25_2_MLID.SRS","ID26_MLID.BWS","ID27_MLID.BWS","ID28_AS.mosaic.MILD","ID29_MLID.BWS","ID30_MLID.BWS","ID31_MLID.SRS","ID32_MLID.BWS","ID33_MLID.BWS","ID34_1_MLID.SRS","ID34_2_MLID.SRS","ID35_1_MLID.SRS","ID35_2_MLID.SRS","ID36_1_MLID.PHPIb","ID36_2_MLID.PHPIb","ID37_1_MLID.SRS.&.BWS","ID37_2_MLID.SRS.&.BWS","CB1","CB2","CB3","CB4","CB5","CB6","CB7","CB8","CB9","CB10","CB11","CB12","CB13","CB14","CB15","CB16","CB17","CB18","CB19","CB20","CB21","CB22","CB23","CB24","CB25","CB26","CB27","CB28","CB29","CB30","CB31","CB32_1","CB32_2","CB33_1","CB33_2","CB33_3","CB34_1","CB34_2","CB35_1","CB35_2","CB36_1","CB36_2","CB37_1","CB37_2","CB38_2","CB38_1","CB39","DMR","Freq")
head(aggregate_bens_rezwan_counted1)

rearranged_rezwan_beta1 <- read.table("rezwan_beta1_rechr_humanICR_shared_Probed_overlappedrezwan.txt", header = FALSE)
head(rearranged_rezwan_beta1)
head(rearranged_rezwan_beta1[,c(20,5:15)])
rearranged_rezwan_beta1 <- rearranged_rezwan_beta1[,c(20,5:15)]
head(rearranged_rezwan_beta1, 10)
Count_rearranged_rezwan_beta1 <- count(rearranged_rezwan_beta1, "V20")
head(Count_rearranged_rezwan_beta1)
aggregate_rezwan_beta1 = aggregate(rearranged_rezwan_beta1[,2:12],by=list(rearranged_rezwan_beta1$V20), median)
head(aggregate_rezwan_beta1, 3)
aggregate_rezwan_beta1_counted <- cbind(aggregate_rezwan_beta1, Count_rearranged_rezwan_beta1)
head(aggregate_rezwan_beta1_counted)
dim(aggregate_rezwan_beta1_counted)
aggregate_rezwan_beta1_counted1 <- aggregate_rezwan_beta1_counted[which(aggregate_rezwan_beta1_counted$freq >= 3),]
head(aggregate_rezwan_beta1_counted1)
dim(aggregate_rezwan_beta1_counted1)
colnames(aggregate_rezwan_beta1_counted1) <- c("DMR","CTRL_1_1","CTRL001b","CTRL003","NS","AO","JL","DB","RM","SF","WA","MW","DMR","Freq")
head(aggregate_rezwan_beta1_counted1)

Merged_LBR <- cbind(aggregate_Locatelli_rezwan_counted1, aggregate_bens_rezwan_counted1, aggregate_rezwan_beta1_counted1)
head(Merged_LBR)
dim(Merged_LBR)
write.table(Merged_LBR , "/media/ankitv/Archivio2/ankit/BWS_array/Median/Merged_LBR.txt", sep="\t", quote = FALSE, append = FALSE)
Merged_LBR1 <- Merged_LBR[,c(1:11, 15:109, 113:123)]
head(Merged_LBR1)
write.table(Merged_LBR1, "/media/ankitv/Archivio2/ankit/BWS_array/Median/Merged_LBR1.txt", sep="\t", quote = FALSE, append = FALSE)

################################### Heatmap human ICR ######################################
head(Merged_LBR1,1)
rownames(Merged_LBR1)
Merged_LBR1[,1]
rownames(Merged_LBR1)=Merged_LBR1[,1]
rownames(Merged_LBR1)
colnames(Merged_LBR1)
Merged_LBR1 = Merged_LBR1[,-1]
head(Merged_LBR1)
colnames(Merged_LBR1)
Merged_LBR1 = as.matrix(Merged_LBR1)
dim(Merged_LBR1)
head(Merged_LBR1,1)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
#Col cluster
svg(filename="/media/ankitv/Archivio2/ankit/BWS_array/Median/Heatmap_Merged_LBR1dedup.reorderICR_heamap.svg", width=10, height=10, pointsize=12)
heatmap.2(Merged_LBR1,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="column", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(0,1, length.out = 100))
dev.off()
Merged_LBR1 <- data.frame(Merged_LBR1)
#Take average by R instead of excel
Merged_LBR2 <- data.frame((cbind(Merged_LBR1[,c(1:11)], 
                                 rowMeans(Merged_LBR1[,12:13]), 
                                 Merged_LBR1[,c(14:15)], 
                                 rowMeans(Merged_LBR1[,16:17]), 
                                 Merged_LBR1[,c(18:30)], 
                                 rowMeans(Merged_LBR1[,31:32]), 
                                 rowMeans(Merged_LBR1[,33:34]), 
                                 Merged_LBR1[,c(35)], 
                                 rowMeans(Merged_LBR1[,36:38]), 
                                 Merged_LBR1[,c(39:40)], 
                                 rowMeans(Merged_LBR1[,41:42]), 
                                 Merged_LBR1[,c(43:50)], 
                                 rowMeans(Merged_LBR1[,51:52]), 
                                 rowMeans(Merged_LBR1[,53:54]), 
                                 rowMeans(Merged_LBR1[,55:56]), 
                                 rowMeans(Merged_LBR1[,57:58]), 
                                 Merged_LBR1[,c(59:89)], 
                                 rowMeans(Merged_LBR1[,90:91]), 
                                 rowMeans(Merged_LBR1[,92:94]), 
                                 rowMeans(Merged_LBR1[,95:96]), 
                                 rowMeans(Merged_LBR1[,97:98]), 
                                 rowMeans(Merged_LBR1[,99:100]), 
                                 rowMeans(Merged_LBR1[,101:102]), 
                                 rowMeans(Merged_LBR1[,103:104]), 
                                 Merged_LBR1[,c(105:116)])))
head(Merged_LBR2)
#write.table(Merged_LBR2, "Merged_LBR2.txt", sep="\t", quote = FALSE, append = FALSE)
colnames(Merged_LBR2) <- c("Lorenzo", "Mario", "Madre", "Padre", "CL1", "CL2", "CL3", "CL4", "CL5", "CL6", "ID1_TNDM","ID2_SRS","ID3_BWS.H19.Hyper","ID4_SRS","ID5_BWS.H19.Hyper","ID6_BWS.Kcnq1ot1.Hypo","ID7_MLID.BWS","ID8_MLID.SRS","ID9_AS","ID10_PWS","ID11_PWS","ID12_PWS","ID13_AS.mosaic","ID14_AS.mosaic","ID15_AS","ID16_AS","ID17_SRS","ID18_BWS.Kcnq1ot1.Hypo","ID19_BWS.Kcnq1ot1.Hypo","ID20_BWS.UPD","ID21_MLID.SRS","ID22_TS.MILD","ID23_MLID.BWS","ID24_MLID.BWS","ID25_MLID.SRS","ID26_MLID.BWS","ID27_MLID.BWS","ID28_AS.mosaic.MILD","ID29_MLID.BWS","ID30_MLID.BWS","ID31_MLID.SRS","ID32_MLID.BWS","ID33_MLID.BWS","ID34_MLID.SRS","ID35_MLID.SRS","ID36_MLID.PHPIb","ID37_MLID.SRS.BWS","CB1","CB2","CB3","CB4","CB5","CB6","CB7","CB8","CB9","CB10","CB11","CB12","CB13","CB14","CB15","CB16","CB17","CB18","CB19","CB20","CB21","CB22","CB23","CB24","CB25","CB26","CB27","CB28","CB29","CB30","CB31","CB32","CB33","CB34","CB35","CB36","CB37","CB38","CB39","CTRL_1_1","CTRL001b","CTRL003","NS","AO","JL","DB","RM","SF","WA","MW")
head(Merged_LBR2)
#write.table(Merged_LBR2, "Merged_LBR2dedup.txt", sep="\t", quote = FALSE, append = FALSE)

#Extract MLID cases from bens_rezwan
head(Merged_LBR2)
Merged_LBR2mlid <- as.matrix(Merged_LBR2[,c(1:4, 5:10, 31:97)])
head(Merged_LBR2mlid)
write.table(Merged_LBR2mlid, "/media/ankitv/Archivio2/ankit/BWS_array/Median/Merged_LBR2mliddedup.txt", sep="\t", quote = FALSE, append = FALSE)
#Heatmap human ICR merged, dedup data
library(pheatmap)
my_sample_col <- data.frame(Annotations= c("A_Lorenzo","A_Mario","A_Mradre","A_Padre","yControl","yControl","yControl","yControl","yControl","yControl","B_MLID.SRS","B_MLID.TS","B_MLID.BWS","B_MLID.BWS","B_MLID.SRS","B_MLID.BWS","B_MLID.BWS","B_MLID.AS.mosaic","B_MLID.BWS","B_MLID.BWS","B_MLID.SRS","B_MLID.BWS","B_MLID.BWS","B_MLID.SRS","B_MLID.SRS","B_MLID.PHPIb","B_MLID.SRS.BWS","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","zControl","xControl","xControl","xControl","C_NS","C_Mackay_BWS","C_Mackay_BWS","C_DB","C_RM","C_SF","C_WA","C_Mackay_BWS"))
row.names(my_sample_col) <- colnames(Merged_LBR2mlid)
#my_colour = list(Annotations = c(A_Lorenzo = "#CC0000", A_Mario = "darkgreen", B_MLID.AS.mosaic = "Red", B_MLID.BWS = "Orange", B_MLID.PHPIb = "#00FFFF", B_MLID.SRS = "#8B0000", B_MLID.SRS.BWS = "#CC99FF", B_MLID.TS = "#666600", C_AO = "#99F38F", C_DB = "#2CC71B", C_JL = "#cccc00", C_MW = "#993399", C_NS = "#cc6600", C_RM = "#CC0000", C_SF = "#006600", C_WA = "navy", xControl = "grey", yControl = "#808080", zControl = "#D3D3D3"))
my_colour = list(Annotations = c(zControl = "lightgray", yControl = "gray", xControl = "dimgray", A_Lorenzo = "darkred", A_Mario = "orange", A_Mradre = "magenta", A_Padre = "blue", B_MLID.AS.mosaic = "#8cff66", B_MLID.BWS = "mediumseagreen", B_MLID.PHPIb = "darkgreen", B_MLID.SRS = "darkkhaki", B_MLID.SRS.BWS = "#666633", B_MLID.TS = "yellowgreen", C_DB = "#c299ff", C_Mackay_BWS = "#b3ccff", C_NS = "cyan", C_RM = "deepskyblue", C_SF = "#884dff", C_WA = "navy"))
breaksList = seq(0, 1, by =0.01)
pheatmap(Merged_LBR2mlid,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)),
         annotation_colors = my_colour,
         fontsize = 8,
         annotation_col = my_sample_col,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)
#################################################### Transform data #####################################à
head(Merged_LBR2)
Merged_LBR3mlid <- as.matrix(Merged_LBR2[,c(1:10, 31:97)])
tMerged_LBR3mlid = t(Merged_LBR3mlid)
head(tMerged_LBR3mlid)
dim(tMerged_LBR3mlid)
rownames(tMerged_LBR3mlid)
tMerged_LBR3mlid = data.frame(tMerged_LBR3mlid)
head(tMerged_LBR3mlid)
#write.table(tMerged_LBR3mlid , "tMerged_LBR3mliddedup.txt", sep="\t", quote = FALSE, append = FALSE)
#tMerged_LBR3mlid["Color"] <-  c("A_Lorenzo","A_Mario","A_Mradre", "A_Padre","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","B_MLID/SRS","B_MLID/TS","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/BWS","B_MLID/BWS","B_MLID/AS.mosaic","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/SRS","B_MLID/PHPIb","B_MLID/SRS.BWS","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_rezwan","A_Control_rezwan","A_Control_rezwan","C_NS","C_AO","C_JL","C_DB","C_RM","C_SF","C_WA","C_MW")
tMerged_LBR3mlid["Color"] <-  c("A_Lorenzo","A_Mario","A_Mradre", "A_Padre","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","A_Control_Loca","B_MLID/SRS","B_MLID/TS","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/BWS","B_MLID/BWS","B_MLID/AS.mosaic","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/BWS","B_MLID/BWS","B_MLID/SRS","B_MLID/SRS","B_MLID/PHPIb","B_MLID/SRS.BWS","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_bens","A_Control_rezwan","A_Control_rezwan","A_Control_rezwan","C_NS","C_Mackay_BWS","C_Mackay_BWS","C_DB","C_RM","C_SF","C_WA","C_Mackay_BWS")

head(tMerged_LBR3mlid)
write.table(tMerged_LBR3mlid , "/media/ankitv/Archivio2/ankit/BWS_array/Median/tMerged_LBR3mliddedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tMerged_LBR3mlid)
dfx <-tMerged_LBR3mlid[c(1:42)]
head(dfx)
PC<-prcomp(dfx)
head(PC)
PCi<-data.frame(PC$x,Color=tMerged_LBR3mlid$Color)
percentageICR <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentageICR <- paste( colnames(PCi), "(", paste( as.character(percentageICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

#df <- as.data.frame(dataMergedLBJM1)
p<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentageICR[1]) + ylab(percentageICR[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("lightgray","gray","dimgray","darkred","orange","magenta","blue","#8cff66","mediumseagreen","darkgreen","darkkhaki","#666633","yellowgreen","#c299ff","#b3ccff","cyan","deepskyblue","#884dff","navy"))+
  scale_shape_manual(values=c(1, 2, 3,19,19,7,7,8,9,10,10,12,16,14,9,16,17,18,21))
p <- p+theme_classic()
p
ggsave("/media/ankitv/Archivio2/ankit/BWS_array/Median/PCA-TMergedLBR.humanICR.averageddedup2.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96)



#Correlation plot, mean of beta values vs PC1
head(Merged_LBR3mlid)
dat=prcomp(t(as.matrix(Merged_LBR3mlid)))
head(dat)
autoplot(dat,  label = TRUE)


Methyl_mean=colMeans(Merged_LBR3mlid)

PC1corMehyl=cor(dat$x[,1], Methyl_mean, method = "pearson")
PC2corMehyl=cor(dat$x[,2], Methyl_mean, method = "pearson")

A=cbind("PC1" =dat$x[,1], "Methyl_mean"=Methyl_mean)
head(A)
A <- data.frame(A)
library(GGally)
g <- ggpairs(A, upper = list(continuous = wrap("cor", size = 7)),
             diag=list(continuous = wrap("barDiag", bins = 20, fill="navy")))
g
ggsave("/media/ankitv/Archivio2/ankit/BWS_array/Median/correlationLBR_PC1_methylmean.svg", width=12*1.25, height=12*1.25, units="cm", dpi=96)

B=cbind("PC1" =dat$x[,1:2], "Methyl_mean"=Methyl_mean)
head(B)
B <- data.frame(B)
library(GGally)
g <- ggpairs(B, upper = list(continuous = wrap("cor", size = 7)),
             diag=list(continuous = wrap("barDiag", bins = 20, fill="navy")))
g
ggsave("/media/ankitv/Archivio2/ankit/BWS_array/Median/correlationLBR_PC1-PC2_methylmean.svg", width=12*1.25, height=12*1.25, units="cm", dpi=96)


###########################################################################################
###########################################################################################
###########################################################################################
################################# Bens  ###################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

#Bens Data was downloaded from GSE78773 one by one from processed table in GSM codes mentioned
#See excel sheet commands.xlsx in this folder to know how the below files were produced
#Only ProbeID and Beta values were extracted Probes were already sorted

library(ggplot2)
library(gplots)
library(plyr)
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Bens/Samples")
GSM2076276 <- read.table("GSM2076276-30182.Probe.Beta-value.txt", header=FALSE)
GSM2076277 <- read.table("GSM2076277-30183.Probe.Beta-value.txt", header=FALSE)
GSM2076278 <- read.table("GSM2076278-30184.Probe.Beta-value.txt", header=FALSE)
GSM2076279 <- read.table("GSM2076279-30185.Probe.Beta-value.txt", header=FALSE)
GSM2076280 <- read.table("GSM2076280-30186.Probe.Beta-value.txt", header=FALSE)
GSM2076281 <- read.table("GSM2076281-30187.Probe.Beta-value.txt", header=FALSE)
GSM2076282 <- read.table("GSM2076282-30188.Probe.Beta-value.txt", header=FALSE)
GSM2076283 <- read.table("GSM2076283-30189.Probe.Beta-value.txt", header=FALSE)
GSM2076284 <- read.table("GSM2076284-30190.Probe.Beta-value.txt", header=FALSE)
GSM2076285 <- read.table("GSM2076285-30191.Probe.Beta-value.txt", header=FALSE)
GSM2076286 <- read.table("GSM2076286-30192.Probe.Beta-value.txt", header=FALSE)
GSM2076287 <- read.table("GSM2076287-30193.Probe.Beta-value.txt", header=FALSE)
GSM2076288 <- read.table("GSM2076288-30194.Probe.Beta-value.txt", header=FALSE)
GSM2076289 <- read.table("GSM2076289-30195.Probe.Beta-value.txt", header=FALSE)
GSM2076290 <- read.table("GSM2076290-30196.Probe.Beta-value.txt", header=FALSE)
GSM2076291 <- read.table("GSM2076291-30197.Probe.Beta-value.txt", header=FALSE)
GSM2076292 <- read.table("GSM2076292-30198.Probe.Beta-value.txt", header=FALSE)
GSM2076293 <- read.table("GSM2076293-30199.Probe.Beta-value.txt", header=FALSE)
GSM2076294 <- read.table("GSM2076294-30200.Probe.Beta-value.txt", header=FALSE)
GSM2076295 <- read.table("GSM2076295-30201.Probe.Beta-value.txt", header=FALSE)
GSM2076296 <- read.table("GSM2076296-30202.Probe.Beta-value.txt", header=FALSE)
GSM2076297 <- read.table("GSM2076297-30203.Probe.Beta-value.txt", header=FALSE)
GSM2076298 <- read.table("GSM2076298-30204.Probe.Beta-value.txt", header=FALSE)
GSM2076299 <- read.table("GSM2076299-30205.Probe.Beta-value.txt", header=FALSE)
GSM2076300 <- read.table("GSM2076300-30206.Probe.Beta-value.txt", header=FALSE)
GSM2076301 <- read.table("GSM2076301-30207.Probe.Beta-value.txt", header=FALSE)
GSM2076302 <- read.table("GSM2076302-30208.Probe.Beta-value.txt", header=FALSE)
GSM2076303 <- read.table("GSM2076303-30209.Probe.Beta-value.txt", header=FALSE)
GSM2076304 <- read.table("GSM2076304-30210.Probe.Beta-value.txt", header=FALSE)
GSM2076305 <- read.table("GSM2076305-30211.Probe.Beta-value.txt", header=FALSE)
GSM2076306 <- read.table("GSM2076306-30212.Probe.Beta-value.txt", header=FALSE)
GSM2076307 <- read.table("GSM2076307-30213.Probe.Beta-value.txt", header=FALSE)
GSM2076308 <- read.table("GSM2076308-30214.Probe.Beta-value.txt", header=FALSE)
GSM2076309 <- read.table("GSM2076309-30215.Probe.Beta-value.txt", header=FALSE)
GSM2076310 <- read.table("GSM2076310-30216.Probe.Beta-value.txt", header=FALSE)
GSM2076311 <- read.table("GSM2076311-30217.Probe.Beta-value.txt", header=FALSE)
GSM2076312 <- read.table("GSM2076312-30218.Probe.Beta-value.txt", header=FALSE)
GSM2076313 <- read.table("GSM2076313-30219.Probe.Beta-value.txt", header=FALSE)
GSM2076314 <- read.table("GSM2076314-30220.Probe.Beta-value.txt", header=FALSE)
GSM2076315 <- read.table("GSM2076315-30221.Probe.Beta-value.txt", header=FALSE)
GSM2076316 <- read.table("GSM2076316-30222.Probe.Beta-value.txt", header=FALSE)
GSM2076317 <- read.table("GSM2076317-30223.Probe.Beta-value.txt", header=FALSE)
GSM2076318 <- read.table("GSM2076318-30224.Probe.Beta-value.txt", header=FALSE)
GSM2076319 <- read.table("GSM2076319-30225.Probe.Beta-value.txt", header=FALSE)
GSM2076320 <- read.table("GSM2076320-30226.Probe.Beta-value.txt", header=FALSE)
GSM2076321 <- read.table("GSM2076321-30227.Probe.Beta-value.txt", header=FALSE)
GSM2076322 <- read.table("GSM2076322-30228.Probe.Beta-value.txt", header=FALSE)
GSM2076323 <- read.table("GSM2076323-30229.Probe.Beta-value.txt", header=FALSE)
GSM2076324 <- read.table("GSM2076324-30230.Probe.Beta-value.txt", header=FALSE)
GSM2076325 <- read.table("GSM2076325-30231.Probe.Beta-value.txt", header=FALSE)
GSM2076326 <- read.table("GSM2076326-30232.Probe.Beta-value.txt", header=FALSE)
GSM2076327 <- read.table("GSM2076327-30233.Probe.Beta-value.txt", header=FALSE)
GSM2076328 <- read.table("GSM2076328-30234.Probe.Beta-value.txt", header=FALSE)
GSM2076329 <- read.table("GSM2076329-30235.Probe.Beta-value.txt", header=FALSE)
GSM2076330 <- read.table("GSM2076330-30236.Probe.Beta-value.txt", header=FALSE)
GSM2076331 <- read.table("GSM2076331-30237.Probe.Beta-value.txt", header=FALSE)
GSM2076332 <- read.table("GSM2076332-30238.Probe.Beta-value.txt", header=FALSE)
GSM2076333 <- read.table("GSM2076333-30239.Probe.Beta-value.txt", header=FALSE)
GSM2076334 <- read.table("GSM2076334-30240.Probe.Beta-value.txt", header=FALSE)
GSM2076335 <- read.table("GSM2076335-30241.Probe.Beta-value.txt", header=FALSE)
GSM2076336 <- read.table("GSM2076336-30242.Probe.Beta-value.txt", header=FALSE)
GSM2076337 <- read.table("GSM2076337-30243.Probe.Beta-value.txt", header=FALSE)
GSM2076338 <- read.table("GSM2076338-30244.Probe.Beta-value.txt", header=FALSE)
GSM2076339 <- read.table("GSM2076339-30245.Probe.Beta-value.txt", header=FALSE)
GSM2076340 <- read.table("GSM2076340-30246.Probe.Beta-value.txt", header=FALSE)
GSM2076341 <- read.table("GSM2076341-30247.Probe.Beta-value.txt", header=FALSE)
GSM2076342 <- read.table("GSM2076342-30248.Probe.Beta-value.txt", header=FALSE)
GSM2076343 <- read.table("GSM2076343-30249.Probe.Beta-value.txt", header=FALSE)
GSM2076344 <- read.table("GSM2076344-30250.Probe.Beta-value.txt", header=FALSE)
GSM2076345 <- read.table("GSM2076345-30251.Probe.Beta-value.txt", header=FALSE)
GSM2076346 <- read.table("GSM2076346-30252.Probe.Beta-value.txt", header=FALSE)
GSM2076347 <- read.table("GSM2076347-30253.Probe.Beta-value.txt", header=FALSE)
GSM2076348 <- read.table("GSM2076348-30254.Probe.Beta-value.txt", header=FALSE)
GSM2076349 <- read.table("GSM2076349-30255.Probe.Beta-value.txt", header=FALSE)
GSM2076350 <- read.table("GSM2076350-30256.Probe.Beta-value.txt", header=FALSE)
GSM2076351 <- read.table("GSM2076351-30257.Probe.Beta-value.txt", header=FALSE)
GSM2076352 <- read.table("GSM2076352-30258.Probe.Beta-value.txt", header=FALSE)
GSM2076353 <- read.table("GSM2076353-30259.Probe.Beta-value.txt", header=FALSE)
GSM2076354 <- read.table("GSM2076354-30260.Probe.Beta-value.txt", header=FALSE)
GSM2076355 <- read.table("GSM2076355-30261.Probe.Beta-value.txt", header=FALSE)
GSM2076356 <- read.table("GSM2076356-30262.Probe.Beta-value.txt", header=FALSE)
GSM2076357 <- read.table("GSM2076357-30263.Probe.Beta-value.txt", header=FALSE)
GSM2076358 <- read.table("GSM2076358-30264.Probe.Beta-value.txt", header=FALSE)
GSM2076359 <- read.table("GSM2076359-30265.Probe.Beta-value.txt", header=FALSE)
GSM2076360 <- read.table("GSM2076360-30266.Probe.Beta-value.txt", header=FALSE)
GSM2076361 <- read.table("GSM2076361-30267.Probe.Beta-value.txt", header=FALSE)
GSM2076362 <- read.table("GSM2076362-30268.Probe.Beta-value.txt", header=FALSE)
GSM2076363 <- read.table("GSM2076363-30269.Probe.Beta-value.txt", header=FALSE)
GSM2076364 <- read.table("GSM2076364-30270.Probe.Beta-value.txt", header=FALSE)
GSM2076365 <- read.table("GSM2076365-30271.Probe.Beta-value.txt", header=FALSE)
GSM2076366 <- read.table("GSM2076366-30272.Probe.Beta-value.txt", header=FALSE)
GSM2076367 <- read.table("GSM2076367-30273.Probe.Beta-value.txt", header=FALSE)
GSM2076368 <- read.table("GSM2076368-30274.Probe.Beta-value.txt", header=FALSE)
GSM2076369 <- read.table("GSM2076369-30275.Probe.Beta-value.txt", header=FALSE)
GSM2076370 <- read.table("GSM2076370-30276.Probe.Beta-value.txt", header=FALSE)

#Total Probes in Bens each data: 482421
GSM_Bens_Betavalues <- cbind(GSM2076276,GSM2076277,GSM2076278,GSM2076279,GSM2076280,GSM2076281,GSM2076282,GSM2076283,GSM2076284,GSM2076285,GSM2076286,GSM2076287,GSM2076288,GSM2076289,GSM2076290,GSM2076291,GSM2076292,GSM2076293,GSM2076294,GSM2076295,GSM2076296,GSM2076297,GSM2076298,GSM2076299,GSM2076300,GSM2076301,GSM2076302,GSM2076303,GSM2076304,GSM2076305,GSM2076306,GSM2076307,GSM2076308,GSM2076309,GSM2076310,GSM2076311,GSM2076312,GSM2076313,GSM2076314,GSM2076315,GSM2076316,GSM2076317,GSM2076318,GSM2076319,GSM2076320,GSM2076321,GSM2076322,GSM2076323,GSM2076324,GSM2076325,GSM2076326,GSM2076327,GSM2076328,GSM2076329,GSM2076330,GSM2076331,GSM2076332,GSM2076333,GSM2076334,GSM2076335,GSM2076336,GSM2076337,GSM2076338,GSM2076339,GSM2076340,GSM2076341,GSM2076342,GSM2076343,GSM2076344,GSM2076345,GSM2076346,GSM2076347,GSM2076348,GSM2076349,GSM2076350,GSM2076351,GSM2076352,GSM2076353,GSM2076354,GSM2076355,GSM2076356,GSM2076357,GSM2076358,GSM2076359,GSM2076360,GSM2076361,GSM2076362,GSM2076363,GSM2076364,GSM2076365,GSM2076366,GSM2076367,GSM2076368,GSM2076369,GSM2076370)
head(GSM_Bens_Betavalues)
tail(GSM_Bens_Betavalues)
colnames(GSM_Bens_Betavalues) <- c("Probe","GSM2076276","Probe","GSM2076277","Probe","GSM2076278","Probe","GSM2076279","Probe","GSM2076280","Probe","GSM2076281","Probe","GSM2076282","Probe","GSM2076283","Probe","GSM2076284","Probe","GSM2076285","Probe","GSM2076286","Probe","GSM2076287","Probe","GSM2076288","Probe","GSM2076289","Probe","GSM2076290","Probe","GSM2076291","Probe","GSM2076292","Probe","GSM2076293","Probe","GSM2076294","Probe","GSM2076295","Probe","GSM2076296","Probe","GSM2076297","Probe","GSM2076298","Probe","GSM2076299","Probe","GSM2076300","Probe","GSM2076301","Probe","GSM2076302","Probe","GSM2076303","Probe","GSM2076304","Probe","GSM2076305","Probe","GSM2076306","Probe","GSM2076307","Probe","GSM2076308","Probe","GSM2076309","Probe","GSM2076310","Probe","GSM2076311","Probe","GSM2076312","Probe","GSM2076313","Probe","GSM2076314","Probe","GSM2076315","Probe","GSM2076316","Probe","GSM2076317","Probe","GSM2076318","Probe","GSM2076319","Probe","GSM2076320","Probe","GSM2076321","Probe","GSM2076322","Probe","GSM2076323","Probe","GSM2076324","Probe","GSM2076325","Probe","GSM2076326","Probe","GSM2076327","Probe","GSM2076328","Probe","GSM2076329","Probe","GSM2076330","Probe","GSM2076331","Probe","GSM2076332","Probe","GSM2076333","Probe","GSM2076334","Probe","GSM2076335","Probe","GSM2076336","Probe","GSM2076337","Probe","GSM2076338","Probe","GSM2076339","Probe","GSM2076340","Probe","GSM2076341","Probe","GSM2076342","Probe","GSM2076343","Probe","GSM2076344","Probe","GSM2076345","Probe","GSM2076346","Probe","GSM2076347","Probe","GSM2076348","Probe","GSM2076349","Probe","GSM2076350","Probe","GSM2076351","Probe","GSM2076352","Probe","GSM2076353","Probe","GSM2076354","Probe","GSM2076355","Probe","GSM2076356","Probe","GSM2076357","Probe","GSM2076358","Probe","GSM2076359","Probe","GSM2076360","Probe","GSM2076361","Probe","GSM2076362","Probe","GSM2076363","Probe","GSM2076364","Probe","GSM2076365","Probe","GSM2076366","Probe","GSM2076367","Probe","GSM2076368","Probe","GSM2076369","Probe","GSM2076370")
head(GSM_Bens_Betavalues)
GSM_Bens_Betavalues_PCA <- GSM_Bens_Betavalues[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190)]
head(GSM_Bens_Betavalues_PCA)

############################################################## PCA Bens #####################################à

#Plot PCA by group color and labelling
mydata1=GSM_Bens_Betavalues_PCA
rownames(mydata1)
mydata1[,1]
rownames(mydata1)=mydata1[,1]
rownames(mydata1)
colnames(mydata1)
mydata1 = mydata1[,-1]
head(mydata1)
str(mydata1)
#write.table(mydata1, "GSM_Bens_Betavalues_PCA.txt", sep="\t", quote = FALSE, append = FALSE, col.names = TRUE)

#df <- as.data.frame(mydata1)
df <- mydata1
head(df)
df_pca <- prcomp(df)
df_out <- as.data.frame(df_pca$rotation)
head(df_out)
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

#Color By hroup as mentioned in file color
z=read.table("color.txt", header = FALSE, stringsAsFactors = FALSE)
z<-z[2:96]
z <- as.list(z)
z <- unlist(z)
z

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=z, label=z ))
p<-p+geom_point()+theme + xlab(percentage[1]) + ylab(percentage[2])
p <- p+theme_bw()
p
ggsave("PCA_GSMBens.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#################################################### PCA Bens Averaged #####################################
#grep "cg" GSM_Bens_Betavalues_PCA.txt > GSM_Bens_Betavalues_PCA-header.txt
probename_chrpos_Infinium <- read.table("/media/ankitv/Archivio2/ankit/BWS_array/Bens/probename_chrpos_Infinium.txt", header = TRUE)
head(probename_chrpos_Infinium)
GSM_Bens_Betavalues_PCA_chr <- cbind(probename_chrpos_Infinium, GSM_Bens_Betavalues_PCA)
GSM_Bens_Betavalues_PCA_chr_re <- GSM_Bens_Betavalues_PCA_chr[,c(2:4, 1, 6:100)]
head(GSM_Bens_Betavalues_PCA_chr_re)
dim(GSM_Bens_Betavalues_PCA_chr_re)
#write.table(GSM_Bens_Betavalues_PCA_chr_re, "GSM_Bens_Betavalues_PCA_chr_re.txt", sep="\t", quote = FALSE, append = FALSE, col.names = FALSE, row.names = F)
#bedtools intersect -wa -wb -a GSM_Bens_Betavalues_PCA_chr_re.txt -b ~/ref_av/human_ICR.bed > GSM_Bens_Betavalues_PCA.sort.chr_humanICR.txt
##Extract Bens Probes lies under human ICR: awk '{print $4}' GSM_Bens_Betavalues_PCA.sort.chr_humanICR.txt > GSM_Bens_Betavalues_PCA.sort.chr_humanICR.ProbeID
#These probes will be used for any further overlappping and extractions
setwd("/media/ankitv/Archivio2/ankit/BWS_array/Bens/Samples")
data1 <- read.table("GSM_Bens_Betavalues_PCA.sort.chr_humanICR.txt", header = FALSE)
rownames(data1)
colnames(data1)
head(data1)
data2 <- data1[,c(103, 5:99)]
colnames(data2) <- c("DMR","GSM2076276","GSM2076277","GSM2076278","GSM2076279","GSM2076280","GSM2076281","GSM2076282","GSM2076283","GSM2076284","GSM2076285","GSM2076286","GSM2076287","GSM2076288","GSM2076289","GSM2076290","GSM2076291","GSM2076292","GSM2076293","GSM2076294","GSM2076295","GSM2076296","GSM2076297","GSM2076298","GSM2076299","GSM2076300","GSM2076301","GSM2076302","GSM2076303","GSM2076304","GSM2076305","GSM2076306","GSM2076307","GSM2076308","GSM2076309","GSM2076310","GSM2076311","GSM2076312","GSM2076313","GSM2076314","GSM2076315","GSM2076316","GSM2076317","GSM2076318","GSM2076319","GSM2076320","GSM2076321","GSM2076322","GSM2076323","GSM2076324","GSM2076325","GSM2076326","GSM2076327","GSM2076328","GSM2076329","GSM2076330","GSM2076331","GSM2076332","GSM2076333","GSM2076334","GSM2076335","GSM2076336","GSM2076337","GSM2076338","GSM2076339","GSM2076340","GSM2076341","GSM2076342","GSM2076343","GSM2076344","GSM2076345","GSM2076346","GSM2076347","GSM2076348","GSM2076349","GSM2076350","GSM2076351","GSM2076352","GSM2076353","GSM2076354","GSM2076355","GSM2076356","GSM2076357","GSM2076358","GSM2076359","GSM2076360","GSM2076361","GSM2076362","GSM2076363","GSM2076364","GSM2076365","GSM2076366","GSM2076367","GSM2076368","GSM2076369","GSM2076370")
head(data2)
aggregate1 = aggregate(data2[,2:96],by=list(data2$DMR), median)
head(aggregate1, 2)
#Plot PCA by group color and labelling
dataICR1=aggregate1
rownames(dataICR1)
dataICR1[,1]
rownames(dataICR1)=dataICR1[,1]
rownames(dataICR1)
colnames(dataICR1)
dataICR1 = dataICR1[,-1]
head(dataICR1)
str(dataICR1)

#df <- as.data.frame(dataICR1)
dfICR <- dataICR1
head(dfICR)
dfICR_pca <- prcomp(dfICR)
dfICR_out <- as.data.frame(dfICR_pca$rotation)
head(dfICR_out)
percentageICR <- round(dfICR_pca$sdev / sum(dfICR_pca$sdev) * 100, 2)
percentageICR <- paste( colnames(dfICR_out), "(", paste( as.character(percentageICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

#Color By hroup as mentioned in file color
z=read.table("color.txt", header = FALSE, stringsAsFactors = FALSE)
z<-z[2:96]
z <- as.list(z)
z <- unlist(z)
z

p<-ggplot(dfICR_out,aes(x=PC1,y=PC2,color=z, label=z ))
p<-p+geom_point()+theme + xlab(percentageICR[1]) + ylab(percentageICR[2])
p <- p+theme_bw()
p
ggsave("PCA_GSMBens.humanICR.averaged.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

################################### Heatmap human ICR ######################################àà

dataICR2 = as.matrix(dataICR1)
head(dataICR2)
dim(dataICR2)
colfunc <- colorRampPalette(c("#291BEB","white", "red"))
#Control Samples
heatmap.2(dataICR2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#Let it cluster
heatmap.2(dataICR2,trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))

############################################à Data M-> Beta ####################################
#https://rdrr.io/bioc/lumi/man/m2beta.html
#Joslo
library(lumi)
setwd("/media/ankitv/Archivio2/ankit/BWS_array/")
Joslo <- read.table("JosLo_pvalue_M_annot.6164647024_R03C02BWS.CTRL50.chTest.ALL.rearranged.onlyMvalues.txt", header = TRUE)
rownames(Joslo)
Joslo[,1]
rownames(Joslo)=Joslo[,1]
rownames(Joslo)
colnames(Joslo)
Joslo = Joslo[,-1]
head(Joslo)
Joslom2beta = m2beta(Joslo)
head(Joslom2beta)
#write.table(Joslom2beta, "JosLo_pvalue_M_annot.6164647024_R03C02BWS.CTRL50.chTest.ALL.rearranged.Betavalues.txt", sep="\t", quote = FALSE, append = FALSE)

#MW
setwd("/media/ankitv/Archivio2/ankit/BWS_array/")
MW <- read.table("MW_IoW7.resultsALLanno.rearranged.onlyMvalues.txt", header = TRUE)
rownames(MW)
MW[,1]
rownames(MW)=MW[,1]
rownames(MW)
colnames(MW)
MW = MW[,-1]
head(MW)
MWm2beta = m2beta(MW)
head(MWm2beta)
dim(MWm2beta)
#write.table(MWm2beta, "MW_IoW7.resultsALLanno.rearranged.Betavalues.txt", sep="\t", quote = FALSE, append = FALSE)

############################################ Data #################################################
rearranged_Locatelli <- read.table("PCA_Locatelli_family_Controls-dec2018.chrPos_shared_Probed_overlapped_humanICR.txt", header = FALSE)
head(rearranged_Locatelli)
head(rearranged_Locatelli[,c(20,6:15)])
rearranged_Locatelli <- rearranged_Locatelli[,c(20,6:15)]
head(rearranged_Locatelli, 10)
aggregate_Locatelli = aggregate(rearranged_Locatelli[,2:11],by=list(rearranged_Locatelli$V20), median)
head(aggregate_Locatelli, 3)
# system("grep DIRAS3:ex2-DMR PCA_Locatelli_family_Controls-dec2018.chrPos_shared_Probed_overlapped_humanICR.txt")

rearranged_bens <- read.table("GSM_Bens_Betavalues_PCA.sort.chr_shared_Probe_overlapped_humanICR.txt", header = FALSE)
head(rearranged_bens)
head(rearranged_bens[,c(104,5:99)])
rearranged_bens <- rearranged_bens[,c(104,5:99)]
head(rearranged_bens, 1)
aggregate_bens = aggregate(rearranged_bens[,2:96],by=list(rearranged_bens$V104), median)
head(aggregate_bens, 1)

rearranged_Joslo <- read.table("JosLo_pvalue_M_annot.6164647024_R03C02BWS.CTRL50.chTest.ALL.rearranged.Betavalues.sorted_overlapped_humanICR.txt", header = FALSE)
head(rearranged_Joslo)
head(rearranged_Joslo[,c(12,4:5)])
rearranged_Joslo <- rearranged_Joslo[,c(12,4:5)]
head(rearranged_Joslo, 1)
aggregate_Joslo = aggregate(rearranged_Joslo[,2:3],by=list(rearranged_Joslo$V12), median)
head(aggregate_Joslo, 1)

rearranged_MW <- read.table("MW_IoW7.resultsALLanno.rearranged.Betavalues.sorted_overlapped_humanICR.txt", header = FALSE)
head(rearranged_MW)
head(rearranged_MW[,c(12,4:5)])
rearranged_MW <- rearranged_MW[,c(12,4:5)]
head(rearranged_MW, 1)
aggregate_MW = aggregate(rearranged_MW[,2:3],by=list(rearranged_MW$V12), median)
head(aggregate_MW, 1)

Merged_LBJM <- cbind(aggregate_Locatelli, aggregate_bens, aggregate_Joslo, aggregate_MW)
head(Merged_LBJM )
#write.table(Merged_LBJM , "Merged_LBJM.txt", sep="\t", quote = FALSE, append = FALSE)

#################################################### Merged data #####################################à
setwd("/media/ankitv/Archivio2/ankit/BWS_array/")
MergedLBJM21 <- read.table("MMergedLBJM1dedup.txt", header = TRUE, stringsAsFactors = FALSE)
rownames(MergedLBJM21)
head(MergedLBJM21)
dataMergedLBJM1=MergedLBJM21
dataMergedLBJM1[,1]
rownames(dataMergedLBJM1)=dataMergedLBJM1[,1]
rownames(dataMergedLBJM1)
colnames(dataMergedLBJM1)
dataMergedLBJM1 = dataMergedLBJM1[,-1]
head(dataMergedLBJM1)

#df <- as.data.frame(dataMergedLBJM1)
dfICRLBJM1 <- dataMergedLBJM1
head(dfICRLBJM1)
colnames(dfICRLBJM1)
dfICRLBJM1_pca <- prcomp(dfICRLBJM1)
dfICRLBJM1_out <- as.data.frame(dfICRLBJM1_pca$rotation, Color=z)
head(dfICRLBJM1_out)
percentageICR <- round(dfICRLBJM1_pca$sdev^2 / sum(dfICRLBJM1_pca$sdev^2) * 100, 2)
percentageICR <- paste( colnames(dfICRLBJM1_out), "(", paste( as.character(percentageICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

#Color By hroup as mentioned in file color
z=read.table("color.txt", header = FALSE, stringsAsFactors = FALSE)
z<-z[2:91]
z <- as.list(z)
z <- unlist(z)
z

p<-ggplot(dfICRLBJM1_out,aes(x=PC1,y=PC2,col=z))
p<-p+geom_point()+theme + xlab(percentageICR[1]) + ylab(percentageICR[2])
p <- p+theme_bw()
p
ggsave("PCA-MergedLBJM.humanICR.averageddedup.svg", width=25*1.25, height=18*1.25, units="cm", dpi=96)

################################### Heatmap human ICR ######################################àà
setwd("/media/ankitv/Archivio2/ankit/BWS_array/")
library(ggplot2)
library(gplots)
dataMergedLBJM2 <- read.table("dataMergedLBJM2dedup.reorderICR.txt", header = TRUE)
head(dataMergedLBJM2,1)
rownames(dataMergedLBJM2)
dataMergedLBJM2[,1]
rownames(dataMergedLBJM2)=dataMergedLBJM2[,1]
rownames(dataMergedLBJM2)
colnames(dataMergedLBJM2)
dataMergedLBJM2 = dataMergedLBJM2[,-1]
head(dataMergedLBJM2)
colnames(dataMergedLBJM2)
dataMergedLBJM2 = as.matrix(dataMergedLBJM2)
dim(dataMergedLBJM2)
head(dataMergedLBJM2,1)
colfunc <- colorRampPalette(c("#291BEB","white", "red"))
#Col cluster
svg(filename="Heatmap_dataMergedLBJM2dedup.reorderICR_heamap.svg", width=10, height=10, pointsize=12)
heatmap.2(dataMergedLBJM2, Rowv = "NA",trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="column", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(0,1, length.out = 100))
dev.off()
#Heatmap rearranged1
dataMergedLBJM2 <- read.table("dataMergedLBJM2dedup.reorderICR-selected.txt", header = TRUE)
head(dataMergedLBJM2,1)
rownames(dataMergedLBJM2)
dataMergedLBJM2[,1]
rownames(dataMergedLBJM2)=dataMergedLBJM2[,1]
rownames(dataMergedLBJM2)
colnames(dataMergedLBJM2)
dataMergedLBJM2 = dataMergedLBJM2[,-1]
head(dataMergedLBJM2)
colnames(dataMergedLBJM2)
dataMergedLBJM2 = as.matrix(dataMergedLBJM2)
dim(dataMergedLBJM2)
head(dataMergedLBJM2,1)
colfunc <- colorRampPalette(c("#291BEB","white", "red"))
#No cluster
heatmap.2(dataMergedLBJM2, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 0.7, margins =c(6,5), breaks = seq(0,1, length.out = 100))
heatmap.2(dataMergedLBJM2, trace = "none", lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2,  density.info=c("none"), cexRow=0.8, font=3, cexCol = 0.7, margins =c(6,5), breaks = seq(0,1, length.out = 100))


#Heatmap rearranged2
dataMergedLBJM2 <- read.table("dataMergedLBJM2dedup.reorderICR-selected_heamap.txt", header = TRUE)
head(dataMergedLBJM2,1)
rownames(dataMergedLBJM2)
dataMergedLBJM2[,1]
rownames(dataMergedLBJM2)=dataMergedLBJM2[,1]
rownames(dataMergedLBJM2)
colnames(dataMergedLBJM2)
dataMergedLBJM2 = dataMergedLBJM2[,-1]
head(dataMergedLBJM2)
colnames(dataMergedLBJM2)
dataMergedLBJM2 = as.matrix(dataMergedLBJM2)
dim(dataMergedLBJM2)
head(dataMergedLBJM2,1)
colfunc <- colorRampPalette(c("#291BEB","white", "red"))
#Col cluster
svg(filename="Heatmap_dataMergedLBJM2dedup.reorderICR-selected_heamap.svg", width=5, height=10, pointsize=12)
heatmap.2(dataMergedLBJM2, Rowv = "NA",trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="column", cexRow=0.8, font=3, cexCol = 0.7, margins =c(8,6), breaks = seq(0,1, length.out = 100))
dev.off()

#Scripts copied from transform_BWS-bens.R and merged in this for easy use
setwd("/media/ankitv/Archivio2/ankit/BWS_array/")
#################################################### Transform data #####################################à

setwd("/media/ankitv/Archivio2/ankit/BWS_array/")
MergedLBJM21 <- read.table("Merged_LBJM_dedup.txt", header = FALSE)
rownames(MergedLBJM21)
colnames(MergedLBJM21)
head(MergedLBJM21)
tMergedLBJM21 = t(MergedLBJM21)
head(tMergedLBJM21)
tMergedLBJM21 = data.frame(tMergedLBJM21)
tMergedLBJM21 = tMergedLBJM21[,c(2:50,1)]
head(tMergedLBJM21)
#write.table(dataMergedLBJM1 , "tMergedLBJM1dedup.txt", sep="\t", quote = FALSE, append = FALSE)

#Select only MLID cases using excel sheet

TMergedLBJM21 <- read.table("tMergedLBJM1dedup.txt", header = TRUE)
dataMergedLBJM1=TMergedLBJM21
rownames(dataMergedLBJM1)
dataMergedLBJM1[,1]
rownames(dataMergedLBJM1)=dataMergedLBJM1[,1]
rownames(dataMergedLBJM1)
colnames(dataMergedLBJM1)
dataMergedLBJM1 = dataMergedLBJM1[,-1]
head(dataMergedLBJM1)

dfx <-dataMergedLBJM1[c(1:48)]
head(dfx)
PC<-prcomp(dfx)
PCi<-data.frame(PC$x,Color=dataMergedLBJM1$Color)
percentageICR <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2)
percentageICR <- paste( colnames(PCi), "(", paste( as.character(percentageICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))


#df <- as.data.frame(dataMergedLBJM1)
p<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentageICR[1]) + ylab(percentageICR[2])+
  geom_point(size=2,alpha=1)+
  scale_color_manual(values = c("#D3D3D3","#808080","#696969","Blue","Magenta","#FA8072","Red","Orange","#00FFFF","#8B0000","#CC99FF","#666600","#99F38F","#2CC71B"))
p <- p+theme_bw()
p
ggsave("PCA-TMergedLBJM.humanICR.averageddedup.svg", width=15*1.25, height=10*1.25, units="cm", dpi=96)
#Checking labels on PCA
#plot(PCi$PC1, PCi$PC2)
#text(PCi$PC1, PCi$PC2, rownames(dataMergedLBJM1) , pos = 2, col = "black")
################################### Heatmap human ICR ######################################àà

library(ggplot2)
library(gplots)
dataMergedLBJM2 <- read.table("dataMergedLBJM2dedup.txt", header = TRUE)
head(dataMergedLBJM2)
rownames(dataMergedLBJM2)
dataMergedLBJM2[,1]
rownames(dataMergedLBJM2)=dataMergedLBJM2[,1]
rownames(dataMergedLBJM2)
colnames(dataMergedLBJM2)
dataMergedLBJM2 = dataMergedLBJM2[,-1]
head(dataMergedLBJM2)
colnames(dataMergedLBJM2)
dataMergedLBJM2 = as.matrix(dataMergedLBJM2)
dim(dataMergedLBJM2)
head(dataMergedLBJM2,1)
colfunc <- colorRampPalette(c("#291BEB","white", "red"))
#Control Samples
heatmap.2(dataMergedLBJM2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#Let it cluster
heatmap.2(dataMergedLBJM2,trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(dataMergedLBJM2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))


