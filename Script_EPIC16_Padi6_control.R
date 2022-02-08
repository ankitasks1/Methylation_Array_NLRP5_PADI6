if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
#install via terminal
git clone https://github.com/YuanTian1991/ChAMP.git
sudo R CMD INSTALL ChAMP
library(ChAMP)

setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control")
dir()
#Note data read by champ.import() can not be used for SWAN normalization and FunctionNormalization in champ.norm() function. If user want to use SWAN, you may still consider champ.load() function, but remember to set "method" parameter as "ChAMP", which is "ChAMP" in default.
myImport <- champ.import(getwd(), arraytype="EPIC")
dim(myImport$beta) 
head(myImport$beta) 
myImport2 <- data.frame(myImport$beta[order(rownames(myImport$beta)),])
head(myImport2,2)
dim(myImport2)
#Get detection p-value
mydetP <- myImport$detP
head(mydetP)
mydetP2 <- data.frame(mydetP[order(rownames(mydetP)),])
colnames(mydetP2) <- c("1-FM_detPvalue","2-FS_detPvalue","3-FM_detPvalue","4-AM_detPvalue","5-AA_detPvalue","6-AM_detPvalue","7-PM_detPvalue","8-PA_detPvalue","9-PG_detPvalue","10-PM_detPvalue","CTRL01_detPvalue","CTRL02_detPvalue","CTRL03_detPvalue","CTRL04_detPvalue","CTRL05_detPvalue","CTRL06_detPvalue","11-Ctrl1_detPvalue","12-Ctrl2_detPvalue","13-Ctrl3_detPvalue","14-Ctrl4_detPvalue","15-Ctrl5_detPvalue","16-Ctrl6_detPvalue")
head(mydetP2)
dim(mydetP2)

#Get intensity unmeth and meth
myIUnmeth <- myImport$UnMeth
head(myIUnmeth)
myIUnmeth2 <- data.frame(myIUnmeth[order(rownames(myIUnmeth)),])
colnames(myIUnmeth2) <- c("1-FM_Unmethylated_Signal","2-FS_Unmethylated_Signal","3-FM_Unmethylated_Signal","4-AM_Unmethylated_Signal","5-AA_Unmethylated_Signal","6-AM_Unmethylated_Signal","7-PM_Unmethylated_Signal","8-PA_Unmethylated_Signal","9-PG_Unmethylated_Signal","10-PM_Unmethylated_Signal","CTRL01_Unmethylated_Signal","CTRL02_Unmethylated_Signal","CTRL03_Unmethylated_Signal","CTRL04_Unmethylated_Signal","CTRL05_Unmethylated_Signal","CTRL06_Unmethylated_Signal","11-Ctrl1_Unmethylated_Signal","12-Ctrl2_Unmethylated_Signal","13-Ctrl3_Unmethylated_Signal","14-Ctrl4_Unmethylated_Signal","15-Ctrl5_Unmethylated_Signal","16-Ctrl6_Unmethylated_Signal")
head(myIUnmeth2)
dim(myIUnmeth2)


myIMeth <- myImport$Meth
head(myIMeth)
myIMeth2 <- data.frame(myIMeth[order(rownames(myIMeth)),])
colnames(myIMeth2) <- c("1-FM_Methylated_Signal","2-FS_Methylated_Signal","3-FM_Methylated_Signal","4-AM_Methylated_Signal","5-AA_Methylated_Signal","6-AM_Methylated_Signal","7-PM_Methylated_Signal","8-PA_Methylated_Signal","9-PG_Methylated_Signal","10-PM_Methylated_Signal","CTRL01_Methylated_Signal","CTRL02_Methylated_Signal","CTRL03_Methylated_Signal","CTRL04_Methylated_Signal","CTRL05_Methylated_Signal","CTRL06_Methylated_Signal","11-Ctrl1_Methylated_Signal","12-Ctrl2_Methylated_Signal","13-Ctrl3_Methylated_Signal","14-Ctrl4_Methylated_Signal","15-Ctrl5_Methylated_Signal","16-Ctrl6_Methylated_Signal")
head(myIMeth2)
dim(myIMeth2)

Matrix_signal_intensities <- cbind.data.frame(myIUnmeth2, myIMeth2, mydetP2)
head(Matrix_signal_intensities,1)
dim(Matrix_signal_intensities)
Matrix_signal_intensitiesre <- Matrix_signal_intensities[,c(1,23,45,2,24,46,3,25,47,4,26,48,5,27,49,6,28,50,7,29,51,8,30,52,9,31,53,10,32,54,11,33,55,12,34,56,13,35,57,14,36,58,15,37,59,16,38,60,17,39,61,18,40,62,19,41,63,20,42,64,21,43,65,22,44,66)]
dim(Matrix_signal_intensitiesre)
head(Matrix_signal_intensitiesre,1)
write.table(Matrix_signal_intensitiesre,"Matrix_signal_intensitiesre.txt",row.names=TRUE,quote=FALSE, append=F,sep = "\t")
#GEO Submission
system("cp Matrix_signal_intensitiesre.txt ./../GEOSubmission/Matrix_signal_intensities.txt")
#edit manually to add GEO header
#write.table(myImport$beta,"myImportbeta.txt",row.names=TRUE,quote=FALSE)
#myfilter <- champ.filter(beta=myImport$beta,M=NULL,pd=myImport$pd,intensity=NULL,Meth=NULL,UnMeth=NULL,detP=myImport$detP,beadcount=myImport$beadcount,autoimpute=TRUE,filterDetP=TRUE,ProbeCutoff=0,SampleCutoff=0.1,detPcut=0.01,filterBeads=TRUE,beadCutoff=0.05,filterNoCG = TRUE,filterSNPs = TRUE,population = NULL,filterMultiHit = TRUE,filterXY = TRUE,fixOutlier = TRUE,arraytype = "EPIC")
#head(myfilter$beta)
#myLoad <- champ.filter(beta=myImport$beta, arraytype = "EPIC", filterDetP=TRUE, detPcut=0.01)
#data(EPIC.manifest.hg38)
#LOAD FILE (method=ChAMP) PER SWAN NORMALIZATION###########
myLoad_2<-champ.load(directory = getwd(),
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     filterBeads=TRUE,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=FALSE,
                     arraytype="EPIC")
#write.table(myLoad_2$beta,"myLoad_2_beta.txt",row.names=TRUE,quote=FALSE)
summary(myLoad_2)
head(myLoad_2$beta)
head(myLoad_2$pd)
dim(myLoad_2$beta)
myLoad_3 <- data.frame(myLoad_2$beta[order(rownames(myLoad_2$beta)),])
head(myLoad_3,2)
dim(myLoad_3)
#QC
#QC.GUI(myLoad_2$beta,
#       pheno=myLoad_2$pd$Sample_Group, 
#       arraytype="EPIC")
head(myLoad_2$beta)
DATA_myLoad2 <- cbind.data.frame(rownames(myLoad_2$beta), myLoad_2$beta)
head(DATA_myLoad2)
colnames(DATA_myLoad2) <- c("TargetID","1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6")
head(DATA_myLoad2)
dim(DATA_myLoad2)
DATA_myLoad2sort <- DATA_myLoad2[order(DATA_myLoad2$TargetID),]
head(DATA_myLoad2sort)
head(mydetP2)
DATA_mydetP2 <- cbind.data.frame(rownames(mydetP2), mydetP2)
head(DATA_mydetP2)
colnames(DATA_mydetP2) <- c("TargetID","1-FM_detPvalue","2-FS_detPvalue","3-FM_detPvalue","4-AM_detPvalue","5-AA_detPvalue","6-AM_detPvalue","7-PM_detPvalue","8-PA_detPvalue","9-PG_detPvalue","10-PM_detPvalue","CTRL01_detPvalue","CTRL02_detPvalue","CTRL03_detPvalue","CTRL04_detPvalue","CTRL05_detPvalue","CTRL06_detPvalue","11-Ctrl1_detPvalue","12-Ctrl2_detPvalue","13-Ctrl3_detPvalue","14-Ctrl4_detPvalue","15-Ctrl5_detPvalue","16-Ctrl6_detPvalue")
head(DATA_mydetP2)
dim(DATA_mydetP2)
DATA_mydetP2sort <- DATA_mydetP2[order(DATA_mydetP2$TargetID),]
head(DATA_mydetP2sort)
MERGE_LoaddetP <- merge(DATA_myLoad2sort,DATA_mydetP2sort,by.x="TargetID",by.y="TargetID")
#check for dimensions
dim(DATA_mydetP2)
dim(MERGE_LoaddetP)
head(MERGE_LoaddetP,2)
MERGE_LoaddetPre <- MERGE_LoaddetP[,c(1,2,24,3,25,4,26,5,27,6,28,7,29,8,30,9,31,10,32,11,33,12,34,13,35,14,36,15,37,16,38,17,39,18,40,19,41,20,42,21,43,22,44,23,45)]
head(MERGE_LoaddetPre,1)
#write files
write.table(MERGE_LoaddetPre,"MERGE_LoaddetPre.txt",row.names=FALSE,quote=FALSE,append = F, sep = "\t")
#GEO submission
system("cp MERGE_LoaddetPre.txt ./../GEOSubmission/Matrix_processed.txt")

# NORMALIZE DATASET (ChAMP load)
#myNorm2<- champ.norm(myLoad_2$beta,
#                    method="BMIQ",
#                    arraytype="EPIC")
dim(myNorm2)
head(myNorm2)
write.table(myNorm2,"myNorm2.txt",row.names=TRUE,quote=FALSE,append = F,sep = "\t")

myNorm2 <- read.table("myNorm2.txt", header = TRUE)
myNorm3 <- data.frame(myNorm2[order(rownames(myNorm2)),])
head(myNorm3)
write.table(myNorm3,"myNorm3.txt",row.names=TRUE,quote=FALSE,append = F,sep = "\t")

#other features not used
#QUALITY CONTROL-TIME CONSUMING BE-CARE
#QC.GUI(myNorm2,
#       pheno=myLoad_2$pd$Sample_Group,
#       arraytype="EPIC")
# Do SVD check on data set
champ.SVD(myNorm2,
          rgSet=NULL,
          pd=myLoad_2$pd)

#save as Batch_effects_SVD.svg
#######################################################
##################### No Batch correction #############
#######################################################

#write.table(myCombat2,"myCombat_2.txt",row.names=TRUE,quote=FALSE)
# If Batch detected, run champ.runCombat() here.
#MERGE myNorm2 to annotations of Illumina Manifest
#before Add TargetID to first header of myNorm2.txt
#reload modified file
#DATA_myNorm_2<-read.delim("myNorm2.txt", sep=" ")
DATA_myNorm_2 <- cbind.data.frame(rownames(myNorm2), myNorm2)
head(DATA_myNorm_2)
colnames(DATA_myNorm_2) <- c("TargetID","1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6")
head(DATA_myNorm_2)
dim(DATA_myNorm_2)
#LOAD MANIFEST hg19 B4 https://bioconductor.org/packages/release/data/annotation/manuals/IlluminaHumanMethylationEPICanno.ilm10b4.hg19/man/IlluminaHumanMethylationEPICanno.ilm10b4.hg19.pdf. Rearranged to csv
MANIFEST <- read.csv("/media/ankitv/Archivio1/testing_array/Loca/test/MethylationEPIC_v-1-0_B4_ext.csv")
#MANIFEST <- EPIC.manifest.hg19
colnames(MANIFEST)
dim(MANIFEST)
head(MANIFEST)
MERGE_myNorm_3 <- merge(DATA_myNorm_2,MANIFEST,by.x="TargetID",by.y="IlmnID")
#check for dimensions
dim(DATA_myNorm_2)
dim(MERGE_myNorm_3)
head(MERGE_myNorm_3)
#write files
write.table(MERGE_myNorm_3,"Merge_myNorm_3x.txt",row.names=FALSE,quote=FALSE)
head(myNorm2)
tmyNorm2 = t(myNorm2)
#head(tmyNorm2)
dim(tmyNorm2)
#rownames(tmyNorm2)
tmyNorm2 = data.frame(tmyNorm2)
#head(tmyNorm2)
#write.table(tmyNorm2 , "tmyNorm2dedup.txt", sep="\t", quote = FALSE, append = FALSE)
#tmyNorm2["Color"] <-  c("cNormalEpic16","bBWSEpic16","bBWSEpic16","cNormalEpic16","bBWSEpic16","cNormalEpic16","cNormalEpic16","bBWSEpic16","cNormalEpic16","cNormalEpic16","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16")
tmyNorm2["Color"] <-  c("a1-FM","a2-FS","a3-FM","a4-AM","a5-AA","a6-AM","a7-PM","a8-PA","a9-PG","b10-PM","cCTRL01","cCTRL02","cCTRL03","cCTRL04","cCTRL05","cCTRL06","cCtrl1","cCtrl2","cCtrl3","cCtrl4","cCtrl5","cCtrl6")
dim(tmyNorm2)
dfx <-tmyNorm2[c(1:736048)]
#head(dfx)
PC<-prcomp(dfx, center = TRUE, scale. = TRUE)
#head(PC)
PCi<-data.frame(PC$x,Color=tmyNorm2$Color)
percentageAll <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentageAll <- paste( colnames(PCi), "(", paste( as.character(percentageAll), "%", ")", sep="") )
library(ggplot2)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p <- ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentageAll[1]) + ylab(percentageAll[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("pink","navy","navy","pink","navy","skyblue","pink","navy","skyblue","skyblue","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  scale_shape_manual(values=c(1,1,1,2,2,2,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0))
p <- p+theme_classic()
#p + xlim(-50,50)+ ylim(-50,50)
p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA-tmyNorm2.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96)
svg(filename="boxplot_myNorm2.svg", width=10, height=5, pointsize=12)
boxplot(myNorm2, main="Normalized beta values", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue"))
dev.off()

#myDMP <- champ.DMP()
myDMP2 <- champ.DMP(beta = myNorm2,
                    pheno=myLoad_2$pd$Sample_Group, 
                    compare.group = NULL, 
                    adjPVal = 0.05, 
                    adjust.method = "BH", 
                    arraytype = "EPIC")
head(myDMP2[[1]])
DMP.GUI(DMP=myDMP2[[1]],
        beta=myNorm2,
        pheno=myLoad_2$pd$Sample_Group) #Must press submit
#myDMR <- champ.DMR()
myDMR2 <- champ.DMR(beta=myNorm2,
                    pheno=myLoad_2$pd$Sample_Group,
                    compare.group=NULL,
                    arraytype="EPIC",
                    method="Bumphunter")
head(myDMR2$BumphunterDMR) #head(myDMR2$METHODNAME)
DMR.GUI(DMR=myDMR2,
        beta=myNorm2,
        pheno=myLoad_2$pd$Sample_Group,
        compare.group= c("A","B"),
        arraytype="EPIC") #Must press submit
myBlock2 <- champ.Block(beta=myNorm2,
                        pheno=myLoad_2$pd$Sample_Group,
                        arraytype="EPIC")
head(myBlock2$Block)
Block.GUI(Block=myBlock2,
          beta=myNorm2,
          pheno=myLoad_2$pd$Sample_Group,
          runDMP=TRUE,
          compare.group=c("A","B"),
          arraytype="EPIC")
myGSEA2 <- champ.GSEA(beta=myNorm2,
                      DMP=myDMP2[[1]], 
                      DMR=myDMR2, 
                      arraytype="EPIC",
                      adjPval=0.05, 
                      method="fisher")
# myDMP and myDMR could (not must) be used directly.
head(myGSEA2$DMP)
head(myGSEA2$DMR)
#Empirical Bayes GSEA method
myebayGSEA2 <- champ.ebGSEA(beta=myNorm2,
                            pheno=myLoad_2$pd$Sample_Group,
                            arraytype="EPIC")
myEpiMod2 <- champ.EpiMod(beta=myNorm2,
                          pheno=myLoad_2$pd$Sample_Group, 
                          arraytype="EPIC")
myCNA2 <- champ.CNA(intensity=myLoad_2$intensity,
                    pheno=myLoad_2$pd$Sample_Group, 
                    arraytype="EPIC")
# If DataSet is Blood samples, run champ.refbase() here.
myRefBase <- champ.refbase(beta=myNorm2,
                           arraytype="EPIC")
library(ggpubr)
library(ggplot2)
library(plyr)


##################################################  Batch correction #################################

##############################################
#########  With Batch Correction #############
#Batch Effect detected#
#Remove batch effect using Combat#
##################
##################
sample_info <- read.table("sample_info.txt", header = TRUE)
head(sample_info)
sample_info[,1]
rownames(sample_info)=sample_info[,1]
sample_info = sample_info[,-1]
head(sample_info)

phenobatch <- read.table("phenobatch.txt", header = TRUE)
rownames(phenobatch)
phenobatch[,1]
rownames(phenobatch)=phenobatch[,1]
rownames(phenobatch)
colnames(phenobatch)
phenobatch = phenobatch[,-1]
head(phenobatch)
colnames(phenobatch)
dim(phenobatch)
#Lab batch correction
Batch = data.frame(phenobatch$Batch)
head(Batch)
rownames(Batch) <- colnames(myNorm2)
head(Batch)
Modcombat = model.matrix(~1, data=Batch)
library(sva)
#Logit transformation is required for beta values as explained in ChAMP manual to keep interval [0,1]
#See also this paper https://www.tandfonline.com/doi/full/10.1080/15592294.2018.1530008 M=logit2(beta)
myNorm2logit = ComBat(dat=as.matrix(logit2(myNorm2)), batch=as.numeric(phenobatch$Batch), mod=Modcombat, par.prior=TRUE, prior.plots=FALSE)
#Reverse logit2 transformation
myNorm2batch <- ilogit2(myNorm2logit)
champ.SVD(myNorm2batch,rgSet=NULL,
          pd=myLoad_2$pd)
summary(myNorm2batch)
#Sex variable correction
Sex = data.frame(phenobatch$Sex)
head(Sex)
rownames(Sex) <- colnames(myNorm2)
head(Sex)
Sexcombat = model.matrix(~1, data=Sex)
#BiocManager::install("sva")
library(sva)

myNorm2batchlogit = ComBat(dat=as.matrix(logit2(myNorm2batch)), batch=as.numeric(phenobatch$Sex), mod=Sexcombat, par.prior=TRUE, prior.plots=FALSE)
myNorm2batchx <- ilogit2(myNorm2batchlogit)

head(myNorm2batchx)
champ.SVD(myNorm2batchx,rgSet=NULL,
          pd=myLoad_2$pd)
summary(myNorm2batchx)
#Condition variable correction (=Sample Group):Do not remove as it is condition specific
#Condition = data.frame(phenobatch$Condition)
#head(Condition)
#rownames(Condition) <- colnames(myNorm2)
#head(Condition)
#Conditioncombat = model.matrix(~1, data=Condition)
#BiocManager::install("sva")
#myNorm2batchxc = ComBat(dat=as.matrix(myNorm2batchx), batch=as.numeric(phenobatch$Condition), mod=Conditioncombat, par.prior=TRUE, prior.plots=FALSE)
#champ.SVD(myNorm2batchxc,rgSet=NULL,
#          pd=myLoad_2$pd)
#save as SVD.myNorm2batchxc.svg
##############################################
#Correct batch effect on "Slide" factor.
#myCombat<-champ.runCombat(myNorm2,
#                          myLoad_2$pd,
#                          variablename="Sample_Group",
#                          batchname=c("Slide"),
#                          logitTrans=TRUE)


myCombat2 <- myNorm2batchx
head(myCombat2)
dim(myCombat2)

# Do SVD check on data set
champ.SVD(myCombat2,rgSet=NULL,
          pd=myLoad_2$pd)
#save as Batch_effect_corrected_SVD.svg
head(myCombat2)
myCombat2sort <- myCombat2[order(rownames(myCombat2)),]
head(myCombat2sort)
dim(myCombat2sort)
write.table(myCombat2sort,"myCombat_2.txt",row.names=TRUE,quote=FALSE, append = F, sep = "\t")
#GEO Submission
system("cp myCombat_2.txt ./../GEOSubmission/Matrix_normalized_batchcorrected.txt")

#QC.GUI(myCombat2,
#       pheno=myLoad_2$pd$Sample_Group,
#       arraytype="EPIC")
head(myCombat2)
tmyCombat2 = t(myCombat2)
#head(tmyCombat2)
dim(tmyCombat2)
rownames(tmyCombat2)
tmyCombat2 = data.frame(tmyCombat2)
#head(tmyCombat2)
#write.table(tmyCombat2 , "tmyCombat2dedup.txt", sep="\t", quote = FALSE, append = FALSE)
#tmyCombat2["Color"] <-  c("cNormalEpic16","bBWSEpic16","bBWSEpic16","cNormalEpic16","bBWSEpic16","cNormalEpic16","cNormalEpic16","bBWSEpic16","cNormalEpic16","cNormalEpic16","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16")
tmyCombat2["Color"] <- c("a1-FM","a2-FS","a3-FM","a4-AM","a5-AA","a6-AM","a7-PM","a8-PA","a9-PG","b10-PM","cCTRL01","cCTRL02","cCTRL03","cCTRL04","cCTRL05","cCTRL06","cCtrl1","cCtrl2","cCtrl3","cCtrl4","cCtrl5","cCtrl6")
library(ggplot2)
dim(tmyCombat2)
cobdfx <-tmyCombat2[c(1:736048)]
#head(dfx)
PCcob<-prcomp(cobdfx, center = TRUE, scale. = TRUE)
#head(PC)
PCcobi<-data.frame(PCcob$x,Color=tmyCombat2$Color)
percentagecobAll <- round(PCcob$sdev^2 / sum(PCcob$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecobAll <- paste( colnames(PCcobi), "(", paste( as.character(percentagecobAll), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pcob1<-ggplot(PCcobi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecobAll[1]) + ylab(percentagecobAll[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("pink","navy","navy","pink","navy","skyblue","pink","navy","skyblue","skyblue","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  scale_shape_manual(values=c(1,1,1,2,2,2,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0))
pcob1 <- pcob1+theme_classic()
#pcob1 + xlim(-40,30)+ ylim(-30,30)
pcob1 + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA-tmyCombat2.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96)


head(myCombat2)
tmyCombat3 = t(myCombat2)
#head(tmyCombat3)
dim(tmyCombat3)
#rownames(tmyCombat3)
tmyCombat3 = data.frame(tmyCombat3)
#head(tmyCombat3)
#write.table(tmyCombat3 , "tmyCombat3dedup.txt", sep="\t", quote = FALSE, append = FALSE)

tmyCombat3["Color"] <-  c("F1_Madre","F1_BWS1","F1_BWS2","F2_Madre","F2_BWS","F2_Sib","F3_Madre","F3_BWS","F3_Sib1","F3_Sib2","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control")
#head(tmyCombat3)
dim(tmyCombat3)
dfxtmyCombat3 <-tmyCombat3[c(1:736048)]
#head(dfx)
PCtmyCombat3 <-prcomp(dfxtmyCombat3, center = TRUE, scale. = TRUE)
#head(PC)
PCtmyCombat3i<-data.frame(PCtmyCombat3$x,Color=tmyCombat3$Color)
percentagetmyCombat3All <- round(PCtmyCombat3$sdev^2 / sum(PCtmyCombat3$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagetmyCombat3All <- paste( colnames(PCtmyCombat3i), "(", paste( as.character(percentagetmyCombat3All), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p33<-ggplot(PCtmyCombat3i,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagetmyCombat3All[1]) + ylab(percentagetmyCombat3All[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("darkgrey","navy","navy","magenta","navy","magenta","skyblue","navy","magenta","skyblue","skyblue"))+
  scale_shape_manual(values=c(9,12,17,17,6,6,6,19,19,15,19))
p33 <- p33 +theme_classic()
p33 + xlim(-50,50)+ ylim(-50,50)
p33 + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA-tmyCombat3.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96)
ggsave("PCA-tmyCombat3.png", width=17*1.25, height=12*1.25, units="cm", dpi=96)
mycomDMP <- champ.DMP()
mycomDMP2 <- champ.DMP(beta = myCombat2,
                       pheno=myLoad_2$pd$Sample_Group, 
                       compare.group = NULL, 
                       adjPVal = 0.05, 
                       adjust.method = "BH", 
                       arraytype = "EPIC")
head(mycomDMP2[[1]])
DMP.GUI(DMP=mycomDMP2[[1]],
        beta=myCombat2,
        pheno=myLoad_2$pd$Sample_Group) #Must press submit
#mycomDMR <- champ.DMR()
mycomDMR2 <- champ.DMR(beta=myCombat2,
                       pheno=myLoad_2$pd$Sample_Group,
                       compare.group=NULL,
                       arraytype="EPIC",
                       method="Bumphunter")
head(mycomDMR2$BumphunterDMR) #head(mycomDMR2$METHODNAME)
DMR.GUI(DMR=mycomDMR2,
        beta=myCombat2,
        pheno=myLoad_2$pd$Sample_Group,
        compare.group= c("Case","Control"),
        arraytype="EPIC") #Must press submit
mycomBlock2 <- champ.Block(beta=myCombat2,
                           pheno=myLoad_2$pd$Sample_Group,
                           arraytype="EPIC")
head(mycomBlock2$Block)
Block.GUI(Block=mycomBlock2,
          beta=myCombat2,
          pheno=myLoad_2$pd$Sample_Group,
          runDMP=TRUE,
          compare.group=c("A","B"),
          arraytype="EPIC")
mycomGSEA2 <- champ.GSEA(beta=myCombat2,
                         DMP=mycomDMP2[[1]], 
                         DMR=mycomDMR2, 
                         arraytype="EPIC",
                         adjPval=0.05, 
                         method="fisher")
# mycomDMP and mycomDMR could (not must) be used directly.
head(mycomGSEA2$DMP)
head(mycomGSEA2$DMR)
#Empirical Bayes GSEA method
myebayGSEA2 <- champ.ebGSEA(beta=myCombat2,
                            pheno=myLoad_2$pd$Sample_Group,
                            arraytype="EPIC")
mycomEpiMod2 <- champ.EpiMod(beta=myCombat2,
                             pheno=myLoad_2$pd$Sample_Group, 
                             arraytype="EPIC")
mycomCNA2 <- champ.CNA(intensity=myLoad_2$intensity,
                       pheno=myLoad_2$pd$Sample_Group, 
                       arraytype="EPIC")
# If DataSet is Blood samples, run champ.refbase() here.
mycomRefBase <- champ.refbase(beta=myCombat2,
                              arraytype="EPIC")


head(myCombat2)
#write.table(myCombat2,"myCombat2.txt",row.names=TRUE,quote=FALSE)
myCombat3 <- data.frame(myCombat2[order(rownames(myCombat2)),])
head(myCombat3)
#Take median of controls
head(myCombat3)
dim(myCombat3)
myCombat3re <- as.matrix(myCombat3)
head(myCombat3re)
dim(myCombat3re)
svg(filename="boxplot_myCombat3re.svg", width=10, height=5, pointsize=12)
boxplot(myCombat3re, main="Normalized beta values", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue"))
dev.off()


dim(myCombat3re)
#Medians of controls
myCombat3reavg <- cbind.data.frame(rownames(myCombat3re), 
                                   myCombat3re[,c(1:10)],
                                   rowMedians(myCombat3re[,c(11:22)]))
head(myCombat3reavg)
##################################     Data analysed with controls median  ################################
colnames(myCombat3reavg) <- c("TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","Controls")
head(myCombat3reavg)
#LOAD MANIFEST
MANIFEST <- read.csv("/media/ankitv/Archivio1/testing_array/Loca/test/MethylationEPIC_v-1-0_B4_ext.csv")
#MANIFEST <- EPIC.manifest.hg19
colnames(MANIFEST)
dim(MANIFEST)
head(MANIFEST)
MERGE_myCombat3reavg <- merge(myCombat3reavg,MANIFEST,by.x="TargetID",by.y="IlmnID")
head(MERGE_myCombat3reavg,1)
MERGE_myCombat3reavg_pos <- MERGE_myCombat3reavg[,c(23,24,24,1:12)]
#check for dimensions
dim(MERGE_myCombat3reavg_pos)
head(MERGE_myCombat3reavg_pos)
tail(MERGE_myCombat3reavg_pos)
write.table(MERGE_myCombat3reavg_pos, "MERGE_myCombat3reavg_pos.txt", col.names = F, quote = F, row.names = F, sep = "\t")

awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' MERGE_myCombat3reavg_pos.txt | sort -k1,1 -k2,2n > MERGE_myCombat3reavg_pos_chr.txt
bedtools intersect -wa -wb -a MERGE_myCombat3reavg_pos_chr.txt -b ~/ref_av/human_ICR.bed > MERGE_myCombat3re.avg_human_ICR.txt
awk '{print $19"\t"$20"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' MERGE_myCombat3re.avg_human_ICR.txt > MERGE_myCombat3re.avg_human_ICR.rearranged.txt


library(ggpubr)
library(ggplot2)
library(plyr)
#Add column names as follows: "DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM", "Controls"
#PCA with 12 Controls Avg Human ICR
setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control")
comdata1 <- read.table("MERGE_myCombat3re.avg_human_ICR.rearranged.txt", header = F)
colnames(comdata1) <- c("DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM", "Controls")
rownames(comdata1)
colnames(comdata1)
head(comdata1)
Count_comdataICR1 <- count(comdata1, "DMR")
head(Count_comdataICR1)
comaggregate1 = aggregate(comdata1[,7:17],by=list(comdata1$DMR), mean)
head(comaggregate1, 2)
#Plot PCA by group color and labelling
comdataICR1=comaggregate1
rownames(comdataICR1)
comdataICR1[,1]
rownames(comdataICR1)=comdataICR1[,1]
rownames(comdataICR1)
colnames(comdataICR1)
comdataICR1 = comdataICR1[,-1]
head(comdataICR1)
dim(comdataICR1)
comdataICR1_counted <- cbind(comdataICR1, Count_comdataICR1)
head(comdataICR1_counted)
dim(comdataICR1_counted)
comdataICR1_counted1 <- comdataICR1_counted[which(comdataICR1_counted$freq >= 3),]
head(comdataICR1_counted1)
dim(comdataICR1_counted1)
#write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
comdataICR1_counted2 <- comdataICR1_counted1[,1:11]
head(comdataICR1_counted2)
summary(comdataICR1_counted2)
#df <- as.data.frame(dataICR1)
combdfICR <- comdataICR1_counted2
head(combdfICR)
dim(combdfICR)
combdfICR = data.frame(combdfICR)
#write.table(dfICR , "dfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
tcombdfICR = t(combdfICR)
tcombdfICR = data.frame(tcombdfICR)
#write.table(tdfICR , "tdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tcombdfICR)
tcombdfICR["Color"] <-  c("a1-FM","a2-FS","a3-FM","a4-AM","a5-AA","a6-AM","a7-PM","a8-PA","a9-PG","b10-PM", "cControls")
#head(tdfICR)
dim(tcombdfICR)
combdfx <-tcombdfICR[c(1:43)]
combPC<-prcomp(combdfx, center = TRUE, scale. = TRUE)
combPCi<-data.frame(combPC$x,Color=tcombdfICR$Color)
percentageICRcomb <- round(combPC$sdev^2 / sum(combPC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentageICRcomb <- paste( colnames(combPCi), "(", paste( as.character(percentageICRcomb), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

rcombp <-ggplot(combPCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentageICRcomb[1]) + ylab(percentageICRcomb[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("magenta","navy","navy","magenta","navy","skyblue","magenta","navy","skyblue","skyblue","darkgrey"))+
  scale_shape_manual(values=c(1,1,1,2,2,2,5,5,5,5,0))
rcombp <- rcombp +theme_classic()
rcombp + xlim(-20,20)+ ylim(-20,20)
#rcombp + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_combptdfICR_mediumCntrl.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

################################### Heatmap human ICR ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(comdataICR1_counted1)
dim(comdataICR1_counted1)
sum(comdataICR1_counted1$freq)

#Heatmap with 12 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
datacomICR2 = as.matrix(comdataICR1_counted2)
head(datacomICR2)
dim(datacomICR2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
datacomICR3 <- datacomICR2[,c(11,1:10)]
head(datacomICR3)
write.table(datacomICR3, "datacomICR3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
#Control Samples
#heatmap.2(dataICR2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#Let it cluster
svg(filename="Heatmap_com_avg_EPIC16.medCntrl_human_ICR.rearranged_SVG.svg", width=5, height=10, pointsize=12)
#heatmap.2(dataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(datacomICR3, col = colfunc, Colv = "NA",dendrogram = c("none"), trace = "none", keysize=0.8, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

svg(filename="Heatmap_com_avg_EPIC16.medcntrl_human_ICR.rearranged_SVG_clusterboth.svg", width=5, height=10, pointsize=12)
#heatmap.2(dataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(datacomICR3, col = colfunc, dendrogram = c("both"), trace = "none", key=TRUE, density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

#Box and Violin Plot human ICR 
dcavgata1 <- read.table("MERGE_myCombat3re.avg_human_ICR.rearranged.txt", header = F)
colnames(dcavgata1) <- c("DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM", "Controls")
rownames(dcavgata1)
colnames(dcavgata1)
head(dcavgata1)
Count_dcavgataICR1 <- count(dcavgata1, "DMR")
head(Count_dcavgataICR1)
dcavgaggregate1 = aggregate(dcavgata1[,7:17],by=list(dcavgata1$DMR), mean)
head(dcavgaggregate1, 2)
dcavgataICR1=dcavgaggregate1
rownames(dcavgataICR1)
dcavgataICR1[,1]
rownames(dcavgataICR1)=dcavgataICR1[,1]
rownames(dcavgataICR1)
colnames(dcavgataICR1)
dcavgataICR1 = dcavgataICR1[,-1]
head(dcavgataICR1)
dim(dcavgataICR1)
write.table(dcavgataICR1, "aggregated_dcavgataICR1.txt", sep="\t", quote = FALSE, append = FALSE)
comdataICR1_counted <- cbind(dcavgataICR1, Count_dcavgataICR1)
head(comdataICR1_counted)
dim(comdataICR1_counted)
comdataICR1_counted1 <- comdataICR1_counted[which(comdataICR1_counted$freq >= 3),]
head(comdataICR1_counted1)
dim(comdataICR1_counted1)
#write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
comdataICR1_counted2 <- comdataICR1_counted1[,1:11]
head(comdataICR1_counted2)
comVioICR <- data.frame(comdataICR1_counted2[,c(11,1:10)])
head(comVioICR)
write.table(comVioICR, "comVioICR.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_comVioICR_medcntrl.svg", width=10, height=5, pointsize=12)
boxplot(comVioICR, main="Average methylation at human ICRs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("grey","pink","red","darkgreen","pink","brown","pink","pink","blue","pink","pink","magenta","cyan","orange","cyan"))
dev.off()
head(comVioICR)
dim(comVioICR)
comVioICR1 <- comVioICR[,c(1,2,5,8,7,10,11,3,4,6,9)]
head(comVioICR1)
comVioICR2 <- stack(comVioICR1)
head(comVioICR2)
colnames(comVioICR2) <- c("Methylation", "Datasets")
head(comVioICR2)
ggplot(comVioICR2, aes(x=Datasets, y=Methylation, color=Datasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","pink","pink","pink","skyblue","skyblue","skyblue","navy","navy","navy","navy"))
ggsave("Violin_complot_EPIC16_humanICR.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#Statistics
summary(comVioICR)[c(2,3,5),]
write.table(summary(comVioICR)[c(2,3,5),], "summarycomVioICR.txt", append = F, quote = F)

wilcox.test(comVioICR$Controls, comVioICR$X1.FM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X2.FS, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X3.FM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X4.AM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X5.AA, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X6.AM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X7.PM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X8.PA, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X9.PG, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(comVioICR$Controls, comVioICR$X10.PM, paired = FALSE, alternative = "two.sided")$p.value

#Box and Violin Plot Global methylation
comdataGlob=myCombat3reavg
rownames(comdataGlob)
comdataGlob[,1]
comdataGlob1 = comdataGlob[,-1]
comdataGlob1 <- data.frame(comdataGlob1)
head(comdataGlob1)
str(comdataGlob1)
ViocomGlob <- comdataGlob1[,c(11,1:10)]
head(ViocomGlob)
write.table(ViocomGlob, "ViocomGlob.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

svg(filename="boxplot_ViocomGlob.svg", width=10, height=5, pointsize=12)
boxplot(ViocomGlob, main="Global average methylation", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("grey","pink","red","darkgreen","pink","brown","pink","pink","blue","pink","pink","magenta","cyan","orange","cyan"))
dev.off()
head(ViocomGlob)
ViocomGlob1 <- ViocomGlob[,c(1,2,5,8,7,10,11,3,4,6,9)]
head(ViocomGlob1)
ViocomGlob2 <- stack(ViocomGlob1)
colnames(ViocomGlob2) <- c("Methylation", "Datasets")
head(ViocomGlob2)
ggplot(ViocomGlob2, aes(x=Datasets, y=Methylation, color=Datasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","pink","pink","pink","skyblue","skyblue","skyblue","navy","navy","navy","navy"))
ggsave("Violincomg_plot_EPIC16_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)


#Statistics
summary(ViocomGlob)
write.table(summary(ViocomGlob)[c(2,3,5),], "summaryViocomGlob.txt", append = F, quote = F)

wilcox.test(ViocomGlob$Controls, ViocomGlob$X1.FM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X2.FS, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X3.FM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X4.AM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X5.AA, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X6.AM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X7.PM, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X8.PA, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X9.PG, paired = FALSE, alternative = "two.sided")$p.value
wilcox.test(ViocomGlob$Controls, ViocomGlob$X10.PM, paired = FALSE, alternative = "two.sided")$p.value

#ggboxplot(ViocomGlob2, x = "Datasets", y = "Methylation",
#          title = "Boxplot_EPIC16_Globalmeth", ylab = "Methylation",
#          palette =c("darkgrey","pink","red","darkgreen","pink","brown","pink","pink","blue","pink","pink","magenta","cyan","orange","cyan"),color = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic())+ylim(0,1)
#ggsave("Boxplot_EPIC16_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)







#*********************************   Individual samples ************************************
head(myCombat3re)
myiCombat3reprobed <- cbind.data.frame(rownames(myCombat3re),
                                       myCombat3re[,1:22])

head(myiCombat3reprobed)
colnames(myiCombat3reprobed) <- c("TargetID","1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6")
MERGE_myiCombat3re <- merge(myiCombat3reprobed,MANIFEST,by.x="TargetID",by.y="IlmnID")
head(MERGE_myiCombat3re)
MERGE_myiCombat3re_pos <- MERGE_myiCombat3re[,c(34,35,35,1:23)]
#check for dimensions
dim(MERGE_myiCombat3re_pos)
head(MERGE_myiCombat3re_pos)
tail(MERGE_myiCombat3re_pos)
write.table(MERGE_myiCombat3re_pos, "MERGE_myiCombat3re_pos.txt", col.names = F, quote = F, row.names = F)

awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26}' MERGE_myiCombat3re_pos.txt | sort -k1,1 -k2,2n > MERGE_myiCombat3re_pos_chr.txt
bedtools intersect -wa -wb -a MERGE_myiCombat3re_pos_chr.txt -b ~/ref_av/human_ICR.bed > MERGE_myiCombat3re_human_ICR.txt
awk '{print $30"\t"$31"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26}' MERGE_myiCombat3re_human_ICR.txt > MERGE_myiCombat3re_human_ICR.rearranged.txt

library(ggpubr)
library(ggplot2)
library(plyr)
#Add column names as follows: "DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM", "CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6"
#PCA with 6 Controls  Human ICR
setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control")
icomdata1 <- read.table("MERGE_myiCombat3re_human_ICR.rearranged.txt", header = F)
colnames(icomdata1) <- c("DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6")
rownames(icomdata1)
colnames(icomdata1)
head(icomdata1)
dim(icomdata1)
Count_icomdataICR1 <- count(icomdata1, "DMR")
head(Count_icomdataICR1)
dim(Count_icomdataICR1)
scomaggregate1 = aggregate(icomdata1[,7:28],by=list(icomdata1$DMR), mean)
head(scomaggregate1, 2)
#Plot PCA by group color and labelling
icomdataICR1=scomaggregate1
rownames(icomdataICR1)
icomdataICR1[,1]
rownames(icomdataICR1)=icomdataICR1[,1]
rownames(icomdataICR1)
colnames(icomdataICR1)
icomdataICR1 = icomdataICR1[,-1]
head(icomdataICR1)
dim(icomdataICR1)
icomdataICR1_counted <- cbind(icomdataICR1, Count_icomdataICR1)
head(icomdataICR1_counted)
dim(icomdataICR1_counted)
icomdataICR1_counted1 <- icomdataICR1_counted[which(icomdataICR1_counted$freq >= 3),]
head(icomdataICR1_counted1)
dim(icomdataICR1_counted1)
write.table(icomdataICR1_counted1, "icomdataICR1_counted1.filtmin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

icomdataICR1_counted2 <- icomdataICR1_counted1[,1:22]
head(icomdataICR1_counted2)
summary(icomdataICR1_counted2)
#df <- as.icomdata.frame(icomdataICR1)
icomdfICR <- icomdataICR1_counted2
head(icomdfICR)
dim(icomdfICR)
icomdfICR = data.frame(icomdfICR)
#write.table(icomdfICR , "icomdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
ticomdfICR = t(icomdfICR)
ticomdfICR = data.frame(ticomdfICR)
#write.table(ticomdfICR , "ticomdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(ticomdfICR)
ticomdfICR["Color"] <-  c("a1-FM","a2-FS","a3-FM","a4-AM","a5-AA","a6-AM","a7-PM","a8-PA","a9-PG","b10-PM","cCTRL01","cCTRL02","cCTRL03","cCTRL04","cCTRL05","cCTRL06","cCtrl1","cCtrl2","cCtrl3","cCtrl4","cCtrl5","cCtrl6")
#head(ticomdfICR)
#write.table(ticomdfICR , "ticomdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(ticomdfICR)
icomdfx <-ticomdfICR[c(1:43)]
PcomC<-prcomp(icomdfx, center = TRUE, scale. = TRUE)
PcomCi<-data.frame(PcomC$x,Color=ticomdfICR$Color)
percentagecomICR <- round(PcomC$sdev^2 / sum(PcomC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomICR <- paste( colnames(PcomCi), "(", paste( as.character(percentagecomICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pcom1 <-ggplot(PcomCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomICR[1]) + ylab(percentagecomICR[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("pink","navy","navy","pink","navy","skyblue","pink","navy","skyblue","skyblue","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  scale_shape_manual(values=c(1,1,1,2,2,2,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0))
pcom1 <- pcom1 +theme_classic()
pcom1 + xlim(-20,20)+ ylim(-20,20)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_ticomdfICR_CntrlIndiv.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

################################### Heatmap_indiv human ICR ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(icomdataICR1_counted1)
dim(icomdataICR1_counted1)
sum(icomdataICR1_counted1$freq)

#Heatmap_indiv with 12 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
icomdataICR2 = as.matrix(icomdataICR1_counted2)
head(icomdataICR2)
dim(icomdataICR2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
icomdataICR3 <- icomdataICR2[,c(11:22,1:10)]
head(icomdataICR3)
write.table(icomdataICR3, "Heatmap_indiv6_EPIC16_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
#Control Samples
#Heatmap_indiv.2(icomdataICR2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#Let it cluster
svg(filename="Heatmap_com_indiv_EPIC16.avg_human_ICR.rearranged_SVG.svg", width=5, height=10, pointsize=12)
#Heatmap_indiv.2(icomdataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(icomdataICR3, col = colfunc, Colv = "NA",dendrogram = c("none"), trace = "none", keysize=0.8, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

svg(filename="Heatmap_com_indiv_EPIC16.avg_human_ICR.rearranged_SVG_clusterboth.svg", width=15, height=15, pointsize=12)
#Heatmap_indiv.2(icomdataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(icomdataICR3, col = colfunc, dendrogram = c("both"), trace = "none", key=TRUE, density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

#Box and Violin Plot human ICR using 6 Controls from Dec 2018
setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control")
icomdata1 <- read.table("MERGE_myiCombat3re_human_ICR.rearranged.txt", header = F)
colnames(icomdata1) <- c("DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6")
rownames(icomdata1)
colnames(icomdata1)
head(icomdata1)
dim(icomdata1)
Count_icomdataICR1 <- count(icomdata1, "DMR")
head(Count_icomdataICR1)
icomagregateicr1 = aggregate(icomdata1[,7:28],by=list(icomdata1$DMR), mean)
head(icomagregateicr1, 2)
icomdataICR1=icomagregateicr1
rownames(icomdataICR1)
icomdataICR1[,1]
rownames(icomdataICR1)=icomdataICR1[,1]
rownames(icomdataICR1)
colnames(icomdataICR1)
icomdataICR1 = icomdataICR1[,-1]
head(icomdataICR1)
dim(icomdataICR1)
write.table(icomdataICR1, "icomaggregatedicr_icomdataICR1.txt", sep="\t", quote = FALSE, append = FALSE)
icomdataICR1_counted <- cbind(icomdataICR1, Count_icomdataICR1)
head(icomdataICR1_counted)
dim(icomdataICR1_counted)
icomdataICR1_counted1 <- icomdataICR1_counted[which(icomdataICR1_counted$freq >= 3),]
head(icomdataICR1_counted1)
dim(icomdataICR1_counted1)
write.table(icomdataICR1_counted1, "icomdataICR1_counted1.filtmin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
icomdataICR1_counted2 <- icomdataICR1_counted1[,1:22]
head(icomdataICR1_counted2)
ViocomIndvCR <- data.frame(icomdataICR1_counted2[,c(11:22,1:10)])
head(ViocomIndvCR)
write.table(ViocomIndvCR, "ViocomIndvCR.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_ViocomIndvICR.svg", width=10, height=5, pointsize=12)
boxplot(ViocomIndvCR, main="Average methylation at human ICRs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","navy","navy","pink","navy","skyblue","pink","navy","skyblue","skyblue"))
dev.off()
dim(ViocomIndvCR)
ViocomIndvCR <- data.frame(ViocomIndvCR)
ViocomIndvCR1 <- ViocomIndvCR[,c(1:12,13,16,19,18,21,22,14,15,17,20)]
head(ViocomIndvCR1)
ViocomIndvCR2 <- stack(ViocomIndvCR1)
head(ViocomIndvCR2)
colnames(ViocomIndvCR2) <- c("Methylation", "icomdatasets")
head(ViocomIndvCR2)
ggplot(ViocomIndvCR2, aes(x=icomdatasets, y=Methylation, color=icomdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","pink","pink","skyblue","skyblue","skyblue","navy","navy","navy","navy"))
ggsave("Violin_plot_ViocomIndvCR2_EPIC16_humanICR.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#Box and Violin Plot Global methylation 
setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control/")
head(myCombat3re)
dim(myCombat3re)
icomcomdataGlob <- myCombat3re
rownames(icomcomdataGlob)
icomcomdataGlob[,1]
icomcomdataGlob1 <- data.frame(icomcomdataGlob)
head(icomcomdataGlob1)
str(icomcomdataGlob1)
ViocomGlobIndv <- icomcomdataGlob1[,c(11:22,1:10)]
head(ViocomGlobIndv)
write.table(ViocomGlobIndv, "ViocomGlobIndv.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

svg(filename="boxplot_myCombat3re.svg", width=10, height=5, pointsize=12)
boxplot(ViocomGlobIndv, main="Global average methylation", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","navy","navy","pink","navy","skyblue","pink","navy","skyblue","skyblue"))
dev.off()
head(ViocomGlobIndv)
ViocomGlobIndv1 <- ViocomGlobIndv[,c(1:12,13,16,19,18,21,22,14,15,17,20)]
head(ViocomGlobIndv1)
ViocomGlobIndv2 <- stack(ViocomGlobIndv1)
colnames(ViocomGlobIndv2) <- c("Methylation", "icomdatasets")
head(ViocomGlobIndv2)
ggplot(ViocomGlobIndv2, aes(x=icomdatasets, y=Methylation, color=icomdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","pink","pink","skyblue","skyblue","skyblue","navy","navy","navy","navy"))
ggsave("Violin_plot_ViocomGlobIndv2_EPIC16_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#ggboxplot(ViocomGlobIndv2, x = "icomdatasets", y = "Methylation",
#          title = "Boxplot_indiv_EPIC16_Globalmeth", ylab = "Methylation",
#         palette =c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","blue", "magenta", "red", "orange","darkred", "darkorange", "darkgreen", "blue","pink","magenta"),color = "icomdatasets", bxp.errorbar = TRUE, ggtheme = theme_classic())+ylim(0,1)
#ggsave("Boxplot_com_indiv_EPIC16_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)


#Heatmap of control normalized and  Log for only padi6 cases
padidata1 <- read.table("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat/dataICR3ratioavg3.txt", header = TRUE)
rownames(padidata1)
colnames(padidata1)
head(padidata1)
padidata1[,1]
head(padidata1)
dim(padidata1)
padidata2 <- as.matrix(padidata1)
padidata2Log <- log2(padidata2)
head(padidata2Log)
write.table(padidata2Log , "/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat/padidata2Log.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
svg(filename="/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat/Heatmap_padidata2Log.normcontrolLog.svg", width=5, height=10, pointsize=12)
#heatmap.2(dataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(padidata2Log[,2:11], col = colfunc, Colv = "NA",dendrogram = c("row"), trace = "none", keysize=0.8, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(-2,2, length.out = 100))
dev.off()
padidata2Logminuscontrol <- padidata2Log[,c(2:11)]
head(padidata2Logminuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= c("a_1.FM","a_2.FS","a_3.FM","a_4.AM","a_5-AA","a_6-AM","a_7-PM","a_8-PA","a_9-PG","a_10-PM"))
row.names(my_sample_col2) <- colnames(padidata2Logminuscontrol)
my_colour2 = list(Annotations = c("a_1.FM" = "pink","a_2.FS" = "navy","a_3.FM"="navy","a_4.AM"="pink","a_5-AA"="navy","a_6-AM"="skyblue","a_7-PM"="pink","a_8-PA"="navy","a_9-PG"="skyblue","a_10-PM"="skyblue"))
breaksList2 = seq(-1, 1, by = 0.01)
pheatmap(padidata2Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList2)),
         breaks = breaksList2,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cutree_cols = 2)
#save as padidata2Logminuscontrol.svg
#Heatmap of control normalized Log other samples merged
head(icomdataICR1_counted2)
dim(icomdataICR1_counted2)
icomdataICR1_counted2ratio <- as.matrix(icomdataICR1_counted2)
head(icomdataICR1_counted2ratio)
dim(icomdataICR1_counted2ratio)
icomdataICR1_counted2ratioavg <- data.frame(cbind((rowMedians(icomdataICR1_counted2ratio[,11:22])),
                                                  (icomdataICR1_counted2ratio[,1:10])))

head(icomdataICR1_counted2ratioavg)
colnames(icomdataICR1_counted2ratioavg) <- c( "padilocacontrol", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM")


icomdataICR1_counted2ratioavg1 <- icomdataICR1_counted2ratioavg
icomdataICR1_counted2ratioavg2 <- as.matrix(icomdataICR1_counted2ratioavg1)
head(icomdataICR1_counted2ratioavg2)
icomdataICR1_counted2ratioavg3 <- data.frame(cbind((icomdataICR1_counted2ratioavg2[,1]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,2]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,3]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,4]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,5]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,6]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,7]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,8]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,9]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,10]/icomdataICR1_counted2ratioavg2[,1]),
                                                   (icomdataICR1_counted2ratioavg2[,11]/icomdataICR1_counted2ratioavg2[,1])))
head(icomdataICR1_counted2ratioavg3)
colnames(icomdataICR1_counted2ratioavg3) <- c("padilocacontrol", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM")

head(icomdataICR1_counted2ratioavg3)
icomdataICR1_counted2ratioavg3 <- as.matrix(icomdataICR1_counted2ratioavg3)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
svg(filename="Heatmap_icomdataICR1_counted2ratioavg.normcontrol.svg", width=15, height=8, pointsize=12)
heatmap.2(icomdataICR1_counted2ratioavg3,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="column", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(0,1, length.out = 100))
dev.off()

icomdataICR1_counted2ratioavg3Log <- log2(icomdataICR1_counted2ratioavg3)
head(icomdataICR1_counted2ratioavg3Log)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
svg(filename="Heatmap_icomdataICR1_counted2ratioavg.normcontrol_Log.svg", width=15, height=8, pointsize=12)
heatmap.2(icomdataICR1_counted2ratioavg3Log,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="both", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(-1,1, length.out = 100))

dev.off()
write.table(icomdataICR1_counted2ratioavg3Log , "Heatmap_Loca_Padi6_normcontrolicomdataICR1_counted2ratioavg2Log.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)

icomdataICR1_counted2ratioavg3Logminuscontrol <- icomdataICR1_counted2ratioavg3Log[,c(2:11)]
head(icomdataICR1_counted2ratioavg3Logminuscontrol)
library(pheatmap)
my_sample_col1 <- data.frame(Annotations= c("a_1.FM","a_2.FS","a_3.FM","a_4.AM","a_5-AA","a_6-AM","a_7-PM","a_8-PA","a_9-PG","a_10-PM"))

row.names(my_sample_col1) <- colnames(icomdataICR1_counted2ratioavg3Logminuscontrol)
my_colour1 = list(Annotations = c("a_1.FM" = "pink","a_2.FS" = "navy","a_3.FM"="navy","a_4.AM"="pink","a_5-AA"="navy","a_6-AM"="#9fc5e8ff","a_7-PM"="pink","a_8-PA"="navy","a_9-PG"="#9fc5e8ff","a_10-PM"="#9fc5e8ff"))
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(icomdataICR1_counted2ratioavg3Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)
colors1 <- seq(-1,1,by=0.01)
my_palette <- colorRampPalette(colors = c("navy", "white", "firebrick3"))(length(colors1))
pheatmap(icomdataICR1_counted2ratioavg3Logminuscontrol,
         color = my_palette, 
         breaks = colors1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cutree_cols = 2)
#save as pHeatmap_icomdataICR1_counted2ratioavg.normcontrol.svg
#Sort chromsome wise
head(icomdataICR1_counted2ratioavg3Logminuscontrol)
icomdataICR1_counted2ratioavg3Logminuscontrol_chr <- data.frame(icomdataICR1_counted2ratioavg3Logminuscontrol)
icomdataICR1_counted2ratioavg3Logminuscontrol_chr["Gene"] <- rownames(icomdataICR1_counted2ratioavg3Logminuscontrol_chr)
human_ICR <- read.table("/home/ankitv/ref_av/human_ICR_.bed", header = F)
head(human_ICR)
colnames(human_ICR) <- c("chr", "start", "end", "Gene", "DMR")
head(human_ICR)
head(icomdataICR1_counted2ratioavg3Logminuscontrol_chr)
icomdataICR1_counted2ratioavg3Logminuscontrol_chrgene = merge(icomdataICR1_counted2ratioavg3Logminuscontrol_chr, human_ICR, by="Gene", all.x=FALSE)
head(icomdataICR1_counted2ratioavg3Logminuscontrol_chrgene)
icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort <- icomdataICR1_counted2ratioavg3Logminuscontrol_chrgene[order(icomdataICR1_counted2ratioavg3Logminuscontrol_chrgene$chr, icomdataICR1_counted2ratioavg3Logminuscontrol_chrgene$start),]
head(icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort)
icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort  <- icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort[c(1:3, 34:43,4:33),]
head(icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort)
icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1 <- icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort[,1:11]
rownames(icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1) <- icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1[,1]
icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1 <- icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1[,-1]
head(icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1)
pheatmap(icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1,
         color = my_palette, 
         breaks = colors1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         cluster_rows = F,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cutree_cols = 2)
#save as pHeatmap_icomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1.svg

#normalize controls also and include in heatmap

head(icomdataICR1_counted2ratio)
dim(icomdataICR1_counted2ratio)
wicomdataICR1_counted2ratioavg <- data.frame(cbind((rowMedians(icomdataICR1_counted2ratio[,11:22])),
                                                   (icomdataICR1_counted2ratio[,1:10]),
                                                   (icomdataICR1_counted2ratio[,11:22])))

head(wicomdataICR1_counted2ratioavg)


colnames(wicomdataICR1_counted2ratioavg) <- c( "padilocacontrol", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","Ctrl1","Ctrl2","Ctrl3","Ctrl4","Ctrl5","Ctrl6")


wicomdataICR1_counted2ratioavg1 <- wicomdataICR1_counted2ratioavg
wicomdataICR1_counted2ratioavg2 <- as.matrix(wicomdataICR1_counted2ratioavg1)
head(wicomdataICR1_counted2ratioavg2)
wicomdataICR1_counted2ratioavg3 <- data.frame(cbind((wicomdataICR1_counted2ratioavg2[,1]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,2]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,3]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,4]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,5]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,6]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,7]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,8]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,9]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,10]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,11]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,12]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,13]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,14]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,15]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,16]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,17]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,18]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,19]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,20]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,21]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,22]/wicomdataICR1_counted2ratioavg2[,1]),
                                                    (wicomdataICR1_counted2ratioavg2[,23]/wicomdataICR1_counted2ratioavg2[,1])))
head(wicomdataICR1_counted2ratioavg3)
colnames(wicomdataICR1_counted2ratioavg3) <- c("padilocacontrol", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","Ctrl1","Ctrl2","Ctrl3","Ctrl4","Ctrl5","Ctrl6")

head(wicomdataICR1_counted2ratioavg3)
wicomdataICR1_counted2ratioavg3 <- as.matrix(wicomdataICR1_counted2ratioavg3)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
svg(filename="Heatmap_wicomdataICR1_counted2ratioavg.normcontrol.svg", width=15, height=8, pointsize=12)
heatmap.2(wicomdataICR1_counted2ratioavg3,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="column", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(0,1, length.out = 100))
dev.off()

wicomdataICR1_counted2ratioavg3Log <- log2(wicomdataICR1_counted2ratioavg3)
head(wicomdataICR1_counted2ratioavg3Log)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
svg(filename="Heatmap_wicomdataICR1_counted2ratioavg.normcontrol_Log.svg", width=15, height=8, pointsize=12)
heatmap.2(wicomdataICR1_counted2ratioavg3Log,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="both", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(-1,1, length.out = 100))

dev.off()
write.table(wicomdataICR1_counted2ratioavg3Log , "Heatmap_Loca_Padi6_normcontrolwicomdataICR1_counted2ratioavg2Log.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)

wicomdataICR1_counted2ratioavg3Logminuscontrol <- wicomdataICR1_counted2ratioavg3Log[,c(2:23)]
head(wicomdataICR1_counted2ratioavg3Logminuscontrol)
library(pheatmap)
my_sample_col1 <- data.frame(Annotations= c("a_1.FM","a_2.FS","a_3.FM","a_4.AM","a_5-AA","a_6-AM","a_7-PM","a_8-PA","a_9-PG","a_10-PM","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl","c_Cntrl"))

row.names(my_sample_col1) <- colnames(wicomdataICR1_counted2ratioavg3Logminuscontrol)
my_colour1 = list(Annotations = c("a_1.FM" = "pink","a_2.FS" = "navy","a_3.FM"="navy","a_4.AM"="pink","a_5-AA"="navy","a_6-AM"="#9fc5e8ff","a_7-PM"="pink","a_8-PA"="navy","a_9-PG"="#9fc5e8ff","a_10-PM"="#9fc5e8ff","c_Cntrl"="darkgrey"))
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(wicomdataICR1_counted2ratioavg3Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)
colors1 <- seq(-1,1,by=0.01)
my_palette <- colorRampPalette(colors = c("navy", "white", "firebrick3"))(length(colors1))
pheatmap(wicomdataICR1_counted2ratioavg3Logminuscontrol,
         color = my_palette, 
         breaks = colors1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cutree_cols = 2)
#save as pHeatmap_wicomdataICR1_counted2ratioavg.normcontrol.svg
#Sort chromsome wise
head(wicomdataICR1_counted2ratioavg3Logminuscontrol)
wicomdataICR1_counted2ratioavg3Logminuscontrol_chr <- data.frame(wicomdataICR1_counted2ratioavg3Logminuscontrol)
wicomdataICR1_counted2ratioavg3Logminuscontrol_chr["Gene"] <- rownames(wicomdataICR1_counted2ratioavg3Logminuscontrol_chr)
human_ICR <- read.table("/home/ankitv/ref_av/human_ICR_.bed", header = F)
head(human_ICR)
colnames(human_ICR) <- c("chr", "start", "end", "Gene", "DMR")
head(human_ICR)
head(wicomdataICR1_counted2ratioavg3Logminuscontrol_chr)
wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgene = merge(wicomdataICR1_counted2ratioavg3Logminuscontrol_chr, human_ICR, by="Gene", all.x=FALSE)
head(wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgene)
wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort <- wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgene[order(wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgene$chr, wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgene$start),]
head(wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort)
wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort  <- wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort[c(1:3, 34:43,4:33),]
head(wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort)
wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1 <- wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort[,1:23]
rownames(wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1) <- wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1[,1]
wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1 <- wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1[,-1]
head(wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1)
pheatmap(wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1,
         color = my_palette, 
         breaks = colors1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         cluster_rows = F,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cutree_cols = 2)
#save as pHeatmap_wicomdataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1.svg

#Selected samples

#After suggestion of Flavia mam filter somed samples X13.cntrl3 and X14.cntrl4
#Medians of controls
head(myCombat3re)
mycomselbat3reavgselected <- cbind.data.frame(rownames(myCombat3re), 
                                              myCombat3re[,c(1:10)],
                                              rowMedians(myCombat3re[,c(11:18, 21:22)]))
head(mycomselbat3reavgselected)
##################################     Data analysed with controls median  ################################
colnames(mycomselbat3reavgselected) <- c("TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","ControlsAC")
head(mycomselbat3reavgselected)
#LOAD MANIFEST
MANIFEST <- read.csv("/media/ankitv/Archivio1/testing_array/Loca/test/MethylationEPIC_v-1-0_B4_ext.csv")
#MANIFEST <- EPIC.manifest.hg19
colnames(MANIFEST)
dim(MANIFEST)
head(MANIFEST)
MERGE_mycomselbat3reavgselected <- merge(mycomselbat3reavgselected,MANIFEST,by.x="TargetID",by.y="IlmnID")
head(MERGE_mycomselbat3reavgselected)
MERGE_mycomselbat3reavgselected_pos <- MERGE_mycomselbat3reavgselected[,c(23,24,24,1:12)]
#check for dimensions
dim(MERGE_mycomselbat3reavgselected_pos)
head(MERGE_mycomselbat3reavgselected_pos)
tail(MERGE_mycomselbat3reavgselected_pos)
write.table(MERGE_mycomselbat3reavgselected_pos, "MERGE_mycomselbat3reavgselected_pos.txt", col.names = F, quote = F, row.names = F)

awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' MERGE_mycomselbat3reavgselected_pos.txt | sort -k1,1 -k2,2n > MERGE_mycomselbat3reavgselected_pos_chr.txt
bedtools intersect -wa -wb -a MERGE_mycomselbat3reavgselected_pos_chr.txt -b ~/ref_av/human_ICR.bed > MERGE_mycomselbat3re.avg_human_ICR.txt
awk '{print $19"\t"$20"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}' MERGE_mycomselbat3re.avg_human_ICR.txt > MERGE_mycomselbat3re.avg_human_ICR.rearranged.txt


library(ggpubr)
library(ggplot2)
library(plyr)
#Add column names as follows: "DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM", "Controls"
#PCA with 10 Controls Avg Human ICR
setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control")
comseldata1 <- read.table("MERGE_mycomselbat3re.avg_human_ICR.rearranged.txt", header = F)
colnames(comseldata1) <- c("DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM", "ControlsAC")
rownames(comseldata1)
colnames(comseldata1)
head(comseldata1)
Count_comseldataICR1 <- count(comseldata1, "DMR")
head(Count_comseldataICR1)
comselaggregate1 = aggregate(comseldata1[,7:17],by=list(comseldata1$DMR), mean)
head(comselaggregate1, 2)
#Plot PCA by group color and labelling
comseldataICR1=comselaggregate1
rownames(comseldataICR1)
comseldataICR1[,1]
rownames(comseldataICR1)=comseldataICR1[,1]
rownames(comseldataICR1)
colnames(comseldataICR1)
comseldataICR1 = comseldataICR1[,-1]
head(comseldataICR1)
dim(comseldataICR1)
comseldataICR1_counted <- cbind(comseldataICR1, Count_comseldataICR1)
head(comseldataICR1_counted)
dim(comseldataICR1_counted)
comseldataICR1_counted1 <- comseldataICR1_counted[which(comseldataICR1_counted$freq >= 3),]
head(comseldataICR1_counted1)
dim(comseldataICR1_counted1)
#write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
comseldataICR1_counted2 <- comseldataICR1_counted1[,1:11]
head(comseldataICR1_counted2)
summary(comseldataICR1_counted2)
#df <- as.data.frame(dataICR1)
comselbdfICR <- comseldataICR1_counted2
head(comselbdfICR)
dim(comselbdfICR)
comselbdfICR = data.frame(comselbdfICR)
#write.table(dfICR , "dfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
tcomselbdfICR = t(comselbdfICR)
tcomselbdfICR = data.frame(tcomselbdfICR)
#write.table(tdfICR , "tdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tcomselbdfICR)
tcomselbdfICR["Color"] <-  c("a1-FM","a2-FS","a3-FM","a4-AM","a5-AA","a6-AM","a7-PM","a8-PA","a9-PG","b10-PM", "cControlsAC")
#head(tdfICR)
dim(tcomselbdfICR)
comselbdfx <-tcomselbdfICR[c(1:43)]
comselbPC<-prcomp(comselbdfx, center = TRUE, scale. = TRUE)
comselbPCi<-data.frame(comselbPC$x,Color=tcomselbdfICR$Color)
percentageICRcomselb <- round(comselbPC$sdev^2 / sum(comselbPC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentageICRcomselb <- paste( colnames(comselbPCi), "(", paste( as.character(percentageICRcomselb), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

rcomselbp <-ggplot(comselbPCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentageICRcomselb[1]) + ylab(percentageICRcomselb[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("pink","navy","navy","pink","navy","#6fa8dcff","pink","navy","#6fa8dcff","#6fa8dcff","grey"))+
  scale_shape_manual(values=c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19))
rcomselbp <- rcomselbp +theme_classic()
rcomselbp + xlim(-15,15)+ ylim(-15,15)
#rcomselbp + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_comselbptdfICR_mediumCntrl.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

################################### Heatmap human ICR ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(comseldataICR1_counted1)
dim(comseldataICR1_counted1)
sum(comseldataICR1_counted1$freq)

#Heatmap with 12 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
datacomselICR2 = as.matrix(comseldataICR1_counted2)
head(datacomselICR2)
dim(datacomselICR2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
datacomselICR3 <- datacomselICR2[,c(11,1:10)]
head(datacomselICR3)
write.table(datacomselICR3, "datacomselICR3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
#Control Samples
#heatmap.2(dataICR2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#Let it cluster
svg(filename="Heatmap_com_avg_EPIC16.medCntrl_human_ICR.rearranged_SVG.svg", width=5, height=10, pointsize=12)
#heatmap.2(dataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(datacomselICR3, col = colfunc, Colv = "NA",dendrogram = c("none"), trace = "none", keysize=0.8, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

svg(filename="Heatmap_comsel_avg_EPIC16.medcntrl_human_ICR.rearranged_SVG_clusterboth.svg", width=5, height=10, pointsize=12)
#heatmap.2(dataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(datacomselICR3, col = colfunc, dendrogram = c("both"), trace = "none", key=TRUE, density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

#Box and Violin Plot human ICR 
dcavgselata1 <- read.table("MERGE_mycomselbat3re.avg_human_ICR.rearranged.txt", header = F)
head(dcavgselata1)
colnames(dcavgselata1) <- c("DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM", "ControlsAC")
rownames(dcavgselata1)
colnames(dcavgselata1)
head(dcavgselata1)
Count_dcavgselataICR1 <- count(dcavgselata1, "DMR")
head(Count_dcavgselataICR1)
dcavgselaggregate1 = aggregate(dcavgselata1[,7:17],by=list(dcavgselata1$DMR), mean)
head(dcavgselaggregate1, 2)
dcavgselataICR1=dcavgselaggregate1
rownames(dcavgselataICR1)
dcavgselataICR1[,1]
rownames(dcavgselataICR1)=dcavgselataICR1[,1]
rownames(dcavgselataICR1)
colnames(dcavgselataICR1)
dcavgselataICR1 = dcavgselataICR1[,-1]
head(dcavgselataICR1)
dim(dcavgselataICR1)
write.table(dcavgselataICR1, "aggregated_dcavgselataICR1.txt", sep="\t", quote = FALSE, append = FALSE)
comseldataICR1_counted <- cbind(dcavgselataICR1, Count_dcavgselataICR1)
head(comseldataICR1_counted)
dim(comseldataICR1_counted)
comseldataICR1_counted1 <- comseldataICR1_counted[which(comseldataICR1_counted$freq >= 3),]
head(comseldataICR1_counted1)
dim(comseldataICR1_counted1)
#write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
comseldataICR1_counted2 <- comseldataICR1_counted1[,1:11]
head(comseldataICR1_counted2)
#Remove L1, L2, L3 and L4
comselVioICR <- data.frame(comseldataICR1_counted2[,c(11,1:10)])
head(comselVioICR)
write.table(comselVioICR, "comselVioICR.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_comselVioICR_medcntrl.svg", width=10, height=5, pointsize=12)
boxplot(comselVioICR, main="Average methylation at human ICRs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("grey","pink","navy","navy","pink","navy","#9fc5e8ff","pink","navy","#9fc5e8ff","#9fc5e8ff"))
dev.off()
head(comselVioICR)
#rearrange the samples together as sugested  by Flavia mam
comselVioICR1  <- comselVioICR[,c(1,2,5,8,7,10,11,3,4,6,9)]
head(comselVioICR1)
comselVioICR2 <- stack(comselVioICR1)
head(comselVioICR2)
colnames(comselVioICR2) <- c("Methylation", "Datasets")
head(comselVioICR2)
ggplot(comselVioICR2, aes(x=Datasets, y=Methylation, color=Datasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","pink","pink","pink","#9fc5e8ff","#9fc5e8ff","#9fc5e8ff","navy","navy","navy","navy"))
ggsave("Violin_comselplot_EPIC16_humanICR.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#Box and Violin Plot Global methylation
comseldataGlob=mycomselbat3reavgselected
rownames(comseldataGlob)
comseldataGlob[,1]
comseldataGlob1 = comseldataGlob[,-1]
comseldataGlob1 <- data.frame(comseldataGlob1)
head(comseldataGlob1)
str(comseldataGlob1)
#Remove L1, L2, L3 and L4
ViocomselGlob <- comseldataGlob1[,c(11,1:10)]
head(ViocomselGlob)
write.table(ViocomselGlob, "ViocomselGlob.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

svg(filename="boxplot_ViocomselGlob.svg", width=10, height=5, pointsize=12)
boxplot(ViocomselGlob, main="Global average methylation", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","pink","navy","navy","pink","navy","#9fc5e8ff","pink","navy","#9fc5e8ff","#9fc5e8ff"))
dev.off()
head(comseldataGlob1)
#rearrange the samples together as sugested  by Flavia mam
ViocomselGlob1  <- ViocomselGlob[,c(1,2,5,8,7,10,11,3,4,6,9)]
head(ViocomselGlob1)
ViocomselGlob2 <- stack(ViocomselGlob1)
colnames(ViocomselGlob2) <- c("Methylation", "Datasets")
head(ViocomselGlob2)
ggplot(ViocomselGlob2, aes(x=Datasets, y=Methylation, color=Datasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() + scale_color_manual(values=c("darkgrey","pink","pink","pink","#9fc5e8ff","#9fc5e8ff","#9fc5e8ff","navy","navy","navy","navy"))
ggsave("Violincomselg_plot_EPIC16_selGlobalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#ggboxplot(ViocomselGlob2, x = "Datasets", y = "Methylation",
#          title = "Boxplot_EPIC16_Globalmeth", ylab = "Methylation",
#          palette =c("darkgrey","pink","red","darkgreen","pink","brown","pink","pink","blue","pink","pink","magenta","cyan","orange","cyan"),color = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic())+ylim(0,1)
#ggsave("Boxplot_EPIC16_selGlobalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)



#*********************************   Individual samples ************************************
head(myCombat3re)
myicomselbat3reprobed <- cbind.data.frame(rownames(myCombat3re),
                                          myCombat3re[,1:22])

head(myicomselbat3reprobed)
colnames(myicomselbat3reprobed) <- c("TargetID","1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6")
MERGE_myicomselbat3re <- merge(myicomselbat3reprobed,MANIFEST,by.x="TargetID",by.y="IlmnID")
head(MERGE_myicomselbat3re)
MERGE_myicomselbat3re_pos <- MERGE_myicomselbat3re[,c(34,35,35,1:23)]
#check for dimensions
dim(MERGE_myicomselbat3re_pos)
head(MERGE_myicomselbat3re_pos)
tail(MERGE_myicomselbat3re_pos)
write.table(MERGE_myicomselbat3re_pos, "MERGE_myicomselbat3re_pos.txt", col.names = F, quote = F, row.names = F)

awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26}' MERGE_myicomselbat3re_pos.txt | sort -k1,1 -k2,2n > MERGE_myicomselbat3re_pos_chr.txt
bedtools intersect -wa -wb -a MERGE_myicomselbat3re_pos_chr.txt -b ~/ref_av/human_ICR.bed > MERGE_myicomselbat3re_human_ICR.txt
awk '{print $30"\t"$31"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26}' MERGE_myicomselbat3re_human_ICR.txt > MERGE_myicomselbat3re_human_ICR.rearranged.txt

library(ggpubr)
library(ggplot2)
library(plyr)
#Add column names as follows: "DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM", "CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6"
#PCA with 10 Controls  Human ICR
setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control")
icomseldata1 <- read.table("MERGE_myicomselbat3re_human_ICR.rearranged.txt", header = F)
colnames(icomseldata1) <- c("DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6")
rownames(icomseldata1)
colnames(icomseldata1)
head(icomseldata1)
dim(icomseldata1)
Count_icomseldataICR1 <- count(icomseldata1, "DMR")
head(Count_icomseldataICR1)
dim(Count_icomseldataICR1)
scomselaggregate1 = aggregate(icomseldata1[,7:28],by=list(icomseldata1$DMR), mean)
head(scomselaggregate1, 2)
#Plot PCA by group color and labelling
icomseldataICR1=scomselaggregate1
rownames(icomseldataICR1)
icomseldataICR1[,1]
rownames(icomseldataICR1)=icomseldataICR1[,1]
rownames(icomseldataICR1)
colnames(icomseldataICR1)
icomseldataICR1 = icomseldataICR1[,-1]
head(icomseldataICR1)
dim(icomseldataICR1)
icomseldataICR1_counted <- cbind(icomseldataICR1, Count_icomseldataICR1)
head(icomseldataICR1_counted)
dim(icomseldataICR1_counted)
icomseldataICR1_counted1 <- icomseldataICR1_counted[which(icomseldataICR1_counted$freq >= 3),]
head(icomseldataICR1_counted1)
dim(icomseldataICR1_counted1)
#Filter samples after suggestions
icomseldataICR1_counted2 <- icomseldataICR1_counted1[,c(1:10, 11:18,21:22)]
head(icomseldataICR1_counted2)
summary(icomseldataICR1_counted2)

#df <- as.icomseldata.frame(icomseldataICR1)
icomseldfICR <- icomseldataICR1_counted2
head(icomseldfICR)
dim(icomseldfICR)
icomseldfICR = data.frame(icomseldfICR)
#write.table(icomseldfICR , "icomseldfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
ticomseldfICR = t(icomseldfICR)
ticomseldfICR = data.frame(ticomseldfICR)
#write.table(ticomseldfICR , "ticomseldfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(ticomseldfICR)
ticomseldfICR["Color"] <-  c("a1-FM","a2-FS","a3-FM","a4-AM","a5-AA","a6-AM","a7-PM","a8-PA","a9-PG","b10-PM", "cCTRL01","cCTRL02","cCTRL03","cCTRL04","cCTRL05","cCTRL06","cCtrl1","cCtrl2","cCtrl5","cCtrl6")
#head(ticomseldfICR)
#write.table(ticomseldfICR , "ticomseldfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(ticomseldfICR)
icomseldfx <-ticomseldfICR[c(1:43)]
PcomselC<-prcomp(icomseldfx, center = TRUE, scale. = TRUE)
PcomselCi<-data.frame(PcomselC$x,Color=ticomseldfICR$Color)
percentagecomselICR <- round(PcomselC$sdev^2 / sum(PcomselC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomselICR <- paste( colnames(PcomselCi), "(", paste( as.character(percentagecomselICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pcomsel1 <-ggplot(PcomselCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomselICR[1]) + ylab(percentagecomselICR[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("pink","navy","navy","pink","navy","skyblue","pink","navy","skyblue","skyblue","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  scale_shape_manual(values=c(1,1,1,2,2,2,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0))
pcomsel1 <- pcomsel1 +theme_classic()
pcomsel1 + xlim(-20,20)+ ylim(-20,20)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_ticomseldfICR_CntrlIndiv.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)


################################### Heatmap_indiv human ICR ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(icomseldataICR1_counted1)
dim(icomseldataICR1_counted1)
sum(icomseldataICR1_counted1$freq)

#Heatmap_indiv with 12 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
icomseldataICR2 = as.matrix(icomseldataICR1_counted2)
head(icomseldataICR2)
dim(icomseldataICR2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
icomseldataICR3 <- icomseldataICR2[,c(11:20,1:10)]
head(icomseldataICR3)
write.table(icomseldataICR3, "Heatmap_indiv6_EPIC16_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
#Control Samples
#Heatmap_indiv.2(icomseldataICR2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#Let it cluster
svg(filename="Heatmap_comsel_indiv_EPIC16.avg_human_ICR.rearranged_SVG.svg", width=5, height=10, pointsize=12)
#Heatmap_indiv.2(icomseldataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(icomseldataICR3, col = colfunc, Colv = "NA",dendrogram = c("none"), trace = "none", keysize=0.8, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

svg(filename="Heatmap_comsel_indiv_EPIC16.avg_human_ICR.rearranged_SVG_clusterboth.svg", width=15, height=15, pointsize=12)
#Heatmap_indiv.2(icomseldataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(icomseldataICR3, col = colfunc, dendrogram = c("both"), trace = "none", key=TRUE, density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

#Box and Violin Plot human ICR using 10 Controls
setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control")
icomseldata1 <- read.table("MERGE_myicomselbat3re_human_ICR.rearranged.txt", header = F)
head(icomseldata1)
colnames(icomseldata1) <- c("DMR", "DMRType","chr","start","end","TargetID", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6")
rownames(icomseldata1)
colnames(icomseldata1)
head(icomseldata1)
dim(icomseldata1)
Count_icomseldataICR1 <- count(icomseldata1, "DMR")
head(Count_icomseldataICR1)
icomselagregateicr1 = aggregate(icomseldata1[,7:28],by=list(icomseldata1$DMR), mean)
head(icomselagregateicr1, 2)
icomseldataICR1=icomselagregateicr1
rownames(icomseldataICR1)
icomseldataICR1[,1]
rownames(icomseldataICR1)=icomseldataICR1[,1]
rownames(icomseldataICR1)
colnames(icomseldataICR1)
icomseldataICR1 = icomseldataICR1[,-1]
head(icomseldataICR1)
dim(icomseldataICR1)
write.table(icomseldataICR1, "icomselaggregatedicr_icomseldataICR1.txt", sep="\t", quote = FALSE, append = FALSE)
icomseldataICR1_counted <- cbind(icomseldataICR1, Count_icomseldataICR1)
head(icomseldataICR1_counted)
dim(icomseldataICR1_counted)
icomseldataICR1_counted1 <- icomseldataICR1_counted[which(icomseldataICR1_counted$freq >= 3),]
head(icomseldataICR1_counted1)
dim(icomseldataICR1_counted1)
#write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
icomseldataICR1_counted2 <- icomseldataICR1_counted1[,1:22]
head(icomseldataICR1_counted2)

#Now filt samples x13.cntrl3 and X14.cntrl4 after Flavia mam suggestion
ViocomselIndvCR <- data.frame(icomseldataICR1_counted2[,c(11:18,21:22,1:10)])
head(ViocomselIndvCR)
write.table(ViocomselIndvCR, "ViocomselIndvCR.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_ViocomselIndvICR.svg", width=10, height=5, pointsize=12)
boxplot(ViocomselIndvCR, main="Average methylation at human ICRs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","navy","navy","pink","navy","skyblue","pink","navy","skyblue","skyblue"))
dev.off()
head(ViocomselIndvCR)
ViocomselIndvCR1 <- ViocomselIndvCR[,c(1:10,11,14,17,16,19,20,12,13,15,18)]
head(ViocomselIndvCR1)
ViocomselIndvCR2 <- stack(ViocomselIndvCR1)
head(ViocomselIndvCR2)
colnames(ViocomselIndvCR2) <- c("Methylation", "icomseldatasets")
head(ViocomselIndvCR2)
ggplot(ViocomselIndvCR2, aes(x=icomseldatasets, y=Methylation, color=icomseldatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","pink","pink","skyblue","skyblue","skyblue","navy","navy","navy","navy"))
ggsave("Violin_plot_ViocomselIndvCR2_EPIC16_humanICR.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#Box and Violin Plot Global methylation 
setwd("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat_control/")
head(myCombat3re)
mycomselbat3re <- myCombat3re
head(mycomselbat3re)
dim(mycomselbat3re)
#Filter samples as per suggestions X13 and X14
icomcomseldataGlob <- mycomselbat3re[,c(1:18,21:22)]
rownames(icomcomseldataGlob)
icomcomseldataGlob[,1]
icomcomseldataGlob1 <- data.frame(icomcomseldataGlob)
head(icomcomseldataGlob1)
str(icomcomseldataGlob1)
ViocomselGlobIndv <- icomcomseldataGlob1[,c(11:20,1:10)]
head(ViocomselGlobIndv)
write.table(ViocomselGlobIndv, "ViocomselGlobIndv.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

svg(filename="boxplot_mycomselbat3re.svg", width=10, height=5, pointsize=12)
boxplot(ViocomselGlobIndv, main="Global average methylation", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","navy","navy","pink","navy","skyblue","pink","navy","skyblue","skyblue"))
dev.off()
head(ViocomselGlobIndv)
ViocomselGlobIndv1 <-ViocomselGlobIndv[,c(1:10,11,14,17,16,19,20,12,13,15,18)]
head(ViocomselGlobIndv1)
ViocomselGlobIndv2 <- stack(ViocomselGlobIndv1)
colnames(ViocomselGlobIndv2) <- c("Methylation", "icomseldatasets")
head(ViocomselGlobIndv2)
ggplot(ViocomselGlobIndv2, aes(x=icomseldatasets, y=Methylation, color=icomseldatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","pink","pink","skyblue","skyblue","skyblue","navy","navy","navy","navy"))
ggsave("Violin_plot_ViocomselGlobIndv2_EPIC16_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)

#ggboxplot(ViocomselGlobIndv2, x = "icomseldatasets", y = "Methylation",
#          title = "Boxplot_indiv_EPIC16_Globalmeth", ylab = "Methylation",
#         palette =c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","blue", "magenta", "red", "orange","darkred", "darkorange", "darkgreen", "blue","pink","magenta"),color = "icomseldatasets", bxp.errorbar = TRUE, ggtheme = theme_classic())+ylim(0,1)
#ggsave("Boxplot_com_indiv_EPIC16_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)


#Heatmap of control normalized and  Log for only padi6 cases
#Nothing to filter, already 6 cntrls from new dataset
padidata1 <- read.table("/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat/dataICR3ratioavg3.txt", header = TRUE)
rownames(padidata1)
colnames(padidata1)
head(padidata1)
padidata1[,1]
head(padidata1)
dim(padidata1)
padidata2 <- as.matrix(padidata1)
padidata2Log <- log2(padidata2)
head(padidata2Log)
write.table(padidata2Log , "/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat/padidata2Logfilter.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
#svg(filename="/media/ankitv/Archivio2/ankit/EPIC_16/processing/ChAMP_idat/Heatmap_padidata2Log.normcontrolLog.svg", width=5, height=10, pointsize=12)
#heatmap.2(dataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#heatmap.2(padidata2Log[,2:11], col = colfunc, Colv = "NA",dendrogram = c("row"), trace = "none", keysize=0.8, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(-2,2, length.out = 100))
#dev.off()
padidata2Logminuscontrol <- padidata2Log[,c(2:11)]
head(padidata2Logminuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= c("a_1.FM","a_2.FS","a_3.FM","a_4.AM","a_5-AA","a_6-AM","a_7-PM","a_8-PA","a_9-PG","a_10-PM"))
row.names(my_sample_col2) <- colnames(padidata2Logminuscontrol)
my_colour2 = list(Annotations = c("a_1.FM" = "pink","a_2.FS" = "navy","a_3.FM"="navy","a_4.AM"="pink","a_5-AA"="navy","a_6-AM"="skyblue","a_7-PM"="pink","a_8-PA"="navy","a_9-PG"="skyblue","a_10-PM"="skyblue"))
breaksList2 = seq(-1, 1, by = 0.01)
pheatmap(padidata2Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList2)),
         breaks = breaksList2,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cutree_cols = 2)
#save as padidata2Logminuscontrol_sel.svg

#Heatmap of control normalized Log other samples merged
head(icomseldataICR1_counted2)
dim(icomseldataICR1_counted2)
#Filter samples as per Flavia mam suggestion:X13 and X14
icomseldataICR1_counted2ratio <- as.matrix(icomseldataICR1_counted2[,c(1:18,21:22)])
head(icomseldataICR1_counted2ratio)
dim(icomseldataICR1_counted2ratio)
icomseldataICR1_counted2ratioavg <- data.frame(cbind((rowMedians(icomseldataICR1_counted2ratio[,c(11:20)])),
                                                     (icomseldataICR1_counted2ratio[,1:10])))

head(icomseldataICR1_counted2ratioavg)
colnames(icomseldataICR1_counted2ratioavg) <- c( "padilocacontrol", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM")


icomseldataICR1_counted2ratioavg1 <- icomseldataICR1_counted2ratioavg
icomseldataICR1_counted2ratioavg2 <- as.matrix(icomseldataICR1_counted2ratioavg1)
head(icomseldataICR1_counted2ratioavg2)
icomseldataICR1_counted2ratioavg3 <- data.frame(cbind((icomseldataICR1_counted2ratioavg2[,1]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,2]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,3]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,4]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,5]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,6]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,7]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,8]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,9]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,10]/icomseldataICR1_counted2ratioavg2[,1]),
                                                      (icomseldataICR1_counted2ratioavg2[,11]/icomseldataICR1_counted2ratioavg2[,1])))
head(icomseldataICR1_counted2ratioavg3)
colnames(icomseldataICR1_counted2ratioavg3) <- c("padilocacontrol", "1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM")

head(icomseldataICR1_counted2ratioavg3)
icomseldataICR1_counted2ratioavg3 <- as.matrix(icomseldataICR1_counted2ratioavg3)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
svg(filename="Heatmap_icomseldataICR1_counted2ratioavg.normcontrol.svg", width=15, height=8, pointsize=12)
heatmap.2(icomseldataICR1_counted2ratioavg3,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="column", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(0,1, length.out = 100))
dev.off()

icomseldataICR1_counted2ratioavg3Log <- log2(icomseldataICR1_counted2ratioavg3)
head(icomseldataICR1_counted2ratioavg3Log)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
svg(filename="Heatmap_icomseldataICR1_counted2ratioavg.normcontrol_Logfilt.svg", width=15, height=8, pointsize=12)
heatmap.2(icomseldataICR1_counted2ratioavg3Log,trace = "none", col = colfunc , keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), dendrogram="both", cexRow=0.7, font=3, cexCol = 0.6, margins =c(8,6), breaks = seq(-1,1, length.out = 100))

dev.off()
write.table(icomseldataICR1_counted2ratioavg3Log , "Heatmap_Loca_Padi6_normcontrolicomseldataICR1_counted2ratioavg2Logfilt.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)

icomseldataICR1_counted2ratioavg3Logminuscontrol <- icomseldataICR1_counted2ratioavg3Log[,c(2:11)]
head(icomseldataICR1_counted2ratioavg3Logminuscontrol)
library(pheatmap)
my_sample_col1 <- data.frame(Annotations= c("a_1.FM","a_2.FS","a_3.FM","a_4.AM","a_5-AA","a_6-AM","a_7-PM","a_8-PA","a_9-PG","a_10-PM"))

row.names(my_sample_col1) <- colnames(icomseldataICR1_counted2ratioavg3Logminuscontrol)
my_colour1 = list(Annotations = c("a_1.FM" = "pink","a_2.FS" = "navy","a_3.FM"="navy","a_4.AM"="pink","a_5-AA"="navy","a_6-AM"="#9fc5e8ff","a_7-PM"="pink","a_8-PA"="navy","a_9-PG"="#9fc5e8ff","a_10-PM"="#9fc5e8ff"))
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(icomseldataICR1_counted2ratioavg3Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)
colors1 <- seq(-1,1,by=0.01)
my_palette <- colorRampPalette(colors = c("navy", "white", "firebrick3"))(length(colors1))
pheatmap(icomseldataICR1_counted2ratioavg3Logminuscontrol,
         color = my_palette, 
         breaks = colors1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cutree_cols = 2)
#save as pHeatmap_icomseldataICR1_counted2ratioavg.normcontrolfilt.svg
#Sort chromsome wise
head(icomseldataICR1_counted2ratioavg3Logminuscontrol)
icomseldataICR1_counted2ratioavg3Logminuscontrol_chr <- data.frame(icomseldataICR1_counted2ratioavg3Logminuscontrol)
icomseldataICR1_counted2ratioavg3Logminuscontrol_chr["Gene"] <- rownames(icomseldataICR1_counted2ratioavg3Logminuscontrol_chr)
human_ICR <- read.table("/home/ankitv/ref_av/human_ICR_.bed", header = F)
head(human_ICR)
colnames(human_ICR) <- c("chr", "start", "end", "Gene", "DMR")
head(human_ICR)
head(icomseldataICR1_counted2ratioavg3Logminuscontrol_chr)
icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgene = merge(icomseldataICR1_counted2ratioavg3Logminuscontrol_chr, human_ICR, by="Gene", all.x=FALSE)
head(icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgene)
icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort <- icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgene[order(icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgene$chr, icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgene$start),]
head(icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort)
icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort  <- icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort[c(1:3, 34:43,4:33),]
head(icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort)
icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1 <- icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort[,1:11]
rownames(icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1) <- icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1[,1]
icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1 <- icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1[,-1]
head(icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1)
pheatmap(icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1,
         color = my_palette, 
         breaks = colors1,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_col = my_sample_col1,
         cluster_rows = F,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cutree_cols = 2)
#save as pHeatmap_icomseldataICR1_counted2ratioavg3Logminuscontrol_chrgenesort1filt.svg


###################### NEW DMR regions in Locatelli data: CH t-test #####################à
#Hypomethylation CH-t-test
#Remove repeated CpGs: 
head(MERGE_myiCombat3re_pos)
dim(MERGE_myiCombat3re_pos)
fgrep -f CpGmayrepeated_cordsmerge.txt MERGE_myiCombat3re_pos_chr.txt -v > MERGE_myiCombat3re_pos_chr.repCpGminus.txt

awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26}' MERGE_myiCombat3re_pos_chr.repCpGminus.txt | grep chr > MERGE_myiCombat3re_pos_chr.repCpGminus.sep.txt

bedtools intersect -wa -wb -a hg19.2kb+1_corrected.txt -b MERGE_myiCombat3re_pos_chr.repCpGminus.sep.txt | awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29}' > MERGE_myiCombat3re_pos_chr.repCpGminus.sep.2kbp.txt

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#
iEPIC16_repCpGminus.2kb_1 <- read.table("MERGE_myiCombat3re_pos_chr.repCpGminus.sep.2kbp.txt", header = FALSE, stringsAsFactors = FALSE)
head(iEPIC16_repCpGminus.2kb_1)
dim(iEPIC16_repCpGminus.2kb_1)
#Calculate mean methylation in each 2kb bins (mean of beta-values)
aggregate_iEPIC16_repCpGminus.2kb_1 = aggregate(iEPIC16_repCpGminus.2kb_1[,6:27],by=list(iEPIC16_repCpGminus.2kb_1$V1), mean)
head(aggregate_iEPIC16_repCpGminus.2kb_1)
dim(aggregate_iEPIC16_repCpGminus.2kb_1)
aggregate_iEPIC16_repCpGminus.2kb_1_sorted = aggregate_iEPIC16_repCpGminus.2kb_1[order(aggregate_iEPIC16_repCpGminus.2kb_1$Group.1),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted)
rownames(aggregate_iEPIC16_repCpGminus.2kb_1_sorted)
rownames(aggregate_iEPIC16_repCpGminus.2kb_1_sorted)=aggregate_iEPIC16_repCpGminus.2kb_1_sorted[,1]
aggregate_iEPIC16_repCpGminus.2kb_1_sorted = aggregate_iEPIC16_repCpGminus.2kb_1_sorted[,-1]
aggregate_iEPIC16_repCpGminus.2kb_1_sorted = as.matrix(aggregate_iEPIC16_repCpGminus.2kb_1_sorted)
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted)

#Count number of CpGs in each regions and combine information with mean methylation
Count_iEPIC16_repCpGminus.2kb_1 <- count(iEPIC16_repCpGminus.2kb_1, "V1") 
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted <- cbind(aggregate_iEPIC16_repCpGminus.2kb_1_sorted, Count_iEPIC16_repCpGminus.2kb_1)
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted)
tail(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted)

#Extract 2kb regions with minimum CpGs
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted$freq >= 10),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff)
#write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff, "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff.txt", quote = FALSE, append = FALSE, sep = "\t", row.names = T)

#Reassign column names
colnames(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff) <- c("1-FM","2-FS","3-FM","4-AM","5-AA","6-AM","7-PM","8-PA","9-PG","10-PM","CTRL01","CTRL02","CTRL03","CTRL04","CTRL05","CTRL06","11-Ctrl1","12-Ctrl2","13-Ctrl3","14-Ctrl4","15-Ctrl5","16-Ctrl6","Regions","freq")
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff)
#write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff, "aggregate_Loca20_repCpGminus.2kb_1_sorted_counted_10cutoff.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff)
# Mean = X, 
#In apply function, 1 is row and 2 is column, This command is validated by simulation data, matrix <- matrix(1:12, nrow = 3)
#Rowwise effect.apply(matrix,1,mean), Columnwise.effect.apply(matrix,2,mean)
Control12mean.2kb <- data.frame(apply(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff[,11:22],1,mean)) #This command is validated by excel sheet manually
colnames(Control12mean.2kb) <- "Control12mean.2kb"
head(Control12mean.2kb)
#Rowwise standard deviation, apply(matrix,1,sd)
Control12sd.2kb <- data.frame(apply(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff[,11:22],1,sd)) # Standard Deviation = s, #This command is validated by excel sheet manually see copy. 1.7.19
colnames(Control12sd.2kb) <- "Control12sd.2kb"
head(Control12sd.2kb)
#Alignment checked
#Alignment_check <- cbind(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff[,1:4],
#                                                                        rownames(Control12mean.2kb),
#                                                                        Control12mean.2kb,
#                                                                        rownames(Control12sd.2kb),
#                                                                        Control12sd.2kb,
#                                                                        aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff[,25:26])
##write.table(Alignment_check, "Alignment_check.txt", quote = FALSE, append = FALSE, sep = "\t", row.names = T)
#awk '{print $1"\t"$6"\t"$8"\t"$10}' Alignment_check.txt | head/tail

aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd <- cbind(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff[,1:10],
                                                                            Control12mean.2kb,
                                                                            Control12sd.2kb,
                                                                            aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff[,23:24])

head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)


aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["2FS.C"] <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`2-FS`/aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$Control12mean.2kb
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`2FS.C` < 0.70),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo[,c(2,11,15)], 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)

aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`2FS.C` >= 2.0),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper[,c(2,11,15)], 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)


aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["3FM.C"] <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`3-FM`/aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$Control12mean.2kb
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`3FM.C` < 0.70),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo[,c(3,11,16)], 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)

aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`3FM.C` >= 2.0),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper[,c(3,11,16)], 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)

aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["5AA.C"] <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`5-AA`/aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$Control12mean.2kb
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`5AA.C` < 0.70),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo[,c(5,11,17)], 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)

aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`5AA.C` >= 2.0),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper[,c(5,11,17)], 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)

aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["8PA.C"] <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`8-PA`/aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$Control12mean.2kb
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`8PA.C` < 0.70),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo[,c(8,11,18)], 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)

aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`8PA.C` >= 2.0),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper[,c(8,11,18)], 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)





write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd, 
            "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.txt",
            quote = F,
            append = F,
            sep = "\t",
            row.names = T)
awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.txt | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.sep.txt


#Check overlap
awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo.txt  | sort -k1,1 -k2,2n | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo_sep.txt
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo_sep.txt -b ~/ref_av/human_ICR.bed 

#Overlap wit Hmole hg19 data (Hmole_CpG-10-less0.65.hg19.txt generated from lift over of /media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/2kbp+1/CpG-10-less0.65.sep.gene.txt)
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hypo_sep.txt -b Hmole_CpG-10-less0.65.hg19.txt 

awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper.txt  | sort -k1,1 -k2,2n | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper_sep.txt
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper_sep.txt -b ~/ref_av/human_ICR.bed 

#Overlap wit Hmole hg19 data (Hmole_CpG-10-less0.65.hg19.txt generated from lift over of /media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/2kbp+1/CpG-10-less0.65.sep.gene.txt)
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.2FS_hyper_sep.txt -b Hmole_CpG-10-less0.65.hg19.txt 

awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo.txt  | sort -k1,1 -k2,2n | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo_sep.txt
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo_sep.txt -b ~/ref_av/human_ICR.bed 

#Overlap wit Hmole hg19 data (Hmole_CpG-10-less0.65.hg19.txt generated from lift over of /media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/2kbp+1/CpG-10-less0.65.sep.gene.txt)
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hypo_sep.txt -b Hmole_CpG-10-less0.65.hg19.txt 


awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper.txt  | sort -k1,1 -k2,2n | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper_sep.txt
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper_sep.txt -b ~/ref_av/human_ICR.bed 

#Overlap wit Hmole hg19 data (Hmole_CpG-10-less0.65.hg19.txt generated from lift over of /media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/2kbp+1/CpG-10-less0.65.sep.gene.txt)
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.3FM_hyper_sep.txt -b Hmole_CpG-10-less0.65.hg19.txt 


awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo.txt  | sort -k1,1 -k2,2n | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo_sep.txt
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo_sep.txt -b ~/ref_av/human_ICR.bed 

#Overlap wit Hmole hg19 data (Hmole_CpG-10-less0.65.hg19.txt generated from lift over of /media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/2kbp+1/CpG-10-less0.65.sep.gene.txt)
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hypo_sep.txt -b Hmole_CpG-10-less0.65.hg19.txt 


awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper.txt  | sort -k1,1 -k2,2n | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper_sep.txt
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper_sep.txt -b ~/ref_av/human_ICR.bed 

#Overlap wit Hmole hg19 data (Hmole_CpG-10-less0.65.hg19.txt generated from lift over of /media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/2kbp+1/CpG-10-less0.65.sep.gene.txt)
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.5AA_hyper_sep.txt -b Hmole_CpG-10-less0.65.hg19.txt 


awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo.txt  | sort -k1,1 -k2,2n | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo_sep.txt
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo_sep.txt -b ~/ref_av/human_ICR.bed 

#Overlap wit Hmole hg19 data (Hmole_CpG-10-less0.65.hg19.txt generated from lift over of /media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/2kbp+1/CpG-10-less0.65.sep.gene.txt)
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hypo_sep.txt -b Hmole_CpG-10-less0.65.hg19.txt 

awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper.txt  | sort -k1,1 -k2,2n | grep chr > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper_sep.txt
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper_sep.txt -b ~/ref_av/human_ICR.bed 

#Overlap wit Hmole hg19 data (Hmole_CpG-10-less0.65.hg19.txt generated from lift over of /media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/2kbp+1/CpG-10-less0.65.sep.gene.txt)
bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd.8PA_hyper_sep.txt -b Hmole_CpG-10-less0.65.hg19.txt 

n <- 12
v <- sqrt((n + 1)/n)
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["DenoBase"] <- Control12sd.2kb * v #Calculate Denominator
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
#Apply CH t-test
# tch = (x* -X) / s * sqrt((n+1)/n) 


#2FS
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["xC2FS"] <- data.frame(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$`2-FS` - aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$Control12mean.2kb)  #x* - X, where x* =Mean methylation of each 2kb regions in the case sample
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["txC2FS"] <- data.frame(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$xC2FS/aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$DenoBase)				#Calculate t-values, txCL/DenoBase
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
#pt(tval), df=degfree)) #one-tailed  #Dr Claudia Suggested and explained to me
#  pval <- 2*(1-pt(abs(tval), df=degfree)) #two-tailed p-value, to identify both hypo and hypermeth regions
#aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["ptxCL"] <- data.frame(pt((aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$txCL), df=n-1))	#Calculate p-values
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["ptxC2FS"] <- data.frame(2*(1-pt(abs(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$txC2FS), df=n-1)))	#Calculate p-values
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd["padjptxC2FS"] <- data.frame(p.adjust(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$ptxC2FS, method = "BH"))		#Calculate p-adj values, BH multiplicity correction
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd)

#Export 2FS p<0.05
aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd_p2FSless0.05 <- aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd[which(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd$padjtxCL < 0.05),]
head(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd_p2FSless0.05)
dim(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd_p2FSless0.05)
write.table(aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd_p2FSless0.05, "aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd_p2FSless0.05.txt", quote = FALSE, append = FALSE, sep = "\t", row.names =F)
#grep chr aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd_p2FSless0.05.txt | awk '{print $7"\t"$13}' | awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n > aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd_p2FSless0.05.bed
#bedtools intersect -wa -wb -a aggregate_iEPIC16_repCpGminus.2kb_1_sorted_counted_10cutoff_meansd_p2FSless0.05.bed -b ~/ref_av/human_ICR.bed -v
#Assign genes to novel DMRs
#bedtools intersect -wa -wb -a novel_DMRat10CpGcutoff.txt -b /home/ankitv/ref_av/gencodes/human/hg19/human_gene_cords.bed 
sessionInfo()
