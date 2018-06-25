###############################################################################################################################
# Annotation parsing

# 120 matched normal + tumor 450k
# 50 WES tumor

clinical = read.csv("/home/rtm/ctrad/GC_Illumina02/DNA_Methylation/GCF_clinical_data_edit.csv",stringsAsFactors=F)
sample = read.table("/home/rtm/ctrad/GC_Illumina02/DNA_Methylation/sample_sheet.txt",sep = "\t",header=T,stringsAsFactors=F)
sample = sample[,c(1,6,7,4)]
colnames(sample) = c("Sample_ID","Sentrix_ID","Sentrix_Position","Sample_Plate")
sampleNames = gsub("-T","",sample[,1])

ix = match(sampleNames,clinical[,1])
Tumor = gsub(".+\\-","",sample[,1],perl="T")
annotation = cbind(sample,clinical[ix,],Tumor)
annotation[annotation=="N.A."]=NA
row.names(annotation) = annotation[,1]

write.csv(annotation,"/home/rtm/ctrad/GC_Illumina02/DNA_Methylation/sample_annotation_edited.csv",row.names=F)
###############################################################################################################################
# Normalization and pre processing
library(RnBeads)

# Multiprocess
num.cores <- 20
parallel.setup(num.cores)

rnb.options(assembly = "hg38")
## prepro
idat.dir <- file.path("/home/rtm/ctrad/GC_Illumina02/DNA_Methylation/idat_files/")
print(idat.dir)

	sample.annotation <- file.path("/home/rtm/ctrad/GC_Illumina02/DNA_Methylation/sample_annotation_edited.csv")
	print(sample.annotation)

        system("rm -fr /home/rtm/ctrad/GC_Illumina02/DNA_Methylation/RnBeads")
        system("mkdir /home/rtm/ctrad/GC_Illuqmina02/DNA_Methylation/RnBeads")
	report.dir <- file.path("/home/rtm/ctrad/GC_Illumina02/DNA_Methylation/RnBeads")
	print(report.dir)

	rnb.options(import.table.separator=",")
	data.source <- c(idat.dir, sample.annotation)
	result <- rnb.run.import(data.source=data.source,data.type="infinium.idat.dir", dir.reports=report.dir)
	rnb.set.norm <- rnb.execute.normalization(result$rnb.set, method="swan",bgcorr.method="methylumi.noob")

	save.rnb.set(rnb.set.norm,path="/home/rtm/ctrad/GC_Illumina02/DNA_Methylation/rnb.set.norm.RData")
print("DONE")
###############################################################################################################################

suppressMessages(library(RnBeads))
rnb.set.norm=load.rnb.set("rnb.set.norm.RData.zip")
rnb.set.norm_no12=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which((rnb.set.norm@pheno$Tumor) %in% 1:2)])
TvsN_dmc <- rnb.execute.computeDiffMeth(rnb.set.norm_no12,pheno.cols=c("Tumor"))
#
comparison <- get.comparisons(TvsN_dmc)[1]
noNormal_dmr_table <-get.table(TvsN_dmc, comparison, "sites", return.data.frame=TRUE)

meth.norm<-meth(rnb.set.norm)

colnames(meth.norm) = as.character(rnb.set.norm@pheno$Tumor)
rownames(meth.norm) = rownames(rnb.set.norm@sites)

meth.norm.sig=meth.norm[noNormal_dmr_table$diffmeth.p.adj.fdr<0.05 & abs(noNormal_dmr_table[,3])>.15,]

options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

track=colnames(meth.norm.sig)
track[track=="N"]=1
track[track=="T"]=2
track[track=="1"]=3
track[track=="2"]=3
track=as.numeric(track)

colores=c("#db4e68","#497bd1","#d1c349")
clab=as.character(colores[track])

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))
meth.norm.sig = meth.norm.sig[complete.cases(meth.norm.sig),]
pdf("heatmap.pdf")
x = heatmap.2(as.matrix(meth.norm.sig),col=colors,scale="none", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,labCol = "",xlab="", ylab="CpGs",key.title="Methylation lvl",ColSideColors=clab)
dev.off()
