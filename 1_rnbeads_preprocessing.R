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
