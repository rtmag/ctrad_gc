library(RnBeads)

clinical = read.csv("GCF_clinical_data_edit.csv",stringsAsFactors=F)
sample = read.table("sample_sheet.txt",sep = "\t",header=T,stringsAsFactors=F)
sample = sample[,c(1,6,7,4)]
colnames(sample) = c("Sample_ID","Sentrix_ID","Sentrix_Position","Sample_Plate")
gsub("-T","",sample[,1])


# Multiprocess
num.cores <- 23
parallel.setup(num.cores)
## prepro
idat.dir <- file.path("/home/rtm/backup/CSI/CSI-Stephanie_Met450K/Raw_Data")
print(idat.dir)

	sample.annotation <- file.path("/home/rtm/backup/CSI/CSI-Stephanie_Met450K/sample_annotation.txt")
	print(sample.annotation)

	report.dir <- file.path("/home/rtm/backup/CSI/CSI-Stephanie_Met450K/RnBeads")
	print(report.dir)

	rnb.options(import.table.separator="\t")
	data.source <- c(idat.dir, sample.annotation)
	result <- rnb.run.import(data.source=data.source,data.type="infinium.idat.dir", dir.reports=report.dir)
	rnb.set.norm <- rnb.execute.normalization(result$rnb.set, method="swan",bgcorr.method="methylumi.noob")

	save.rnb.set(rnb.set.norm,path="/home/rtm/backup/CSI/CSI-Stephanie_Met450K/rnb.set.norm.RData")
print("DONE")
