#!/usr/local/bin/Rscript --vanilla

#
# Usage:     expression_plots.R list.file outdir LFC_cutoff FDR_cutoff analysis_type(deseq or cuffdiff)
#
#
# Example:   expression_plots.R list.file test 2 0.06 cuffdiff
# 
#
# pkumari@som.umaryland.edu 
#
###############################################################################

# Read command line options
cat("\nReading Command-line arguments .....\n")
aCmdLineOptions <- commandArgs(TRUE)

stopifnot(length(aCmdLineOptions) > 1)

cat("\nReading List file .....\n")
sSampleInfo = aCmdLineOptions[1]

cat("\nReading output directory .....\n")
sOutDir = aCmdLineOptions[2]

cat("\nReading Log Fold Change Cutoff .....\n")
LFC_cutoff <- as.numeric(aCmdLineOptions[3])

cat("\nReading FDR Cutoff .....\n")
FDR_cutoff <- as.numeric(aCmdLineOptions[4])

analysis <- as.character(aCmdLineOptions[5])

if (is.na(aCmdLineOptions[3]) ) {LFC_cutoff <- 1}
if (is.na(aCmdLineOptions[4]) ) {FDR_cutoff <- 0.05}

print(paste("LFC Cutoff",LFC_cutoff,sep="="))
print(paste("FDR Cutoff",FDR_cutoff,sep="="))

require(ggplot2)



DET_plot <- function(oD, samples, sOutDir, LFC_cutoff, FDR_cutoff, analysis) {
	#png.file <- paste(samples[1], samples[2],"DET", "png", sep=".")
	#png(png.file, width=800, height=800)
	heading1<- paste(samples[1], samples[2], sep=" vs ")
	heading2<- paste(heading1," (Abs(LFC) >",LFC_cutoff," and FDR <",FDR_cutoff,")", sep="")
	
	for(i in 1:nrow(oD)){ 
		if (oD$LFC[i] > LFC_cutoff & oD$FDR[i] < FDR_cutoff ){oD$Significant[i] <- "True up"} 
		else if (oD$LFC[i] < -(LFC_cutoff) & oD$FDR[i] < FDR_cutoff) {oD$Significant[i] <- "True down"} 
		else {oD$Significant[i] <- "False"}
	}
	if (analysis == "cuffdiff") {
		xlable <- paste(samples[2],"(log10(FPKM))",sep= ":")
		ylable <- paste(samples[1],"(log10(FPKM))",sep= ":")
	}
	if (analysis == "deseq") {
		xlable <- paste(samples[2],"(log10(RC))",sep= ":")
		ylable <- paste(samples[1],"(log10(RC))",sep= ":")
	}
	
	deg <-  ggplot(data=oD, aes(x=log10(value_1) , y = log10(value_2),colour=Significant)) + geom_point(alpha=0.5,size=1.00)+ xlab(xlable) + ylab(ylable)+ labs(title=heading2) + scale_colour_manual(values = c("red","blue", "green"))
	
	print(deg)
	
}

logFoldvsFDR <- function(oD, samples, sOutDir, LFC_cutoff, FDR_cutoff) {
	#png.file <- paste(samples[1], samples[2], "png", sep=".")
	#png(png.file, width=800, height=800)
	heading1<- paste(samples[1], samples[2], sep=" vs ")	
	heading2<- paste(heading1," (Abs(LFC) >",LFC_cutoff," and FDR <",FDR_cutoff,")", sep="")

	oD$Significant = as.factor(abs(oD$LFC) > LFC_cutoff & oD$FDR < FDR_cutoff)
	vol_plot <- ggplot(data=oD, aes(x=LFC , y = -log10(FDR),colour=Significant)) + geom_point(alpha=0.5,size=1.00)+ xlim(c(-10, 10)) + ylim(c(0, 15))+ xlab("log2 fold change") + ylab("-log10 FDR") + labs(title=heading2)		
	
	print(vol_plot)
	
	
}



cat("\nReading contents of input file", sSampleInfo, ".....\n")
oS = read.delim(sSampleInfo, sep="", header=F)
nI = 1
pdf(paste(sOutDir,"/","Expression_Plots.pdf",sep=""))
while (nI <= nrow(oS)) {
        sFile = as.character(oS[nI,1])

	cat("\nReading File content for", sFile, ".....\n")
	oD = read.delim(sFile, sep="\t", header=TRUE)
	cat("\n")

	# determine sample names from the file name
	tmp <- gsub(".*/", "", as.character(sFile), perl=TRUE)
	if (analysis == "cuffdiff") {
		base.name <- gsub(".sample.*", "", tmp, perl=TRUE)
		samples <- unlist( strsplit(base.name, "vs" ))
		colnames(oD)[10] <- "LFC"
		colnames(oD)[13] <- "FDR"
	} else if (analysis == "deseq") {
		colnames(oD)[6] <- "LFC"
		colnames(oD)[3] <- "value_1"
		colnames(oD)[4] <- "value_2"
		base.name <- gsub(".de_genes.*", "", tmp, perl=TRUE)
		samples <- unlist( strsplit(base.name, "_vs_" ))
	}


	#Plot volcano plot
	cat("\nVolcano plot for log10(FDR) vs log fold change...\n")
	logFoldvsFDR(oD, samples, sOutDir, LFC_cutoff, FDR_cutoff)

	#Plot for up/down/non DET
	cat("\nPlot for log10(sample 1) vs log10(sample 2)...\n")
	DET_plot(oD, samples, sOutDir, LFC_cutoff, FDR_cutoff, analysis)
	
	nI = nI + 1
}
dev.off()
quit()




