#!/usr/local/bin/Rscript --vanilla

#
# Usage: cummerbund_isoform_analysis.R [data directory to all cuffdiff runs]
#
#
# Example:  cummerbund_isoform_analysis.R /local/projects-dr/ifx_core/tcreasy/XGWIE_ioannis/figures_for_paper/cummerbund/raw_data output_dir
# 
#
# tcreasy@som.umaryland.edu 
#
###############################################################################


rm(list=ls(all=T))

# input arguments
args <- commandArgs(TRUE)

stopifnot(length(args) > 1)

#Reading input data directory
data.dir <- args[1]

# Initialize outout directory
sOutDir <- args[2]

dir.list <- list.dirs(data.dir, recursive=F)
#print(dir.list)

cat("\nLoading cummeRbund...\n")
suppressMessages(library(cummeRbund))

#
# Function to plot number of reads per transcript
#
readsPerTranscripts <- function(cuff, samples, sOutDir) {

	pdf.file <- paste(sOutDir,"/" ,samples[1],".",samples[2],".reads_per_transcript." ,"pdf", sep="")
	pdf(pdf.file)

	sigIsoformIds <- getSig(cuff, alpha=1.0, level="isoforms")

	# for debugging for faster processing (only grab 1000 isoforms)
	#sigIsoformIds <- sigIsoformIds[1:1000]

	print("Number of isoforms to plot: ")
	print(length(sigIsoformIds))
	
	sigIsoforms <- getFeatures(cuff, sigIsoformIds, level="isoforms")
	isoform.counts <- count(sigIsoforms)
	
    # isoforms from sample 1 (q1)
	q1.isoforms <- subset(isoform.counts, sample_name == samples[2])
	q1.isoforms <- q1.isoforms[order(q1.isoforms$count, decreasing=TRUE),] 
	q1.matrix <- as.matrix(q1.isoforms$count)
	print(paste("Maximum read count",samples[2],sep=" "))
	print(head(q1.matrix))

    # isoforms from sample 2 (q2)
	q2.isoforms <- subset(isoform.counts, sample_name == samples[1])
	q2.isoforms <- q2.isoforms[order(q2.isoforms$count, decreasing=TRUE),] 
	q2.matrix <- as.matrix(q2.isoforms$count)
	print(paste("Maximum read count",samples[1],sep=" "))
	print(head(q2.matrix))

	# layout for plots
	layout(matrix(c(1,2), nrow=2, ncol=1))
	par(pty="s", mar=c(2,2,2,2), new=TRUE)

    # plot for transcripts from sample 1
	q1.title = paste("Distribution of Reads per Transcript: ", as.character(samples[1]), sep="")
	plot(q1.matrix, col="blue", log="y", pch=20, cex=.5, main=q1.title, ylab="Normalized Read Counts", xlab="Transcripts", ylim=c(1, 50000), add=TRUE)
	legend.txt <- paste("# of Isoforms:", length(q1.isoforms$count), sep=" ")
	mtext(legend.txt, line=-5, adj=0.7)

    # plot for transcripts from sample 2
	q2.title = paste("Distribution of Reads per Transcript: ", as.character(samples[2]), sep="")
	plot(q2.matrix, col="red", log="y", pch=20, cex=.5, main=q2.title, ylab="Normalized Read Counts", xlab="Transcripts", ylim=c(1, 50000), add=TRUE)
	legend.txt <- paste("# of Isoforms:", length(q2.isoforms$count), sep=" ")
	mtext(legend.txt, line=-5, adj=0.7)

	dev.off()
}


#
# Function to produce cummeRbund plots
#
runCummerbund <- function(cuff, samples, sOutDir) {

	pdf.file <- paste(sOutDir,"/" ,samples[1],".",samples[2],"." ,"pdf", sep="")
	print(pdf.file)
	pdf(pdf.file)
	
	# all isoforms
	allIsoforms <- cuff@isoforms

	# significant isoforms (FDR not considered, only removing test != OK)
	sigIsoformIds <- getSig(cuff, alpha=1.0, level="isoforms")

	# for debugging for faster processing (only grab 1000 isoforms)
	sigIsoformIds <- sigIsoformIds[1:1000]

	sigIsoforms <- getFeatures(cuff, sigIsoformIds, level="isoforms")
	
	print( csDensity(allIsoforms, replicates=FALSE ) )
	print( csBoxplot(allIsoforms, replicates=FALSE ) )
	#print( csDendro(allIsoforms, replicates=FALSE ) )
	print( csScatterMatrix(allIsoforms, replicates=FALSE) )
	print( MAplot(allIsoforms, samples[2], samples[1], logMode=TRUE, pseudocount=1.0) )
	print( csVolcano(allIsoforms, samples[2], samples[1], alpha=0.05, showSignificant=TRUE ) )	
	#print( PCAplot(allIsoforms, samples[2], samples[1]) )


    #
    # Heatmap plots useless until we can narrow down the # of isoforms
    #
    #print( csHeatmap(sigIsoforms, cluster='both') )
    #print( expressionBarplot(sigIsoforms) )


	
	dev.off()
}


for( cuff.dir in dir.list ) {

	cat("\nProcessing cuffdiff data from: ")
	cat(cuff.dir)
	cat("\n")

    # get sample names from directory name (this has expectations that the directory name is of the form "sample1vssample2")
	tmp <- gsub(".*/", "", as.character(cuff.dir), perl=TRUE)

	# hacky?
	base.name <- gsub(".sample", "", tmp, perl=TRUE)
	samples <- unlist( strsplit(base.name, "vs" ) )

	cat("\nReading cuffdiff data...\n")
	#cuff <- readCufflinks(dir=cuff.dir, rebuild=TRUE)
	cuff <- readCufflinks(dir=cuff.dir, rebuild=FALSE)

	cat("\nPlot isoform data...\n")
	
	# produce cummerbund plots
	runCummerbund(cuff, samples, sOutDir)

	#
	# CUSTOM PLOT FUNCTIONS
	#

	# plot the number of reads / transcript
	readsPerTranscripts(cuff, samples,sOutDir)
}

quit()


