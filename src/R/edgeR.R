#!/usr/local/packages/R-2.14.1/bin/Rscript --vanilla

# arguments:
# 1: sample table
# 2: file output path
# 3: FDR cutoff (DEFAULT: < 0.05)

###############################################################################

cat("\n***** Starting edgeR (v3.0.6) analysis ******\n\n")

suppressMessages( library(edgeR) )

## #########################
## data loading function
loadCountsTable <- function(sample.table,count.colnum=2){

  ## ################
  ## read all files
  raw.counts.files <- lapply(1:nrow(sample.table),function(x){
    read.table( sample.table$path[x],stringsAsFactors=FALSE)
  })

  ## #########################################
  ## check to make sure the same genes are in the files
  allgenes <- raw.counts.files[[1]][,1]
  for( i in 2:length(raw.counts.files) ){
    if( !all(raw.counts.files[[i]][,1]==allgenes) ){
      stop("The list of genes in the first column of file ",i,
           "does not exactly match that of file 1.\n" )
    }
  }

  ## ################################
  ## make a nice table
  counts.table <- sapply(raw.counts.files,function(x) x[,count.colnum])
  rownames(counts.table) <- allgenes
  colnames(counts.table) <- sample.table$SID

  return(counts.table)
}


## object prep functions
seqprep <- function(ctable, sampgroups, normmethod="TMM") {	
    seqobj <- edgeRprep(ctable, sampgroups, normmethod=normmethod)
	return(seqobj)
}


edgeRprep <- function(ctable, sampgroups, normmethod="TMM"){
	### BEGIN edgeR methods ###
	
	# create DGEList object
	suppressMessages( dge <- DGEList(counts=ctable, group=sampgroups) )
	dge <- calcNormFactors(dge, method="TMM")
	return(dge)
}

## ###########################################
## fit functions

## general function
seqfit <- function(seqobj, compgroups=NULL){
	
	result <- edgeRfit(seqobj,compgroups=compgroups)
		
    result$logFC <- ifelse( is.infinite(result$logFC), sign(result$logFC), ifelse( is.na(result$logFC), 0, result$logFC ) )
	result$p.value <- ifelse( is.na(result$p.value) | is.infinite(result$p.value), 1, result$p.value )
	result$FDR <- p.adjust(result$p.value, method="fdr")

	## make output names more similar to the DESeq component
	colnames(result)[1] <- c("ID")
  
	return(result)
}


edgeRfit <- function(dge, compgroups=NULL) {

  # estimate dispersions
  dge <- estimateCommonDisp(dge)
  dge <- estimateTagwiseDisp(dge)
  
  # test the results
  et <- exactTest(dge, pair=compgroups)
 
  # get all of the genes selected by the exactTest
  top <- topTags(et, n=nrow(dge$counts))$table
  
  # get total genes
  total.genes <- as.numeric(length(rownames(top)))

  # significant DEGs
  sig.DGEs <- rownames(top)[top$FDR < FDR.cutoff & abs(top$logFC) > logFC.cutoff]
  sig.upDGEs <- rownames(top)[top$FDR < FDR.cutoff & top$logFC > logFC.cutoff]
  sig.downDGEs <- rownames(top)[top$FDR < FDR.cutoff & top$logFC < -logFC.cutoff]
  
  # counts
  total.sigDGEs <- as.numeric(length(sig.DGEs))
  total.upDGEs <- length(sig.upDGEs)
  total.downDGEs <- length(sig.downDGEs)

    cutoffs.str <- paste("(FDR: ", FDR.cutoff, ", logFC: ", logFC.cutoff, ")", sep="")

  # print out totals
  cat("Summary:\n---------------------------------------------------\n")
  cat( paste("* Total Number of Genes: ", total.genes, "\n", sep="") )
  cat( paste("* Total Number of Significant Genes ", cutoffs.str, ": ", total.sigDGEs, "\n", sep="") )
  cat( paste("* Number of UP-regulated Significant Genes ", cutoffs.str, ": ", total.upDGEs, "\n", sep="") )	
  cat( paste("* Number of DOWN-regulated Significant Genes ", cutoffs.str, ": ", total.downDGEs, "\n\n", sep="") )	

  result <- et$table
  result$gene <- rownames(result)

  # comparison groups
  ugroups <- et$comparison
  
  ret <- result[,c("gene","logFC","PValue")]
  detags <- rownames(top)
  # ret[,4:5] <- dge$pseudo.counts
  ret <- cbind(ret, cpm(dge)[detags,])

  # assign headers
  gr <- paste(ugroups[2], ugroups[1], sep=" vs ")
  colnames(ret) <- c("gene", "logFC", "p.value", colnames(cpm(dge)[detags,]))

  # open PDF file for writing
  pdf.file <- paste(args[2], "/", ugroups[2], "_vs_", ugroups[1], ".FDR-", FDR.cutoff, "._MAplot.pdf", sep="")
  pdf(pdf.file)

  plot.title <- paste("DEG MA Plot ", cutoffs.str, "\n", sep="")

  plotSmear(dge, de.tags=sig.DGEs, main=plot.title, pair=compgroups, pch=16, cex=0.3, xlab="Log Read-Counts", ylab="Log Fold-Change")
  abline( h=c(-2,2) , col="dodgerblue" )
  
  total.txt <- paste("TOTAL: ", total.sigDGEs, "")
  up.txt <- paste("UP: ", total.upDGEs, sep="")
  down.txt <- paste("DOWN: ", total.downDGEs, sep="")

  legend(14, 10, c(total.txt, up.txt, down.txt, "\n"), border="white", cex=0.7, box.col="white")

  dev.off() 
  return(ret)
}


#### MAIN ####
args <- commandArgs(TRUE)

sample.table <- read.table(args[1], header=TRUE, stringsAsFactors=FALSE)

# input FDR cutoff as argument (default: 0.05)
FDR.cutoff = as.numeric(args[3])
if( is.na(FDR.cutoff) ) {
	FDR.cutoff = 0.05
}

logFC.cutoff = as.numeric(args[4])
if( is.na(logFC.cutoff) ) {
	logFC.cutoff = 1.0
}

cat("DEG Filtering Cutoffs:\n---------------------------------------------------\n")
cat(paste("* FDR: ", FDR.cutoff, "\n", sep=""))
cat(paste("* logFC: ", logFC.cutoff, "\n\n", sep=""))

sample.groups <- sample.table$CID
compgroups <- unique(sample.groups)[1:2]

# load the data based on *sample.table*
counts.table <- loadCountsTable(sample.table,count.colnum=2)

# remove non-gene rows from counts table
idx <- which(rownames(counts.table)=="no_feature")
counts.table <- counts.table[1:(idx-1),]
rownames(counts.table) <- c(rownames(counts.table)[1:(idx-1)])
counts.table <- counts.table[which(rowSums(counts.table[,1:ncol(counts.table)])>0),]
out <- cbind(rownames(counts.table), counts.table)
colnames(out) = c("Gene.ID", colnames(counts.table))

output.filename.txt <- paste(args[2],"/edgeR_",compgroups[2],"_vs_",compgroups[1],".all_counts_noZero.txt",sep="")
write.table(out, file=output.filename.txt, row.names=FALSE, col.names=TRUE, quote=FALSE,sep="\t")

# counts.table <- read.delim(file=output.filename.txt, header=T, sep="\t" )
counts.table <- sapply(2:ncol(out), function(x) as.numeric(out[,x]))
# counts.table <- counts.table[,c(2:ncol(counts.table))]

rownames(counts.table) <- out[,1]
colnames(counts.table) <- colnames(out)[2:ncol(out)]
counts.table <- data.matrix(counts.table)
cat("\n* Results Snippet:\n")
print(head(counts.table))

# create the initial edgeR object and normalize
seqobj <- seqprep(counts.table, sample.groups, normmethod="TMM")

# do the model fit, dispersions, etc
ret <- seqfit(seqobj, compgroups=compgroups)

# cat("\n* Results Snippet:\n")
# print(head(ret))

# sort the results table
ret.sort <- ret[order(ret$FDR),]
# cat("\n* Results Snippet:\n")
# print(head(ret.sort))

ret.sort <- ret.sort[,c(1,4:(ncol(ret.sort) - 1),2,3,ncol(ret.sort))]
# cat("\n* Results Snippet:\n")
# print(head(ret.sort))

# write the results table
output.filename.txt <- paste(args[2],"/edgeR_",compgroups[2],"_vs_",compgroups[1],".de_genes.txt",sep="")

write.table(ret.sort, file=output.filename.txt, row.names=FALSE, quote=FALSE,sep="\t")
