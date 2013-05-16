#! /usr/local/packages/R-2.15.2/bin/Rscript --slave --vanilla

#
# Usage: deseq.R $sample_info $out_dir 
#
# tcreasy@som.umaryland.edu 
#
###############################################################################

##
## sample_info: absolute path of sample info file. Sample info file format: sample_name\tphenotype\tabsolute path of read count file, no header.
## out_dir: absolute path of output directory
## annotation_file: optional. Absolute path of annotation file, first column of annotation file must be gene ID used for read counting, tab-delimited.
## 
## tcreasy

rm(list=ls(all=T))
cat("\n***** Starting DESeq (v1.10.1) analysis ******\n\n")

# input arguments
args <- commandArgs(TRUE)
print(length(args))
stopifnot(length(args) > 1 )

# sample comparison input file
sample.info <- read.delim(args[1], header=FALSE, sep="\t")

# output directory
output.dir <- args[2]

# set the working directory
setwd(output.dir)

## combine single file into a data matrix and output counting statistics
data.tab <- read.delim( as.character(sample.info[1,3]), header=F, sep="\t" )
idx <- which(data.tab=="no_feature")
count.stat <- matrix(c("sum counted", "total", sum(data.tab[,2][1:(idx-1)]), sum(data.tab[,2])), nrow=2, ncol=2)
count.stat <- rbind(data.tab[idx:nrow(data.tab),], count.stat)
colnames(count.stat) <- c("Stat", as.character(sample.info[1,1]))
data.tab <- data.tab[1:(idx-1),]

for (i in 2:nrow(sample.info)) {
	d <- read.delim(as.character(sample.info[i,3]), header=F, sep="\t")
	tmp <- matrix(c(sum(d[,2][1:(idx-1)]), sum(d[,2])), nrow=2, ncol=1)
	tmp <- data.frame(x=c(d[idx:nrow(d),2], tmp))
	n <- colnames(count.stat)
	count.stat <- cbind(count.stat, tmp)
	colnames(count.stat) <-c (n, as.character(sample.info[i,1]))
	d <- d[1:(idx-1),]
	data.tab <- merge(data.tab, d, by=1)
}


colnames(data.tab) <- c("ID", as.character(sample.info[,1]))
write.table(data.tab, file.path("all_counts"), na="", quote=F, row.names=F, sep="\t")
write.table(count.stat, file.path("count_stat"), na="", quote=F, row.names=F, sep="\t")


# remove genes without reads for all samples
data.tab <- data.tab[which(rowSums(data.tab[,2:ncol(data.tab)])>0),]
write.table(data.tab, file.path("all_counts_noZero"), na="", quote=F, row.names=F, sep="\t")
cat(paste("\n* There are ", ncol(data.tab)-1, " samples and ", nrow(data.tab), " genes used for DE analysis\n", sep=""))


# import DESeq and gplots
suppressMessages( library("DESeq") )
suppressMessages( library("RColorBrewer") )
suppressMessages( library("gplots") )


# phenotypes to compare
pheno <- unique(as.character(sample.info[,2]))
cat(paste("\n* Phenotypes found: ", toString(pheno), "\n", sep=""))
d <- data.tab[,2:ncol(data.tab)]
rownames(d) <- data.tab$ID


# condition
condition <- factor(as.character(sample.info[,2]))
cds <- newCountDataSet(d, condition)

cds <- estimateSizeFactors(cds)

#
# Estimate Variance (i.e. dispersions)
#
cat("\n* Estimating dispersions...\n")

if (sum(ifelse((data.frame(table(as.character(sample.info[,2])))$Freq < 2), 1, 0)) > 0) {
	# WITH NO replicates
	#
	# Options: 
	#  - method = "blind"
	#  - sharingMode = "fit-only"
	#  - fitType = "local"
	#
	cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only", fitType="local")
} else {
	# WITH FEW replicates (2 or less replicates)
	#
	# Options: 
	#  - method = "pooled"
	#  - sharingMode = "fit-only"
	#  - fitType = "local"
	#
	cds <- estimateDispersions(cds, method="pooled", sharingMode="fit-only", fitType="local")
}

# WITH replicates (as least 3 each)
#
# Options: 
#  - method = "per-condition"
#  - sharingMode = "maximum"
#  - fitType = "parametric"
#
#cds <- estimateDispersions(cds, method="per-condition", sharingMode="maximum", fitType="parametric")

# create an output file name for the output PDF
pdf.name <- paste(pheno[1], "-", pheno[2], ".pdf", sep="")
pdf(pdf.name)


# variance testing
cat("\n* Estimating variance...\n")
vsd <- varianceStabilizingTransformation(cds)
#select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:20]

# color palette for plots
hmcol = colorRampPalette(brewer.pal(9, "RdBu"))(100)

# heatmap of variance stabilized transformed data
#var.title <- paste("Top 30 DEGs: ", pheno[1], " vs ", pheno[2], sep="")
#heatmap.2(exprs(vsd)[select,], col=hmcol, trace="none", main=var.title, margin=c(13,13), cexRow=0.8, cexCol=0.8, keysize=1.0)

# cat("\n* Results Snippet:\n")
# print(head(exprs(vsd)))

# output to tab file
out <- cbind(rownames(exprs(vsd)), exprs(vsd))
colnames(out) <- c("ID", colnames(exprs(vsd)))
write.table(out, file.path("all_counts_noZero_normalized"), na="", sep="\t", quote=F, row.names=F)

# Heatmap showing clustering of samples
dists = dist( t( exprs(vsd) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat)
dist.title <- paste("Sample Clustering", sep="")
heatmap.2(mat, col=rev(hmcol), trace="none", main=dist.title, margin=c(13,13), cexRow=0.8, cexCol=0.8, keysize=1.0)


for (k in 1:(length(pheno)-1)) {
	
	for (m in (k+1):length(pheno)) {
		
		cat(paste("\n* Running nbinomTest for: ", pheno[k], " vs ", pheno[m], "\n", sep=""))

        # write png image of variance estimation
	    #png(file=paste(pheno[k], "_vs_", pheno[m], ".variance_estimation.png", sep=""), width=720, height=720)
		#plotDispEsts(cds)
		#dev.off()

		# run the negative binomial DEG test
		res <- nbinomTest(cds, pheno[m], pheno[k])
		#res <- cbind(res, NA, NA)

		# order output by FDR
		res <- res[order(res$padj),]
		
		cat("\n* Results Snippet: res\n")
		print(head(res))
		
		# plot the results using FDR=0.05 as the cutoff
		ma.title <- paste("DEG MA Plot", " (FDR < 0.05)", sep="")
		plotMA(res, main=ma.title, cex=0.3, col=ifelse(res$padj>=0.05, "black", "red"), linecol="#ff000080", xlab="Mean of Normalized Counts", ylab=paste("LFC: ", pheno[k], " VS ", pheno[m], sep=""))

		# plot the results using abs(LFC)>=1.0 as the cutoff
		#plotMA(res, main="DEG MA Plot", cex=.2, col=ifelse( abs(res$log2FoldChange)>=1.0, "black", "red"), linecol="black", xlab="Mean of Normalized Counts", ylab=paste("LFC: ", pheno[k], " VS ", pheno[m], sep=""))
		
		
		# get read counts for each group for the top 30 most significant DEGs
		# order output by absolute LFC
		res <- res[order(-abs(res$log2FoldChange)),]
		sig.genes <- res[res$padj<=0.05,]
		print(dim(sig.genes))
		
		if(nrow(sig.genes) < 2) {
			sig.genes <- res[res$pval<=0.05,]
			print(dim(sig.genes))
		}
		
		if(nrow(sig.genes) > 30) {
			sig.genes <- sig.genes[1:30,]
			print(dim(sig.genes))
		}
		
		cat("\n* Results Snippet: sig.genes\n")
		print(head(sig.genes))
		
		read.counts <- cbind(as.numeric(sig.genes[,3]), as.numeric(sig.genes[,4]))

		colnames(read.counts) <- c(pheno[m], pheno[k])
		rownames(read.counts) <- c(sig.genes[,1])
		
		cat("\n* Results Snippet: read.counts.1\n")
		print(head(read.counts))

		# draw heatmap of normalized read counts for the significant genes of each sample
		sig.title <- paste("Top Significant DEGs", " (per condition)", sep="")		
		heatmap.2(read.counts, col=hmcol, trace="none", main=sig.title, margin=c(13,13), cexRow=0.8, cexCol=0.8, keysize=1.0)
		
		read.counts <- sig.genes[,c(1,6)]
		colnames(read.counts) <- c("ID", "LFC")
		read.counts <- merge(read.counts, out, by="ID", x.all=TRUE)
		read.counts <- read.counts[order(-abs(read.counts$LFC)),]
		
		cat("\n* Results Snippet: read.counts.2\n")
		print(head(read.counts))
		
		write.table(read.counts, file.path(paste(pheno[k], "_vs_", pheno[m], ".top30.counts.txt", sep="")), na="", quote=F, row.names=F, sep="\t")
		
		hmap <- read.delim(file.path(paste(pheno[k], "_vs_", pheno[m], ".top30.counts.txt", sep="")), header=T, sep="\t" )
		
		cat("\n* Results Snippet: hmap.1\n")
		print(head(hmap))
		
		hmap <- hmap[,c(3:ncol(hmap))]
		colnames(hmap) <- c(colnames(read.counts)[3:ncol(read.counts)])
		rownames(hmap) <- c(read.counts[,1])
		hmap <- data.matrix(hmap)
		
		cat("\n* Results Snippet: hmap.2\n")
		print(head(hmap))
		
		# draw heatmap of normalized read counts for the significant genes of each sample
		sig.title <- paste("Top Significant DEGs",  " (per sample)", sep="\n")		
		heatmap.2(hmap, col=hmcol, trace="none", main=sig.title, margin=c(13,13), cexRow=0.8, cexCol=0.8, keysize=1.0)
		
		# Change column names for clarity and brevity
		#colnames(res) <- c("Feature.ID", "Read.Count.All", paste("Read.Count.", pheno[m], sep=""), paste("Read.Count.", pheno[k], sep=""), "FC", paste("LFC(", pheno[k], "/", pheno[m], ")", sep=""), "p.Value", "FDR", "NA", "NA")
		colnames(res) <- c("Feature.ID", "Read.Count.All", paste("Read.Count.", pheno[m], sep=""), paste("Read.Count.", pheno[k], sep=""), "FC", paste("LFC(", pheno[k], "/", pheno[m], ")", sep=""), "p.Value", "FDR")	
		# write data to tsv file
		write.table(res, file.path(paste(pheno[k], "_vs_", pheno[m], ".de_genes.txt", sep="")), na="", quote=F, row.names=F, sep="\t")
	}
}

data.tab <- NULL
d <- NULL

cat("\n\n* Garbage Collection Information\n\n")
gc()

cat("\n\n***** DEG Analysis Complete *****\n\n")
