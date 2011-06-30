## Usage: /usr/local/bin/R --slave --vanilla --args $sample_info $out_dir $annotation_file < run_deseq.R
##
## sample_info: absolute path of sample info file. Sample info file format: sample_name\tphenotype\tabsolute path of read count file, no header.
## out_dir: absolute path of output directory
## annotation_file: optional. Absolute path of annotation file, first column of annotation file must be gene ID used for read counting, tab-delimited.
## 
## by ysun

rm(list=ls(all=T))
print("Start DE analysis...")

arg<-commandArgs()
x<-arg[5]
if (is.na(x)) {stop("No sample info file!")}
sample.info<-read.delim(x, header=F, sep="\t")
print(paste("Sample info file: ", x, sep=""))

x<-arg[6]
if (is.na(x)) {
	stop("No output directory!")
} else if (!file.exists(x)) {stop("Ouput directory doesn't exist!")}
print(paste("Output directory: ", x, sep=""))
setwd(x)

x<-arg[7]
noAnn<-is.na(x)
if (noAnn) {
	print("No annotation file.")
} else {
	print(paste("Annotation file: ", x, sep=""))
	ann<-read.delim(x, header=T, sep="\t")
}

## combine single file into a data matrix and output counting statistics

data.tab<-read.delim(as.character(sample.info[1,3]), header=F, sep="\t")
idx<-which(data.tab=="no_feature")

count.stat<-matrix(c("sum counted", "total", sum(data.tab[,2][1:(idx-1)]), sum(data.tab[,2])), nrow=2, ncol=2)
count.stat<-rbind(data.tab[idx:nrow(data.tab),], count.stat)
colnames(count.stat)<-c("Stat", as.character(sample.info[1,1]))
data.tab<-data.tab[1:(idx-1),]

for (i in 2:nrow(sample.info)){
	d<-read.delim(as.character(sample.info[i,3]), header=F, sep="\t")
	tmp<-matrix(c(sum(d[,2][1:(idx-1)]), sum(d[,2])), nrow=2, ncol=1)
	tmp<-data.frame(x=c(d[idx:nrow(d),2], tmp))
	n<-colnames(count.stat)
	count.stat<-cbind(count.stat, tmp)
	colnames(count.stat)<-c(n, as.character(sample.info[i,1]))
	d<-d[1:(idx-1),]
	data.tab<-merge(data.tab, d, by=1)
}

colnames(data.tab)<-c("ID", as.character(sample.info[,1]))
write.table(data.tab, file.path("all_counts"), na="", quote=F, row.names=F, sep="\t")
write.table(count.stat, file.path("count_stat"), na="", quote=F, row.names=F, sep="\t")

## remove gene without read for all samples

data.tab<-data.tab[which(rowSums(data.tab[,2:ncol(data.tab)])>0),]
write.table(data.tab, file.path("all_counts_noZero"), na="", quote=F, row.names=F, sep="\t")
print(paste("There are ", ncol(data.tab)-1, " samples and ", nrow(data.tab), " genes used for DE analysis.", sep=""))

## DESeq procedure

if(!require(DESeq)) {stop("Can't load DESeq library!")}

pheno<-unique(as.character(sample.info[,2]))
print(paste("Phenotypes found: ", toString(pheno), sep=""))

d<-data.tab[,2:ncol(data.tab)]
rownames(d)<-data.tab$ID
cond<-factor(as.character(sample.info[,2]))
cds<-newCountDataSet(d, cond)
cds<-estimateSizeFactors(cds)

if (ncol(d)<3) {cds<-estimateVarianceFunctions(cds, method="blind")} else {cds<-estimateVarianceFunctions(cds)}

tiff(file="variance_estimation.tiff", compression=c("none"), width=720, height=720)
scvPlot(cds, ylim=c(0,2))
dev.off()

vsd<-getVarianceStabilizedData(cds)
dists<-dist(t(vsd))
tiff(file="sample_clustering.tiff", compression=c("none"), width=720, height=720)
heatmap(as.matrix(dists), symm=T, margins=c(10,10))
dev.off()
out<-cbind(rownames(vsd), vsd)
colnames(out)<-c("ID", colnames(vsd))
write.table(out, file.path("all_counts_noZero_normalized"), na="", sep="\t", quote=F, row.names=F)

for (k in 1:(length(pheno)-1)){
	for (m in (k+1):length(pheno)){
		print(paste("Running ", pheno[k], " vs ", pheno[m], sep=""))
		res<-nbinomTest(cds, pheno[m], pheno[k])

		tiff(file=paste(pheno[k], "_vs_", pheno[m], ".tiff", sep=""), compression=c("none"), width=720, height=720)
		plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1, main=paste(pheno[k], " vs ", pheno[m], sep=""), ylab=paste("log2(", pheno[k], "/", pheno[m], ")", sep=""), xlab="Average counts across 2 conditions", cex.axis=1.2, cex.lab=1.5, cex.main=1.5)		
		dev.off()
		
		if(!noAnn){
			res<-merge(ann, res, by=1, all.y=T)
			colnames(res)<-c(colnames(ann), "Avg.normalized.read.count.all", paste("Avg.normalized.read.count.", pheno[m], sep=""), paste("Avg.normalized.read.count.", pheno[k], sep=""), "Fold.change", paste("log2(", pheno[k], "/", pheno[m], ")", sep=""), "Raw.p.value", "FDR", paste("resVar.", pheno[m], sep=""), paste("resVar.", pheno[k], sep=""))			
		}
		else {
			colnames(res)<-c("ID", "Avg.normalized.read.count.all", paste("Avg.normalized.read.count.", pheno[m], sep=""), paste("Avg.normalized.read.count.", pheno[k], sep=""), "Fold.change", paste("log2(", pheno[k], "/", pheno[m], ")", sep=""), "Raw.p.value", "FDR", paste("resVar.", pheno[m], sep=""), paste("resVar.", pheno[k], sep=""))
		}
		res<-res[order(res$FDR),]
		write.table(res, file.path(paste(pheno[k], "_vs_", pheno[m], sep="")), na="", quote=F, row.names=F, sep="\t")
	}
}

data.tab<-NULL
d<-NULL
gc()
print("DE analysis succeed!")
