#! /usr/local/packages/R-2.14.1/bin/Rscript --vanilla


###############################################################################

## arguments:
## 1: sample table
## 2: file output path


##library(DESeq)
##library(baySeq)
library(edgeR)
##library(DSS)


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


## ###########################################
## object prep functions

## general function
seqprep <- function(ctable,
                    sampgroups,
                    use="deseq",
                    normmethod="none"){ ## "TMM", "deseq", or anything else for no norm factors
  if( use=="deseq" ){
    seqobj <- deseqprep(ctable,sampgroups,normmethod=normmethod)
  }else if( use %in% c("edgeR","ttest") ){
    seqobj <- edgeRprep(ctable,sampgroups,normmethod=normmethod)
  }else if( use=="bayseq" ){
    seqobj <- bayseqprep(ctable,sampgroups,normmethod=normmethod)
  }else if( use=="DSS" ){
    seqobj <- DSSprep(ctable,sampgroups,normmethod=normmethod)
  }
  return(seqobj)
}


edgeRprep <- function(ctable,
                      sampgroups,
                      normmethod="TMM"){
  ## #################
  ## BEGIN edgeR methods
  
  ## #####
  ## create object
  dge <- DGEList(counts=ctable, group=sampgroups)
  
  ## #####
  ## estimate norm.factors
  if( normmethod=="TMM" ){
    dge <- calcNormFactors(dge,method="TMM")
  }else if(normmethod=="deseq"){
    cds <- newCountDataSet( ctable, sampgroups )
    dge$samples$norm.factors <- sizeFactors(estimateSizeFactors(cds))
  }else{
    ## this probably isn't right, but gives reproducible rankings
    dge$samples$norm.factors <- apply(dge$counts,2,sum) / mean(apply(dge$counts,2,sum))
  }
  ##effective.library.size <- dge$samples$lib.size * dge$samples$norm.factor

  return(dge)
}

deseqprep <- function(ctable,
                      sampgroups,
                      normmethod="TMM"){
  ## create object
  cds <- newCountDataSet( ctable, sampgroups )

  ## #########
  ## estimate library sizes
  if( normmethod=="TMM" ){
    ## to use the TMM sizes from baySeq
    CD <- new("countData",
              data=as.matrix(ctable), #simData,
              replicates=sampgroups) 
                                        #groups=groups)
    sizeFactors(cds) <- getLibsizes(CD,estimationType="edgeR")
    ## ############
  }else if(normmethod=="deseq"){
    cds <- estimateSizeFactors(cds)
  }else{
    sizeFactors(cds) <- apply(counts(cds),2,sum) / mean(apply(counts(cds),2,sum))
  }

  return(cds)
}


bayseqprep <- function(ctable,
                       sampgroups,
                       normmethod="TMM"){
  repls <- sampgroups
  
  grps <- list(NDE=rep(1,length(repls)),
                 DE=sapply(repls,function(x)
                   which( unique(repls)==x)) )
  
  ## construct ’countData’ object
  CD <- new("countData",
            data=as.matrix(ctable), #simData,
            replicates=repls,
            groups=grps)
  
  
  ## ###########################
  ## estimating library sizes
  if( normmethod=="TMM" ){
    ##CD@libsizes <- getLibsizes(CD)
    ##CD@libsizes <- getLibsizes(CD,estimationType="total")
    ## CD@libsizes <- getLibsizes(CD,estimationType="quantile",
    ##                            quantile=0.75) ## this is the default
    CD@libsizes <- getLibsizes(CD,estimationType="edgeR") *
      mean(getLibsizes(CD,estimationType="total"))
  }else if(normmethod=="deseq"){
    cds <- newCountDataSet( CD@data, CD@replicates ) ## using the DESeq method
    CD@libsizes <- sizeFactors(estimateSizeFactors(cds)) *
      mean(getLibsizes(CD,estimationType="total"))
  }else{
    CD@libsizes <- getLibsizes(CD,estimationType="total")
  }

  return(CD)
}


DSSprep <- function(ctable,
                    sampgroups,
                    normmethod="TMM"){

  cmat <- as.matrix(ctable)
  colnames(cmat) <- NULL
  
  seqData <- newSeqCountSet(cmat, sampgroups)

  seqData <- estNormFactors(seqData,method="total")

  ## ###########################
  ## estimating library sizes
  if( normmethod=="TMM" ){
    ## to use the TMM sizes from baySeq
    CD <- new("countData",
              data=as.matrix(ctable),
              replicates=sampgroups) 
    seqData@normalizationFactor <- getLibsizes(CD,estimationType="edgeR")
  }else if(normmethod=="deseq"){
    cds <- newCountDataSet( ctable, sampgroups )
    seqData@normalizationFactor <- sizeFactors(estimateSizeFactors(cds))
  }else{
    seqData <- estNormFactors(seqData,method="total")
  }

  return(seqData)
}




## ###########################################
## fit functions

## general function
seqfit <- function(seqobj,use="deseq",compgroups=NULL){
  if( class(seqobj)=="DGEList" ){
    if( use=="edgeR" ){
      result <- edgeRfit(seqobj,compgroups=compgroups)
    }else if( use=="ttest" ){
      result <- ttestfit(seqobj,compgroups=compgroups)
    }
  }else if( use=="deseq" && class(seqobj)=="CountDataSet" ){
    result <- deseqfit(seqobj,compgroups=compgroups,method="pooled",sharingMode="fit-only")
  }else if( use=="bayseq" && class(seqobj)=="countData" ){
    result <- bayseqfit(seqobj,compgroups=compgroups,estimation="QL")
  }else if( use=="DSS" && class(seqobj)=="SeqCountSet" ){
    result <- DSSfit(seqobj,compgroups=compgroups)
  }

  
  
  ## result$logFC <- ifelse( is.na(result$logFC) |
  ##                        is.infinite(result$logFC), 0, result$logFC )
  ## result$p.value <- ifelse( is.na(result$p.value) |
  ##                          is.infinite(result$p.value), 1, result$p.value )

  result$logFC <- ifelse( is.infinite(result$logFC),
                    sign(result$logFC), ifelse( is.na(result$logFC), 0, result$logFC ))
  result$p.value <- ifelse( is.na(result$p.value) |
                      is.infinite(result$p.value), 1, result$p.value )
  result$FDR <- p.adjust(result$p.value,method="fdr")
  ## result$signlogpval <- log10(result$p.value) * -sign(result$logFC)

  ## make output names more similar to the DESeq component
  colnames(result)[1] <- c("ID")
  
  return(result)
}


edgeRfit <- function(dge,compgroups=NULL){
  ## #####
  ## estimate dispersions
  dge <- estimateCommonDisp(dge)
  dge <- estimateTagwiseDisp(dge)

  ## #####
  ## test the results
  et <- exactTest(dge,pair=compgroups)
 
  ##top <- topTags(et, n=nrow(dge$counts))

  result <- et$table
  ##result$FDR <- p.adjust(result$PValue)
  result$gene <- rownames(result)
  ## result <- result[, c(ncol(result), 1:(ncol(result) - 1))]
  
  ##result$DE.call.1 <- rank(result$PValue) <= 0.1 * nrow(result)
  ##result$DE.call.1 <- result$PValue <= 0.1
  ##sum(result$DE.call.1)

  ##ugroups <- levels(dge@.Data[[1]]$group)
  ugroups <- et$comparison
  
  ret <- result[,c("gene","logFC","PValue")]
  ret[,4:5] <- dge$conc$conc.group
  colnames(ret) <- c("gene","logFC","p.value",paste("avg",ugroups,sep="_"))
  return(ret)
}


deseqfit <- function(cds,method="pooled",
                     sharingMode="fit-only",
                     compgroups=NULL){
  
  ## #####
  ## estimate dispersions
 
  ##cds.disp <- estimateDispersions( cds )
  cds.disp <- estimateDispersions( cds,
                                  method=method, #c( # "pooled"
                                  ##"pooled-CR" # implemented by edgeR, 2010
                                  ## "per-condition" # used to be called "normal"
                                  ## "blind" # ignores conditions
                                  ##  ),
                                  sharingMode=sharingMode, ##c( #"maximum"
                                  ##"fit-only"
                                  ##"gene-est-only"
                                  ##  ),
                                  fitType=c( #"parametric"
                                    "local"
                                    )
                                  )
  
  ## estimateDispersions( object,
  ##                     method = c( "pooled", "pooled-CR", "per-condition", "blind" ),
  ##                     sharingMode = c( "maximum", "fit-only", "gene-est-only" ),
  ##                     fitType = c("parametric", "local"),
  ##                     locfit_extra_args=list(), lp_extra_args=list(),
  ##                     modelFrame = NULL, modelFormula = count ~ condition, ... )
  
  
  ## str( fitInfo(cds,name="ERY") )
  ## str( fitInfo(cds,name="uninduced") )
  
  ## fData(cds) these are the estimated, fitted dispersion values
  
  ## #####
  ## test the results
  ##result <- nbinomTest(cds.disp,"uninduced","ERY")
  result <- nbinomTest(cds.disp,
                       compgroups[1],
                       compgroups[2] )
  
  ##head(result[order(result$padj),])

  ugroups <- compgroups

  ret <- result[,c("id","foldChange","pval","baseMeanA","baseMeanB")]
  colnames(ret) <- c("gene","logFC","p.value",paste("avg",ugroups,sep="_"))
  ret[,"logFC"] <- log(result[,"foldChange"])
  
  return(ret)
}


bayseqfit <- function(CD,estimation="QL",
                      compgroups=NULL){

  ## estimate prior distributions on ’countData’ object using negative binomial
  ## method. Other methods are available - see getPriors
  CDPriors <- getPriors.NB(CD,
                           ##estimation=estimation,
                           estimation="QL", # default
                           ##estimation="ML", 
                           ## estimation="edgeR", 
                           cl = NULL)
  
  ## estimate posterior likelihoods for each row of data belonging to each hypothesis
  CDPost <- getLikelihoods(CDPriors, cl = NULL)
  
  ## ## display the rows of data showing greatest association with the second
  ## ## hypothesis (differential expression)
  ## topCounts(CDPost, group = "DE", number = 1 )

  ugroups <- unique(names(CDPost@groups$DE))

  ret <- data.frame(gene=rownames(CDPost@data),
                    logFC=getRawLogFoldChange(CDPost@data,
                      compgroups),
                    p.value=exp(CDPost@posteriors[,"NDE"]),
                    row.names=NULL)
  ret[,4:5] <- getMeanCounts(CDPost@data,
                             names(CDPost@groups$DE))
names(ret)[4:5] <- paste("avg",ugroups,sep="_")

  return(ret)
}



DSSfit <- function(seqData,compgroups=NULL){

  seqData <- estDispersion(seqData)

  ugroups <- levels( seqData@phenoData@data$designs )
  result <- waldTest(seqData, compgroups[1], compgroups[2])

  result$gene <- rownames(result)
  sampmeans <- result[,c("muA","muB")]
  
  ret <- result[,c("gene","lfc","pval")]
  ret[,4:5] <- sampmeans
  colnames(ret) <- c("gene","logFC","p.value",paste("avg",ugroups,sep="_"))
  
  return(ret)
}


ttestfit <- function(dge,compgroups=NULL){
  ## dge is an object from edgeR

  ugroups <- compgroups ##levels(dge$samples$group)

  ## using limma
  xpr <- NULL
  ## xpr$E <- normalizeBetweenArrays(as.matrix(expertab.list),
  ##                                              method="quantile")
  xpr$E <- log(dge$counts+1)
  xpr$genes <- data.frame(GeneName=rownames(dge$counts),
                          stringsAsFactors=FALSE)
  
  xpr$targets <- data.frame(groups=dge$samples$group,
                            stringsAsFactors=FALSE )
  colnames(xpr$targets) <- "groups"
  xpr <- new("EList",xpr)
  design <- model.matrix( ~ 0 + xpr$targets$group )
  colnames(design) <- gsub("xpr$targets$group","",
                           colnames(design),fixed=TRUE)
  
  fit <- lmFit(xpr, design)

  contrasts <- paste(compgroups,sep="",collapse="-")
  contrmat <- makeContrasts( contrasts=contrasts,levels=design )


  ## for genes
  fit2 <- contrasts.fit(fit, contrmat)
  fit3 <- eBayes(fit2)

  ttab <- topTable(fit3, number=100000, adjust.method="BH",
                   sort.by="p" )

  ret <- data.frame(gene=ttab$GeneName,
                    logFC=ttab$logFC,
                    p.value=ttab$P.Value,
                    row.names=NULL)


  ret[,4:5] <- getMeanCounts(dge$counts,
                      dge$samples$group)
  names(ret)[4:5] <- paste("avg",ugroups,sep="_")

  return(ret)
}



## #####################
## utility functions

getRawLogFoldChange <- function(ctable,sampgroups){
  xtab <- getMeanCounts(ctable,sampgroups)
  ret <- log(xtab[,2])-log(xtab[,1])
  return(ret)
}

getMeanCounts <- function(ctable,sampgroups){
  ugrps <- unique(sampgroups)
  ret <- t( apply(ctable,1,function(xrow){
    logrow <- log(xrow+1)
    ret2 <- c(exp(mean(logrow[sampgroups==ugrps[1]])),
              exp(mean(logrow[sampgroups==ugrps[2]])))-1
    return(ret2)
  }) )
  return(ret)
}



##
args <- commandArgs(TRUE)


sample.table <- read.table(args[1],
                           header=TRUE,
                           stringsAsFactors=FALSE)



## load the data based on *sample.table*
counts.table <- loadCountsTable(sample.table,count.colnum=2)
sample.groups <- sample.table$CID
compgroups <- unique(sample.groups)[1:2]


## create the initial edgeR object and normalize
seqobj <- seqprep(counts.table,
                  sample.groups,
                  use="edgeR",
                  normmethod="TMM")

## do the model fit, dispersions, etc
ret <- seqfit(seqobj,
              use="edgeR",
              compgroups=compgroups)

## sort the results table
ret.sort <- ret[order(ret$FDR),]


## write the results table
output.filename.txt <- paste(args[2],"/edgeR_",compgroups[2],"_vs_",compgroups[1],".de_genes.txt",sep="")

write.table(ret.sort,
            file=output.filename.txt,
            row.names=FALSE,
            quote=FALSE,sep="\t")


## create an MA plot
output.filename.pdf <- paste(args[2],"/edgeR_",compgroups[2],"_vs_",compgroups[1],"_MAplot.pdf",sep="")

##x11()
pdf(output.filename.pdf)
maPlot(ret[,4],ret[,5],pch="x",
       main=paste("MA plot: ",compgroups[2],"_vs_",compgroups[1],sep=""))
dev.off()
