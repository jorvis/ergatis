invisible(options(echo = TRUE))

## read in data
pangenome <- read.table("###input_file###", header=FALSE)
genome_count <- max(pangenome$V8)
genomes <- (pangenome$V9[1:genome_count])
print(genomes)
pangenome <- pangenome[ pangenome$V1 > 1, ]
attach(pangenome)


## Calculate the means
v4allmeans <- as.vector(tapply(V4,V1,FUN=mean))
print(v4allmeans)
v1allmeans <- as.vector(tapply(V1,V1,FUN=mean))
print(v1allmeans)

## Calculate the medians
v4allmedians <- as.vector(tapply(V4,V1,FUN=median))
print(v4allmedians)
v1allmedians <- as.vector(tapply(V1,V1,FUN=median))
print(v1allmedians)

## plot points from each new comparison genome in its own color 
row_count <- length(V1)
source_colors <- rainbow(genome_count)
p_color <- c()
for ( ii in c(1:row_count) ) {
    p_color[ii] <- source_colors[V8[ii]]
#    points(temp_v1, temp_v4, pch=17, col=p_color)
}
## end of color block

# Hard-coded start and end values for the iterations. Note that the set needs 
# to have at least 5 genomes in it.
start_count <- max(floor(genome_count/2 + 0.5),5)
end_count <- 4

# Loop through from the start count to the end count (descending)
for(gcnt in start_count:end_count) {

# Generate a subset of the data
pansub <- pangenome[ pangenome$V1 >= gcnt, ]

## Calculate the means
v4submeans <- as.vector(tapply(pansub$V4,pansub$V1,FUN=mean))
v1submeans <- as.vector(tapply(pansub$V1,pansub$V1,FUN=mean))

## Calculate the medians
v4submedians <- as.vector(tapply(pansub$V4,pansub$V1,FUN=median))
v1submedians <- as.vector(tapply(pansub$V1,pansub$V1,FUN=median))

## exponential model based on medianss
nlmodel_pot <- nls(v4submedians ~ th1* v1submedians^(- th2), data=pansub,
start=list(th1=476, th2=0.5))


# Open up the output file for the log graph
postscript(file=paste(sep="", "###output_path###new_genes_power_law2_medians_log_min_",gcnt,"_error_bar.ps"), width=11, height=8.5, paper='special')
layout(matrix(c(1,2),byrow=TRUE), heights=c(7.5,1))

# Draw the axis
plot(V1,V4, xlab="number of genomes", ylab="new genes", main=paste("###TITLE### new genes power law 2 log axis ",gcnt," or more genomes"), col=p_color, cex=0.5, log="xy", type="n")

superpose.eb <-
function (x,y, ...) {
    sum = summary(y)
    q1 <- quantile(y, names=FALSE)[2]
    print(q1)
    med <- quantile(y,names=FALSE)[3]
    q4 <- quantile(y, names=FALSE)[4]
    print(q4)
    print(med)
    print(x)
   arrows(as.integer(x), sum[[2]],as.integer(x) , sum[[5]], angle = 90, code = 3,
   length = 0.08, ...)
}

m <- tapply(pangenome$V4,pangenome$V1,c)

for( x in names(m) ) {
    superpose.eb(x, m[[x]])
}

# plot the medians
points(tapply(pangenome$V4,pangenome$V1,FUN=median)~tapply(pangenome$V1,pangenome$V1,FUN=median),pch=5,col='black')

# plot the means
#points(tapply(V4,V1,FUN=mean)~tapply(V1,V1,FUN=mean),pch=6,col='black')

# plot the regression
# plot the regression
x <- seq(par()$xaxp[1]-1,as.integer(0.5 + 10^par()$usr[[2]]))
lines(x, predict(nlmodel_pot,data.frame(v1submedians=x)), lwd=2, col="black")



expr_pot <- substitute(
                expression(y == th1 * x^(-th2)), 
                list(
                    th1=round(nlmodel_pot$m$getPars()[1], digit=2),
                    th1err = round(summary(nlmodel_pot)[10][[1]][3], digit=2),
                    th2=round(nlmodel_pot$m$getPars()[2], digit=2),
                    th2err = round(summary(nlmodel_pot)[10][[1]][4], digit=2)
                    )
                      )

par(mai=c(.2,0,0,0))
height<- (10^(par()$usr[4]) - 10^(par()$usr[3]))
width<- (10^(par()$usr[2]) - 10^(par()$usr[1]))
plot.new()
par(xpd=T)
legend("top", c(eval(expr_pot)), lwd=c(2,2), yjust=0.5,xjust=0)
#legend(10^(par()$usr[2])+(0.01*width),10^(par()$usr[3]) + height/2, c(eval(expr_pot)), lwd=c(2,2), yjust=0.5,xjust=0)

dev.off()
}
