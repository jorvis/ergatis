invisible(options(echo = TRUE))

## read in data
pangenome <- read.table("###input_file###", header=FALSE)
genome_count <- max(pangenome$V8)
genomes <- (pangenome$V9[1:genome_count])
print(genomes)
pangenome <- pangenome[ pangenome$V1 > 1, ]
attach(pangenome)


## Calculate the means
v2means <- as.vector(tapply(V2,V1,FUN=mean))
print(v2means)
v1means <- as.vector(tapply(V1,V1,FUN=mean))
print(v1means)

## Calculate the medians
v2allmedians <- as.vector(tapply(V2,V1,FUN=median))
print(v2allmedians)
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

## power law model based on medianss
nlmodel_pot <- nls(v2allmedians ~ th1 + th2* v1allmedians^(th3), data=pangenome, 
start=list(th1=4000, th2=-100, th3=.5))

# Open up the output file for the log graph
postscript(file="###output_path###core_genes_power_law_medians_log_error_bar.ps", width=11, height=8.5, paper='special'))
layout(matrix(c(1,2),byrow=TRUE), heights=c(7.5,1))

# Draw the axis
plot(V1,V2, xlab="number of genomes", ylab="new genes", main="###TITLE### core genes power law log axis", col=p_color, cex=0.5, log="xy", type="n")


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

m <- tapply(pangenome$V2,pangenome$V1,c)

for( x in names(m) ) {
    superpose.eb(x, m[[x]])
}

# plot the medians
points(tapply(pangenome$V2,pangenome$V1,FUN=median)~tapply(pangenome$V1,pangenome$V1,FUN=median),pch=5,col='black')

# plot the means
points(tapply(V2,V1,FUN=mean)~tapply(V1,V1,FUN=mean),pch=6,col='black')

# plot the regression
x <- seq(par()$xaxp[1]-1,as.integer(1.0 + 10^par()$usr[[2]]))
lines(x, predict(nlmodel_pot, data.frame(v1allmedians=x)), lwd=2, col="black")

expr_pot <- substitute(
                expression(y == th1err %+-% th1 + th2err %+-% th2 * x^(th3 %+-% th3err)), 
                list(
                    th1=round(nlmodel_pot$m$getPars()[1], digit=2),
                    th1err = round(summary(nlmodel_pot)[10][[1]][3], digit=2),
                    th2=round(nlmodel_pot$m$getPars()[2], digit=2),
                    th2err = round(summary(nlmodel_pot)[10][[1]][4], digit=2),
                    th3=round(nlmodel_pot$m$getPars()[3], digit=2),
                    th3err = round(summary(nlmodel_pot)[10][[1]][5], digit=2)
                    )
                )


par(mai=c(.2,0,0,0))
height<- (10^(par()$usr[4]) - 10^(par()$usr[3]))
width<- (10^(par()$usr[2]) - 10^(par()$usr[1]))
plot.new()
legend("top", c(eval(expr_pot)), lwd=c(2,2), yjust=0.5,xjust=0)
#legend(10^(par()$usr[2])+(0.01*width),10^(par()$usr[3]) + height/2, c(eval(expr_pot)), lwd=c(2,2), yjust=0.5,xjust=0)