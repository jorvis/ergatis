invisible(options(echo = TRUE))

## read in data
pangenome <- read.table("###input_file###", header=FALSE)
#genome_count <- max(pangenome$V8)
#genomes <- (pangenome$V9[1:genome_count])
#print(genomes)
pangenome <- pangenome[ pangenome$V1 > 1, ]
attach(pangenome)


## Calculate the means
v2means <- as.vector(tapply(V2,V1,FUN=mean))
v1means <- as.vector(tapply(V1,V1,FUN=mean))

## Calculate the medians
v2allmedians <- as.vector(tapply(V2,V1,FUN=median))
v1allmedians <- as.vector(tapply(V1,V1,FUN=median))

## exponential model based on medianss
nlmodel_exp <- nls(v2allmedians ~ th1 + th2* exp(-v1allmedians / th3), data=pangenome,
start=list(th1=2000, th2=-200, th3=2))
#summary(nlmodel_exp)

# Open up the output file for the log graph
postscript(file="###output_path###pangenome_exponential_medians_log_error_bar_with_points.ps", width=11, height=8.5, paper='special')
layout(matrix(c(1,2),byrow=TRUE), heights=c(7.5,1))

# Draw the axis
plot(V1,V2, ylim=c(min(V2),nlmodel_exp$m$getPars()[1]), xlab="number of genomes", ylab="pan-genome size", main="###TITLE### pangenome exponential log axis", cex=0.5, log="xy")


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
points(tapply(pangenome$V2,pangenome$V1,FUN=median)~tapply(pangenome$V1,pangenome$V1,FUN=median),pch=5,col='red')

# plot the means
#points(tapply(V2,V1,FUN=mean)~tapply(V1,V1,FUN=mean),pch=6,col='red')

# plot the regression
x <- seq(par()$xaxp[1]-1,as.integer(1.0 + 10^par()$usr[[2]]))
lines(x, predict(nlmodel_exp, data.frame(v1allmedians=x)), lwd=2, col="red")
abline(h=nlmodel_exp$m$getPars()[1], lty=2, lwd=2,col="red")

expr_exp <- substitute(
                expression(y == th1 + th2 * italic(e)^(-x / th3)), 
                list(
                    th1 = round(nlmodel_exp$m$getPars()[1], digit=2),
                    th1err = round(summary(nlmodel_exp)[10][[1]][3], digit=2),
                    th2 = round(nlmodel_exp$m$getPars()[2], digit=2),
                    th2err = round(summary(nlmodel_exp)[10][[1]][4], digit=2),
                    th3 = round(nlmodel_exp$m$getPars()[3], digit=2),
                    th3err = round(summary(nlmodel_exp)[10][[1]][5], digit=2)
                    )
                )

par(mai=c(.2,0,0,0))
height<- (10^(par()$usr[4]) - 10^(par()$usr[3]))
width<- (10^(par()$usr[2]) - 10^(par()$usr[1]))
plot.new()
legend("top", c(eval(expr_exp)), lwd=c(2,2), yjust=0.5,xjust=0)
#legend(10^(par()$usr[2])+(0.01*width),10^(par()$usr[3]) + height/2, c(eval(expr_exp)), lwd=c(2,2), yjust=0.5,xjust=0)