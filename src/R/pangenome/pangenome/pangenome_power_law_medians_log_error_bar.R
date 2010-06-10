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
nlmodel_pot <- nls(v2allmedians ~ th1*v1allmedians^(th2), data=pangenome,
start=list(th1=2000, th2=0))
#summary(nlmodel_pot)

# Open up the output file for the log graph
postscript(file="###output_path###pangenome_power_law_medians_log_error_bar.ps", width=11, height=8.5, paper='special')
layout(matrix(c(1,2),byrow=TRUE), heights=c(7.5,1))

# Draw the axis
plot(V1,V2, ylim=c(min(V2),max(V2)), xlab="number of genomes", ylab="pan-genome size", main="###TITLE### pan-genome power law log axis", cex=0.5, log="xy",type="n")

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
lines(x, predict(nlmodel_pot, data.frame(v1allmedians=x)), lwd=2, col="red")


#abline(h=nlmodel_pot$m$getPars()[1], lty=2, lwd=2,col="red")
expr_pot <- substitute(
                expression(y == th1 %+-% th1err * x^(th2 %+-% th2err)),
                list(
                    th1=round(nlmodel_pot$m$getPars()[1], digit=2),
                    th1err = round(summary(nlmodel_pot)[10][[1]][3], digit=2),
                    th2=round(nlmodel_pot$m$getPars()[2], digit=2),
                    th2err = round(summary(nlmodel_pot)[10][[1]][4], digit=4)
                    )
                      )

par(mai=c(.2,0,0,0))
height<- (10^(par()$usr[4]) - 10^(par()$usr[3]))
width<- (10^(par()$usr[2]) - 10^(par()$usr[1]))
plot.new()
par(xpd=T)
legend("top", c(eval(expr_pot)), lwd=c(2,2), yjust=0.5,xjust=0)
#legend(10^(par()$usr[2])+(0.01*width),10^(par()$usr[3]) + height/2, c(eval(expr_pot)), lwd=c(2,2), yjust=0.5,xjust=0)
