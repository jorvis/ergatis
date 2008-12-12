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
v1means <- as.vector(tapply(V1,V1,FUN=mean))

## Calculate the medians
v2allmedians <- as.vector(tapply(V2,V1,FUN=median))
v1allmedians <- as.vector(tapply(V1,V1,FUN=median))

point_colors <- c(
                    "palegreen4",
                    "paleturquoise4",
                    "palevioletred4",
                    "peachpuff4",
                    "plum4",
                    "purple4",
                    "red4",
                    "salmon",
                    "yellow4",
                    "snow4",
                    "steelblue4",
                    "wheat3",
                    "yellowgreen",
                    "rosybrown3",
                    "orangered3",
                    "hotpink3",
                    "khaki4",
                    "orange3",
                    "yellow3",
                    "lawngreen"
                 )

# plot points from each new comparison genome in its own color 
row_count <- length(V1)
source_colors <- rainbow(genome_count)
p_color <- c()
for ( ii in c(1:row_count) ) {
    p_color[ii] <- source_colors[V8[ii]]
#    points(temp_v1, temp_v4, pch=17, col=p_color)
}
## end of color block

## exponential model based on medianss
nlmodel_exp <- nls(v2allmedians ~ th1 + th2* exp(-v1allmedians / th3), data=pangenome,
start=list(th1=33, th2=476, th3=1.5))
#summary(nlmodel_exp)

# Open up the output file for the log graph
postscript(file="###output_path###core_genes_exponential_medians_log.ps")

# Add some space on the right for the legend(s)
par(mar=par()$mar+c(0,0,0,14))

# Draw the axis
plot(V1,V2, xlab="number of genomes", ylab="new genes", main="###TITLE### core genes exponential log axis", cex=0.5, log="xy", col=p_color)

# plot the medians
points(tapply(pangenome$V2,pangenome$V1,FUN=median)~tapply(pangenome$V1,pangenome$V1,FUN=median),pch=5,col='black')

# plot the means
points(tapply(V2,V1,FUN=mean)~tapply(V1,V1,FUN=mean),pch=6,col='black')

# plot the regression
x <- seq(par()$xaxp[1]-1,par()$xaxp[2]+1)
lines(x, predict(nlmodel_exp, data.frame(v1allmedians=x)), lwd=2, col="black")
abline(h=nlmodel_exp$m$getPars()[1], lty=2, lwd=2,col="black")

expr_exp <- substitute(
                expression(y == th1 + th2 * italic(e)^(-x / th3)), 
                list(
                    th1 = round(nlmodel_exp$m$getPars()[1], digit=2),
                    th2 = round(nlmodel_exp$m$getPars()[2], digit=2),
                    th3 = round(nlmodel_exp$m$getPars()[3], digit=2)
                    )
                )

height<- (10^(par()$usr[4]) - 10^(par()$usr[3]))
width<- (10^(par()$usr[2]) - 10^(par()$usr[1]))
par(xpd=T)
legend(10^(par()$usr[2])+(0.01*width),10^(par()$usr[3]) + height/2, c(eval(expr_exp)), lwd=c(2,2), yjust=0.5,xjust=0)
