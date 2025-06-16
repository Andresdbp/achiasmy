## CDS and Equations

gff <- read.csv2("data/gencode.v47.annotation.gff3", skip = 7, header = F, sep = "\t")

# Gene sites
gff_gene <- gff[gff$V3=="gene",]
gff_gene_X <- gff_gene[gff_gene$V1=="chrX",]
size_genes_x <- abs(
  gff_gene_X$V5 - gff_gene_X$V4
)
sum(size_genes_x)

# CDS sites

gff_cds <- gff[gff$V3=="CDS",]
gff_cds_X <- gff_cds[gff_cds$V1=="chrX",]
size_cds_x <- abs(
  gff_cds_X$V5 - gff_cds_X$V4
)
sum(size_cds_x)





me <- 0.03
u <- 1e-09
acs <- 2e+07
gen <- 1000

fit <- c()
for (t in 1:gen) {
  ml <- (1 - me)^(u * acs * t)
  fit[t] <- ml
}
plot(fit, type = "l", ylim = c(0,1), col = "blue",
     ylab = "Fitness", xlab = "Generations", lwd = 1.5)

fit <- c()
ml <- 0
for (t in 1:gen) {
  ml <- ml + me*u*acs
  w <- 1-ml
  fit[t] <- w
}
lines(fit, col = "red", lwd = 1.5)

points(0, 0.2, col = "red", pch = "-", cex = 2)
text(0, 0.2, pos = 4, labels = "Linear")
points(0, 0.1, col = "blue", pch = "-", cex = 2)
text(0, 0.1, pos = 4, labels = "Non-linear")



ml <- c(0)
for (i in 1:1000) {
  ml[i] <- (u * acs * i)
}
fitness_decline_log <- (1 - me) ^ ml
plot(fitness_decline_log, type = "l", ylim = c(0,1), col = "red")

ml <- c(0)
for (i in 1:1000) {
  ml[i+1] <- ml[i] + (u * acs * me)
}
fitness_decline_mul <- (1 - ml)
lines(fitness_decline_mul, type = "l", col = "blue")

# K-class fitness function

f <- function(k) {
  (1+0.1)*((1-0.03)^k)
}
k <- c(0:99)
y <- f(k)
plot(k, y, type = "l", 
     main = "Plot of (1-S)^k", 
     xlab = "x", ylab = "y", ylim = c(0,1.5))


par(mfrow=c(1,1))

# Define the function
f <- function(x) {
  sqrt( 2 * ( -((1e-8) * (x)) * log(1 - 0.03) ) * log(1 / 0.0001) )
}

# Create a sequence of x values from 0 to 0.99
x <- seq(0, 5e+06, length.out = 1000)
y <- f(x)

image(resY, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "Coding Sites Transitioning to Achiasmy (Mb)",
      main = "Frequency on Y",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
lines(x = x/(5e+06), y = y/(0.5), col = "red",
      lwd = 2)

# Plot the function
plot(x, y, type = "l", 
     main = "Plot of y = sqrt(2*(-((1e-09)*(2e+06))*ln(1-x))*ln(0.99/0.01))", 
     xlab = "x", ylab = "y", ylim = c(0,0.5))

