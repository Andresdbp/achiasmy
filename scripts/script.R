# Andres Barboza
# Model of achiasmatic meiosis with a sexually antagonistic allele
# andresdbp@tamu.edu

library(viridis)
library(fields)

setwd("~/Documents/GitHub/achiasmy")

# Variables needed
## h - dominance factor of the male benefit allele
## s - selection coefficient of the SAL
## r1 - recombination frequency between the SDR - SAL
## r2 - recombination freq between SAL - Achiasmy locus
## gen - generations
## me - mutation effect
## u - mutation rate

## Main Function
simLife <- function(h = 1, s = 0, r1 = 0, r2 = 0,
                    gen = 1000, me = 0.001, u = 0.001,
                    aneu=0, option="Y"){
  #TODO
  # mutational load
  ml <- 0
  
  # fitness matrices
  wf <- c(1, 1/(1+h*s), 1/(1+s))
  wm <- matrix(c(1-aneu, 1-ml, (1+h*s)*(1-aneu), (1+h*s)*(1-ml), (1+s)*(1-aneu), (1+s)*(1-ml)), 2, 3)
  
  # Initial gamete frequencies
  if (option == "Y") {
    # eggs
    xfr.e <- 0.5
    xfa.e <- 0
    xmr.e <- 0.5
    xma.e <- 0
    # sperm
    xfr.s <- 0.25
    xfa.s <- 0
    xmr.s <- 0.25
    xma.s <- 0
    yfr.s <- 0.24995
    yfa.s <- 0.00005
    ymr.s <- 0.24995
    yma.s <- 0.00005
  }
  if (option == "X") {
    # eggs
    xfr.e <- 0.499925
    xfa.e <- 0.000075
    xmr.e <- 0.499925
    xma.e <- 0.000075
    # sperm
    xfr.s <- 0.249925
    xfa.s <- 0.000025
    xmr.s <- 0.249925
    xma.s <- 0.000025
    yfr.s <- 0.25
    yfa.s <- 0
    ymr.s <- 0.25
    yma.s <- 0
  }
  if (option == "A") {
    # eggs
    xfr.e <- 0.49995
    xfa.e <- 0.00005
    xmr.e <- 0.49995
    xma.e <- 0.00005
    # sperm
    xfr.s <- 0.249975
    xfa.s <- 0.000025
    xmr.s <- 0.249975
    xma.s <- 0.000025
    yfr.s <- 0.249975
    yfa.s <- 0.000025
    ymr.s <- 0.249975
    yma.s <- 0.000025
  }

  results <- matrix(NA,gen,8)
  colnames(results) <- c("Xfr", "Xfa", "Xmr", "Xma","Yfr", "Yfa", "Ymr", "Yma")
  results[1,1] <- (xfr.e+xfr.s)/1.5
  results[1,2] <- (xfa.e+xfa.s)/1.5
  results[1,3] <- (xmr.e+xmr.s)/1.5
  results[1,4] <- (xma.e+xma.s)/1.5
  results[1,5] <- yfr.s/0.5
  results[1,6] <- yfa.s/0.5
  results[1,7] <- ymr.s/0.5
  results[1,8] <- yma.s/0.5
  
  for (i in 2:gen) {
    if (((option == "Y") & (sum(results[i-1,5:8] <= 0.99)>3)) | ((option == "X") & (sum(results[i-1,1:4] <= 0.99)>3)) | (option == "A")) {
      #mean fitnesses (weighted sums, which is a mean)
      wbarf <- xfr.e*xfr.s*wf[1] + xfr.e*xfa.s*wf[1] + xfr.e*xmr.s*wf[2] + xfr.e*xma.s*wf[2] + xfa.e*xfr.s*wf[1] + xfa.e*xfa.s*wf[1] + xfa.e*xmr.s*wf[2] + xfa.e*xma.s*wf[2] + xmr.e*xfr.s*wf[2] + xmr.e*xfa.s*wf[2] + xmr.e*xmr.s*wf[3] + xmr.e*xma.s*wf[3] + xma.e*xfr.s*wf[2] + xma.e*xfa.s*wf[2] + xma.e*xmr.s*wf[3] + xma.e*xma.s*wf[3]         
      wbarm <- xfr.e*yfr.s*wm[1,1] + xfr.e*yfa.s*wm[2,1] + xfr.e*ymr.s*wm[1,2] + xfr.e*yma.s*wm[2,2] + xfa.e*yfr.s*wm[2,1] + xfa.e*yfa.s*wm[2,1] + xfa.e*ymr.s*wm[2,2] + xfa.e*yma.s*wm[2,2] + xmr.e*yfr.s*wm[1,2] + xmr.e*yfa.s*wm[2,2] + xmr.e*ymr.s*wm[1,3] + xmr.e*yma.s*wm[2,3] + xma.e*yfr.s*wm[2,2] + xma.e*yfa.s*wm[2,2] + xma.e*ymr.s*wm[2,3] + xma.e*yma.s*wm[2,3]       
      
      #female gamete equations
      xfr.ep <- (xfr.e*xfr.s*wf[1] + 0.5*xfr.e*xfa.s*wf[1] + 0.5*xfr.e*xmr.s*wf[2] + 0.5*xfr.e*xma.s*wf[2]*(1-r2) + 0.5*xfa.e*xfr.s*wf[1] + 0.5*xfa.e*xmr.s*wf[2]*r2 + 0.5*xmr.e*xfr.s*wf[2] + 0.5*xmr.e*xfa.s*wf[2]*r2 + 0.5*xma.e*xfr.s*wf[2]*(1-r2))/wbarf
      xfa.ep <- (0.5*xfr.e*xfa.s*wf[1] + 0.5*xfr.e*xma.s*wf[2]*r2 + 0.5*xfa.e*xfr.s*wf[1] + xfa.e*xfa.s*wf[1] + 0.5*xfa.e*xmr.s*wf[2]*(1-r2) + 0.5*xfa.e*xma.s*wf[2] + 0.5*xmr.e*xfa.s*wf[2]*(1-r2) + 0.5*xma.e*xfr.s*wf[2]*r2 + 0.5*xma.e*xfa.s*wf[2])/wbarf
      xmr.ep <- (0.5*xfr.e*xmr.s*wf[2] + 0.5*xfr.e*xma.s*wf[2]*r2 + 0.5*xfa.e*xmr.s*wf[2]*(1-r2) + 0.5*xmr.e*xfr.s*wf[2] + 0.5*xmr.e*xfa.s*wf[2]*(1-r2) + xmr.e*xmr.s*wf[3] + 0.5*xmr.e*xma.s*wf[3] + 0.5*xma.e*xfr.s*wf[2]*r2 + 0.5*xma.e*xmr.s*wf[3])/wbarf
      xma.ep <- (0.5*xfr.e*xma.s*wf[2]*(1-r2) + 0.5*xfa.e*xmr.s*wf[2]*r2 + 0.5*xfa.e*xma.s*wf[2] + 0.5*xmr.e*xfa.s*wf[2]*r2 + 0.5*xmr.e*xma.s*wf[3] + 0.5*xma.e*xfr.s*wf[2]*(1-r2) + 0.5*xma.e*xfa.s*wf[2] + 0.5*xma.e*xmr.s*wf[3] + xma.e*xma.s*wf[3])/wbarf
      
      if (option == "Y" | option == "X") {
        xfr.sp <- (0.5*xfr.e*yfr.s*wm[1,1] + 0.5*xfr.e*yfa.s*wm[2,1] + 0.5*xfr.e*ymr.s*wm[1,2]*(1-r1) + 0.5*xfr.e*yma.s*wm[2,2] + 0.5*xmr.e*yfr.s*wm[1,2]*r1)/wbarm
        xfa.sp <- (0.5*xfa.e*yfr.s*wm[2,1] + 0.5*xfa.e*yfa.s*wm[2,1] + 0.5*xfa.e*ymr.s*wm[2,2] + 0.5*xfa.e*yma.s*wm[2,2])/wbarm
        xmr.sp <- (0.5*xmr.e*yfa.s*wm[2,2] + 0.5*xmr.e*ymr.s*wm[1,3] + 0.5*xmr.e*yfr.s*wm[1,2]*(1-r1) + 0.5*xmr.e*yma.s*wm[2,3] + 0.5*xfr.e*ymr.s*wm[1,2]*r1)/wbarm
        xma.sp <- (0.5*xma.e*yfr.s*wm[2,2] + 0.5*xma.e*yfa.s*wm[2,2] + 0.5*xma.e*ymr.s*wm[2,3] + 0.5*xma.e*yma.s*wm[2,3])/wbarm
        
        yfr.sp <- (0.5*xfr.e*yfr.s*wm[1,1] + 0.5*xfr.e*ymr.s*wm[1,2]*r1 + 0.5*xfa.e*yfr.s*wm[2,1] + 0.5*xmr.e*yfr.s*wm[1,2]*(1-r1) + 0.5*xma.e*yfr.s*wm[2,2])/wbarm
        yfa.sp <- (0.5*xfr.e*yfa.s*wm[2,1] + 0.5*xfa.e*yfa.s*wm[2,1] + 0.5*xmr.e*yfa.s*wm[2,2] + 0.5*xma.e*yfa.s*wm[2,2])/wbarm
        ymr.sp <- (0.5*xfr.e*ymr.s*wm[1,2]*(1-r1) + 0.5*xfa.e*ymr.s*wm[2,2] + 0.5*xmr.e*yfr.s*wm[1,2]*r1 + 0.5*xmr.e*ymr.s*wm[1,3] + 0.5*xma.e*ymr.s*wm[2,3])/wbarm
        yma.sp <- (0.5*xfr.e*yma.s*wm[2,2] + 0.5*xfa.e*yma.s*wm[2,2] + 0.5*xmr.e*yma.s*wm[2,3] + 0.5*xma.e*yma.s*wm[2,3])/wbarm
      }
      if (option == "A") {
        xfr.sp <- (0.5*xfr.e*yfr.s*wm[1,1] + 0.25*xfr.e*yfa.s*wm[2,1] + 0.5*xfr.e*ymr.s*wm[1,2]*(1-r1) + 0.25*xfr.e*yma.s*wm[2,2] + 0.25*xfa.e*yfr.s*wm[2,1] + 0.25*xfa.e*ymr.s*wm[2,2] + 0.5*xmr.e*yfr.s*wm[1,2]*r1)/wbarm
        xfa.sp <- (0.25*xfr.e*yfa.s*wm[2,1] + 0.25*xfr.e*yma.s*wm[2,2] + 0.25*xfa.e*yfr.s*wm[2,1] + 0.5*xfa.e*yfa.s*wm[2,1] + 0.25*xfa.e*ymr.s*wm[2,2] + 0.5*xfa.e*yma.s*wm[2,2])/wbarm
        xmr.sp <- (0.5*xfr.e*ymr.s*wm[1,2]*r1 + 0.5*xmr.e*yfr.s*wm[1,2]*(1-r1) + 0.25*xmr.e*yfa.s*wm[2,2] + 0.5*xmr.e*ymr.s*wm[1,3] + 0.25*xmr.e*yma.s*wm[2,3] + 0.25*xma.e*yfr.s*wm[2,2] + 0.25*xma.e*ymr.s*wm[2,3])/wbarm
        xma.sp <- (0.25*xmr.e*yfa.s*wm[2,2] + 0.25*xmr.e*yma.s*wm[2,3] + 0.25*xma.e*yfr.s*wm[2,2] + 0.25*xma.e*yfa.s*wm[2,2] + 0.25*xma.e*ymr.s*wm[2,3] + 0.5*xma.e*yma.s*wm[2,3])/wbarm
        
        yfr.sp <- (0.5*xfr.e*yfr.s*wm[1,1] + 0.25*xfr.e*yfa.s*wm[2,1] + 0.5*xfr.e*ymr.s*wm[1,2]*r1 + 0.25*xfa.e*yfr.s*wm[2,1] + 0.5*xmr.e*yfr.s*wm[1,2]*(1-r1) + 0.25*xmr.e*yfa.s*wm[2,2] + 0.25*xma.e*yfr.s*wm[2,2])/wbarm
        yfa.sp <- (0.25*xfr.e*yfa.s*wm[2,1] + 0.25*xfa.e*yfr.s*wm[2,1] + 0.5*xfa.e*yfa.s*wm[2,1] + 0.25*xmr.e*yfa.s*wm[2,2] + 0.25*xma.e*yfr.s*wm[2,2] + 0.25*xma.e*yfa.s*wm[2,2])/wbarm
        ymr.sp <- (0.5*xfr.e*ymr.s*wm[1,2]*(1-r1) + 0.25*xfr.e*yma.s*wm[2,2] + 0.25*xfa.e*ymr.s*wm[2,2] + 0.5*xmr.e*yfr.s*wm[1,2]*r1 + 0.5*xmr.e*ymr.s*wm[1,3] + 0.25*xmr.e*yma.s*wm[2,3] + 0.25*xma.e*ymr.s*wm[2,3])/wbarm
        yma.sp <- (0.25*xfr.e*yma.s*wm[2,2] + 0.25*xfa.e*ymr.s*wm[2,2] + 0.5*xfa.e*yma.s*wm[2,2] + 0.25*xmr.e*yma.s*wm[2,3] + 0.25*xma.e*ymr.s*wm[2,3] + 0.5*xma.e*yma.s*wm[2,3])/wbarm
      }
      
      #fill in new values
      xfr.e <- xfr.ep
      xfa.e <- xfa.ep
      xmr.e <- xmr.ep
      xma.e <- xma.ep
      xfr.s <- xfr.sp
      xfa.s <- xfa.sp
      xmr.s <- xmr.sp
      xma.s <- xma.sp
      yfr.s <- yfr.sp
      yfa.s <- yfa.sp
      ymr.s <- ymr.sp
      yma.s <- yma.sp
      
      #TODO select frequencies to track
      results[i,1] <- (xfr.e+xfr.s)/1.5
      results[i,2] <- (xfa.e+xfa.s)/1.5
      results[i,3] <- (xmr.e+xmr.s)/1.5
      results[i,4] <- (xma.e+xma.s)/1.5
      results[i,5] <- yfr.s/0.5
      results[i,6] <- yfa.s/0.5
      results[i,7] <- ymr.s/0.5
      results[i,8] <- yma.s/0.5
      
      #increase mutational load
      ml <- ml + me*u
      
      #recalculate male fitness after mutational load changes
      wm <- matrix(c(1-aneu, 1-ml, (1+h*s)*(1-aneu), (1+h*s)*(1-ml), (1+s)*(1-aneu), (1+s)*(1-ml)), 2, 3)
    }
    if (((option == "Y") & (sum(results[i-1,5:8] > 0.99)>0)) | ((option == "X") & (sum(results[i-1,1:4] > 0.99)>0))) {
      results[i,1] <- (xfr.e+xfr.s)/1.5
      results[i,2] <- (xfa.e+xfa.s)/1.5
      results[i,3] <- (xmr.e+xmr.s)/1.5
      results[i,4] <- (xma.e+xma.s)/1.5
      results[i,5] <- yfr.s/0.5
      results[i,6] <- yfa.s/0.5
      results[i,7] <- ymr.s/0.5
      results[i,8] <- yma.s/0.5
    }
  }
  return(results)
}

### Individual Runs ####
# Y chromosome
res <- simLife(h=1,r1=0.1,r2=0.1,s=0.5,me=0,gen=1000,aneu=0,option = "A")

#Plot X-carrying gametes
par(mfrow=c(1,2))
plot(res[,1], type = "l", col = viridis(1, begin = 0.25), main = "Gamete Frequency",
     ylim = c(0, 1), xlim = c(0, length(res[,1])),
     xlab = "Generation", ylab = "Frequency", lwd = 3, lty = 2)
lines(res[,2], col = viridis(1, begin = 0.25), lwd = 3)
lines(res[,3], col = viridis(1, begin = 0.75), lwd = 3, lty = 2)
lines(res[,4], col = viridis(1, begin = 0.75), lwd = 3)
abline(h=0, col = "black")
abline(h=1, col = "black")
# points(x=0, y=0.95, pch = 16, col = "red")
# points(x=0, y=0.9, pch = 16, col = "blue")
# points(x=0, y=0.85, pch = 16, col = "green")
# points(x=0, y=0.8, pch = 16, col = "orange")
# text(x=0, y=0.95,
#      paste("Female Benefit"), pos=4, cex=.7)
# text(x=0, y=0.9,
#      paste("Male Benefit"), pos=4, cex=.7)
# text(x=0, y=0.85,
#      paste("Achiasmy - Female Benefit"), pos=4, cex=.7)
# text(x=0, y=0.8,
#      paste("Achiasmy - Male Benefit"), pos=4, cex=.7)

#Plot Y-carrying gametes
plot(res[,5], type = "l", col = viridis(1, begin = 0.25), main = "Gamete Frequency",
     ylim = c(0, 1), xlim = c(0, length(res[,1])),
     xlab = "Generation", ylab = "Frequency", lwd = 3, lty = 2)
lines(res[,6], col = viridis(1, begin = 0.25), lwd = 3)
lines(res[,7], col = viridis(1, begin = 0.75), lwd = 3, lty = 2)
lines(res[,8], col = viridis(1, begin = 0.75), lwd = 3)
abline(h=0, col = "black")
abline(h=1, col = "black")
# points(x=0, y=0.95, pch = 16, col = "red")
# points(x=0, y=0.9, pch = 16, col = "blue")
# points(x=0, y=0.85, pch = 16, col = "green")
# points(x=0, y=0.8, pch = 16, col = "orange")
# text(x=0, y=0.95,
#      paste("Chiasmatic - Female Benefit"), pos=4, cex=.7)
# text(x=0, y=0.9,
#      paste("Chiasmatic - Male Benefit"), pos=4, cex=.7)
# text(x=0, y=0.85,
#      paste("Achiasmatic - Male Benefit"), pos=4, cex=.7)
# text(x=0, y=0.8,
#      paste("Achiasmatic - Female Benefit"), pos=4, cex=.7)

# X chromosome
res <- simLife(h=0,r1=0.1,r2=0.1,s=0.4,me=0.001,gen=5000,aneu=0,option = "X")
#Plot X-carrying gametes
par(mfrow=c(1,2))
plot(res[,1], type = "l", col = viridis(1, begin = 0.25), main = "Gamete Frequency",
     ylim = c(0, 1), xlim = c(0, length(res[,1])),
     xlab = "Generation", ylab = "Frequency", lwd = 3, lty = 2)
lines(res[,2], col = viridis(1, begin = 0.25), lwd = 3)
lines(res[,3], col = viridis(1, begin = 0.75), lwd = 3, lty = 2)
lines(res[,4], col = viridis(1, begin = 0.75), lwd = 3)
abline(h=0, col = "black")
abline(h=1, col = "black")
#Plot Y-carrying gametes
plot(res[,5], type = "l", col = viridis(1, begin = 0.25), main = "Gamete Frequency",
     ylim = c(0, 1), xlim = c(0, length(res[,1])),
     xlab = "Generation", ylab = "Frequency", lwd = 3, lty = 2)
lines(res[,6], col = viridis(1, begin = 0.25), lwd = 3)
lines(res[,7], col = viridis(1, begin = 0.75), lwd = 3, lty = 2)
lines(res[,8], col = viridis(1, begin = 0.75), lwd = 3)
abline(h=0, col = "black")
abline(h=1, col = "black")

#####

### Comparing ME and S values ####
steps <- 200
mes <- seq(from = 0, to = 0.2, length.out=steps)
svals <- seq(from = 0, to = 0.5, length.out=steps)
gen <- 1000

resY <- resX <- resA <- matrix(NA, steps, steps)
row.names(resY) <- row.names(resX) <- row.names(resA) <- paste("me.", round(mes, 3), sep="")
colnames(resY) <- colnames(resX)<- colnames(resA) <- paste("sval.", round(svals, 3), sep="")

# Y mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on me", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1,
                   s = svals[j], me = mes[i],
                   gen = gen, option = "Y")
    resY[i, j] <- (res[gen,6]+res[gen,8])
  }
}
# X mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on me", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1,
                   s = svals[j], me = mes[i],
                   gen = gen, option = "X")
    resX[i, j] <- (res[gen,2]+res[gen,4])
  }
}
# Autosome Mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on me", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1,
                   s = svals[j], me = mes[i],
                   gen = gen, option = "A")
    resA[i, j] <- (((res[gen,2]+res[gen,4])*1.5)+((res[gen,6]+res[gen,8])*0.5))/2
  }
}
write.csv(resY, file = "data/resY-me-s.csv", row.names = F)
write.csv(resX, file = "data/resX-me-s.csv", row.names = F)
write.csv(resA, file = "data/resA-me-s.csv", row.names = F)


# Plotting
resY <- as.matrix(read.csv("data/resY-me-s.csv"))
resX <- as.matrix(read.csv("data/resX-me-s.csv"))
resA <- as.matrix(read.csv("data/resA-me-s.csv"))

layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,
       heights = c(5, 5, 5, 5), 
       widths = c(5, 5, 5, 1))
par(mar = c(5.1,4.6,4.1,2.1))
image(resY, col = viridis(100, begin = 0.25, end = 1), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "Deleterious Mutation Effect",
      main = "Frequency on Y",
      cex.lab = 1.6,
      cex.main = 1.8)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.2, length.out=6),
     cex.axis = 1.3)
image(resX, col = viridis(100, begin = 0.25), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "Deleterious Mutation Effect",
      main = "Frequency on X",
      cex.lab = 1.6,
      cex.main = 1.8)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.2, length.out=6),
     cex.axis = 1.3)
image(resA, col = viridis(100, begin = 0.25), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "Deleterious Mutation Effect",
      main = "Frequency on Autosome",
      cex.lab = 1.6,
      cex.main = 1.8)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.2, length.out=6),
     cex.axis = 1.3)

image.plot(resY, col = viridis(100, begin = 0.2), legend.only = T,
           horizontal = F, legend.width = 5, legend.mar = 0,
           legend.shrink = 0.8)
#####

### ME vs Aneu ####
steps <- 200
mes <- seq(from = 0, to = 0.5, length.out=steps)
kvals <- seq(from = 0, to = 0.5, length.out=steps)
gen <- 1000

resY <- resX <- resA <- matrix(NA, steps, steps)
row.names(resY) <- row.names(resX) <- row.names(resA) <- paste("me.", round(mes, 3), sep="")
colnames(resY) <- colnames(resX)<- colnames(resA) <- paste("kval.", round(kvals, 3), sep="")

# Y mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on me", i))
  for (j in 1:steps) { # cycles through aneuploidy rates (k)
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1, s=0,
                   aneu = kvals[j], me = mes[i],
                   gen = gen, option = "Y")
    resY[i, j] <- (res[gen,6]+res[gen,8])
  }
}
# X mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on me", i))
  for (j in 1:steps) { # cycles through aneuploidy rates (k)
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1, s=0,
                   aneu = kvals[j], me = mes[i],
                   gen = gen, option = "X")
    resX[i, j] <- (res[gen,2]+res[gen,4])
  }
}
# A mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on me", i))
  for (j in 1:steps) { # cycles through aneuploidy rates (k)
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1, s=0,
                   aneu = kvals[j], me = mes[i],
                   gen = gen, option = "A")
    resA[i, j] <- (((res[gen,2]+res[gen,4])*1.5)+((res[gen,6]+res[gen,8])*0.5))/2
  }
}
write.csv(resY, file = "data/resY-me-k.csv", row.names = F)
write.csv(resX, file = "data/resX-me-k.csv", row.names = F)
write.csv(resA, file = "data/resA-me-k.csv", row.names = F)

# Plotting
resY <- as.matrix(read.csv("data/resY-me-k.csv"))
resX <- as.matrix(read.csv("data/resX-me-k.csv"))
resA <- as.matrix(read.csv("data/resA-me-k.csv"))

layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,
       heights = c(5, 5, 5, 5),
       widths = c(5, 5, 5, 1))
par(mar = c(5.1,4.6,4.1,2.1))
image(resY, col = viridis(100, begin = 0.25, end = 1), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Aneuploidy Rate",
      xlab = "Deleterious Mutation Effect",
      main = "Frequency on Y",
      cex.lab = 1.6,
      cex.main = 1.8)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)
image(resX, col = viridis(100, begin = 0.25), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Aneuploidy Rate",
      xlab = "Deleterious Mutation Effect",
      main = "Frequency on X",
      cex.lab = 1.6,
      cex.main = 1.8)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)
image(resA, col = viridis(100, begin = 0.25), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Aneuploidy Rate",
      xlab = "Deleterious Mutation Effect",
      main = "Frequency on Autosome",
      cex.lab = 1.6,
      cex.main = 1.8)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 1.3)

image.plot(resY, col = viridis(100, begin = 0.2), legend.only = T,
           horizontal = F, legend.width = 5, legend.mar = 0,
           legend.shrink = 0.8)
#####

### Comparing across different recombination frequencies ####
steps <- 200
r2s <- seq(from = 0, to = 0.5, length.out=steps)
r1s <- seq(from = 0, to = 0.5, length.out=steps)
gen <- 1000

resY <- resX <- resA <- matrix(NA, steps, steps)
row.names(resXr) <- row.names(resXa) <- row.names(resYr) <- row.names(resYa) <- paste("R2.", round(r2s, 3), sep="")
colnames(resXr) <- colnames(resXa)<- colnames(resYr) <- colnames(resYa) <- paste("R1.", round(r1s, 3), sep="")

# Y mutation
for (i in 1:steps) { # cycles through r2s
  print(paste("working on R2", i))
  for (j in 1:steps) { # cycles through r1s
    res <- simLife(r2 = r2s[i], r1 = r1s[j], h=0,
                   s = 0.1, me = 0.001,
                   gen = gen, option = "Y")
    resY[i, j] <- (res[gen,6]+res[gen,8])
  }
}
# X mutation
for (i in 1:steps) { # cycles through r2s
  print(paste("working on R2", i))
  for (j in 1:steps) { # cycles through r1s
    res <- simLife(r2 = r2s[i], r1 = r1s[j], h=0,
                   s = 0.1, me = 0.001,
                   gen = gen, option = "X")
    resX[i, j] <- (res[gen,2]+res[gen,4])
  }
}
# Autosome Mutation
### TODO

#Plotting
par(mfrow=c(1,3))
### TODO Fix labels
image(resY, col = viridis(100), yaxt='n', xaxt='n', ylab = "R1 (SDR - SAL)",
      xlab = "R2 (SAL - Achiasmy)", zlim = c(0,1), main = "Ya Frequency")
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6))
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6))
### TODO Fix labels
image(resX, col = viridis(100), yaxt='n', xaxt='n', ylab = "R1 (SDR - SAL)",
      xlab = "R2 (SAL - Achiasmy)", zlim = c(0,1), main = "Xa Frequency")
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6))
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6))
### TODO Plot Autosome

#####

