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
## acs - achiasmatic coding sites, sites that transition to achiasmy

## Main Function
simLife <- function(h = 1, s = 0.1, r1 = 0.1, r2 = 0.1,
                    gen = 1000, me = 0.003, u = 1e-09, acs = 5e+06,
                    aneu=0, option="Y"){
  #TODO
  # mutational load
  ml <- 0
  
  # fitness matrices
  wf <- c(1, 1/(1+h*s), 1/(1+s))
  wm <- matrix(c(1-aneu, 1-ml, (1+h*s)*(1-aneu), (1+h*s)*(1-ml), (1+s)*(1-aneu), (1+s)*(1-ml)), 2, 3)
  
  # Initial gamete frequencies
  # Autosomes: new mutation frequency = 1/(2N)
  # X chromosome: new mutation frequency = 2/(3N)
  # Y chromosome: new mutation frequency = 2/(N)
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
    xfa.e <- 0.0001
    xfr.e <- 0.4999
    xma.e <- 0.0001
    xmr.e <- 0.4999
    # sperm
    xfa.s <- 0.00005
    xfr.s <- 0.24995
    xma.s <- 0.00005
    xmr.s <- 0.24995
    yfr.s <- 0.25
    yfa.s <- 0
    ymr.s <- 0.25
    yma.s <- 0
  }
  
  if (option == "A") {
    # eggs
    xfa.e <- 0.0001
    xfr.e <- 0.4999
    xma.e <- 0.0001
    xmr.e <- 0.4999
    # sperm
    xfa.s <- 0.00005
    xfr.s <- 0.24995
    xma.s <- 0.00005
    xmr.s <- 0.24995
    yfa.s <- 0.00005
    yfr.s <- 0.24995
    yma.s <- 0.00005
    ymr.s <- 0.24995
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
        xma.sp <- (0.25*xmr.e*yfa.s*wm[2,2] + 0.25*xmr.e*yma.s*wm[2,3] + 0.25*xma.e*yfr.s*wm[2,2] + 0.5*xma.e*yfa.s*wm[2,2] + 0.25*xma.e*ymr.s*wm[2,3] + 0.5*xma.e*yma.s*wm[2,3])/wbarm
        
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
      #ml <- ml + me*u*acs
      ml <- 1 - ( (1 - me)^(u * acs * i) )
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
par(mfrow=c(1,1), mar=c(6, 6, 4, 1), oma = c(2, 2, 0, 5))

res <- simLife(h = 1, s = 0.01, r1 = 0.1, r2 = 0.1,
               gen = 5000, me = 0.001, u = 1e-09, acs = 2e+07,
               aneu=0, option="Y")

plot((res[,6]+res[,8]), type = "l", col = "#317ec2", main = "Gamete Frequency\n Sexual Antagonism",
     ylim = c(0, 1), xlim = c(0, length(res[,1])),
     xlab = "", ylab = "Frequency", lwd = 3, cex.lab = 1.3)

#Comparison across chromosomes
## SA
#par(mfrow=c(2,1), mar=c(6, 6, 4, 1), oma = c(2, 2, 0, 5))
res <- simLife(h=1,r1=0.1,r2=0.1,s=0.5,gen=3000,aneu=0,option = "Y")
plot((res[,6]+res[,8]), type = "l", col = "#317ec2", main = "Gamete Frequency\n Sexual Antagonism",
     ylim = c(0, 1), xlim = c(0, length(res[,1])),
     xlab = "", ylab = "Frequency", lwd = 3, cex.lab = 1.3)
abline(v=150, lty=3, lwd = 2, col = "#317ec2")
res <- simLife(h=1,r1=0.1,r2=0.1,s=0.5,gen=3000,aneu=0,option = "X")
lines((res[,2]+res[,4]), col = "#c03830", lwd = 3)
abline(v=650, lty=3, lwd = 2, col = "#c03830")
res <- simLife(h=1,r1=0.1,r2=0.1,s=0.5,gen=3000,aneu=0,option = "A")
lines((res[,2]+res[,4]+res[,6]+res[,8])/2, col = "#5aaa46", lwd = 3)
abline(v=2700, lty=3, lwd = 2, col = "#5aaa46")
## Legend
legend(x=1200, y=0.3, legend = c("Y Chromosome", "X Chromosome", "Autosomes", "Dotted Lines: Equilbrium Point"),
       col = c("#317ec2", "#c03830", "#5aaa46", "black"), lty = c(1, 1, 1, 3), lwd = c(6, 6, 6, 3), bty ="n", cex = 1.2)
## Aneu
res <- simLife(h=1,r1=0.1,r2=0.1,s=0,gen=3000,aneu=0.1,option = "Y")
plot((res[,6]+res[,8]), type = "l", col = "#317ec2", main = "Gamete Frequency\n Heteromorphy-Dependant Aneuploidy",
     ylim = c(0, 1), xlim = c(0, length(res[,1])),
     xlab = "Generation", ylab = "Frequency", lwd = 3, cex.lab = 1.3)
abline(v=150, lty=3, lwd = 2, col = "#317ec2")
res <- simLife(h=1,r1=0.1,r2=0.1,s=0,gen=3000,aneu=0.1,option = "X")
lines((res[,2]+res[,4]), col = "#c03830", lwd = 3)
abline(v=450, lty=3, lwd = 2, col = "#c03830")
res <- simLife(h=1,r1=0.1,r2=0.1,s=0,gen=3000,aneu=0.1,option = "A")
lines((res[,2]+res[,4]+res[,6]+res[,8])/2, col = "#5aaa46", lwd = 3)
abline(v=1950, lty=3, lwd = 2, col = "#5aaa46")


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
acss <- seq(from = 0, to = 5e+06, length.out=steps)
svals <- seq(from = 0, to = 0.5, length.out=steps)
gen <- 1000

resY <- resX <- resA <- matrix(NA, steps, steps)
row.names(resY) <- row.names(resX) <- row.names(resA) <- paste("acs.", round(acss, 3), sep="")
colnames(resY) <- colnames(resX)<- colnames(resA) <- paste("sval.", round(svals, 3), sep="")

## h1-r01
# Y mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "Y")
    resY[i, j] <- (res[gen,6]+res[gen,8])
  }
}
# X mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "X")
    resX[i, j] <- (res[gen,2]+res[gen,4])
  }
}
# Autosome Mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "A")
    resA[i, j] <- (((res[gen,2]+res[gen,4])*1.5)+((res[gen,6]+res[gen,8])*0.5))/2
  }
}
write.csv(resY, file = "data/resY-h1-r01.csv", row.names = F)
write.csv(resX, file = "data/resX-h1-r01.csv", row.names = F)
write.csv(resA, file = "data/resA-h1-r01.csv", row.names = F)

## h05-r01
# Y mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=0.5,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "Y")
    resY[i, j] <- (res[gen,6]+res[gen,8])
  }
}
# X mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=0.5,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "X")
    resX[i, j] <- (res[gen,2]+res[gen,4])
  }
}
# Autosome Mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=0.5,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "A")
    resA[i, j] <- (((res[gen,2]+res[gen,4])*1.5)+((res[gen,6]+res[gen,8])*0.5))/2
  }
}
write.csv(resY, file = "data/resY-h05-r01.csv", row.names = F)
write.csv(resX, file = "data/resX-h05-r01.csv", row.names = F)
write.csv(resA, file = "data/resA-h05-r01.csv", row.names = F)

## h0-r01
# Y mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=0,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "Y")
    resY[i, j] <- (res[gen,6]+res[gen,8])
  }
}
# X mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=0,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "X")
    resX[i, j] <- (res[gen,2]+res[gen,4])
  }
}
# Autosome Mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.1, h=0,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "A")
    resA[i, j] <- (((res[gen,2]+res[gen,4])*1.5)+((res[gen,6]+res[gen,8])*0.5))/2
  }
}
write.csv(resY, file = "data/resY-h0-r01.csv", row.names = F)
write.csv(resX, file = "data/resX-h0-r01.csv", row.names = F)
write.csv(resA, file = "data/resA-h0-r01.csv", row.names = F)

## h1-r03
# Y mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.3, h=1,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "Y")
    resY[i, j] <- (res[gen,6]+res[gen,8])
  }
}
# X mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.3, h=1,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "X")
    resX[i, j] <- (res[gen,2]+res[gen,4])
  }
}
# Autosome Mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through the selection coef.
    res <- simLife(r2 = 0.1, r1 = 0.3, h=1,
                   s = svals[j], acs = acss[i],
                   gen = gen, option = "A")
    resA[i, j] <- (((res[gen,2]+res[gen,4])*1.5)+((res[gen,6]+res[gen,8])*0.5))/2
  }
}
write.csv(resY, file = "data/resY-h1-r03.csv", row.names = F)
write.csv(resX, file = "data/resX-h1-r03.csv", row.names = F)
write.csv(resA, file = "data/resA-h1-r03.csv", row.names = F)

#### h1-r01 plot ####
resY <- as.matrix(read.csv("data/resY-h1-r01.csv"))
resX <- as.matrix(read.csv("data/resX-h1-r01.csv"))
resA <- as.matrix(read.csv("data/resA-h1-r01.csv"))

layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,
       heights = c(5, 5, 5, 5), 
       widths = c(5, 5, 5, 1))
par(mar = c(10.2,9.2,8.2,4.2),
    mgp = c(5,1.5,0))
image(resY, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "",
      main = "Frequency on Y",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resX, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "Coding Sites Transitioning to Achiasmy (Mb)",
      main = "Frequency on X",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resA, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "",
      main = "Frequency on Autosome",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)

# image.plot(resY, col = viridis::mako(4000), legend.only = T,
#            horizontal = F, legend.width = 4, legend.mar = 0,
#            legend.shrink = 0.8)
#####
#### h05-r01 plot ####
resY <- as.matrix(read.csv("data/resY-h05-r01.csv"))
resX <- as.matrix(read.csv("data/resX-h05-r01.csv"))
resA <- as.matrix(read.csv("data/resA-h05-r01.csv"))

layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,
       heights = c(5, 5, 5, 5), 
       widths = c(5, 5, 5, 1))
par(mar = c(10.2,9.2,8.2,4.2),
    mgp = c(5,1.5,0))
image(resY, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "",
      main = "Frequency on Y",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resX, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "Coding Sites Transitioning to Achiasmy (Mb)",
      main = "Frequency on X",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resA, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "",
      main = "Frequency on Autosome",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)

# image.plot(resY, col = viridis::mako(4000), legend.only = T,
#            horizontal = F, legend.width = 5, legend.mar = 0,
#            legend.shrink = 0.8)
#####
#### h0-r01 plot ####
resY <- as.matrix(read.csv("data/resY-h0-r01.csv"))
resX <- as.matrix(read.csv("data/resX-h0-r01.csv"))
resA <- as.matrix(read.csv("data/resA-h0-r01.csv"))

layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,
       heights = c(5, 5, 5, 5), 
       widths = c(5, 5, 5, 1))
par(mar = c(10.2,9.2,8.2,4.2),
    mgp = c(5,1.5,0))
image(resY, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "",
      main = "Frequency on Y",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resX, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "Coding Sites Transitioning to Achiasmy (Mb)",
      main = "Frequency on X",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resA, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "",
      main = "Frequency on Autosome",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)

# image.plot(resY, col = viridis::mako(4000), legend.only = T,
#            horizontal = F, legend.width = 5, legend.mar = 0,
#            legend.shrink = 0.8)
#####
#### h1-r03 plot ####
resY <- as.matrix(read.csv("data/resY-h1-r03.csv"))
resX <- as.matrix(read.csv("data/resX-h1-r03.csv"))
resA <- as.matrix(read.csv("data/resA-h1-r03.csv"))

layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,
       heights = c(5, 5, 5, 5), 
       widths = c(5, 5, 5, 1))
par(mar = c(10.2,9.2,8.2,4.2),
    mgp = c(5,1.5,0))
image(resY, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "",
      main = "Frequency on Y",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resX, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "Coding Sites Transitioning to Achiasmy (Mb)",
      main = "Frequency on X",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resA, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Selection Coefficient",
      xlab = "",
      main = "Frequency on Autosome",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.5, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)

# image.plot(resY, col = viridis::mako(4000), legend.only = T,
#            horizontal = F, legend.width = 5, legend.mar = 0,
#            legend.shrink = 0.8)
#####

### ME vs Aneu ####
steps <- 200
acss <- seq(from = 0, to = 5e+06, length.out=steps)
kvals <- seq(from = 0, to = 0.05, length.out=steps)
gen <- 1000

resY <- resX <- resA <- matrix(NA, steps, steps)
row.names(resY) <- row.names(resX) <- row.names(resA) <- paste("acs.", round(acss, 3), sep="")
colnames(resY) <- colnames(resX)<- colnames(resA) <- paste("kval.", round(kvals, 3), sep="")

# Y mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through aneuploidy rates (k)
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1, s=0,
                   aneu = kvals[j], acs = acss[i],
                   gen = gen, option = "Y")
    resY[i, j] <- (res[gen,6]+res[gen,8])
  }
}
# X mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through aneuploidy rates (k)
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1, s=0,
                   aneu = kvals[j], acs = acss[i],
                   gen = gen, option = "X")
    resX[i, j] <- (res[gen,2]+res[gen,4])
  }
}
# A mutation
for (i in 1:steps) { # cycles through the mutational effects factors
  print(paste("working on acs", i))
  for (j in 1:steps) { # cycles through aneuploidy rates (k)
    res <- simLife(r2 = 0.1, r1 = 0.1, h=1, s=0,
                   aneu = kvals[j], acs = acss[i],
                   gen = gen, option = "A")
    resA[i, j] <- (((res[gen,2]+res[gen,4])*1.5)+((res[gen,6]+res[gen,8])*0.5))/2
  }
}
write.csv(resY, file = "data/resY-k.csv", row.names = F)
write.csv(resX, file = "data/resX-k.csv", row.names = F)
write.csv(resA, file = "data/resA-k.csv", row.names = F)

# Plotting
resY <- as.matrix(read.csv("data/resY-k.csv"))
resX <- as.matrix(read.csv("data/resX-k.csv"))
resA <- as.matrix(read.csv("data/resA-k.csv"))

layout.matrix <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,
       heights = c(5, 5, 5, 5), 
       widths = c(5, 5, 5, 1))
par(mar = c(10.2,9.2,8.2,4.2),
    mgp = c(5,1.5,0))
image(resY, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Aneuploidy Rate",
      xlab = "",
      main = "Frequency on Y",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.05, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resX, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Aneuploidy Rate",
      xlab = "Coding Sites Transitioning to Achiasmy (Mb)",
      main = "Frequency on X",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.05, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)
image(resA, col = viridis::mako(4000), yaxt='n', xaxt='n', zlim = c(0,1),
      ylab = "Aneuploidy Rate",
      xlab = "",
      main = "Frequency on Autosome",
      cex.lab = 3.2,
      cex.main = 3.6)
axis(2, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 0.05, length.out=6),
     cex.axis = 2.2)
axis(1, at = seq(from = 0, to = 1, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 2.2)

# par(mfrow=c(1,1))
# image.plot(resY, col = viridis(100, begin = 0.25), legend.only = T,
#            horizontal = T, legend.width = 1, legend.mar = 5,
#            legend.shrink = 0.8)
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
