## Time-stamp: "Last modified 2022-03-15 18:04:27 delucia"
library(Rcpp)
library(RcppEigen)
library(ReacTran)
library(deSolve)

options(width=110)

setwd("app")

## This creates the "diff1D" function with our BTCSdiffusion code
sourceCpp("Rcpp-interface.cpp")

### FTCS explicit (same name)
sourceCpp("RcppFTCS.cpp")



## Grid 101
## Set initial conditions
N <- 501
D.coeff <- 1E-3
C0 <- 1             ## Initial concentration (mg/L)
X0 <- 0             ## Location of initial concentration (m)
## Yini <- c(C0, rep(0,N-1))

## Ode1d solution
xgrid <- setup.grid.1D(x.up = 0, x.down = 1, N = N)
x <- xgrid$x.mid
Diffusion <- function (t, Y, parms){
  tran <- tran.1D(C = Y, C.up = 0, C.down = 0, D = parms$D, dx = xgrid)
  return(list(tran$dC))
}


## gaussian pulse as initial condition
sigma <- 0.01
Yini <- exp(-0.5*((x-1/2.0)**2)/sigma**2)

## plot(x, Yini, type="l")

parms1 <- list(D=D.coeff)
# 1 timestep, 10 s
times <- seq(from = 0, to = 1, by = 1)

system.time({
    out1 <- ode.1D(y = Yini, times = times, func = Diffusion,
                  parms = parms1, dimens = N)[2,-1]
})

## Now with BTCS 
alpha <- rep(D.coeff, N)

system.time({
    out2 <- diff1D(n=N, length=1, field=Yini, alpha=alpha, timestep = 1, 0, 0, iterations = 1)
})

## plot(out1, out2)
## abline(0,1)

## matplot(cbind(out1,out2),type="l", col=c("black","red"),lty="solid", lwd=2,
##         xlab="grid element", ylab="Concentration", las=1)
## legend("topright", c("ReacTran ode1D", "BTCS 1d"), text.col=c("black","red"), bty = "n")





system.time({
    out3 <- RcppFTCS(n=N, length=1, field=Yini, alpha=1E-3, bc_left = 0, bc_right = 0, timestep = 1)
})

## Poor man's 
mm <- colMeans(rbind(out2,out3))



matplot(cbind(Yini, out1, out2, out3, mm),type="l", col=c("black","grey","red","blue","green4"), lty=c("dashed", rep("solid",4)), lwd=2,
        xlab="grid element", ylab="Concentration", las=1)
legend("topright", c("Init","ReacTran ode1D", "BTCS 1d", "FTCS", "poor man's CN"), text.col=c("black","grey","red","blue","green4"), bty = "n")


