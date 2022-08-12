## Time-stamp: "Last modified 2022-03-16 14:01:11 delucia"
library(Rcpp)
library(RcppEigen)
library(ReacTran)
library(deSolve)

options(width=110)

setwd("app")

## This creates the "diff1D" function with our BTCSdiffusion code
sourceCpp("Rcpp-BTCS-1d.cpp")

### FTCS explicit (same name)
sourceCpp("RcppFTCS.cpp")



## Grid 101
## Set initial conditions
N <- 1001
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
sigma <- 0.02
Yini <- 0.5*exp(-0.5*((x-1/2.0)**2)/sigma**2)

## plot(x, Yini, type="l")

parms1 <- list(D=D.coeff)
# 1 timestep, 10 s
times <- seq(from = 0, to = 1, by = 0.1)

system.time({
    out1 <- ode.1D(y = Yini, times = times, func = Diffusion,
                  parms = parms1, dimens = N)[11,-1]
})

## Now with BTCS 
alpha <- rep(D.coeff, N)

system.time({
    out2 <- diff1D(n=N, length=1, field=Yini, alpha=alpha, timestep = 0.1, 0, 0, iterations = 10)
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



matplot(cbind(Yini,out1, out2, out3, mm),type="l", col=c("grey","black","red","blue","green4"), lty="solid", lwd=2,
        xlab="grid element", ylab="Concentration", las=1)
legend("topright", c("init","ReacTran ode1D", "BTCS 1d", "FTCS", "poor man's CN"), text.col=c("grey","black","red","blue","green4"), bty = "n")


sum(Yini)
sum(out1)
sum(out2)
sum(out3)
sum(mm)

## Yini <- 0.2*sin(pi/0.1*x)+0.2
## plot(Yini)

## plot(out3)

Fun <- function(dx) {
    tmp <- diff1D(n=N, length=1, field=Yini, alpha=alpha, timestep = dx, 0, 0, iterations = floor(1/dx))
    sqrt(sum({out1-tmp}^2))
}

reso <- optimise(f=Fun, interval=c(1E-5, 1E-1), maximum = FALSE)


dx <- 0.0006038284
floor(1/dx)

1/dx

system.time({
    out2o <- diff1D(n=N, length=1, field=Yini, alpha=alpha, timestep = dx, 0, 0, iterations = 1656)
})

matplot(cbind(out1, out2o),type="l", col=c("black","red"), lty="solid", lwd=2,
        xlab="grid element", ylab="Concentration", las=1)
legend("topright", c("ReacTran ode1D", "BTCS 1d dx=0.0006"), text.col=c("black","red"), bty = "n")


dx <- 0.05

system.time({
    out2o <- diff1D(n=N, length=1, field=Yini, alpha=alpha, timestep = dx, 0, 0, iterations = 1/dx)
})

matplot(cbind(out1, out2o),type="l", col=c("black","red"), lty="solid", lwd=2,
        xlab="grid element", ylab="Concentration", las=1)
legend("topright", c("ReacTran ode1D", "BTCS 1d dx=0.0006"), text.col=c("black","red"), bty = "n")

Matplot




## This creates the "diff1D" function with our BTCSdiffusion code
sourceCpp("Rcpp-BTCS-2d.cpp")

n <- 256
a2d <- rep(1E-3, n^2)

init2d <- readRDS("gs1.rds")

ll <- {init2d - min(init2d)}/diff(range(init2d))

system.time({
    res1 <- diff2D(nx=N, ny=N, lenx=1, leny=1, field=ll, alpha=a2d, timestep = 0.1, iterations = 10)
})

hist(ll,32)














