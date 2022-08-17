## Time-stamp: "Last modified 2022-01-25 11:22:30 delucia"
library(ReacTran)
library(deSolve)
options(width=114)


Diffusion <- function (t, Y, parms){
  tran <- tran.1D(C = Y, C.up = 1, C.down = 0, D = parms$D, dx = xgrid)
  return(list(tran$dC))
}


## sol 100
## Set initial conditions
N <- 100
xgrid <- setup.grid.1D(x.up = 0, x.down = 1, N = N)
x <- xgrid$x.mid
D.coeff <- 1E-3
C0 <- 1             ## Initial concentration (mg/L)
X0 <- 0             ## Location of initial concentration (m)
Yini <- c(C0, rep(0,N-1))
parms1 <- list(D=D.coeff)
times <- seq(from = 0, to = 5, by = 1)
system.time({
    out100 <- ode.1D(y = Yini, times = times, func = Diffusion,
                     parms = parms1, dimens = N)
})


## sol 1000
## Set initial conditions
N <- 1000
xgrid <- setup.grid.1D(x.up = 0, x.down = 1, N = N)
x <- xgrid$x.mid
D.coeff <- 1E-3
C0 <- 1             ## Initial concentration (mg/L)
X0 <- 0             ## Location of initial concentration (m)
Yini <- c(C0, rep(0,N-1))
parms1 <- list(D=D.coeff)
times <- seq(from = 0, to = 5, by = 1)

system.time({
    out1000 <- ode.1D(y = Yini, times = times, func = Diffusion,
                      parms = parms1, dimens = N)
})

## sol 5000
## Set initial conditions
N <- 5000
xgrid <- setup.grid.1D(x.up = 0, x.down = 1, N = N)
x <- xgrid$x.mid
D.coeff <- 1E-3
C0 <- 1             ## Initial concentration (mg/L)
X0 <- 0             ## Location of initial concentration (m)
Yini <- c(C0, rep(0,N-1))
parms1 <- list(D=D.coeff)
times <- seq(from = 0, to = 5, by = 1)
system.time({
    out5000 <-        
        ode.1D(y = Yini, times = times, func = Diffusion, parms = parms1,
               dimens = N)
})

## ./mdl_main 100 > Sol100.dat
## ./mdl_main 1000 > Sol1000.dat
## ./mdl_main 5000 > Sol5000.dat
a100 <- data.table::fread("build/src/Sol100.dat")
a1000 <- data.table::fread("build/src/Sol1000.dat")
a5000 <- data.table::fread("build/src/Sol5000.dat")


## run: phreeqc pqc.in pqc.out /PATH_TO_PHREEQC.DAT
p100 <- data.table::fread("./selected_output_1.sel")

pqc <- subset(p100, subset = time==5) |>
    subset(subset = soln > 0 & soln < 101) |>
    subset(select = Mg)
    
cairo_pdf("test100.pdf", family="serif")
matplot.1D(out100, type = "l", subset = time == 4, xaxs="i", xlab="grid element", ylab="C",
           main="Comparison ode.1D vs btcs, discretization 100", lwd=2)
lines(a100$inp1, lwd=2, col="blue")
lines(a100$inp2, lwd=2, col="red")
lines(a100$inp3, lwd=2, col="green4")
lines(pqc$Mg, lwd=2, col="orange")
legend("topright", c("R deSolve ode.1D reference, ts=1",
                     "btcs, ts=1", "btcs, ts=0.1","btcs, ts=5","PHREEQC ts=1"), lty="solid",
       lwd=2, col=c("black","blue","red", "green4", "orange"),
       bty="n")
dev.off()

cairo_pdf("test1000.pdf", family="serif")
matplot.1D(out1000, type = "l", subset = time == 4, xaxs="i", xlab="grid element", ylab="C",
           main="Comparison ode.1D vs btcs, discretization 1000", lwd=2)
lines(a1000$inp1, lwd=2, col="blue")
lines(a1000$inp2, lwd=2, col="red")
lines(a1000$inp3, lwd=2, col="green4")
legend("topright", c("R deSolve ode.1D reference, ts=1",
                     "btcs, ts=1", "btcs, ts=0.1","btcs, ts=5"), lty="solid",
       lwd=2, col=c("black","blue","red", "green4"),
       bty="n")
dev.off()

cairo_pdf("test5000.pdf", family="serif")
matplot.1D(out5000, type = "l", subset = time == 4, xaxs="i", xlab="grid element", ylab="C",
           main="Comparison ode.1D vs btcs, discretization 5000", lwd=2)
lines(a5000$inp1, lwd=2, col="blue")
lines(a5000$inp2, lwd=2, col="red")
lines(a5000$inp3, lwd=2, col="green4")
legend("topright", c("R deSolve ode.1D reference, ts=1",
                     "btcs, ts=1", "btcs, ts=0.1","btcs, ts=5"), lty="solid",
       lwd=2, col=c("black","blue","red", "green4"),
       bty="n")
dev.off()


############## old stuff

a <- data.table::fread("build/src/Sol1000.dat")

erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)


Analytical <- function(x, t, a, v) {
    erfc
}

cairo_pdf("test.pdf", family="serif")
matplot.1D(out, type = "l", subset = time == 4, xaxs="i", xlab="grid element", ylab="C",
           main="Comparison ode.1D vs btcs", lwd=2)
lines(a$inp1, lwd=2, col="blue")
lines(a$inp2, lwd=2, col="red")
lines(a$inp3, lwd=2, col="green4")
legend("topright", c("R deSolve ode.1D reference, ts=1",
                     "btcs, ts=1", "btcs, ts=0.1","btcs, ts=5"), lty="solid",
       lwd=2, col=c("black","blue","red", "green4"),
       bty="n")
dev.off()


## sol 5000
## Set initial conditions
N <- 5000
xgrid <- setup.grid.1D(x.up = 0, x.down = 1, N = N)
x <- xgrid$x.mid
D.coeff <- 1E-3
C0 <- 1             ## Initial concentration (mg/L)
X0 <- 0             ## Location of initial concentration (m)
Yini <- c(C0, rep(0,N-1))

parms1 <- list(D=D.coeff)

times <- seq(from = 0, to = 5, by = 1)
system.time(
    out5000 <- ode.1D(y = Yini, times = times, func = Diffusion,
                      parms = parms1, dimens = N))

a <- data.table::fread("build/src/Sol5000.dat")

cairo_pdf("test.pdf", family="serif")
matplot.1D(out5000, type = "l", subset = time == 4, xaxs="i", xlab="grid element", ylab="C",
           main="Comparison ode.1D vs btcs", lwd=2)
lines(a$inp1, lwd=2, col="blue")
lines(a$inp2, lwd=2, col="red")
lines(a$inp3, lwd=2, col="green4")
legend("topright", c("R deSolve ode.1D reference, ts=1",
                     "btcs, ts=1", "btcs, ts=0.1","btcs, ts=5"), lty="solid",
       lwd=2, col=c("black","blue","red", "green4"),
       bty="n")
dev.off()



### ##
## diffusR <- function(t, C, parameters) {
##   deltax  <- c (0.5*delx, rep(delx, numboxes-1), 0.5*delx)
##   Flux    <- -D*diff(c(0, C, 0))/deltax
##   dC <- -diff(Flux)/delx 

##   list(dC)   # the output
## }
  
## ## ==================
## ## Model application
## ## ==================

## ## the model parameters:
## D         <- 0.3    # m2/day  diffusion rate
## r         <- 0.01   # /day    net growth rate
## delx      <- 1      # m       thickness of boxes
## numboxes  <- 50 

## ## distance of boxes on plant, m, 1 m intervals
## Distance  <- seq(from = 0.5, by = delx, length.out = numboxes)

## ## Initial conditions, ind/m2
## ## aphids present only on two central boxes
## C        <- rep(0, times = numboxes)
## C[30:31] <- 1
## state         <- c(C = C)      # initialise state variables 
                  
## ## RUNNING the model:
## times <- seq(0, 200, by = 1)   # output wanted at these time intervals
## out   <- ode.band(state, times, diffusR, parms = 0, 
##                   nspec = 1, names = "C")

## ## ================
## ## Plotting output
## ## ================
## image(out, grid = Distance, method = "filled.contour", 
##       xlab = "time, days", ylab = "Distance on plant, m",
##       main = "Aphid density on a row of plants")

## matplot.1D(out, grid = Distance, type = "l", 
##    subset = time %in% seq(0, 200, by = 10))

## # add an observed dataset to 1-D plot (make sure to use correct name):
## data <- cbind(dist  = c(0,10, 20,  30,  40, 50, 60), 
##               Aphid = c(0,0.1,0.25,0.5,0.25,0.1,0))

## matplot.1D(out, grid = Distance, type = "l", 
##    subset = time %in% seq(0, 200, by = 10), 
##    obs = data, obspar = list(pch = 18, cex = 2, col="red"))

## ## Not run: 
## plot.1D(out, grid = Distance, type = "l")

