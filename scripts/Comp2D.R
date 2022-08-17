## Time-stamp: "Last modified 2022-03-01 15:41:58 delucia"
library(ReacTran)
library(deSolve)
options(width=114)

N     <- 51          # number of grid cells
XX    <- 1           # total size
dy    <- dx <- XX/N  # grid size
Dy    <- Dx <- 0.1   # diffusion coeff, X- and Y-direction

N2  <- ceiling(N/2)

# The model equations

Diff2D <- function (t, y, parms)  {
 CONC  <- matrix(nrow = N, ncol = N, y)
 # CONC[25, 25] <- 1
 dCONC <- tran.2D(CONC, D.x = Dx, D.y = Dy, dx = dx, dy = dy)$dC

 return (list(dCONC))

}

# initial condition: 0 everywhere, except in central point
y <- matrix(nrow = N, ncol = N, data = 0)
y[N2, N2] <- 1  # initial concentration in the central point...

as.numeric(y)[1301]

# solve for 8 time units
times <- 0:5
out <- ode.2D (y = y, func = Diff2D, t = times, parms = NULL,
                dim = c(N,N), lrw=155412)


image(out, ask = FALSE, mfrow = c(2, 3), main = paste("time", times), asp=1)

str(out)

## adi <- data.matrix(data.table::fread("./bu2d/src/test2D", header=FALSE))
## str(adi)

adi <- data.table::fread("./bu2d/src/tt", header=FALSE)

attributes(adi) <- attributes(out)

madi1 <- matrix(adi[1,-1], 51,51,byrow = TRUE)
madi2 <- matrix(adi[2,-1], 51,51,byrow = TRUE)
madi6 <- matrix(adi[6,-1], 51,51,byrow = TRUE)
mout6 <- matrix(out[6,-1], 51,51,byrow = TRUE)

madi6[1300]


madi6 <- matrix(unlist(adi[6,-1]), 501,501,byrow = TRUE)
class(madi6)
image(madi6, asp=1)
contour(madi6, asp=1)

par(mfrow = c(1,3))
image(madi1, asp=1)
image(madi2, asp=1)
image(madi5, asp=1)

range(out[6, -1])

x11()

plot(adi[6, -1], out[6, -1], pch=4)
abline(0,1)

range(o6 <- out[6, -1])
range(o2 <- out[2, -1])

image(madi6)

contour(madi6)

contour(mout6)


library(ReacTran)

#####################################################
## FIRST SIMPLEST MODEL:  Fick second law for diffusion
# d^2(C(x,t))/dx^2 = D. d^2(C(x,t))/dx^2
# x: position
# t: time
# D: Diffusion coefficient
# C(x,t): Concentration at position x and time t

# Model definition

Fick_model=function(t=0,y,parms=NULL)
{
	
	tran= tran.1D(C=y,D=1000,dx=grid)$dC
	return(list(dC=tran))
}

# Grid definition

C <- dnorm(seq(-10,10,.1))
L <- 200   # length of the diffusion space
N <- length(C)  # number of grid layers
grid <- setup.grid.1D(x.up = 0,L = L, N = N)


# Initial conditions + simulation 

Fick_solved_20 <- ode.1D(y = C, times=seq(0,20),func=Fick_model,nspec=1,method="lsoda")
Fick_solved_20000 <- ode.1D(y = C, times=seq(0,20000),func=Fick_model,nspec=1,method="lsoda")


# Plot Results

plot(Fick_solved_20[1,2:dim(Fick_solved_20)[2]],type="l",xlab="x",ylab="Concentration",lwd=2,ylim=c(0,.5))
par(new=T)
plot(Fick_solved_20[2,2:dim(Fick_solved_20)[2]],type="l",xlab="x",ylab="Concentration",lwd=2,ylim=c(0,.5),col="blue")
par(new=T)
plot(Fick_solved_20000[2,2:dim(Fick_solved_20000)[2]],type="l",xlab="x",ylab="Concentration",lwd=2,ylim=c(0,.5),col="red")

Fick_solved_20[2,2:dim(Fick_solved_20)[2]]==Fick_solved_20000[2,2:dim(Fick_solved_20000)[2]]

mass <- function(time, state, params) {
     with(as.list(c(state, params)), {
             H2O <- K3 * H * O
             dH <- -K1 * H2O
             dO <- -K2 * H2O
             dH2O <- H2O
             list(c(dH, dO, dH2O))
     })
}

params <- c(K1 = 2,
             K2 = 1,
             K3 = 0.005)
state <- c(H = 200,
            O = 100,
            H2O = 0)
time <- seq(0,5)
out <- ode(state, time, mass, params)
plot(out)
