
## Brutal implementation of 2D ADI scheme
## Square NxN grid with dx=dy=1
ADI <- function(n, dt, iter, alpha) {

    nx <- ny <- n
    dx <- dy <- 1
    field <- matrix(0, nx, ny)
    ## find out the center of the grid to apply conc=1
    cen <- ceiling(n/2)
    field[cen, cen] <- 1

    ## prepare containers for computations and outputs
    tmpX <- tmpY <- res <- field
    out <- vector(mode="list", length=iter)
    
    for (it in seq(1, iter)) {
        for (i in seq(1, ny))
            tmpX[i,] <- SweepByRow(i, res, dt=dt, alpha=alpha)

        resY <- t(tmpX)
        for (i in seq(1, nx))
            tmpY[i,] <- SweepByRow(i, resY, dt=dt, alpha=alpha)

        res <- t(tmpY)
        out[[it]] <- res }

    return(out)
}


## Workhorse function to fill A, B and solve for a given *row* of the
## grid matrix
SweepByRow <- function(i, field, dt, alpha) {
    dx <- 1 ## fixed in our test
    A <- matrix(0, nrow(field), ncol(field))
    Sx <- Sy <- alpha*dt/2/dx/dx

    ## diagonal of A at once
    diag(A) <- -1-2*Sx

    ## adjacent diagonals "Sx"
    for (ii in seq(1, nrow(field)-1)){
        A[ii+1, ii] <- Sx
        A[ii, ii+1] <- Sx
    }

    B <- numeric(ncol(field))

    ## We now distinguish the top and bottom rows
    if (i == 1){
        ## top boundary, "i-1" doesn't exist or is at a ghost
        ## node/cell boundary (TODO)
        for (ii in seq_along(B))
            B[ii] <- (-1 +2*Sy)*field[i,ii] - Sy*field[i+1,ii]
    } else if (i == nrow(field)){ 
        ## bottom boundary, "i+1" doesn't exist or is at a ghost
        ## node/cell boundary (TODO)
        for (ii in seq_along(B))
            B[ii] <- -Sy*field[i-1, ii] + (-1 +2*Sy)*field[i,ii]
        
    } else {
        ## inner grid row, full expression
        for (ii in seq_along(B))
            B[ii] <- -Sy*field[i-1, ii] + (-1 +2*Sy)*field[i,ii] - Sy*field[i+1,ii]
    }
    
    x <- solve(A, B)
    x
}



DoRef <- function(n, alpha, dt, iter) {
    require(ReacTran)
    require(deSolve)
    
    N     <- n          # number of grid cells
    XX    <- n          # total size
    dy    <- dx <- XX/N  # grid size
    Dy    <- Dx <- alpha  # diffusion coeff, X- and Y-direction
    
    ## The model equations
    
    Diff2D <- function (t, y, parms)  {
        CONC  <- matrix(nrow = N, ncol = N, y)
        dCONC <- tran.2D(CONC, D.x = Dx, D.y = Dy, dx = dx, dy = dy)$dC
        return (list(dCONC))
    }
    
    ## initial condition: 0 everywhere, except in central point
    y <- matrix(nrow = N, ncol = N, data = 0)
    cen <- ceiling(N/2)
    y[cen, cen] <- 1  ## initial concentration in the central point...
    

    ## solve for  time units
    times <- seq(0,iter)*dt
    
    out <- ode.2D (y = y, func = Diff2D, t = times, parms = NULL,
                   dim = c(N,N), lrw=155412)

    ref <- matrix(out[length(times),-1], N, N)
    return(ref)
}

## test number 1
adi1 <- ADI(n=25, dt=100, iter=50, alpha=1E-3)
ref1 <- DoRef(n=25, alpha=1E-3, dt=100, iter=50)

plot(adi1[[length(adi1)]], ref1, log="xy", xlab="ADI", ylab="ode.2D (reference)",
     las=1, xlim=c(1E-15, 1), ylim=c(1E-15, 1))
abline(0,1)

sapply(adi1, sum)

## test number 2
adi2 <- ADI(n=51, dt=10, iter=200, alpha=1E-3)
ref2 <- DoRef(n=51, alpha=1E-3, dt=10, iter=200)

plot(adi2[[length(adi2)]], ref2, log="xy", xlab="ADI", ylab="ode.2D (reference)",
     las=1, xlim=c(1E-15, 1), ylim=c(1E-15, 1))
abline(0,1)

