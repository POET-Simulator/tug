## Time-stamp: "Last modified 2023-01-05 17:52:55 delucia"

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
        out[[it]] <- res
    }

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
    if (i == 1) {
        ## top boundary, "i-1" doesn't exist or is at a ghost
        ## node/cell boundary (TODO)
        for (ii in seq_along(B))
            B[ii] <- (-1 +2*Sy)*field[i,ii] - Sy*field[i+1,ii]
    } else if (i == nrow(field)) { 
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


## Test heterogeneous scheme, chain rule
ADIHet <- function(field, dt, iter, alpha) {

    if (!all.equal(dim(field), dim(alpha)))
        stop("field and alpha are not matrix")
    
    ## now both field and alpha must be nx*ny matrices
    nx <- ncol(field)
    ny <- nrow(field)
    dx <- dy <- 1

    ## find out the center of the grid to apply conc=1
    cenx <- ceiling(nx/2)
    ceny <- ceiling(ny/2)
    field[cenx, ceny] <- 1

    Aij <- Bij <- alpha

    for (i in seq(2,ncol(field)-1)) {
        for (j in seq(2,nrow(field)-1)) {
            Aij[i,j] <- (alpha[i+1,j]-alpha[i-1,j])/4 + alpha[i,j]
            Bij[i,j] <- (alpha[i,j+1]-alpha[i,j-1])/4 + alpha[i,j]
        }
    }

    if (any(Aij<0) || any(Bij<0))
        stop("Aij or Bij are negative!")

    ## prepare containers for computations and outputs
    tmpX <- tmpY <- res <- field
    out <- vector(mode="list", length=iter)
    
    for (it in seq(1, iter)) {
        for (i in seq(1, ny))
            tmpX[i,] <- SweepByRowHet(i, res, dt=dt, alpha=alpha, Aij, Bij)

        resY <- t(tmpX)
        for (i in seq(1, nx))
            tmpY[i,] <- SweepByRowHet(i, resY, dt=dt, alpha=alpha, Bij, Aij)

        res <- t(tmpY)
        out[[it]] <- res
    }

    return(out)
}


## Workhorse function to fill A, B and solve for a given *row* of the
## grid matrix
SweepByRowHet <- function(i, field, dt, alpha, Aij, Bij) {
    dx <- 1 ## fixed in our test
    Sx <- Sy <- dt/2/dx/dx
    
    ## diagonal of A at once
    A <- matrix(0, nrow(field), ncol(field))
    diag(A) <- 1+2*Sx*diag(alpha)

    ## adjacent diagonals "Sx"
    for (ii in seq(1, nrow(field)-1)) {
        A[ii+1, ii] <- -Sx*Aij[ii+1,ii]
        A[ii, ii+1] <- -Sx*Aij[ii,ii+1]
    }

    B <- numeric(ncol(field))

    ## We now distinguish the top and bottom rows
    if (i == 1) {
        ## top boundary, "i-1" doesn't exist or is at a ghost
        ## node/cell boundary (TODO)
        for (ii in seq_along(B))
            B[ii] <- Sy*Bij[i+1,ii]*field[i+1,ii] + (1-2*Sy*Bij[i,ii])*field[i, ii]
    } else if (i == nrow(field)) { 
        ## bottom boundary, "i+1" doesn't exist or is at a ghost
        ## node/cell boundary (TODO)
        for (ii in seq_along(B))
            B[ii] <- (1-2*Sy*Bij[i,ii])*field[i, ii] + Sy*Bij[i-1,ii]*field[i-1,ii]
        
    } else {
        ## inner grid row, full expression
        for (ii in seq_along(B))
            B[ii] <- Sy*Bij[i+1,ii]*field[i+1,ii] + (1-2*Sy*Bij[i,ii])*field[i, ii] + Sy*Bij[i-1,ii]*field[i-1,ii]
    }
    
    x <- solve(A, B)
    x
}

## adi2 <- ADI(n=51, dt=10, iter=200, alpha=1E-3)
## ref2 <- DoRef(n=51, alpha=1E-3, dt=10, iter=200)

n <- 51
field <- matrix(0, n, n)
alphas <- matrix(1E-3*runif(n*n, 1,1.2), n, n)

## for (i in seq(1,nrow(alphas)))
##     alphas[i,] <- seq(1E-7,1E-3, length=n)

#diag(alphas) <- rep(1E-2, n)

adih1 <- ADIHet(field=field, dt=10, iter=100, alpha=alphas)
adi2  <- ADI(n=n, dt=10, iter=100, alpha=1E-3)


par(mfrow=c(1,3))
image(adi2[[length(adi2)]])
image(adih1[[length(adih1)]])
points(0.5,0.5, col="red",pch=4)
plot(adih1[[length(adih1)]], adi2[[length(adi2)]], pch=4, log="xy")
abline(0,1)


sapply(adih1, sum)
sapply(adi2, sum)

adi2


par(mfrow=c(1,2))
image(alphas)
image(adih1[[length(adih1)]])
points(0.5,0.5, col="red",pch=4)




## Test heterogeneous scheme, direct discretization
ADIHetDir <- function(field, dt, iter, alpha) {

    if (!all.equal(dim(field), dim(alpha)))
        stop("field and alpha are not matrix")
    
    ## now both field and alpha must be nx*ny matrices
    nx <- ncol(field)
    ny <- nrow(field)
    dx <- dy <- 1

    ## find out the center of the grid to apply conc=1
    cenx <- ceiling(nx/2)
    ceny <- ceiling(ny/2)
    field[cenx, ceny] <- 1

    ## prepare containers for computations and outputs
    tmpX <- tmpY <- res <- field
    out <- vector(mode="list", length=iter)
    
    for (it in seq(1, iter)) {
        for (i in seq(2, ny-1)) {
            Aij <- cbind(harm(alpha[i,], alpha[i-1,]), harm(alpha[i,], alpha[i+1,]))
            Bij <- cbind(harm(alpha[,i], alpha[,i-1]), harm(alpha[,i], alpha[,i+1]))
            tmpX[i,] <- SweepByRowHetDir(i, res, dt=dt, Aij, Bij)
        }
        resY <- t(tmpX)
        for (i in seq(2, nx-1))
            tmpY[i,] <- SweepByRowHetDir(i, resY, dt=dt, Bij, Aij)
        res <- t(tmpY)
        out[[it]] <- res
    }

    return(out)
}

harm <- function(x,y) 1/(1/x+1/y)

harm(1,4)


## Direct discretization, Workhorse function to fill A, B and solve
## for a given *row* of the grid matrix
SweepByRowHetDir <- function(i, field, dt, Aij, Bij) {
    dx <- 1 ## fixed in our test
    Sx <- Sy <- dt/2/dx/dx
    
    ## diagonal of A at once
    A <- matrix(0, nrow(field), ncol(field))
    diag(A) <- 1 + Sx*(Aij[,1]+Aij[,2])
    
    ## adjacent diagonals "Sx"
    for (ii in seq(1, nrow(field)-1)) {
        A[ii+1, ii] <- -Sx*Aij[ii,2] # i-1/2
        A[ii, ii+1] <- -Sx*Aij[ii,1] # i+1/2
    }

    B <- numeric(ncol(field))

    for (ii in seq_along(B))
        B[ii] <- Sy*Bij[ii,2]*field[i+1,ii] + (1 - Sy*(Bij[ii,1]+Bij[ii,2]))*field[i, ii] + Sy*Bij[ii,1]*field[i-1,ii]

    lastA <<- A
    lastB <<- B
    x <- solve(A, B)
    x
}

## adi2 <- ADI(n=51, dt=10, iter=200, alpha=1E-3)
## ref2 <- DoRef(n=51, alpha=1E-3, dt=10, iter=200)

n <- 5
field <- matrix(0, n, n)
alphas <- matrix(1E-5*runif(n*n, 1,2), n, n)

alphas1 <- matrix(3E-5, n, 25)
alphas2 <- matrix(1E-5, n, 26)

alphas <- cbind(alphas1, alphas2)

## for (i in seq(1,nrow(alphas)))
##     alphas[i,] <- seq(1E-7,1E-3, length=n)

#diag(alphas) <- rep(1E-2, n)

adih  <- ADIHetDir(field=field, dt=10, iter=200, alpha=alphas)
adi2  <- ADI(n=n, dt=10, iter=200, alpha=1E-5)


par(mfrow=c(1,3))
image(adi2[[length(adi2)]])
image(adih[[length(adih)]])
points(0.5,0.5, col="red",pch=4)
plot(adih[[length(adih)]], adi2[[length(adi2)]], pch=4, log="xy")
abline(0,1)


sapply(adih, sum)
sapply(adi2, sum)

adi2


par(mfrow=c(1,2))
image(alphas)
points(0.5,0.5, col="red",pch=4)
image(adih[[length(adih)]])
points(0.5,0.5, col="red",pch=4)

options(width=110)
