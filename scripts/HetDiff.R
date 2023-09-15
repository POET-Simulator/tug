## Time-stamp: "Last modified 2023-07-31 16:28:48 delucia"

library(ReacTran)
library(deSolve)
options(width=114)

## harmonic mean
harm <- function(x,y) {
    if (length(x) != 1 || length(y) != 1)
        stop("x & y have different lengths")
    2/(1/x+1/y)
}

## harm(0, 1) ## 0
## harm(1, 2) ## 0


############# Providing coeffs on the interfaces
N     <- 11          # number of grid cells
ini   <- 1           # initial value at x=0
N2    <- ceiling(N/2)
L     <- 10          # domain side

## Define diff.coeff per cell, in  4 quadrants
alphas <- matrix(0, N, N) 
alphas[1:N2, 1:N2] <- 1
alphas[1:N2, seq(N2+1,N)] <- 0.1
alphas[seq(N2+1,N), 1:N2] <- 0.01
alphas[seq(N2+1,N), seq(N2+1,N)] <- 0.001

image(log10(alphas), col=heat.colors(4))

r180 <- function(x) {
    xx <- rev(x)
    dim(xx) <- dim(x)
    xx
}
mirror <- function (x) {
    xx <- as.data.frame(x)
    xx <- rev(xx)
    xx <- as.matrix(xx)
    xx
}

array_to_LaTeX <- function(arr) {
    rows <- apply(arr, MARGIN=1, paste, collapse = " & ")
    matrix_string <- paste(rows, collapse = " \\\\ ")
    return(paste("\\begin{bmatrix}", matrix_string, "\\end{bmatrix}"))
}


cat(array_to_LaTeX(mirror(r180(alphas))))



r180(alphas)
  
filled.contour(log10(alphas), col=terrain.colors(4), nlevels=4)

cmpharm <- function(x) {
    y <- c(0, x, 0)
    ret <- numeric(length(x)+1)
    for (i in seq(2, length(y))) {
        ret[i-1] <- harm(y[i], y[i-1])
    }
    ret
}

## Construction of the 2D grid
x.grid <- setup.grid.1D(x.up = 0, L = L, N = N)
y.grid <- setup.grid.1D(x.up = 0, L = L, N = N)
grid2D <- setup.grid.2D(x.grid, y.grid)

D.grid <- list()

# Diffusion on x-interfaces
D.grid$x.int <- apply(alphas, 1, cmpharm)

# Diffusion on y-interfaces
## matrix(nrow = N, ncol = N+1, data = rep(c(rep(1E-1, 50),5.E-1,rep(1., 50)), 100) )
D.grid$y.int <- t(apply(alphas, 2, cmpharm))

dx <- L/N
dy <- L/N

# The model equations
Diff2Dc <- function(t, y, parms) {
    CONC  <- matrix(nrow = N, ncol = N, data = y)
    dCONC <- tran.2D(CONC, dx = dx, dy = dy, D.grid = D.grid)$dC
    return(list(dCONC))
}

## initial condition: 0 everywhere, except in central point
y <- matrix(nrow = N, ncol = N, data = 0)
y[N2, N2] <- ini  # initial concentration in the central point...

## solve for 10 time units
times <- 0:10
outc <- ode.2D(y = y, func = Diff2Dc, t = times, parms = NULL,
                dim = c(N, N), lrw = 1860000)

outtimes <- c(0, 4, 7, 10)

## NB: assuming current working dir is "tug"
cairo_pdf("doc/images/deSolve_AlphaHet1.pdf", family="serif", width=12, height=12)
image(outc, ask = FALSE, mfrow = c(2, 2), main = paste("time", outtimes),
      legend = TRUE, add.contour = FALSE, subset = time %in% outtimes,
      xlab="",ylab="", axes=FALSE, asp=1)
dev.off()

## outc is a matrix with 11 rows and 122 columns (first column is
## simulation time);
str(outc)

## extract only the results and transpose the matrix for storage
ret <- data.matrix(t(outc[ , -1]))
rownames(ret) <- NULL

## NB: assuming current working dir is "tug"
data.table::fwrite(ret, file="scripts/gold/HetDiff1.csv", col.names=FALSE)




#################### 2D visualization

## Version of Rmufits::PlotCartCellData with the ability to fix the
## "breaks" for color coding of 2D simulations
Plot2DCellData <- function(data, grid, nx, ny, contour = TRUE,
                           nlevels = 12, breaks, palette = "heat.colors",
                           rev.palette = TRUE, scale = TRUE, plot.axes=TRUE, ...) {
    if (!missing(grid)) {
        xc <- unique(sort(grid$cell$XCOORD))
        yc <- unique(sort(grid$cell$YCOORD))
        nx <- length(xc)
        ny <- length(yc)
        if (!length(data) == nx * ny) 
            stop("Wrong nx, ny or grid")
    } else {
        xc <- seq(1, nx)
        yc <- seq(1, ny)
    }
    z <- matrix(round(data, 6), ncol = nx, nrow = ny, byrow = TRUE)
    pp <- t(z[rev(seq(1, nrow(z))), ])

    if (missing(breaks)) {
        breaks <- pretty(data, n = nlevels)
    }
    
    breakslen <- length(breaks)
    colors <- do.call(palette, list(n = breakslen - 1))
    if (rev.palette) 
        colors <- rev(colors)
    if (scale) {
        par(mfrow = c(1, 2))
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(4, 
            1))
    }
    par(las = 1, mar = c(5, 5, 3, 1))
    image(xc, yc, pp, xlab = "X [m]", ylab = "Y[m]", las = 1, asp = 1,
          breaks = breaks, col = colors, axes = FALSE, ann=plot.axes,
          ...)

    if (plot.axes) {
        axis(1)
        axis(2)
    }
    if (contour) 
        contour(unique(sort(xc)), unique(sort(yc)), pp, breaks = breaks, 
            add = TRUE)
    if (scale) {
        par(las = 1, mar = c(5, 1, 5, 5))
        PlotImageScale(data, breaks = breaks, add.axis = FALSE, 
            axis.pos = 4, col = colors)
        axis(4, at = breaks)
    }
    invisible(pp)
}

PlotImageScale <- function(z, zlim, col = grDevices::heat.colors(12), breaks, 
                           axis.pos = 1, add.axis = TRUE, ...) {
    if (!missing(breaks)) {
        if (length(breaks) != (length(col) + 1)) {
            stop("must have one more break than colour")
        }
    }
    if (missing(breaks) & !missing(zlim)) {
        breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 1))
    }
    if (missing(breaks) & missing(zlim)) {
        zlim <- range(z, na.rm = TRUE)
        zlim[2] <- zlim[2] + c(zlim[2] - zlim[1]) * (0.001)
        zlim[1] <- zlim[1] - c(zlim[2] - zlim[1]) * (0.001)
        breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 1))
    }
    
    poly <- vector(mode = "list", length(col))
    for (i in seq(poly)) {
        poly[[i]] <- c(breaks[i], breaks[i + 1], breaks[i + 1], 
            breaks[i])
    }
    if (axis.pos %in% c(1, 3)) {
        ylim <- c(0, 1)
        xlim <- range(breaks)
    }
    if (axis.pos %in% c(2, 4)) {
        ylim <- range(breaks)
        xlim <- c(0, 1)
    }
    plot(1, 1, t = "n", ylim = ylim, xlim = xlim, axes = FALSE, 
        xlab = "", ylab = "", xaxs = "i", yaxs = "i", ...)
    for (i in seq(poly)) {
        if (axis.pos %in% c(1, 3)) {
            polygon(poly[[i]], c(0, 0, 1, 1), col = col[i], border = NA)
        }
        if (axis.pos %in% c(2, 4)) {
            polygon(c(0, 0, 1, 1), poly[[i]], col = col[i], border = NA)
        }
    }
    box()
    if (add.axis) {
        axis(axis.pos)
    }
}

cairo_pdf("AlphaHet1.pdf", family="serif", width=8)
par(mar = c(1,1,1,1))
Plot2DCellData(log10(mirror(alphas)), nx=N, ny=N, nlevels=5, palette = terrain.colors, contour=FALSE, plot.axes=FALSE,
               scale = F,
               main=expression(log[10](alpha)))
text(3,8,"1")
text(8,8,"0.1")
text(3,3,"0.01")
text(8,3,"0.001")
# title("Diff. Coefficients (log10)")
dev.off()

