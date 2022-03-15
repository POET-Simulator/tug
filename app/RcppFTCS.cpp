//  Time-stamp: "Last modified 2022-03-15 18:15:39 delucia"
#include <Rcpp.h>
#include <iostream> // for std
#include <vector>   // for vector

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector RcppFTCS(int n,
                       double length,    
                       NumericVector & field,
                       double alpha,
                       double bc_left,
                       double bc_right,
                       double timestep)
{
    // dimension of grid
    NumericVector ext (clone(field));
    double dx = length / ((double) n - 1.);
    double dt = 0.25*dx*dx/alpha;

    
    double afac = alpha*dt/dx/dx;
    int iter = (int) (timestep/dt);

    Rcout << "dt: " << dt << "; inner iterations: " << iter << endl;


    for (int it = 0; it < iter; it++){
        for (int i = 1; i < ext.size()-1; i++) {
            ext[i] = (1. - 2*afac)*ext[i] + afac*(ext[i+1]+ext[i-1]);
        }
        ext[0] = bc_left;
        ext[n-1] = bc_right;
    }
    return(ext);
}

