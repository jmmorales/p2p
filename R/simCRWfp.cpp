#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame simCRWfp(int nind, double m, double xp, double yp, double rs, double per, double mscale, double mshape, 
double kappa, double meanspeed){ 
  
  double pi = 3.14159265358979323846;
  Function rweibull("rweibull");
  Function rvonmises("rvonmises");
  Function rexp("rexp");
  Function runif("runif");
  RNGScope scope;    // ensure RNG gets set/reset

// xp = 0.2   # patch coordinates
// yp = 0
// rs <- 0.1  # patch radius

NumericVector res(nind);
NumericVector resT(nind);
// m <- 0.00001   # death rate

NumericVector nm(nind);
nm = ceil(rexp(nind,m)) + 1;   // number of moves


for(int i=0; i<nind; i++){
  NumericVector R(nm[i]);
  NumericVector Tn(nm[i]);
  R = rweibull(nm[i], mshape, mscale);
  Tn = rvonmises(nm[i], 0, kappa);
  NumericVector bearing = runif(nm[i],0,2*pi);             // initial movement direction
  NumericVector x(nm[i]);
  NumericVector y(nm[i]);
  double cT = 0;
  double dst;
  //res[i] =0;
  
  for(int j=1; j<nm[i]; j++){
    bearing[j] = bearing[j-1] + Tn[j]; // fmod((bearing + Tn[j]), (2*pi)); // # keep directions between 0 and 2*pi;
    x[j] = x[j-1] + R[j]*cos(bearing[j]);
    y[j] = y[j-1] + R[j]*sin(bearing[j]);
    cT = cT + meanspeed*R[j];
    dst = sqrt(pow(x[j]-xp,2) + pow(y[j]-yp,2)) - (rs + per);
    //Rcout << dst;
    if(dst < 0){
      res[i] = 1;
      resT[i] = cT;
      break;
      }
  }
}

return DataFrame::create(Named("res")= res, Named("resT")= resT, Named("nm")=nm);
}