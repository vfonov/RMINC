library(RMINC)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)

getRMINCTestData("./")

gf <- read.csv("rminctestdata/test_data_set.csv")

voxel_left <- mincGetVoxel(gf$jacobians_fixed_2[1:10], 0,0,0)
voxel_right <- mincGetVoxel(gf$jacobians_fixed_2[11:20], 0,0,0)
Sex <- gf$Sex[1:10]
Scale <- gf$scale[1:10]
Coil <- as.factor(gf$coil[1:10])

gf$coil <- as.factor(gf$coil)
gftest <- gf[1:10,]

src <- '
Rcpp::List fLmSEXP(SEXP Xs, SEXP ys) {
    Rcpp::NumericMatrix Xr(Xs);
    Rcpp::NumericVector yr(ys);
    int n = Xr.nrow(), k = Xr.ncol();
    arma::mat X(Xr.begin(), n, k, false);
    arma::colvec y(yr.begin(), yr.size(), false);
    int df = n - k;
    // fit model y ~ X, extract residuals
    arma::colvec coef = arma::solve(X, y);
    arma::colvec res  = y - X*coef;
    double s2 = std::inner_product(res.begin(), res.end(),
                                   res.begin(), 0.0)/df;
    // std.errors of coefficients
    arma::colvec sderr = arma::sqrt(s2 *
       arma::diagvec(arma::pinv(arma::trans(X)*X)));
    return Rcpp::List::create(Rcpp::Named("coefficients")=coef,
                              Rcpp::Named("stderr")      =sderr,
                              Rcpp::Named("df")          =df);
}
'
cppFunction(code=src, depends="RcppArmadillo") #creates fLmSEXP

pred_mat <- as.matrix(as.numeric(gftest$Sex))


benches <- benchmark(mincLm = mincLm(jacobians_fixed_2 ~ Sex, data = gftest)
                    , fLmSEXP = mincApplyRCPP(gftest$jacobians_fixed_2
                                                , function(x) fLmSEXP(pred_mat, as.matrix(x))
                                                , slab_sizes = c(10,10,1))
                    , fLmSEXP2 = mincApplyRCPP2(gftest$jacobians_fixed_2
                                                , function(x) fLmSEXP(pred_mat, as.matrix(x))
                                                , slab_sizes = c(10,10,1))
                    , lm.fit = mincApplyRCPP(gftest$jacobians_fixed_2
                                            , function(x) lm.fit(pred_mat, as.matrix(x))
                                            , slab_sizes = c(10,10,1))
                    , columns = c("test", "replications", "relative", "elapsed")
                    , order = "relative")
benches