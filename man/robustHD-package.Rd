\name{robustHD-package}
\alias{robustHD-package}
\alias{robustHD}
\docType{package}
\title{
Robust methods for high-dimensional data
}
\description{
Robust methods for high-dimensional data, in particular linear
    model selection techniques based on least angle regression and sparse
    regression.
}
\details{
\tabular{ll}{
Package: \tab robustHD\cr
Type: \tab Package\cr
Version: \tab 0.5.1\cr
Date: \tab 2014-12-16\cr
Depends: \tab R (>= 3.0.2), ggplot2 (>= 0.9.2), perry (>= 0.2.0), robustbase
(>= 0.9-5)\cr
Imports: \tab MASS, parallel, stats\cr
LinkingTo: \tab Rcpp (>= 0.9.10), RcppArmadillo (>= 0.3.0)\cr
Suggests: \tab mvtnorm\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
AIC.seqModel            Information criteria for a sequence of
                        regression models
TopGear                 Top Gear car data
coef.seqModel           Extract coefficients from a sequence of
                        regression models
coefPlot                Coefficient plot of a sequence of regression
                        models
corHuber                Robust correlation based on winsorization.
critPlot                Optimality criterion plot of a sequence of
                        regression models
diagnosticPlot          Diagnostic plots for a sequence of regression
                        models
fitted.seqModel         Extract fitted values from a sequence of
                        regression models
fortify.seqModel        Convert a sequence of regression models into a
                        data frame for plotting
getScale                Extract the residual scale of a robust
                        regression model
grplars                 (Robust) groupwise least angle regression
lambda0                 Penalty parameter for sparse LTS regression
perry.seqModel          Resampling-based prediction error for a
                        sequential regression model
plot.seqModel           Plot a sequence of regression models
predict.seqModel        Predict from a sequence of regression models
residuals.seqModel      Extract residuals from a sequence of regression
                        models
rlars                   Robust least angle regression
robustHD-deprecated     Deprecated functions in package 'robustHD'
robustHD-package        Robust methods for high-dimensional data
sparseLTS               Sparse least trimmed squares regression
standardize             Data standardization
tsBlocks                Construct predictor blocks for time series
                        models
tslars                  (Robust) least angle regression for time series
                        data
tslarsP                 (Robust) least angle regression for time series
                        data with fixed lag length
winsorize               Data cleaning by winsorization
wt                      Extract outlier weights from sparse LTS
                        regression models
}
}
\author{
Andreas Alfons [aut, cre]

Maintainer: Andreas Alfons <alfons@ese.eur.nl>
}
\keyword{package}
