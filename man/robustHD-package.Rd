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
Version: \tab 0.2.2\cr
Date: \tab 2012-09-12\cr
Depends: \tab R (>= 2.14.1), Rcpp (>= 0.9.10), ggplot2 (>= 0.9.2), MASS,
parallel, pcaPP, perry, robustbase (>= 0.6-0)\cr
Imports: \tab RcppArmadillo (>= 0.3.0), stats\cr
LinkingTo: \tab Rcpp, RcppArmadillo\cr
Suggests: \tab mvtnorm\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
AIC.sparseLTS           Information criteria for sparse LTS regression
                        models
coef.seqModel           Extract coefficients from a sequence of
                        regression models
coef.sparseLTS          Extract coefficients from sparse LTS regression
                        models
coefPlot                Coefficient plot of a sequence of regression
                        models
corHuber                Robust correlation based on winsorization.
critPlot                Optimality criterion plot of a sequence of
                        regression models
diagnosticPlot          Diagnostic plots for sparse LTS regression
                        models
fitted.seqModel         Extract fitted values from a sequence of
                        regression models
fitted.sparseLTS        Extract fitted values from sparse LTS
                        regression models
fortify.sparseLTS       Convert a sparse LTS regression model into a
                        data frame for plotting
lambda0                 Penalty parameter for sparse LTS regression
perry.seqModel          Resampling-based prediction error for a
                        sequential regression model
perry.sparseLTS         Resampling-based prediction error for sparse
                        LTS regression models
plot.seqModel           Plot a sequence of regression models
predict.seqModel        Predict from a sequence of regression models
predict.sparseLTS       Predict from sparse LTS regression models
residuals.seqModel      Extract residuals from a sequence of regression
                        models
residuals.sparseLTS     Extract residuals from sparse LTS regression
                        models
rlars                   Robust least angle regression
robustHD-package        Robust methods for high-dimensional data
sparseLTS               Sparse least trimmed squares regression
sparseLTSGrid           Sparse least trimmed squares regression
standardize             Data standardization
winsorize               Data cleaning by winsorization
wt                      Extract outlier weights from sparse LTS
                        regression models
}
}
\author{
Andreas Alfons [aut, cre]

Maintainer: Andreas Alfons <andreas.alfons@kuleuven.be>
}
\keyword{package}
