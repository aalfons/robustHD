# robustHD: Robust Methods for High-Dimensional Data


## Summary

In regression analysis with high-dimensional data, variable selection is an important step to (i) overcome computational problems, (ii) improve prediction performance by variance reduction, and (iii) increase interpretability of the resulting models due to the smaller number of variables.  However, robust methods are necessary to prevent outlying data points from distorting the results.  The add-on package `robustHD` for the statistical computing environment `R` provides functionality for robust linear model selection with high-dimensional data.  More specifically, the implemented functionality includes robust least angle regression ([Khan et al., 2007](https://doi.org/10.1198/016214507000000950)), robust groupwise least angle regression ([Alfons et al., 2016](https://doi.org/10.1016/j.csda.2015.02.007)), as well as sparse least trimmed squares regression [Alfons et al., 2013](https://doi.org/10.1214/12-AOAS575). The latter can be seen as a trimmed version of the popular lasso regression estimator ([Tibshirani, 1996](https://doi.org/10.1111/j.2517-6161.1996.tb02080.x)).  Selecting the optimal model can be done via cross-validation or an information criterion, and various plots are available to illustrate model selection and to evaluate the final model estimates.  Furthermore, the package includes functionality for pre-processing and cleaning the data, such as robust standardization and winsorization.  Finally, `robustHD` follows a clear object-oriented design and takes advantage of C++ code and parallel computing to reduce computing time. 


## Main functionality

 * `sparseLTS()`: Sparse least trimmed squares regression.
 
 * `rlars()`: Robust least angle regression.
 
 * `grplars()` and `rgrplars()`: (Robust) groupwise least angle regression.
 
 * `tslars()` and `rtslars()`: (Robust) least angle regression for time series data.
 
 * `corHuber()`: Robust correlation based on winsorization.
 
 * `winsorize()`: Data cleaning by winsorization.
 
 * `robStandardize()`: Data standardization with given functions for computing center and scale. By default, the median and MAD are used.


## Installation

Package `robustHD` is on CRAN (The Comprehensive R Archive Network), hence the latest release can be easily installed from the `R` command line via

```
install.packages("robustHD")
```


## Building from source

To install the latest (possibly unstable) development version from GitHub, you can pull this repository and install it from the `R` command line via

```
install.packages("devtools")
devtools::install_github("aalfons/robustHD")
```

If you already have package `devtools` installed, you can skip the first line.  Moreover, package `robustHD` contains `C++` code that needs to be compiled, so you may need to download and install the [necessary tools for MacOS](https://cran.r-project.org/bin/macosx/tools/) or the [necessary tools for Windows](https://cran.r-project.org/bin/windows/Rtools/).


## Community guidelines

### Report issues and request features

If you experience any bugs or issues or if you have any suggestions for additional features, please submit an issue via the *Issues* tab of this repository.  Please have a look at existing issues first to see if your problem for feature request has already been discussed.

### Contribute to the package

For modifying the code in accordance with the GNU GPL, you can fork this repository.  In order to contribute to the package, you can create a pull request after forking the repository and implementing the relevant functionality.

### Ask for help

If you need help using the package, or if you are interested in collaborations related to this project, please get in touch: alfons at ese dot eur dot nl


## References

Alfons, A., Croux, C. and Gelper, S. (2013) Sparse least trimmed squares regression for analyzing high-dimensional large data sets. The Annals of Applied Statistics, 7(1), 226–248. DOI [10.1214/12-AOAS575](https://doi.org/10.1214/12-AOAS575).

Alfons, A., Croux, C. and Gelper, S. (2016) Robust groupwise least angle regression. Computational Statistics & Data Analysis, 93, 421–435. DOI [10.1016/j.csda.2015.02.007](https://doi.org/10.1016/j.csda.2015.02.007).

Khan, J.A., Van Aelst, S. and Zamar, R.H. (2007) Robust linear model selection based on least angle regression. Journal of the American Statistical Association, 102(480), 1289–1299. DOI [10.1198/016214507000000950](https://doi.org/10.1198/016214507000000950).

Tibshirani, R. (1996) Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society, Series B, 58(1), 267–288. DOI [10.1111/j.2517-6161.1996.tb02080.x](https://doi.org/10.1111/j.2517-6161.1996.tb02080.x).
