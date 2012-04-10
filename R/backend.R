# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

# access environment containing function to call back end
backendEnv <- function() {
    pos <-  match(".BackendEnv", search())
    if(is.na(pos)) { # must create it
        attach(list(), name=".BackendEnv")
        pos <- match(".BackendEnv", search())
    }
    return(pos.to.env(pos))
}

# put in back end environment
putBackendEnv <- function(x, value) {
    assign(x, value, envir=backendEnv())
}

# get from back end environment
getBackendEnv <- function (x, mode = "any") { 
    get(x, envir=backendEnv(), mode=mode, inherits=FALSE)
}

# does object exist in back end environment?
existsBackendEnv <- function (x, mode = "any") { 
    exists(x, envir=backendEnv(), mode=mode, inherits=FALSE)
}

# function to call built-in RcppArmadillo back end
.CallRobustHD <- function(.NAME, ..., PACKAGE) {
    .Call(.NAME, ..., PACKAGE="robustHD")
}

# register built-in RcppArmadillo back end
registerBackend <- function() {
    putBackendEnv("armadillo", .CallRobustHD)
}

## get function to call C++ back end
#getBackend <- function(which = c("eigen", "armadillo")) {
#    which <- match.arg(which)
#    if(which == "eigen") {
#        ## use RcppEigen back end if possible
#        # check if RcppEigen back end is installed
#        ok <- isTRUE("sparseLTSEigen" %in% .packages(all.available=TRUE))
#        # check platform (RcppEigen back end does not work with 32-bit Windows)
#        ok <- ok && !isTRUE(.Platform$OS.type == "windows" && .Platform$r_arch == "i386")
#        # if everything is ok so far, check if RcppEigen back end is loaded
#        if(ok) {
#            ok <- isTRUE("sparseLTSEigen" %in% .packages())
#            if(!ok) {
#                warning("RcppEigen back end is not loaded; ", 
#                    "using built-in RcppArmadillo back end")
#            }
#        } else {
#            warning("RcppEigen back end is not available; ", 
#                "using built-in RcppArmadillo back end")
#        }
#        # use RcppArmadillo back end if RcppEigen back end cannot be used
#        if(!ok) which <- "armadillo"
#    }
#    # set built-in Armadillo back end
#    getBackendEnv(which)
#}
getBackend <- function() {
    # if the RcppEigen back end is registered, check the platform
    # (RcppEigen back end does not work with 32-bit Windows)
    ok <- existsBackendEnv("eigen") && 
        !isTRUE(.Platform$OS.type == "windows" && .Platform$r_arch == "i386")
    # use the RcppEigen back end if everything is ok, otherwise use the 
    # built-in RcppArmadillo back end
    if(ok) {
        getBackendEnv("eigen")
    } else getBackendEnv("armadillo")
}
