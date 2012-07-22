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

# get function to call C++ back end
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
