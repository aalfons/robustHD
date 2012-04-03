## generate data
library("mvtnorm")
set.seed(1234)  # for reproducibility
n <- 100  # number of observations
p <- 200  # number of variables
beta <- rep.int(c(1, 0), c(5, p-5))  # coefficients
sigma <- 0.5      # controls signal-to-noise ratio
epsilon <- 0.1    # contamination level
Sigma <- 0.5^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
x <- rmvnorm(n, sigma=Sigma)    # predictor matrix
e <- rnorm(n)                   # error terms
i <- 1:ceiling(epsilon*n)       # observations to be contaminated
e[i] <- e[i] + 5                # vertical outliers
y <- c(x %*% beta + sigma * e)  # response
x[i,] <- x[i,] + 5              # bad leverage points

## sparse LTS
# fit model
fit <- sparseLTS(x, y, lambda = 0.05, mode = "fraction")
# extract coefficients
coef(fit, zeros = FALSE)
coef(fit, fit = "both", zeros = FALSE)

## sparse LTS over a grid of values for lambda
# fit model
frac <- seq(0.25, 0.05, by = -0.05)
fitGrid <- sparseLTSGrid(x, y, lambda = frac, mode = "fraction")
# extract coefficients
coef(fitGrid, zeros = FALSE)
coef(fitGrid, fit = "both", zeros = FALSE)
coef(fitGrid, s = NULL, zeros = FALSE)
coef(fitGrid, fit = "both", s = NULL, zeros = FALSE)
