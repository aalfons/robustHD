# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## prepare BIC values for plotting
# model .... object of class "bicSelect"
# data ..... data frame containing additional information such as step numbers 
#            or values of a tuning parameter
# select ... indicates columns of a BIC matrix to keep
fortify.bicSelect <- function(model, data = NULL, select = NULL, ...) {
  # extract BIC values and make sure they are in matrix form
  bic <- model$values
  if(!is.null(dim(bic)) && !is.null(select)) bic <- bic[, select, drop=FALSE]
  bic <- as.data.frame(bic)
  # define index of the submodels
  d <- dim(bic)
  # check data frames for BIC and additional information
  if(any(d == 0)) stop("BIC data has no rows or columns")
  if(is.null(data)) data <- bic[, 0, drop=FALSE]  # NULL data frame
  else if(!is.data.frame(data) || nrow(data) != d[1]) {
    stop(sprintf("'data' must be a data frame with %d rows", d[1]))
  }
  # combine all information into one data frame
  index <- seq_len(d[1])
  if(d[2] == 1) {
    bic <- cbind(index, data, bic)
    names(bic) <- c("index", names(data), "BIC")
  } else {
    # reshape BIC matrix and add additional information (works if data is NULL)
    bic <- mapply(function(val, nam) {
      cbind(fit=rep.int(nam, d[1]), index=index, data, BIC=val)
    }, val=bic, nam=names(bic), SIMPLIFY=FALSE, USE.NAMES=FALSE)
    bic <- do.call(rbind, bic)
    attr(bic, "facets") <- . ~ fit
  }
  # return data frame
  bic
}


# fortify.bicSelect <- function(model, data = NULL, select = NULL, ...) {
#   haveVector <- is.null(dim(bic <- model$values))
#   bic <- if(haveVector) data.frame(BIC=bic) else as.data.frame(bic)
#   n <- nrow(bic)
#   if(is.null(data)) data <- data.frame(Index=seq_len(n))
#   else {
#     data <- if(is.null(dim(data))) data.frame(x=data) else as.data.frame(data)
#     if(nrow(data) != n) stop(sprintf("'data' must have %d rows", n))
#   }
#   bic <- mapply(function(v, n) {
#     data.frame(data, BIC=v, Fit=rep.int(n, length(v)))
#   }, bic, names(bic), SIMPLIFY=FALSE, USE.NAMES=FALSE)
#   bic <- do.call(rbind, bic)
#   
#   
#   
#   bic <- model$values
#   bic <- if(is.null(dim(bic))) data.frame(BIC=bic) else as.data.frame(bic)
#   if(is.null(data) {
#     data <- 
#     } else data <- as.data.frame(data)
#   bic <- mapply(function(v, n) {
#     data.frame(Fit=rep.int(n, length(v)), BIC=v, data)
#     }, bic, names(bic), SIMPLIFY=FALSE, USE.NAMES=FALSE)
#   bic <- do.call(rbind, bic)
#   if(!is.null(data)) 
# }
# 
# 
# bicPlot <- function(object, x = NULL, size = c(0.5, 2), main = NULL, 
#                     xlab, ylab, ..., mapping, data) {
#   
# }
# 
# 
# plot.bicSelect <- function()
# 
# #' @S3method plot bicSelect
# plot.bicSelect <- function(critData, abscissa = c("step", "lambda"), 
#                            size = c(0.5, 2), main = NULL, xlab, ylab, ..., mapping, data) {
#   # initializations
#   abscissa <- match.arg(abscissa)
#   crit <- setdiff(names(critData), c("fit", "step", "lambda"))
#   size <- as.numeric(size)
#   size <- c(size, rep.int(NA, max(0, 2-length(size))))[1:2]  # ensure length 2
#   size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
#   # define default axis labels
#   if(missing(xlab)) xlab <- switch(abscissa, step="Step", lambda="lambda")
#   if(missing(ylab)) ylab <- crit
#   # define aesthetic mapping for plotting coefficients
#   mapping <- aes_string(x=abscissa, y=crit)
#   # draw minor grid lines for each step, but leave 
#   # major grid lines and tick marks pretty
#   gridX <- unique(critData[, abscissa])
#   # create plot
#   ggplot(critData, mapping) + 
#     geom_line(size=size[1], ...) + 
#     geom_point(size=size[2], ...) + 
#     scale_x_continuous(minor_breaks=gridX) + 
#     labs(title=main, x=xlab, y=ylab)
# }
