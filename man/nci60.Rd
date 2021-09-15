\name{nci60}
\alias{nci60}
\alias{protein}
\alias{X}
\docType{data}
\title{
NCI-60 cancer cell panel
}
\description{
The data set is a pre-processed version of the NCI-60 cancer cell panel as used
in Alfons, Croux & Gelper (2013).  One observation was removed since all values
in the gene expression data were missing.
}
\usage{data("nci60")}
\format{
  Protein and gene expression data on 59 observations are stored in two
  separate matrices:
  \describe{
    \item{\code{protein}}{a matrix containing protein expressions based on
    antibodies (162 variables), acquired via reverse-phase protein lysate
    arrays and log2 transformed.}
    \item{\code{X}}{a matrix containing gene expression data (22283 variables),
    obtained with an Affymetrix HG-U133A chip and normalized with the GCRMA
    method.}
  }
}
\source{
The original data were downloaded from
\url{https://discover.nci.nih.gov/cellminer/} on 2012-01-27.  They can be
obtained from \url{https://github.com/aalfons/nci60}, together with our script
for pre-processing.
}
\examples{
# load data
data("nci60")
# define response variable
y <- protein[, 92]
# screen most correlated predictor variables
correlations <- apply(X, 2, corHuber, y)
keep <- partialOrder(abs(correlations), 100, decreasing = TRUE)
X <- X[, keep]
}
\references{
Reinhold, W.C., Sunshine, M., Liu, H., Varma, S., Kohn, K.W., Morris, J.,
Doroshow, J. and Pommier, Y. (2012) CellMiner: A Web-Based Suite of Genomic and
Pharmacologic Tools to Explore Transcript and Drug Patterns in the NCI-60 Cell
Line Set. \emph{Cancer Research}, \bold{72}(14), 3499--3511.

Alfons, A., Croux, C. and Gelper, S. (2013) Sparse least trimmed squares
regression for analyzing high-dimensional large data sets. \emph{The Annals
of Applied Statistics}, \bold{7}(1), 226--248.
}
\keyword{datasets}
