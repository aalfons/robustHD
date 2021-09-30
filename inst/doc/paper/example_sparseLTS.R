# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# load package and data
library("robustHD")
data("nci60")  # contains matrices 'protein' and 'gene'

# define response variable
y <- protein[, 92]
# screen most correlated predictor variables
correlations <- apply(gene, 2, corHuber, y)
keep <- partialOrder(abs(correlations), 100, decreasing = TRUE)
X <- gene[, keep]

# fit sparse least trimmed squares regression and print results
lambda <- seq(0.01, 0.5, length.out = 10)
fit <- sparseLTS(X, y, lambda = lambda, mode = "fraction", crit = "PE",
                 splits = foldControl(K = 5, R = 1), seed = 20210507)
fit


# create optimality criterion plot
p1 <- critPlot(fit) +
  labs(title = "Optimality criterion plot")

# create diagnostic plot of optimal model fit
p2 <- diagnosticPlot(fit, which = "rdiag", id.n = 0) +
  labs(title = "Regression diagnostic plot") +
  theme(legend.position = "top", legend.title = element_blank())


## arrange plots in file
library("gridExtra")
# pdf file
pdf(file = "figure_sparseLTS.pdf", width = 6.5, height = 3.5)
grid.arrange(p1, p2, nrow = 1)
dev.off()
# png file
png(file = "figure_sparseLTS.png", width = 6.5*72, height = 3.5*72)
grid.arrange(p1, p2, nrow = 1)
dev.off()
