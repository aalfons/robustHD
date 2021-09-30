# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# load package and data
library("robustHD")
data("TopGear")

# keep complete observations and remove information on car model
keep <- complete.cases(TopGear)
TopGear <- TopGear[keep, -(1:3)]
# log-transform price
TopGear$Price <- log(TopGear$Price)

# fit robust groupwise least angle regression and print results
fit <- rgrplars(MPG ~ ., data = TopGear, sMax = 15,
                crit = "BIC", seed = 20210507)
fit


# create coefficient plot of sequence of model fits
p1 <- coefPlot(fit) +
  labs(title = "Coefficient plot") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.16)))

# create optimality criterion plot
p2 <- critPlot(fit) +
  labs(title = "Optimality criterion plot")

# create diagnostic plot of optimal model fit
p3 <- diagnosticPlot(fit, covArgs = list(alpha = 0.8),
                     which = "rdiag", id.n = 0) +
  labs(title = "Regression diagnostic plot") +
  theme(legend.position = "top", legend.title = element_blank())


## arrange plots in file
library("gridExtra")
library("svglite")
# # pdf file
pdf(file = "figure_rgrplars.pdf", width = 9.75, height = 3.5)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()
# svg file
svglite(file = "figure_rgrplars.svg", width = 9.75, height = 3.5)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()
