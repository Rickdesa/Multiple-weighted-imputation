# Multiple weighted imputation

The R code contains functions which permit to do multiple imputations on single missing variables by predictive mean matching, single regression tree, bootstrap regression trees and Bayesian regression trees with weighted units.
The role of all variables must be specified with column indexes.

## Reproducible Example

Load the functions with `source("Code for weighted multiple imputations.R")` before running this example.

```r
set.seed(1234)
y <- rnorm(25)
x <- rnorm(25)
z <- rnorm(25, 0.5)
v <- runif(25)
x[sample(1:25, 3)] <- NA
dat.frame <- data.frame(y, x, z, v)

pmm_imputation(
  ampdf = dat.frame,
  donors.pmm = 4,
  no.imputs = 5,
  seed = 123,
  miss.ind = 2,
  cov.ind = c(1, 3),
  weights.ind = 4
)

cart_imputation_bayes_boot(
  ampdf = dat.frame,
  minbucket.cart = 3,
  no.imputs = 5,
  seed = 123,
  miss.ind = 2,
  cov.ind = c(1, 3),
  weights.ind = 4
)

cart_imputation_boot(
  ampdf = dat.frame,
  minbucket.cart = 3,
  no.imputs = 5,
  seed = 123,
  miss.ind = 2,
  cov.ind = c(1, 3),
  weights.ind = 4
)

cart_imputation(
  ampdf = dat.frame,
  minbucket.cart = 3,
  no.imputs = 5,
  seed = 123,
  miss.ind = 2,
  cov.ind = c(1, 3),
  weights.ind = 4
)
```
