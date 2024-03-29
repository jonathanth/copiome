% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_utilities.R
\name{myadonis}
\alias{myadonis}
\title{PERMANOVA analysis easy set-up for block and non-block designs.}
\usage{
myadonis(
  dist,
  phy,
  var,
  covars = c(),
  interaction = NA,
  formula = NA,
  nproc = 1,
  seed = 123,
  nperm = 999,
  verbose = TRUE,
  by = "margin",
  logvars = c(),
  block = NA,
  transient = TRUE
)
}
\arguments{
\item{dist}{Distance matrix.}

\item{phy}{Phyloseq object.}

\item{var}{Variable outcome of interest.}

\item{covars}{Vector of covariates to adjust for.}

\item{interaction}{Interaction term.}

\item{formula}{User-specified formula. If not provided, `dist ~ var (+ covars)`.}

\item{nproc}{Number of parallel processors available to use.}

\item{seed}{Random seed. Default is `123`.}

\item{nperm}{Number of permutations required.}

\item{verbose}{Printed output. Default is TRUE.}

\item{by}{by = "terms" will assess significance for each term (sequentially from first to last), by = "margin" will assess the marginal effects of the terms (each marginal term analysed in a model with all other variables), and by = NULL will assess the overall significance of all terms together. Default is "margin".}

\item{logvars}{Vector of (co)variables to log-transform. Must be included in `vars` or `covars`, too.}

\item{block}{Name of group (strata) within which to constrain permutations.}

\item{transient}{Default TRUE. If the variable outcome is a transient phenotype, when `TRUE` the groups evaluated are transient vs. cross-sectional diseased; when `FALSE` the groups evaluated are transient vs. healthy.}
}
\value{
Returns an anova.cca result as data.frame object with a new column for total observations (n) and missing (n_missing).
}
\description{
PERMANOVA analysis easy set-up for block and non-block designs.
}
\examples{
library(phyloseq)
data(GlobalPatterns)
sample_data(GlobalPatterns)$contaminated <- sample(c(rep(0, length = 8), rep(1, length = 18)))
dist <- distance(transform_phy(GlobalPatterns, transform = "log", pseudocount = 1), method = "bray")

## non-block design
myadonis(dist, GlobalPatterns, var = "contaminated", covars = c("SampleType"))

## with block design
block <- "sampletype"
myadonis(dist, GlobalPatterns, var = "contaminated", block = "SampleType")
}
