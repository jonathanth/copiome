% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_utilities.R
\name{prevalence_df}
\alias{prevalence_df}
\title{Taxa data.frame prevalence}
\usage{
prevalence_df(phy, taxa = FALSE)
}
\arguments{
\item{phy}{A phyloseq object.}

\item{taxa}{FALSE/TRUE for taxonomic information. Default is FALSE.}
}
\value{
Prevalence data.frame in long format. freq = frequency (range 0-1) and prevalence (%). Each row is a sample feature.
}
\description{
Taxa data.frame prevalence
}
\examples{
library(phyloseq)
data(GlobalPatterns)
prev <- prevalence_df(GlobalPatterns, tax = TRUE)
head(prev)
}
