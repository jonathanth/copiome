% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_utilities.R
\name{tidycox}
\alias{tidycox}
\title{Tidy function for cox regression
Adds n and n_missing}
\usage{
tidycox(model, ...)
}
\arguments{
\item{model}{The coxph model object}

\item{...}{Other arguments passed on to tidy()}
}
\value{
A tidy data.frame
}
\description{
Tidy function for cox regression
Adds n and n_missing
}
\examples{
library(survival)
data("cancer")
coxmodel <- coxph(Surv(time, status) ~ sex + age, data = cancer)
tidycox(coxmodel)
}
