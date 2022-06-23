
#' Compute relative abundances
#'
#' @param x A numeric vector of counts
#'
#' @return A numeric vector of relative abundances
#' @export
#'
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' gp_ra <- transform_sample_counts(GlobalPatterns, relabu)
#' sample_sums(gp_ra)
relabu <- function(x){
  x/sum(x)
}

#' Tidy function for linear or logistic regression
#' Adds n and n_missing
#'
#' @param model The lm/glm model object
#' @param ... Other arguments passed on to tidy()
#'
#' @return A data.frame with model parameters, similar to tidy()
#' @export
#'
#' @examples
#' library(dplyr)
#' library(magrittr)
#' glm(factor(am) ~ disp, data = mtcars, family = binomial) %>% tidylog(exp = TRUE, conf.int = TRUE)
#'
tidylog <- function(model, ...){
  dplyr::mutate(broom::tidy(model, ...),
                n = stats::nobs(model),
                n_missing = length(model$na.action))
}



#' Tidy function for cox regression
#' Adds n and n_missing
#'
#' @param model The coxph model object
#' @param ... Other arguments passed on to tidy()
#'
#' @return A tidy data.frame
#' @export
#'
#' @examples
#' library(survival)
#' data("cancer")
#' coxmodel <- coxph(Surv(time, status) ~ sex + age, data = cancer)
#' tidycox(coxmodel)
tidycox <- function(model, ...){
  tab <- broom::tidy(model, ...)
  data.frame(tab, n = model$n, n_missing = length(model$na.action))
}

#' Round to exact number of digits
#' Useful when plotting, since round() cuts away zeroes at the end
#'
#' @param vec Numeric vector that you want to round
#' @param n number of decimal places
#'
#' @return Character vector with exact rounding.
#' @export
#'
#' @examples
#' round(3.9999, 3)
#' roundex(3.9999, 3)
roundex <- function(vec, n){
  sprintf(paste0("%.", n, "f"), vec)
}


#' Remove samples from distance object based on the samples of a phyloseq object
#'
#' @param d A numeric matrix, data frame or "dist" object.
#' @param phy A phyloseq object.
#'
#' @return A subsetted "dist" object.
#' @export
#' @importFrom stats as.dist
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' dist.bray <- distance(GlobalPatterns, "bray")
#' subset_dist(dist.bray, subset_samples(GlobalPatterns, SampleType == "Feces"))
subset_dist <- function(d, phy) {
  keep <- phyloseq::sample_names(phy)
  return(as.dist(as.matrix(d)[keep, keep]))
}
