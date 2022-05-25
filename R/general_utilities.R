
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

#' Tidy function for logistic regression
#' Adds n and n_missing
#'
#' @param model The glm model object
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
                n_missing = nrow(model$data) - n)
}
