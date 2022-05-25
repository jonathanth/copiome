
#' Compute relative abundances
#'
#' @param x A numeric vector of counts
#'
#' @return A numeric vector of relative abundances
#' @export
#'
#' @examples
#' transform_sample_counts(phy, relabu)
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
#' library(tidyverse)
#' glm(factor(am) ~ disp, data = mtcars, family = binomial) %>% tidylog(exp = TRUE, conf.int = TRUE)
#'
tidylog <- function(model, ...){
  require(broom)
  tidy(model, ...) %>%
    mutate(n = nobs(model), n_missing = nrow(model$data) - n)
}
