#' Make a data.frame with handy taxon-wise summaries
#'
#'
#' @param phy A phyloseq object
#'
#' @return a data.frame which contains the folowing variables: tax (taxa names), mra (mean relative abundance), prevalence, and a copy of the tax table appended.
#' @export
#'
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' mradat <- make_mradat(GlobalPatterns)
#' head(mradat)
make_mradat <- function(phy){
  require(phyloseq)
  require(tidyverse)
  data.frame(tax = taxa_names(phy),
             mra = rowMeans(phy %>% transform_sample_counts(function(x) x/sum(x)) %>%
                              otu_table),
             prevalence = rowMeans(otu_table(phy) > 0),
             tax_table(phy) %>% as("matrix") %>% as.data.frame)
}
