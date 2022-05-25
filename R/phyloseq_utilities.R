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
  data.frame(tax = phyloseq::taxa_names(phy),
             mra = rowMeans(phyloseq::otu_table(phyloseq::transform_sample_counts(phy, function(x) x/sum(x)))),
             prevalence = rowMeans(phyloseq::otu_table(phy) > 0),
             as.data.frame(methods::as(phyloseq::tax_table(phy), "matrix")))
}
