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

#' Return OTU table from phyloseq object as data.frame
#'
#' @param phy A phyloseq object.
#'
#' @return a data.frame containing abundances with samples (rows) and taxa(columns).
#' @export
#'
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' otumat <- get_otu(GlobalPatterns)
#' head(otumat)[,1:5]
get_otu <- function(phy) {
  otu <- phyloseq::otu_table(phy)
  if (phyloseq::taxa_are_rows(otu)) {
    otu <- t(otu)
  }
  return(data.frame(methods::as(otu, "matrix")))
}



#' Return sample data from phyloseq object as data.frame.
#'
#' @param phy A phyloseq object.
#'
#' @return a data.frame containing sample data with samples (rows) and metadata (columns).
#' @export
#'
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' samdat <- sample_df(GlobalPatterns)
#' head(samdat)[,1:5]
sample_df <- function(phy) {
  return(methods::as(phyloseq::sample_data(phy), "data.frame"))
}
