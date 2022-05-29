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
#' otumat <- otu_df(GlobalPatterns)
#' head(otumat)[,1:5]
otu_df <- function(phy) {
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



#' Return tax data from phyloseq object as data.frame
#'
#' @param phy A phyloseq object.
#'
#' @return a data.frame containing tax data with taxa (rows) and ranks (columns). If the phyloseq object is
#' agglomerated to lower taxonomic resolution, the lower ranks are not shown.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' taxmat <- tax_df(GlobalPatterns)
#' head(taxmat)
tax_df = function(phy) {
  taxmat <- data.frame(phy@tax_table@.Data)
  group <- taxmat %>%
    dplyr::select_if(~ !all(is.na(.))) %>% colnames %>% utils::tail(1)
  if (tolower(group) != "species") {
    taxmat <- taxmat %>%
      .data[1:which(group == colnames(.data))[[1]]]
  }
  return(taxmat)
}
