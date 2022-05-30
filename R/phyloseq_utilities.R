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


#' Return abundance long-formatted data.frame from phyloseq object.
#'
#' @param phy A phyloseq object
#' @param tax FALSE/TRUE for taxonomic information. Default is FALSE.
#' @param sample FALSE/TRUE for sample metadata. Default is FALSE.
#' @param id sample or study identifier. Default is "abcno".
#' @param vars List of variables of interest from the metadata (sample data).
#'
#' @return Abundance data.frame in long format. Each row is a sample feature.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' ab <- abundance_df(GlobalPatterns, tax = TRUE, sample = TRUE, vars = c("SampleType"))
#' head(ab)
abundance_df = function(phy, tax=FALSE, sample=FALSE, id=character(0),
                        vars=c()) {
  . <- NULL
  ab <- otu_df(phy) %>%
    tidyr::pivot_longer(cols = colnames(.),
                        names_to="taxid",
                        values_to="abundance") %>%
    dplyr::mutate(id = rep(phyloseq::sample_names(phy),
                           each=phyloseq::ntaxa(phy)),
                  taxid = rep(phyloseq::taxa_names(phy),
                              each=phyloseq::nsamples(phy))) %>%
    dplyr::select(id, dplyr::everything())

  if (!rlang::is_empty(id)) {
    ab <- ab %>%
      dplyr::mutate(abcno = rep(sample_df(phy)[[id]],
                                each=phyloseq::ntaxa(phy)))
  }
  if (tax) {
    group <- utils::tail(colnames(tax_df(phy)), 1)
    if (tolower(group) != "species") {
      tax <- tax_df(phy) %>%
        .data[1:which(group == colnames(.))[[1]]]
    } else {
      tax <- tax_df(phy)
    }
    ab <- ab %>%
      dplyr::left_join(tax %>% dplyr::mutate(taxid = phyloseq::taxa_names(phy)))
  }
  if (sample) {
    sam <- sample_df(phy) %>%
      dplyr::mutate(id = phyloseq::sample_names(phy))
    if (length(vars)>0) {
     sam <- sam %>%
        dplyr::select_at(.vars=vars) %>%
        dplyr::mutate(id = phyloseq::sample_names(phy))
    }
    ab <- ab %>%
      dplyr::left_join(sam)
  }
  return(ab)
}



#' Taxa data.frame prevalence
#'
#' @param ps A phyloseq object.
#' @param tax FALSE/TRUE for taxonomic information. Default is FALSE.
#'
#' @return Prevalence data.frame in long format. freq = frequency (range 0-1) and prevalence (%). Each row is a sample feature.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' prev <- prevalence_df(GlobalPatterns, tax = TRUE)
#' head(prev)
prevalence_df = function(phy, tax = FALSE) {
  . <- NULL
  prev <- phy %>%
    transform_ps(., transform="pa") %>%
    otu_df(.)

  prevdf = data.frame(tax = rep(phyloseq::taxa_names(phy),
                                each = phyloseq::nsamples(phy)),
                      n = unname(rowSums(prev)),
                      N = ncol(prev)) %>%
    dplyr::mutate(., freq = n/N,
           prevalence = freq*100)

  if (tax) {
    prevdf <- prevdf %>%
      dplyr::left_join(., tax_df(phy) %>%
                         dplyr::mutate(tax = phyloseq::taxa_names(phy)))
  }

  return(prevdf)
}
