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
  otu <- as.data.frame(methods::as(otu, "matrix"))
  return(otu)
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
    dplyr::select_if(~ !all(is.na(.))) %>% colnames %>% utils::tail(1) %>%
    dplyr::mutate(tax = phyloseq::taxa_names(phy))
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
                        names_to="tax",
                        values_to="abundance") %>%
    dplyr::mutate(id = rep(phyloseq::sample_names(phy),
                           each=phyloseq::ntaxa(phy)),
                  tax = rep(phyloseq::taxa_names(phy),
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
      dplyr::left_join(tax %>% dplyr::mutate(tax = phyloseq::taxa_names(phy)))
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
#' @param phy A phyloseq object.
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
prevalence_df <-  function(phy, taxa = FALSE) {
  . <- NULL
  prev <- phy %>%
    transform_phy(., transform="pa") %>%
    otu_df(.)

  prevdf = data.frame(tax = rep(phyloseq::taxa_names(phy),
                                each = phyloseq::nsamples(phy)),
                      n = unname(rowSums(prev)),
                      N = ncol(prev)) %>%
    dplyr::mutate(., freq = n/N,
           prevalence = freq*100)

  if (taxa) {
    prevdf <- prevdf %>%
      dplyr::left_join(., tax_df(phy) %>%
                         dplyr::mutate(tax = phyloseq::taxa_names(phy)))
  }

  return(prevdf)
}



#' Common transformations for phyloseq count data.
#'
#' @param phy A phyloseq object.
#' @param transform Transformation to be applied: "compositional" (i.e. relative abundance), "clr" (centered log ratio), "log" (log10) and "pa" (presence/absence).
#' By default "clr" applies a pseudocount of min(relative abundance)/2 and "log" applies a pseudocount of min/2 before taking logs. The user can also define a certain value
#' for the pseudocount.
#' @param pseudocount Pseudocount to be used in "clr" or "log" transformations. Default is NULL
#'
#' @return Transformed phyloseq object
#' @export
#'
#' @examples
#' library(phyloseq)
#' data(esophagus)
#'
#' # No transformation
#' head(otu_df(esophagus))
#'
#' # Relative abundance
#' head(otu_df(transform_phy(esophagus, transform = "compositional")))
#'
#' # Centered log-ratio (clr) transformation
#' head(otu_df(transform_phy(esophagus, transform = "clr")))
#' head(otu_df(transform_phy(esophagus, transform = "clr", pseudocount = 0.1)))
#'
#' # Log10 transformation
#' head(otu_df(transform_phy(esophagus, transform = "log")))
#' head(otu_df(transform_phy(esophagus, transform = "log", pseudocount = 1)))
#'
#' # Presence/absence transformation
#' head(otu_df(transform_phy(esophagus, transform = "log")))
transform_phy <-  function(phy, transform="compositional",
                          pseudocount=NULL) {
  x <- otu_df(phy) %>% t()
  if (transform == "compositional") {
    print("x/max(sum(x))")
      xt <- apply(x, 2, function(x) {
        x/max(sum(x), 1e-32)})
  }
  else if (transform == "clr") {
    if (any(x < 0)) {
      stop("Non-negative counts needed for clr transformation.")
    }
    xt <- apply(x, 2, function(x) {
      x/max(sum(x), 1e-32)})
    colnames(xt) <- colnames(x)
    if (any(xt == 0)) {
      v <- as.vector(xt)
      if (is.null(pseudocount)) {
        print("log(x + min/2) - mean(log(x + min/2))")
        minval <- min(v[v > 0])/2
        xt <- xt + minval
      } else {
        print(paste0("log(x + ", pseudocount, ") - mean(log(x + ", pseudocount, "))"))
        xt <- xt + pseudocount
      }
    }
    d <- t(apply(xt, 2, function(x) {
      log(x) - mean(log(x))
    }))
    if (nrow(d) == ncol(xt)) {
      rownames(d) <- colnames(xt)
      colnames(d) <- rownames(xt)
    }
    else {
      colnames(d) <- colnames(xt)
      rownames(d) <- rownames(xt)
    }
    xt <- t(d)
  }
  else if (transform == "log") {
    if (is.null(pseudocount)) {
      print("log10(x + min/2)")
      xt <- x %>% apply(2, function(x) log(x + min(x[x > 0])/2))
    } else {
      print(paste0("log10(x + ", pseudocount, ")"))
      xt <- x %>% apply(2, function(x) log(x + pseudocount))
    }
  }
  else if (transform == "pa") {
    print("presence/absence")
    xt <- x
    xt[xt > 0] <- 1
  }
  otu_table(phy) <- phyloseq::otu_table(xt, taxa_are_rows = TRUE)
  return(phy)
}
