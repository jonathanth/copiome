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
#' agglomerated to lower taxonomic resolution, the higher ranks are not shown.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' taxmat <- tax_df(GlobalPatterns)
#' head(taxmat)
tax_df = function(phy) {
  . <- NULL
  taxmat <- data.frame(phy@tax_table@.Data)
  group <- taxmat %>%
    dplyr::select_if(~ !all(is.na(.))) %>% colnames %>% utils::tail(1)
  if (tolower(group) != "species") {
    taxmat <- taxmat[1:which(group == colnames(taxmat))[[1]]]
  }
  taxmat <- taxmat %>%
    dplyr::mutate(tax = phyloseq::taxa_names(phy))
  return(taxmat)
}


#' Return abundance long-formatted data.frame from phyloseq object.
#'
#' @param phy A phyloseq object
#' @param taxa FALSE/TRUE for taxonomic information. Default is FALSE.
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
#' ab <- abundance_df(GlobalPatterns, taxa = TRUE, sample = TRUE, vars = c("SampleType"))
#' head(ab)
abundance_df = function(phy, taxa=FALSE, sample=FALSE, id=character(0),
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
  if (taxa) {
    group <- utils::tail(colnames(tax_df(phy))[colnames(tax_df(phy)) != "tax"], 1)
    if (tolower(group) != "species") {
      taxmat <- tax_df(phy) %>%
        .data[1:which(group == colnames(.))[[1]]]
    } else {
      taxmat <- tax_df(phy)
    }
    ab <- ab %>%
      dplyr::left_join(taxmat %>% dplyr::mutate(tax = phyloseq::taxa_names(phy)))
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
#' @param taxa FALSE/TRUE for taxonomic information. Default is FALSE.
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
  prevdf <- data.frame(tax = phyloseq::taxa_names(phy),
                     n = otu_df(phy %>%
                                  transform_phy(., transform="pa")) %>%
                       colSums(),
                     N = phyloseq::nsamples(phy)) %>%
    dplyr::mutate(freq = n / N, prevalence = freq * 100)

  if (taxa) {
    prevdf <- prevdf %>%
      dplyr::left_join(., tax_df(phy))
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
#' # Log transformation
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
      print("log(x + min/2)")
      xt <- x %>% apply(2, function(x) log(x + min(x[x > 0])/2))
    } else {
      print(paste0("log(x + ", pseudocount, ")"))
      xt <- x %>% apply(2, function(x) log(x + pseudocount))
    }
  }
  else if (transform == "pa") {
    print("presence/absence")
    xt <- x
    xt[xt > 0] <- 1
  }
  phyloseq::otu_table(phy) <- phyloseq::otu_table(xt, taxa_are_rows = TRUE)
  return(phy)
}



#' Filter phyloseq object based on abundance and prevalence.
#'
#' @param phy A phyloseq object
#' @param abundance Minimum abundance (in relative abundance) to retain. Default is NULL.
#' @param prevalence Minimum prevalence (in percent) to retain. Default is NULL.
#' @param compositional FALSE/TRUE if count table is relative (compositional). Default is TRUE. If FALSE, absolute counts are returned but filtering is done in relative abundances.
#' @param tss FALSE/TRUE if relative abundance is re-calculated after filtering taxa. Default is TRUE.
#'
#' @return A filtered phyloseq object.
#' @export
#'
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' GlobalPatterns
#' filter_phy(GlobalPatterns, abundance = 0.01, prevalence = 5)
filter_phy <-  function(phy, abundance=NULL, prevalence=NULL, compositional=TRUE,
                        tss = TRUE) {
  if (!compositional) {
    phyt <- transform_phy(phy, transform = "compositional")
  } else {
    phyt <- phy
  }
  if (!is.null(abundance)) {
    phyt <- phyloseq::filter_taxa(phyt, function(x) mean(x) > abundance, TRUE)
  }
  if (!is.na(prevalence)) {
    phyt <- phyloseq::prune_taxa(phyloseq::taxa_sums(transform_phy(phyt, transform="pa"))
                     >= ceiling(prevalence * phyloseq::nsamples(phy) / 100), phyt)
  }

  if (compositional & tss) {
    phyt <- transform_phy(phyt, transform = "compositional")
    return(phyt)
  }

  if (compositional & !tss) {
    return(phyt)
  }

  if (!compositional) {
    phyt <- phyloseq::prune_taxa(phyloseq::taxa_names(phyt), phy)
    return(phyt)
  }
}



#' Return virome phyloseq object formatted ready for agglomerating bacteria host taxa.
#'
#' @param phy A virome phyloseq object with host information.
#'
#' @return A virome phyloseq object with taxonomic information related to host, i.e. easy to agglomerate.
#' @export
#' @importFrom magrittr %>%
get_virome_host_phy <- function(phy) {
  taxmat <- tax_df(phy) %>%
    dplyr::select_at(.vars=c("hostKingdom", "hostPhylum", "hostClass",
                             "hostOrder", "hostFamily", "hostGenus", "hostSpecies",
                             "species", "OTU"))
  phyloseq::tax_table(phy) <- phyloseq::tax_table(as.matrix(taxmat))
  return(phy)
}



#' Return number of NAs in sample data
#'
#' @param phy A phyloseq object
#'
#' @return returns a named vector of the number of NAs in sample variables
#' @export
#'
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' na_phy(GlobalPatterns)
na_phy <- function(phy) {
  apply(apply(sample_df(phy), 2, is.na), 2, sum)
}



#' PERMANOVA analysis easy set-up for block and non-block designs.
#'
#' @param dist Distance matrix.
#' @param phy Phyloseq object.
#' @param var Variable outcome of interest.
#' @param covars Vector of covariates to adjust for.
#' @param interaction Interaction term.
#' @param formula User-specified formula. If not provided, `dist ~ var (+ covars)`.
#' @param nproc Number of parallel processors available to use.
#' @param seed Random seed. Default is `123`.
#' @param nperm Number of permutations required.
#' @param verbose Printed output. Default is TRUE.
#' @param by by = "terms" will assess significance for each term (sequentially from first to last), by = "margin" will assess the marginal effects of the terms (each marginal term analysed in a model with all other variables), and by = NULL will assess the overall significance of all terms together. Default is "margin".
#' @param logvars Vector of (co)variables to log-transform. Must be included in `vars` or `covars`, too.
#' @param block Name of group (strata) within which to constrain permutations.
#' @param transient Default TRUE. If the variable outcome is a transient phenotype, when `TRUE` the groups evaluated are transient vs. cross-sectional diseased; when `FALSE` the groups evaluated are transient vs. healthy.
#'
#' @return Returns an anova.cca result as data.frame object with a new column for total observations (n) and missing (n_missing).

#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' sample_data(GlobalPatterns)$contaminated <- sample(c(rep(0, length = 8), rep(1, length = 18)))
#' dist <- distance(transform_phy(GlobalPatterns, transform = "log", pseudocount = 1), method = "bray")
#'
#' ## non-block design
#' myadonis(dist, GlobalPatterns, var = "contaminated", covars = c("SampleType"))
#'
#' ## with block design
#' block <- "sampletype"
#' myadonis(dist, GlobalPatterns, var = "contaminated", block = "SampleType")
myadonis <- function(dist,
                     phy,
                     var,
                     covars=c(),
                     interaction=NA,
                     formula = NA,
                     nproc = 1,
                     seed = 123,
                     nperm = 999,
                     verbose = TRUE,
                     by = "margin",
                     logvars = c(),
                     block = NA,
                     transient = TRUE) {
  phenodata <- sample_df(phy) %>%
    dplyr::filter_at(dplyr::all_of(c(var, covars)),
                     dplyr::all_vars(!is.na(.))) %>%
    dplyr::mutate_at(.vars = logvars, .funs = log)
  if (grepl("transient", var)) {
    if (transient) {
      phenodata <- phenodata[phenodata[[stringr::str_replace(var, "transient", "ever")]] == 1 | phenodata[[var]] == 1, ]
    } else {
      phenodata <- phenodata[phenodata[[stringr::str_replace(var, "transient", "ever")]] == 0 | phenodata[[var]] == 1, ]
    }
  }
  phy <- phyloseq::prune_samples(rownames(phenodata), phy)
  dist <-  subset_dist(dist, phy)
  if (!is.na(interaction)) {var <- paste0(var, "*", interaction)}
  if (is.na(formula)) {
    if (length(covars)==0) {
      formula <- stats::as.formula(paste0("dist ~ ", var))
    } else {
      formula <- stats::as.formula(paste0("dist ~ ", var,
                                          " + ",
                                          paste(covars, collapse = " + ")))
    }
  } else {
    formula <- stats::as.formula(formula)
  }
  if (verbose) {
    print(table(phenodata[[var]]))
    print(formula)
  }
  doParallel::registerDoParallel(nproc)
  if (is.na(block)) {
    set.seed(seed)
    res <- vegan::adonis2(formula,
                          data = phenodata,
                          by = by,
                          permutations = nperm,
                          parallel = nproc)
  } else { # block design, allows to constrain permutations within each block/strata
    if (verbose) {print(paste0("Correcting for strata using block design for ", block))}
    perm <- permute::how(nperm = nperm)
    set.seed(seed)
    assign(block, block)
    permute::setBlocks(perm) <- with(phenodata, get(block))
    res <- vegan::adonis2(formula,
                          data = phenodata,
                          by = by,
                          permutations = perm,
                          parallel = nproc)
  }
  if (verbose) {
    print(res)
  }
  res_df <- data.frame(res) %>%
    dplyr::mutate(variable = rownames(.)) %>%
    dplyr::rename(pval = colnames(.)[5]) %>%
    dplyr::mutate(n = nrow(phenodata),
                  n_missing = phyloseq::nsamples(phy) - nrow(phenodata))
  return(res_df)
}

