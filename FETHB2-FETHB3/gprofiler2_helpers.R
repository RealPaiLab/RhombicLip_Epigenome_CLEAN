library(dplyr)

#' Run g:GOSt gene set enrichment analysis from g:Profiler and save results.
#'
#' @description For more information on the input parameters, see
#'   `gprofiler2::gost()`.
#'
#' @param query Gene list (as vector) or multiple gene lists (as list of
#'   vectors).
#' @param organism Unique ID for custom gene sets from a GMT file that was
#'   uploaded using `gprofiler2::upload_GMT_file()`. Defaults to
#'   `"gp__CEZi_x2oO_5Lg"` which is a [custom gene set from Bader
#'   lab](https://download.baderlab.org/EM_Genesets/) (GOBP, all pathways, no
#'   IEA, 2023-08-08). See `<species_dge software directory>/README.md` for more
#'   details.
#' @param significant If `TRUE` (default), returns only significant results.
#' @param evcodes If `TRUE` (default), includes evidence codes in the results.
#' @param correction_method Defaults to `"fdr"`, see `gprofiler2` documentation
#'   for additional options.
#' @param custom_bg Optional gene list to be used as background for statistical
#'   testing.
#' @param filename File name to write results to (without extension). E.g., if
#'   `filename = "results"`, then output files will be `"results.csv"`,
#'   `"results.rds"`, and `"results.html"`.
#'
#' @return Output of `gprofiler2::gost()` (named list). If `filename` is not
#'   `NULL`, then the full output will be saved to `*.rds`, the table of
#'   enriched pathways will be saved to `*.csv`, and the Manhattan plot will be
#'   saved to `*.html`.
#' 
run_gost <- function(
  query,
  organism = "gp__CEZi_x2oO_5Lg",
  significant = TRUE,
  evcodes = TRUE,
  correction_method = "fdr",
  custom_bg = NULL,
  filename = NULL,
  ...
) {
  gost_res <- gost(
    query = query,
    organism = organism,
    significant = significant,
    evcodes = evcodes,
    correction_method = correction_method,
    custom_bg = custom_bg,
    ...
  )

  if (is.null(gost_res)) {
    warning("***WARNING: `gost` did not return any results***")
    return(gost_res)
  }
  
  # save GOSt results
  if (!is.null(filename)) {
    message(sprintf("Saving %s.rds", filename))
    readr::write_rds(
      x = gost_res,
      file = paste0(filename, ".rds")
    )
    
    message(sprintf("Saving %s.csv", filename))
    readr::write_csv(
      x = gost_res$result,
      file = paste0(filename, ".csv")
    )
    
    # save plot
    message(sprintf("Saving %s.html", filename))
    gostplot(gost_res, pal = "#1e3765") %>% 
      htmlwidgets::saveWidget(
        widget = .,
        file = paste0(filename, ".html"),
        title = basename(filename)
      )
  }
  
  return(gost_res)
}


#' Make Generic Enrichment Map from g:GOSt results.
#'
#' @description Takes the output from runnng `gost` and converts the results
#'   into a dataframe that conforms to the Generic Enrichment Map format for
#'   downstream analysis with EnrichmentMap.
#'
#' @param gost_res Output from running `gost`, i.e., a named list of length 2
#'   (`result` and `meta`).
#' @param min_max_terms Set minimum and maximum (inclusive) term size of the
#'   pathways. For example, to filter out pathways with < 10 genes or > 250
#'   genes, set this argument to `c(10, 250)`. Default is `c(1, Inf)` which
#'   doesn't filter out the pathways. Note that the q-values are not adjusted
#'   for the filtering and remain the same as the unfiltered q-values.
#' @param phenotype Either a number (`"+1"` or `"-1"`) or a list of length 2
#'   named `down` and `up` with the elements being the query names. For example,
#'   `phenotype = list(down = "query_down", up = "query_up")` where
#'   `"query_down"` and `"query_up"` are in `gost_res$result$query`. `down` will
#'   be mapped to `"-1"` while `up` will be mapped to `"+1"`. If `NULL`, the
#'   phenotype column will not be added.
#' @param keep_query If TRUE, the dataframe includes the `query` column from the
#'   gost results. Useful if the dataframe will be saved to multiple files
#'   (e.g., one for each query). The `query` column must be removed before
#'   before saving to a file.
#'
#' @returns A dataframe in the Generic Enrichment Map format.
#' 
gost_res2gem <- function(
  gost_res,
  min_max_terms = c(1, Inf),
  phenotype = NULL,
  keep_query = FALSE
) {
  if (min_max_terms[1] != 1 | !is.infinite(min_max_terms[2])) {
    gost_res$result <- gost_res$result %>% 
      dplyr::filter(
        term_size >= min_max_terms[1] & term_size <= min_max_terms[2]
      )
  }
  
  gem <- gost_res$result %>% 
    dplyr::select(query, term_id, term_name, p_value, intersection) %>% 
    dplyr::rename(
      GO.ID = term_id,
      Description = term_name,
      p.Val = p_value,
      Genes = intersection
    ) %>% 
    mutate(FDR = p.Val, .after = p.Val)
  
  if (length(phenotype) == 1) {
    gem <- mutate(
      gem,
      Phenotype = phenotype,
      .before = Genes
    )
  } else if (length(phenotype) == 2) {
    gem <- mutate(
      gem,
      Phenotype = case_when(
        query %in% phenotype$down ~ "-1",
        query %in% phenotype$up ~ "+1",
        TRUE ~ NA_character_
      ),
      .before = Genes
    )
  }
  
  if (!keep_query) {
    gem <- dplyr::select(gem, !query)
  }
  
  return(gem)
}


#' Wrapper to save a Generic Enrichment Map file.
#'
#' @description Saves a dataframe as a tab-delimited text file for downstream
#'   analysis with EnrichmentMap.
#'
#' @param gem_df A dataframe conforming to the Generic Enrichment Map format.
#'   Generated using `gost_res2gem`.
#' @param file File to write to. Extension should be `.txt`.
#'
#' @returns None.
#'
write_gem <- function(
  gem_df,
  file
) {
  write_tsv(x = gem_df, file = file)
}
