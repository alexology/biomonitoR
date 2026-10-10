#' Compute a Biological Trait Index
#'
#' This function computes a biological trait-based index (e.g., trophic index)
#' based on trait composition and fuzzy coding. It supports both predefined and custom trait groupings.
#'
#' @param x Results of `aggregate_taxa()`.
#' @param tax_lev Character. The taxonomic level to use for trait assignment. Default is `"Taxa"`.
#' @param trait_db Optional. A custom trait database to use. If `NULL`, internal trait values will be assigned.
#' @param type Character. Data type: `"ab"` (abundance), `"po"` (presence only). Default is `"ab"`.
#' @param trans Optional The function used to transform the abundances, by default log1p.
#' @param index Character vector. Index/Indices to calculate. The option are: trophic_t, saprobity_t, temperature_t, and custom. Default trophic_t.
#' @param rm_taxa Logical. If TRUE, taxa with all trait values equal to 0 are removed. If FALSE, taxa with a 0 in one or more traits are retained. This applies only when more than one trait is used to calculate the index.
#' @param group Character vector. Set the biotic group of interest. Possible values are "mi" for macroinvertebrates, "mf" for macrophytes and "fi" for fish. This option will not be considered if a custom reference database is provided. Default to mi.
#' @param num Character vector. Trait(s) to be used in the numerator for a custom index. Required if `index = "custom"`.
#' @param den Character vector. Trait(s) to be used in the denominator for a custom index. Required if `index = "custom"`.
#' @param col_blocks Integer vector. A vector containing the number of modalities for each trait. IF NULL, the default Tachet et al. (2010) col_blocks will be assigned.
#'
#' @return A numeric vector or data frame with the trait-based index per sample.
#' 
#' @details
#'
#' - When `index = "custom"`, you must provide both `num` and `den`, defining the trait columns for the ratio.
#'
#' The function automatically filters taxa with zero values in the selected trait columns and averages traits across blocks.
#'
#' @examples
#' \dontrun{
#' 
#' Example 1
#'  data(mi_prin)
#'  mi_prin <- mi_prin[1:10, 1:5]
#'  mi_prin
#'  data_bio <- as_biomonitor(mi_prin)
#'  data_agr <- aggregate_taxa(data_bio)
#'  bio_index(x = data_agr)
#'  
#'  Example 2
#'  set.seed(16)
#'  t_db <- data.frame(Taxa = data_ts$Taxa,
#'                     a_1 = runif(length(data_ts$Taxa)),
#'                     b_1 = runif(length(data_ts$Taxa)),
#'                     c_1 = runif(length(data_ts$Taxa)),
#'                     a_2 = runif(length(data_ts$Taxa)),
#'                     b_2 = runif(length(data_ts$Taxa)),
#'                     c_2 = runif(length(data_ts$Taxa)),
#'                     d_2 = runif(length(data_ts$Taxa)))
#'  bio_index(x = data_agr, num = "a_1", den = c("b_1", "c_1"), 
#'            index = "custom", col_blocks = 3, trait_db = t_db)
#' }
#'
#' @seealso \code{\link{assign_traits}}, \code{\link{average_traits}}, \code{\link{fuzzy_trait_ratio}}
#' @export

# Set blocks for traits assignation
select_blocks <- function(df, blocks) {
  starts <- cumsum(c(1, blocks[-length(blocks)]))
  ends <- cumsum(blocks)
  
  blocks_list <- mapply(function(start, end) df[, start:end, drop = FALSE],
                        starts, ends, SIMPLIFY = FALSE)
  
  return(blocks_list)
}


# Main function
bio_index <- function(x, tax_lev = "Taxa", trait_db = NULL, type = "ab", trans = NULL, 
                      index = c("trophic_t"), rm_taxa = TRUE, group = "mi", 
                      num = NULL, den = NULL, col_blocks = NULL, col_order = NULL){
  
  index_values <- c("trophic_t", "saprobity_t", "temperature_t", "custom")
  
  # Add a check for index
  if (!all(index %in% index_values)) {
    stop("Invalid 'index' value(s). Must be one or more of: ", paste(index_values, collapse = ", "))
  }
  
  if (any(index != "custom") && (!is.null(num) || !is.null(den))) {
    stop("You can set numerator and denominator ONLY if all indices are 'custom'")
  }
  if ("custom" %in% index && (is.null(num) || is.null(den))) {
    stop("For custom index, both numerator and denominator must be provided")
  }
  
  compute_single_index <- function(idx) {
    
    # Select the numerator/denominator for the index calculation
    if (idx == "trophic_t") {
      numerator <- c("TROPHIC_1")
      denominator <- c("TROPHIC_1", "TROPHIC_2", "TROPHIC_3")
      col_blocks <- 3
      col_order <- c("TROPHIC_1", "TROPHIC_2", "TROPHIC_3")
    }
    
    if (idx == "saprobity_t") {
      numerator <- c("SAPROBITY_1", "SAPROBITY_2")
      denominator <- c("SAPROBITY_1", "SAPROBITY_2","SAPROBITY_3", "SAPROBITY_4", "SAPROBITY_5")
      col_blocks <- 5
      col_order <- c("SAPROBITY_1", "SAPROBITY_2","SAPROBITY_3", "SAPROBITY_4", "SAPROBITY_5")
    }
    
    if (idx == "temperature_t") {
      numerator <- c("TEMPERATURE_1")
      denominator <- c("TEMPERATURE_1", "TEMPERATURE_2", "TEMPERATURE_3")
      col_blocks <- 3
      col_order <- c("TEMPERATURE_1", "TEMPERATURE_2", "TEMPERATURE_3")
    }
    
    if (idx == "custom") {
      numerator <- num
      denominator <- den
    }
    
    # Assign traits
    data_ts <- if (is.null(trait_db)) {
      assign_traits(x, tax_lev = tax_lev, group = group)
    } else {
      assign_traits(x, tax_lev = tax_lev, trait_db = trait_db)
    }
    
    data_ts_temp <- data_ts[, c("Taxa",
                                "Taxa_db",
                                "Traits_taxlev",
                                "Traits_real",
                                "Taxonomic_distance",
                                col_order)]
    
    db <- data_ts_temp[, col_order]
    db_block <- select_blocks(db, col_blocks)
    db_block <- lapply(db_block, function(block) rowSums(block) > 0)
    
    if (rm_taxa) {
      non_zero_rows <- Reduce("&", db_block)
    } else {
      non_zero_rows <- Reduce("|", db_block)
    }
    
    data_ts <- data_ts_temp[non_zero_rows, , drop = FALSE]
    data_ts_av <- average_traits(data_ts, col_blocks = col_blocks)
    data_ts_av[is.na(data_ts_av)] <- 0
    
    res <- fuzzy_trait_ratio(x, tax_lev = tax_lev, trait_db = data_ts_av, type = type, 
                             numerator = numerator, denominator = denominator, trans = trans)
    
    return(res)
  }
  
  
  results <- lapply(index, compute_single_index)
  names(results) <- index
  
  return(results)
  
}
