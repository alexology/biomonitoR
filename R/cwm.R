#' Community-Weighted Mean values of traits
#'
#' @description
#' This function calculates the community-weighted means of trait categories.
#'
#' @details
#' This function first takes the abundance table corresponding to the desired
#' taxonomic level from the `x` aggregatoR object.
#'
#' Then it searches from the trait data base all the information available at
#' the desired level and, if required, calculates the corresponding averaged
#' trait values (e.g. the family trait values are obtained by averaging all the
#' trait values from taxa with trait information within this family).
#'
#' Finally, the community mean trait values are calculated using the transformed
#' abundances (using the `trans` function) as weigths
#'
#' @param x Results of the function `aggregate_taxa()`.
#' @param trait_db A trait data base with a column `Taxa` and the other columns
#'   containing the traits.
#'   By default, the data base used is the one from Tachet *et al* (2010) that
#'   can be retrieved from
#'   [freshwaterecology.info](https://www.freshwaterecology.info/) website
#'   (Schmidt-Kloiber & Hering, 2015).
#' @param tax_lev Character string giving the taxonomic level used to retrieve
#'   trait information. Possible levels are `"Taxa"`, `"Species"`, `"Genus"`,
#'   `"Family"` as returned by the [aggregatoR] function.
#' @param trans The function used to transform the abundances, by default
#'   `log1p`.
#' @param traceB When set to TRUE returns a list with the results of the cwm function
#'   and the traits value used for the calculation.
#'
#' @return a table with the CWM values of each trait (trait modality)
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#'
#' @examples
#' data(macro_ex)
#'
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' data_ts <- assign_traits( data_agr )
#'
#' # averaging
#' data_ts_av <- average_traits( data_ts )
#'
#' # community specialization index
#' cwm(x = data_agr, trait_db = data_ts_av, tax_lev = "Taxa", trans = log1p)
#' cwm(x = data_agr, trait_db = data_ts_av, tax_lev = "Taxa",
#'     trans = function(x) {
#'         ifelse(x > 0, 1, 0)
#'     })
#' cwm(x = data_agr, trait_db = data_ts_av, tax_lev = "Genus", trans = log1p)
#'
#' @seealso [aggregate_taxa]
#'
#' @references Tachet, H., Richoux, P., Bournaud, M., & Usseglio-Polatera, P.
#'   (2010). Invertebres d'eau douce: systematique, biologie, ecologie. Paris:
#'   CNRS editions.
#' @references Schmidt-Kloiber, A., & Hering, D. (2015). www. freshwaterecology.
#'   info-An online tool that unifies, standardises and codifies more than
#'   20,000 European freshwater organisms and their ecological preferences.
#'   Ecological indicators, 53, 271-282.
#'
#' @export

cwm <- function( x , trait_db = NULL, tax_lev = "Taxa" , trans = log1p , traceB = FALSE ) {

  classCheck( x )

  if( is.null( trait_db )){
    trait_db = traitsTachet
  } else {
    trait_db = trait_db
  }


  if (! tax_lev %in% c("Family", "Genus", "Species", "Taxa")) {
    return("tax_lev should be one of the following: Family, Genus, Species or Taxa")
  }

  # create dummy variables to avoid R CMD check NOTES

  traitsTachet <- Taxa <- modality <- affinity <- Phylum <- Subspecies <-
    Abundance <- Sample <- Weight <- Affinity <- totWeight <-
    weightedAffinity <- Category <- . <- NULL

  abundances <- x[[tax_lev]]
  colnames(abundances)[1] <- "Taxa"

  if( inherits( x , "bin" ) ){
    abundances <- to_bin( abundances )
  }

  # remove unassigned taxa from abundances
  if("unassigned" %in% abundances[ , "Taxa"]){
    z <- which(abundances[ , "Taxa" ] == "unassigned")
    abundances <- abundances[ -z ,] # remove unassigned row from the species count
  }

  abundances <- merge( abundances , trait_db[ , "Taxa" , drop = FALSE ] )
  trait_db <- merge( trait_db , abundances[ , "Taxa" , drop = FALSE ] )

  abundances$Taxa <- as.character( abundances$Taxa )
  trait_db <- trait_db %>% mutate( Taxa = as.character( Taxa ) )
  res <- abundances                                       %>%
    gather(key = Sample, value = Abundance, -Taxa) %>%
    mutate(Sample = factor(Sample,
                           levels = colnames(abundances)[-1]),
           Taxa   = as.character(Taxa),
           Weight = trans(Abundance))              %>%
    left_join(group_by(., Sample)                  %>%
                summarise(totWeight = sum(Weight)),
              by = "Sample")                       %>%
    left_join(gather(trait_db,
                     key   = Category,
                     value = Affinity, -Taxa),
              by = "Taxa")                         %>%
    mutate(weightedAffinity = (Weight * Affinity) /
             totWeight)                            %>%
    group_by(Sample, Category)                     %>%
    summarise(Affinity = sum(weightedAffinity,
                             na.rm = TRUE))        %>%
    spread(key = Category, value = Affinity)       %>%
    as.data.frame()


  if( traceB == FALSE ){
    res
  } else {

    df1 <- trait_db %>% semi_join( abundances , "Taxa" )
    df2 <- abundances %>% semi_join( trait_db , "Taxa" )
    df3 <- abundances[ ! abundances$Taxa %in% df2$Taxa , "Taxa" ]
    if( length( df3) == 0) ( df3 <- NA )
    res.list <- list( res , df1 , df2 , df3 )
    names( res.list ) <- c( "results" , "traits" , "taxa" , "taxa_excluded" )
    res.list

  }
}
