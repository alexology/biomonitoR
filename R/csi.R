#' Community trait specialization
#'
#' This function calculates the community trait specialization.
#' @param x results of function aggregatoR
#' @param traitDB a trait data base with a column `Taxa` and the other columns
#'   containing the traits.Please check [traitScaling] for building the trait database.
#' @param taxLev character string giving the taxonomic level used to retrieve
#'   trait information. Possible levels are `"Taxa"`, `"Species"`, `"Genus"`,
#'   `"Family"` as returned by the [aggregatoR] function.
#' @param trans the function used to transform the abundances, by default
#'   `log1p`.
#' @param traceB when set to TRUE returns a list with the results of the cwm function
#'   and the traits value used for the calculation.
#'
#' @details This function first takes the abundance table corresponding to the desired
#' taxonomic level from the `x` aggregatoR object.
#'
#' For each taxon and each trait, a taxon specialization index is calculated
#' using the following formula:
#'
#' \deqn{TSI = (sum(c^{2}_tik) - 1/K) / (1 - 1/K)}
#'
#' with \eqn{c^{2}_tik} being the affinity of taxon `t` for the modality `k` (among
#' `K`) of trait `i`.
#'
#' Finally, the community trait specialization is calculated for each trait by
#' averaging the TSIs using the transformed abundances (using the `trans`
#' function) as weigths. Only the taxa having traits contributes to the calculation of
#' index.
#'
#' @note USE WITH CAUTION, STILL IN DEVELOPMENT.
#'
#' @return a table with the CSI values for each trait
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup n_distinct
#' @importFrom tidyr gather spread
#'
#' @examples
#' data(macro_ex)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.ts <- traitScaling( data.agR )
#'
#' # averaging
#' data.ts.av <- traitsMean( data.ts )
#'
#' # taxon specialization index
#'
#' tsi(x = data.agR, traitDB = data.ts.av, taxLev = "Taxa")
#'
#' # community specialization index
#' csi(x = data.agR, traitDB = data.ts.av, taxLev = "Taxa", trans = log1p)
#' csi(x = data.agR, traitDB = data.ts.av, taxLev = "Taxa",
#'     trans = function(x) {
#'         ifelse(x > 0, 1, 0)
#'     })
#' csi(x = data.agR, traitDB = data.ts.av, taxLev = "Genus", trans = log1p)
#'
#' @seealso [aggregatoR]
#'
#' @references Mondy, C. P., & Usseglio-Polatera P. (2013) Using Fuzzy-Coded
#'   Traits to Elucidate the Non-Random Role of Anthropogenic Stress in the
#'   Functional Homogenisation of Invertebrate Assemblages. Freshwater Biology,
#'   59 (3), 584-600. https://doi.org/10.1111/fwb.12289.
#' @references Tachet, H., Richoux, P., Bournaud, M., & Usseglio-Polatera, P.
#'   (2010). Invertebres d'eau douce: systematique, biologie, ecologie. Paris:
#'   CNRS editions.
#' @references Schmidt-Kloiber, A., & Hering, D. (2015). www. freshwaterecology.
#'   info-An online tool that unifies, standardises and codifies more than
#'   20,000 European freshwater organisms and their ecological preferences.
#'   Ecological indicators, 53, 271-282.
#'
#' @export
#' @export tsi



csi <- function( x, traitDB = NULL, taxLev = "Taxa", trans = log1p, traceB = FALSE ) {

  classCheck( x )

  if( is.null( traitDB )){
    trait_db <- traitsTachet
  } else {
    trait_db = traitDB
  }

  # create dummy variables to avoid R CMD check NOTES

  Taxa <- Trait <- Modality <- modality <- affinity <- Affinity <- Phylum <- Subspecies <-
    Abundance <- Sample <- Weight <- totWeight <- k <-
    TSI <- CSI <- weightedTSI <-  . <- NULL

  abundances <- x[[ taxLev ]]
  colnames(abundances)[1] <- "Taxa"

  if( inherits( x , "bin" ) ){
    abundances <- to_bin( abundances )
  }

  # remove unassigned taxa from abundances
  if("unassigned" %in% abundances[ , "Taxa" ] ){
    z <- which( abundances[ , "Taxa" ] == "unassigned")
    abundances <- abundances[ -z ,] # remove unassigned row from the species count
  }

  abundances$Taxa <- as.character( abundances$Taxa )
  trait_db <- trait_db %>% mutate( Taxa = as.character( Taxa ) )

  abundances <- merge( abundances , trait_db[ , "Taxa" , drop = FALSE ] )
  trait_db <- merge( trait_db , abundances[ , "Taxa" , drop = FALSE ] )

  tsi <- trait_db                                       %>%
    semi_join( abundances , "Taxa" )                    %>%
    gather( key = Modality , value = Affinity , -Taxa ) %>%
    mutate(Trait = strsplit( Modality, split = "_" )    %>%
             sapply(FUN = '[[', 1))                     %>%
    group_by(Taxa, Trait)                               %>%
    mutate(k = n_distinct(Modality))                    %>%
    summarise(TSI = (sum(Affinity^2) - 1 / unique(k)) /
                (1 - 1 / unique(k)))

  res <- abundances                                 %>%
    semi_join( trait_db , "Taxa" )                  %>%
    gather( key = Sample , value = Abundance , -Taxa )  %>%
    mutate(Sample = factor(Sample,
                           levels = colnames(abundances)[ -1 ] ),
           Taxa   = as.character( Taxa ) ,
           Weight = trans( Abundance ) )               %>%
    left_join( group_by(., Sample)                   %>%
                summarise( totWeight = sum(Weight) ),
              by = "Sample" )                       %>%
    left_join( tsi , by = "Taxa")                   %>%
    mutate(weightedTSI = (Weight * TSI) /
             totWeight)                             %>%
    group_by( Sample, Trait )                         %>%
    summarise( CSI = sum( weightedTSI, na.rm = TRUE ) ) %>%
    spread( key = Trait, value = CSI )                %>%
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
