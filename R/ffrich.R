#' Fuzzy coded functional richness
#'
#' This function calculates the functional richness based on trait category.
#'
#' This function first takes the abundance table corresponding to the desired
#' taxonomic level from the `x` aggregatoR object.
#'
#' Functional richness (FRic) represents the amount of functional space filled by
#' the community (Villeger et al., 2008) and it is related to the community use of
#' resources and productivity (Mason et al., 2005). FRic is defined by the trait
#' extremes and thus reflects the potential maximum functional dissimilarity.
#' FRic is calculated as the hypervolume enclosing the functional space filled
#' by the community. For this the convex hull approach is applied (Cornwell et al.,
#' 2006) using the Quickhull algorithm (Barber et al., 1996) that estimates the minimum
#'  convex hull which includes all the species considered in the previously defined
#'  functional space. Basically, this algorithm determines the taxa in the most extreme
#'  points of the functional space, links them to build the convex hull in order to
#'  calculate the volume inside it. In particular, the convex hull of a set of points
#'  S in n dimensions is the intersection of all convex sets containing S.
#'  For N points , ..., , the convex hull C is then given by the expression:
#'
#'
#' The functional T dimensional space is built using a certain number of dimensions
#' (T) determined by the axes of a principal component analysis based on the trait
#'  dissimilarity matrix. Using a poor quality functional space could led to
#'  a biased assessment of FRic and false ecological conclusions so the number of
#'  axes retained for its estimation is case-specific and decided following the
#'  method proposed in Maire et al., (2015): a pragmatic approach consisting of
#'  computing all the possible functional spaces and selecting the most
#'  parsimonious one. This method uses the mSD index (mean squared deviation
#'  between the initial functional distance and the scaled distance in the
#'  functional space), which accounts explicitly for the deviation between
#'  the initial and final distance and penalizes the strong deviation. The mSD
#'  index has been widely used in statistics to assess errors and it has been
#'  demonstrated to work in different contexts and situations (Maire et al., 2015).
#'  In addition, when using Gower's distance the mSD ranges from 0 and 1,
#'  which helps to interpret quality. Finally, the resulting FRic variable
#'  is standardized by its maximum, ranging from 0 to 1. In addition, it must
#'  be considered that the number of taxa must be higher than the number of traits
#'  to have reliable FRic values (Villeger et al., 2008).
#'
#' @param x results of function aggregatoR
#' @param traitDB a trait data base with a column `Taxa` and the other columns
#'   containing the traits. If a trait has several modalities they should be
#'   named as follow: TRAIT_MODALITY.
#'
#'   By default, the data base used is the one from Tachet *et al* (2010) that
#'   can be retrieved from
#'   [freshwaterecology.info](https://www.freshwaterecology.info/) website
#'   (Schmidt-Kloiber & Hering, 2015).
#' @param agg should ffrich aggregate user's traitDB of higher taxonomic level? TRUE to aggregate, otherwise FALSE.
#'    For instance, if user's traitDB has both Halesus and Limnephilidae, ffrich will aggregate traits value if ADD = TRUE.
#' @param traitSel interactively select traits.
#' @param taxLev character string giving the taxonomic level used to retrieve
#' trait information. Possible levels are `"Taxa"`, `"Species"`, `"Genus"`,
#' `"Family"` as returned by the [aggregatoR] function.
#' @param traceB if TRUE ffrich will return a list with 2 elements, the first being the ffrich values and the second the database used for the calculation. Useful to check missing taxa.
#' @param nbdim maximum number of dimensions for multidimensional functional spaces. By default, nbdim=7                   ##
#'      		Final number of dimensions depends on the number of positive eigenvalues (after correction) obtained with PCoA
#' @param metric to be used to compute functional distance, "Euclidean" or "Gower" (=default)
#' @param corr_method TRUE will aplly the Cailliez correction.
#'
#' @note USE WITH CAUTION, STILL IN DEVELOPMENT.
#'
#' @return a vector with fuzzy functional richness results.
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#' @importFrom stats complete.cases
#'
#' @examples
#' data(macro_ex)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#'
#' cwm(x = data.agR, taxLev = "Taxa")
#'
#' @seealso [aggregatoR]
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

ffrich <- function(x, traitDB = NULL, agg = FALSE, traitSel = FALSE, taxLev = "Family", traceB = FALSE, nbdim = 7, metric = "Gower", corr_method = FALSE){

  # check if user provided a trait database, otherwise use traitsTachet
  # if traitsTachet has to be used check for class biomonitoR and "mi"
  if( is.null( traitDB )){
    # check if the object d is of class "biomonitoR" & "mi"
    classCheck(x, group = "mi")
    traitDB = traitsTachet
    # useful for the if condition later
    mes <- TRUE
  } else {
    # check if the object d is of class "biomonitoR"
    classCheck(x)
    traitDB = traitDB
    # trim and capitalise the first letter in the DB provided by the user
    traitDB[ , "Taxa"] <- as.factor( sapply( trim( traitDB[ , "Taxa"] ), capWords, USE.NAMES = F ) )
    mes <- FALSE
  }

  if( traitSel == TRUE){
    rma <- select.list( names( traitDB[ -which( names( traitDB ) %in% "Taxa")] ) , title = "Traits selection"  , graphics = TRUE , multiple = T )
    traitDB <- traitDB %>%
      select( c("Taxa", rma) )
    # trim and capitalise the column Taxa of the user' trait database
    traitDB$Taxa <- apply( as.data.frame( trim( traitDB$Taxa ) ) , 1 , capWords)
  }


  abundances <- x[[taxLev]]
  colnames(abundances)[1] <- "Taxa"
  taxa <- as.character(abundances$Taxa)


  if( isTRUE(mes) | agg == TRUE ){
    # create dummy variables to avoid R CMD check NOTES
    traitsTachet <- Taxa <- modality <- affinity <- Phylum <- Subspecies <-
      Abundance <- Sample <- Weight <- Affinity <- totWeight <-
      weightedAffinity <- Category <- . <- NULL

    if (! taxLev %in% c("Family", "Genus", "Species", "Taxa")) {
      return("taxLev should be one of the following: Family, Genus, Species or Taxa")
    }

    # I included also Ad., Gen. and Lv. in gsub pattern. Currently biomonitoR does not handle adults and larvae at the same time
    trait_db <- traitDB                               %>%
      (function(df) {
        mutate(df,
               Taxa = gsub(pattern     = "sp.|Ad.|Lv.|Gen." ,
                           replacement = "",
                           x           = Taxa))
      })                                              %>%
      gather(key = modality, value = affinity, -Taxa) %>%
      group_by(Taxa, modality)                        %>%
      summarise(affinity = mean(affinity))            %>%
      spread(key = modality, value = affinity)        %>%
      ungroup()
      trait_db$Taxa <- trimws(trait_db$Taxa)


    if (length(taxa[taxa != "unassigned"]) == 0) {
      return("At least one taxa should be identified at a level compatible with the indicated taxLev")
    }

    if (taxLev == "Taxa") {
      level <- sapply(select(x$Tree, Phylum:Subspecies),
                      function(i) {
                        as.character(i) == as.character(x$Tree$Taxa)
                      }) %>%
        (function(df) {
          colnames(df)[apply(df, MARGIN = 1, which)]
        })
    } else {
      level <- rep(taxLev, length(taxa))
    }

    ref <- select(mi_ref, Phylum:Taxa)

    taxa_traits <- mutate(ref, Taxa = as.character(Taxa)) %>%
      left_join(mutate(trait_db, Taxa = as.character(Taxa)),
                by = "Taxa")                              %>%
      (function(df) {
        lapply(seq(length(taxa)),
               function(i) {
                 df[df[, level[i]] == taxa[i],] %>%
                   select(-(Phylum:Taxa))       %>%
                   colMeans(na.rm = TRUE)
               })               %>%
          do.call(what = rbind) %>%
          data.frame(Taxa = taxa, ., stringsAsFactors = FALSE)
      })

  }

  if( mes == FALSE & agg == FALSE ){
    taxa_traits <- traitDB
  }

  taxa_traits <- as.data.frame(taxa_traits)
  # be sure that taxa_traits contains only the Taxa present in the user's community data
  taxa_traits <- taxa_traits[ taxa_traits[, "Taxa"]  %in% abundances[, "Taxa"] , ]
  # remove the column called Taxa. I prefer this way and not traitDB[ , -1, drop = F] because numbers can always change
  taxa_trace <- taxa_traits
  taxa_traits <- taxa_traits[ , -which( names( taxa_traits ) %in% "Taxa") , drop = F]
  # remove NA
  taxa_traits <- taxa_traits[complete.cases(taxa_traits), ]
  # remove categories with sum = 0
  taxa_traits <- taxa_traits[ , colSums(taxa_traits) > 0 , drop = F ]
  # remove rows also in abundances
  abundances <- abundances[ rownames(taxa_traits) , ]


  qual_fs <- qfs(taxa_traits , nbdim = nbdim, metric = metric, corr_method = corr_method)
  m <- qual_fs$meanSD<0.01

  if(!any( m == TRUE)) { stop("there is no optimal number of dimension, please check your data for possible problems")}
  m <- min( which( m == TRUE ) ) + 1

  fric <- fric_3d( t(abundances[, -which( names( abundances ) %in% "Taxa")]), fpc = qual_fs$fpc, m = m)
  fric[ which( is.na( fric ) ) ]<-0
  if( traceB == F ){
    fric
  }
  else{
    list( fric, taxa_trace )
  }
}

