#' Fuzzy coded functional richness
#'
#' This function calculates the functional richness based on trait categories.
#'
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
#'  \deqn{C = \sum_{j=1}^{N} \lambda_j \ p_j  : \ \lambda_j \geq \ for \ all \ j \ and \ \sum_{j=1}^{N} \lambda_j = 1}
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
#' @param dfref reference database to be used when a custom trait database is provided and agg equals to TRUE.
#' @param traitSel interactively select traits.
#' @param colB A vector that contains the number of modalities for each trait
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
#' \enumerate{
#'  \item **results**: results of the ffred function;
#'  \item **traits**: a data.frame containing the traits used for the calculations;
#'  \item **taxa**: a data.frame conaining the taxa used for th calculations;
#'  \item **nbdim**: number of dimensions used after calculatin the quality of functional spaces according to Maire et al., (2015);
#'  \item **taxa_not_used**: a vector conaining the names of the taxa exluded from the calculations. Taxa are excluded from the calculations if they have NA values at least in one of the trait modalities;
#'  \item **problematic traits**: traits containing NAs for *taxa_not_used*
#'  \item **traits_not_used**: trait modalities excluded from the calculations because they have 0 value for each species.
#' }
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
#' ffrich(data.agR, taxLev = "Taxa")
#'
#' @seealso [aggregatoR]
#'
#' @references Barber, C. B., Dobkin, D. P., & Huhdanpaa, H. (1996).
#'   The quickhull algorithm for convex hulls. ACM Transactions on
#'   Mathematical Software (TOMS), 22(4), 469-483.
#' @references Cornwell, W. K., Schwilk, D. W., & Ackerly, D. D. (2006).
#'   A trait-based test for habitat filtering: convex hull volume.
#'   Ecology, 87(6), 1465-1471
#' @references Maire, E., Grenouillet, G., Brosse, S., & Villeger, S. (2015).
#'   How many dimensions are needed to accurately assess functional diversity?
#'   A pragmatic approach for assessing the quality of functional spaces. Global
#'   Ecology and Biogeography, 24(6), 728-740.
#' @references Mason, N. W., Mouillot, D., Lee, W. G., and
#'   Wilson, J. B. (2005). Functional richness, functional evenness and functional
#'   divergence: the primary components of functional diversity. Oikos, 111(1),
#'   112-118.
#' @references Villeger, S., Mason, N. W., & Mouillot, D.
#'   (2008). New multidimensional functional diversity indices for a
#'   multifaceted framework in functional ecology. Ecology, 89(8), 2290-2301.
#'
#' @export

ffrich <- function(x, traitDB = NULL, agg = FALSE,  dfref = NULL, traitSel = FALSE, colB = NULL, taxLev = "Family", traceB = FALSE, nbdim = 7, metric = "Gower", corr_method = FALSE){

  # check if user provided a trait database, otherwise use traitsTachet
  # if traitsTachet has to be used check for class biomonitoR and "mi"
  if( is.null( traitDB )){
    # check if the object d is of class "biomonitoR" & "mi"
    classCheck(x, group = "mi")
    traitDB = traitsTachet
    if( is.null( colB ) ) { colB = c( 8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8 ) }
    else { stop("You must not set colB when traitDB = NULL") }
    # useful for the if condition later
    mes <- "no"
  } else {
    # check if the object d is of class "biomonitoR"
    classCheck(x)
    traitDB = traitDB
    # trim and capitalise the first letter in the DB provided by the user
    traitDB[ , "Taxa"] <- as.factor( sapply( trim( traitDB[ , "Taxa"] ), capWords, USE.NAMES = F ) )
    mes <- "yes"
    if( is.null( colB ) ) ( stop("Please provide colB") )
    # check if the number of traits in traitDB equals the sum of colB, otherwise stop
    if( ( nrow( traitDB ) - 1 ) != sum( colB ) ) ( stop("The number of traits in traitDB is not equal to the sum of colB") )
  }

  if( traitSel == TRUE){
    Index <- rep( 1:length( colB ), colB)
    rma <- select.list( names( traitDB[ -which( names( traitDB ) %in% "Taxa")] ) , title = "Traits selection"  , graphics = TRUE , multiple = T )
    # new colB based on user trait selection, -1 because there is the column called Taxa
    colB <- as.vector( table( Index[ which( names( traitDB ) %in% rma ) - 1 ] ) )

    #  trait must have at least two modalities
    if( any( colB < 2 ) ) ( stop( "a trait must have at least two modalities" ) )

    traitDB <- traitDB %>%
      select( c("Taxa", rma) )
    # trim and capitalise the column Taxa of the user' trait database
    traitDB$Taxa <- apply( as.data.frame( trim( traitDB$Taxa ) ) , 1 , capWords)

  }


  # check for taxLev: it needs to be Taxa, Species, Genus or Family
  if (! taxLev %in% c("Family", "Genus", "Species", "Taxa")) {
    return("taxLev should be one of the following: Family, Genus, Species or Taxa")
  }

  abundances <- x[[taxLev]]
  colnames(abundances)[1] <- "Taxa"
  st.names <- names( abundances[ , -which( "Taxa" %in% names( abundances ) ), drop = F ] )


  # remove unassigned taxa from abundances
  if("unassigned" %in% abundances[ , "Taxa"]){
    z <- which(abundances[ , "Taxa" ] == "unassigned")
    abundances <- abundances[ -z , ] # remove unassigned row from the species count
  }

  taxa <- as.character(abundances$Taxa)

  if (length(taxa[taxa != "unassigned"]) == 0) {
    return("At least one taxa should be identified at a level compatible with the indicated taxLev")
  }


  # create dummy variables to avoid R CMD check NOTES
  traitsTachet <- Taxa <- modality <- affinity <- Phylum <- Subspecies <-
    Abundance <- Sample <- Weight <- Affinity <- totWeight <-
    weightedAffinity <- Category <- . <- NULL

  # prepare the taxa trait database
  trait_db <- traitDB                               %>%
    (function(df) {
      mutate(df,
             Taxa = gsub(pattern     = "sp.|Ad.|Lv.|Gen.",
                         replacement = "",
                         x           = Taxa))
    })                                              %>%
    gather(key = modality, value = affinity, -Taxa) %>%
    group_by(Taxa, modality)                        %>%
    summarise(affinity = mean(affinity))            %>%
    spread(key = modality, value = affinity)        %>%
    ungroup()
  trait_db$Taxa <- trimws(trait_db$Taxa)

  if(mes == "no"){
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

    # merge reference database

    ref <- select(mi_ref, Phylum:Taxa)
  }

  if(mes == "yes"){
    if(agg == TRUE){
      if( is.null( dfref ) == TRUE) ( stop("Reference database is needed when agg = TRUE") )
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

      # merge reference database

      ref <- select(dfref, Phylum:Taxa)
    } else {
      ref <- select(x$Tree, Phylum:taxLev)
      ref <- ref[ !duplicated(ref) , ]
      ref <- ref[ref[taxLev] != "", ]
      ref$Taxa <- ref[ , taxLev]
      level <-  rep(taxLev, nrow( ref ) )
    }
  }

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
  # order the traits as in the original data.frame (tratiDB)
  taxa_traits <- taxa_traits[ , match( names( traitDB ), names( taxa_traits ) ) ]
  taxa_traits <- as.data.frame(taxa_traits)
  # be sure that taxa_traits contains only the Taxa present in the user's community data
  taxa_traits <- taxa_traits[ taxa_traits[, "Taxa"]  %in% abundances[, "Taxa"] , ]
  taxa_trace <- taxa_traits
  # remove NA from the rows. Sometimes happens that a trait for a species is set to NA
  taxa_traits <- taxa_traits[complete.cases(taxa_traits[ , -which( names( taxa_traits ) %in% "Taxa") , drop = F]), ]
  taxa_traits_name <- as.character(taxa_traits$Taxa)
  taxa_traits <- taxa_traits[ , -which( names( taxa_traits ) %in% "Taxa") , drop = F]

  # remove categories with sum = 0, we don't want traits equals to zero
  cl.rm <- colSums(taxa_traits) > 0
  taxa_traits <- taxa_traits[ , cl.rm , drop = F ]

  #remove traits with incomplete cases and sum = 0

  Index <- rep( 1:length( colB ), colB)
  colB <- as.vector( table( Index[ cl.rm ] ) )

  if( any( colB < 2 ) ) ( stop( "a trait must have at least two modalities" ) )

  tr_prep <- prep.fuzzy( taxa_traits, col.blocks = colB)
  # check for the problematic traits
  pr.tr <- names(tr_prep[ , is.na( colSums( tr_prep ) )])
  # remove rows also in abundances
  abundances <- abundances[ as.character(abundances$Taxa) %in% taxa_traits_name[ complete.cases( tr_prep ) ], ]
  abu.names <- abundances[ , "Taxa" ]
  tr_prep <- tr_prep[ complete.cases( tr_prep ), ]
  abundances <- abundances[ , -which( names( abundances ) %in% "Taxa") , drop = F ]
  qual_fs <- qfs( tr_prep , nbdim = nbdim, metric = metric, corr_method = corr_method )
  m <- qual_fs$meanSD<0.01

  if(!any( m == TRUE)) { stop("there is no optimal number of dimension, please check your data for possible problems")}
  m <- min( which( m == TRUE ) ) + 1

  fric <- fric_3d( t(abundances), fpc = qual_fs$fpc, m = m )
  fric[ which( is.na( fric ) ) ] <- 0
  names(fric) <- st.names
  if( traceB == F ){
    return( fric )
  }
  else{
    # list the organisms that have not been used for the calculation
    if( length( abu.names ) == length( taxa ) ) { ta.miss <- "none"} else {
      ta.miss <- taxa[ ! taxa %in% abu.names ]
    }
    # list the traits that have not been used for the calculation, because they summed to 0
    if( length( names(tr_prep) == length( names( traitDB  ) ) ) ) { tr.miss <- "none" } else {
      tr.n <- names( traitDB  ) # trait names in traitDB
      tr.s <- names( tr_prep ) # trait names used for the calculation
      tr.miss <- tr.n[ ! tr.n %in% tr.s ]
    }
    if( length( pr.tr ) == 0 ) { tr.pr <- "none" } else { tr.pr <- pr.tr }
    fric.list <- list( fric, data.frame( Taxa = abu.names, tr_prep  ) , data.frame( Taxa= abu.names, abundances ), m , ta.miss, tr.pr , tr.miss )
    names( fric.list ) <- c( "results" , "traits" , "taxa", "nbdim" , "taxa_not_used", "problematic_traits", "traits_not_used" )
    return ( fric.list )
  }
}
