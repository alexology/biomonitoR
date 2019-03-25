#' Fuzzy coded functional dispersion
#'
#' FDis is defined as the mean distance in multidimensional trait space of individual species to the centroid of all species (Laliberte & Legendre, 2010).
#'
#' @param x results of function aggregatoR
#' @param traitDB a trait data base with a column `Taxa` and the other columns
#'   containing the traits. Needed only if users use their own trait database.  If a trait has several modalities they should be
#'   named as follow: TRAIT_MODALITY.
#'
#'   By default, the data base used is the one from Tachet *et al* (2010) that
#'   can be retrieved from
#'   [freshwaterecology.info](https://www.freshwaterecology.info/) website
#'   (Schmidt-Kloiber & Hering, 2015).
#' @param agg should ffrich aggregate user's traitDB of higher taxonomic level? TRUE to aggregate, otherwise FALSE.
#'    For instance, if user's traitDB has both Halesus and Limnephilidae, ffred will aggregate traits value if agg = TRUE.
#' @param traitSel interactively select traits. See details for more info.
#' @param colB A vector that contains the number of modalities for each trait. Needed only if users use their own trait database.
#' @param taxLev character string giving the taxonomic level used to retrieve
#' trait information. Possible levels are `"Taxa"`, `"Species"`, `"Genus"`,
#' `"Family"` as returned by the [aggregatoR] function.
#' @param dfref reference database to be used when a custom trait database is provided and agg = TRUE.
#'    It should be the same reference database used in [asBiomonitor] when dfref = TRUE.
#' @param traceB if TRUE ffred will return a list with 2 elements, the first being the ffrich values and the second the database used for the calculation. Useful to check missing taxa.
#'
#' @details Taxa without traits assigned in the trait database are removed from both the trait and abundance databases.
#' @note USE WITH CAUTION, STILL IN DEVELOPMENT.
#' @return A vector with functional dispersion results.\cr
#' If traceB is set to TRUE a list is provided with:
#' \enumerate{
#'  \item **results**: results of the ffred function;
#'  \item **traits**: a data.frame containing the traits used for the calculations;
#'  \item **taxa**: a data.frame conaining the taxa used for th calculations;
#'  \item **taxa_not_used**: a vector conaining the names of the taxa exluded from the calculations. Taxa are excluded from the calculations if they have NA values at least in one of the trait modalities;
#'  \item **problematic traits**: traits containing NAs for *taxa_not_used*
#'  \item **traits_not_used**: trait modalities excluded from the calculations because they have 0 value for each species.
#' }
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#' @importFrom stats complete.cases na.omit
#' @importFrom ade4 ktab.list.df dist.ktab prep.fuzzy divc quasieuclid is.euclid
#' @importFrom FD fdisp
#'
#' @examples
#' data(macro_ex)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' ffdisp(data.agR)
#'
#' @seealso [aggregatoR]
#' @references Laliberte, E., & Legendre, P. (2010). A distance-based framework for measuring functional diversity from multiple traits. Ecology, 91(1), 299-305.
#' @export

ffdisp <- function(x, traitDB = NULL, agg = FALSE, dfref = NULL, traitSel = FALSE, colB = NULL, taxLev = "Taxa", traceB = FALSE){


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
    if( ( ncol( traitDB ) - 1 ) != sum( colB ) ) ( stop("The number of traits in traitDB is not equal to the sum of colB") )
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
             Taxa = gsub(pattern     = "sp[.]|Ad[.]|Lv[.]|Gen[.]|lv[.]|ad[.]|gen[.]",
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
  tr_ktab <- ktab.list.df( list( tr_prep ) )
  dist_tr <- dist.ktab( tr_ktab, "F" )

  abundances <- abundances[ , -which( names( abundances ) %in% "Taxa") , drop = FALSE ]

  res <- fdisp( d= dist_tr , a = t( as.matrix( abundances ) ) )$FDis

  names( res ) <- st.names
  if( traceB == FALSE ){
    return( res )
  }
  if(traceB == TRUE){
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

    res.list <- list( res, data.frame( Taxa = abu.names, tr_prep ),
                      data.frame( Taxa = abu.names, abundances ), ta.miss, tr.pr , tr.miss  )
    names( res.list ) <- c( "results" , "traits" , "taxa", "taxa_not_used", "problematic_traits", "traits_not_used" )
    return( res.list )
  }

}
