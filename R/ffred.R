#' Fuzzy coded functional redundancy
#'
#' This function calculates the functional richness based on trait categories.
#'
#' Functional redundancy (FR) is measured as the difference between taxonomic
#' diversity and functional diversity (de Bello et al., 2007). It relates positively
#' to ecosystem stability, resistance and resilience (Hooper et al. 2005; Guillemot
#' et al., 2011).
#'
#' \deqn{ FR = D - Q}
#'
#' The Gini-Simpson index is used to quantify Taxonomic Diversity (D, which ranges
#' from 0 to 1) where pi is the proportion of the abundance of taxa i in a
#' biological community.
#'
#' \deqn{D = 1 - \sum_{i=1}^{S} p_i^{2}}
#'
#' Rao quadratic entropy (Q; Rao, 1982) was used to estimate Functional Diversity
#' because it has been considered more appropriate than other indices (Botta-Dukat,
#' 2005; Ricotta, 2005). In this formula (Q) dij is the dissimilarity (ranging
#' from 0 to 1), between species i and j based on a set of specified functional
#' traits. This index is standardized by the maximum value to constrain the values
#' within the range of 0-1. Rao index is estimated using presence or abundance data
#' and the Euclidean transformed version of the traits-based Gower dissimilarity
#' matrix. For this, Gower's dissimilarity index (which ranges from 0 to 1) is used
#' because it can deal with traits of different nature and measuring scales
#' (continuous, nominal, binary, ordinal, etc.; see Podani 1999 for more information)
#'
#' \deqn{Q = \sum_{i=1}^{S} \sum_{j=1}^{S} d_{ij} \ p_i \ p_j}
#'
#' Given that the concept of FR was originally developed to represent the number
#' of taxa contributing similarly to an ecosystem function (Walker 1992; Lawton and
#' Brown 1993; Rosenfeld, 2002), FR and therefore functional diversity should be
#' calculated using only effect traits. Effect traits are those biological
#' features that directly influence a specific function of the ecosystem
#' (e.g. productivity, nutrient cycling). See Schmera et al., 2017 for more
#' information about effect traits in aquatic invertebrate communities. Regarding
#' the interpretation of the results, when taxa within a community differ
#' completely in their functional traits, then Q = D and thus FR = 0. On the
#' other hand, when all taxa have identical functional traits, then Q = 0 and
#' FR = D, and when in addition the number of taxa is very large and they are
#' equally abundant, then D (and in this case FR) approaches 1 (Pillar et al., 2013).
#' Although the concept of FR could suggest that functionally similar species may
#' compensate for the loss or failure of others, there is evidence that ecosystems
#' need such redundancy to perform their functions efficiently and stably over time
#' (Rosenfeld 2002; Biggs et al. 2012). In fact, a decrease in FR could be dramatic
#' in non-redundant communities since the loss or replacement of one species could
#' lead to loss of unique traits or functions (Hooper et al. 2005), increasing
#' ecosystem vulnerability (Elmqvist et al. 2003).
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
#' @param colB A vector that contains the number of modalities for each trait
#' @param taxLev character string giving the taxonomic level used to retrieve
#' trait information. Possible levels are `"Taxa"`, `"Species"`, `"Genus"`,
#' `"Family"` as returned by the [aggregatoR] function.
#' @param dfref reference database to be used when a custom trait database is provided and agg equals to TRUE.
#' @param traceB if TRUE ffrich will return a list with 2 elements, the first being the ffrich values and the second the database used for the calculation. Useful to check missing taxa.
#'
#' @details Taxa with no traits are removed from both the trait and abundance databases.
#' @note USE WITH CAUTION, STILL IN DEVELOPMENT.
#' @return a data.frame with 3 columns: Gini-Simpson richness, rao quadratic entropy and functional redundancy.
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#' @importFrom stats complete.cases na.omit
#' @importFrom ade4 ktab.list.df dist.ktab prep.fuzzy divc quasieuclid is.euclid
#'
#' @examples
#' data(macro_ex)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' ffred(data.agR)
#'
#' @seealso [aggregatoR]
#'
#' @references Biggs, R., Schluter, M., Biggs, D., Bohensky, E. L., BurnSilver, S.,
#'   Cundill, G., ... & Leitch, A. M. (2012). Toward principles for enhancing the
#'   resilience of ecosystem services. Annual Review of Environment and Resources,
#'   37, 421-448.
#' @references Botta-Dukat, Z. (2005). Rao's quadratic entropy as a measure of
#'   functional diversity based on multiple traits. Journal of Vegetation Science,
#'   16(5), 533-540.
#' @references de Bello, F., Leps, J., Lavorel, S., & Moretti, M. (2007).
#'   Importance of species abundance for assessment of trait composition:
#'   an example based on pollinator communities. Community Ecology, 8(2), 163-170.
#' @references Elmqvist, T., Folke, C., Nystrom, M., Peterson, G.,
#'   Bengtsson, J., Walker, B., & Norberg, J. (2003). Response diversity, ecosystem
#'   change, and resilience. Frontiers in Ecology and the Environment, 1(9), 488-494.
#' @references Guillemot, N., Kulbicki, M., Chabanet, P., & Vigliola, L. (2011).
#'   Functional redundancy patterns reveal non-random assembly rules in a
#'  species-rich marine assemblage. PLoS One, 6(10), e26735.
#' @references Hooper, D. U., Chapin, F. S., Ewel, J. J., Hector, A., Inchausti,
#'   P., Lavorel, S., et al. (2005). Effects of biodiversity on ecosystem
#'   functioning: a consensus of current knowledge. Ecological Monographs,
#'   75(1), 3-35.
#' @references Lawton, J.H. & Brown, V.K. (1993) Redundancy in ecosystems.
#'   Biodiversity and Ecosystem Function (eds E.-D. Schulze & H.A. Mooney),
#'   pp. 255-270. Springer-Verlag, Berlin.
#' @references Pillar, V. D., Blanco, C. C., Muller, S. C., Sosinski, E. E.,
#'   Joner, F., & Duarte, L. D. (2013). Functional redundancy and stability
#'   in plant communities. Journal of Vegetation Science, 24(5), 963-974.
#' @references Podani, J. (1999). Extending Gower's general coefficient of
#'   similarity to ordinal characters. Taxon, 331-340.
#' @references Rao, C. R. (1982). Diversity and dissimilarity coefficients:
#'   a unified approach. Theoretical population biology, 21(1), 24-43.
#' @references Ricotta, C. (2005). A note on functional diversity measures.
#'   Basic and Applied Ecology, 6(5), 479-486.
#' @references Rosenfeld, J. S. (2002). Functional redundancy in ecology
#'   and conservation. Oikos, 98(1), 156-162.
#' @references Schmera, D., Heino, J., Podani, J., Eros, T., & Doledec, S. (2017).
#'   Functional diversity: a review of methodology and current knowledge in
#'   freshwater macroinvertebrate research. Hydrobiologia, 787(1), 27-44.
#' @references Walker, B. H. (1992). Biodiversity and ecological redundancy.
#'   Conservation biology, 6(1), 18-23.
#'
#' @export

ffred <- function(x, traitDB = NULL, agg = FALSE, dfref = NULL, traitSel = FALSE, colB = NULL, taxLev = "Family", traceB = FALSE){

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
  
  # since traitDB is sorted in alphabetical order we need to be sure that colB is sorted too
  colB.names <- data.frame( colbN = names( traitDB[ , -1 ] ) , colB = rep( 1:length( colB ), colB) )
  colB.temp <- colB.names[ match( names( trait_db[ , -1 ] ) , colB.names[ , "colbN" ] ), "colB" ]
  
  
  
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
  taxa_traits <- as.data.frame(taxa_traits)
  taxa_traits_name <- as.character(taxa_traits$Taxa)
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
  
  colB.table <- table( colB.temp[ cl.rm ] )
  colB <- as.vector( colB.table[ match( unique( colB.temp ), names( colB.table ) ) ] )

  if( any( colB < 2 ) ) ( stop( "a trait must have at least two modalities" ) )
  
  
  tr_prep <- prep.fuzzy( taxa_traits, col.blocks = colB)
  abundances <- abundances[ as.character(abundances$Taxa) %in% taxa_traits_name[ complete.cases( tr_prep ) ], ]
  tr_prep <- tr_prep[ complete.cases( tr_prep ), ]
  tr_ktab <- ktab.list.df( list( tr_prep ) )
  dist_tr <- dist.ktab( tr_ktab, "F", taxa_traits )
  
  abundances <- abundances[ , -which( names( abundances ) %in% "Taxa") , drop = F]
  
  tax_sim <- divc( abundances )$diversity
  if( is.euclid( dist_tr ) ){
    raoQ <- divc( abundances, dist_tr , scale = T)$diversity
  } else {
    raoQ <- divc( abundances, quasieuclid( dist_tr ), scale = T)$diversity
  }
  FRed <- tax_sim - raoQ
  
  res <- data.frame(GS_rich = tax_sim, raoQ = raoQ, fred = FRed)
  rownames( res ) <- st.names
  if( traceB == FALSE ){
    return( res )
  }
  if(traceB == TRUE){
    res.list <- list( res, data.frame( Taxa = taxa_traits_name, tr_prep ),
                      data.frame( Taxa = taxa_traits_name, abundances ) )
    names( res.list ) <- c( "results" , "traits" , "taxa" )
    return( res.list )
  }
  
}