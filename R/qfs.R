#' @importFrom stats cmdscale dist
#' @importFrom FD gowdis
#' @importFrom geometry convhulln
 

qfs <- function(x, nbdim = nbdim, metric = metric, corr_method = corr_method){
  if( anyNA(x) ) {  warning(" NA detected in 'x' ")     }
  if ( nbdim < 2)   {  stop(" 'nbdim' must be higher than 1")     }
  if ( nrow( x ) < 3 )   {  stop(" there must be at least 3 species in 'x' ")     }
  if ( metric == "Euclidean" & ! any( apply( x , 2 , is.numeric ) ) )  { stop("using Euclidean distance requires that all traits are continuous")   }
  if ( nbdim > ncol(x) )  { stop( paste( "using", metric, "distance requires less dimensions than number of traits" ) )   }
  
  # computing functional dissimilarity between species given their traits values
  if (metric == "Gower") mat_dissim <- gowdis( x , ord = "classic" )
  if (metric == "Euclidean") mat_dissim <- dist( scale( x ) ) # scaling if continuous traits
  
  # lists to store distances matrices
  dist_raw<-list()
  dist_st<-list()
  
  ################################
  # computing PCoA using Caillez correction
  mat_cmd <- cmdscale(mat_dissim , k = nbdim, add = corr_method, eig = TRUE)

  # changing number of dimensions given number of positive eigenvalues
  nbdim <- min(nbdim, sum( mat_cmd$eig > 0) )
  
  # keeping species coordoinates on the 'nbdim' axes
  mat_coord <- mat_cmd$points[ , 1:nbdim ]
  row.names( mat_coord ) <- row.names( x )
  colnames( mat_coord ) <- paste( "PC" , 1:nbdim, sep = "" )
  
  # computing Euclidean distances between species in the (nbdim-1) multidimensionnal functional spaces 
  for ( k in 2:nbdim ) {
    eval( parse( text = paste( "dist_", k , "D<-dist(mat_coord[,1:", k ,"], method='euclidean')", sep = "" ) ) )
    eval( parse( text = paste( "dist_raw$m_", k , "D<-dist_", k , "D" , sep = "" ) ) )
  } # end of k
  
  
  ################################
  # computing mean squared deviation between initial distance and standardized final distance in the functional space
  meanSD <- rep( NA , nbdim-1 ) ; names( meanSD ) <- c(paste(paste("m_", 2:nbdim , "D" , sep = "" ) ) )
  
  z <- mat_dissim # initial distance
  S <- nrow(x) # species richness
  
  # for muldimensionnal spaces
  for ( k in 2:nbdim )  {
    eval( parse( text = paste( "y<-dist_" , k , "D" , sep="") ) )
    yst <- y / max( y ) * max( z )
    eval( parse( text = paste( "dist_st$m_",k , "D<-dist_" , k , "D" , sep = "" ) ) )
    meanSD[paste("m_",k,"D",sep="")]<-round( ( (sum((z-yst)^2)) / (S*(S-1)/2) ) ,6)
  }  # end of k
  
  # list of outputs
  res <- list( meanSD = meanSD , mat_dissim = mat_dissim, fpc = mat_coord )
  
  invisible( res )

}

fric_3d <- function( taxa , fpc , m){
  fric.3d <- rep( NA , nrow( taxa ) )
  convhulln( fpc[ , 1:m ], "FA" )$vol-> fric.3d.max
  apply( taxa > 0, 1, sum ) -> ric
  for ( com in 1:nrow( taxa ) ){
    fpc[ which( unlist( rep( taxa[ com ,] ) ) > 0 ), 1:m ] -> tr.com
    if (ric[ com ]> m + 1 ) convhulln( tr.com , "FA" )$vol/fric.3d.max -> fric.3d[ com ] else NA -> fric.3d[ com ]
  }
  return( fric.3d )
}