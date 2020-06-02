#################################################################################################################################
## quality_funct_space_fromdist: R function for computing the quality of functional multidimensional functional spaces
##
##                  Given a species x trait matrix, the function computes the quality (i.e. mean squared-deviation between
##                  initial functional distance and standardized distance in the functional space) for all the
##                  multidimensional functional spaces from 2 to N dimensions (N selected by the user).
##
##  Original code by Eva Maire & Sébastien Villéger (sebastien.villeger@univ-montp2.fr) in:
##
##    Maire, E., Grenouillet, G., Brosse, S., & Villéger, S. (2015) How many dimensions are needed to
##    accurately assess functional diversity? A pragmatic approach for assessing the quality of functional
##    spaces. Global Ecology and Biogeography, 24, 728–740.
##
##
##  Code modified and adapted by Cayetano Gutiérrez-Cánovas (cayeguti@um.es)
##
##
##     INPUTS:
##
##      - "dist_funct" : Euclidean distance or Gower dissimilarity matrix based on a species x trait matrix                                                                                                                             ##
##      - "nbdim" : maximum number of dimensions for multidimensional functional spaces. By default, nbdim=7
##      		Final number of dimensions depends on the number of positive eigenvalues (after correction) obtained with PCoA
##
##      - "plot" : character string to set the name of the jpeg file for plots illustrating the quality of functional spaces
##                  NA means no plot
##
##     NB: 1/ high value for 'nbdim' can increase computation time
##         2/ if at least one trait is not numeric, 'metric' should equal 'gower'
##         3/ if metric=Euclidean, functional traits are scaled (mean=0, sd=1) before computing functional distances
##         4/ R libraries 'ape', 'clue', 'cluster', 'FD', 'geometry', 'gtools' are required
##
##
##      OUTPUTS: a list with
##
##      - $ meanSD : a vector with mean squared deviation values for all the functional spaces tested
##            names are 'm_kD' (with k=2:nbdim) for multidimensional spaces
##
##      - $ mat_coord : coordinates of species in the nbdim multidimensional functional space (PCoA)
##
##      - $ dist_raw : list of raw functional distance matrices based on multidimensional functional spaces (PCoA)
##            of different dimensions (2D until nbdim)
##
##      - $ dist_st : list of standardized functional distance matrices based on multidimensional functional spaces (PCoA)
##            of different dimensions (2D until nbdim)
##
#################################################################################################################################

# Modified version to evaluate the quality of the functional space from a functional distance

qfs <- function( dist_funct ,  nbdim = nbdim, plot = "quality_funct_space" )
{


  #################################################################################################################################

  # checking data
  if ( ! ( "dist" %in% class( dist_funct ) ) )   {  stop( " 'dist_funct' must be of class 'dist'" )     }
  if ( length( dist_funct ) < 3 )   {  stop( " there must be at least 3 species in 'dist_funct' ")     }
  if ( sum( is.na( dist_funct ) ) != 0 )   {  stop( " NA are not allowed in 'dist_funct' ")     }
  if ( nbdim < 2 )   {  stop( " 'nbdim' must be higher than 1" )     }

  # to avoid RCMD notes
  y <- NULL

  # functional distance
  mat_dissim <- dist_funct

  # species names
  nm_sp<-row.names( as.matrix(dist_funct))

  ################################################################

  suppressWarnings( mat_pcoa<- dudi.pco( mat_dissim , scannf = F , nf = nbdim ) )

  # changing number of dimensions given number of positive eigenvalues
  nbdim <-min( nbdim , ncol( mat_pcoa$li ) )

  # keeping species coordoinates on the 'nbdim' axes
  mat_coord <- mat_pcoa$li[ , 1:nbdim ]
  row.names( mat_coord ) <- nm_sp
  colnames( mat_coord ) <- paste( "PC" , 1:nbdim ,sep = "" )

  # lists to store distance matrices
  dist_raw<-list()
  dist_st<-list()

  # computing Euclidean distances between species in the (nbdim-1) multidimensionnal functional spaces
  for (k in 2:nbdim)
  {
    eval(parse(text=paste("dist_",k,"D<-dist(mat_coord[,1:",k,"],method='euclidean')", sep="")))
    eval(parse(text=paste("dist_raw$m_",k,"D<-dist_",k,"D", sep="")))
  } # end of k

  ################################################################
  # computing mean squared deviation between initial distance and standardized final distance in the functional space
  meanSD<-rep(NA,nbdim) ; names(meanSD)<-paste("m_",2:nbdim,"D",sep="")

  x<-mat_dissim # initial distance
  S<-length(nm_sp) # species richness

  # for muldimensionnal spaces
  for (k in 2:nbdim)
  {
    eval(parse(text=paste("y<-dist_",k,"D",sep="")))
    yst<- y/max(y) * max(x)
    eval(parse(text=paste("dist_st$m_",k,"D<-dist_",k,"D", sep="")))
    meanSD[paste("m_",k,"D",sep="")]<-round( ( (sum((x-yst)^2)) / (S*(S-1)/2) ) ,6)
  }  # end of k

  # plotting change in meanSD with increasing number of dimensions
  par(mfrow=c(1,1),mar=c(5,5,1,1))
  barplot(height=meanSD,names.arg=names(meanSD), xlab="Functional space", ylab= "Quality (Mean SD)",
          space=0, cex.names=0.7, col=c("red", rep("blue",nbdim-1) ) )

  abline(h=0.01,lwd=3,lty=2)

  # plotting quality of each functional space

  # list of outputs

  res <- list( meanSD = meanSD , mat_coord = mat_coord , dist_raw = dist_raw, dist_st = dist_st )

  ################################################################################################################################
  ################################################################################################################################

  invisible(res)

} # end of function quality_funct_space_fromdist

################################################################################################################################
################################################################################################################################


fric_3d <- function( taxa , fpc , m , prec = c( "Qt" , "QJ" ) ){
  fric.3d <- rep( NA ,nrow( taxa ) )
  convhulln( fpc[ , 1:m ] , c( "FA" , prec ) )$vol -> fric.3d.max
  apply( taxa , 1 , FUN = function( x ) sum( x > 0 ) ) -> ric
  for ( com in 1:nrow( taxa ) ){
    fpc[ which( unlist( rep( taxa[ com , ] ) ) > 0 ) , 1:m ] -> tr.com
    if( ric[ com ] >= m + 1 ) convhulln( tr.com , c( "FA" , prec ) )$vol / fric.3d.max -> fric.3d[ com ] else NA -> fric.3d[ com ]
  }
  return( fric.3d )
}



# fdisp_k() estimates the Functional Dispersion of a set of communties,
# based on a subset of taxa within the original space
#
# This function computes weithed-mean distance to the centroid
# of how each community in the functional space
#
# Modified from dbFD() function in FD package by E. Laliberté, P. Legendre,
# B. Shipley
#
# Laliberté, E. and P. Legendre (2010) A distance-based framework for
# measuring functional diversity from multiple traits. Ecology 91:299-305.
#
#
# Inputs:
# d: trait distance or dissimilarity matrix
# a: community data
# tax_sub= subset of taxa to be selected from the funcional space
# m: number of axes to select
#
# Output:
# a vector with the Funcional Dispersion of each community


fdisp_k<-function (d, a, m )

{
  tol = 1e-07
  if (!inherits(d, "dist"))
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels")))
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (ncol(a) != n)
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a)))
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a))
    stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).",
         "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0))
    warning("At least one community has zero-sum abundances (no species).",
            "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0))
    warning("At least one species does not occur in any community (zero total abundance across all communities).",
            "\n")
  if (any(is.na(d)))
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol)
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])),
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  if (m>0) pos<-c(pos[1:m],rep(F,length(pos)-m))

  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(avg.dist.cent)
}


# feve_k() estimates the Functional Eveness of a set of communties
# This function computes weithed-mean regularity in taxon the functional space
#
# Modified from dbFD() function in FD package by E. Laliberté, P. Legendre,
# B. Shipley
#
# Villéger, S., N. W. H. Mason and D. Mouillot (2008) New multidimensional
# functional diversity indices for a multifaceted framework in functional
# ecology. Ecology 89:2290-2301.
#
#
# Inputs:
# fpc: matrix with the functiona space
# taxa: community data
# m: number of axes to select
#
# Output:
# a vector with the Funcional Evenness of each community

feve_k <- function(fpc,taxa,m){

  rdt=1

  # Creating variables
  nrow(taxa)->c
  FEve <- rep(NA, c) ; names(FEve) <- row.names(taxa)
  nb.sp <- apply( taxa , 1 , FUN = function( x ) sum( x > 0 ) )

  # generating the taxonomic matrix arranged according to the replicated trait values
  tax.pool<-ncol(taxa)
  taxa.rep<-data.frame(matrix(NA,nrow(taxa),tax.pool*rdt))
  spp.list<-c(1:(tax.pool*rdt))

  for (spp in 1:tax.pool) {paste(rep("spp",rdt),spp,sep="")->spp.list[((spp-1)*rdt+1):(spp*rdt)]}

  colnames(taxa.rep)<-spp.list

  for (spp in 1:tax.pool){taxa.rep[,((spp-1)*rdt+1):(spp*rdt)]<-taxa[,spp]/rdt}

  # Estimating Functional Evenness for each community

  for (i in 1:c) {
    sppres <- which(taxa.rep[i, ] > 0)
    # number of species in the community
    S <- length(sppres)
    ab <- as.matrix(taxa.rep[i, sppres])
    # scaling of abundances
    abundrel <- ab / sum(ab)

    # selecting the c
    tr <- data.frame(fpc[sppres,1:m ])

    if (nb.sp[i] > 2) {
      tr.dist <- dist(tr)
      linkmst <- mst(tr.dist)
      mstvect <- as.dist(linkmst)
      abund2 <- matrix(0, nrow = S, ncol = S)
      for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] +  abundrel[r]
      abund2vect <- as.dist(abund2)
      EW <- rep(0, S - 1)
      flag <- 1
      for (mv in 1:((S - 1) * S/2)) {
        if (mstvect[mv] != 0) {
          EW[flag] <- tr.dist[mv]/(abund2vect[mv])
          flag <- flag + 1
        }
      }
      minPEW <- rep(0, S - 1)
      OdSmO <- 1/(S - 1)
      for (l in 1:(S - 1)) minPEW[l] <- min((EW[l]/sum(EW)),
                                            OdSmO)
      FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    } else FEve[i] <- NA
  }
  return(FEve)
}


