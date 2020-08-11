#' @importFrom ade4 mantel.randtest is.euclid
#' @importFrom stats cmdscale

pcoaQuality <- function( x , type = "none" , method = NULL , tresh = 0.7 , nbdim = 15 ){

  if( identical( method , "cor" ) ){
    Order <- attr( x ,"Size" )

    if( is.euclid( x ) ){
      res.pco <- suppressWarnings( cmdscale( x , k = Order - 1 , eig = TRUE ) )$points
      res.man.r <- res.man.p <- c()

      for( i in 1:ncol( res.pco ) ){
        res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
        res.man.r[ i ] <- res.man$obs
        res.man.p[ i ] <- res.man$pvalue
      }

      r2.ax <- length( which( res.man.r  <= tresh ) )
      res <- data.frame( Transformation = "none" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] ,  p_value = res.man.p[ r2.ax ]  )
    }


    if( ! is.euclid( x ) ){
      cai.dist <- cailliez( x )
      lin.dist <- lingoes( x )
      sqr.dist <- sqrt( x )
      qua.dist <- quasieuclid( x )

      if( identical( type , "none" ) ){
        res.pco <- suppressWarnings( cmdscale( x , k = Order - 1 , eig = TRUE ) )$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- length( which( res.man.r  <= tresh ) )
        res <- data.frame( Transformation = "none" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )
      }

      if( identical( type , "cailliez" ) ){
        res.pco <- suppressWarnings( cmdscale( cai.dist , k = Order - 1 , eig = TRUE ) )$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- length( which( res.man.r  <= tresh ) )
        if( is.euclid( cai.dist ) ){
          res <- data.frame( Transformation = "cailliez" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )

        } else {
          res <- data.frame( Transformation = "cailliez" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )

        }
      }

      if( identical( type , "lingoes" ) ){
        res.pco <- suppressWarnings( cmdscale( lin.dist , k = Order - 1 , eig = TRUE ) )$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- length( which( res.man.r  <= tresh ) )
        if( is.euclid( lin.dist ) ){
          res <- data.frame( Transformation = "lingoes" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )

        } else {
          res <- data.frame( Transformation = "lingoes" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )

        }
      }

      if( identical( type , "sqrt" ) ){
        res.pco <- suppressWarnings( cmdscale( sqr.dist , k = Order - 1 , eig = TRUE ) )$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- length( which( res.man.r  <= tresh ) )
        if( is.euclid( sqr.dist ) ){
          res <- data.frame( Transformation = "sqrt" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )

        } else {
          res <- data.frame( Transformation = "sqrt" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )

        }
      }

      if( identical( type , "quasi" ) ){
        res.pco <- suppressWarnings( cmdscale( qua.dist , k = Order - 1 , eig = TRUE ) )$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- length( which( res.man.r  <= tresh ) )
        if( is.euclid( qua.dist ) ){
          res <- data.frame( Transformation = "quasi" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )

        } else {
          res <- data.frame( Transformation = "quasi" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  )

        }
      }

    }

  }




  if( identical( method , "legendre" ) ){
    Order <- attr( x ,"Size" )

    if( is.euclid( x ) ){
      res.pco <- suppressWarnings( cmdscale( x , k = Order - 1 , eig = TRUE ) )
      res.eig <- res.pco$eig
      r2.eig <- cumsum( res.eig ) / sum( res.eig )
      r2.ax <- length( which( r2.eig  <= tresh ) )
      res <- data.frame( Transformation = "none" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ] )
    }


    if( ! is.euclid( x ) ){
      cai.dist <- cailliez( x )
      lin.dist <- lingoes( x )
      sqr.dist <- sqrt( x )
      qua.dist <- quasieuclid( x )

      if( identical( type , "none" ) ){
        res.eig <- res.pco$eig
        r2.eig <- cumsum( ( res.eig ) )
        neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
        r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
        r2.ax <- length( which( r2.eig  <= tresh ) )
        res <- data.frame( Transformation = "none" , Euclidean = "no" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]   )
      }

      if( identical( type , "cailliez" ) ){
        if( is.euclid( cai.dist ) ){
          res.pco <- suppressWarnings( cmdscale( cai.dist , k = Order - 1 , eig = TRUE ) )
          res.eig <- res.pco$eig
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          r2.ax <- length( which( r2.eig  <= tresh ) )
          res <- data.frame( Transformation = "cailliez" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]   )
        } else {
          res <- data.frame( Transformation = "cailliez" , Euclidean = "yes" , nbdim = NA  )
        }
      }

      if( identical( type , "lingoes" ) ){
        if( is.euclid( lin.dist ) ){
          res.pco <- suppressWarnings( cmdscale( lin.dist , k = Order - 1 , eig = TRUE ) )
          res.eig <- res.pco$eig
          res.eig <- res.pco$eig
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          r2.ax <- length( which( r2.eig  <= tresh ) )
          res <- data.frame( Transformation = "lingoes" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]   )
        } else {
          res <- data.frame( Transformation = "lingoes" , Euclidean = "no" , nbdim = NA  )
        }
      }

      if( identical( type , "sqrt" ) ){
        if( is.euclid( sqr.dist ) ){
          res.pco <- suppressWarnings( cmdscale( sqr.dist , k = Order - 1 , eig = TRUE ) )
          res.eig <- res.pco$eig
          res.eig <- res.pco$eig
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          r2.ax <- length( which( r2.eig  <= tresh ) )
          res <- data.frame( Transformation = "sqrt" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]   )
        } else {
          res <- data.frame( Transformation = "sqrt" , Euclidean = "no" , nbdim = NA  )
        }
      }

      if( identical( type , "quasi" ) ){
        if( is.euclid( qua.dist ) ){
          res.pco <- suppressWarnings( cmdscale( qua.dist , k = Order - 1 , eig = TRUE ) )
          res.eig <- res.pco$eig
          res.eig <- res.pco$eig
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          r2.ax <- length( which( r2.eig  <= tresh ) )
          res <- data.frame( Transformation = "quasi" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]   )
        } else {
          res <- data.frame( Transformation = "quasi" , Euclidean = "no" , nbdim = NA  )
        }
      }

    }

  }

  if( identical( method , "maire" ) ){

    if( is.euclid( x ) ){
      qual_fs <- qfs( x , nbdim = nbdim )
      m.check <- qual_fs$meanSD < tresh

      if( ! any( na.omit( m.check ) ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
      r2.ax <- min( which( m.check ) )
      res <- data.frame( Transformation = "none" , Euclidean = "yes" , nbdim = r2.ax  )
    }


    if( ! is.euclid( x ) ){
      cai.dist <- cailliez( x )
      lin.dist <- lingoes( x )
      sqr.dist <- sqrt( x )
      qua.dist <- quasieuclid( x )

      if( identical( type , "none" ) ){
        qual_fs <- qfs( x , nbdim = nbdim )
        m.check <- qual_fs$meanSD < tresh
        if( ! any( na.omit( m.check ) ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
        r2.ax <- min( which( m.check ) )
        res <- data.frame( Transformation = "none" , Euclidean = "no" , nbdim = r2.ax  )
      }

      if( identical( type , "cailliez" ) ){
        if( is.euclid( cai.dist ) ){
          qual_fs <- qfs( cai.dist , nbdim = nbdim )
          m.check <- qual_fs$meanSD < tresh
          if( ! any( na.omit( m.check ) ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
          r2.ax <- min( which( m.check ) )
          res <- data.frame( Transformation = "cailliez" , Euclidean = "yes" , nbdim = r2.ax  )
        } else {
          res <- data.frame( Transformation = "cailliez" , Euclidean = "no" , nbdim = NA  )
        }
      }

      if( identical( type , "lingoes" ) ){
        if( is.euclid( lin.dist ) ){
          qual_fs <- qfs( lin.dist , nbdim = nbdim )
          m.check <- qual_fs$meanSD < tresh
          if( ! any( na.omit( m.check ) ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
          r2.ax <- min( which( m.check ) )
          res <- data.frame( Transformation = "lingoes" , Euclidean = "yes" , nbdim = r2.ax  )
        } else {
          res <- data.frame( Transformation = "lingoes" , Euclidean = "no" , nbdim = NA  )
        }
      }

      if( identical( type , "sqrt" ) ){
        if( is.euclid( sqr.dist ) ){
          qual_fs <- qfs( sqr.dist , nbdim = nbdim )
          m.check <- qual_fs$meanSD < tresh
          if( ! any( na.omit( m.check ) ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
          r2.ax <- min( which( m.check ) )
          res <- data.frame( Transformation = "sqrt" , Euclidean = "yes" , nbdim = r2.ax  )
        } else {
          res <- data.frame( Transformation = "sqrt" , Euclidean = "no" , nbdim = NA  )
        }
      }

      if( identical( type , "quasi" ) ){
        if( is.euclid( qua.dist ) ){
          qual_fs <- qfs( sqr.dist , nbdim = nbdim )
          m.check <- qual_fs$meanSD < tresh
          if( ! any( na.omit( m.check ) ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
          r2.ax <- min( which( m.check ) )
          res <- data.frame( Transformation = "quasi" , Euclidean = "yes" , nbdim = r2.ax  )
        } else {
          res <- data.frame( Transformation = "quasi" , Euclidean = "no" , nbdim = NA  )
        }
      }

    }

  }




  res

}
