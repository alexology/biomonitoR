#' @importFrom ade4 mantel.randtest is.euclid
#' @importFrom stats cmdscale

pcoaQuality <- function( x , type = "none" , method = NULL , tresh = 0.7 , nbdim = 15 ){

  if( identical( method , "cor" ) ){
    Order <- attr( x ,"Size" )

    if( is.euclid( x ) ){
      res.pco.c <- suppressWarnings( cmdscale( x , k = Order - 1 , eig = TRUE ) )
      res.eig <- res.pco.c$eig
      r2.eig <- cumsum( res.eig ) / sum( res.eig )
      res.pco <- res.pco.c$points
      res.man.r <- res.man.p <- c()

      for( i in 1:ncol( res.pco ) ){
        res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
        res.man.r[ i ] <- res.man$obs
        res.man.p[ i ] <- res.man$pvalue
      }

      r2.ax <- which( res.man.r > tresh )[ 1 ]
      res <- data.frame( Transformation = "none" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] ,  p_value = res.man.p[ r2.ax ] , r2 =  r2.eig[ r2.ax ] )
    }


    if( ! is.euclid( x ) ){
      cai.dist <- cailliez( x )
      lin.dist <- lingoes( x )
      sqr.dist <- sqrt( x )
      qua.dist <- quasieuclid( x )

      if( identical( type , "none" ) ){
        res.pco.c <- suppressWarnings( cmdscale( x , k = Order - 1 , eig = TRUE ) )
        res.pco <- res.pco.c$points
        res.eig <- res.pco.c$eig
        r2.eig <- cumsum( ( res.eig ) )
        neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
        r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- which( res.man.r > tresh )[ 1 ]
        res <- data.frame( Transformation = "none" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ] , r2 =  r2.eig[ r2.ax ]  )

      }

      if( identical( type , "cailliez" ) ){
        res.pco.c <- suppressWarnings( cmdscale( cai.dist , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco.c$eig
        res.pco <- res.pco.c$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        tresh.max <- max( res.man.r[ 1 ] , tresh )
        r2.ax <- length( which( res.man.r  <= tresh.max ) )

        if( is.euclid( cai.dist ) ){

          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          res <- data.frame( Transformation = "cailliez" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ] , r2 = r2.eig[ r2.ax ] )

        } else {

          res.eig <- res.pco.c$eig
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          res <- data.frame( Transformation = "cailliez" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ]  , r2 = r2.eig[ r2.ax ] )

        }
      }

      if( identical( type , "lingoes" ) ){
        res.pco.c <- suppressWarnings( cmdscale( lin.dist , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco.c$eig
        res.pco <- res.pco.c$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- which( res.man.r > tresh )[ 1 ]

        if( is.euclid( lin.dist ) ){

          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          res <- data.frame( Transformation = "lingoes" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ] , r2 = r2.eig[ r2.ax ]  )

        } else {

          res.eig <- res.pco.c$eig
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          res <- data.frame( Transformation = "lingoes" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ] , r2 = r2.eig[ r2.ax ] )

        }
      }

      if( identical( type , "sqrt" ) ){
        res.pco.c <- suppressWarnings( cmdscale( sqr.dist , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco.c$eig
        res.pco <- res.pco.c$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- which( res.man.r > tresh )[ 1 ]

        if( is.euclid( sqr.dist ) ){

          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          res <- data.frame( Transformation = "sqrt" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ] , r2 = r2.eig[ r2.ax ]   )

        } else {

          res.eig <- res.pco.c$eig
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          res <- data.frame( Transformation = "sqrt" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ] , r2 = r2.eig[ r2.ax ] )

        }
      }

      if( identical( type , "quasi" ) ){
        res.pco.c <- suppressWarnings( cmdscale( qua.dist , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco.c$eig
        res.pco <- res.pco.c$points
        res.man.r <- res.man.p <- c()

        for( i in 1:ncol( res.pco ) ){
          res.man <- mantel.randtest( dist( res.pco[ , 1:i ] ) , x )
          res.man.r[ i ] <- res.man$obs
          res.man.p[ i ] <- res.man$pvalue
        }

        r2.ax <- which( res.man.r > tresh )[ 1 ]

        if( is.euclid( qua.dist ) ){

          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          res <- data.frame( Transformation = "quasi" , Euclidean = "yes" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ] , r2 = r2.eig[ r2.ax ]  )

        } else {

          res.eig <- res.pco.c$eig
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          res <- data.frame( Transformation = "quasi" , Euclidean = "no" , nbdim = r2.ax , cor = res.man.r[ r2.ax ] , p_value = res.man.p[ r2.ax ] , r2 = r2.eig[ r2.ax ]  )

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
      r2.ax <- which( r2.eig > tresh )[ 1 ]
      res <- data.frame( Transformation = "none" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ] )
    }


    if( ! is.euclid( x ) ){
      cai.dist <- cailliez( x )
      lin.dist <- lingoes( x )
      sqr.dist <- sqrt( x )
      qua.dist <- quasieuclid( x )

      if( identical( type , "none" ) ){
        res.pco <- suppressWarnings( cmdscale( x , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco$eig
        r2.eig <- cumsum( ( res.eig ) )
        neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
        r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
        r2.ax <- which( r2.eig > tresh )[ 1 ]
        res <- data.frame( Transformation = "none" , Euclidean = "no" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]   )
      }

      if( identical( type , "cailliez" ) ){
        res.pco <- suppressWarnings( cmdscale( cai.dist , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco$eig

        if( is.euclid( cai.dist ) ){

          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          r2.ax <- which( r2.eig > tresh )[ 1 ]
          res <- data.frame( Transformation = "cailliez" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ] )

        } else {

          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          r2.ax <- which( r2.eig > tresh )[ 1 ]
          res <- data.frame( Transformation = "cailliez" , Euclidean = "no" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ] )
        }
      }

      if( identical( type , "lingoes" ) ){
        res.pco <- suppressWarnings( cmdscale( lin.dist , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco$eig

        if( is.euclid( lin.dist ) ){
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          r2.ax <- which( r2.eig > tresh )[ 1 ]
          res <- data.frame( Transformation = "lingoes" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ] )
        } else {
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          r2.ax <- which( r2.eig > tresh )[ 1 ]
          res <- data.frame( Transformation = "lingoes" , Euclidean = "no" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ] )
        }
      }

      if( identical( type , "sqrt" ) ){
        res.pco <- suppressWarnings( cmdscale( sqr.dist , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco$eig

        if( is.euclid( sqr.dist ) ){
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          r2.ax <- which( r2.eig > tresh )[ 1 ]
          res <- data.frame( Transformation = "sqrt" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]   )
        } else {
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          r2.ax <- which( r2.eig > tresh )[ 1 ]
          res <- data.frame( Transformation = "sqrt" , Euclidean = "no" , nbdim = NA , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]  )
        }
      }

      if( identical( type , "quasi" ) ){
        res.pco <- suppressWarnings( cmdscale( qua.dist , k = Order - 1 , eig = TRUE ) )
        res.eig <- res.pco$eig

        if( is.euclid( qua.dist ) ){
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          r2.ax <- which( r2.eig > tresh )[ 1 ]
          res <- data.frame( Transformation = "quasi" , Euclidean = "yes" , nbdim = r2.ax , r2 = r2.eig[ r2.ax ]   )
        } else {
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          r2.ax <- which( r2.eig > tresh )[ 1 ]
          res <- data.frame( Transformation = "quasi" , Euclidean = "no" , nbdim = NA , nbdim = r2.ax , r2 = r2.eig[ r2.ax ] )
        }
      }

    }

  }

  if( identical( method , "maire" ) ){
    Order <- attr( x ,"Size" )

    if( is.euclid( x ) ){
      qual_fs <- qfs( x , nbdim = nbdim )
      r2.ax <- which( qual_fs$meanSD < tresh )[ 1 ] + 1

      if( is.na( r2.ax ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
      res.eig <- qual_fs$mat_eig
      r2.eig <- cumsum( res.eig ) / sum( res.eig )
      res <- data.frame( Transformation = "none" , Euclidean = "yes" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ] )
      rownames( res ) <- NULL
    }


    if( ! is.euclid( x ) ){
      cai.dist <- cailliez( x )
      lin.dist <- lingoes( x )
      sqr.dist <- sqrt( x )
      qua.dist <- quasieuclid( x )

      if( identical( type , "none" ) ){
        qual_fs <- qfs( x , nbdim = nbdim )
        r2.ax <- which( qual_fs$meanSD < tresh )[ 1 ] + 1

        if( is.na( r2.ax ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
        res.eig <- suppressWarnings( cmdscale( x , k = Order - 1 , eig = TRUE ) )$eig
        r2.eig <- cumsum( ( res.eig ) )
        neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
        r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
        res <- data.frame( Transformation = "none" , Euclidean = "no" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax ] , r2 = r2.eig[ r2.ax ]  )
        rownames( res ) <- NULL
      }

      if( identical( type , "cailliez" ) ){

        qual_fs <- qfs( cai.dist , nbdim = nbdim )
        r2.ax <- which( qual_fs$meanSD < tresh )[ 1 ] + 1
        if( is.na( r2.ax ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
        res.eig <- suppressWarnings( cmdscale( cai.dist , k = Order - 1 , eig = TRUE ) )$eig

        if( is.euclid( cai.dist ) ){
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          res <- data.frame( Transformation = "cailliez" , Euclidean = "yes" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ]   )
          rownames( res ) <- NULL

        } else {
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          res <- data.frame( Transformation = "cailliez" , Euclidean = "no" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ]   )
          rownames( res ) <- NULL

        }
      }

      if( identical( type , "lingoes" ) ){

        qual_fs <- qfs( lin.dist , nbdim = nbdim )
        r2.ax <- which( qual_fs$meanSD < tresh )[ 1 ] + 1
        if( is.na( r2.ax ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
        res.eig <- suppressWarnings( cmdscale( lin.dist , k = Order - 1 , eig = TRUE ) )$eig


        if( is.euclid( lin.dist ) ){
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          res <- data.frame( Transformation = "lingoes" , Euclidean = "yes" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ]   )
          rownames( res ) <- NULL
        } else {
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          res <- data.frame( Transformation = "lingoes" , Euclidean = "no" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ]  )
          rownames( res ) <- NULL
        }
      }

      if( identical( type , "sqrt" ) ){
        qual_fs <- qfs( sqr.dist , nbdim = nbdim )
        r2.ax <- which( qual_fs$meanSD < tresh )[ 1 ] + 1
        if( is.na( r2.ax ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
        res.eig <- suppressWarnings( cmdscale( sqr.dist , k = Order - 1 , eig = TRUE ) )$eig

        if( is.euclid( sqr.dist ) ){
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          res <- data.frame( Transformation = "sqrt" , Euclidean = "yes" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ]   )
          rownames( res ) <- NULL
        } else {
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          res <- data.frame( Transformation = "sqrt" , Euclidean = "no" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ]  )
          rownames( res ) <- NULL
        }
      }

      if( identical( type , "quasi" ) ){

        qual_fs <- qfs( qua.dist , nbdim = nbdim )
        r2.ax <- which( qual_fs$meanSD < tresh )[ 1 ] + 1
        if( is.na( r2.ax ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
        res.eig <- suppressWarnings( cmdscale( qua.dist , k = Order - 1 , eig = TRUE ) )$eig


        if( is.euclid( qua.dist ) ){
          r2.eig <- cumsum( res.eig ) / sum( res.eig )
          res <- data.frame( Transformation = "quasi" , Euclidean = "yes" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ]   )
          rownames( res ) <- NULL
        } else {
          r2.eig <- cumsum( ( res.eig ) )
          neg.eig <- max( abs( res.eig[ which( res.eig < 0 ) ] ) )
          r2.eig  <-  ( r2.eig + ( 1:Order ) * neg.eig ) / ( sum( res.eig ) + ( Order - 1  ) * neg.eig )
          res <- data.frame( Transformation = "quasi" , Euclidean = "no" , nbdim = r2.ax , SD = qual_fs$meanSD[ r2.ax - 1 ] , r2 = r2.eig[ r2.ax ]   )
          rownames( res ) <- NULL
        }
      }

    }

  }

  res

}
