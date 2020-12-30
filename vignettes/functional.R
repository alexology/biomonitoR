## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(biomonitoR)

# import the example dataset in biomonitoR
# macro_ex includes taxa identified at species, genus and family level

data(macro_ex)

data.bio <- asBiomonitor(macro_ex)
data.agR <- aggregatoR(data.bio)


# use the traitScaling function

data.ts <- traitScaling( data.agR , taxLev = "Taxa" )

# we can use traitScaling at the desired taxonomic level
# for example setting taxLev = "Family"
# traitscaling does not allow taxonomic levels higher than families



## ----head traitScaling--------------------------------------------------------
# data.ts has more rows because more entries have been found for one or more taxa.
# only the first two traits are showed to not fill the page with numbers

data.ts[ 1:5 , 1:7 ]



## ----dim traitScaling---------------------------------------------------------
# data.ts has more rows because more entries have been found for one or more taxa.

dim( macro_ex )
dim( data.ts )

# look at Beraeidae for which five taxa of the trait dataset have been associated
# to one taxon of the taxonomic dataset

data.ts[ data.ts$Taxa %in% "Beraeidae" , 1:5 ]



## ----filter traitScaling------------------------------------------------------
traitScaling( data.agR , filter_by_distance = 0 )[, 1:5 ]

data.ts[ data.ts$Taxonomic_distance == 0 , 1:5 ]


## ----manage traitScaling------------------------------------------------------
# select the traits belonging to the nearest taxa that have a taxonomic level equal or finer than the target one

data.ts_plus <- manageTraits( data.ts , method = "nearest+" , traceB = TRUE )

# with the traceB option set to TRUE you can check the taxa removed because they did not meet the requirments
# Here Serratella ignita is at species level while the trait at genus level (Serratella)

data.ts_plus$taxa_excluded



## ----mean traitScaling--------------------------------------------------------
# take the average. If colB is set to NULL the default values will be used.
# Please consider that the default value is correct only for the default trait dataset!
# only the first 5 rows and columns will be showed

data.ts.av <- traitsMean( data.ts )
data.ts.av[ 1:5 , 1:5 ]

# alternatively a random selection can be done
# only the first 5 rows and  columns will be showed

sampleTraits( data.ts )[ 1:5 , 1:5 ]

## ----functional indices-------------------------------------------------------
colB <- c( 8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8 )

# functional dispersion
f_disp( data.agR , traitDB = data.ts.av , nbdim = 2 , type = "F" , colB = colB )

# functional diversity 
f_divs( data.agR , traitDB = data.ts.av , type = "F" , colB = colB )

# functional evenness
f_eve( data.agR , traitDB = data.ts.av , type = "F" , colB = colB )

# functional redundancy
f_red( data.agR , traitDB = data.ts.av , type = "F" , colB = colB )

# functional richness on 2 dimensions
f_rich( data.agR , traitDB = data.ts.av , nbdim = 2 , type = "F" , colB = colB )

## ----traceB-------------------------------------------------------------------
frich <- f_rich( data.agR , traitDB = data.ts.av , nbdim = 2 , type = "F" , colB = colB , traceB = TRUE )

# look at the traits with NA values
head( frich$NA_detection )

## ----cwm----------------------------------------------------------------------

# community weighted mean and community specialization index
# only the first 5 columns are showed

cwm(x = data.agR, traitDB = data.ts.av, taxLev = "Taxa", trans = log1p )[ , 1:5 ]
csi(x = data.agR, traitDB = data.ts.av, taxLev = "Taxa", trans = log1p )[ , 1:5 ]


## ----dist ma------------------------------------------------------------------
library( ade4 )

rownames( data.ts.av ) <- data.ts.av$Taxa
traits.prep <- prep.fuzzy( data.ts.av[ , -1 ], col.blocks = colB )

traits.dist <- ktab.list.df( list( traits.prep ) )
traits.dist <- dist.ktab( traits.dist , type = "F" )

f_rich( data.agR , traitDB = traits.dist , nbdim = 2 )

## ----pcoaAxes-----------------------------------------------------------------
selectPcoaAxes( traits.dist , method = "cor" , tresh = 0.7)
selectPcoaAxes( traits.dist , method = "legendre" , tresh = 0.65)
selectPcoaAxes( traits.dist , method = "maire" , tresh = 0.01)

## ----family-------------------------------------------------------------------
data_fam.ts <- traitScaling( data.agR , taxLev = "Family" )

data_fam.ts.av <- traitsMean( data_fam.ts )

# functional richness
f_rich( data.agR , traitDB = data_fam.ts.av , taxLev = "Family" , nbdim = 2 , type = "F" , colB = colB )


## ----macroph------------------------------------------------------------------
# load macrophyte data
data( oglio )

# transform to presence-absence data for simplicity
# ignore the warning
oglio.bin <- oglio
oglio.bin[ oglio.bin > 0 ] <- 1

# importing data in the biomonitoR format
oglio.asb <- asBiomonitor( oglio.bin , group = "mf" , FUN = bin )
oglio.agg <- aggregatoR( oglio.asb )


