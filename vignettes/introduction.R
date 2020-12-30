## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(biomonitoR)

# built-in example

data( macro_ex )

head( macro_ex )


## ----import-------------------------------------------------------------------
# import macroinvertebrate data
data.bio <- asBiomonitor( macro_ex , group = "mi" )

data.bio

## ----import bin---------------------------------------------------------------
# import macroinvertebrate data
data.bio.bin <- asBiomonitor( macro_ex , group = "mi" , FUN = bin )

data.bio.bin

## ----spell--------------------------------------------------------------------
# introduce an error into taxonomic data

macro_ex_wrong <- macro_ex
macro_ex_wrong$Taxa <- as.character( macro_ex_wrong$Taxa )
macro_ex_wrong[ macro_ex_wrong$Taxa %in% "Acentrella" , "Taxa" ] <- "Acentrela"

data.bio.wrong <- asBiomonitor( macro_ex_wrong , group = "mi" , traceB = TRUE )

# "Acentrela" has been removed and a suggestion is proposed
data.bio.wrong

## ----as.data------------------------------------------------------------------
# export taxonomic dataset and suggestions as data.frame

macro_ex.df <- as.data.frame( data.bio.wrong )
macro_ex_correction.df <- as.data.frame( data.bio.wrong , object = 2 )
macro_ex_correction.df


## ----aggrega------------------------------------------------------------------
data.agg <- aggregatoR( data.bio )


## ----subset-------------------------------------------------------------------
# select EPT Taxa (Ephemeroptera, Plecoptera and Trichoptera)

subset(data.bio, taxa = c( "Ephemeroptera" , "Plecoptera" , "Trichoptera" ) )

# select Trichoptera excluding the trichopteran family Hydropsychidae

tricho <- subset(data.bio, taxa = "Trichoptera" , exclude = "Hydropsychidae" )
tricho

# import tricho with the aggregatoR function
tricho.agg <- aggregatoR( tricho )


## ----refe---------------------------------------------------------------------
# import data

data( Tree )

# Tree is a standard taxonomic tree that needs to be transformed into a reference database
# We can use the function refFromTree

ref_custom <- refFromTree(Tree)

ref_custom



