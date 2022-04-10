# biomonitoR

A package for managing taxonomic and functional information and for calculating indices for biomonitoring of running waters, with a focus on the macroinvertebrate community.

[![codecov](https://codecov.io/gh/alexology/biomonitoR/branch/main/graph/badge.svg?token=Ix3zzcWgko)](https://codecov.io/gh/alexology/biomonitoR)
[![Build Status](https://travis-ci.org/alexology/biomonitoR.svg?branch=main)](https://travis-ci.org/alexology/biomonitoR)



# Installation

```R
# install the devtools package and then
# The devtools package requires and updated version of Rtools.
library(devtools)

install_github("alexology/biomonitoR", ref = "main", build_vignettes = TRUE)
```

To install the newest version:

```R
# install the devtools package and then
library(devtools)

install_github("alexology/biomonitoR", ref = "develop", build_vignettes = TRUE)
```

To install old and deprecated version:

```R
# install the devtools package and then
library(devtools)

install_github("alexology/biomonitoR", ref = "old_version")
```


# Basic usage

```R
library(biomonitoR)

# load example dataset. biomonitoR package needs a dataset with taxa names in the first column called "Taxa" and samples on the columns. Take a look to macro_ex for an example:

data(macro_ex)

# Prepare data for the analysis.

data_bio <- as_biomonitor(macro_ex)
data_agr <- aggregate_taxa(data_bio)

# calculate genus and family richness
richness(data_agr, tax_lev = "Genus")
richness(data_agr, tax_lev = "Family")

# calculate shannon index
shannon(data_agr, tax_lev = "Family")

# calculate italian bmwp and aspt
bmwp(data_agr, method = "ita", agg = TRUE)
aspt(data_agr, method = "ita", agg = TRUE)

# calculate iberian bmwp and aspt
bmwp(data_agr, method = "spa")
aspt(data_agr, method = "spa")

```

## Vignettes

To have some more information on how to use biomonitoR take a look to the vignettes:

```R
# import data in biomonitoR
vignette("introduction", package = "biomonitoR")

# trait-based analysis
vignette("functional", package = "biomonitoR")

```


## Acknowledgments
This package is based upon work from COST Action CA15113 (SMIRES, Science and Management of Intermittent Rivers and Ephemeral Streams,[www.smires.eu](http://www.smires.eu/)), supported by COST (European Cooperation in Science and Technology).

## References
Reference database is from [freshwaterecology.info](http://www.freshwaterecology.info/).

Schmidt-Kloiber, A. & Hering D. (2015): www.freshwaterecology.info - an online tool that unifies, standardises and codifies more than 20,000 European freshwater organisms and their ecological preferences. Ecological Indicators 53: 271-282. doi: 10.1016/j.ecolind.2015.02.007
