# biomonitoR
A package for calculating indices for biomonitoring of running water with a focus on macroinvertebrate community. Still in development.

# Installation

```R
# install the devtools package and then
library(devtools)

install_github("alexology/biomonitoR")
```

To install the newest version (recommended):

```R
# install the devtools package and then
library(devtools)

install_github("alexology/biomonitoR" , ref = "develop")
```


# Usage

```R
library(biomonitoR)

# load example dataset. biomonitoR package needs a dataset with taxa names in the first column called "Taxa" and samples on the columns. Take a look to macro_ex for an example:

data(macro_ex)

# Prepare data for the analysis.
data(macro_ex)
data.bio <- asBiomonitor(macro_ex)
data.agR <- aggregatoR(data.bio)

# calculate genus and family richness
richness(data.agR, "Genus")
famNumb(data.agR, "Family")

# calculate shannon index
shannon(data.agR, taxLev = "Family")

# calculate italian bmwp and aspt
bmwp(data.agR, method="ita_agg")
aspt(data.agR, method="ita_agg")

# calculate iberian bmwp and aspt
bmwp(data.agR, method="spa")
aspt(data.agR, method="spa")

```

## Acknowledgments
This package is based upon work from COST Action CA15113 (SMIRES, Science and Management of Intermittent Rivers and Ephemeral Streams,[www.smires.eu](http://www.smires.eu/)), supported by COST (European Cooperation in Science and Technology).

## References
Reference database is from [freshwaterecology.info](http://www.freshwaterecology.info/).

Schmidt-Kloiber, A. & Hering D. (2015): www.freshwaterecology.info - an online tool that unifies, standardises and codifies more than 20,000 European freshwater organisms and their ecological preferences. Ecological Indicators 53: 271-282. doi: 10.1016/j.ecolind.2015.02.007
