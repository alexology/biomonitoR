# biomonitoR
A package for calculating indices for biomonitoring of running water with a focus on macroinvertebrate community. 
Still in development. See www.biomonitor.it for tips.


# Installation

```R
# install the devtools package and then
library(devtools)

install_github("alexology/biomonitoR")
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
genNumb(data.agR)
famNumb(data.agR)

# calculate shannon index
shannon(data.agR, taxLev = "Family")

# calculate italian bmwp and aspt
bmwp(data.agR, method="ita_agg")
aspt(data.agR, method="ita_agg")

# calculate iberian bmwp and aspt
bmwp(data.agR, method="spa")
aspt(data.agR, method="spa")

```

# To do list
1. Check and fill missing data in reference taxa list;
2. Update taxonomy of the reference datasets (e.g. those for the calculation of BMWP, LIFE, etc.) according to the [freshwaterecology.info](http://www.freshwaterecology.info/) list;
3. Implement other indices (e.g. metrics of I2M2, PSI, SPEAR, etc.);
4. Implement combTaxa, a function to combine taxa of a particular taxonomic level in order to find the best combination of taxa that fit with a target environmental variable;
5. Improve help;
6. Test biomonitoR extensively.

## Acknowledgments
This package is based upon work from COST Action CA15113 (SMIRES, Science and Management of Intermittent Rivers and Ephemeral Streams,[www.smires.eu](http://www.freshwaterecology.info/)), supported by COST (European Cooperation in Science and Technology).

## References
Reference database is from [freshwaterecology.info](http://www.smires.eu/).

Schmidt-Kloiber, A. & Hering D. (2015): www.freshwaterecology.info - an online tool that unifies, standardises and codifies more than 20,000 European freshwater organisms and their ecological preferences. Ecological Indicators 53: 271-282. doi: 10.1016/j.ecolind.2015.02.007
