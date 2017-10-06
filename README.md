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
bmwp(data.agR, method="b")
aspt(data.agR, method="b")

# calculate iberian bmwp and aspt
bmwp(data.agR, method="i")
aspt(data.agR, method="i")

```

## References
Reference database is from [freshwaterecology.info](http://www.freshwaterecology.info/).

Schmidt-Kloiber, A. & Hering D. (2015): www.freshwaterecology.info - an online tool that unifies, standardises and codifies more than 20,000 European freshwater organisms and their ecological preferences. Ecological Indicators 53: 271-282. doi: 10.1016/j.ecolind.2015.02.007
