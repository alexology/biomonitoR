---
title: "Trait-based analysis with biomonitoR"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{functional}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Trait-based analysis is a standard procedure in community ecology and represents a good complement to taxonomic information. Some efforts have also been made to include trait-based metrics into current biomonitoring practices. The `biomonitoR` package provides functions to manipulate trait data and to calculate functional metrics (e.g. richness, diversity). Some examples on how to deal with functional analysis in `biomonitoR` will be provided.

## Macroinvertebrates

For macroinvertebrates, `biomonitoR` provides by default the Tachet database (Tachet et al. 2010) available at [freshwaterecology.info](https://www.freshwaterecology.info/) website. Briefly, the Tachet database is a fuzzy-coded database with 21 categories and a total of 113 modalities. One of the first steps in trait-based analysis is to retrieve trait information for the taxa present in the taxonomic dataset. There are some time-consuming steps when retrieving traits for macroinvertebrates because of inconsistencies between taxonomic resolution of the taxonomic dataset and those of traits. `biomonitoR` accomplish this task with the function `traitScaling`.


```r
library(biomonitoR)

# import the example dataset in biomonitoR
# macro_ex includes taxa identified at species, genus and family level

data(macro_ex)

data_bio <- as_biomonitor(macro_ex)
data_agr <- aggregate_taxa(data_bio)


# use the traitScaling function

data_ts <- assign_traits( data_agr , tax_lev = "Taxa" )

# we can use traitScaling at the desired taxonomic level
# for example setting tax_lev = "Family"
# assign_traits does not allow taxonomic levels higher than families

```


`traitScaling` returns a `data.frame` reporting taxa names of the taxonomic dataset, the taxonomic level of each taxon of the taxonomic dataset, the name of the taxa whose traits have been associated to a taxon of the taxonomic dataset and the corresponding taxonomic level, the taxonomic distance between taxa of the taxonomic dataset and the one from trait dataset and traits. Taxonomic distance is calculated based on the taxonomic levels implemented in `biomonitoR`. There are some rules of aggregation that can be consulted on the help page of `traitScaling`.



```r
# data.ts has more rows because more entries have been found for one or more taxa.
# only the first two traits are showed to not fill the page with numbers

data_ts[ 1:5 , 1:7 ]
#>         Taxa Taxa_db       Traits_taxlev Traits_real Taxonomic_distance
#> 1 Acentrella   Genus  Acentrella sinaica     Species                  1
#> 2    Ancylus   Genus Ancylus fluviatilis     Species                  1
#> 3     Baetis   Genus              Baetis       Genus                  0
#> 4 Beraeamyia   Genus Beraeamyia squamosa     Species                  1
#> 5  Beraeidae  Family              Beraea       Genus                  3
#>   LONGITUDINAL_1 LONGITUDINAL_2
#> 1      0.0000000      0.2857143
#> 2      0.1538462      0.2307692
#> 3      0.1428571      0.1904762
#> 4      0.0000000      0.0000000
#> 5      0.2500000      0.2500000
```


The `traitScaling` function links available traits to a taxon by looking at the taxonomic tree. Multiple assignements for the same taxon can thus be obtained.



```r
# data_ts has more rows because more entries have been found for one or more taxa.

dim( macro_ex )
#> [1] 34  3
dim( data_ts )
#> [1]  71 118

# look at Beraeidae for which five taxa of the trait dataset have been associated
# to one taxon of the taxonomic dataset

data_ts[ data_ts$Taxa %in% "Beraeidae" , 1:5 ]
#>        Taxa Taxa_db       Traits_taxlev Traits_real Taxonomic_distance
#> 5 Beraeidae  Family              Beraea       Genus                  3
#> 6 Beraeidae  Family Beraeamyia squamosa     Species                  4
#> 7 Beraeidae  Family   Beraeodes minutus     Species                  4
#> 8 Beraeidae  Family Beraeodina palpalis     Species                  4
#> 9 Beraeidae  Family             Ernodes       Genus                  3
```

Taxonomic distance can be useful, for example when only exact matches are needed. The `filter_by_distance` option of the `traitScaling` function or simple subsets can be used.



```r
assign_traits( data_agr , filter_by_distance = 0 )[, 1:5 ]
#>              Taxa   Taxa_db  Traits_taxlev Traits_real Taxonomic_distance
#> 3          Baetis     Genus         Baetis       Genus                  0
#> 10         Caenis     Genus         Caenis       Genus                  0
#> 21         Cyphon     Genus         Cyphon       Genus                  0
#> 22    Electrogena     Genus    Electrogena       Genus                  0
#> 23          Elmis     Genus          Elmis       Genus                  0
#> 26         Esolus     Genus         Esolus       Genus                  0
#> 27        Gyrinus     Genus        Gyrinus       Genus                  0
#> 28 Habroleptoides     Genus Habroleptoides       Genus                  0
#> 29   Habrophlebia     Genus   Habrophlebia       Genus                  0
#> 30  Holocentropus     Genus  Holocentropus       Genus                  0
#> 31       Hydraena     Genus       Hydraena       Genus                  0
#> 32    Hydropsyche     Genus    Hydropsyche       Genus                  0
#> 33     Hydroptila     Genus     Hydroptila       Genus                  0
#> 43      Laccobius     Genus      Laccobius       Genus                  0
#> 45        Leuctra     Genus        Leuctra       Genus                  0
#> 62    Nebrioporus     Genus    Nebrioporus       Genus                  0
#> 63  Onychogomphus     Genus  Onychogomphus       Genus                  0
#> 64 Orthocladiinae Subfamily Orthocladiinae   Subfamily                  0
#> 66    Rhithrogena     Genus    Rhithrogena       Genus                  0
#> 68        Setodes     Genus        Setodes       Genus                  0
#> 71    Tanypodinae Subfamily    Tanypodinae   Subfamily                  0

data_ts[ data_ts$Taxonomic_distance == 0 , 1:5 ]
#>              Taxa   Taxa_db  Traits_taxlev Traits_real Taxonomic_distance
#> 3          Baetis     Genus         Baetis       Genus                  0
#> 10         Caenis     Genus         Caenis       Genus                  0
#> 21         Cyphon     Genus         Cyphon       Genus                  0
#> 22    Electrogena     Genus    Electrogena       Genus                  0
#> 23          Elmis     Genus          Elmis       Genus                  0
#> 26         Esolus     Genus         Esolus       Genus                  0
#> 27        Gyrinus     Genus        Gyrinus       Genus                  0
#> 28 Habroleptoides     Genus Habroleptoides       Genus                  0
#> 29   Habrophlebia     Genus   Habrophlebia       Genus                  0
#> 30  Holocentropus     Genus  Holocentropus       Genus                  0
#> 31       Hydraena     Genus       Hydraena       Genus                  0
#> 32    Hydropsyche     Genus    Hydropsyche       Genus                  0
#> 33     Hydroptila     Genus     Hydroptila       Genus                  0
#> 43      Laccobius     Genus      Laccobius       Genus                  0
#> 45        Leuctra     Genus        Leuctra       Genus                  0
#> 62    Nebrioporus     Genus    Nebrioporus       Genus                  0
#> 63  Onychogomphus     Genus  Onychogomphus       Genus                  0
#> 64 Orthocladiinae Subfamily Orthocladiinae   Subfamily                  0
#> 66    Rhithrogena     Genus    Rhithrogena       Genus                  0
#> 68        Setodes     Genus        Setodes       Genus                  0
#> 71    Tanypodinae Subfamily    Tanypodinae   Subfamily                  0
```

What to do next? `biomonitoR` offers a more advanced subsetting tool with the function `manageTraits`.


```r
# select the traits belonging to the nearest taxa that have a taxonomic level equal or finer than the target one

data_ts_plus <- manageTraits( data_ts , method = "nearest+" , traceB = TRUE )

# with the traceB option set to TRUE you can check the taxa removed because they did not meet the requirments
# Here Serratella ignita is at species level while the trait at genus level (Serratella)

data_ts_plus$taxa_excluded
#> [1] "Serratella ignita"
```

When you are satisfied with the list of taxa you can finalize your trait dataset with two strategies. The function `traitsMean` allow the computation of the average traits for each Taxon. It is currently designed only for fuzzy data. For continuos data the default `aggregate` function can be used. Alternatively, taxa can be chosen randomly when multiple matches are available.



```r
# take the average. If col_blocks is set to NULL the default values will be used.
# Please consider that the default value is correct only for the default trait dataset!
# only the first 5 rows and columns will be showed

data_ts_av <- traitsMean( data_ts )
data_ts_av[ 1:5 , 1:5 ]
#>         Taxa LONGITUDINAL_1 LONGITUDINAL_2 LONGITUDINAL_3 LONGITUDINAL_4
#> 1 Acentrella      0.0000000      0.2857143      0.4285714     0.28571429
#> 2    Ancylus      0.1538462      0.2307692      0.3076923     0.07692308
#> 3     Baetis      0.1428571      0.1904762      0.1904762     0.19047619
#> 4 Beraeamyia      0.0000000      0.0000000      0.0000000     0.00000000
#> 5  Beraeidae      0.2291667      0.2708333      0.1145833     0.03125000

# alternatively a random selection can be done
# only the first 5 rows and  columns will be showed

sampleTraits( data_ts )[ 1:5 , 1:5 ]
#> Warning: 'sampleTraits' è deprecata
#> Usare 'sample_traits' al suo posto.
#> Si veda help("Deprecated") e help("biomonitoR-deprecated").
#>         Taxa LONGITUDINAL_1 LONGITUDINAL_2 LONGITUDINAL_3 LONGITUDINAL_4
#> 1 Acentrella      0.0000000      0.2857143      0.4285714     0.28571429
#> 2    Ancylus      0.1538462      0.2307692      0.3076923     0.07692308
#> 3     Baetis      0.1428571      0.1904762      0.1904762     0.19047619
#> 4 Beraeamyia      0.0000000      0.0000000      0.0000000     0.00000000
#> 5  Beraeidae      0.2500000      0.2500000      0.1250000     0.12500000
```

Now we have obtained a trait dataset that can be used for calculating functional metrics.

### Functional indices

The calculation of functional dispersion, diversity, evenness, redundancy and richness are straightforward.















