2# biomonitoR 0.9.3

* Implemented 2 plot functions.
* Added the function `abundance_classes` to transform abundance-coverage data into abundance-coverage classes.
* New macroinvertebrate-based indices! Functions `dehli` and `flow_T` added.
* New functions for building reference datasets from online resources! Function `get_gbif_taxa_tree`, `get_iucn_taxa_tree`, `get_nbn_taxa_tree` and `get_worms_taxa_tree` added.
* New macrophyte-based index! Function `ibmr` added.

2# biomonitoR 0.9.2

* Function deprecated in version 0.9.1 are now defunct.
* Code coverage added.
* Bug fixed. `eptd_families` argument of `eptd` incorrectly managed taxa names.
* Bug fixed. `select_pcoa_axes` incorrectly calculated the SD value for `none` correction option of the method method `maire`. It overestimated the SD value by assigning the value of the true number of axes plus 1. The estimated number of axes was right.
* Citation file added.
* Added a `plot` utility for objects generated with `as_biomonitor`.
* Bug fixed. The option `zerodist_rm` of the `f_rich` function did not worked as expected when set to `TRUE`.

2# biomonitoR 0.9.1

* Function names changed according to the tidyverse style guide.
* Argument `fuzzy` added to the `add_bias_to_traits()` function. This option allows to add bias also to non-fuzzy data. Argument `traceB` removed.
* Code styled according to the tidyverse style guide with the `styler` package.
* Tests added with the testthat package.
* Bug fixed. `as_biomonitor` always imported binary data as an object of class `bin`. Users are free to do whatever they want with their data and this behavior has been removed.
* Bug fixed. `psi` incorrectly assigned site names.

2# biomonitoR 0.9

* Added a `NEWS.md` file to track changes to the package.

