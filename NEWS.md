2# biomonitoR 0.9.2

* Function deprecated in version 0.9.1 are now defunct.
* Code coverage added.
* Bug fixed. `eptd_families` argument of `eptd` incorrectly managed taxa names.

2# biomonitoR 0.9.1

* Function names changed according to the tidiverse style guide.
* Argument `fuzzy` added to the `add_bias_to_traits()` function. This option allows to add bias also to non-fuzzy data. Argument `traceB` removed.
* Code styled according to the tidyverse style guide with the `styler` package.
* Tests added with the testthat package.
* Bug fixed. `as_biomonitor` always imported binary data as an object of class `bin`. Users are free to do whatever they want with their data and this behaviour has been removed.
* Bug fixed. `psi` incorrectly assigned site names.

2# biomonitoR 0.9

* Added a `NEWS.md` file to track changes to the package.

