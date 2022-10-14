
<!-- badges:start -->

[![R-CMD-check](https://github.com/Rekyt/fddimensionality_ms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Rekyt/fddimensionality_ms/actions/workflows/R-CMD-check.yaml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh///master?urlpath=rstudio)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5570256.svg)](https://doi.org/10.5281/zenodo.5570256)
<!-- badges:end -->

# fd-dimensionality: Effect of the number of traits on fd-env relationships

This is the companion repository for the following paper:

> Grenié M., Munoz F., Violle C. (2022). *When more is less: adding more
> traits dilutes the functional diversity-environment relationship*.
> *Submitted*

### How to cite

Please cite this compendium as:

> Grenié M., Munoz F., Violle C., (2022). *Compendium of R code and data
> for When more is less: adding more traits dilutes the functional
> diversity-environment relationship*. Accessed 14 Oct 2022. Online at
> <https://doi.org/10.5281/zenodo.5570256>

### How to download or install

You can install this compendium as an R package, `fddimensionality`,
from GitHub with:

``` r
# install.packages("devtools")
remotes::install_github("Rekyt/fddimensionality_ms")
```

### How to run

The analyses of the compendium heavily rely on the use the of workflow
management package [`drake`](https://cran.r-project.org/package=drake).
All analyses are contained within the `R` directory.

To run them run the following lines:

``` r
# Beware this may take a long time (>10 hours on a single computer)
pkgload::load_all()
drake::make(global_workflow())
```
