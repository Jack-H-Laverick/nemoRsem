
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nemoRsem <a href='https://jack-h-laverick.github.io/nemoRsem'><img src='man/figures/logo.svg' align="right" height="139" /></a>

<!-- badges: start -->

[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R build
status](https://github.com/Jack-H-Laverick/nemoRsem/workflows/R-CMD-check/badge.svg)](https://github.com/Jack-H-Laverick/nemoRsem/actions)
[![Codecov test
coverage](https://codecov.io/gh/Jack-H-Laverick/nemoRsem/branch/main/graph/badge.svg)](https://codecov.io/gh/Jack-H-Laverick/nemoRsem?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/Jack-H-Laverick/nemoRsem/badge)](https://www.codefactor.io/repository/github/Jack-H-Laverick/nemoRsem)
[![](https://img.shields.io/github/last-commit/Jack-H-Laverick/nemoRsem.svg)](https://github.com/Jack-H-Laverick/nemoRsem/commits/main)
<!-- badges: end -->

nemoRsem is an R package for creating common spatial and temporal
summaries of NEMO-ERSEM model output. The functions contained in this
package streamline the process of sampling along transects, within
volumes, and across depths. The functions balance processing speed
against memory load so a standard laptop should be able to perform these
operations.

<br/>

Typically the amount of data in memory at any one time is small enough
to allow files to be processed in parallel. When tested at the
University of Strathclyde the limit on speed was accessing data files
from the server, instead of processing time of each file. The largest
effect you can have on speed as an end user is likely where you choose
to store the NEMO-ERSEM model output relative to the machine running R.

## Installation

The internal facing / development version can be downloaded from
[GitHub](https://github.com/) with:

``` r
remotes::install_github("Jack-H-Laverick/nemoRsem")
```

## Project Documentation

Head to
[GitHub.io](https://jack-h-laverick.github.io/nemoRsem/index.html) to
view the full documentation.

<br/>

The navbar contains links to an index of documented functions and a
change log.
