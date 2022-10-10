
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phyf

<!-- badges: start -->

[![R-CMD-check](https://github.com/rdinnager/phyf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rdinnager/phyf/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of phyf is to implement a `tibble` subclass useful for
statistical modelling on phylogenetic trees. It mainly implements an
phylogenetic flow (`pf`) object that is essentially a tibble with one or
more phylogenetic flow collection (`pfc`) columns. Phylogenetic flow
collection columns are collections of phylogenetic flow paths (`pfp`),
which store data on how information flows through a phylogeny from its
root node (phylogenies must be rooted to work with `phyf`) to its tips
(and its internal nodes). This allows for easy manipulation of the
phylogeny and associated data. The objects are used in the package
`fibre` for phylogenetic branch regression models, a highly felxible
framework for comparative analysis and modelling trait evolution across
a phylogeny.

## Installation

You can install the development version of phyf from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rdinnager/phyf")
```

## Example

..Nothing here yet..
