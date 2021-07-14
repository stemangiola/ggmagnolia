ggmagnolia
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02472/status.svg)](https://doi.org/10.21105/joss.02472)
<!-- badges: end -->

## Citation

\[to come\]

Please have a look also to

-   [tidyHeatmap](https://github.com/stemangiola/tidyHeatmap/) for tidy
    high-level data analysis and manipulation
-   [tidygate](https://github.com/stemangiola/tidygate/) for adding
    custom gate information to your tibble
-   [tidySingleCellExperiment](https://stemangiola.github.io/tidySingleCellExperiment/)
    for tidy manipulation of Seurat objects
-   [tidyseurat](https://stemangiola.github.io/tidyseurat/) for tidy
    manipulation of Seurat objects
-   [tidybulk](https://stemangiola.github.io/tidybulk/) for tidy
    high-level data analysis and manipulation
-   [tidySummarizedExperiment](https://stemangiola.github.io/tidySummarizedExperiment/)
    for heatmaps produced with tidy principles

website:
[stemangiola.github.io/ggmagnolia](https://stemangiola.github.io/ggmagnolia/)

`ggmagnolia` is a package that introduces tidy principles to the
creation of information-rich heatmaps. This package uses
[ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
as graphical engine.

## Functions/utilities available

| Function   | Description                        |
|------------|------------------------------------|
| `heatmap`  | Plot base heatmap                  |
| `add_tile` | Add tile annotation to the heatmap |

## Installation

To install the most up-to-date version

``` r
devtools::install_github("stemangiola/ggmagnolia")
```

To install the most stable version (however please keep in mind that
this package is under a maturing lifecycle stage)

``` r
install.packages("ggmagnolia")
```

## Contribution

If you want to contribute to the software, report issues or problems
with the software or seek support please open an issue
[here](https://github.com/stemangiola/ggmagnolia/issues)

## Input data frame

The heatmaps visualise a multi-element, multi-feature dataset, annotated
with independent variables. Each observation is a element-feature pair
(e.g., person-physical characteristics).

| element         | feature         | value     | independent\_variables |
|-----------------|-----------------|-----------|------------------------|
| `chr` or `fctr` | `chr` or `fctr` | `numeric` | …                      |

Let’s transform the mtcars dataset into a tidy
“element-feature-independent variables” data frame. Where the
independent variables in this case are ‘hp’ and ‘vs’.

``` r
 magnolia_input |> 
        plot_polar()
```

    ## Warning: Ignoring unknown aesthetics: width

    ## Warning: Ignoring unknown aesthetics: width

    ## Warning: Ignoring unknown aesthetics: width

    ## Warning: Ignoring unknown aesthetics: width

![](man/figures/unnamed-chunk-5-1.png)<!-- -->
