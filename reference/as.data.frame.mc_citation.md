# Data Frame Method For An mc_citation Object

Converts the
[`cite_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/cite_clocks.md)
output to a data.frame.

## Usage

``` r
# S3 method for class 'mc_citation'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `mc_citation` object. The value returned by
  [`cite_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/cite_clocks.md).

- row.names:

  A character vector. Not used by this method. Default is `NULL`.

- optional:

  A boolean. Not used by this method. Default is `FALSE`.

- ...:

  Not used.

## Value

A data.frame. One row for each clock, with the bib key of the paper it
cites and publication details such as the title, the authors, and the
DOI.

## Examples

``` r
cites <- cite_clocks(c("Horvath1", "Hannum"))
as.data.frame(cites)
#>   clock_id     pmid    role               bib_key
#> 1 Horvath1 24138928 primary Horvath_2013_24138928
#> 2   Hannum 23177740 primary  Hannum_2013_23177740
#>                                                                             title
#> 1                             DNA methylation age of human tissues and cell types
#> 2 Genome-wide methylation profiles reveal quantitative views of human aging rates
#>                                                                                                                                                                                                                                                                                         author
#> 1                                                                                                                                                                                                                                                                               Horvath, Steve
#> 2 Hannum, Gregory and Guinney, Justin and Zhao, Ling and Zhang, Li and Hughes, Guy and Sadda, SriniVas and Klotzle, Brandy and Bibikova, Marina and Fan, Jian-Bing and Gao, Yuan and Deconde, Rob and Chen, Menzies and Rajapakse, Indika and Friend, Stephen and Ideker, Trey and Zhang, Kang
#>   year        journal volume number    pages                          doi
#> 1 2013 Genome biology     14     10     R115   10.1186/gb-2013-14-10-r115
#> 2 2013 Molecular cell     49      2 359--367 10.1016/j.molcel.2012.10.016
#>                                            url
#> 1   https://doi.org/10.1186/gb-2013-14-10-r115
#> 2 https://doi.org/10.1016/j.molcel.2012.10.016
```
