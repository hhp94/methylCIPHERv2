# Clock Citations

Builds the citations for a set of clocks, or for the clocks scored in an
`mc_result` object.

## Usage

``` r
cite_clocks(x, ...)

# S3 method for class 'character'
cite_clocks(x, ...)

# S3 method for class 'mc_result'
cite_clocks(x, ...)

# Default S3 method
cite_clocks(x, ...)
```

## Arguments

- x:

  A character vector. The clock ids or group ids to cite, or an
  `mc_result` object from
  [`calc_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/calc_clocks.md).

- ...:

  Not used.

## Value

An `mc_citation` object. It holds the clock-to-paper links and the
bibtex text for each paper.

## Details

A character vector cites the clocks it names, and a group id cites every
clock in the group. An `mc_result` object cites the clocks in its
scores. Any other class raises an error that names the two accepted
types.

## Examples

``` r
cite_clocks(c("Horvath1", "Hannum"))
#> <mc_citation> 2 clock(s) x 2 paper(s)
#> 
#> $bibtex [2 paper(s)]
#> @article{Horvath_2013_24138928,
#>   title = {{DNA} methylation age of human tissues and cell types},
#>   author = {Horvath, Steve},
#>   year = {2013},
#>   journal = {Genome biology},
#>   volume = {14},
#>   number = {10},
#>   pages = {R115},
#>   doi = {10.1186/gb-2013-14-10-r115},
#>   pmid = {24138928},
#>   url = {https://doi.org/10.1186/gb-2013-14-10-r115}
#> }
#> 
#> @article{Hannum_2013_23177740,
#>   title = {Genome-wide methylation profiles reveal quantitative views of human aging rates},
#>   author = {Hannum, Gregory and Guinney, Justin and Zhao, Ling and Zhang, Li and Hughes, Guy and Sadda, SriniVas and Klotzle, Brandy and Bibikova, Marina and Fan, Jian-Bing and Gao, Yuan and Deconde, Rob and Chen, Menzies and Rajapakse, Indika and Friend, Stephen and Ideker, Trey and Zhang, Kang},
#>   year = {2013},
#>   journal = {Molecular cell},
#>   volume = {49},
#>   number = {2},
#>   pages = {359--367},
#>   doi = {10.1016/j.molcel.2012.10.016},
#>   pmid = {23177740},
#>   url = {https://doi.org/10.1016/j.molcel.2012.10.016}
#> }
#> 
#> ℹ `as.data.frame(x)` returns the clock-to-paper table.
#> ℹ `writeLines(toBibtex(x), "refs.bib")` writes the bibtex to a file.
#> ℹ `citation("methylCIPHERv2")` cites the package itself.

clocks <- c("Horvath1", "Hannum")
sim <- sim_DNAm(clocks, n = 5)
res <- calc_clocks(sim[["DNAm"]], clocks)
cite_clocks(res)
#> <mc_citation> 2 clock(s) x 2 paper(s)
#> 
#> $bibtex [2 paper(s)]
#> @article{Horvath_2013_24138928,
#>   title = {{DNA} methylation age of human tissues and cell types},
#>   author = {Horvath, Steve},
#>   year = {2013},
#>   journal = {Genome biology},
#>   volume = {14},
#>   number = {10},
#>   pages = {R115},
#>   doi = {10.1186/gb-2013-14-10-r115},
#>   pmid = {24138928},
#>   url = {https://doi.org/10.1186/gb-2013-14-10-r115}
#> }
#> 
#> @article{Hannum_2013_23177740,
#>   title = {Genome-wide methylation profiles reveal quantitative views of human aging rates},
#>   author = {Hannum, Gregory and Guinney, Justin and Zhao, Ling and Zhang, Li and Hughes, Guy and Sadda, SriniVas and Klotzle, Brandy and Bibikova, Marina and Fan, Jian-Bing and Gao, Yuan and Deconde, Rob and Chen, Menzies and Rajapakse, Indika and Friend, Stephen and Ideker, Trey and Zhang, Kang},
#>   year = {2013},
#>   journal = {Molecular cell},
#>   volume = {49},
#>   number = {2},
#>   pages = {359--367},
#>   doi = {10.1016/j.molcel.2012.10.016},
#>   pmid = {23177740},
#>   url = {https://doi.org/10.1016/j.molcel.2012.10.016}
#> }
#> 
#> ℹ `as.data.frame(x)` returns the clock-to-paper table.
#> ℹ `writeLines(toBibtex(x), "refs.bib")` writes the bibtex to a file.
#> ℹ `citation("methylCIPHERv2")` cites the package itself.
```
