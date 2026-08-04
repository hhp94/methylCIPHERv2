# Print Method For An mc_citation Object

Prints the clock and paper counts, then the bibtex text, for an
`mc_citation` object.

## Usage

``` r
# S3 method for class 'mc_citation'
print(x, ...)
```

## Arguments

- x:

  An `mc_citation` object. The value returned by
  [`cite_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/cite_clocks.md).

- ...:

  Not used.

## Value

An `mc_citation` object. Returns `x`, invisibly, after printing it.

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
```
