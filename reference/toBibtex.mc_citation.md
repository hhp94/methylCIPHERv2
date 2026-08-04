# Bibtex Method For An mc_citation Object

Returns the bibtex text for the papers cited in an `mc_citation` object.

## Usage

``` r
# S3 method for class 'mc_citation'
toBibtex(object, ...)
```

## Arguments

- object:

  An `mc_citation` object. The value returned by
  [`cite_clocks()`](https://hhp94.github.io/methylCIPHERv2/reference/cite_clocks.md).

- ...:

  Not used.

## Value

A character vector of class `Bibtex`. One bibtex entry for each cited
paper, as lines of text.

## Examples

``` r
cites <- cite_clocks(c("Horvath1", "Hannum"))
toBibtex(cites)
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
```
