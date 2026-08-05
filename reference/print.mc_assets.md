# Print Method For An mc_assets Object

Prints the clock count and the CpG count of every group in an
`mc_assets` object.

## Usage

``` r
# S3 method for class 'mc_assets'
print(x, ...)
```

## Arguments

- x:

  An `mc_assets` object. The value returned by
  [`load_mc_assets()`](https://hhp94.github.io/methylCIPHERv2/reference/load_mc_assets.md).

- ...:

  Not used.

## Value

An `mc_assets` object. Returns `x`, invisibly, after printing it.

## Examples

``` r
if (FALSE) { # interactive()
assets <- load_mc_assets("PCBrainAge")
print(assets)
}
```
