# Visualize multiple dose realizations

Create a descriptive figure to visualize the distribution of dose and
its uncertainty.

## Usage

``` r
ecdfplot(data, dosevars, xlab="Dose", ylab="Cumulative distribution")
```

## Arguments

- data:

  data frame containing columns with dose vectors

- dosevars:

  names or column indices of dose vectors.

- xlab:

  label for the x-axis, default `"Dose"`

- ylab:

  label for the y-axis, default `"Cumulative distribution"`

## Details

In the left panel, the empirical cumulative distribution function (ECDF)
is plotted for each dose realization. In other words, each curve shows
one distribution of dose across individuals. The spread within
individual curves reflects the dose range across individuals, while the
spread between curves reflects between-realization variation on the
cohort level.

In the right panel, ECDFs are plotted for each individual, showing
distributions within individuals. A wide spread within individual curves
is indicative of large within-individual variation, while the spread
between curves reflects between-individual variation.

Any zeros in the doses are excluded while plotting.

## Examples

``` r
   data(data, package="ameras")
   ecdfplot(data,  dosevars=paste0("V", 1:10))
```
