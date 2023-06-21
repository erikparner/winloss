
<!-- README.md is generated from README.Rmd. Please edit that file -->

# winloss

<!-- badges: start -->
<!-- badges: end -->

Inference for win-loss parameters using censored event data. The event
data can be single or recurrent events. Only the event type that is
prioritized last is allow to be recurrent, in which case the number of
recurrent event is compared.

## Installation

You can install the development version of winloss from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("erikparner/winloss")
```

## Example

The package uses the R package prodlim.

``` r
library(winloss)
library(prodlim)
## basic example code
fit <- winloss(id=hf_action$id,
               time=hf_action$time, 
               status=hf_action$status, 
               group=hf_action$group, 
               type=c(1,2), 
               at=47)
```

The fit contain a vector of the win-loss parameter *wl* ordered as:
first win, first loss, second win, second loss, ect. The asumptotic
variance *sigma*. For all win-loss parameters, stub say, the fit
contains the standard error, *se_stub*, lower and upper 95% confidence
interval, *l_stub* and *u_stub*. The coverage of the confidence interval
can be changed in the win-loss function. The parameters are

- win-loss of each event type: *wl*
- total win ratio and log win ratio: *wr* and *logwr*
- total win difference: *wd*
- event specific win ratio and log win ratio: *wrk* and *logwrk*
- event specific win difference: *wrk* and *logwrk*
- total win and loss: *w* and *l*
- event specific ranking probability: *rankedk*
- total ranking probability: *ranked*

For example, display the win-ratio with 95% confidence interval:

``` r
c(fit$wr,fit$l_wr,fit$u_wr)
#> [1] 1.3249359 0.8872591 1.9785146
```
