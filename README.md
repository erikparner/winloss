
<!-- README.md is generated from README.Rmd. Please edit that file -->

# winloss

<!-- badges: start -->
<!-- badges: end -->

Inference for win-loss parameters using censored event data. The event
data can in general be single event. However, we allow the last event to
be recurrent events, in which case the number of recurrent event is
compared.

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

# The HF-ACTION study

The HF-ACTION data comes from the WR package. We performed some initial
data cleaning:

``` r
library(WR)
library(tidyverse)
library(magrittr)

data <- hfaction_cpx9

# Rename and recode id.
data$patid <- as.numeric(substr(data$patid, 6, 12))
names(data) <- c("id","time","status","group","age60")

# Time=0 problem (id==1359)
data <- data %>%
  mutate(time=ifelse(time<0.00001, 0.01, time)) 

# time=lag(time) problem (id=662).
data <- data %>%
  group_by(id) %>%
  mutate(time=ifelse(time==lag(time, default=0),time+0.001,time))

head(data)
#> # A tibble: 6 × 5
#> # Groups:   id [2]
#>      id   time status group age60
#>   <dbl>  <dbl>  <int> <int> <int>
#> 1     1  7.25       2     0     1
#> 2     1 12.6        0     0     1
#> 3     2  0.754      2     0     1
#> 4     2  4.30       2     0     1
#> 5     2  4.75       2     0     1
#> 6     2 45.9        0     0     1

# hf_action <- data
```
