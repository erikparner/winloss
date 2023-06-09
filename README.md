
<!-- README.md is generated from README.Rmd. Please edit that file -->

# winloss

<!-- badges: start -->
<!-- badges: end -->

Inference for win-loss parameters using censored event data. The
implementation focus on two prioritized events, where the first is death
and the second is either a single event or a recurrent event. The
methods are described in the manuscripts:

- Parner and Overgaard (2023). Estimation of win, loss probabilities and
  win ratio based on right-censored event data. Submitted.
- Parner and Overgaard (2023). Win-loss parameters for right-censored
  event data, with application to recurrent events. Submitted.

The method seem to generalize to multiple single events, Further, the
events may be a combination of single events and recurrent events as
long as the recurrent events are the last prioritized event.

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

For example, displaying the win-ratio with 95% confidence interval:

``` r
c(fit$wr,fit$l_wr,fit$u_wr)
#> [1] 1.3249359 0.8872591 1.9785146
```

## The HF-ACTION study

The HF-ACTION data comes from the WR package. We performed some initial
data cleaning for the version appearing in this package:

``` r
library(WR)
library(tidyverse)
library(magrittr)

hf_action <- hfaction_cpx9

# Rename and recode id.
hf_action$patid <- as.numeric(substr(hf_action$patid, 6, 12))
names(hf_action) <- c("id","time","status","group","age60")

# Time=0 problem (id==1359)
hf_action <- hf_action %>%
  mutate(time=ifelse(time<0.00001, 0.01, time)) 

# time=lag(time) problem (id=662).
hf_action <- hf_action %>%
  group_by(id) %>%
  mutate(time=ifelse(time==lag(time, default=0),time+0.001,time))
```
