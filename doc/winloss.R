## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE-----------------------------------------------------
library(prodlim)
library(winloss)
library(tidyverse)
library(magrittr)

## -----------------------------------------------------------------------------
library(WR)
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

## -----------------------------------------------------------------------------
fit <- winloss(id=data$id,
               time=data$time, 
               status=data$status, 
               group=data$group, 
               type=c(1,2), 
               at=47)

## -----------------------------------------------------------------------------
c(fit$wr,fit$l_wr,fit$u_wr)

