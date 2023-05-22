## code to prepare `hf_action` dataset goes here

library(tidyverse)
library(WR)

data(hfaction_cpx9)
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

usethis::use_data(hf_action, overwrite = TRUE)
