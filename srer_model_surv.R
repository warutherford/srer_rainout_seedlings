# USDA Rainout - Santa Rita Experimental Range
# Modeling PRVE germination and survival with long-term precipitation
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2022-01-24

# Load packages
library(tidyverse)
library(ggpubr)
library(vroom)
library(psych)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(multcomp)
library(MuMIn)
library(emmeans)
library(car)
library(ggeffects)

### Santa Rita Desert Grassland Station Rain Can/DESGR
ppt_srer <- vroom("Data/site-env-data/precip_srer_DESGR.csv",
                  col_types = c(.default = "n",
                                STATION = "c"))

glimpse(ppt_srer)

ppt_srer_avg_all <- ppt_srer %>%
  group_by(YEAR) %>% 
  rowwise() %>% 
  dplyr::summarise(ppt_sum = sum(c_across(JAN:DEC), na.rm = TRUE)) %>% 
  mutate(ppt_yr_in = ppt_sum/100) %>% 
  mutate(ppt_yr_mm = ppt_yr_in*25.4)

summary(ppt_srer_avg_all)

ppt_srer_avg_monsoon <- ppt_srer %>%
  group_by(YEAR) %>% 
  rowwise() %>% 
  dplyr::summarise(ppt_sum = sum(c_across(JUN:SEP), na.rm = TRUE)) %>% 
  mutate(ppt_yr_in = ppt_sum/100) %>% 
  mutate(ppt_yr_mm = ppt_yr_in*25.4)

summary(ppt_srer_avg_monsoon)

# arrange annual ppt data by year

ppt_arrange <- ppt_srer_avg_all %>% 
  arrange(ppt_yr_mm)

# pull out 10th and 90th percentile of annual ppt
quantile(ppt_arrange$ppt_yr_mm, c(.10, .90))

ppt_drought <- ppt_arrange %>% filter(ppt_yr_mm <= 290.3728) %>% 
  mutate(precip = "RO")

ppt_wet <- ppt_arrange %>% filter(ppt_yr_mm >= 534.5684) %>% 
  mutate(precip = "IR")

ppt_ambient <- ppt_arrange %>% filter(ppt_yr_mm > 290 & ppt_yr_mm < 534) %>% 
  mutate(precip = "CO")

ppt_srer_new <- rbind(ppt_ambient, ppt_drought, ppt_wet) %>% mutate(precip_cont = (precip), precip_cont = ppt_yr_mm)

# create a new data frame for prediction
nd_ann <- ppt_srer_new %>%
  mutate(cohort = as.factor(YEAR), sampID = as.factor(YEAR)) %>%
  dplyr::select(precip_cont, precip, cohort, sampID) %>% 
  arrange(YEAR)

nd_comp_ann <- rep(nd_ann, 4) %>% group_by(YEAR) %>% 
  mutate(excl = factor(c("Control", "Ants", "Rodents", "Total")),
         precip = as.factor(precip))

# use original survival model for annual survival prediction
zi.srer.surv.cont.pred <- glmmTMB(survival ~ precip_cont  + excl + (1|sampID),
                             data = precip_cont_surv_df,
                             ziformula = ~.,
                             family = poisson(link = 'log'))

summary(zi.srer.surv.cont.pred)

# create prediction data set
ppt_pred <- predict(zi.srer.surv.cont.pred, newdata = nd_comp_ann, type = "response", allow.new.levels = T)

ppt_pred <- as.data.frame(ppt_pred)

# combine prediction data with new data/years
srer_comb_ppt_surv <- cbind(ppt_pred, nd_comp_ann)

srer_mod_surv <- srer_comb_ppt_surv %>% 
  mutate(surv = 10*ppt_pred,
         year = YEAR,
         precip = precip_cont) %>% 
  dplyr::select(surv, year,precip, excl)

srer_mod_surv %>% filter(precip <= 290.3728) %>% 
  mutate(pptx = "Drought") %>% 
  group_by(pptx, excl) %>% 
  summarise(surv_mean = mean(surv),
            sd_surv = sd(surv),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts)))

srer_mod_surv %>% filter(precip >= 534.5684) %>% 
  mutate(pptx = "Wet") %>% 
  group_by(pptx, excl) %>% 
  summarise(surv_mean = mean(surv),
            sd_surv = sd(surv),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts)))

srer_mod_surv %>% filter(precip > 290 & precip < 534) %>% 
  mutate(pptx = "Ambient")%>% 
  group_by(pptx, excl) %>% 
  summarise(surv_mean = mean(surv),
            sd_surv = sd(surv),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts)))

srer_mod_surv %>% 
  ggplot(aes(x = precip, y = surv, alpha = year)) +
  geom_point()+
  ylim(0,30)

# Monsoon season only
ppt_arrange_mon <- ppt_srer_avg_monsoon %>% 
  arrange(ppt_yr_mm)

quantile(ppt_arrange_mon$ppt_yr_mm, c(.10, .90))

ppt_drought_mon <- ppt_arrange_mon %>% filter(ppt_yr_mm <= 151.0284) %>% 
  mutate(precip = "RO")

ppt_wet_mon <- ppt_arrange_mon %>% filter(ppt_yr_mm >= 332.1812) %>% 
  mutate(precip = "IR")

ppt_ambient_mon <- ppt_arrange_mon %>% filter(ppt_yr_mm > 151.0284 & ppt_yr_mm < 332.1812) %>% 
  mutate(precip = "CO")

ppt_srer_new_mon <- rbind(ppt_ambient_mon, ppt_drought_mon, ppt_wet_mon) %>% mutate(precip_cont = ppt_yr_mm)

nd_mon <- ppt_srer_new_mon %>%
  mutate(cohort = as.factor(YEAR), sampID = as.factor(YEAR)) %>%
  dplyr::select(precip_cont,precip, cohort, sampID) %>% 
  arrange(YEAR)

nd_comp <- rep(nd_mon, 4) %>% group_by(YEAR) %>% 
  mutate(excl = factor(c("Control", "Ants", "Rodents", "Total")),
         precip = as.factor(precip))

# germination model
zi.srer.germ.cont.pred <- glmmTMB(tot_germination ~ as.integer(as.character(precip_cont)) + excl + (1|cohort) + (1|sampID),
                                  data = precip_cont_germ_df,
                                  family = binomial(link = "logit"))
zi.srer.germ.cont.pred.sum <- summary(zi.srer.germ.cont.pred)
zi.srer.germ.cont.pred.sum

ppt_pred_mon <- predict(zi.srer.germ.cont.pred, newdata = nd_comp, type = "response", allow.new.levels = T, re.form =NA)

ppt_pred_mon <- as.data.frame(ppt_pred_mon)

srer_comb_ppt_germ <- cbind(ppt_pred_mon, nd_comp)

srer_mod_germ <- srer_comb_ppt_germ %>% 
  mutate(germ = 100*ppt_pred_mon,
         year = YEAR,
         precip = precip_cont) %>% 
  dplyr::select(germ, year,precip, excl)

srer_mod_germ %>% filter(precip <= 151.0284) %>% 
  mutate(pptx = "Drought") %>% 
  group_by(pptx, excl) %>% 
  summarise(germ_mean = mean(germ),
            sd_germ = sd(germ),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts)))

srer_mod_germ %>% filter(precip >= 332.1812) %>% 
  mutate(pptx = "Wet") %>% 
  group_by(pptx, excl) %>% 
  summarise(germ_mean = mean(germ),
            sd_germ = sd(germ),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts)))

srer_mod_germ %>% filter(precip > 151 & precip < 332) %>% 
  mutate(pptx = "Ambient")%>% 
  group_by(pptx, excl) %>% 
  summarise(germ_mean = mean(germ),
            sd_germ = sd(germ),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts)))

srer_mod_germ %>% 
  ggplot(aes(x = precip, y = germ, alpha = excl)) +
  geom_point()+
  ylim(0,30)



