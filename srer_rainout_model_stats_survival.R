# USDA Rainout - Santa Rita Experimental Range
# Stats and Modeling-Survival
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-01-25

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
library(ggpmisc)


### Read in seedlings data (cohorts 1-3), make all columns factor and date a date
seedlings <- vroom("Data/seedlings_combined.csv",
                   col_select = -c(1),
                   col_types = c(.default = "f",
                                 date = "D"))
str(seedlings)

# Fate counts
seedlings_fate_full <- seedlings %>% 
  group_by(block, precip, clip, excl, date, cohort) %>%
  count(fate) %>% 
  pivot_wider(names_from = fate,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_germ = "0",
         survival = "1",
         died = "2") %>% 
  mutate(tot_germination = survival + died,
         surv_perc = (survival/tot_germination)) %>% 
         mutate(surv_perc = replace_na(surv_perc, 0))

# Herbivory counts
seedlings_herb_full <- seedlings %>% 
  group_by(block, precip, clip, excl, date, cohort) %>%
  count(herbivory) %>% 
  pivot_wider(names_from = herbivory,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_herb = "0",
         herb_lived = "1",
         herb_died = "2") %>% 
  mutate(tot_herbivory = herb_lived + herb_died)

# Granivory counts
seedlings_gran_full <- seedlings %>% 
  group_by(block, precip, clip, excl, date, cohort) %>%
  count(granivory) %>% 
  pivot_wider(names_from = granivory,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_gran = "0",
         granivory = "1")

# Joining fate, herbivory, and granivory variables of interest to one table
seedlings_all_full <- seedlings_fate_full %>% 
  left_join(seedlings_herb_full) %>% 
  left_join(seedlings_gran_full) 
# %>% 
#   dplyr::select(-c(no_germ, no_herb, no_gran))

## percent germ (relative to lab conditions), lost to abiotic constraints, lost to preds
germ1 <- seedlings_all_full %>% 
  filter(cohort == 1) %>% 
  mutate(rel_germ = ((100*(tot_germination/10))*0.99)) %>% 
  mutate(abiotic = 10*(no_germ-granivory)) %>% 
  group_by(cohort, precip, clip, excl) %>% 
  summarise(mean_gran = 10*mean(granivory),
            sd_gran = 10*sd(granivory),
            counts = n(),
            se_gran = (sd_gran/sqrt(counts)),
            mean_ab = mean(abiotic),
            sd_ab = sd(abiotic),
            se_ab = (sd_ab/sqrt(counts)),
            mean_rel = mean(rel_germ),
            sd_rel = sd(rel_germ),
            se_rel = (sd_rel/sqrt(counts)))


germ2 <- seedlings_all_full %>% 
  filter(cohort == 2) %>% 
  mutate(rel_germ = ((100*(tot_germination/10))*0.99)) %>% 
  mutate(abiotic = 10*(no_germ-granivory)) %>% 
  group_by(cohort, precip, clip, excl) %>% 
  summarise(mean_gran = 10*mean(granivory),
            sd_gran = 10*sd(granivory),
            counts = n(),
            se_gran = (sd_gran/sqrt(counts)),
            mean_ab = mean(abiotic),
            sd_ab = sd(abiotic),
            se_ab = (sd_ab/sqrt(counts)),
            mean_rel = mean(rel_germ),
            sd_rel = sd(rel_germ),
            se_rel = (sd_rel/sqrt(counts)))

germ3 <- seedlings_all_full %>% 
  filter(cohort == 3) %>% 
  mutate(rel_germ = ((100*(tot_germination/10))*0.99)) %>% 
  mutate(abiotic = 10*(no_germ-granivory)) %>% 
  group_by(cohort, precip, clip, excl) %>% 
  summarise(mean_gran = 10*mean(granivory),
            sd_gran = 10*sd(granivory),
            counts = n(),
            se_gran = (sd_gran/sqrt(counts)),
            mean_ab = mean(abiotic),
            sd_ab = sd(abiotic),
            se_ab = (sd_ab/sqrt(counts)),
            mean_rel = mean(rel_germ),
            sd_rel = sd(rel_germ),
            se_rel = (sd_rel/sqrt(counts)))

germ_full <- rbind(germ1,germ2,germ3)

germ_clean <- germ_full %>% dplyr::select(-sd_gran, -counts, -sd_ab, -mean_rel, -sd_rel,-se_rel) %>%
  mutate(excl = recode_factor(excl, # coded based on access, control = all excl, total = no excl
                              "Control" = "Total",
                              "Rodents" = "Ants Excl",
                              "Ants" = "Rodents Excl",
                              "Total" = "None"))
  

# Histogram of variables, all zero-inflated poisson except for tot_germination
seedlings_all_full %>% 
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram(bins = 100, binwidth = 1)

glimpse(seedlings_all_full)
summary(seedlings_all_full)

# descriptive stats for each cohort 1-3
describeBy(seedlings_all_full, group = "cohort")

# make a new ID column for each observation and plot
seedlings_obs <- seedlings_all_full %>% 
  rowid_to_column("ObsID") %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("sampID", block:excl, sep = "_", remove = FALSE) %>%
  unite("plotID", block:precip, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         date = as.factor(date),
         ObsID = as.factor(ObsID),
         sampID = as.factor(sampID),
         plotID = as.factor(plotID))

# remove two weeks after sowing
seedlings_obs_year <- seedling_fate_year %>% 
  rowid_to_column("ObsID") %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("plotID",block:precip, sep = "_", remove = FALSE) %>%
  unite("sampID", block:excl, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         plotID = as.factor(plotID),
         date = as.factor(date),
         ObsID = as.factor(ObsID),
         sampID = as.factor(sampID))

# mixed effects model with nesting

hist(log(seedlings_obs_year$survival)) # still zero-inflated
hist(sqrt(seedlings_obs_year$survival)) # still zero-inflated

# count data...poisson and zero-inflated (41% of data are zeros)
num_obs <- seedlings_obs_year %>% ungroup() %>% summarise(obs = n())

zeros <- seedlings_obs_year %>%
  ungroup() %>% 
  filter(survival == "0") %>% 
  summarise(zero = 100*(n()/num_obs$obs))


#glmmTMB function to build zero-inflated model, poisson dist., random factors

# accounting for sample independence and each year and date
# are the main random factors to include in the models
# (e.g, (1|sampID/date))
# Thus, plot ID and date and the date-by-plot interaction are used in building final model

# test if survival sig diff across all blocks?
zi.block <- glmmTMB(survival ~ block + (1|sampID/date),
                    data = seedlings_obs_year,
                    ziformula = ~block,
                    family = poisson())
zi.block.sum <- summary(zi.block)
zi.block.sum

# not all blocks sig, so exclude from potential models

# start with full model, and make simpler
# precipitation, clipping, and exclusion total interactions as fixed factors
seed_test_1<-seedlings_obs_year %>% filter(year == 1)
seed_test_2<-seedlings_obs_year %>% filter(year == 2)
seed_test_3<-seedlings_obs_year %>% filter(year == 3)


zi.srer.surv.full <- glmmTMB(survival ~ precip + precip/clip+ precip/excl + excl/clip + excl + clip + precip/clip/excl + year +
                               precip/year + precip/clip/year + precip/excl/year + excl/clip/year + excl/year + clip/year +
                               precip/clip/excl/year + (1|sampID/date),
                             data = seedlings_obs_year,
                             ziformula = ~1,
                             family = poisson())
zi.srer.surv.full.sum <- summary(zi.srer.surv.full)
zi.srer.surv.full.sum

# model for each year of survival
zi.srer.surv.1 <- glmmTMB(survival ~ precip + precip/clip+ precip/excl + clip/excl + excl + clip + precip/clip/excl + cohort +
                               precip/cohort + precip/clip/cohort + precip/excl/cohort + excl/clip/cohort + excl/cohort + clip/cohort +
                               precip/clip/excl/cohort +(1|sampID) + ar1(date + 0|cohort),
                             data = seed_test_1,
                             ziformula = ~precip+excl+clip+cohort,
                             family = poisson())
zi.srer.surv.1.sum <- summary(zi.srer.surv.1)
zi.srer.surv.1.sum


zi.srer.surv.2 <- glmmTMB(survival ~ precip + precip/clip+ precip/excl + clip/excl + excl + clip + precip/clip/excl + cohort +
                            precip/cohort + precip/clip/cohort + precip/excl/cohort + excl/clip/cohort + excl/cohort + clip/cohort +
                            precip/clip/excl/cohort + (1|sampID) + ar1(date + 0|cohort),
                             data = seed_test_2,
                             ziformula = ~precip+excl+clip+cohort,
                             family = poisson())
zi.srer.surv.2.sum <- summary(zi.srer.surv.2)
zi.srer.surv.2.sum

zi.srer.surv.3 <- glmmTMB(survival ~ precip + precip/clip+ precip/excl +
                            clip/excl + excl + clip + precip/clip/excl +(1|sampID) + ar1(date + 0|cohort),
                             data = seed_test_3,
                             ziformula = ~precip+excl+clip,
                             family = poisson())
zi.srer.surv.3.sum <- summary(zi.srer.surv.3)
zi.srer.surv.3.sum


# type II wald's for fixed effect significance
Anova(zi.srer.surv.full)
Anova(zi.srer.surv.1)
Anova(zi.srer.surv.2)
Anova(zi.srer.surv.3)

post.hoc.full <- emmeans::emmeans(zi.srer.surv.full, specs = ~clip*year, type = "response")

# get lettering report on post-hoc test
post.hoc.letters.full <- cld(post.hoc.full, Letters = letters, covar = T, level = 0.0001)
post.hoc.letters.full

post.hoc.surv <- as.data.frame(post.hoc.full) %>% mutate(surv_per = 10*rate)

# use first year for ppt modeling
# precipitation, clipping, and exclusion fixed factors
# poisson model convergence issue with cohort and sample with temp autocorrelation, use neg binomial
zi.srer.surv.pce <- glmmTMB(survival ~ precip + clip + excl + (1|sampID) + ar1(date + 0|cohort),
                            data = seed_test_1,
                            ziformula = ~precip + clip + excl,
                            family = poisson())
zi.srer.surv.pce.sum <- summary(zi.srer.surv.pce)
zi.srer.surv.pce.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pce)

# clipping not sig, remove from model
# precipitation and exclusion fixed factors
zi.srer.surv.pe <- glmmTMB(survival ~ precip + excl + (1|sampID) + ar1(date + 0|cohort),
                           data = seed_test_1,
                           ziformula = ~precip+excl,
                           family = poisson())
zi.srer.surv.pe.sum <- summary(zi.srer.surv.pe)
zi.srer.surv.pe.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pe)

# precipitation only fixed factor
zi.srer.surv.p <- glmmTMB(survival ~ precip + (1|sampID) + ar1(date + 0|cohort),
                          data = seed_test_1,
                          ziformula = ~precip,
                          family = poisson())
zi.srer.surv.p.sum <- summary(zi.srer.surv.p)
zi.srer.surv.p.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.p)

# precipitation, exclusion, and precipitation and exclusion interaction fixed factors
zi.srer.surv.pe.int <- glmmTMB(survival ~ precip + precip/excl + (1|sampID) + ar1(date + 0|cohort),
                               data = seed_test_1,
                               ziformula = ~precip+excl,
                               family = poisson())
zi.srer.surv.pe.int.sum <- summary(zi.srer.surv.pe.int)
zi.srer.surv.pe.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pe.int)

# precipitation, exclusion, and precipitation and exclusion interaction fixed factors
zi.srer.surv.pe.int.1 <- glmmTMB(survival ~ precip/excl + (1|sampID) + ar1(date + 0|cohort),
                               data = seed_test_1,
                               ziformula = ~precip+excl,
                               family = poisson())
zi.srer.surv.pe.int.1.sum <- summary(zi.srer.surv.pe.int.1)
zi.srer.surv.pe.int.1.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pe.int.1)

# precipitation, clip, and precipitation and clip interaction fixed factors
zi.srer.surv.pc.int <- glmmTMB(survival ~ precip + precip/clip + (1|sampID) + ar1(date + 0|cohort),
                               data = seed_test_1,
                               ziformula = ~precip+clip,
                               family = poisson())
zi.srer.surv.pc.int.sum <- summary(zi.srer.surv.pc.int)
zi.srer.surv.pc.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pc.int)

# add year
zi.srer.surv.pcey <- glmmTMB(survival ~ precip + clip + excl + cohort + (1|sampID) + ar1(date + 0|cohort),
                             data = seed_test_1,
                             ziformula = ~precip + clip + excl + cohort,
                             family = poisson())
zi.srer.surv.pcey.sum <- summary(zi.srer.surv.pcey)
zi.srer.surv.pcey.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pcey)

# clipping not sig, remove from model
# precipitation and exclusion fixed factors
zi.srer.surv.pey <- glmmTMB(survival ~ precip + excl + cohort +  (1|sampID) + ar1(date + 0|cohort),
                            data = seed_test_1,
                            ziformula = ~precip + excl + cohort,
                            family = poisson())
zi.srer.surv.pey.sum <- summary(zi.srer.surv.pey)
zi.srer.surv.pey.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pey)

# precipitation only fixed factor
zi.srer.surv.py <- glmmTMB(survival ~ precip + cohort + (1|sampID) + ar1(date + 0|cohort),
                           data = seed_test_1,
                           ziformula = ~precip+cohort,
                           family = poisson())
zi.srer.surv.py.sum <- summary(zi.srer.surv.py)
zi.srer.surv.py.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.py)

# precipitation, exclusion, and precipitation and exclusion interaction fixed factors
zi.srer.surv.pe.int.yr <- glmmTMB(survival ~ precip + precip/excl + cohort + (1|sampID) + ar1(date + 0|cohort),
                                  data = seed_test_1,
                                  ziformula = ~precip+excl+cohort,
                                  family = poisson())
zi.srer.surv.pe.int.yr.sum <- summary(zi.srer.surv.pe.int.yr)
zi.srer.surv.pe.int.yr.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pe.int.yr)

# precipitation, clip, and precipitation and clip interaction fixed factors
zi.srer.surv.pc.int.yr <- glmmTMB(survival ~ precip + precip/clip + cohort + (1|sampID) + ar1(date + 0|cohort),
                                  data = seed_test_1,
                                  ziformula = ~precip+clip+cohort,
                                  family = poisson())
zi.srer.surv.pc.int.yr.sum <- summary(zi.srer.surv.pc.int.yr)
zi.srer.surv.pc.int.yr.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pc.int.yr)

# add year interaction
zi.srer.surv.py.int <- glmmTMB(survival ~ precip/cohort + (1|sampID) + ar1(date + 0|cohort),
                               data = seed_test_1,
                               ziformula = ~precip+cohort,
                               family = poisson())
zi.srer.surv.py.int.sum <- summary(zi.srer.surv.py.int)
zi.srer.surv.py.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.py)

# clipping not sig, remove from model
# precipitation and exclusion fixed factors
zi.srer.surv.ey.int <- glmmTMB(survival ~ excl/cohort +(1|sampID) + ar1(date + 0|cohort),
                               data = seed_test_1,
                               ziformula = ~excl+cohort,
                               family = poisson())
zi.srer.surv.ey.int.sum <- summary(zi.srer.surv.ey.int)
zi.srer.surv.ey.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.ey.int)

# precipitation only fixed factor
zi.srer.surv.cyr.int <- glmmTMB(survival ~ clip/cohort + (1|sampID) + ar1(date + 0|cohort),
                                data = seed_test_1,
                                ziformula = ~clip+cohort,
                                family = poisson())
zi.srer.surv.cyr.int.sum <- summary(zi.srer.surv.cyr.int)
zi.srer.surv.cyr.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.cyr.int)

# precipitation, exclusion, and precipitation and exclusion interaction fixed factors
zi.srer.surv.pe.int.yr.int <- glmmTMB(survival ~ precip/excl/cohort + (1|sampID) + ar1(date + 0|cohort),
                                      data = seed_test_1,
                                      ziformula = ~precip+excl+cohort,
                                      family = poisson())
zi.srer.surv.pe.int.yr.int.sum <- summary(zi.srer.surv.pe.int.yr.int)
zi.srer.surv.pe.int.yr.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pe.int.yr.int)

# precipitation, clip, and precipitation and clip interaction fixed factors
zi.srer.surv.pc.int.yr.int <- glmmTMB(survival ~ precip/clip/cohort + (1|sampID) + ar1(date + 0|cohort),
                                      data = seed_test_1,
                                      ziformula = ~precip+clip+cohort,
                                      family = poisson())
zi.srer.surv.pc.int.yr.int.sum <- summary(zi.srer.surv.pc.int.yr.int)
zi.srer.surv.pc.int.yr.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.surv.pc.int.yr.int)

# compare AIC scores of all potential models for model selection
aic.compare.final.surv <- AICtab(zi.srer.surv.p,
                                 zi.srer.surv.pe,
                                 zi.srer.surv.pce,
                                 #zi.srer.surv.1,
                                 zi.srer.surv.pe.int,
                                 zi.srer.surv.pe.int.1,
                                 zi.srer.surv.pc.int,
                                 zi.srer.surv.pcey,
                                 zi.srer.surv.pey,
                                 zi.srer.surv.py,
                                 zi.srer.surv.pe.int.yr,
                                 zi.srer.surv.pc.int.yr,
                                 zi.srer.surv.py.int,
                                 zi.srer.surv.ey.int,
                                 zi.srer.surv.cyr.int,
                                 zi.srer.surv.pe.int.yr.int,
                                 zi.srer.surv.pc.int.yr.int,
                                 logLik = TRUE)

aic.compare.final.surv

# does removing clipping sig improve model?
anova(zi.srer.surv.pce, zi.srer.surv.pe) # yes

# final PRVE seedlings survival model will precipitation and exclusion as only fixed factors
zi.srer.surv.final <- glmmTMB(survival ~ precip + excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                              data = seedlings_obs_year,
                              ziformula = ~.,
                              family = poisson())
zi.srer.surv.final.sum <- summary(zi.srer.surv.final)
zi.srer.surv.final.sum

# define new variable to be best model for diagnostics
best.model.surv <- zi.srer.surv.final

# extract only the variance correlation matrix
VarCorr(best.model.surv)

# examine model residuals
seedlings_obs$res_surv <- residuals(best.model.surv, quantileFunction = qpois)

# examine histogram of model residuals
hist(seedlings_obs$res_surv) # look normal

# extract model co-variance matrix
vcov.matrix.surv <- vcov(best.model.surv)

# extract conditional model coefficients
cond.cov.surv <- as.data.frame(vcov.matrix.surv$cond)

# extract zero inflation model coefficients (probably not needed)
zi.cov.surv <- as.data.frame(vcov.matrix.surv$zi)

# use Dharma package to examine model fit
seedlings.best.simres.surv <- simulateResiduals(best.model.surv, n = 1000, plot = TRUE, integerResponse = TRUE)

# save scaled residuals
seedlings_obs$sim_surv <- seedlings.best.simres.surv$scaledResiduals

# save model fitted residuals
seedlings_obs$sim_fit_surv <- seedlings.best.simres.surv$fittedResiduals

# look at simulated/model residuals
hist(seedlings_obs$sim_surv) # pretty level, slight inc at end but makes sense with slight over prediction

hist(seedlings_obs$sim_fit_surv) # normal residuals

# are there outliers?
testOutliers(seedlings.best.simres.surv, type = 'bootstrap') # no outliers

# is the model overly dispersed
testDispersion(seedlings.best.simres.surv) # not overly dispersed

# KS test
testUniformity(seedlings.best.simres.surv) # slight hump in middle, observed survival slightly higher than model expected

# model account for temp autocorrelation
testTemporalAutocorrelation(seedlings.best.simres.surv) # temporal autocorrelation is accounted for in model

# model account for zero inflation?
testZeroInflation(seedlings.best.simres.surv) # zero inflation accounted for in the model

# does the model mean prediction fit?
means <- function(x) mean(x) 
testGeneric(seedlings.best.simres.surv, summary = means) # yes, does not sig deviate (p = 0.77)

# does the model standard deviation prediction fit?
spread <- function(x) sd(x)
testGeneric(seedlings.best.simres.surv, summary = spread) # yes, does not sig deviate (p = 0.22)

# take a look at model residuals
plotResiduals(seedlings.best.simres.surv, seedlings_obs$precip, quantreg = T) # no strange residual patterns

# save model predicted responses as a new column
seedlings_obs$pred_surv <- predict(best.model.surv, type = "response")

# look histograms of model predicted survival and original survival data 
hist(seedlings_obs$pred_surv)

hist(seedlings_obs$survival)

# another method for post-hoc test of precip and exclusion 
post.hoc <- emmeans::emmeans(zi.srer.surv.pe.int.yr.int, specs = ~precip*excl*year, type = "response")
post.hoc

# get lettering report on post-hoc test
post.hoc.letters <- cld(post.hoc, Letters = letters, covar = T)
post.hoc.letters

post.hoc.surv <- as.data.frame(post.hoc) %>% mutate(surv_per = 10*response)

# in case random effect coeff and model coeff are needed
ranef(best.model.surv)

coef(best.model.surv)

# write seedlings_obs data frame to csv for use/graphing
write.csv(seedlings_obs, file = "Data/seedlings_obs.csv", row.names = FALSE)

# write model outputs to text file

# sink("Model_Outputs/model-outputs-survival.txt")
# 
# zi.srer.surv.p.sum
# zi.srer.surv.pe.sum
# zi.srer.surv.pce.sum
# zi.srer.surv.full.sum
# zi.srer.surv.pe.int.sum
# aic.compare.final
# zi.srer.surv.final
# Anova(zi.srer.surv.full)
# Anova(zi.srer.surv.pce)
# Anova(zi.srer.surv.pe)
# Anova(zi.srer.surv.p)
# Anova(zi.srer.surv.pe.int)
# Anova(zi.srer.surv.pc.int)
# post.hoc
# post.hoc.letters
# 
# sink()

# Visualization: treatment ppt values with ggeffects
precip_cont_surv_df_1 <- seedlings_obs_year %>% 
  group_by(cohort) %>%
  filter(cohort == "1") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "339",
                                     "IR" = "559",
                                     "RO" = "119")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

precip_cont_surv_df_2 <-seedlings_obs_year %>% 
  group_by(cohort) %>%
  filter(cohort == "2") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "536",
                                     "IR" = "884",
                                     "RO" = "188")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

precip_cont_surv_df_3 <-seedlings_obs_year %>% 
  group_by(cohort) %>%
  filter(cohort == "3") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "569",
                                     "IR" = "939",
                                     "RO" = "199")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

# combine into one dataframe
precip_cont_surv_df <- rbind(precip_cont_surv_df_1, precip_cont_surv_df_2, precip_cont_surv_df_3)

# make ppt a factor for model
precip_cont_surv_df_1y <- precip_cont_surv_df %>% 
  filter(year == 1) %>% 
  mutate(precip_cont = as.factor(precip_cont))

zi.srer.surv.cont <- glmmTMB(survival ~ precip_cont+excl+cohort+ (1|sampID) + ar1(date + 0|cohort),
                             data = precip_cont_surv_df_1y,
                             ziformula = ~precip_cont+excl+cohort,
                             family = poisson())
zi.srer.surv.cont.sum <- summary(zi.srer.surv.cont)
zi.srer.surv.cont.sum

# get predictions of model
mydf <- ggpredict(zi.srer.surv.cont, type = "simulate_random", terms = c("precip_cont", "excl", "cohort"))

# create graph
ggeff_pred_ppt_fig <- as.data.frame(mydf) %>%
  mutate(excl = group, #cohort = facet,
         x = as.numeric(as.character(x))) %>%
  mutate(excl = recode_factor(excl, # coded based on access
                              "Control" = "All Excl",
                              "Ants" = "Rodents Excl",
                              "Rodents" = "Ants Excl",
                              "Total" = "None")) %>%  
  ggplot(aes(x = x, y = predicted*10)) +
  geom_point(aes(color = excl)) +
  #geom_pointrange(aes(ymin = 10*conf.low, ymax = 10*conf.high, color = excl), size = 0.5) +
  geom_smooth(method = "glm", formula = y ~ log(x) + x, se = T, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Precipitation (mm)",
       color = "Exclusion") +
  scale_x_continuous(breaks = c(0,100,200,300,400,500,600,700,800,900,1000), limits = c(0, 1000))+
  ylim(0,110) +
  theme_pubr(legend = "right")+
  labs_pubr(base_size = 24)

ggeff_pred_ppt_fig

ggsave(filename = "Figures_Tables/pred_surv_cont.tiff",
       plot = ggeff_pred_ppt_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# by exclusion tx
ggeff_excl_ppt_fig <- as.data.frame(mydf) %>%
  mutate(excl = group, cohort = facet,
         x = as.numeric(as.character(x))) %>%
  mutate(excl = recode_factor(excl, 
                              "Control" = "All Excl",
                              "Ants" = "RodentsExcl",
                              "Rodents" = "Ants Excl",
                              "Total" = "None")) %>%  
  ggplot(aes(x = x, y = predicted*10, group = excl, color = excl)) +
  geom_smooth(method = "glm", formula = y ~ log(x) + x, se = F, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Precipitation (mm)",
       color = "Exclusion") +
  scale_x_continuous(breaks = c(0,100,200,300,400,500,600,700,800,900,1000), limits = c(0, 1000))+
  ylim(0, 110) +
  theme_pubr(legend = "right")+
  labs_pubr(base_size = 24)

ggeff_excl_ppt_fig

ggsave(filename = "Figures_Tables/pred_surv_excl_cont.tiff",
       plot = ggeff_excl_ppt_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

