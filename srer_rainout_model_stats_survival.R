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
  mutate(tot_germination = survival + died)

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
  left_join(seedlings_gran_full) %>% 
  dplyr::select(-c(no_germ, no_herb, no_gran))

# Histogram of variables, all zero-inflated poisson except for tot_germination
seedlings_all_full %>% 
  keep(is.integer) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram(bins = 50, binwidth = 1)

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
  unite("plotID", block:clip, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         date = as.factor(date),
         ObsID = as.factor(ObsID),
         sampID = as.factor(sampID),
         plotID = as.factor(plotID))

# mixed effects model with nesting, sampID as random and date (cohort) within year for temporal autocorrelation

hist(log(seedlings_obs$survival)) # still zero-inflated
hist(sqrt(seedlings_obs$survival)) # still zero-inflated

# count data...poisson and zero-inflated (41% of data are zeros)
num_obs <- seedlings_obs %>% ungroup() %>% summarise(obs = n())

zeros <- seedlings_obs %>%
  ungroup() %>% 
  filter(survival == "0") %>% 
  summarise(zero = 100*(n()/num_obs$obs))


#glmmTMB function to build zero-inflated model, poisson dist., random factors


# examine possible random factors

# intercept only, starting with sample ID within each cohort as random factor
zi.srer.fix <- glmmTMB(survival ~ 1 + (1|cohort/sampID) + ar1(date + 0|cohort),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum <- summary(zi.srer.fix)
zi.srer.fix.sum

zi.srer.fix.0 <- glmmTMB(survival ~ 1 + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                       data = seedlings_obs,
                       ziformula = ~.,
                       family = poisson)
zi.srer.fix.sum.0 <- summary(zi.srer.fix.0)
zi.srer.fix.sum.0

# intercept only, sample ID to account for measurement independence
zi.srer.fix.1 <- glmmTMB(survival ~ 1 + (1|sampID) + ar1(date + 0|cohort),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.1 <- summary(zi.srer.fix.1)
zi.srer.fix.sum.1

# intercept only, each plot ID within each cohort as random factor
zi.srer.fix.2 <- glmmTMB(survival ~ 1 + (1|cohort/plotID) + ar1(date + 0|cohort),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.2 <- summary(zi.srer.fix.2)
zi.srer.fix.sum.2

# intercept only, only plot ID
zi.srer.fix.3 <- glmmTMB(survival ~ 1 + (1|plotID) + ar1(date + 0|cohort),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.3 <- summary(zi.srer.fix.3)
zi.srer.fix.sum.3

# calculate differences b/w AIC for a number of models
aic.compare.random <- AICtab(zi.srer.fix,
                             zi.srer.fix.0,
                             zi.srer.fix.1,
                             zi.srer.fix.2,
                             zi.srer.fix.3,
                      logLik = TRUE)

aic.compare.random

# Based on AIC and Log Likelihood, and model convergence/overdispersion issues,
# accounting for sample independence and each year plus temp 
# are the main random factors to include in the models
# (e.g, (1|cohort) + (1|sampID) + ar1(date + 0|cohort))
# Thus, cohort, sample ID and temp autocorrelation are used in building final model

# test if survival sig diff across all blocks?
zi.block <- glmmTMB(survival ~ block + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                    data = seedlings_obs,
                    ziformula = ~.,
                    family = poisson)
zi.block.sum <- summary(zi.block)
zi.block.sum

# not all blocks sig, so exclude from potential models

# start with full model, and make simpler
# precipitation, clipping, and exclusion total interactions as fixed factors
zi.srer.surv.full <- glmmTMB(survival ~ precip/clip/excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.surv.full.sum <- summary(zi.srer.surv.full)
zi.srer.surv.full.sum

# interactions not significant and not informative, use each tx individually
# precipitation, clipping, and exclusion fixed factors
# model convergence issue with cohort and sample with temp autocorrelation
# try neg binomial dist instead of poisson
zi.srer.surv.pce <- glmmTMB(survival ~ precip + clip + excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.surv.pce.sum <- summary(zi.srer.surv.pce)
zi.srer.surv.pce.sum

# clipping not sig, remove from model
# precipitation and exclusion fixed factors
zi.srer.surv.pe <- glmmTMB(survival ~ precip + excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                           data = seedlings_obs,
                           ziformula = ~.,
                           family = poisson)
zi.srer.surv.pe.sum <- summary(zi.srer.surv.pe)
zi.srer.surv.pe.sum

# precipitation only fixed factor
zi.srer.surv.p <- glmmTMB(survival ~ precip + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.surv.p.sum <- summary(zi.srer.surv.p)
zi.srer.surv.p.sum

# precipitation, exclusion, and precipitation and exclusion interaction fixed factors
zi.srer.surv.pe.int <- glmmTMB(survival ~ precip + precip/excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs,
                               ziformula = ~.,
                               family = poisson)
zi.srer.surv.pe.int.sum <- summary(zi.srer.surv.pe.int)
zi.srer.surv.pe.int.sum

# compare AIC scores of all potential models for model selection
aic.compare.final <- AICctab(zi.srer.surv.p,
                            zi.srer.surv.pe,
                            zi.srer.surv.pce,
                            zi.srer.surv.full,
                            zi.srer.surv.pe.int,
                            logLik = TRUE)

aic.compare.final

# does removing clipping sig improve model?
anova(zi.srer.surv.pce, zi.srer.surv.pe) # yes

# final PRVE seedlings survival model will precipitation and exclusion as only fixed factors
zi.srer.surv.final <- glmmTMB(survival ~ precip + excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                        data = seedlings_obs,
                        ziformula = ~.,
                        family = poisson)
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

# post-hoc test of within precip treatment differences
g.precip <- glht(best.model.surv, linfct = mcp(precip = "Tukey", interaction_average = T))
summary(g.precip) # drought treatments sig different from control and wet

# post-hoc test of within exclusion treatment differences
g.excl <- glht(best.model.surv, linfct = mcp(excl = "Tukey", interaction_average = T))
summary(g.excl) # total treatments sig different from control, ant and rodents not a big influence

# another method for post-hoc test of precip and exclusion 
post.hoc <- emmeans::emmeans(best.model.surv, specs = ~precip*excl)
post.hoc

# get lettering report on post-hoc test
post.hoc.letters <- cld(post.hoc, Letters = letters, covar = T)
post.hoc.letters

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
# post.hoc
# post.hoc.letters
# 
# sink()




