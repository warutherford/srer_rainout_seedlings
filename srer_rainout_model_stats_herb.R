# USDA Rainout - Santa Rita Experimental Range
# Stats and Modeling-herbivory
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-01-26

# Load packages
library(tidyverse)
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
  rowid_to_column("Obs_ID") %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("sampID", block:excl, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         date = as.factor(date),
         Obs_ID = as.factor(Obs_ID),
         sampID = as.factor(sampID))

#glmmTMB function to build zero-inflated model, poisson dist., random factors
# intercept only, starting with Observation ID within each cohort as random factor
zi.srer.fix <- glmmTMB(tot_herbivory ~ 1 + (1|cohort/Obs_ID),
                       data = seedlings_obs,
                       ziformula = ~.,
                       family = poisson)
zi.srer.fix.sum <- summary(zi.srer.fix)
zi.srer.fix.sum

# intercept only, sample ID to account for measurement independence
zi.srer.fix.1 <- glmmTMB(tot_herbivory ~ 1 + (1|sampID),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.1 <- summary(zi.srer.fix.1)
zi.srer.fix.sum.1

# intercept only, seedling cohort/year
zi.srer.fix.2 <- glmmTMB(tot_herbivory ~ 1 + (1|cohort),
                         data = seedlings_obs,
                         ziformula = ~1,
                         family = poisson)
zi.srer.fix.sum.2 <- summary(zi.srer.fix.2)
zi.srer.fix.sum.2

# intercept only, sample ID and cohort
zi.srer.fix.3 <- glmmTMB(tot_herbivory ~ 1 + (1|sampID) + (1|cohort),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.3 <- summary(zi.srer.fix.3)
zi.srer.fix.sum.3

# intercept only, sample ID within each cohort/year (cohort/sample ID interaction)
zi.srer.fix.4 <- glmmTMB(tot_herbivory ~ 1 + (1|cohort/sampID),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.4 <- summary(zi.srer.fix.4)
zi.srer.fix.sum.4

# intercept only, measurement date within each cohort/year
zi.srer.fix.5 <- glmmTMB(tot_herbivory ~ 1 + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.5 <- summary(zi.srer.fix.5)
zi.srer.fix.sum.5

# intercept only, sample ID and measurement date within each cohort
zi.srer.fix.6 <- glmmTMB(tot_herbivory ~ 1 + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.6 <- summary(zi.srer.fix.6)
zi.srer.fix.sum.6

# intercept only, temporal autocorrelation for measurement date based on each cohort
zi.srer.fix.auto <- glmmTMB(tot_herbivory ~ 1 + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.fix.sum.auto <- summary(zi.srer.fix.auto)
zi.srer.fix.sum.auto

# intercept only, sample ID and temporal autocorrelation
zi.srer.fix.auto1 <- glmmTMB(tot_herbivory ~ 1 + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.fix.sum.auto1 <- summary(zi.srer.fix.auto1)
zi.srer.fix.sum.auto1

# intercept only, sample/measurement for split based on cohort and temporal autocorrelation
# when adding fixed factors, unable to get model to converge and fill covariance matrix
zi.srer.fix.auto2 <- glmmTMB(tot_herbivory ~ 1 + (1|cohort/sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.fix.sum.auto2 <- summary(zi.srer.fix.auto2)
zi.srer.fix.sum.auto2

# calculate differences b/w AIC for a number of models
aic.compare.random <- AICtab(zi.srer.fix,
                             zi.srer.fix.1,
                             zi.srer.fix.2,
                             zi.srer.fix.3,
                             zi.srer.fix.4,
                             zi.srer.fix.5,
                             zi.srer.fix.6,
                             zi.srer.fix.auto,
                             zi.srer.fix.auto1,
                             zi.srer.fix.auto2,
                             logLik = TRUE)

aic.compare.random

# Based on AIC and Log Likelihood, sample ID and temp autocorrelation random factors accounts for most variation
# Thus, sample ID and temp autocorrelation is used in building final model
# See script "srer_rainout_model_stats_survival.R" for code to justify random factors

# look at distribution of germination counts
hist(seedlings_obs$tot_herbivory)

# test if survival sig diff across all blocks?
zi.block <- glmmTMB(tot_herbivory ~ block + (1|sampID) + ar1(date + 0|cohort),
                    data = seedlings_obs,
                    ziformula = ~.,
                    family = poisson)
zi.block.sum <- summary(zi.block)
zi.block.sum

#not all blocks sig, so exclude from potential models

# individual models with treatments using same random effects
# precipitation only
zi.srer.herb.p <- glmmTMB(tot_herbivory ~ precip + (1|sampID) + ar1(date + 0|cohort),
                          data = seedlings_obs,
                          ziformula = ~.,
                          family = poisson)
zi.srer.herb.p.sum <- summary(zi.srer.herb.p)
zi.srer.herb.p.sum

# precipitation and clipping only
zi.srer.herb.pc <- glmmTMB(tot_herbivory ~ precip + clip + (1|sampID) + ar1(date + 0|cohort),
                           data = seedlings_obs,
                           ziformula = ~.,
                           family = poisson)
zi.srer.herb.pc.sum <- summary(zi.srer.herb.pc)
zi.srer.herb.pc.sum

# precipitation and exclusion only
zi.srer.herb.pe <- glmmTMB(tot_herbivory ~ precip + excl + (1|sampID) + ar1(date + 0|cohort),
                           data = seedlings_obs,
                           ziformula = ~.,
                           family = poisson)
zi.srer.herb.pe.sum <- summary(zi.srer.herb.pe)
zi.srer.herb.pe.sum

# precipitation, clipping, and exclusion
zi.srer.herb.pce <- glmmTMB(tot_herbivory ~ precip + clip + excl + (1|sampID) + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.herb.pce.sum <- summary(zi.srer.herb.pce)
zi.srer.herb.pce.sum

# all interactions b/w precipitation, clipping, and exclusion 
zi.srer.herb.full <- glmmTMB(tot_herbivory ~ precip/clip/excl + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.herb.full.sum <- summary(zi.srer.herb.full)
zi.srer.herb.full.sum

# precipitation, clipping and precipitation interaction, and exclusion
zi.srer.herb.pc.int <- glmmTMB(tot_herbivory ~ precip + precip/clip + excl + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs,
                               ziformula = ~.,
                               family = poisson)
zi.srer.herb.pc.int.sum <- summary(zi.srer.herb.pc.int)
zi.srer.herb.pc.int.sum

# precipitation, clipping, and precipitation and exclusion interaction
zi.srer.herb.pce.int <- glmmTMB(tot_herbivory ~ precip + clip + precip/excl + (1|sampID) + ar1(date + 0|cohort),
                                data = seedlings_obs,
                                ziformula = ~.,
                                family = poisson)
zi.srer.herb.pce.int.sum <- summary(zi.srer.herb.pce.int)
zi.srer.herb.pce.int.sum

# precipitation, clipping, and clipping and exclusion interaction
zi.srer.herb.ce.int <- glmmTMB(tot_herbivory ~ precip + clip + clip/excl + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs,
                               ziformula = ~.,
                               family = poisson)
zi.srer.herb.ce.int.sum <- summary(zi.srer.herb.ce.int)
zi.srer.herb.ce.int.sum

# precipitation, exclusion, and precipiration and exlclusion interaction
zi.srer.herb.pe.int <- glmmTMB(tot_herbivory ~ precip + excl + precip/excl + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs,
                               ziformula = ~.,
                               family = poisson)
zi.srer.herb.pe.int.sum <- summary(zi.srer.herb.pe.int)
zi.srer.herb.pe.int.sum

# compare all AIC scores for model selection
aic.compare.final <- AICtab(zi.srer.herb.p,
                            zi.srer.herb.pc,
                            zi.srer.herb.pe,
                            zi.srer.herb.pce,
                            zi.srer.herb.full,
                            zi.srer.herb.pc.int,
                            zi.srer.herb.pce.int,
                            zi.srer.herb.pe.int,
                            zi.srer.herb.ce.int,
                            logLik = TRUE)

aic.compare.final

# best model based on AIC score
zi.srer.herb.aic <- glmmTMB(tot_herbivory ~ precip + excl + (1|sampID) + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.herb.aic.sum <- summary(zi.srer.herb.aic)
zi.srer.herb.aic.sum

# define new variable to be best model for diagnostics
best.model.herb <- zi.srer.herb.aic

### below is MuMIn model selection procedure with dredge function

zi.srer.herb.all <- glmmTMB(tot_herbivory ~ precip*clip*excl + (1|sampID) + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.herb.all.sum <- summary(zi.srer.herb.all)
zi.srer.herb.all.sum

# model comparisons/ranks/weights by AIC
modelcomp.herb <- MuMIn::dredge(zi.srer.herb.all, rank = "AIC")

# examine model averages
MuMIn::model.avg(modelcomp.germ)

# extract only the variance correlation matrix
VarCorr(best.model.germ)

# examine model residuals
seedlings_obs$res_germ <- residuals(best.model.germ, quantileFunction = qpois)

# examine histogram of model residuals
hist(seedlings_obs$res_germ)

# extract model co-variance matrix
vcov.matrix.germ <- vcov(best.model.germ)

# extract conditional model coefficients
cond.cov.germ <- as.data.frame(vcov.matrix.germ$cond)

# extract zero inflation model coefficients (probably not needed)
zi.cov.germ <- as.data.frame(vcov.matrix.germ$zi)

# use Dharma package to examine model fit
seedlings.best.simres.germ <- simulateResiduals(best.model.germ, n = 1000, plot = TRUE, integerResponse = TRUE)

seedlings_obs$sim_germ <- seedlings.best.simres.germ$scaledResiduals

hist(seedlings_obs$sim_germ)

testOutliers(seedlings.best.simres.germ, type = 'bootstrap')

testDispersion(seedlings.best.simres.germ)

testUniformity(seedlings.best.simres.germ)

testResiduals(seedlings.best.simres.germ)

testTemporalAutocorrelation(seedlings.best.simres.germ)

testZeroInflation(seedlings.best.simres.germ)

plotResiduals(seedlings.best.simres.germ, seedlings_obs$excl, quantreg = T)

plot(seedlings.best.simres.germ)

hist(seedlings.best.simres.germ)

seedlings_obs$pred_germ <- predict(best.model.germ, type = "response")

hist(seedlings_obs$pred_germ)

# post-hoc test of within precip treatment differences
g.precip <- glht(best.model.germ, linfct = mcp(precip = "Tukey", interaction_average = T))
summary(g.precip)


g2 <- car::Anova(best.model.germ, type = "III")
summary(g2)

# another method for post-hoc test of precip and exclusion 
post.hoc <- emmeans::emmeans(best.model.germ, specs = ~precip*excl)

# get lettering report on post-hoc test
post.hoc.letters <- cld(post.hoc, Letters = letters, covar = T)

# write model outputs

# sink("Model_Outputs/model-outputs-germination.txt")
# 
# zi.block.sum
# zi.srer.germ.p.sum
# zi.srer.germ.pc.sum
# zi.srer.germ.pe.sum
# zi.srer.germ.pce.sum
# zi.srer.germ.full.sum
# zi.srer.germ.pc.int.sum
# zi.srer.germ.pce.int.sum
# zi.srer.germ.pe.int.sum
# zi.srer.germ.ce.int.sum
# aic.compare.final
# modelcomp.germ
# zi.srer.germ.aic.sum
# post.hoc
# post.hoc.letters
# 
# sink()

# plot observed germination values and predicted germination values
seedlings_obs %>%
  ggplot(aes(x = 10*tot_germination, y = 10*pred_germ)) +
  geom_point() +
  geom_smooth(method = "lm", formula = (y) ~ (x)) +
  facet_wrap(~precip)

ranef(best.model.germ)

coef(best.model.germ)


