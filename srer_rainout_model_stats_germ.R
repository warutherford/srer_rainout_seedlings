# USDA Rainout - Santa Rita Experimental Range
# Stats and Modeling-Germination
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

# Based on AIC and Log Likelihood, sample ID and temp autocorrelation random factors accounts for most variation
# Thus, sample ID and temp autocorrelation is used in building final model
# See script "srer_rainout_model_stats_survival.R" for code to justify random factors

# look at distribution of germination counts
hist(seedlings_obs$tot_germination)

# test if survival sig diff across all blocks?
zi.block <- glmmTMB(tot_germination ~ block + (1|sampID) + ar1(date + 0|cohort),
                    data = seedlings_obs,
                    ziformula = ~.,
                    family = poisson)
zi.block.sum <- summary(zi.block)
zi.block.sum

#not all blocks sig, so exclude from potential models

# individual models with treatments using same random effects
# precipitation only
zi.srer.germ.p <- glmmTMB(tot_germination ~ precip + (1|sampID) + ar1(date + 0|cohort),
                          data = seedlings_obs,
                          ziformula = ~.,
                          family = poisson)
zi.srer.germ.p.sum <- summary(zi.srer.germ.p)
zi.srer.germ.p.sum

# precipitation and clipping only
zi.srer.germ.pc <- glmmTMB(tot_germination ~ precip + clip + (1|sampID) + ar1(date + 0|cohort),
                           data = seedlings_obs,
                           ziformula = ~.,
                           family = poisson)
zi.srer.germ.pc.sum <- summary(zi.srer.germ.pc)
zi.srer.germ.pc.sum

# precipitation and exclusion only
zi.srer.germ.pe <- glmmTMB(tot_germination ~ precip + excl + (1|sampID) + ar1(date + 0|cohort),
                           data = seedlings_obs,
                           ziformula = ~.,
                           family = poisson)
zi.srer.germ.pe.sum <- summary(zi.srer.germ.pe)
zi.srer.germ.pe.sum

# precipitation, clipping, and exclusion
zi.srer.germ.pce <- glmmTMB(tot_germination ~ precip + clip + excl + (1|sampID) + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.germ.pce.sum <- summary(zi.srer.germ.pce)
zi.srer.germ.pce.sum

# all interactions b/w precipitation, clipping, and exclusion 
zi.srer.germ.full <- glmmTMB(tot_germination ~ precip/clip/excl + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.germ.full.sum <- summary(zi.srer.germ.full)
zi.srer.germ.full.sum

# precipitation, clipping and precipitation interaction, and exclusion
zi.srer.germ.pc.int <- glmmTMB(tot_germination ~ precip + precip/clip + excl + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs,
                               ziformula = ~.,
                               family = poisson)
zi.srer.germ.pc.int.sum <- summary(zi.srer.germ.pc.int)
zi.srer.germ.pc.int.sum

# precipitation, clipping, and precipitation and exclusion interaction
zi.srer.germ.pce.int <- glmmTMB(tot_germination ~ precip + clip + precip/excl + (1|sampID) + ar1(date + 0|cohort),
                                data = seedlings_obs,
                                ziformula = ~.,
                                family = poisson)
zi.srer.germ.pce.int.sum <- summary(zi.srer.germ.pce.int)
zi.srer.germ.pce.int.sum

# precipitation, clipping, and clipping and exclusion interaction
zi.srer.germ.ce.int <- glmmTMB(tot_germination ~ precip + clip + clip/excl + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs,
                               ziformula = ~.,
                               family = poisson)
zi.srer.germ.ce.int.sum <- summary(zi.srer.germ.ce.int)
zi.srer.germ.ce.int.sum

# precipitation, exclusion, and precipiration and exlclusion interaction
zi.srer.germ.pe.int <- glmmTMB(tot_germination ~ precip + excl + precip/excl + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs,
                               ziformula = ~.,
                               family = poisson)
zi.srer.germ.pe.int.sum <- summary(zi.srer.germ.pe.int)
zi.srer.germ.pe.int.sum

# compare all AIC scores for model selection
aic.compare.final <- AICtab(zi.srer.germ.p,
                            zi.srer.germ.pc,
                            zi.srer.germ.pe,
                            zi.srer.germ.pce,
                            zi.srer.germ.full,
                            zi.srer.germ.pc.int,
                            zi.srer.germ.pce.int,
                            zi.srer.germ.pe.int,
                            zi.srer.germ.ce.int,
                            logLik = TRUE)

aic.compare.final

# best model based on AIC score
zi.srer.germ.aic <- glmmTMB(tot_germination ~ precip + excl + (1|sampID) + ar1(date + 0|cohort),
                        data = seedlings_obs,
                        ziformula = ~.,
                        family = poisson)
zi.srer.germ.aic.sum <- summary(zi.srer.germ.aic)
zi.srer.germ.aic.sum

# define new variable to be best model for diagnostics
best.model.germ <- zi.srer.germ.aic

### below is MuMIn model selection procedure with dredge function

zi.srer.germ.all <- glmmTMB(tot_germination ~ precip*clip*excl + (1|sampID) + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.germ.all.sum <- summary(zi.srer.germ.all)
zi.srer.germ.all.sum

# model comparisons/ranks/weights by AIC
modelcomp.germ <- MuMIn::dredge(zi.srer.germ.all, rank = "AIC")

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


