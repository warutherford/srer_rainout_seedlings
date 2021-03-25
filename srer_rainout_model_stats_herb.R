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

# Herbivory counts
seedlings_herb_full <- seedlings %>% 
  group_by(block, precip, clip, excl, side, rep, date, cohort) %>%
  count(herbivory) %>% 
  pivot_wider(names_from = herbivory,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_herb = "0",
         herb_lived = "1",
         herb_died = "2") %>% 
  mutate(tot_herbivory = herb_lived + herb_died)

glimpse(seedlings_herb_full)
summary(seedlings_herb_full)

# descriptive stats for each cohort 1-3
describeBy(seedlings_herb_full, group = "cohort")

# make a new ID column for each observation and plot
seedlings_obs_herb <- seedlings_herb_full %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("sampID", block:rep, sep = "_", remove = FALSE) %>%
  unite("plotID", block:precip, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         date = as.factor(date),
         sampID = as.factor(sampID),
         plotID = as.factor(plotID))

# mixed effects model with nesting, sampID as random and date (cohort) within year for temporal autocorrelation
hist(seedlings_obs_herb$tot_herbivory)
hist(log(seedlings_obs_herb$tot_herbivory)) # still zero-inflated
hist(sqrt(seedlings_obs_herb$tot_herbivory)) # still zero-inflated

#glmmTMB function to build zero-inflated model, binomial dist., random factors

# Based on AIC and Log Likelihood, and model convergence/overdispersion issues,
# accounting for sample independence and each year plus temp 
# are the main random factors to include in the models
# (e.g, (1|cohort) + (1|sampID) + ar1(date + 0|cohort))
# Thus, cohort, sample ID and temp autocorrelation are used in building final model
# Keep same random factors as with survival model

# test if tot_herbivory sig diff across all blocks?
zi.block <- glmmTMB(tot_herbivory ~ block + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                    data = seedlings_obs_herb,
                    family = binomial)
zi.block.sum <- summary(zi.block)
zi.block.sum

# not all blocks sig, so exclude from potential models

# start with full model, and make simpler
# precipitation, clipping, and exclusion total interactions as fixed factors
zi.srer.herb.full <- glmmTMB(tot_herbivory ~ precip/clip/excl + (1|cohort) + (1|sampID),# + ar1(date + 0|cohort),
                             data = seedlings_obs_herb,
                             family = binomial)
zi.srer.herb.full.sum <- summary(zi.srer.herb.full)
zi.srer.herb.full.sum

# interactions not significant and not informative, use each tx individually
# precipitation, clipping, and exclusion fixed factors
# model convergence issue with cohort and sample with temp autocorrelation
# try neg binomial dist instead of poisson
zi.srer.herb.pce <- glmmTMB(tot_herbivory ~ precip + clip + excl + (1|cohort) + (1|sampID),# + ar1(date + 0|cohort),
                            data = seedlings_obs_herb,
                            family = binomial)
zi.srer.herb.pce.sum <- summary(zi.srer.herb.pce)
zi.srer.herb.pce.sum

# clipping not sig, remove from model
# precipitation and exclusion fixed factors
zi.srer.herb.pe <- glmmTMB(tot_herbivory ~ precip + excl + (1|cohort) + (1|sampID),# + ar1(date + 0|cohort),
                           data = seedlings_obs_herb,
                           family = binomial)
zi.srer.herb.pe.sum <- summary(zi.srer.herb.pe)
zi.srer.herb.pe.sum

# precipitation only fixed factor
zi.srer.herb.p <- glmmTMB(tot_herbivory ~ precip + (1|cohort) + (1|sampID),# + ar1(date + 0|cohort),
                          data = seedlings_obs_herb,
                          family = binomial)
zi.srer.herb.p.sum <- summary(zi.srer.herb.p)
zi.srer.herb.p.sum

# precipitation, exclusion, and precipitation and exclusion interaction fixed factors
zi.srer.herb.pe.int <- glmmTMB(tot_herbivory ~ precip + precip/excl + (1|cohort) + (1|sampID),# + ar1(date + 0|cohort),
                               data = seedlings_obs_herb,
                               family = binomial)
zi.srer.herb.pe.int.sum <- summary(zi.srer.herb.pe.int)
zi.srer.herb.pe.int.sum

# compare AIC scores of all potential models for model selection
aic.compare.final <- AICtab(zi.srer.herb.p,
                            zi.srer.herb.pe,
                            zi.srer.herb.pce,
                            zi.srer.herb.full,
                            zi.srer.herb.pe.int,
                            logLik = TRUE)

aic.compare.final

# does removing clipping sig improve model?
anova(zi.srer.herb.pce, zi.srer.herb.pe) # yes

# final PRVE seedlings tot_herbivory model will precipitation and exclusion as only fixed factors
zi.srer.herb.final <- glmmTMB(tot_herbivory ~ precip + excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                              data = seedlings_obs_herb,
                              family = binomial)
zi.srer.herb.final.sum <- summary(zi.srer.herb.final)
zi.srer.herb.final.sum

# define new variable to be best model for diagnostics
best.model.herb <- zi.srer.herb.final

# extract only the variance correlation matrix
VarCorr(best.model.herb)

# examine model residuals
seedlings_obs$res_herb <- residuals(best.model.herb, quantileFunction = qpois)

# examine histogram of model residuals
hist(seedlings_obs$res_herb) # look normal

# extract model co-variance matrix
vcov.matrix.herb <- vcov(best.model.herb)

# extract conditional model coefficients
cond.cov.herb <- as.data.frame(vcov.matrix.herb$cond)

# extract zero inflation model coefficients (probably not needed)
zi.cov.herb <- as.data.frame(vcov.matrix.herb$zi)

# use Dharma package to examine model fit
seedlings.best.simres.herb <- simulateResiduals(best.model.herb, n = 1000, plot = TRUE, integerResponse = TRUE)

# save scaled residuals
seedlings_obs$sim_herb <- seedlings.best.simres.herb$scaledResiduals

# save model fitted residuals
seedlings_obs$sim_fit_herb <- seedlings.best.simres.herb$fittedResiduals

# look at simulated/model residuals
hist(seedlings_obs$sim_herb) # pretty level, slight inc at end but makes sense with slight over prediction

hist(seedlings_obs$sim_fit_herb) # normal residuals

# are there outliers?
testOutliers(seedlings.best.simres.herb, type = 'bootstrap') # no outliers

# is the model overly dispersed
testDispersion(seedlings.best.simres.herb) # not overly dispersed

# KS test
testUniformity(seedlings.best.simres.herb) # slight hump in middle, observed tot_herbivory slightly higher than model expected

# model account for temp autocorrelation
testTemporalAutocorrelation(seedlings.best.simres.herb) # temporal autocorrelation is accounted for in model

# model account for zero inflation?
testZeroInflation(seedlings.best.simres.herb) # zero inflation accounted for in the model

# does the model mean prediction fit?
means <- function(x) mean(x) 
testGeneric(seedlings.best.simres.herb, summary = means) # yes, does not sig deviate (p = 0.77)

# does the model standard deviation prediction fit?
spread <- function(x) sd(x)
testGeneric(seedlings.best.simres.herb, summary = spread) # yes, does not sig deviate (p = 0.22)

# take a look at model residuals
plotResiduals(seedlings.best.simres.herb, seedlings_obs$precip, quantreg = T) # no strange residual patterns

# save model predicted responses as a new column
seedlings_obs$pred_herb <- predict(best.model.herb, type = "response")

# look histograms of model predicted tot_herbivory and original tot_herbivory data 
hist(seedlings_obs$pred_herb)

hist(seedlings_obs$tot_herbivory)

# post-hoc test of within precip treatment differences
g.precip <- glht(best.model.herb, linfct = mcp(precip = "Tukey", interaction_average = T))
summary(g.precip) # drought treatments sig different from control and wet

# post-hoc test of within exclusion treatment differences
g.excl <- glht(best.model.herb, linfct = mcp(excl = "Tukey", interaction_average = T))
summary(g.excl) # total treatments sig different from control, ant and rodents not a big influence

# another method for post-hoc test of precip and exclusion 
post.hoc <- emmeans::emmeans(best.model.herb, specs = ~precip*excl)
post.hoc

# get lettering report on post-hoc test
post.hoc.letters <- cld(post.hoc, Letters = letters, covar = T)
post.hoc.letters

# in case random effect coeff and model coeff are needed
ranef(best.model.herb)

coef(best.model.herb)

# write seedlings_obs data frame to csv for use/graphing
write.csv(seedlings_obs, file = "Data/seedlings_obs_herb.csv", row.names = FALSE)

# write model outputs to text file

# sink("Model_Outputs/model-outputs-tot_herbivory.txt")
# 
# zi.srer.herb.p.sum
# zi.srer.herb.pe.sum
# zi.srer.herb.pce.sum
# zi.srer.herb.full.sum
# zi.srer.herb.pe.int.sum
# aic.compare.final
# zi.srer.herb.final
# post.hoc
# post.hoc.letters
# 
# sink()

