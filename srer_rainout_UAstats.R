# USDA Rainout - Santa Rita Experimental Range
# Stats and Modeling for UA consults
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-12-14

# Load packages
library(tidyverse)
library(vroom)
library(psych)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(multcomp)

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
  unite("plotID",block:precip, sep = "_", remove = FALSE) %>%
  unite("sampID", block:excl, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         plotID = as.factor(plotID),
         date = as.factor(date),
         Obs_ID = as.factor(Obs_ID),
         sampID = as.factor(sampID))

# Fate counts for complete dataset (include side/rep)
seedlings_fate_full_comp <- seedlings %>% 
  group_by(block, precip, clip, excl, side, date, cohort) %>%
  count(fate) %>% 
  pivot_wider(names_from = fate,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_germ = "0",
         survival = "1",
         died = "2") %>% 
  mutate(tot_germination = survival + died)

# make a new ID column for each observation and plot for the full dataset
seedlings_obs_comp <- seedlings_fate_full_comp %>% 
  rowid_to_column("Obs_ID") %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("plotID",block:precip, sep = "_", remove = FALSE) %>%
  unite("sampID", block:side, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         plotID = as.factor(plotID),
         date = as.factor(date),
         Obs_ID = as.factor(Obs_ID),
         sampID = as.factor(sampID))

#glmmTMB function to build zero-inflated model, poisson dist., random factors

# test if survival sig diff across all blocks?
zi.block <- glmmTMB(survival ~ block + (1|sampID) + (1|cohort/date),
                    data = seedlings_obs,
                    ziformula = ~1,
                    family = poisson)
zi.block.sum <- summary(zi.block)
zi.block.sum

#not all blocks sig, so exclude from model

# recreate model similar to colleagues, no nesting but plot and year random
zi.jrn <- glmmTMB(survival ~ precip + clip + excl + (1|plotID) + (1|cohort),
                  data = seedlings_obs,
                  ziformula = ~1,
                  family = poisson)
zi.jrn.sum <- summary(zi.jrn)
zi.jrn.sum

# step through random factors

# intercept only, sample ID as random factor
zi.srer.fix.1 <- glmmTMB(survival ~ 1 + (1|sampID),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.1 <- summary(zi.srer.fix.1)
zi.srer.fix.sum.1

# compare with zi poisson
zi.srer.fix.1.nbi <- glmmTMB(survival ~ 1 + (1|sampID),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = nbinom2)
zi.srer.fix.sum.1.nbi <- summary(zi.srer.fix.1.nbi)
zi.srer.fix.sum.1.nbi

zi.srer.fix.2 <- glmmTMB(survival ~ 1 + (1|cohort),
                         data = seedlings_obs,
                         ziformula = ~1,
                         family = poisson)
zi.srer.fix.sum.2 <- summary(zi.srer.fix.2)
zi.srer.fix.sum.2

zi.srer.fix.3 <- glmmTMB(survival ~ 1 + (1|sampID) + (1|cohort),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.3 <- summary(zi.srer.fix.3)
zi.srer.fix.sum.3

zi.srer.fix.4 <- glmmTMB(survival ~ 1 + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.4 <- summary(zi.srer.fix.4)
zi.srer.fix.sum.4

zi.srer.fix.5 <- glmmTMB(survival ~ 1 + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.5 <- summary(zi.srer.fix.5)
zi.srer.fix.sum.5

zi.srer.fix.full.nbi <- glmmTMB(survival ~ precip+clip+excl + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~precip+clip+excl,
                         family = nbinom2)
zi.srer.fix.sum.full.nbi <- summary(zi.srer.fix.full.nbi)
zi.srer.fix.sum.full.nbi

zi.srer.fix.auto <- glmmTMB(survival ~ 1 + (1|sampID) + ar1(date + 0|cohort),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.auto <- summary(zi.srer.fix.auto)
zi.srer.fix.sum.auto

zi.srer.fix.auto1 <- glmmTMB(survival ~ 1 + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.fix.sum.auto1 <- summary(zi.srer.fix.auto1)
zi.srer.fix.sum.auto1


zi.srer.fix.auto2 <- glmmTMB(survival ~ precip + clip + excl + (1|sampID) + ar1(date + 0|cohort),
                            data = seedlings_obs,
                            ziformula = ~.,
                            family = poisson)
zi.srer.fix.sum.auto2 <- summary(zi.srer.fix.auto2)
zi.srer.fix.sum.auto2

zi.srer.fix.auto3 <- glmmTMB(survival ~ precip + precip*clip + excl + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.fix.sum.auto3 <- summary(zi.srer.fix.auto3)
zi.srer.fix.sum.auto3

zi.srer.fix.auto4 <- glmmTMB(survival ~ precip + clip + precip*excl + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.fix.sum.auto4 <- summary(zi.srer.fix.auto4)
zi.srer.fix.sum.auto4

zi.srer.fix.auto5 <- glmmTMB(survival ~ precip*clip*excl + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.fix.sum.auto5 <- summary(zi.srer.fix.auto5)
zi.srer.fix.sum.auto5

zi.srer.fix.auto6 <- glmmTMB(survival ~ precip*clip*excl + ar1(date + 0|cohort),
                             data = seedlings_obs,
                             ziformula = ~.,
                             family = poisson)
zi.srer.fix.sum.auto6 <- summary(zi.srer.fix.auto6)
zi.srer.fix.sum.auto6

zi.srer.fix.full <- glmmTMB(survival ~ precip*clip*excl + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.full <- summary(zi.srer.fix.full)
zi.srer.fix.sum.full

zi.srer.fix.6 <- glmmTMB(survival ~ precip + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.6 <- summary(zi.srer.fix.6)
zi.srer.fix.sum.6

zi.srer.fix.7 <- glmmTMB(survival ~ precip + precip*clip + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.7 <- summary(zi.srer.fix.7)
zi.srer.fix.sum.7

zi.srer.fix.8 <- glmmTMB(survival ~ precip + clip + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.8 <- summary(zi.srer.fix.8)
zi.srer.fix.sum.8

zi.srer.fix.9 <- glmmTMB(survival ~ precip + clip + excl + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = poisson)
zi.srer.fix.sum.9 <- summary(zi.srer.fix.9)
zi.srer.fix.sum.9

zi.srer.fix.10.nbi <- glmmTMB(survival ~ precip + clip + precip*excl + (1|sampID) + (1|cohort/date),
                         data = seedlings_obs,
                         ziformula = ~.,
                         family = nbinom2)
zi.srer.fix.sum.10.nbi <- summary(zi.srer.fix.10.nbi)
zi.srer.fix.sum.10.nbi

zi.srer.fix.10 <- glmmTMB(survival ~ precip + clip + precip*excl + (1|sampID) + (1|cohort/date),
                          data = seedlings_obs,
                          ziformula = ~.,
                          family = poisson)
zi.srer.fix.sum.10 <- summary(zi.srer.fix.10)
zi.srer.fix.sum.10

anova(zi.srer.fix.1,
      zi.srer.fix.2,
      zi.srer.fix.3,
      zi.srer.fix.4,
      zi.srer.fix.5,
      zi.srer.fix.full,
      zi.srer.fix.full.nbi,
      zi.srer.fix.auto,
      zi.srer.fix.auto1,
      zi.srer.fix.auto2,
      zi.srer.fix.auto3,
      zi.srer.fix.auto4,
      zi.srer.fix.auto5,
      zi.srer.fix.auto6,
      zi.srer.fix.6,
      zi.srer.fix.7,
      zi.srer.fix.8,
      zi.srer.fix.9,
      zi.srer.fix.10,
      zi.srer.fix.10.nbi,
      test = "Chisq")


buildglmmTMB(survival ~ precip + clip + excl + block +
               precip*clip + precip*excl + excl*clip + precip*clip*excl + cohort +
               precip*clip*block*cohort +
               (1|sampID) + (1|cohort/sampID) + (1|cohort/date) + ar1(date + 0|cohort),
             data = seedlings_obs,
             ziformula = ~.,
             family = poisson,
             crit = "AIC",
             elim = "AIC")

zi.srer.surv <- glmmTMB(survival ~ precip + excl + (1|sampID) + ar1(date + 0|cohort),
                       data = seedlings_obs,
                       ziformula = ~.,
                       family = poisson)
zi.srer.surv.sum <- summary(zi.srer.surv)
zi.srer.surv.sum


zi.srer.cohort <- glmmTMB(survival ~ precip + excl + (1|cohort/sampID),
                        data = seedlings_obs,
                        ziformula = ~.,
                        family = poisson)
zi.srer.cohort.sum <- summary(zi.srer.cohort)
zi.srer.cohort.sum

zi.srer.cohort.full <- glmmTMB(survival ~ precip + clip + excl + (1|cohort/sampID),
                          data = seedlings_obs,
                          ziformula = ~.,
                          family = poisson)
zi.srer.cohort.full.sum <- summary(zi.srer.cohort.full)
zi.srer.cohort.full.sum


zi.srer.surv1 <- glmmTMB(survival ~ precip + clip + excl + (1|sampID) + ar1(date + 0|cohort),
                        data = seedlings_obs,
                        ziformula = ~.,
                        family = poisson)
zi.srer.surv1.sum <- summary(zi.srer.surv1)
zi.srer.surv1.sum


buildglmmTMB(tot_germination ~ precip + clip + excl + 
               precip*clip + precip*excl + excl*clip + precip*clip*excl + 
               (1|sampID) + (1|cohort/sampID) + (1|cohort/date) + ar1(date + 0|cohort),
             data = seedlings_obs,
             ziformula = ~1,
             family = poisson,
             crit = "AIC",
             elim = "AIC")

zi.srer.germ <- glmmTMB(tot_germination ~ precip + excl + clip + precip:clip + excl:clip + (1|sampID) + ar1(date + 0|cohort),
                        data = seedlings_obs,
                        ziformula = ~.,
                        family = poisson)
zi.srer.germ.sum <- summary(zi.srer.germ)
zi.srer.germ.sum


zi.srer.fin2 <- glmmTMB(survival ~ precip + excl + clip + (1|sampID) + ar1(date + 0|cohort),
                       data = seedlings_obs,
                       ziformula = ~.,
                       family = poisson(link = "log"))
summary(zi.srer.fin2)

# calculate differences b/w AIC for a number of models
aic.compare <- AICtab(zi.jrn,
       zi.srer.1,
       zi.srer.1.nbi,
       zi.srer.2,
       zi.srer.3,
       zi.srer.3.nbi,
       zi.srer.4,
       zi.srer.5,
       zi.srer.6,
       zi.srer.7,
       zi.srer.8,
       zi.srer.9,
       logLik = TRUE)

sink("model-outputs.txt")

zi.jrn.sum
zi.srer.sum.1
zi.srer.sum.1.nbi
zi.srer.sum.2
zi.srer.sum.3
zi.srer.sum.3.nbi
zi.srer.sum.4
zi.srer.sum.5
zi.srer.sum.6
aic.compare

sink()
unlink("model-outputs.txt")

best.model <- zi.srer.cohort.full

# with nesting, also get interactions? summary looks like yes...

modelcomp <- MuMIn::dredge(best.model)

MuMIn::model.avg(modelcomp)

VarCorr(best.model)
seedlings_obs$res <- residuals(best.model, quantileFunction = qpois)

hist(seedlings_obs$res)

vcov.matrix <- vcov(best.model)
cond.cov <- as.data.frame(vcov.matrix$cond)
zi.cov <- as.data.frame(vcov.matrix$zi)


seedlings.best.simres <- simulateResiduals(best.model, n = 1000, plot = TRUE, integerResponse = TRUE)

seedlings_obs$sim <- seedlings.best.simres$scaledResiduals

hist(seedlings_obs$sim)


# group.sim <- recalculateResiduals(seedlings.best.simres, group = seedlings_obs$cohort)
# 
# plot(group.sim)

testOutliers(seedlings.best.simres, type = 'bootstrap')

testDispersion(seedlings.best.simres)

testUniformity(seedlings.best.simres)

testResiduals(seedlings.best.simres)

testTemporalAutocorrelation(seedlings.best.simres)

testZeroInflation(seedlings.best.simres)

plotResiduals(seedlings.best.simres, seedlings_obs$excl)

plot(seedlings.best.simres)

hist(seedlings.best.simres)

seedlings_obs$pred <- predict(best.model, type = "response")

hist(seedlings_obs$pred)

g.precip <- glht(best.model, linfct = mcp(precip = "Tukey", interaction_average = T))
summary(g.precip)


g.excl <- glht(best.model, linfct = mcp(excl = "Tukey", interaction_average = T))
summary(g.excl)


g2 <- car::Anova(best.model, type = "III")
summary(g2)


post.hoc <- emmeans::emmeans(best.model, specs = ~precip*excl)

cld(post.hoc, Letters = letters, covar = T)

seedlings_obs %>%
  ggplot(aes(x = 10*survival, y = 10*pred, color = excl)) +
  geom_point() +
  geom_smooth(method = "glm", formula = (y) ~ sqrt(x)) +
  facet_wrap(~precip)

ranef(best.model)
coef(best.model)
