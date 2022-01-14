# USDA Rainout - Santa Rita Experimental Range
# Stats and Modeling-Germination
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-01-26

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

### Read in seedlings data (cohorts 1-3), make all columns factor and date a date
seedlings <- vroom("Data/seedlings_combined.csv",
                   col_select = -c(1),
                   col_types = c(.default = "f",
                                 date = "D"))
str(seedlings)

# For germination data, either the seed germinated or it didn't
# Create data frame for presence/absence (1 or 0) of germination, binomial distribution
seedlings_germ_full <- seedlings %>% 
  group_by(block, precip, clip, excl, side, rep, date, cohort) %>%
  count(fate) %>% 
  pivot_wider(names_from = fate,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_germ = "0",
         survival = "1",
         died = "2") %>% 
  mutate(tot_germination = (survival + died))

# Create possible random factor variables and fix any data structures needed for modeling
seedlings_obs_germ <- seedlings_germ_full %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("plotID",block:precip, sep = "_", remove = FALSE) %>%
  unite("sampID", block:rep, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         plotID = as.factor(plotID),
         date = as.factor(date),
         sampID = as.factor(sampID))

# descriptive stats for each cohort 1-3
describeBy(seedlings_obs_germ , group = "cohort")

# mixed effects model with nesting, sampID as random and date (cohort) within year for temporal autocorrelation
hist(seedlings_obs_germ$tot_germination)# binomial (0 or 1)

#glmmTMB function to build zero-inflated model, binomial dist., random factors

# Based on AIC and Log Likelihood, and model convergence/overdispersion issues,
# accounting for sample independence and each year plus temp 
# are the main random factors to include in the models
# (e.g, (1|cohort) + (1|sampID) + ar1(date + 0|cohort))
# Thus, cohort, sample ID and temp autocorrelation are used in building final model
# Keep same random factors as with survival model

# test if tot_germination sig diff across all blocks?
zi.block <- glmmTMB(tot_germination ~ block + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                    data = seedlings_obs_germ,
                    family = binomial(link = "logit"))
zi.block.sum <- summary(zi.block)
zi.block.sum

# not all blocks sig, so exclude from potential models

# start with full model, and make simpler
# precipitation, clipping, and exclusion total interactions as fixed factors
zi.srer.germ.full <- glmmTMB(tot_germination ~ precip + precip/clip+ precip/excl + excl + clip + precip/clip/excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                             data = seedlings_obs_germ,
                             family = binomial(link = "logit"))
zi.srer.germ.full.sum <- summary(zi.srer.germ.full)
zi.srer.germ.full.sum

# type II wald's for fixed effect significance
Anova(zi.srer.germ.full)

# interactions not all significant and not informative, use each tx individually
# precipitation, clipping, and exclusion fixed factors
# model convergence issue with cohort and sample with temp autocorrelation
zi.srer.germ.pce <- glmmTMB(tot_germination ~ precip + clip + excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                            data = seedlings_obs_germ,
                            family = binomial(link = "logit"))
zi.srer.germ.pce.sum <- summary(zi.srer.germ.pce)
zi.srer.germ.pce.sum

# type II wald's for fixed effect significance
Anova(zi.srer.germ.pce)

# clipping not sig, remove from model
# precipitation and exclusion fixed factors
zi.srer.germ.pe <- glmmTMB(tot_germination ~ precip + excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                           data = seedlings_obs_germ,
                           family = binomial(link = "logit"))
zi.srer.germ.pe.sum <- summary(zi.srer.germ.pe)
zi.srer.germ.pe.sum

# type II wald's for fixed effect significance
Anova(zi.srer.germ.pe)

# precipitation only fixed factor
zi.srer.germ.p <- glmmTMB(tot_germination ~ precip + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                          data = seedlings_obs_germ,
                          family = binomial(link = "logit"))
zi.srer.germ.p.sum <- summary(zi.srer.germ.p)
zi.srer.germ.p.sum

# type II wald's for fixed effect significance
Anova(zi.srer.germ.p)

# precipitation, exclusion, and precipitation and exclusion interaction fixed factors
zi.srer.germ.pe.int <- glmmTMB(tot_germination ~ precip + precip/excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs_germ,
                               family = binomial(link = "logit"))
zi.srer.germ.pe.int.sum <- summary(zi.srer.germ.pe.int)
zi.srer.germ.pe.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.germ.pe.int)

# precipitation, clip, and precipitation and clip interaction fixed factors
zi.srer.germ.pc.int <- glmmTMB(tot_germination ~ precip + precip/clip + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                               data = seedlings_obs_germ,
                               family = binomial(link = "logit"))
zi.srer.germ.pc.int.sum <- summary(zi.srer.germ.pc.int)
zi.srer.germ.pc.int.sum

# type II wald's for fixed effect significance
Anova(zi.srer.germ.pc.int)

# compare AIC scores of all potential models for model selection
aic.compare.final <- AICtab(zi.srer.germ.p,
                            zi.srer.germ.pe,
                            zi.srer.germ.pce,
                            zi.srer.germ.full,
                            zi.srer.germ.pe.int,
                            zi.srer.germ.pc.int,
                            logLik = TRUE)

aic.compare.final

# does removing clipping sig improve model?
anova(zi.srer.germ.pce, zi.srer.germ.pe) # not sig, but just pe is a more parsimonious model

# final PRVE seedlings tot_germination model will precipitation and exclusion as only fixed factors
zi.srer.germ.final <- glmmTMB(tot_germination ~ precip + excl + (1|cohort) + (1|sampID) + ar1(date + 0|cohort),
                              data = seedlings_obs_germ,
                              family = binomial(link = "logit"))
zi.srer.germ.final.sum <- summary(zi.srer.germ.final)
zi.srer.germ.final.sum

# define new variable to be best model for diagnostics
best.model.germ <- zi.srer.germ.final

# extract only the variance correlation matrix
VarCorr(best.model.germ)

# examine model residuals
seedlings_obs_germ$res_germ <- residuals(best.model.germ, quantileFunction = qpois)

# examine histogram of model residuals
hist(seedlings_obs_germ$res_germ) # look normal

# extract model co-variance matrix
vcov.matrix.germ <- vcov(best.model.germ)

# extract conditional model coefficients
cond.cov.germ <- as.data.frame(vcov.matrix.germ$cond)

# extract zero inflation model coefficients (probably not needed)
zi.cov.germ <- as.data.frame(vcov.matrix.germ$zi)

# use Dharma package to examine model fit
seedlings.best.simres.germ <- simulateResiduals(best.model.germ, n = 1000, plot = TRUE, integerResponse = TRUE)

# save scaled residuals
seedlings_obs_germ$sim_germ <- seedlings.best.simres.germ$scaledResiduals

# save model fitted residuals
seedlings_obs_germ$sim_fit_germ <- seedlings.best.simres.germ$fittedResiduals

# look at simulated/model residuals
hist(seedlings_obs_germ$sim_germ) # pretty level, slight inc at end but makes sense with slight over prediction

hist(seedlings_obs_germ$sim_fit_germ) # normal residuals

# are there outliers?
testOutliers(seedlings.best.simres.germ, type = 'bootstrap') # no outliers

# is the model overly dispersed
testDispersion(seedlings.best.simres.germ) # not overly dispersed

# KS test
testUniformity(seedlings.best.simres.germ) # slight hump in middle, observed tot_germination slightly higher than model expected

# model account for temp autocorrelation
testTemporalAutocorrelation(seedlings.best.simres.germ) # temporal autocorrelation is accounted for in model

# model account for zero inflation?
testZeroInflation(seedlings.best.simres.germ) # zero inflation accounted for in the model

# does the model mean prediction fit?
means <- function(x) mean(x) 
testGeneric(seedlings.best.simres.germ, summary = means) # yes, does not sig deviate (p = 0.37)

# does the model standard deviation prediction fit?
spread <- function(x) sd(x)
testGeneric(seedlings.best.simres.germ, summary = spread) # yes, does not sig deviate (p = 0.67)

# take a look at model residuals
plotResiduals(seedlings.best.simres.germ, seedlings_obs$precip, quantreg = T) # no strange residual patterns

# save model predicted responses as a new column
seedlings_obs_germ$pred_germ <- predict(best.model.germ, type = "response")

# look histograms of model predicted tot_germination and original tot_germination data 
hist(seedlings_obs_germ$pred_germ)

hist(seedlings_obs_germ$tot_germination)

# another method for post-hoc test of precip and exclusion 
post.hoc <- emmeans::emmeans(best.model.germ, specs = ~precip*excl)
post.hoc

# get lettering report on post-hoc test
post.hoc.letters <- cld(post.hoc, Letters = letters, covar = T)
post.hoc.letters

# in case random effect coeff and model coeff are needed
ranef(best.model.germ)

coef(best.model.germ)

# write seedlings_obs data frame to csv for use/graphing
write.csv(seedlings_obs_germ, file = "Data/seedlings_obs_germ.csv", row.names = FALSE)

# write model outputs to text file

# sink("Model_Outputs/model-outputs-tot_germination.txt")
# 
# zi.srer.germ.p.sum
# zi.srer.germ.pe.sum
# zi.srer.germ.pce.sum
# zi.srer.germ.full.sum
# zi.srer.germ.pe.int.sum
# aic.compare.final
# zi.srer.germ.final
# Anova(zi.srer.germ.full)
# Anova(zi.srer.germ.pce)
# Anova(zi.srer.germ.pe)
# Anova(zi.srer.germ.p)
# Anova(zi.srer.germ.pe.int)
# Anova(zi.srer.germ.pc.int)
# post.hoc
# post.hoc.letters
# 
# sink()

### if want to look at effects treating monsoon precip as continuous
precip_cont_germ_df_1 <- seedlings_obs_germ %>% 
  group_by(cohort) %>%
  filter(cohort == "1") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "280",
                                     "IR" = "462",
                                     "RO" = "98")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

precip_cont_germ_df_2 <-seedlings_obs_germ %>% 
  group_by(cohort) %>%
  filter(cohort == "2") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "330",
                                     "IR" = "545",
                                     "RO" = "115")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

precip_cont_germ_df_3 <-seedlings_obs_germ %>% 
  group_by(cohort) %>%
  filter(cohort == "3") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "292",
                                     "IR" = "482",
                                     "RO" = "102")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

precip_cont_germ_df <- rbind(precip_cont_germ_df_1, precip_cont_germ_df_2, precip_cont_germ_df_3)

precip_cont_germ_df <- precip_cont_germ_df %>% mutate(precip_cont = as.factor(precip_cont))
                                                      
# germination model
zi.srer.germ.cont <- glmmTMB(tot_germination ~ precip_cont + excl + (1|cohort) + (1|sampID),# + ar1(date + 0|cohort),
                              data = precip_cont_germ_df,
                              family = binomial(link = "logit"))
zi.srer.germ.cont.sum <- summary(zi.srer.germ.cont)
zi.srer.germ.cont.sum

# get predictions of model
mydf_germ <- ggpredict(zi.srer.germ.cont, type = "re", terms = c("precip_cont", "excl", "cohort"))

# create graph
ggeff_germ_ppt_fig <- as.data.frame(mydf_germ) %>%
  mutate(excl = group,
         cohort = facet,
         precip = as.integer(as.character(x))) %>%
  mutate(excl = recode_factor(excl, 
                              "Control" = "None",
                              "Ants" = "Ants Excl",
                              "Rodents" = "Rodents Excl",
                              "Total" = "All Excl")) %>% 
  ggplot(aes(x = precip, y = predicted*100)) +
  geom_point(aes(color = excl)) +
  geom_pointrange(aes(ymin = ((100*predicted) - std.error), ymax = ((100*predicted) + std.error), color = excl), size = 0.5) +
  geom_smooth(method = "glm", formula = y ~ log(x) + x, se = F, size = 2)+
  labs(y = "Seed Germination (%)",
       x = "Precipitation (mm)",
       color = "Exclusion") +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550), limits = c(0, 600))+
  ylim(NA, 100) +
  theme_pubr(legend = "right")+
  labs_pubr(base_size = 24)

ggeff_germ_ppt_fig

ggsave(filename = "Figures_Tables/pred_germ_cont.tiff",
       plot = ggeff_germ_ppt_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# by exclusion tx
ggeff_excl_ppt_fig <- as.data.frame(mydf_germ) %>%
  mutate(excl = group,
         cohort = facet,
         precip = as.integer(as.character(x))) %>%
  mutate(excl = recode_factor(excl, 
                              "Control" = "None",
                              "Ants" = "Ants Excl",
                              "Rodents" = "Rodents Excl",
                              "Total" = "All Excl")) %>%  
  ggplot(aes(x = precip, y = predicted*100, group = excl, color = excl)) +
  #geom_point() +
  #geom_pointrange(aes(ymin = 10*lower, ymax = 10*upper, color = excl), size = 0.5) +
  geom_smooth(method = "glm", formula = y ~ log(x)+x, se = F, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Precipitation (mm)",
       color = "Exclusion") +
  scale_x_continuous(breaks = c(0,50, 100,150, 200,250, 300,350, 400,450, 500, 550), limits = c(0, 550))+
  ylim(0, 100) +
  theme_pubr(legend = "right")+
  labs_pubr(base_size = 24)

ggeff_excl_ppt_fig




