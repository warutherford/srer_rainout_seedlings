

ppt_arrange <- ppt_srer_avg_all %>% 
  arrange(ppt_yr_mm)

quantile(ppt_arrange$ppt_yr_mm, c(.10, .90))

ppt_drought <- ppt_arrange %>% filter(ppt_yr_mm <= 290.3728) %>% 
  mutate(precip = "Drought")

ppt_wet <- ppt_arrange %>% filter(ppt_yr_mm >= 534.5684) %>% 
  mutate(precip = "Wet")

ppt_ambient <- ppt_arrange %>% filter(ppt_yr_mm > 290 & ppt_yr_mm < 534) %>% 
  mutate(precip = "Ambient")

ppt_srer_new <- rbind(ppt_ambient, ppt_drought, ppt_wet) %>% mutate(precip_cont = (precip), precip_cont = ppt_yr_mm) %>% 
  mutate(excl = NA)


nd <- ppt_srer_new %>%  mutate(cohort = as.factor(YEAR)) %>% dplyr::select(precip_cont, cohort) %>% 
  arrange(YEAR)

zi.srer.surv.cont.pred <- glmmTMB(survival/10 ~ as.numeric(as.character(precip_cont)) + log(as.numeric(as.character(precip_cont))) + (1|cohort),
                             data = precip_cont_surv_df,
                             family = gaussian())

summary(zi.srer.surv.cont.pred)

ppt_pred <- predict(zi.srer.surv.cont.pred, newdata = nd, allow.new.levels = T, re.form = NULL)

ppt_pred <- as.data.frame(ppt_pred)

srer_comb_ppt_surv <- cbind(ppt_pred, nd)

srer_mod_surv <- srer_comb_ppt_surv %>% 
  mutate(surv = 100*ppt_pred,
         year = YEAR,
         precip = precip_cont) %>% 
  dplyr::select(surv, year,precip)

srer_mod_surv %>% filter(precip <= 290) %>% 
  mutate(pptx = "Drought") %>% 
  group_by(pptx) %>% 
  summarise(surv_mean = mean(surv),
            sd_surv = sd(surv),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts)))

srer_mod_surv %>% filter(precip >= 534) %>% 
  mutate(pptx = "Wet") %>% 
  group_by(pptx) %>% 
  summarise(surv_mean = mean(surv),
            sd_surv = sd(surv),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts)))

srer_mod_surv %>% filter(precip > 290 & precip < 534) %>% 
  mutate(pptx = "Ambient")%>% 
  group_by(pptx) %>% 
  summarise(surv_mean = mean(surv),
            sd_surv = sd(surv),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts)))

srer_mod_surv %>% 
  ggplot(aes(x = precip, y = surv, alpha = year)) +
  geom_point()+
  ylim(0,30)
