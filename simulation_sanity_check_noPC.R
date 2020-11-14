
# Average # of competitors by t_LOE:
MEPS_noPIV_mean <- MEPS_noPIV %>%
  filter(year > 2012) %>%
  group_by(t_LOE) %>%
  summarise(mean_competitor = mean(competitor),
            mean_P_b = mean(P_b, na.rm = T),
            mean_P_g = mean(P_g, na.rm = T)) %>%
  mutate(simulation = 0)


# With PC only
# Varying increase in hazard
model_noPIV_PWPGT_post <- coxph(Surv(genericnoPIV_postGDUFA$gaptime_start, genericnoPIV_postGDUFA$gaptime, entry1) ~ strata(ncompetitor) + route + AG + ETASU + guidance_before + indexyear + cluster(index), method = "breslow", data = genericnoPIV_postGDUFA)
summary(model_noPIV_PWPGT_post)

genericnoPIV_postGDUFA$predicted_risk <- predict(model_noPIV_PWPGT_post, type = "risk")

test <- genericnoPIV_postGDUFA %>%
  dplyr::select(index, order, gaptime, predicted_risk) 

#hazard_sim <- seq(2, 2, 1)
hazard_sim <- 1
n_simulation <- 1000

h_noPIV <- basehaz(model_noPIV_PWPGT_post, centered = T)
h_noPIV <- h_noPIV %>% 
  mutate(h2 = hazard*hazard_sim)

h_noPIV_1 <- h_noPIV %>%
  filter(strata == "ncompetitor=0")

h_noPIV_2 <- h_noPIV %>%
  filter(strata == "ncompetitor=1")

h_noPIV_3 <- h_noPIV %>%
  filter(strata == "ncompetitor=2")

E <- rep(NA, n_simulation)
E_simulated <- rep(NA, n_simulation)

MEPS_noPIV_1 <- MEPS_noPIV %>%
  filter(competitor >= 1,
         !is.na(P_g),
         !is.na(P_b))

MEPS_noPIV_id <- MEPS_noPIV_1 %>%
  distinct(index)

genericnoPIV_MEPS <- genericnoPIV_postGDUFA %>%
  inner_join(MEPS_noPIV_id, by = "index")

genericnoPIV_MEPS_1 <- genericnoPIV_MEPS %>%
  filter(ncompetitor == 0,
         gaptime > 0.001,
         entry1 == 1)

genericnoPIV_MEPS_2 <- genericnoPIV_MEPS %>%
  filter(ncompetitor == 1,
         gaptime > 0.001,
         entry1 == 1)

genericnoPIV_MEPS_3 <- genericnoPIV_MEPS %>%
  filter(ncompetitor == 2,
         gaptime > 0.001,
         entry1 == 1)

fit_noPIV_MEPS_1 <- survfit(Surv(gaptime, entry1) ~ 1,
                            data = genericnoPIV_MEPS_1)

ggsurvplot(fit_noPIV_MEPS_1)

fit_noPIV_MEPS_2 <- survfit(Surv(gaptime, entry1) ~ 1,
                          data = genericnoPIV_MEPS_2)

ggsurvplot(fit_noPIV_MEPS_2)

fit_noPIV_MEPS_3 <- survfit(Surv(gaptime, entry1) ~ 1,
                          data = genericnoPIV_MEPS_3)

ggsurvplot(fit_noPIV_MEPS_3)

genericnoPIV_MEPS_1_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_1, type="risk",se.fit=TRUE, reference = "strata")
genericnoPIV_MEPS_1$predicted <- genericnoPIV_MEPS_1_predicted$fit
genericnoPIV_MEPS_1$runif <- runif(nrow(genericnoPIV_MEPS_1), min=0, max=1)
genericnoPIV_MEPS_1$predicted_h <- - (log(genericnoPIV_MEPS_1$runif) / exp(genericnoPIV_MEPS_1$predicted + log(hazard_sim)))

genericnoPIV_MEPS_1$predicted_t <- sapply(genericnoPIV_MEPS_1$predicted_h, get_time_H0, data = h_noPIV_1)

summary(genericnoPIV_MEPS_1$predicted_t)

summary(genericnoPIV_MEPS_1$gaptime)

genericnoPIV_MEPS_2_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_2, type="risk",se.fit=TRUE, reference = "strata")
genericnoPIV_MEPS_2$predicted <- genericnoPIV_MEPS_2_predicted$fit
genericnoPIV_MEPS_2$runif <- runif(nrow(genericnoPIV_MEPS_2), min=0, max=1)
genericnoPIV_MEPS_2$predicted_h <- - (log(genericnoPIV_MEPS_2$runif) / exp(genericnoPIV_MEPS_2$predicted + log(hazard_sim)))

genericnoPIV_MEPS_2$predicted_t <- sapply(genericnoPIV_MEPS_2$predicted_h, get_time_H0, data = h_noPIV_2)

summary(genericnoPIV_MEPS_2$predicted_t)

summary(genericnoPIV_MEPS_2$gaptime)

genericnoPIV_MEPS_3_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_3, type="risk",se.fit=TRUE, reference = "strata")
genericnoPIV_MEPS_3$predicted <- genericnoPIV_MEPS_3_predicted$fit
genericnoPIV_MEPS_3$runif <- runif(nrow(genericnoPIV_MEPS_3), min=0, max=1)
genericnoPIV_MEPS_3$predicted_h <- - (log(genericnoPIV_MEPS_3$runif) / exp(genericnoPIV_MEPS_3$predicted + log(hazard_sim)))

genericnoPIV_MEPS_3$predicted_t <- sapply(genericnoPIV_MEPS_3$predicted_h, get_time_H0, data = h_noPIV_3)

genericnoPIV_MEPS %>%
  group_by(ncompetitor) %>%
  #filter(gaptime > 0.001) %>%
  summarise(mean = mean(gaptime),
            median = median(gaptime),
            min = min(gaptime),
            max = max(gaptime),
            n = n())

genericnoPIV_postGDUFA %>%
  group_by(ncompetitor) %>%
  summarise(mean = mean(gaptime),
            median = median(gaptime),
            max = max(gaptime),
            n = n())

summary(genericnoPIV_MEPS_3$predicted_t)

summary(genericnoPIV_MEPS_3$gaptime)

genericnoPIV_MEPS_1_simulated <- genericnoPIV_MEPS_1 %>%
  dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime) %>%
  mutate(sim = 0)

genericnoPIV_MEPS_2_simulated <- genericnoPIV_MEPS_2 %>%
  dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime) %>%
  mutate(sim = 0)

genericnoPIV_MEPS_3_simulated <- genericnoPIV_MEPS_3 %>%
  dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime) %>%
  mutate(sim = 0)

for (i in 1:n_simulation){
  MEPS_noPIV_1 <- MEPS_noPIV %>%
    filter(competitor >= 1,
           !is.na(P_g),
           !is.na(P_b))
  
  MEPS_noPIV_id <- MEPS_noPIV_1 %>%
    distinct(index)
  
  genericnoPIV_MEPS <- genericnoPIV_postGDUFA %>%
    inner_join(MEPS_noPIV_id, by = "index")
  
  genericnoPIV_MEPS_1 <- genericnoPIV_MEPS %>%
    filter(ncompetitor == 0,
           gaptime > 0.001,
           entry1 == 1)
  
  genericnoPIV_MEPS_2 <- genericnoPIV_MEPS %>%
    filter(ncompetitor == 1,
           gaptime > 0.001,
           entry1 == 1)
  
  genericnoPIV_MEPS_3 <- genericnoPIV_MEPS %>%
    filter(ncompetitor == 2,
           gaptime > 0.001,
           entry1 == 1)
  
  enericnoPIV_MEPS_1_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_1, type="risk",se.fit=TRUE, reference = "strata")
  genericnoPIV_MEPS_1$predicted <- genericnoPIV_MEPS_1_predicted$fit
  genericnoPIV_MEPS_1$runif <- runif(nrow(genericnoPIV_MEPS_1), min=0, max=1)
  genericnoPIV_MEPS_1$predicted_h <- - (log(genericnoPIV_MEPS_1$runif) / exp(genericnoPIV_MEPS_1$predicted + log(hazard_sim)))
  genericnoPIV_MEPS_1$predicted_t <- sapply(genericnoPIV_MEPS_1$predicted_h, get_time_H0, data = h_noPIV_1)
  
  
  genericnoPIV_MEPS_2_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_2, type="risk",se.fit=TRUE, reference = "strata")
  genericnoPIV_MEPS_2$predicted <- genericnoPIV_MEPS_2_predicted$fit
  genericnoPIV_MEPS_2$runif <- runif(nrow(genericnoPIV_MEPS_2), min=0, max=1)
  genericnoPIV_MEPS_2$predicted_h <- - (log(genericnoPIV_MEPS_2$runif) / exp(genericnoPIV_MEPS_2$predicted + log(hazard_sim)))
  genericnoPIV_MEPS_2$predicted_t <- sapply(genericnoPIV_MEPS_2$predicted_h, get_time_H0, data = h_noPIV_2)
  
  
  genericnoPIV_MEPS_3_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_3, type="risk",se.fit=TRUE, reference = "strata")
  genericnoPIV_MEPS_3$predicted <- genericnoPIV_MEPS_3_predicted$fit
  genericnoPIV_MEPS_3$runif <- runif(nrow(genericnoPIV_MEPS_3), min=0, max=1)
  genericnoPIV_MEPS_3$predicted_h <- - (log(genericnoPIV_MEPS_3$runif) / exp(genericnoPIV_MEPS_3$predicted + log(hazard_sim)))
  
  genericnoPIV_MEPS_3$predicted_t <- sapply(genericnoPIV_MEPS_3$predicted_h, get_time_H0, data = h_noPIV_3)

  genericnoPIV_MEPS_1 <- genericnoPIV_MEPS_1 %>%
    dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)
  
  genericnoPIV_MEPS_2 <- genericnoPIV_MEPS_2 %>%
    dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)
  
  genericnoPIV_MEPS_3 <- genericnoPIV_MEPS_3 %>%
    dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)
  
  genericnoPIV_MEPS_simulated <-genericnoPIV_MEPS_1 %>%
    bind_rows(genericnoPIV_MEPS_2) %>%
    bind_rows(genericnoPIV_MEPS_3)
  
  genericnoPIV_MEPS_1_simulated <- genericnoPIV_MEPS_1 %>%
    dplyr::select(-gaptime) %>%
    rename(gaptime = predicted_t) %>%
    mutate(sim = i) %>%
    bind_rows(genericnoPIV_MEPS_1_simulated)
  
  genericnoPIV_MEPS_2_simulated <- genericnoPIV_MEPS_2 %>%
    dplyr::select(-gaptime) %>%
    rename(gaptime = predicted_t) %>%
    mutate(sim = i) %>%
    bind_rows(genericnoPIV_MEPS_2_simulated)
  
  genericnoPIV_MEPS_3_simulated <- genericnoPIV_MEPS_3 %>%
    dplyr::select(-gaptime) %>%
    rename(gaptime = predicted_t) %>%
    mutate(sim = i) %>%
    bind_rows(genericnoPIV_MEPS_3_simulated)
  
  # update new gap time for second and third entrants only (with simulated)
  genericnoPIV_MEPS <- genericnoPIV_MEPS %>%
    left_join(genericnoPIV_MEPS_simulated, by = c("index", "Appl_No", "Product_No", "Strength", "gaptime", "order")) %>%
    distinct() %>%
    mutate(simulated_t = ifelse(is.na(predicted_t), gaptime, predicted_t))
  
  # get simulated approval date and approval year
  genericnoPIV_MEPS <- genericnoPIV_MEPS %>% 
    group_by(index) %>% 
    mutate(simulated_t_total = cumsum(simulated_t),
           simulated_approval_date = as.Date(exclusivity) + simulated_t,
           simulated_approval_year = year(simulated_approval_date))
}

genericnoPIV_MEPS_1_simulated %>%
  filter(sim != 0) %>%
  ungroup() %>%
  summarise(mean_gaptime = mean(gaptime),
            median_gaptime = median(gaptime)) %>%
  summarise(mean = mean(mean_gaptime),
            median = median(median_gaptime))

genericnoPIV_MEPS_2_simulated %>%
  filter(sim != 0) %>%
  ungroup() %>%
  summarise(mean_gaptime = mean(gaptime),
            median_gaptime = median(gaptime)) %>%
  summarise(mean = mean(mean_gaptime),
            median = median(median_gaptime))

genericnoPIV_MEPS_3_simulated %>%
  filter(sim != 0) %>%
  ungroup() %>%
  summarise(mean = mean(gaptime),
            median = median(gaptime))

genericnoPIV_MEPS_3_simulated %>%
  filter(sim != 0) %>%
  ungroup() %>%
  summarise(mean_gaptime = mean(gaptime),
            median_gaptime = median(gaptime)) %>%
  summarise(mean = mean(mean_gaptime),
            median = median(median_gaptime))


# plot simulated survival

genericnoPIV_MEPS_1_simulated <- genericnoPIV_MEPS_1_simulated %>%
  mutate(entry1 = 1)

genericnoPIV_MEPS_1_simulated <- genericnoPIV_MEPS_1_simulated %>%
  group_by(sim) %>%
  mutate(sim_rank = rank(gaptime, ties.method = "random"))

genericnoPIV_MEPS_1_simulated_sim025 <- genericnoPIV_MEPS_1_simulated %>%
  filter(sim > 0) %>%
  group_by(sim_rank) %>%
  summarise(gaptime = quantile(gaptime, 0.025)) %>%
  mutate(group = "2.5p") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_1_simulated_sim975 <- genericnoPIV_MEPS_1_simulated %>%
  filter(sim > 0) %>%
  group_by(sim_rank) %>%
  summarise(gaptime = quantile(gaptime, 0.975)) %>%
  mutate(group = "97.5p") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_1_simulated_origin <- genericnoPIV_MEPS_1_simulated %>%
  filter(sim == 0) %>%
  ungroup() %>%
  mutate(group = "data") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_1_simulated_95ci <- genericnoPIV_MEPS_1_simulated_origin %>%
  rbind(genericnoPIV_MEPS_1_simulated_sim025) %>%
  rbind(genericnoPIV_MEPS_1_simulated_sim975) %>%
  mutate(entry1 = 1)

fit_genericnoPIV_MEPS_1_simulated_95ci <- survfit(Surv(gaptime, entry1) ~ group,
                                                  data = genericnoPIV_MEPS_1_simulated_95ci)

ggsurvplot(fit_genericnoPIV_MEPS_1_simulated_95ci)



genericnoPIV_MEPS_2_simulated <- genericnoPIV_MEPS_2_simulated %>%
  mutate(entry1 = 1)

genericnoPIV_MEPS_2_simulated <- genericnoPIV_MEPS_2_simulated %>%
  group_by(sim) %>%
  mutate(sim_rank = rank(gaptime, ties.method = "random"))

genericnoPIV_MEPS_2_simulated_sim025 <- genericnoPIV_MEPS_2_simulated %>%
  filter(sim > 0) %>%
  group_by(sim_rank) %>%
  summarise(gaptime = quantile(gaptime, 0.025)) %>%
  mutate(group = "2.5p") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_2_simulated_sim975 <- genericnoPIV_MEPS_2_simulated %>%
  filter(sim > 0) %>%
  group_by(sim_rank) %>%
  summarise(gaptime = quantile(gaptime, 0.975)) %>%
  mutate(group = "97.5p") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_2_simulated_origin <- genericnoPIV_MEPS_2_simulated %>%
  filter(sim == 0) %>%
  ungroup() %>%
  mutate(group = "data") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_2_simulated_95ci <- genericnoPIV_MEPS_2_simulated_origin %>%
  rbind(genericnoPIV_MEPS_2_simulated_sim025) %>%
  rbind(genericnoPIV_MEPS_2_simulated_sim975) %>%
  mutate(entry1 = 1)

fit_genericnoPIV_MEPS_2_simulated_95ci <- survfit(Surv(gaptime, entry1) ~ group,
                                                data = genericnoPIV_MEPS_2_simulated_95ci)

ggsurvplot(fit_genericnoPIV_MEPS_2_simulated_95ci)


genericnoPIV_MEPS_3_simulated <- genericnoPIV_MEPS_3_simulated %>%
  mutate(entry1 = 1)

genericnoPIV_MEPS_3_simulated <- genericnoPIV_MEPS_3_simulated %>%
  group_by(sim) %>%
  mutate(sim_rank = rank(gaptime, ties.method = "random"))

genericnoPIV_MEPS_3_simulated_sim025 <- genericnoPIV_MEPS_3_simulated %>%
  filter(sim > 0) %>%
  group_by(sim_rank) %>%
  summarise(gaptime = quantile(gaptime, 0.025)) %>%
  mutate(group = "2.5p") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_3_simulated_sim975 <- genericnoPIV_MEPS_3_simulated %>%
  filter(sim > 0) %>%
  group_by(sim_rank) %>%
  summarise(gaptime = quantile(gaptime, 0.975)) %>%
  mutate(group = "97.5p") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_3_simulated_origin <- genericnoPIV_MEPS_3_simulated %>%
  filter(sim == 0) %>%
  ungroup() %>%
  mutate(group = "data") %>%
  dplyr::select(group, gaptime)

genericnoPIV_MEPS_3_simulated_95ci <- genericnoPIV_MEPS_3_simulated_origin %>%
  rbind(genericnoPIV_MEPS_3_simulated_sim025) %>%
  rbind(genericnoPIV_MEPS_3_simulated_sim975) %>%
  mutate(entry1 = 1)

fit_genericnoPIV_MEPS_3_simulated_95ci <- survfit(Surv(gaptime, entry1) ~ group,
                                                  data = genericnoPIV_MEPS_3_simulated_95ci)

ggsurvplot(fit_genericnoPIV_MEPS_3_simulated_95ci)

