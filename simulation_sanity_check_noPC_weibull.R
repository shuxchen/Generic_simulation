# Average # of competitors by t_LOE:
MEPS_noPIV_mean <- MEPS_noPIV %>%
  filter(year > 2012) %>%
  group_by(t_LOE) %>%
  summarise(mean_competitor = mean(competitor),
            mean_P_b = mean(P_b, na.rm = T),
            mean_P_g = mean(P_g, na.rm = T)) %>%
  mutate(simulation = 0)


genericnoPIV_postGDUFA_first <- genericnoPIV_postGDUFA %>%
  filter(ncompetitor == 0,
         gaptime > 0.001)

genericnoPIV_postGDUFA_second <- genericnoPIV_postGDUFA %>%
  filter(ncompetitor == 1,
         gaptime > 0.001)

genericnoPIV_postGDUFA_third <- genericnoPIV_postGDUFA %>%
  filter(ncompetitor == 2,
         gaptime > 0.001)

model_noPIV_first <- flexsurvreg(Surv(gaptime, entry1) ~ AG + guidance_before + indexyear, data = genericnoPIV_postGDUFA_first, dist = "weibullPH")
plot(model_noPIV_first, type = "cumhaz")

model_noPIV_second <- flexsurvreg(Surv(gaptime, entry1) ~ AG + guidance_before + indexyear, data = genericnoPIV_postGDUFA_second, dist = "weibullPH")

model_noPIV_third <- flexsurvreg(Surv(gaptime, entry1) ~ AG + guidance_before + indexyear, data = genericnoPIV_postGDUFA_third, dist = "weibullPH")

#hazard_sim <- seq(2, 2, 1)
hazard_sim <- 5
n_simulation <- 100



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
         gaptime > 0.001)

genericnoPIV_MEPS_2 <- genericnoPIV_MEPS %>%
  filter(ncompetitor == 1,
         gaptime > 0.001)

genericnoPIV_MEPS_3 <- genericnoPIV_MEPS %>%
  filter(ncompetitor == 2,
         gaptime > 0.001)

genericnoPIV_MEPS_1$predicted <- genericnoPIV_MEPS_1$AG*model_noPIV_first$res[3,"est"] + genericnoPIV_MEPS_1$guidance_before*model_noPIV_first$res[4,"est"] + genericnoPIV_MEPS_1$indexyear*model_noPIV_first$res[5,"est"]
genericnoPIV_MEPS_1$runif <- runif(nrow(genericnoPIV_MEPS_1), min=0, max=1)
genericnoPIV_MEPS_1$predicted_t <- (1/(model_noPIV_first$res[2,"est"])) * (-log(genericnoPIV_MEPS_1$runif) / exp(genericnoPIV_MEPS_1$predicted + log(hazard_sim)))^(1/model_noPIV_first$res[1,"est"])  
genericnoPIV_MEPS_1$predicted_h <- (-log(genericnoPIV_MEPS_1$runif) / exp(genericnoPIV_MEPS_1$predicted + log(hazard_sim))) 


genericnoPIV_MEPS_2$predicted <- genericnoPIV_MEPS_2$AG*model_noPIV_second$res[3,"est"] + genericnoPIV_MEPS_2$guidance_before*model_noPIV_second$res[4,"est"] + genericnoPIV_MEPS_2$indexyear*model_noPIV_second$res[5,"est"]
genericnoPIV_MEPS_2$runif <- runif(nrow(genericnoPIV_MEPS_2), min=0, max=1)
genericnoPIV_MEPS_2$predicted_t <- (1/(model_noPIV_second$res[2,"est"])) * (-log(genericnoPIV_MEPS_2$runif) / exp(genericnoPIV_MEPS_2$predicted + log(hazard_sim)))^(1/model_noPIV_second$res[1,"est"])  
genericnoPIV_MEPS_2$predicted_h <- (-log(genericnoPIV_MEPS_2$runif) / exp(genericnoPIV_MEPS_2$predicted + log(hazard_sim))) 

genericnoPIV_MEPS_3$predicted <- genericnoPIV_MEPS_3$AG*model_noPIV_third$res[3,"est"] + genericnoPIV_MEPS_3$guidance_before*model_noPIV_third$res[4,"est"] + genericnoPIV_MEPS_3$indexyear*model_noPIV_third$res[5,"est"]
genericnoPIV_MEPS_3$runif <- runif(nrow(genericnoPIV_MEPS_3), min=0, max=1)
genericnoPIV_MEPS_3$predicted_t <- (1/(model_noPIV_third$res[2,"est"])) * (-log(genericnoPIV_MEPS_3$runif) / exp(genericnoPIV_MEPS_3$predicted + log(hazard_sim)))^(1/model_noPIV_third$res[1,"est"])  

genericnonoPIV_MEPS %>%
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
           gaptime > 0.001)
  
  genericnoPIV_MEPS_2 <- genericnoPIV_MEPS %>%
    filter(ncompetitor == 1,
           gaptime > 0.001)
  
  genericnoPIV_MEPS_3 <- genericnoPIV_MEPS %>%
    filter(ncompetitor == 2,
           gaptime > 0.001)
  
  genericnoPIV_MEPS_1$predicted <- genericnoPIV_MEPS_1$AG*model_noPIV_first$res[3,"est"] + genericnoPIV_MEPS_1$guidance_before*model_noPIV_first$res[4,"est"] + genericnoPIV_MEPS_1$indexyear*model_noPIV_first$res[5,"est"]
  genericnoPIV_MEPS_1$runif <- runif(nrow(genericnoPIV_MEPS_1), min=0, max=1)
  genericnoPIV_MEPS_1$predicted_t <- (1/(model_noPIV_first$res[2,"est"])) * (-log(genericnoPIV_MEPS_1$runif) / exp(genericnoPIV_MEPS_1$predicted + log(hazard_sim)))^(1/model_noPIV_first$res[1,"est"])  
  genericnoPIV_MEPS_1$predicted_h <- (-log(genericnoPIV_MEPS_1$runif) / exp(genericnoPIV_MEPS_1$predicted + log(hazard_sim))) 
  
  
  genericnoPIV_MEPS_2$predicted <- genericnoPIV_MEPS_2$AG*model_noPIV_second$res[3,"est"] + genericnoPIV_MEPS_2$guidance_before*model_noPIV_second$res[4,"est"] + genericnoPIV_MEPS_2$indexyear*model_noPIV_second$res[5,"est"]
  genericnoPIV_MEPS_2$runif <- runif(nrow(genericnoPIV_MEPS_2), min=0, max=1)
  genericnoPIV_MEPS_2$predicted_t <- (1/(model_noPIV_second$res[2,"est"])) * (-log(genericnoPIV_MEPS_2$runif) / exp(genericnoPIV_MEPS_2$predicted + log(hazard_sim)))^(1/model_noPIV_second$res[1,"est"])  
  genericnoPIV_MEPS_2$predicted_h <- (-log(genericnoPIV_MEPS_2$runif) / exp(genericnoPIV_MEPS_2$predicted + log(hazard_sim))) 
  
  genericnoPIV_MEPS_3$predicted <- genericnoPIV_MEPS_3$AG*model_noPIV_third$res[3,"est"] + genericnoPIV_MEPS_3$guidance_before*model_noPIV_third$res[4,"est"] + genericnoPIV_MEPS_3$indexyear*model_noPIV_third$res[5,"est"]
  genericnoPIV_MEPS_3$runif <- runif(nrow(genericnoPIV_MEPS_3), min=0, max=1)
  genericnoPIV_MEPS_3$predicted_t <- (1/(model_noPIV_third$res[2,"est"])) * (-log(genericnoPIV_MEPS_3$runif) / exp(genericnoPIV_MEPS_3$predicted + log(hazard_sim)))^(1/model_noPIV_third$res[1,"est"])  
  
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