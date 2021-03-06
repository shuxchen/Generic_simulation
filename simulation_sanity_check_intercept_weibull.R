
# Average # of competitors by t_LOE:
MEPS_PIV_mean <- MEPS_PIV %>%
  filter(year > 2012) %>%
  group_by(t_LOE) %>%
  summarise(mean_competitor = mean(competitor),
            mean_P_b = mean(P_b, na.rm = T),
            mean_P_g = mean(P_g, na.rm = T)) %>%
  mutate(simulation = 0)


# With PC only
genericPIV_postGDUFA_second <- genericPIV_postGDUFA %>%
  filter(ncompetitor == 1,
         gaptime > 0.001)

genericPIV_postGDUFA_third <- genericPIV_postGDUFA %>%
  filter(ncompetitor == 2,
         gaptime > 0.001)

model_PIV_weibull_second <- flexsurvreg(Surv(gaptime, entry2) ~ AG + guidance_before + indexyear, data = genericPIV_postGDUFA_second, dist = "weibull")

model_PIV_weibull_second$res[-(1:2),"est"]
model_PIV_weibull_second$res

test <- genericPIV_postGDUFA %>%
  dplyr::select(index, order, gaptime, predicted_risk) 

#hazard_sim <- seq(2, 2, 1)
hazard_sim <- 1
n_simulation <- 100

h_PIV <- basehaz(model_PIV_PWPGT_post, centered = T)
#h_PIV <- h_PIV %>% 
#  mutate(h2 = hazard*hazard_sim)

h_PIV_2 <- h_PIV %>%
  filter(strata == "ncompetitor=1")

h_PIV_3 <- h_PIV %>%
  filter(strata == "ncompetitor=2")

E <- rep(NA, n_simulation)
E_simulated <- rep(NA, n_simulation)

MEPS_PIV_2 <- MEPS_PIV %>%
  filter(competitor >= 2,
         !is.na(P_g),
         !is.na(P_b))

MEPS_PIV_id <- MEPS_PIV_2 %>%
  distinct(index)

genericPIV_MEPS <- genericPIV_postGDUFA %>%
  inner_join(MEPS_PIV_id, by = "index")

genericPIV_MEPS_2 <- genericPIV_MEPS %>%
  filter(ncompetitor == 1,
         gaptime > 0.001)

genericPIV_MEPS_3 <- genericPIV_MEPS %>%
  filter(ncompetitor == 2,
         gaptime > 0.001)

genericPIV_MEPS_2_predicted <- predict(model_PIV_PWPGT_post, genericPIV_MEPS_2, type="lp",se.fit=TRUE, reference = "strata")
genericPIV_MEPS_2$predicted <- genericPIV_MEPS_2_predicted$fit
genericPIV_MEPS_2$runif <- runif(nrow(genericPIV_MEPS_2), min=0, max=1)
genericPIV_MEPS_2$predicted_h <- - (log(genericPIV_MEPS_2$runif) / exp(genericPIV_MEPS_2$predicted + log(hazard_sim)))
genericPIV_MEPS_2$predicted_t <- sapply(genericPIV_MEPS_2$predicted_h, get_time_H0, data = h_PIV_2)



genericPIV_MEPS_3_predicted <- predict(model_PIV_PWPGT_post, genericPIV_MEPS_3, type="lp",se.fit=TRUE, reference = "strata")
genericPIV_MEPS_3$predicted <- genericPIV_MEPS_3_predicted$fit
genericPIV_MEPS_3$runif <- runif(nrow(genericPIV_MEPS_3), min=0, max=1)
genericPIV_MEPS_3$predicted_h <- - (log(genericPIV_MEPS_3$runif) / exp(genericPIV_MEPS_3$predicted + log(hazard_sim)))

genericPIV_MEPS_3$predicted_t <- sapply(genericPIV_MEPS_3$predicted_h, get_time_H0, data = h_PIV_3)

genericPIV_MEPS_2 <- genericPIV_MEPS_2 %>%
  dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)

genericPIV_MEPS_3 <- genericPIV_MEPS_3 %>%
  dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)

genericPIV_MEPS_simulated <-genericPIV_MEPS_2 %>%
  bind_rows(genericPIV_MEPS_3)

genericPIV_MEPS_2_simulated <- genericPIV_MEPS_2 %>%
  dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime) %>%
  mutate(sim = 0)

genericPIV_MEPS_3_simulated <- genericPIV_MEPS_3 %>%
  dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime) %>%
  mutate(sim = 0)

for (i in 1:n_simulation){
  MEPS_PIV_2 <- MEPS_PIV %>%
    filter(competitor >= 2,
           !is.na(P_g),
           !is.na(P_b))
  
  MEPS_PIV_id <- MEPS_PIV_2 %>%
    distinct(index)
  
  genericPIV_MEPS <- genericPIV_postGDUFA %>%
    inner_join(MEPS_PIV_id, by = "index")
  
  genericPIV_MEPS_2 <- genericPIV_MEPS %>%
    filter(ncompetitor == 1,
           gaptime > 0.001)
  
  genericPIV_MEPS_3 <- genericPIV_MEPS %>%
    filter(ncompetitor == 2,
           gaptime > 0.001)
  
  genericPIV_MEPS_2_predicted <- predict(model_PIV_PWPGT_post, genericPIV_MEPS_2, type="lp",se.fit=TRUE, reference = "strata")
  genericPIV_MEPS_2$predicted <- genericPIV_MEPS_2_predicted$fit
  genericPIV_MEPS_2$runif <- runif(nrow(genericPIV_MEPS_2), min=0, max=1)
  genericPIV_MEPS_2$predicted_h <- - (log(genericPIV_MEPS_2$runif) / exp(genericPIV_MEPS_2$predicted + log(hazard_sim)))
  genericPIV_MEPS_2$predicted_t <- sapply(genericPIV_MEPS_2$predicted_h, get_time_H0, data = h_PIV_2)
  
  
  
  genericPIV_MEPS_3_predicted <- predict(model_PIV_PWPGT_post, genericPIV_MEPS_3, type="lp",se.fit=TRUE, reference = "strata")
  genericPIV_MEPS_3$predicted <- genericPIV_MEPS_3_predicted$fit
  genericPIV_MEPS_3$runif <- runif(nrow(genericPIV_MEPS_3), min=0, max=1)
  genericPIV_MEPS_3$predicted_h <- - (log(genericPIV_MEPS_3$runif) / exp(genericPIV_MEPS_3$predicted + log(hazard_sim)))
  
  genericPIV_MEPS_3$predicted_t <- sapply(genericPIV_MEPS_3$predicted_h, get_time_H0, data = h_PIV_3)
  
  genericPIV_MEPS_2 <- genericPIV_MEPS_2 %>%
    dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)
  
  genericPIV_MEPS_3 <- genericPIV_MEPS_3 %>%
    dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)
  
  genericPIV_MEPS_simulated <-genericPIV_MEPS_2 %>%
    bind_rows(genericPIV_MEPS_3)
  
  genericPIV_MEPS_2_simulated <- genericPIV_MEPS_2 %>%
    dplyr::select(-gaptime) %>%
    rename(gaptime = predicted_t) %>%
    mutate(sim = i) %>%
    bind_rows(genericPIV_MEPS_2_simulated)
  
  genericPIV_MEPS_3_simulated <- genericPIV_MEPS_3 %>%
    dplyr::select(-gaptime) %>%
    rename(gaptime = predicted_t) %>%
    mutate(sim = i) %>%
    bind_rows(genericPIV_MEPS_3_simulated)
  
  # update new gap time for second and third entrants only (with simulated)
  genericPIV_MEPS <- genericPIV_MEPS %>%
    left_join(genericPIV_MEPS_simulated, by = c("index", "Appl_No", "Product_No", "Strength", "gaptime", "order")) %>%
    distinct() %>%
    mutate(simulated_t = ifelse(is.na(predicted_t), gaptime, predicted_t))
  
  # get simulated approval date and approval year
  genericPIV_MEPS <- genericPIV_MEPS %>% 
    group_by(index) %>% 
    mutate(simulated_t_total = cumsum(simulated_t),
           simulated_approval_date = as.Date(exclusivity) + simulated_t,
           simulated_approval_year = year(simulated_approval_date))
  
  genericPIV_MEPS_compare_ncompetitor <- genericPIV_MEPS %>% 
    dplyr::select(approveyear, simulated_approval_year) 
  
  competitor <- rep(0, 11)
  
  index <- (genericPIV_MEPS %>%
              distinct(index))$index
  
  MEPS_PIV_simulated <- map_dfr(index, function(id){
    year <- c(2012:2017)
    competitor <- rep(0, 6)
    
    genericPIV_MEPS <- genericPIV_MEPS %>%
      filter(index == id,
             entry2 == 1)
    
    for (j in 2012:2017){
      for (i in 1:nrow(genericPIV_MEPS)){
        if(genericPIV_MEPS$simulated_approval_year[i] == j){
          competitor[j-2011]= competitor[j-2011] + 1
        }
      }
      if(j<2017){
        competitor[(j-2012):6]=competitor[j-2011]
      }
    }
    
    output <- data.frame(year, competitor)
    output <- output %>%
      mutate(index = id) %>%
      dplyr::select(index, year, competitor) %>%
      rename(competitor_simulated = competitor) %>%
      data.table()
    
    return(output)
  })
  
  MEPS_PIV_simulated <- MEPS_PIV_simulated %>%
    filter(year > 2012) %>%
    mutate(year = as.numeric(year),
           competitor_simulated = competitor_simulated + 1) %>%
    as_tibble()
  
  MEPS_PIV_simulated_v1 <- MEPS_PIV %>%
    left_join(MEPS_PIV_simulated, by = c("index", "year"))
  
  # keep drugs when simulated result is available
  MEPS_PIV_simulated_v1 <- MEPS_PIV_simulated_v1 %>%
    filter(!is.na(competitor_simulated)) %>%
    mutate(competitor_diff = competitor_simulated - competitor)
  
  MEPS_PIV_simulated_v1 <- MEPS_PIV_simulated_v1 %>%
    filter(!is.na(P_g) & !is.na(P_b)) %>% ## comment this if fill in previous years!!!
    mutate(P_b_simulated = P_b * (1 + rnorm(1, 0.0097, 0.005) * competitor_diff),
           P_g_simulated = P_g * (1 + rnorm(1, -0.077, 0.025) * competitor_diff),
           N_g_simulated = N_g * (1 + -0.16 * (P_g_simulated - P_g)/P_g), 
           E = P_b * N_b + P_g * N_g,
           E_simulated = P_b_simulated * (N_b + N_g - N_g_simulated) + P_g_simulated * N_g_simulated)
  
  E[i] = sum(MEPS_PIV_simulated_v1$E)
  E_simulated[i] = sum(MEPS_PIV_simulated_v1$E_simulated)
  
  simulation_mean <- MEPS_PIV_simulated_v1 %>%
    group_by(t_LOE) %>%
    summarise(mean_competitor = mean(competitor_simulated),
              mean_P_b = mean(P_b_simulated, na.rm = T),
              mean_P_g = mean(P_g_simulated, na.rm = T)) %>%
    mutate(simulation = i)
  
  MEPS_PIV_mean <- MEPS_PIV_mean %>%
    rbind(simulation_mean)
  
}

genericPIV_MEPS_2_simulated %>%
  filter(sim != 0) %>%
  ungroup() %>%
  summarise(mean_gaptime = mean(gaptime),
            median_gaptime = median(gaptime)) %>%
  summarise(mean = mean(mean_gaptime),
            median = median(median_gaptime))

genericPIV_MEPS_3_simulated %>%
  filter(sim != 0) %>%
  ungroup() %>%
  summarise(mean = mean(gaptime),
            median = median(gaptime))