
# No PC only
# Varying increase in hazard
model_noPIV_PWPGT_post <- coxph(Surv(genericnoPIV_postGDUFA$gaptime_start, genericnoPIV_postGDUFA$gaptime, entry1) ~ strata(ncompetitor) + route + AG + ETASU + guidance_before + indexyear + cluster(index), method = "breslow", data = genericnoPIV_postGDUFA)
summary(model_noPIV_PWPGT_post)

#hazard_sim <- seq(1, 5, 0.1)
hazard_sim <- seq(1, 10, 1)
E_diff_mean <- rep(NA, length(hazard_sim))
E_diff_025 <- rep(NA, length(hazard_sim))
E_diff_975 <- rep(NA, length(hazard_sim))
E_ratio_mean <- rep(NA, length(hazard_sim))
E_ratio_025 <- rep(NA, length(hazard_sim))
E_ratio_975 <- rep(NA, length(hazard_sim))
E_ratio_median <- rep(NA, length(hazard_sim))

n_simulation = 100

for (j in 1:length(hazard_sim)){
  
  h_noPIV <- basehaz(model_noPIV_PWPGT_post, centered = T)
  h_noPIV <- h_noPIV %>% 
    mutate(h2 = hazard*hazard_sim[j])
  
  h_noPIV_1 <- h_noPIV %>%
    filter(strata == "ncompetitor=0")
  
  h_noPIV_2 <- h_noPIV %>%
    filter(strata == "ncompetitor=1")
  
  h_noPIV_3 <- h_noPIV %>%
    filter(strata == "ncompetitor=2")
  
  E <- rep(NA, n_simulation)
  E_simulated <- rep(NA, n_simulation)
  
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
      filter(order == 1,
             gaptime > 0.001,
             entry1 == 1)
    
    genericnoPIV_MEPS_2 <- genericnoPIV_MEPS %>%
      filter(order == 2,
             gaptime > 0.001,
             entry1 == 1)
    
    genericnoPIV_MEPS_3 <- genericnoPIV_MEPS %>%
      filter(order == 3,
             gaptime > 0.001,
             entry1 == 1)
    
    genericnoPIV_MEPS_1_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_1, type="lp",se.fit=TRUE, reference = "strata")
    genericnoPIV_MEPS_1$predicted <- genericnoPIV_MEPS_1_predicted$fit
    genericnoPIV_MEPS_1$runif <- runif(nrow(genericnoPIV_MEPS_1), min=0, max=1)
    genericnoPIV_MEPS_1$predicted_h <- - (log(genericnoPIV_MEPS_1$runif) / exp(genericnoPIV_MEPS_1$predicted + log(hazard_sim[i])))
    
    genericnoPIV_MEPS_1$predicted_t <- sapply(genericnoPIV_MEPS_1$predicted_h, get_time_H0, data = h_noPIV_1)
    
    genericnoPIV_MEPS_2_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_2, type="lp",se.fit=TRUE, reference = "strata")
    genericnoPIV_MEPS_2$predicted <- genericnoPIV_MEPS_2_predicted$fit
    genericnoPIV_MEPS_2$runif <- runif(nrow(genericnoPIV_MEPS_2), min=0, max=1)
    genericnoPIV_MEPS_2$predicted_h <- - (log(genericnoPIV_MEPS_2$runif) / exp(genericnoPIV_MEPS_2$predicted + log(hazard_sim[i])))
    genericnoPIV_MEPS_2$predicted_t <- sapply(genericnoPIV_MEPS_2$predicted_h, get_time_H0, data = h_noPIV_2)
    
    genericnoPIV_MEPS_3_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_3, type="lp",se.fit=TRUE, reference = "strata")
    genericnoPIV_MEPS_3$predicted <- genericnoPIV_MEPS_3_predicted$fit
    genericnoPIV_MEPS_3$runif <- runif(nrow(genericnoPIV_MEPS_3), min=0, max=1)
    genericnoPIV_MEPS_3$predicted_h <- - (log(genericnoPIV_MEPS_3$runif) / exp(genericnoPIV_MEPS_3$predicted + log(hazard_sim[i])))
    
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
    
    # update new gap time for second and third entrants only (with simulated)
    genericnoPIV_MEPS <- genericnoPIV_MEPS %>%
      left_join(genericnoPIV_MEPS_simulated, by = c("index", "Appl_No", "Product_No", "Strength", "gaptime", "order")) %>%
      distinct() %>%
      mutate(simulated_t = ifelse(is.na(predicted_t), gaptime, predicted_t))
    
    # get simulated approval date and approval year
    genericnoPIV_MEPS <- genericnoPIV_MEPS %>% 
      group_by(index) %>% 
      mutate(simulated_t_total = cumsum(simulated_t),
             simulated_approval_date = as.Date(min) + simulated_t,
             simulated_approval_year = year(simulated_approval_date))
    
    genericnoPIV_MEPS_compare_ncompetitor <- genericnoPIV_MEPS %>% 
      dplyr::select(approveyear, simulated_approval_year) 
    
    competitor <- rep(0, 11)
    
    index <- (genericnoPIV_MEPS %>%
                distinct(index))$index
    
    MEPS_noPIV_simulated <- map_dfr(index, function(id){
      year <- c(2012:2017)
      competitor <- rep(0, 6)
      
      genericnoPIV_MEPS <- genericnoPIV_MEPS %>%
        filter(index == id,
               entry1 == 1)
      
      for (j in 2012:2017){
        for (i in 1:nrow(genericnoPIV_MEPS)){
          if(genericnoPIV_MEPS$simulated_approval_year[i] == j){
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
    
    MEPS_noPIV_simulated <- MEPS_noPIV_simulated %>%
      filter(year > 2012) %>%
      mutate(year = as.numeric(year),
             competitor_simulated = competitor_simulated + 1) %>%
      as_tibble()
    
    MEPS_noPIV_simulated_v1 <- MEPS_noPIV %>%
      left_join(MEPS_noPIV_simulated, by = c("index", "year"))
    
    # keep drugs when simulated result is available
    MEPS_noPIV_simulated_v1 <- MEPS_noPIV_simulated_v1 %>%
      filter(!is.na(competitor_simulated)) %>%
      mutate(competitor_diff = competitor_simulated - competitor - 1)
    
    MEPS_noPIV_simulated_v1 <- MEPS_noPIV_simulated_v1 %>%
      filter(!is.na(P_g) & !is.na(P_b)) %>% ## comment this if fill in previous years!!!
      mutate(P_b_simulated = P_b * (1 + rnorm(1, 0.0097, 0.005) * competitor_diff),
             P_g_simulated = P_g * (1 + rnorm(1, -0.077, 0.025) * competitor_diff),
             N_g_simulated = N_g * (1 + -0.16 * (P_g_simulated - P_g)/P_g), 
             E = P_b * N_b + P_g * N_g,
             E_simulated = P_b_simulated * (N_b + N_g - N_g_simulated) + P_g_simulated * N_g_simulated)
    
    #MEPS_noPIV_simulated_v1 <- MEPS_noPIV_simulated_v1 %>%
    #  filter(!is.na(P_g) & !is.na(P_b)) %>% ## comment this if fill in previous years!!!
    #  mutate(P_b_simulated = P_b * (1 + 0.01 * competitor_diff),
    #         P_g_simulated = P_g * (1 - 0.08 * competitor_diff),
    #         E = P_b * N_b + P_g * N_g,
    #         E_simulated = P_b_simulated * N_b + P_g_simulated * N_g)
    
    E[i] = sum(MEPS_noPIV_simulated_v1$E)
    E_simulated[i] = sum(MEPS_noPIV_simulated_v1$E_simulated)
    
  }
  
  E_df <- cbind(data.frame(E), data.frame(E_simulated))
  E_df <- E_df %>%
    mutate(E_diff = E - E_simulated,
           E_ratio = E_diff/E)
  
  E_diff_mean[j] <- mean(E_df$E_diff)
  E_diff_025[j] <- unname(quantile(E_df$E_diff, 0.025))
  E_diff_975[j] <- unname(quantile(E_df$E_diff, 0.975))
  
  E_ratio_mean[j] <- mean(E_df$E_ratio)
  E_ratio_median[j] <- median(E_df$E_ratio)
  E_ratio_025[j] <- unname(quantile(E_df$E_ratio, 0.025))
  E_ratio_975[j] <- unname(quantile(E_df$E_ratio, 0.975))
}

E_overall <- cbind(data.frame(hazard_sim), data.frame(E_diff_mean), data.frame(E_diff_025), data.frame(E_diff_975), data.frame(E_ratio_mean), data.frame(E_ratio_025), data.frame(E_ratio_975))

E_overall <- cbind(data.frame(hazard_sim), data.frame(E_diff_mean), data.frame(E_ratio_mean), data.frame(E_ratio_median))

E_noPC <- E_overall 

ggplot(data=E_noPC, aes(x=hazard_sim, y=E_ratio_median)) +
  geom_line()+
  geom_point() +
  ylim(c(0, 0.01))
