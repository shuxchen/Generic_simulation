
model_noPIV_PWPGT_post <- coxph(Surv(genericnoPIV_postGDUFA$gaptime_start, genericnoPIV_postGDUFA$gaptime, entry1) ~ strata(ncompetitor) + route + AG + ETASU + guidance_before + indexyear + cluster(index), method = "breslow", data = genericnoPIV_postGDUFA)
summary(model_noPIV_PWPGT_post)

model_PIV_PWPGT_post <- coxph(Surv(genericPIV_postGDUFA$gaptime_start, genericPIV_postGDUFA$gaptime, entry2) ~ strata(ncompetitor) + route + AG + ETASU + guidance_before + indexyear + cluster(index), method = "breslow", data = genericPIV_postGDUFA)
summary(model_PIV_PWPGT_post)

#hazard_sim <- seq(1, 5, 0.1)
hazard_sim <- seq(0.1, 5, 0.2)
E_diff_mean <- rep(NA, length(hazard_sim))
E_diff_025 <- rep(NA, length(hazard_sim))
E_diff_975 <- rep(NA, length(hazard_sim))
E_ratio_mean <- rep(NA, length(hazard_sim))
E_ratio_025 <- rep(NA, length(hazard_sim))
E_ratio_975 <- rep(NA, length(hazard_sim))
E_ratio_median <- rep(NA, length(hazard_sim))

n_simulation = 500

for (j in 1:length(hazard_sim)){
  
  h_PIV <- basehaz(model_PIV_PWPGT_post, centered = T)
  h_PIV <- h_PIV %>% 
    mutate(h2 = hazard*hazard_sim[j])
  
  h_PIV_2 <- h_PIV %>%
    filter(strata == "ncompetitor=1")
  
  h_PIV_3 <- h_PIV %>%
    filter(strata == "ncompetitor=2")
  
  
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
    
    #1. With PC
    MEPS_PIV_2 <- MEPS_PIV %>%
      filter(competitor >= 2,
             !is.na(P_g),
             !is.na(P_b))
    
    MEPS_PIV_id <- MEPS_PIV_2 %>%
      distinct(index)
    
    genericPIV_MEPS <- genericPIV_postGDUFA %>%
      filter(entry2 == 1) %>%
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
    genericPIV_MEPS_2$predicted_h <- - (log(genericPIV_MEPS_2$runif) / exp(genericPIV_MEPS_2$predicted + log(hazard_sim[j])))
    genericPIV_MEPS_2$predicted_t <- sapply(genericPIV_MEPS_2$predicted_h, get_time_H0, data = h_PIV_2)
    
    genericPIV_MEPS_3_predicted <- predict(model_PIV_PWPGT_post, genericPIV_MEPS_3, type="lp",se.fit=TRUE, reference = "strata")
    genericPIV_MEPS_3$predicted <- genericPIV_MEPS_3_predicted$fit
    genericPIV_MEPS_3$runif <- runif(nrow(genericPIV_MEPS_3), min=0, max=1)
    genericPIV_MEPS_3$predicted_h <- - (log(genericPIV_MEPS_3$runif) / exp(genericPIV_MEPS_3$predicted + log(hazard_sim[j])))
    genericPIV_MEPS_3$predicted_t <- sapply(genericPIV_MEPS_3$predicted_h, get_time_H0, data = h_PIV_3)
    
    genericPIV_MEPS_2 <- genericPIV_MEPS_2 %>%
      dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)
    
    genericPIV_MEPS_3 <- genericPIV_MEPS_3 %>%
      dplyr::select(index, Appl_No, Product_No, Strength, order, gaptime, predicted_t)
    
    genericPIV_MEPS_simulated <-genericPIV_MEPS_2 %>%
      bind_rows(genericPIV_MEPS_3)
    
    # update new gap time for second and third entrants only (with simulated)
    genericPIV_MEPS <- genericPIV_MEPS %>%
      left_join(genericPIV_MEPS_simulated, by = c("index", "Appl_No", "Product_No", "Strength", "gaptime", "order")) %>%
      distinct() %>%
      mutate(simulated_t = ifelse(is.na(predicted_t), gaptime, predicted_t))
    
    # get simulated approval date and approval year
    genericPIV_MEPS <- genericPIV_MEPS %>% 
      group_by(index) %>% 
      mutate(simulated_t_total = cumsum(simulated_t),
             simulated_approval_date = as.Date(exclusivity) + simulated_t_total,
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
      #mutate(competitor_diff = ifelse(year == 2017, 0, competitor_simulated - competitor))
      mutate(competitor_diff = competitor_simulated - competitor)
    
    MEPS_PIV_simulated_v1 <- MEPS_PIV_simulated_v1 %>%
      filter(!is.na(P_g) & !is.na(P_b)) %>% ## comment this if fill in previous years!!!
      mutate(#branded_change = rnorm(n(), exp(-0.04370)-1, exp(0.01100) - 1),
             #branded_change = rnorm(n(), exp(0.00099)-1, exp(0.00458) - 1),
             branded_change = rnorm(n(), exp(-0.02874)-1, exp(0.00970) - 1),
             #generic_change = rnorm(n(), exp(-0.06644)-1, exp(0.01007) - 1),
             #generic_change = rnorm(n(), exp(-0.06929)-1, exp(0.00897) - 1),
             generic_change = rnorm(n(), exp(-0.06728)-1, exp(0.01010) - 1),
             P_b_simulated = P_b * (1 + branded_change * competitor_diff),
             P_g_simulated = P_g * (1 + generic_change * competitor_diff),
             N_g_simulated = N_g * (1 + -0.16 * (P_g_simulated - P_g)/P_g), 
             E = P_b * N_b + P_g * N_g,
             E_simulated = P_b_simulated * (N_b + N_g - N_g_simulated) + P_g_simulated * N_g_simulated)
    
    
    #2. Without PC
    MEPS_noPIV_1 <- MEPS_noPIV %>%
      filter(competitor >= 1,
             !is.na(P_g),
             !is.na(P_b))
    
    MEPS_noPIV_id <- MEPS_noPIV_1 %>%
      distinct(index)
    
    genericnoPIV_MEPS <- genericnoPIV_postGDUFA %>%
      filter(entry1 == 1) %>%
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
    genericnoPIV_MEPS_1$predicted_h <- - (log(genericnoPIV_MEPS_1$runif) / exp(genericnoPIV_MEPS_1$predicted + log(hazard_sim[j])))
    
    genericnoPIV_MEPS_1$predicted_t <- sapply(genericnoPIV_MEPS_1$predicted_h, get_time_H0, data = h_noPIV_1)
    
    genericnoPIV_MEPS_2_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_2, type="lp",se.fit=TRUE, reference = "strata")
    genericnoPIV_MEPS_2$predicted <- genericnoPIV_MEPS_2_predicted$fit
    genericnoPIV_MEPS_2$runif <- runif(nrow(genericnoPIV_MEPS_2), min=0, max=1)
    genericnoPIV_MEPS_2$predicted_h <- - (log(genericnoPIV_MEPS_2$runif) / exp(genericnoPIV_MEPS_2$predicted + log(hazard_sim[j])))
    genericnoPIV_MEPS_2$predicted_t <- sapply(genericnoPIV_MEPS_2$predicted_h, get_time_H0, data = h_noPIV_2)
    
    genericnoPIV_MEPS_3_predicted <- predict(model_noPIV_PWPGT_post, genericnoPIV_MEPS_3, type="lp",se.fit=TRUE, reference = "strata")
    genericnoPIV_MEPS_3$predicted <- genericnoPIV_MEPS_3_predicted$fit
    genericnoPIV_MEPS_3$runif <- runif(nrow(genericnoPIV_MEPS_3), min=0, max=1)
    genericnoPIV_MEPS_3$predicted_h <- - (log(genericnoPIV_MEPS_3$runif) / exp(genericnoPIV_MEPS_3$predicted + log(hazard_sim[j])))
    
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
             simulated_approval_date = as.Date(min) + simulated_t_total,
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
             competitor_simulated = competitor_simulated) %>%
      as_tibble()
    
    MEPS_noPIV_simulated_v1 <- MEPS_noPIV %>%
      left_join(MEPS_noPIV_simulated, by = c("index", "year"))
    
    # keep drugs when simulated result is available
    MEPS_noPIV_simulated_v1 <- MEPS_noPIV_simulated_v1 %>%
      filter(!is.na(competitor_simulated)) %>%
      #mutate(competitor_diff = ifelse(year == 2017, 0, competitor_simulated - competitor))
      mutate(competitor_diff = competitor_simulated - competitor)

    MEPS_noPIV_simulated_v1 <- MEPS_noPIV_simulated_v1 %>%
      filter(!is.na(P_g) & !is.na(P_b)) %>% ## comment this if fill in previous years!!!
      mutate(#branded_change = rnorm(n(), exp(-0.04370)-1, exp(0.01100) - 1),
        #branded_change = rnorm(n(), exp(0.00099)-1, exp(0.00458) - 1),
        branded_change = rnorm(n(), exp(-0.02874)-1, exp(0.00970) - 1),
        #generic_change = rnorm(n(), exp(-0.06644)-1, exp(0.01007) - 1),
        #generic_change = rnorm(n(), exp(-0.06929)-1, exp(0.00897) - 1),
        generic_change = rnorm(n(), exp(-0.06728)-1, exp(0.01010) - 1),
        P_b_simulated = P_b * (1 + branded_change * competitor_diff),
        P_g_simulated = P_g * (1 + generic_change * competitor_diff),
        N_g_simulated = N_g * (1 + -0.16 * (P_g_simulated - P_g)/P_g), 
        E = P_b * N_b + P_g * N_g,
        E_simulated = P_b_simulated * (N_b + N_g - N_g_simulated) + P_g_simulated * N_g_simulated)
    
    #3. Combine both
    MEPS_simulated_v1 <- MEPS_noPIV_simulated_v1 %>%
      ungroup() %>%
      bind_rows(MEPS_PIV_simulated_v1 %>% ungroup())

    E[i] = sum(MEPS_simulated_v1$E)
    E_simulated[i] = sum(MEPS_simulated_v1$E_simulated)
    
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

E_both <- E_overall 

E_both$E_ratio_mean <- 100*E_both$E_ratio_mean
E_both$E_ratio_025 <- 100*E_both$E_ratio_025
E_both$E_ratio_975 <- 100*E_both$E_ratio_975

#E_both$hazard_sim <- as.factor(E_both$hazard_sim)

ggplot(data=E_both, aes(x=hazard_sim, y=E_ratio_mean, group = 1)) +
  geom_line(color = "#9966FF")+
  geom_point(color = "#9966FF") +
  geom_ribbon(aes(ymin=E_both$E_ratio_025, ymax=E_both$E_ratio_975), linetype=2, alpha=0.25, fill="blue") + 
  scale_colour_manual("",values="blue")+
  scale_fill_manual("",values="#CCCCFF") +
  ylim(c(-10, 25)) +
  xlab("Policy shock k") +
  ylab("Change in expenditure reduction (%)") 


E_both_mean <- E_both %>%
  dplyr::select(hazard_sim, E_ratio_mean) %>%
  rename(E_ratio = E_ratio_mean) %>%
  mutate(group = "Mean")

E_both_025 <- E_both %>%
  dplyr::select(hazard_sim, E_ratio_025) %>%
  rename(E_ratio = E_ratio_025) %>%
  mutate(group = "Lower bound")

E_both_975 <- E_both %>%
  dplyr::select(hazard_sim, E_ratio_975) %>%
  rename(E_ratio = E_ratio_975) %>%
  mutate(group = "Upper bound")

E_both <- E_both_mean %>%
  rbind(E_both_025) %>%
  rbind(E_both_975)

ggplot(data=E_both, aes(x=hazard_sim, y=E_ratio, group = group, color = group)) +
  geom_line()+
  geom_point() +
  ylim(c(-1, 25)) +
  xlab("Policy shock k") +
  ylab("Change in expenditure reduction (%)")
