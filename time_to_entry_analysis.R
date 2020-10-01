#load lib and datasets
setwd("~/Dropbox/Advanced Method Project/Data")
source("Aim2/load_packages.R")
source("Aim3/Generic_simulation/load_functions.R")
load("genericPIV_postGDUFA.RData")
load("genericnoPIV_postGDUFA.RData")
load("MEPS_summary_weighted.Rdata")

# With PC
model_PIV_PWPGT_post <- coxph(Surv(genericPIV_postGDUFA$gaptime_start, genericPIV_postGDUFA$gaptime, entry2) ~ strata(ncompetitor) + route + AG + ETASU + guidance_before + indexyear + cluster(index), method = "breslow", data = genericPIV_postGDUFA)
summary(model_PIV_PWPGT_post)

h_PIV <- basehaz(model_PIV_PWPGT_post, centered = T)
h_PIV <- h_PIV %>% 
  mutate(h2 = hazard*2)

# With PC, for second entrants:
h_PIV_1 <- h_PIV %>% 
  filter(strata == "ncompetitor=1")

h_PIV_1_long <- melt(h_PIV_1, id = "time", measure = c("hazard", "h2"))

h_PIV_1_long$tend <- c(h_PIV_1_long$time[2:nrow(h_PIV_1)], NA)

p1 <- (ggplot(h_PIV_1_long, aes(x=time, y=value, xend=tend, yend=value, color = variable)) +
         geom_vline(aes(xintercept=time), linetype=2, color="grey") +
         geom_point() +  # Solid points to left
         geom_point(aes(x=tend, y=value), shape=1) +  # Open points to right
         geom_segment())  # Horizontal line segments
p1


# MEPS data (add PIV info; include post-GDUFA only)
# first fill in missing prices 
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

MEPS_summary_weighted[is.nan(MEPS_summary_weighted )] <- NA

MEPS_summary_weighted$N_b[MEPS_summary_weighted$N_b == 0] <- NA
MEPS_summary_weighted$N_g[MEPS_summary_weighted$N_g == 0] <- NA

MEPS_summary_weighted <- MEPS_summary_weighted %>%
  group_by(index) %>%
  fill(P_b, P_g, N_b, N_g)

genericPIV_postGDUFA_id <- genericPIV_postGDUFA %>%
  distinct(index)

MEPS_PIV <- MEPS_summary_weighted %>%
  inner_join(genericPIV_postGDUFA_id, by = "index") %>%
  dplyr::select(-(15:28)) %>%
  filter(year > 2012) %>%
  distinct() %>%
  ungroup()

genericnoPIV_postGDUFA_id <- genericnoPIV_postGDUFA %>%
  distinct(index)
MEPS_noPIV <- MEPS_summary_weighted %>%
  inner_join(genericnoPIV_postGDUFA_id, by = "index") %>%
  dplyr::select(-(15:28)) %>%
  filter(year > 2012) %>%
  distinct()

# Missingness too high

## Keep drugs that can be affected:
# With PC: keep N >= 2
MEPS_PIV_2 <- MEPS_PIV %>%
  filter(competitor >= 2,
         !is.na(P_g),
         !is.na(P_b))

# Without PC: keep N >= 2
MEPS_noPIV_1 <- MEPS_noPIV %>%
  filter(competitor >= 1,
         !is.na(P_g),
         !is.na(P_b))
  
# get index for corresponding drugs in MEPS (with and without PC)
MEPS_PIV_id <- MEPS_PIV_2 %>%
  distinct(index)

MEPS_noPIV_id <- MEPS_noPIV_1 %>%
  distinct(index)

# get OB subsets
genericPIV_MEPS <- genericPIV_postGDUFA %>%
  inner_join(MEPS_PIV_id)

genericnoPIV_MEPS <- genericnoPIV_postGDUFA %>%
  inner_join(MEPS_noPIV_id)

# generate one set of simulated time, for second and third entrants with PC
genericPIV_MEPS_2 <- genericPIV_MEPS %>%
  filter(order == 2)

genericPIV_MEPS_3 <- genericPIV_MEPS %>%
  filter(order == 3)

h_PIV_2 <- h_PIV %>%
  filter(strata == "ncompetitor=1")

h_PIV_3 <- h_PIV %>%
  filter(strata == "ncompetitor=2")

genericPIV_MEPS_2_predicted <- predict(model_PIV_PWPGT_post, genericPIV_MEPS_2, type="risk",se.fit=TRUE, reference = "strata")
genericPIV_MEPS_2$predicted <- genericPIV_MEPS_2_predicted$fit
genericPIV_MEPS_2$runif <- runif(nrow(genericPIV_MEPS_2), min=0, max=1)
genericPIV_MEPS_2$predicted_h <- - (genericPIV_MEPS_2$predicted * log(genericPIV_MEPS_2$runif))

genericPIV_MEPS_2$predicted_t <- sapply(genericPIV_MEPS_2$predicted_h, get_time, data = h_PIV_1)

genericPIV_MEPS_3_predicted <- predict(model_PIV_PWPGT_post, genericPIV_MEPS_3, type="risk",se.fit=TRUE, reference = "strata")
genericPIV_MEPS_3$predicted <- genericPIV_MEPS_3_predicted$fit
genericPIV_MEPS_3$runif <- runif(nrow(genericPIV_MEPS_3), min=0, max=1)
genericPIV_MEPS_3$predicted_h <- - (genericPIV_MEPS_3$predicted * log(genericPIV_MEPS_3$runif))

genericPIV_MEPS_3$predicted_t <- sapply(genericPIV_MEPS_3$predicted_h, get_time, data = h_PIV_1)

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
  mutate(P_b_simulated = P_b * (1 + 0.01 * competitor_diff),
         P_g_simulated = P_g * (1 - 0.08 * competitor_diff),
         E = P_b * N_b + P_g * N_g,
         E_simulated = P_b_simulated * N_b + P_g_simulated * N_g)

E_diff <- sum(MEPS_PIV_simulated_v1$E) - sum(MEPS_PIV_simulated_v1$E_simulated)

E_ratio <- E_diff/sum(MEPS_PIV_simulated_v1$E)
