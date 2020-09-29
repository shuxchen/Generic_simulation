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
  distinct()

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

# generate one set of simulated time, for second entrants with PC
genericPIV_MEPS_2 <- genericPIV_MEPS %>%
  filter(order == 2)

h_PIV_2 <- h_PIV %>%
  filter(strata == "ncompetitor=1")

genericPIV_MEPS_2_predicted <- predict(model_PIV_PWPGT_post, genericPIV_MEPS_2, type="risk",se.fit=TRUE, reference = "strata")
genericPIV_MEPS_2$predicted <- genericPIV_MEPS_2_predicted$fit
genericPIV_MEPS_2$runif <- runif(nrow(genericPIV_MEPS_2), min=0, max=1)
genericPIV_MEPS_2$predicted_h <- - (genericPIV_MEPS_2$predicted * log(genericPIV_MEPS_2$runif))

genericPIV_MEPS_2$predicted_t <- sapply(genericPIV_MEPS_2$predicted_h, get_time, data = h_PIV_1)
