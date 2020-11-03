#load lib and datasets
setwd("~/Dropbox/Advanced Method Project/Data")
source("Aim2/load_packages.R")
source("Aim3/Generic_simulation/load_functions.R")
load("genericPIV_postGDUFA.RData")
load("genericnoPIV_postGDUFA.RData")
load("MEPS_summary_weighted.Rdata")

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


genericPIV_postGDUFA <- genericPIV_postGDUFA %>%
  mutate(same_day = ifelse(gaptime <= 0.001, 1, 0))


genericPIV_postGDUFA_second <- genericPIV_postGDUFA %>%
  filter(ncompetitor == 1)

genericPIV_postGDUFA_third <- genericPIV_postGDUFA %>%
  filter(ncompetitor == 2)

#first-part: logistic
logit_PIV_second <- glm(same_day ~ route + AG + ETASU + guidance_before + indexyear, data = genericPIV_postGDUFA_second, family = "binomial")
summary(logit_PIV_second)

logit_PIV_third <- glm(same_day ~ route + AG + guidance_before + indexyear, data = genericPIV_postGDUFA_third, family = "binomial")
summary(logit_PIV_third)

genericPIV_postGDUFA_second <- genericPIV_postGDUFA %>%
  filter(ncompetitor == 1,
         gaptime > 0.001)

genericPIV_postGDUFA_third <- genericPIV_postGDUFA %>%
  filter(ncompetitor == 2,
         gaptime > 0.001)

#model_PIV_second <- flexsurvreg(Surv(gaptime, entry2) ~ AG + guidance_before + indexyear, data = genericPIV_postGDUFA_second, dist = "gengamma")
#summary(model_PIV_second, type = "link")

#flexsurvreg(Surv(gaptime, entry2) ~ AG + guidance_before + indexyear, data = genericPIV_postGDUFA_second, dist = "exp")
model_PIV_second <- flexsurvreg(Surv(gaptime, entry2) ~ AG + guidance_before + indexyear, data = genericPIV_postGDUFA_second, dist = "weibull")
#flexsurvreg(Surv(gaptime, entry2) ~ AG + guidance_before + indexyear, data = genericPIV_postGDUFA_second, dist = "lognormal")

#model_PIV_second <- survreg(Surv(gaptime, entry2) ~ AG + guidance_before + indexyear, data = genericPIV_postGDUFA_second, dist = "lognormal")

model_PIV_third <- flexsurvreg(Surv(gaptime, entry2) ~ route + AG + guidance_before + indexyear, data = genericPIV_postGDUFA_third, dist = "weibull")
#exclude ETASU as no ETASU = 1


dists <- c("exp", "weibull", "gamma", 
           "lognormal", "gengamma")
dists_long <- c("Exponential", "Weibull",
                "Gamma", "Lognormal", 
                "Generalized gamma")
parametric_haz <- vector(mode = "list", length = length(dists))
for (i in 1:length(dists)){
  fit <- flexsurvreg(Surv(gaptime, entry2) ~ AG + guidance_before + indexyear, data = genericPIV_postGDUFA_second, dist = dists[i]) 
  parametric_haz[[i]] <- summary(fit, type = "hazard", 
                                 ci = FALSE, tidy = TRUE)
  parametric_haz[[i]]$method <- dists_long[i]
}
parametric_haz <- rbindlist(parametric_haz)
haz <- parametric_haz
n_dists <- length(dists) 
ggplot(haz, aes(x = time, y = est, col = method, linetype = method)) +
  geom_line() +
  xlab("Days") + ylab("Hazard") + 
  scale_colour_manual(name = "", 
                      values = c("black", rainbow(n_dists))) +
  scale_linetype_manual(name = "",
                        values = c(1, rep_len(2:6, n_dists)))


model.matrix(model_PIV_second)
model_PIV_second$res[-(1:2),"est"]

exp(model.matrix(model_PIV_second) %*% model_PIV_second$res[-(1:2),"est"])



m1 <- summary(model_PIV_second)[[1]]
plot(est~time, data=m1)



