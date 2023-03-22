# Interactive effects of long-term wave exposure and other stressors on fish 
#   welfare on Atlantic salmon (Salmo salar) farms
# Szewczyk et al. 2023
#
# Mortality model
#
# This script fits the model for weekly mortality rate


# switches ----------------------------------------------------------------

nChains <- 3 # 6
iter <- 100 # 4000
warmup <- 50 # 3000



# set up ------------------------------------------------------------------

library(tidyverse); library(glue); library(brms)
wk.df <- read_csv("data/pub_dataset.csv") 



# formula -----------------------------------------------------------------

pred.mort <- list(
  waveHt="waveHeight", 
  density="biomass_density_z",
  time="days_since_start_z",
  env=paste0(c("temperature", "waterSpeed", "windSpeed"), "_z"),
  env_sq=paste0("I(", c("temperature"), "_z^2)"),
  trt=paste0("wksSince", c("AntiAGD", "Mech", "Medi"), "_z"),
  trt_sq=paste0("I(wksSince", c("AntiAGD", "Mech", "Medi"), "_z^2)"),
  disease=paste0(c("lice_total", "agd_score"), "_z")
)

mort.terms <- c(paste(unlist(pred.mort)), 
                paste(pred.mort$waveHt, 
                      c(pred.mort$density, 
                        pred.mort$time, 
                        pred.mort$env, pred.mort$env_sq, 
                        pred.mort$trt, pred.mort$trt_sq, 
                        pred.mort$disease), 
                      sep=":"),
                paste("temperature_z", 
                      c(pred.mort$density,
                        pred.mort$trt, 
                        pred.mort$trt_sq, 
                        pred.mort$disease), 
                      sep=":"),
                paste("lice_total_z",
                      c(pred.mort$density, 
                        pred.mort$trt, 
                        pred.mort$trt_sq, 
                        "waterSpeed_z"),
                      sep=":"),
                paste("agd_score_z",
                      c(pred.mort$density, 
                        pred.mort$trt[1:2], 
                        pred.mort$trt_sq[1:2], 
                        "waterSpeed_z"),
                      sep=":")
)

mort.rand.terms <- grep("waveHeight", mort.terms, invert=T, value=T)

mort.df <- wk.df %>%
  mutate(lice_total=log1p(lice_total),
         agd_score=log1p(agd_score)) %>% 
  mutate(across(starts_with("wksSince"), ~log1p(.x))) 

mort.df <- mort.df %>%
  select(farm, cage, day, 
         lmortRate, biomass_density,
         days_since_start,
         waveHeight, 
         temperature, waterSpeed, windSpeed,
         lice_total, agd_score,
         wksSinceAntiAGD, wksSinceMech, wksSinceMedi) %>%
  filter(complete.cases(.)) %>%
  mutate(across(where(is.numeric), ~c(scale(.x)), .names="{.col}_z"))
write_csv(mort.df, "out/mort_data.csv")



# run model ---------------------------------------------------------------

mort.fit <- brm(bf(glue("lmortRate_z ~ ", 
                        "{paste(mort.terms, collapse='+')}",
                        "+ (1+{paste(mort.rand.terms, collapse='+')}|farm/cage)")), 
                family=gaussian(), 
                data=mort.df, cores=nChains, chains=nChains, init=0, 
                iter=iter, warmup=warmup, refresh=5,
                prior=c(prior(horseshoe(3, par_ratio=0.3), "b"),
                        prior(normal(0, 1), "Intercept"),
                        prior(normal(0, 1), "sigma"),
                        prior(normal(0, .1), "sd")), 
                control=list(adapt_delta=0.99, max_treedepth=20))
saveRDS(mort.fit, "out/mort_fit.rds")

