# Interactive effects of long-term wave exposure and other stressors on fish 
#   welfare on Atlantic salmon (Salmo salar) farms
# Szewczyk et al. 2023
#
# Sea lice model
#
# This script fits the model for the weekly sea lice infection rate


# switches ----------------------------------------------------------------

nChains <- 3 # 6
iter <- 100 # 4000
warmup <- 50 # 3000



# set up ------------------------------------------------------------------

library(tidyverse); library(glue); library(brms)
wk.df <- read_csv("data/pub_dataset.csv") 



# formula -----------------------------------------------------------------

pred.lice <- list(
  waveHt="waveHeight", 
  density="biomass_density_z",
  time="days_since_start_z",
  env=paste0(c("temperature", "waterSpeed", "windSpeed"), "_z"),
  env_sq=paste0("I(", c("temperature"), "_z^2)"),
  trt=paste0("wksSince", c("AntiAGD", "Mech", "Medi"), "_z"),
  trt_sq=paste0("I(wksSince", c("AntiAGD", "Mech", "Medi"), "_z^2)")
)

lice.terms <- c(paste(unlist(pred.lice)), 
                paste(pred.lice$waveHt, 
                      c(pred.lice$density, 
                        pred.lice$time, 
                        pred.lice$env, pred.lice$env_sq, 
                        pred.lice$trt, pred.lice$trt_sq,
                        "agd_score_z"), 
                      sep=":"),
                paste("temperature_z", 
                      c(pred.lice$density,
                        pred.lice$trt,
                        pred.lice$trt_sq), 
                      sep=":"),
                paste("agd_score_z", 
                      c(pred.lice$density, 
                        "waterSpeed_z"),
                      sep=":")
)

lice.rand.terms <- grep("waveHeight", lice.terms, invert=T, value=T)

lice.df <- wk.df %>%
  mutate(lep_fe=log(lep_fe+1),
         agd_score=log(agd_score+1)) %>% 
  mutate(across(starts_with("wksSince"), ~log(.x+1))) 

lice.df <- lice.df %>%
  select(farm, cage, day,  
         biomass_density,
         days_since_start,
         waveHeight, 
         temperature, waterSpeed, windSpeed,
         lep_fe, agd_score,
         wksSinceAntiAGD, wksSinceMech, wksSinceMedi) %>%
  filter(complete.cases(.)) %>%
  mutate(across(where(is.numeric), ~c(scale(.x)), .names="{.col}_z"))
write_csv(lice.df, "out/lepFe_data.csv")



# run model ---------------------------------------------------------------

lice.fit <- brm(bf(glue("lep_fe ~",
                        "{paste(lice.terms, collapse='+')}",
                        "+ (1+{paste(lice.rand.terms, collapse='+')}|farm/cage)"), 
                   glue("hu ~",
                        "{paste(lice.terms, collapse='+')}",
                        "+ (1+{paste(lice.rand.terms, collapse='+')}|farm/cage)")),
                family=hurdle_lognormal(), 
                data=lice.df, cores=nChains, chains=nChains, init=0, 
                iter=iter, warmup=warmup, refresh=5,
                prior=c(prior(horseshoe(3, par_ratio=0.3), "b"),
                        prior(horseshoe(3, par_ratio=0.3), "b", dpar="hu"),
                        prior(normal(0, 1), "Intercept"),
                        prior(normal(0, 0.5), "Intercept", dpar="hu"),
                        prior(normal(0, 0.1), "sd"),
                        prior(normal(0, 0.1), "sd", dpar="hu"),
                        prior(normal(0, 1), "sigma")), 
                control=list(adapt_delta=0.99, max_treedepth=20))
saveRDS(lice.fit, "out/lepFe_fit.rds")
