# Interactive effects of long-term wave exposure and other stressors on fish 
#   welfare on Atlantic salmon (Salmo salar) farms
# Szewczyk et al. 2023
#
# Helper functions
#
# This script contains functions used to process the fitted models




#' Convert log1p-transformed weeks  to natural-scale days
#' 
#' Used internally in `conditional_effects2()`
#'
#' @param x Numeric vector
#' 
expm1_wkToDay <- function(x) {
  expm1(x) * 7
}






#' Convert logit-scale value to percentage
#' 
#' Used internally in `conditional_effects2()`
#'
#' @param x Numeric vector
#' 
logit_prop_to_pct <- function(x) {
  boot::inv.logit(x)*100
}






#' Summarise posterior effects
#'
#' @param out.ls List of brms output
#' @param farm.df Dataframe of farm information
#' @param nGrpThresh Threshold for including (number of farms with nonzero 95% CIs)
#' @param probs Vector of CI probabilities to calculate
#' 
summarise_effects <- function(out.ls, farm.df, nGrpThresh=1,
                              probs=c(0.5, 0.8, 0.9, 0.95)) {
  library(tidybayes)
  # population-level effects
  eff.pop <- imap_dfr(out.ls, 
                      ~fixef(.x, summary=F) %>% 
                        as_tibble() %>% 
                        gather_variables() %>%
                        mean_qi(.width=probs) %>%
                        mutate(nonZero=sign(.lower)==sign(.upper),
                               model=.y)) %>%
    mutate(term=.variable %>%
             str_remove_all(., "I|E2|hu_"))
  
  # group-level effects
  eff.grp.df <- imap_dfr(out.ls, 
                         ~as_tibble(coef(.x, summary=F)$farm) %>%
                           gather_variables() %>%
                           mean_qi(.width=probs) %>%
                           mutate(farm=str_split_fixed(.variable, "\\.", 2)[,1],
                                  .variable=str_split_fixed(.variable, "\\.", 2)[,2],
                                  nonZero=sign(.lower)==sign(.upper),
                                  model=.y)) %>%
    left_join(., farm.df) %>%
    mutate(term=.variable %>%
             str_remove_all(., "I|E2|hu_") %>%
             gsub("(exposure+)(:)(.+)", "\\3\\2\\1", .))
  
  # non-zero
  sig.eff <- bind_rows(
    eff.grp.df %>% 
      group_by(model, .variable, .width) %>%
      summarise(nSig=sum(nonZero)) %>%
      ungroup %>%
      mutate(level="farm") %>%
      filter(nSig >= nGrpThresh),
    eff.pop %>%
      filter(nonZero) %>%
      select(model, .variable, .width) %>%
      mutate(nSig=8, level="population")
  ) %>%
    mutate(term=.variable %>% str_remove_all(., "I|E2|hu_"))
  
  eff.names <- map(names(out.ls), 
                   ~grep("Intercept", 
                         bind_rows(eff.grp.df %>% 
                                     select(.variable, model) %>%
                                     filter(model==.x),
                                   eff.pop %>% 
                                     select(.variable, model) %>%
                                     filter(model==.x))$.variable,
                         value=T, invert=T) %>%
                     str_remove_all(., "I|E2|hu_") %>%
                     unique()) %>% 
    setNames(names(out.ls))
  eff.single <- map(eff.names, ~grep(":", .x, invert=T, value=T))
  
  sig.eff.names <- sig.eff %>% filter(!grepl("Intercept", .variable)) %>%
    group_by(model, .width) %>% distinct(term) %>% 
    group_by(model) %>% group_split(.keep=F) %>%
    setNames(names(out.ls))
  sig.eff.single <- map(sig.eff.names, ~filter(.x, !grepl(":", term)))
  
  eff.clean <- tribble(~og, ~clean, ~unit,
                       "agd_score_z", "AGD score", "mean",
                       "lice_total_z", "Lice/fish", "mean",
                       "biomass_density_z", "Mass", "kg/m3",
                       "days_since_start_z", "Time since stocking", "d",
                       "temperature_z", "Weekly temperature", "\u00B0C",
                       "wksSinceAntiAGD_z", "Time since anti-AGD trt.", "d",
                       "wksSinceMedi_z", "Time since anti-lice medi. trt.", "d",
                       "wksSinceMech_z", "Time since anti-lice mech. trt.", "d",
                       "waterSpeed_z", "Weekly water speed", "m/s",
                       "windSpeed_z", "Weekly wind speed", "m/s",
                       "waveHeight", "50-yr sig. wave ht.", "m",
                       "lmortRatez", "Mortality rate", "%")
  
  return(list(eff.pop=eff.pop, 
              eff.grp=eff.grp.df, 
              eff.names=eff.names,
              eff.single=eff.single,
              sig.eff=sig.eff, 
              sig.eff.names=sig.eff.names, 
              sig.eff.single=sig.eff.single,
              eff.clean=eff.clean))
}







#' Customized version of brms::conditional_effects()
#' 
#' Returns a list of dataframes for use with ggplot2
#'
#' @param x Output from brms
#' @param effects Character vector of effects to calculate
#' @param conditions # Character vector of conditions to calculate across
#' @param re Logical; Include random effects? F = population-level
#' @param probs Vector of quantiles to calculate
#' @param npoints Points in the x and y direction for rasterised output
#' @param interaction Representation of interaction by "raster" or "line" 
#' @param lb Lower bound in predictions; currently ignored
#' @param ub Upper bound in predictions
#' @param modType "gaussian" or "hurdle" 
#' @param z_scale Logical; plot on z-transformed scale? F = natural scale
#' @param dat.df Dataframe with data
#' @param names Effect names
#'
conditional_effects2 <- function(x, effects, conditions=NULL, re=F, 
                                 probs=c(0.025, 0.975), 
                                 npoints=c(10, 10),
                                 interaction="raster",
                                 lb=-Inf, ub=Inf, modType="hurdle",
                                 z_scale=T, dat.df=NULL, names=NULL) {
  # set up list of dataframes, using *actual* farm:cage
  base.df <- x$data %>%
    summarise(across(c(-starts_with("I("), -farm, -cage, -`farm:cage`), mean)) 
  if(re) {
    base.df <- base.df %>%
      bind_cols(x$data %>%
                  select(farm, cage, `farm:cage`) %>%
                  group_by(farm, cage, `farm:cage`) %>%
                  slice_head(n=1)) %>%
      select(-waveHeight) %>%
      left_join(., conditions) %>% 
      group_by(farm, cage, `farm:cage`)
  }
  
  eff.mx <- str_split_fixed(effects, ":", 2)
  df.ls <- gg.ls <- vector("list", nrow(eff.mx)) %>% setNames(effects)
  for(i in seq_along(df.ls)) {
    nEffects <- sum(eff.mx[i,] != "")
    effect1 <- seq(min(x$data[[eff.mx[i,1]]]),
                   max(x$data[[eff.mx[i,1]]]),
                   length.out=npoints[1])
    if(nEffects==1) {
      df.ls[[i]] <- base.df %>%
        rowwise() %>%
        mutate(effect1__=list(effect1)) %>%
        unnest(effect1__)
      df.ls[[i]][[eff.mx[i,1]]] <- df.ls[[i]]$effect1__
    } else {
      effect2 <- seq(min(x$data[[eff.mx[i,2]]]),
                     max(x$data[[eff.mx[i,2]]]),
                     length.out=npoints[2])
      eff.grd <- expand_grid(effect1=effect1,
                             effect2=effect2)
      df.ls[[i]] <- base.df %>%
        rowwise() %>%
        mutate(effect1__=list(eff.grd$effect1),
               effect2__=list(eff.grd$effect2)) %>%
        unnest(c(effect1__, effect2__))
      df.ls[[i]][[eff.mx[i,1]]] <- df.ls[[i]]$effect1__
      df.ls[[i]][[eff.mx[i,2]]] <- df.ls[[i]]$effect2__
    }
    
    # generate predictions for each dataframe and each dpar, using *actual* farm:cage
    if(modType=="hurdle") {
      if(re) {
        pred.mu <- posterior_epred(x, newdata=df.ls[[i]], dpar="mu", re_formula=NULL)
        pred.hu <- posterior_epred(x, newdata=df.ls[[i]], dpar="hu", re_formula=NULL)
      } else {
        pred.mu <- posterior_epred(x, newdata=df.ls[[i]], dpar="mu", re_formula=NA)
        pred.hu <- posterior_epred(x, newdata=df.ls[[i]], dpar="hu", re_formula=NA)
      }
      pred <- pmin((1-pred.hu) * exp(pred.mu), ub) 
    } else if(modType=="gaussian") {
      if(re) {
        pred <- posterior_epred(x, newdata=df.ls[[i]], re_formula=NULL) 
      } else {
        pred <- posterior_epred(x, newdata=df.ls[[i]], re_formula=NA) 
      }
    }
    
    pred.df <- as_tibble(t(pred)) %>% 
      mutate(obs.id=row_number()) %>%
      pivot_longer(starts_with("V"), names_to="iter", values_to="value")
    
    # rescale from z if needed
    if(!z_scale) {
      trans <- list(lmortRatez=logit_prop_to_pct,
                    agd_score_z=expm1,
                    lice_total_z=expm1,
                    licetotal=expm1,
                    wksSinceAntiAGD_z=expm1_wkToDay,
                    wksSinceMedi_z=expm1_wkToDay,
                    wksSinceMech_z=expm1_wkToDay)
      resp <- ifelse(x$formula$resp=="lmortRatez", "lmortRate", "lepFe")
      if(resp=="lmortRate") {
        pred_unZ <- pred.df$value * sd(dat.df[[resp]]) + mean(dat.df[[resp]])
      } else {
        pred_unZ <- pred.df$value
      }
      if(x$formula$resp %in% names(trans)) {
        pred_unZ <- trans[[x$formula$resp]](pred_unZ)
      }
      pred.df$value <- pred_unZ
      
      eff.mx_unZ <- gsub("_z", "", eff.mx)
      if(eff.mx[i,1] != "waveHeight") {
        eff1_unZ <- df.ls[[i]]$effect1__ * sd(dat.df[[eff.mx_unZ[i,1]]]) + 
          mean(dat.df[[eff.mx_unZ[i,1]]])
        if(eff.mx[i,1] %in% names(trans)) {
          eff1_unZ <- trans[[eff.mx[i,1]]](eff1_unZ)
        }
        df.ls[[i]]$effect1__ <- eff1_unZ 
      }
      
      if(nEffects>1) {
        if(eff.mx[i,2] != "waveHeight") {
          eff2_unZ <- df.ls[[i]]$effect2__ * sd(dat.df[[eff.mx_unZ[i,2]]]) + 
            mean(dat.df[[eff.mx_unZ[i,2]]])
          if(eff.mx[i,2] %in% names(trans)) {
            eff2_unZ <- trans[[eff.mx[i,2]]](eff2_unZ)
          } 
          df.ls[[i]]$effect2__ <- eff2_unZ 
        }
      }
    }
    
    # if re==TRUE, summarise by farm, using *actual* farm:cage
    if(re) {
      if(nEffects==1) {
        df.ls[[i]] <- df.ls[[i]] %>% ungroup %>% 
          mutate(obs.id=row_number()) %>%
          full_join(pred.df) %>%
          group_by(farm, waveHeight, effect1__) %>%
          summarise(estimate__=median(value),
                    lower__=quantile(value, probs[1]),
                    upper__=quantile(value, probs[2]))  
      } else {
        df.ls[[i]] <- df.ls[[i]] %>% ungroup %>% 
          mutate(obs.id=row_number()) %>%
          full_join(pred.df) %>%
          group_by(farm, waveHeight, effect1__, effect2__) %>%
          summarise(estimate__=median(value),
                    lower__=quantile(value, probs[1]),
                    upper__=quantile(value, probs[2])) 
      }
    } else {
      if(nEffects==1) {
        df.ls[[i]] <- df.ls[[i]] %>% ungroup %>% 
          mutate(obs.id=row_number()) %>%
          full_join(pred.df) %>%
          group_by(effect1__) %>%
          summarise(estimate__=median(value),
                    lower__=quantile(value, probs[1]),
                    upper__=quantile(value, probs[2]))
      } else {
        df.ls[[i]] <- df.ls[[i]] %>% ungroup %>% 
          mutate(obs.id=row_number()) %>%
          full_join(pred.df) %>%
          group_by(effect1__, effect2__) %>%
          summarise(estimate__=median(value),
                    lower__=quantile(value, probs[1]),
                    upper__=quantile(value, probs[2]))
      }
    }
    
    # make ggplots
    if(is.null(names)) {
      eff.labs <- list(resp=x$formula$resp, 
                       eff1=eff.mx[i,1],
                       eff2=eff.mx[i,2])
    } else {
      eff.labs <- list(resp=with(filter(names, og==x$formula$resp), 
                                 glue("{clean} ({unit})")),
                       eff1=with(filter(names, og==eff.mx[i,1]), 
                                 glue("{clean} ({unit})")),
                       eff2=with(filter(names, og==eff.mx[i,2]), 
                                 glue("{clean} ({unit})")))
    }
    if(nEffects==1) {
      gg.ls[[i]] <- ggplot(df.ls[[i]], aes(effect1__, estimate__)) + 
        {if(re) geom_line(aes(group=farm, colour=waveHeight))} +
        {if(!re) geom_line(colour="blue")} +
        labs(x=eff.labs$eff1, y=eff.labs$resp) 
    } else {
      if(interaction=="raster") {
        gg.ls[[i]] <- ggplot(df.ls[[i]], aes(effect1__, effect2__, fill=estimate__)) +
          geom_raster() +
          {if(re) facet_wrap(~farm)} +
          labs(x=eff.labs$eff1, y=eff.labs$eff2) +
          scale_fill_viridis_c(eff.labs$resp)
      } else {
        gg.ls[[i]] <- ggplot(df.ls[[i]], aes(effect1__, estimate__, 
                                             colour=effect2__, group=effect2__)) +
          geom_line() +
          {if(re) facet_wrap(~farm)} +
          {if(grepl("trt", eff.labs$eff1)) scale_x_continuous(limits=c(0,100), breaks=seq(0,100,by=25))} +
          labs(x=eff.labs$eff1, y=eff.labs$resp) +
          scale_colour_viridis_c(eff.labs$eff2)
      }
    }
  }
  
  # return list of dataframes
  return(gg.ls)
}
