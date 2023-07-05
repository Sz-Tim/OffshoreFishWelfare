# Interactive effects of long-term wave exposure and other stressors on fish 
#   welfare on Atlantic salmon (Salmo salar) farms
# Szewczyk et al. 2023
#
# Summaries and figures
#
# This script summarizes model output, produces figures, and extracts summaries
# referenced in the manuscript


# setup -------------------------------------------------------------------

library(tidyverse)
library(brms)
library(glue)
library(bayesplot)
library(ggpubr)
library(lubridate)
library(HDInterval)

source("code/00_fn.R")
cm <- readRDS("data/cmr_cmaps.RDS")  # CMasher color schemes: https://cmasher.readthedocs.io/
theme_set(theme_classic())
farm_i <- read_csv("data/farm_i.csv")

# colour schemes and axes
wave_col <- scale_colour_gradientn("50-yr sig.\nwave ht. (m)", colours=cm$amber[1:(0.9*length(cm$amber))])
temp_col <- scale_colour_viridis_c("Temp. (\u00B0C)", option="plasma", end=0.9)
mass_col <- scale_colour_viridis_c("Mass (kg/m3)", option="cividis", end=0.9)
agd_col <- scale_colour_distiller("AGD score", type="seq", palette="PuRd", direction=1, limits=c(0,3))
lice_col <- scale_colour_distiller("Lice/fish\n(all spp.)", type="seq", palette="PuBuGn", direction=-1)
mort_y <- scale_y_continuous("Weekly mortality rate (%)", limits=c(0, 0.3), breaks=c(0, 0.1, 0.2))
lice_y <- scale_y_continuous("Mean lice per fish", limits=c(0, 1), breaks=seq(0, 1, by=0.5))
conf_cols <- c("80-90%"="#67a9cf", "90-95%"="#1c9099", "\u2265 95%"="#016c59")
pub_theme <- theme(plot.title=element_text(size=9), 
                   axis.title=element_text(size=9),
                   axis.text=element_text(size=8),
                   legend.text=element_text(size=7),
                   legend.title=element_text(size=8),
                   strip.text=element_text(size=9))




# mortality effects -------------------------------------------------------

# load data and output
mort.dat <- read_csv("out/mort_data.csv")
mort.fit <- list(mort=readRDS("out/mort_fit.rds"))


# summarise effects and calculate conditional effects
mort.eff <- summarise_effects(mort.fit, farm_i, nGrpThresh=8)
ce.1D <- conditional_effects2(mort.fit[[1]], 
                              effects=grep("waveHeight",
                                           mort.eff$eff.single[[1]], 
                                           invert=T, value=T),
                              conditions=farm_i %>% select(farm, waveHeight),
                              re=T, npoints=c(100, 100), 
                              modType="gaussian", 
                              z_scale=F, dat.df=mort.dat, names=mort.eff$eff.clean)
ce <- conditional_effects2(mort.fit[[1]],
                           effects=mort.eff$eff.names[[1]] %>% 
                             gsub("(biomass_density_z+)(:)(.+)", "\\3\\2\\1", .) %>%
                             gsub("(temperature_z+)(:)(.+)", "\\3\\2\\1", .) %>%
                             gsub("(waveHeight+)(:)(.+)", "\\3\\2\\1", .),
                           re=F, npoints=c(100, 10), interaction="lines",
                           modType="gaussian", 
                           z_scale=F, dat.df=mort.dat, names=mort.eff$eff.clean)
saveRDS(ce.1D, "out/mort_ce1D.rds")
saveRDS(ce, "out/mort_ceAll.rds")



# * Fig 4 -----------------------------------------------------------------
# create plot of effects with >80% confidence
ci_level <- 0.8
ce.1D <- readRDS("out/mort_ce1D.rds")
ce <- readRDS("out/mort_ceAll.rds")
eff.sig <- filter(mort.eff$sig.eff.names[[1]], .width==ci_level)$term %>% 
  gsub("(temperature_z+)(:)(.+)", "\\3\\2\\1", .) %>%
  gsub("(waveHeight+)(:)(.+)", "\\3\\2\\1", .)
eff.sigLevel <- mort.eff$sig.eff.names[[1]] %>%
  arrange(desc(.width)) %>%
  group_by(term) %>%
  slice_head(n=1) %>%
  filter(.width > ci_level) %>%
  mutate(term=term %>% 
           gsub("(temperature_z+)(:)(.+)", "\\3\\2\\1", .) %>%
           gsub("(waveHeight+)(:)(.+)", "\\3\\2\\1", .),
         conf=case_when(.width==0.8 ~ "80-90%",
                        .width==0.9 ~ "90-95%",
                        .width==0.95 ~ "\u2265 95%"),
         col=conf_cols[conf]) 
ce.sig <- ce[eff.sigLevel$term]
plot.ls <- c(map(seq_along(ce.sig), ~ce.sig[[.x]] + 
                   annotate("text", hjust=0, vjust=0.5, size=3,
                            x=ce.sig[[.x]]$data$effect1__[1], y=mort_y$limits[2], 
                            label=paste(eff.sigLevel$conf[.x], "conf."), 
                            colour=eff.sigLevel$col[.x]) +
                   mort_y + wave_col + pub_theme + 
                   theme(legend.key.height=unit(0.1, "cm"))))
mort.80.p <- ggarrange(plotlist=plot.ls, ncol=2, nrow=1, 
                       common.legend=T, legend="bottom", labels="auto")
ggsave("figs/fig_4.jpg", mort.80.p, width=140, height=70, units="mm", dpi=300)



# * Summary values --------------------------------------------------------
# Conditional and marginal R2
bayes_R2(mort.fit[[1]])
bayes_R2(mort.fit[[1]], re_formula=NA)

# Time: Mortality at 1 year
ce$`days_since_start_z:waveHeight`$data %>% 
  ungroup %>% 
  filter(abs(effect1__ - 365) < 2) %>%
  mutate(across(one_of("estimate__", "lower__", "upper__"), 
                ~.x/first(.x), 
                .names="propDiff_{.col}"))

# Temp: Mortality at mean August temp
mnAugTemp <- mort.dat %>% 
  mutate(month=lubridate::month(day)) %>% 
  group_by(month) %>% 
  summarise(temperature=mean(temperature)) %>%
  filter(month==8) %>% 
  select(temperature) %>%
  unlist
ce$`temperature_z:waveHeight`$data %>% 
  ungroup %>% 
  filter(abs(effect1__ - mnAugTemp) < 0.05) %>%
  mutate(across(one_of("estimate__", "lower__", "upper__"), 
                ~.x/first(.x), 
                .names="propDiff_{.col}"))





# sea lice effects --------------------------------------------------------

# load data and output
lice.dat <- read_csv("out/lepFe_data.csv")
lice.fit <- list(lepFe=readRDS("out/lepFe_fit.rds"))
lepFe.png <- png::readPNG("data/L_salmonis_gravid_female.png")


# summarise effects and calculate conditional effects
lice.eff <- summarise_effects(lice.fit, farm_i, nGrpThresh=8)
ce.1D <- conditional_effects2(lice.fit[[1]], 
                              effects=grep("waveHeight",
                                           lice.eff$eff.single[[1]], 
                                           invert=T, value=T),
                              conditions=farm_i %>% select(farm, waveHeight),
                              re=T, npoints=c(100, 100), 
                              modType="hurdle", ub=15,
                              z_scale=F, dat.df=lice.dat, names=lice.eff$eff.clean)
ce <- conditional_effects2(lice.fit[[1]],
                           effects=lice.eff$eff.names[[1]] %>% 
                             gsub("(biomass_density_z+)(:)(.+)", "\\3\\2\\1", .) %>%
                             gsub("(temperature_z+)(:)(.+)", "\\3\\2\\1", .) %>%
                             gsub("(waveHeight+)(:)(.+)", "\\3\\2\\1", .),
                           re=F, npoints=c(100, 10), interaction="lines",
                           modType="hurdle", ub=15,
                           z_scale=F, dat.df=lice.dat, names=lice.eff$eff.clean)
saveRDS(ce.1D, "out/lepFe_ce1D.rds")
saveRDS(ce, "out/lepFe_ceAll.rds")



# * Fig 6 -----------------------------------------------------------------
# create plot of effects with >80% confidence
ci_level <- 0.8
ce.1D <- readRDS("out/lepFe_ce1D.rds")
ce <- readRDS("out/lepFe_ceAll.rds")
eff.sigLevel <- lice.eff$sig.eff.names[[1]] %>%
  arrange(desc(.width)) %>%
  group_by(term) %>%
  slice_head(n=1) %>%
  filter(.width >= ci_level) %>%
  mutate(term=term %>% 
           gsub("(biomass_density_z+)(:)(.+)", "\\3\\2\\1", .) %>%
           gsub("(temperature_z+)(:)(.+)", "\\3\\2\\1", .) %>%
           gsub("(waveHeight+)(:)(.+)", "\\3\\2\\1", .),
         conf=case_when(.width==0.8 ~ "80-90%",
                        .width==0.9 ~ "90-95%",
                        .width==0.95 ~ "\u2265 95%"),
         col=conf_cols[conf]) 

ce.1D.sig <- ce.1D[grep(":|^waveHeight", eff.sigLevel$term, invert=T, value=T)]
ce.sig <- ce[eff.sigLevel$term]
plot.ls <- c(map(names(ce.1D.sig),
                 ~ce.sig[[.x]] +
                   geom_point(data=ce.sig[[.x]]$data[1,] %>% mutate(fakeColumn="A"),
                              aes(colour=fakeColumn)) +
                   annotate("text", hjust=0, vjust=0.5, size=3,
                            x=ce.sig[[.x]]$data$effect1__[1], y=lice_y$limits[2], 
                            label=paste(filter(eff.sigLevel, term==.x)$conf, "conf."), 
                            colour=filter(eff.sigLevel, term==.x)$col) +
                   lice_y + pub_theme +
                   scale_colour_manual(values="white") +
                   theme(legend.title=element_text(colour="white"),
                         legend.key=element_blank(),
                         legend.text=element_text(colour="white")) +
                   geom_ribbon(aes(ymin=pmax(lower__, 0), ymax=pmin(upper__, 1)), fill="grey", alpha=0.5) +
                   geom_line(colour="blue")),
             map(ce.sig["agd_score_z:waveHeight"], ~.x + 
                   annotate("text", hjust=0, vjust=0.5, size=3,
                            x=ce.sig[["agd_score_z:waveHeight"]]$data$effect1__[1], y=lice_y$limits[2], 
                            label=paste(filter(eff.sigLevel, term=="agd_score_z:waveHeight")$conf, "conf."), 
                            colour=filter(eff.sigLevel, term=="agd_score_z:waveHeight")$col) +
                   lice_y + wave_col + pub_theme + 
                   theme(legend.key.width=unit(0.1, "cm"))),
             map(ce.sig["waterSpeed_z:agd_score_z"], ~.x + 
                   annotate("text", hjust=0, vjust=0.5, size=3,
                            x=ce.sig[["waterSpeed_z:agd_score_z"]]$data$effect1__[1], y=lice_y$limits[2], 
                            label=paste(filter(eff.sigLevel, term=="waterSpeed_z:agd_score_z")$conf, "conf."), 
                            colour=filter(eff.sigLevel, term=="waterSpeed_z:agd_score_z")$col) +
                   lice_y + agd_col + pub_theme + 
                   theme(legend.key.width=unit(0.1, "cm"))),
             map(ce.sig["temperature_z:waveHeight"], ~.x + 
                   annotate("text", hjust=0, vjust=0.5, size=3,
                            x=ce.sig[["temperature_z:waveHeight"]]$data$effect1__[1], y=lice_y$limits[2], 
                            label=paste(filter(eff.sigLevel, term=="temperature_z:waveHeight")$conf, "conf."), 
                            colour=filter(eff.sigLevel, term=="temperature_z:waveHeight")$col) +
                   lice_y + wave_col + pub_theme + 
                   theme(legend.key.width=unit(0.1, "cm"))),
             map(ce.sig["wksSinceAntiAGD_z:waveHeight"], ~.x + 
                   annotate("text", hjust=0, vjust=0.5, size=3,
                            x=ce.sig[["wksSinceAntiAGD_z:waveHeight"]]$data$effect1__[1], y=lice_y$limits[2], 
                            label=paste(filter(eff.sigLevel, term=="wksSinceAntiAGD_z:waveHeight")$conf, "conf."), 
                            colour=filter(eff.sigLevel, term=="wksSinceAntiAGD_z:waveHeight")$col) +
                   lice_y + wave_col + pub_theme + 
                   theme(legend.key.width=unit(0.1, "cm"))),
             map(ce.sig["wksSinceAntiAGD_z:temperature_z"], ~.x + 
                   annotate("text", hjust=0, vjust=0.5, size=3,
                            x=ce.sig[["wksSinceAntiAGD_z:temperature_z"]]$data$effect1__[1], y=lice_y$limits[2], 
                            label=paste(filter(eff.sigLevel, term=="wksSinceAntiAGD_z:temperature_z")$conf, "conf."), 
                            colour=filter(eff.sigLevel, term=="wksSinceAntiAGD_z:temperature_z")$col) +
                   lice_y + temp_col + pub_theme + 
                   theme(legend.key.width=unit(0.1, "cm"))),
             map(ce.sig["wksSinceMedi_z:waveHeight"], ~.x + 
                   annotate("text", hjust=0, vjust=0.5, size=3,
                            x=ce.sig[["wksSinceMedi_z:waveHeight"]]$data$effect1__[1], y=lice_y$limits[2], 
                            label=paste(filter(eff.sigLevel, term=="wksSinceMedi_z:waveHeight")$conf, "conf."), 
                            colour=filter(eff.sigLevel, term=="wksSinceMedi_z:waveHeight")$col) +
                   lice_y + wave_col + pub_theme + 
                   theme(legend.key.width=unit(0.1, "cm"))), 
             list(ggplot() + 
                    annotation_custom(grid::rasterGrob(lepFe.png)) + 
                    ggtitle("") +
                    theme(axis.line=element_blank(),
                          axis.title=element_blank(),
                          axis.text=element_blank())))
lice.80.p <- ggarrange(plotlist=plot.ls, ncol=2, nrow=4, 
                       common.legend=F, labels=c(letters[1:7], ""))
ggsave("figs/fig_6.jpg", lice.80.p, width=140, height=200, units="mm", dpi=300)



# * Summary values --------------------------------------------------------
# Conditional and marginal R2
bayes_R2(lice.fit[[1]])
bayes_R2(lice.fit[[1]], re_formula=NA)

# Time: Lice at 1 year
ce$days_since_start_z$data %>% 
  ungroup %>% 
  filter(effect1__ == min(effect1__) | abs(effect1__ - 365) < 2) %>%
  mutate(across(one_of("estimate__", "lower__", "upper__"), 
                ~.x/first(.x), 
                .names="propDiff_{.col}"))

# AGD: Lice at low vs high AGD
ce$`agd_score_z:waveHeight`$data %>% 
  ungroup %>% 
  filter(round(effect1__,2)==1.4) %>%  # expm1(1.4) = approx 3
  mutate(across(one_of("estimate__", "lower__", "upper__"), 
                ~.x/first(.x), 
                .names="propDiff_{.col}"))

# AGD*waterspeed: Lice at low vs high AGD
# Slow water current; AGD=effect2__
ce$`waterSpeed_z:agd_score_z`$data %>% 
  ungroup %>% 
  filter(effect1__==min(effect1__)) %>%
  mutate(across(one_of("estimate__", "lower__", "upper__"), 
                ~.x/first(.x), 
                .names="propDiff_{.col}"))
# Fast water current; AGD=effect2__
ce$`waterSpeed_z:agd_score_z`$data %>% 
  ungroup %>% 
  filter(effect1__==max(effect1__)) %>%
  mutate(across(one_of("estimate__", "lower__", "upper__"), 
                ~.x/first(.x), 
                .names="propDiff_{.col}"))



# fitted lines ------------------------------------------------------------

pred.df <- mort.dat %>%
  left_join(lice.dat %>% select(farm, cage, day, lep_fe))
pred.mort <- posterior_epred(mort.fit[[1]], newdata=pred.df)
pred.lice <- posterior_epred(lice.fit[[1]], newdata=pred.df)
pred.df <- pred.df %>%
  mutate(p_mort_mn=colMeans(pred.mort),
         p_mort_95l=HDInterval::hdi(pred.mort)[1,],
         p_mort_95h=HDInterval::hdi(pred.mort)[2,],
         p_lice_mn=colMeans(pred.lice),
         p_lice_95l=HDInterval::hdi(pred.lice)[1,],
         p_lice_95h=HDInterval::hdi(pred.lice)[2,]) %>%
  mutate(p_mort_mn_natural=boot::inv.logit(p_mort_mn*sd(mort.dat$lmortRate) + 
                                             mean(mort.dat$lmortRate)))


# * Fig 3 -----------------------------------------------------------------
# mortality
ggplot(pred.df, aes(day)) +
  geom_point(alpha=0.5, size=0.5, shape=1, aes(y=boot::inv.logit(lmortRate))) +
  geom_line(alpha=0.75, size=0.25, aes(y=p_mort_mn_natural, group=cage), colour="dodgerblue") +
  scale_x_date(date_breaks="6 months", date_labels="%Y-%b") +
  scale_y_continuous("Weekly mortality rate (%)", trans="logit", 
                     breaks=c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
                     labels=c(0.001, 0.01, 0.1, 1, 10),
                     limits=c(3e-6, 1e-1)) +
  facet_wrap(~farm, ncol=2) +
  theme_classic() + 
  theme(legend.key.height=unit(0.1, "cm"),
        legend.position="bottom", 
        axis.title.x=element_blank(),
        axis.title=element_text(size=9),
        axis.text=element_text(size=8),
        strip.text=element_text(size=9))
ggsave("figs/fig_3.jpg", width=140, height=140, units="mm", dpi=300)



# * Fig 5 -----------------------------------------------------------------
# sea lice
ggplot(pred.df, aes(day)) +
  geom_point(alpha=0.5, size=0.5, shape=1, aes(y=expm1(lep_fe))) +
  geom_line(alpha=0.75, size=0.25, aes(y=expm1(p_lice_mn), group=cage), colour="dodgerblue") +
  scale_x_date(date_breaks="6 months", date_labels="%Y-%b") +
  scale_y_continuous("Mean lice per fish", trans="log1p", 
                     breaks=c(0, 1, 5, 10), limits=c(-0.2, 15)) +
  facet_wrap(~farm, ncol=2) +
  theme_classic() + 
  theme(legend.key.height=unit(0.1, "cm"),
        legend.position="bottom", 
        axis.title.x=element_blank(),
        axis.title=element_text(size=9),
        axis.text=element_text(size=8),
        strip.text=element_text(size=9))
ggsave("figs/fig_5.jpg", width=140, height=140, units="mm", dpi=300)

