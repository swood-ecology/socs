# LOAD PACKAGES
library(tidyverse)    # For data manipulation
library(readxl)       # For reading in .xls files
library(VSURF)        # For random forest model selection
library(doBy)
library(GGally)       # For ggcorr
library(ggthemes)     # For alternate ggplot themes
library(grid)         # For manipulating details of ggplot
library(lme4)         # For mixed-effects regression
library(lmerTest)     # For p-value calculations of lmer models
library(MuMIn)
library(sjstats)
library(stargazer)    # For generating HTML tables for export
library(jtools)       # For interaction plots
source("fig4_function.R")


# READ IN DATA
usda_soil_data <- read_excel("data/usda-soil-data.xlsx") %>%
  mutate(
    MAOM_CN = `MAOM C` / `MAOM N`,
    POM_CN = `POM C` / `POM N`
  )

bd <- read_excel("data/site-data.xlsx",
  sheet = "Bulk Density Data",
  col_types = c("skip", "skip", "numeric", "numeric", "text", "text", "numeric", "numeric", "numeric")
) %>%
  filter(calyr == "2011") %>%
  mutate(
    rep = paste0("Rep ", rep)
  )

inputs <- read_excel("data/site-data.xlsx",
  sheet = "Org Matter and N inputs Summary",
  col_types = c(
    "skip", "skip", "text", "numeric", "numeric", "numeric",
    "numeric", "numeric", "numeric", "numeric", "numeric",
    "numeric", "numeric", "numeric", "numeric", "numeric",
    "numeric", "numeric", "numeric", "numeric"
  )
) %>%
  mutate(
    rep = paste0("Rep ", rep)
  )

icp <- read_excel("data/icp.xlsx")


# DATA MANIPULATION
usda_soil_data$organic_matter <- usda_soil_data$organic_matter*0.58
inputsLabels <- colnames(inputs)
colnames(inputs)[4:18] <- c(
  "total_compost",
  "annual_cc_shoot",
  "total_cc_shoot",
  "annual_legume_cc_shoot",
  "total_legume_cc_shoot",
  "n_fixation",
  "compost_n",
  "fertilizer_n",
  "total_n",
  "perc_n_from_fixation",
  "total_veg_residue_shoot",
  "total_om",
  "fresh_om",
  "fresh_om_perc",
  "total_C_no_roots_no_exudates"
)

usda_soil_data$slope <- ifelse(usda_soil_data$Replicate == "Rep 1" | usda_soil_data$Replicate == "Rep 2", "Upslope", "Downslope")

inputs$RepTreat <- paste(inputs$trt, inputs$rep)
usda_soil_data$RepTreat <- paste(usda_soil_data$Trt_id, usda_soil_data$Replicate)
bd$RepTreat <- paste(bd$trt, bd$rep)
icp$RepTreat <- paste(icp$Trt_id, icp$Replicate)

all_data <- merge(usda_soil_data, bd, by = "RepTreat") %>%
  select(RepTreat, System, Replicate, Trt_id, cover_crop, Compost:slope, blkden) %>%
  mutate(
    pom.stock = (`POM C`*blkden*30)/10,
    pom.n.stock = (`POM N`*blkden*30)/10,
    maom.stock = (`MAOM C`*blkden*30)/10,
    maom.n.stock = (`MAOM N`*blkden*30)/10,
    totC.stock = pom.stock + maom.stock
  )

all_data <- all_data %>%
  merge(inputs, by = "RepTreat") %>%
  select(RepTreat,System,Replicate:totC.stock, total_compost:total_C_no_roots_no_exudates) 

all_data$CC <- recode(all_data$Trt_id,
                      leg3x="Legume-rye every year",
                      mus1x="Mustard every year",
                      nocc="Legume-rye every four years",
                      noccnocp="Legume-rye every four years",
                      rye1x="Rye every year")

all_data <- all_data %>%
  merge(icp %>% select(EC:RepTreat), by="RepTreat")

rm(bd); rm(inputs); rm(usda_soil_data); rm(icp)


# DATA EXPLORATION AND PLOTTING
## Table 1. Summary stats
all_data %>%
  select(cover_crop:POM_CN,blkden,totC.stock,maom.stock,pom.stock,maom.n.stock,pom.n.stock) %>%
  group_by(cover_crop,Compost,Cover_crop_freq) %>%
  summarise(
    mean_sand=mean(perc_sand),
    sd_sand=sd(perc_sand),
    mean_silt=mean(perc_silt),
    sd_silt=sd(perc_silt),
    mean_clay=mean(perc_clay),
    sd_clay=sd(perc_clay),
    mean_pH=mean(pH),
    sd_pH=sd(pH),
    mean_blkden=mean(blkden),
    sd_blkden=sd(blkden),
    mean_om=mean(organic_matter),
    sd_om=sd(organic_matter),
    mean_totC=mean(totC.stock),
    sd_totC=sd(totC.stock),
    mean_POMC=mean(`POM C`),
    sd_POMC=sd(`POM C`),
    mean_POMC_stock=mean(pom.stock),
    sd_POMC_stock=sd(pom.stock),
    mean_POMN=mean(`POM N`),
    sd_POMN=sd(`POM N`),
    mean_POMN_stock=mean(pom.n.stock),
    sd_POMN_stock=sd(pom.n.stock),
    mean_POMCN=mean(`POM_CN`),
    sd_POMCN=sd(`POM_CN`),
    mean_MAOMC=mean(`MAOM C`),
    sd_MAOMC=sd(`MAOM C`),
    mean_MAOMC_stock=mean(`maom.stock`),
    sd_MAOMC_stock=sd(`maom.stock`),
    mean_MAOMN=mean(`MAOM N`),
    sd_MAOMN=sd(`MAOM N`),
    mean_MAOMN_stock=mean(`maom.n.stock`),
    sd_MAOMN_stock=sd(`maom.n.stock`),
    mean_MAOMCN=mean(`MAOM_CN`),
    sd_MAOMCN=sd(`MAOM_CN`),
    mean_WHC=mean(water_holding_capacity, na.rm=T),
    sd_WHC=sd(water_holding_capacity, na.rm=T),
    mean_SIR=mean(substrate_induced_respiration),
    sd_SIR=sd(substrate_induced_respiration)
  ) %>%
  write_excel_csv("tables/table1.csv")

all_data %>%
  filter(cover_crop=="legume-rye") %>%
  select(cover_crop:POM_CN,blkden,totC.stock,maom.stock,pom.stock,maom.n.stock,pom.n.stock) %>%
  group_by(Compost,Cover_crop_freq) %>%
  summarise(
    mean_sand=mean(perc_sand),
    sd_sand=sd(perc_sand),
    mean_silt=mean(perc_silt),
    sd_silt=sd(perc_silt),
    mean_clay=mean(perc_clay),
    sd_clay=sd(perc_clay),
    mean_pH=mean(pH),
    sd_pH=sd(pH),
    mean_blkden=mean(blkden),
    sd_blkden=sd(blkden),
    mean_om=mean(organic_matter),
    sd_om=sd(organic_matter),
    mean_totC=mean(totC.stock),
    sd_totC=sd(totC.stock),
    mean_POMC=mean(`POM C`),
    sd_POMC=sd(`POM C`),
    mean_POMC_stock=mean(pom.stock),
    sd_POMC_stock=sd(pom.stock),
    mean_POMN=mean(`POM N`),
    sd_POMN=sd(`POM N`),
    mean_POMN_stock=mean(pom.n.stock),
    sd_POMN_stock=sd(pom.n.stock),
    mean_POMCN=mean(`POM_CN`),
    sd_POMCN=sd(`POM_CN`),
    mean_MAOMC=mean(`MAOM C`),
    sd_MAOMC=sd(`MAOM C`),
    mean_MAOMC_stock=mean(`maom.stock`),
    sd_MAOMC_stock=sd(`maom.stock`),
    mean_MAOMN=mean(`MAOM N`),
    sd_MAOMN=sd(`MAOM N`),
    mean_MAOMN_stock=mean(`maom.n.stock`),
    sd_MAOMN_stock=sd(`maom.n.stock`),
    mean_MAOMCN=mean(`MAOM_CN`),
    sd_MAOMCN=sd(`MAOM_CN`),
    mean_WHC=mean(water_holding_capacity, na.rm=T),
    sd_WHC=sd(water_holding_capacity, na.rm=T),
    mean_SIR=mean(substrate_induced_respiration),
    sd_SIR=sd(substrate_induced_respiration)
  ) %>% select(mean_om, sd_om) %>%
  write_excel_csv("tables/table1_filtered.csv")


## Fig. 1. Cumulative inputs
input.colors <- c("#fff7bc","#662506")

all_data %>%
  select(System, Compost, 
         "Nitrogen" = total_n,
         "Organic matter" = total_om,
         "Carbon" = total_C_no_roots_no_exudates) %>%
  gather(-System, -Compost, key="var", value="value") -> fig1.data

error <- summaryBy(value ~ var + Compost + System,
                   fig1.data,
                   FUN = c(mean, se)
)

ggplot(fig1.data, aes(x=System,y=value,fill=Compost)) +
  geom_bar(position = position_dodge(preserve = "single"), 
           stat = "summary", fun.y = "mean") +
  scale_fill_manual(values=input.colors) +
  geom_errorbar(data=error, mapping = aes(
    x = System,
    ymin = value.mean - value.se,
    ymax = value.mean + value.se  
  ), width = .25,
  position = position_dodge(0.9, preserve = "single"),
  inherit.aes = FALSE) +
  facet_wrap(~var, scales = "free",
             strip.position = "left",
             labeller = as_labeller(c(
               Carbon="Total C input, excl. roots and exudates (kg ha-1)",
               Nitrogen="Total N input (kg ha-1)",
               `Organic matter`="Total organic inputs (kg ha-1)"
             ))
  ) +
  theme_light() + xlab("") + ylab("") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color="black"),
    strip.placement = "outside",
    axis.text.x = element_text(size=8,angle=90)
  ) -> fig1

# Get gtable object
g <- ggplotGrob(fig1)
# Add an extra top row (make some space)
g <- gtable::gtable_add_rows(x = g, heights = unit(0.65, 'cm'), pos = 2)
# First strip
g <- gtable::gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = "grey70")),
                                  textGrob(label = "Carbon", 
                                           gp = gpar(col = "white"))),
                     t = 3, l = 7, b = 3, r = 7, 
                     name = c("strip-top-1-rectg", "strip-top-1-text"))
# Second strip
g <- gtable::gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = "grey70")),
                                  textGrob(label = "Nitrogen", 
                                           gp = gpar(col = "white"))),
                     t = 3, l = 13, b = 3, r = 13, 
                     name = c("strip-top-2-rectg", "strip-top-2-text"))
# Third strip
g <- gtable::gtable_add_grob(x = g,
                             grobs = list(rectGrob(gp = gpar(col = NA, 
                                                             fill = "grey70")),
                                          textGrob(label = "Organic matter", 
                                                   gp = gpar(col = "white"))),
                             t = 3, l = 19, b = 3, r = 19, 
                             name = c("strip-top-3-rectg", "strip-top-3-text"))

# Draw the edited plot
grid.newpage()
grid.draw(g)

rm(g); rm(fig1)


## Fig 2 model
fig2model <- all_data %>%
  filter(cover_crop=="legume-rye") %>%
  select(
    Replicate, cover_crop, Compost, Cover_crop_freq, organic_matter, water_holding_capacity,
    substrate_induced_respiration:POM_CN, `POM C`, `MAOM C`, slope, pom.stock, pom.n.stock, maom.stock, maom.n.stock, 
    total_C_no_roots_no_exudates, CEC, pH, `X-Ca`
  )

lmer(organic_matter ~ Compost + Cover_crop_freq + (1|Replicate),
  data = fig2model
) -> om.model
summary(om.model)
r.squaredGLMM(om.model)

lmer(pom.stock ~ Compost + Cover_crop_freq + (1|Replicate),
     data = fig2model
) -> pom.model
summary(pom.model)
r.squaredGLMM(pom.model)

lmer(pom.n.stock ~ Compost + Cover_crop_freq + (1|Replicate),
     data = fig2model
) -> pom.n.model
summary(pom.n.model)
r.squaredGLMM(pom.n.model)

lmer(maom.stock ~ Compost + Cover_crop_freq + (1|Replicate),
     data = fig2model
) -> maom.model
summary(maom.model)
r.squaredGLMM(maom.model)

lmer(maom.n.stock ~ Compost + Cover_crop_freq + (1|Replicate),
     data = fig2model
) -> maom.n.model
summary(maom.n.model)
r.squaredGLMM(maom.n.model)

all_data %>%
  select(
    Replicate, cover_crop, Compost, Cover_crop_freq, organic_matter, water_holding_capacity,
    substrate_induced_respiration:POM_CN, `POM C`, `MAOM C`, slope, pom.stock, pom.n.stock, maom.stock, maom.n.stock, 
    total_C_no_roots_no_exudates, CEC, pH, `X-Ca`
  ) -> fig2model.sir

lmer(substrate_induced_respiration ~ cover_crop + Compost + Cover_crop_freq + (1|Replicate),
     data = fig2model.sir
) -> sir.model
summary(sir.model)
r.squaredGLMM(sir.model)

lmer(substrate_induced_respiration ~ cover_crop + Compost + Cover_crop_freq + pom.n.stock + (1|Replicate),
     data = fig2model.sir
) -> sir.model.cn
car::vif(sir.model.cn)
summary(sir.model.cn)
r.squaredGLMM(sir.model.cn)

class(om.model) <- "lmerMod"
class(pom.model) <- "lmerMod"
class(pom.n.model) <- "lmerMod"
class(maom.model) <- "lmerMod"
class(maom.n.model) <- "lmerMod"
class(sir.model) <- "lmerMod"
class(sir.model.cn) <- "lmerMod"

stargazer(om.model, maom.model, maom.n.model, pom.model, pom.n.model, sir.model, sir.model.cn, out="regression.htm")


## Fig. 4 MAOM by soil properties
# Determine best management predictors of MAOM
mgmt <- VSURF(
  y = all_data$maom.stock,
  x = all_data %>%
    select(total_compost:total_C_no_roots_no_exudates),
  parallel = TRUE
)

# Fit model with best management predictors
mgmt.model <- lm(maom.stock ~ Compost+Cover_crop_freq+total_veg_residue_shoot + annual_legume_cc_shoot, data=all_data)
mgmt.model %>% summary()

# Find best non-management predictors of model residuals
soil.prop <- VSURF(
  y = mgmt.model$residuals,
  x = all_data %>%
    select(perc_sand:perc_clay,EC:`Fe (DTPA)`),
  parallel = TRUE
)

# Create data frame of resulting variables, also including clay
resid.data <- data.frame(mgmt.model$residuals,all_data$perc_clay,all_data$`X-K`,all_data$`Ca (SP)`,all_data$`Mn (DTPA)`)
names(resid.data) <- c('model_residuals','perc_clay','K','Ca','Mn')

resid.model <- lm(model_residuals ~ perc_clay + K + Ca + Mn, data=resid.data)
resid.model %>% summary()

effect_plot(model=resid.model, pred = Ca, colors=all_data$System, 
            shape=all_data$Replicate, 
            x.label = "Calcium", y.label = "First-stage model residuals")
effect_plot(model=resid.model, pred = K, colors=all_data$System, 
            shape=all_data$Replicate, 
            x.label = "Potassium", y.label = "First-stage model residuals")
effect_plot(model=resid.model, pred = perc_clay, colors=all_data$System, 
            shape=all_data$Replicate, 
            x.label = "Percent clay", y.label = "First-stage model residuals")
effect_plot(model=resid.model, pred = Mn, colors=all_data$System, 
            shape=all_data$Replicate, 
            x.label = "Manganese", y.label = "First-stage model residuals")

stargazer(mgmt.model, resid.model, out="maom.htm")

# Fig. 5 Soil properties by quantitative inputs
## Random forest models of quantitative predictors
rf.res <- list(
  VSURF(
    y = all_data$organic_matter,
    x = all_data %>%
      select(total_compost:total_C_no_roots_no_exudates),
    parallel = TRUE
  ), # total_C_no_roots_no_exudates, fresh_om, total_om, total_cc_shoot, fresh_om_perc, total_compost
  VSURF(
    y = all_data$pom.stock,
    x = all_data %>%
      select(total_compost:total_C_no_roots_no_exudates),
    parallel = TRUE
  ), # fresh_om
  VSURF(
    y = all_data$maom.stock,
    x = all_data %>%
      select(total_compost:total_C_no_roots_no_exudates),
    parallel = TRUE
  ), # total_veg_residue_shoot, annual_legume_cc_shoot
  VSURF(
    y = all_data$substrate_induced_respiration,
    x = all_data %>%
      select(total_compost:total_C_no_roots_no_exudates),
    parallel = TRUE
  ) # annual_cc_shoot
)

### organic_matter
# total_C_no_roots_no_exudates, total_cc_shoot, fresh_om_perc, total_compost
lm.om <- lmer(organic_matter ~ total_C_no_roots_no_exudates + total_om + total_cc_shoot + fresh_om_perc + (1|Replicate), data=all_data)
summary(lm.om)
r.squaredGLMM(lm.om)
lm.om <- lm(organic_matter ~ total_C_no_roots_no_exudates + total_om + total_cc_shoot + fresh_om_perc, data=all_data)
summary(lm.om)

### pom.stock x fresh_om
lm.pom <- lmer(pom.stock ~ fresh_om + (1|Replicate), data=all_data)
summary(lm.pom)
r.squaredGLMM(lm.pom)
lm.pom <- lm(pom.stock ~ fresh_om, data=all_data)
summary(lm.pom)

### maom.stock 
# total_veg_residue_shoot, annual_legume_cc_shoot
lm.maom <- lmer(maom.stock ~ total_veg_residue_shoot + annual_legume_cc_shoot + (1|Replicate), data=all_data)
r.squaredGLMM(lm.maom)
lm.maom <- lm(maom.stock ~ total_veg_residue_shoot + annual_legume_cc_shoot, data=all_data)
summary(lm.maom)
car::vif(lm.maom)

### SIR x annual_cc_shoot
lm.sir <- lm(substrate_induced_respiration ~ annual_cc_shoot, data=all_data)
summary(lm.sir)

class(lm.om) <- "lmerMod"
class(lm.pom) <- "lmerMod"
class(lm.maom) <- "lmerMod"
class(lm.sir) <- "lmerMod"

stargazer(lm.om, lm.maom, lm.pom, lm.sir, out="regression2.htm")



###################################################################################

###########################
##### DEPRECATED CODE #####
### MAINLY FOR PLOTTING ###
###########################

# ## Fig. S1 POM & MAOM C by cover crop type
# all_data %>%
#   filter(Cover_crop_freq=="Annually") %>%
#   select(cover_crop,pom.stock,pom.n.stock,maom.stock,maom.n.stock,`substrate_induced_respiration`) %>%
#   gather(-cover_crop,key="var",value="value")-> figs1.data
# 
# maom <- figs1.data %>% 
#   filter(var=="maom.stock") %>%
#   ggplot(aes(x = cover_crop, y = value)) +
#   geom_boxplot() +
#   geom_point(size = 1.5) +
#   xlab("") + ylab("Mineral-associated C (Mg ha-1)") +
#   theme_light() +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# pom <- figs1.data %>% 
#   filter(var=="pom.stock") %>%
#   ggplot(aes(x = cover_crop, y = value)) +
#   geom_boxplot() +
#   geom_point(size = 1.5) +
#   xlab("") + ylab("Particulate C (Mg ha-1)") +
#   theme_light() +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# sir <- figs1.data %>% 
#   filter(var=="substrate_induced_respiration") %>%
#   ggplot(aes(x = cover_crop, y = value)) +
#   geom_boxplot() +
#   geom_point(size = 1.5) +
#   xlab("") + ylab("SIR (ug C g soil-1 h-1)") +
#   theme_light() +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )



# ## Fig. S2 Multi-panel of POM & MAOM N by all treatments
# all_data %>%
#   mutate(
#     "Particulate nitrogen" = pom.n.stock,
#     "Mineral-associated nitrogen" = maom.n.stock
#   ) -> figs1.data
# figs2.data$Compost <-
#   recode(figs1.data$Compost,
#          "No" = "No compost",
#          "Yes" = "Compost")
# figs2.data %>%  
#   select(Compost, cover_crop, Cover_crop_freq, `Particulate nitrogen`:`Mineral-associated nitrogen`) %>%
#   gather(-Compost, -cover_crop, -Cover_crop_freq, key = "var", value = "value") -> figs1.data
# 
# bar.colors <- c("#d9f0d3", "#f6e8c3")
# point.colors <- c("#1b7837", "#8c510a")
# 
# pomN <- figs1.data %>% 
#   filter(var == "Particulate nitrogen") %>%
#   ggplot(aes(x = Compost, y = value, fill = Cover_crop_freq)) +
#   geom_boxplot(size=0.25,outlier.shape=NA,position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, name = "Cover Crop Frequency") +
#   geom_point(aes(color = Cover_crop_freq, shape = cover_crop),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Particulate nitrogen") +
#   scale_color_manual(values = point.colors, name = "Cover Crop Frequency") +
#   xlab("") + ylab("Particulate N (Mg ha-1 to 30 cm)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title = "Cover Crop Frequency",
#       title.position = "top",
#       title.hjust = 0.5
#     ),
#     shape = guide_legend(
#       title = "Cover Crop Type",
#       title.position = "top",
#       title.hjust = 0.5
#     )
#   ) +
#   scale_shape(guide = FALSE) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# 
# maomN <- figs1.data %>% 
#   filter(var == "Mineral-associated nitrogen") %>%
#   ggplot(aes(x = Compost, y = value, fill = Cover_crop_freq)) +
#   geom_boxplot(size=0.25,outlier.shape=NA,position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, name = "Cover Crop Frequency") +
#   geom_point(aes(color = Cover_crop_freq, shape = cover_crop),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Mineral-associated nitrogen") +
#   scale_color_manual(values = point.colors, name = "Cover Crop Frequency") +
#   xlab("") + ylab("Mineral-associated N (Mg ha-1 to 30 cm)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title = "Cover Crop Frequency",
#       title.position = "top",
#       title.hjust = 0.5
#     ),
#     shape = guide_legend(
#       title = "Cover Crop Type",
#       title.position = "top",
#       title.hjust = 0.5
#     )
#   ) +
#   scale_shape(guide = FALSE) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )

# ## Fig. S3 Multi-panel of POM, MAOM, organic_matter, and SIR by all treatments
# all_data %>%
#   mutate(
#     "Particulate carbon" = pom.stock,
#     "Mineral-associated carbon" = maom.stock,
#     "Total organic matter" = organic_matter,
#     "Substrate-induced respiration" = substrate_induced_respiration
#   ) -> figs3.data
# figs3.data$Compost <-
#   recode(figs1.data$Compost,
#          "No" = "No compost",
#          "Yes" = "Compost")
# figs3.data %>%  
#   select(Compost, cover_crop, Cover_crop_freq, `Particulate carbon`:`Substrate-induced respiration`) %>%
#   gather(-Compost, -cover_crop, -Cover_crop_freq, key = "var", value = "value") -> figs1.data
# 
# bar.colors <- c("#d9f0d3", "#f6e8c3")
# point.colors <- c("#1b7837", "#8c510a")
# 
# pom <- figs3.data %>% 
#   filter(var == "Particulate carbon") %>%
#   ggplot(aes(x = cover_crop, y = value, fill = Compost)) +
#   geom_boxplot(size=0.25,
#                outlier.shape=NA,
#                position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, 
#                     name = "Compost") +
#   geom_point(aes(color = Compost, shape = Cover_crop_freq),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Particulate carbon") +
#   scale_color_manual(values = point.colors, 
#                      name = "Compost") +
#   scale_shape(name="Cover crop frequency") +
#   xlab("") + ylab("Particulate C (Mg ha-1 to 30 cm)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title.position = "top",
#       title.hjust = 0.5
#     ),
#     shape = guide_legend(
#       title.position = "top",
#       title.hjust = 0.5
#     )
#   ) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# 
# maom <- figs3.data %>% 
#   filter(var == "Mineral-associated carbon") %>%
#   ggplot(aes(x = cover_crop, y = value, fill = Compost)) +
#   geom_boxplot(size=0.25,
#                outlier.shape=NA,
#                position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, 
#                     name = "Compost") +
#   geom_point(aes(color = Compost, shape = Cover_crop_freq),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Mineral-associated carbon") +
#   scale_color_manual(values = point.colors, 
#                      name = "Compost") +
#   scale_shape(name="Cover crop frequency") +
#   xlab("") + ylab("Mineral-associated C (Mg ha-1 to 30 cm)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title.position = "top",
#       title.hjust = 0.5
#     ),
#     shape = guide_legend(
#       title.position = "top",
#       title.hjust = 0.5
#     )
#   ) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# 
# tom <- figs3.data %>% 
#   filter(var == "Total organic matter") %>%
#   ggplot(aes(x = cover_crop, y = value, fill = Compost)) +
#   geom_boxplot(size=0.25,
#                outlier.shape=NA,
#                position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, 
#                     name = "Compost") +
#   geom_point(aes(color = Compost, shape = Cover_crop_freq),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Total organic matter") +
#   scale_color_manual(values = point.colors, 
#                      name = "Compost") +
#   scale_shape(name="Cover crop frequency") +
#   xlab("") + ylab("Total organic matter (%)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title.position = "top",
#       title.hjust = 0.5
#     ),
#     shape = guide_legend(
#       title.position = "top",
#       title.hjust = 0.5
#     )
#   ) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# 
# sir <- figs3.data %>% 
#   filter(var == "Substrate-induced respiration") %>%
#   ggplot(aes(x = cover_crop, y = value, fill = Compost)) +
#   geom_boxplot(size=0.25,
#                outlier.shape=NA,
#                position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, 
#                     name = "Compost") +
#   geom_point(aes(color = Compost, shape = Cover_crop_freq),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Substrate-induced respiration") +
#   scale_color_manual(values = point.colors, 
#                      name = "Compost") +
#   scale_shape(name="Cover crop frequency") +
#   xlab("") + ylab("Substrate-induced respiration (ug C g soil-1 h-1)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title.position = "top",
#       title.hjust = 0.5
#     ),
#     shape = guide_legend(
#       title.position = "top",
#       title.hjust = 0.5
#     )
#   ) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )


# ## Fig. S4 Plot correlations among quantitative soil variables
# ## Define the correlation matrix
# fig6.data <- all_data %>% 
#   select(organic_matter,perc_sand,perc_clay,pH,pom.stock,maom.stock,EC,`NO3-N`:`X-K`,`X-Ca`:`Zn (DTPA)`,CEC)
# 
# names(fig6.data) <- c("Organic matter","Sand","Clay","pH","Particulate C",
#                       "Mineral C","Electric conductivity","Nitrate","Olsen's P","K","Ca","Mg","Zn","CEC")
# 
# ggcorr(fig6.data,
#        nbreaks=6,
#        hjust = 0.75, 
#        size = 3, 
#        label=T,
#        label_size=4,
#        label_alpha=T,
#        layout.exp = 3) +
#   theme(legend.title = element_text(size = 7),
#         legend.text = element_text(size = 7)) 
# 
## OTHER CORRELATION FIGURE TYPES
# ggcorr(fig6.data, nbreaks=5)
# ggcorr(fig6.data, geom = "circle", nbreaks = 5,min_size=0, max_size=6)
# ggcorr(fig6.data,
#        hjust = 0.75, 
#        size = 4, 
#        label=T,
#        label_size=2,
#        label_alpha=T,
#        geom="circle", 
#        min_size=0, max_size=6,
#        name = expression(rho), 
#        legend.position = "bottom", 
#        legend.size = 12,
#        layout.exp = 1) +
#   guides(fill = guide_colorbar(barwidth = 18, title.vjust = 0.75)) +
#   theme(legend.title = element_text(size = 14))
# 
# ggcorr(fig6.data, 
#        geom = "blank", 
#        label = TRUE, 
#        hjust = 0.75) +
#   geom_point(size = 10, aes(color = coefficient < 0, alpha = abs(coefficient) > 0.25)) +
#   scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
#   guides(color = FALSE, alpha = FALSE)

# ## Fig. 2 Multi-panel of POM, MAOM, organic_matter, and SIR by all treatments
# all_data %>%
#   filter(cover_crop == "legume-rye") %>%
#   mutate(
#     "Particulate carbon" = pom.stock,
#     "Mineral-associated carbon" = maom.stock,
#     "Total organic matter" = organic_matter,
#     "Substrate-induced respiration" = substrate_induced_respiration
#   ) -> fig2.data
# fig2.data$Compost <-
#   recode(fig2.data$Compost,
#          "No" = "No compost",
#          "Yes" = "Compost")
# fig2.data %>%  
#   select(Compost, Cover_crop_freq, `Particulate carbon`:`Substrate-induced respiration`) %>%
#   gather(-Compost, -Cover_crop_freq, key = "var", value = "value") -> fig2.data
# 
# bar.colors <- c("#d9f0d3", "#f6e8c3")
# point.colors <- c("#1b7837", "#8c510a")
# 
# pom <- fig2.data %>% 
#   filter(var == "Particulate carbon") %>%
#   ggplot(aes(x = Compost, y = value, fill = Cover_crop_freq)) +
#   geom_boxplot(size=0.25,outlier.shape=NA,position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, name = "Cover Crop Frequency") +
#   geom_point(aes(color = Cover_crop_freq, shape = Compost),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Particulate carbon") +
#   scale_color_manual(values = point.colors, name = "Cover Crop Frequency") +
#   xlab("") + ylab("Particulate C (Mg ha-1 to 30 cm)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title = "Cover Crop Frequency",
#       title.position = "top",
#       title.hjust = 0.5,
#       override.aes = list(shape = c(16, 17))
#     )
#   ) +
#   scale_shape(guide = FALSE) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# 
# maom <- fig2.data %>% 
#   filter(var == "Mineral-associated carbon") %>%
#   ggplot(aes(x = Compost, y = value, fill = Cover_crop_freq)) +
#   geom_boxplot(size=0.25,outlier.shape=NA,position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, name = "Cover Crop Frequency") +
#   geom_point(aes(color = Cover_crop_freq, shape = Compost),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Mineral-associated carbon") +
#   scale_color_manual(values = point.colors, name = "Cover Crop Frequency") +
#   xlab("") + ylab("Mineral-associated C (Mg ha-1 to 30 cm)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title = "Cover Crop Frequency",
#       title.position = "top",
#       title.hjust = 0.5,
#       override.aes = list(shape = c(16, 17))
#     )
#   ) +
#   scale_shape(guide = FALSE) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# 
# tom <- fig2.data %>% 
#   filter(var == "Total organic matter") %>%
#   ggplot(aes(x = Compost, y = value, fill = Cover_crop_freq)) +
#   geom_boxplot(size=0.25,outlier.shape=NA,position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, name = "Cover Crop Frequency") +
#   geom_point(aes(color = Cover_crop_freq, shape = Compost),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Total organic matter") +
#   scale_color_manual(values = point.colors, name = "Cover Crop Frequency") +
#   xlab("") + ylab("Total organic matter (%)\n") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title = "Cover Crop Frequency",
#       title.position = "top",
#       title.hjust = 0.5,
#       override.aes = list(shape = c(16, 17))
#     )
#   ) +
#   scale_shape(guide = FALSE) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   ) 
# 
# sir <- fig2.data %>% filter(var == "Substrate-induced respiration") %>% 
#   ggplot(aes(x = Compost, y = value, fill = Cover_crop_freq)) +
#   geom_boxplot(size=0.25,outlier.shape=NA,position = position_dodge2(preserve = "single")) +  
#   scale_fill_manual(values = bar.colors, name = "Cover Crop Frequency") +
#   geom_point(aes(color = Cover_crop_freq, shape = Compost),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   facet_grid(. ~ "Substrate-induced respiration") +  
#   scale_color_manual(values = point.colors, name = "Cover Crop Frequency") +
#   xlab("") + ylab("Substrate-induced respiration (ug C h-1 g soil-1)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title = "Cover Crop Frequency",
#       title.position = "top",
#       title.hjust = 0.5,
#       override.aes = list(shape = c(16, 17))
#     )
#   ) +
#   scale_shape(guide = FALSE) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# 
# # Fig 3: MAOM C & POM C on same scale
# all_data %>%
#   filter(cover_crop == "legume-rye") %>%
#   mutate(
#     "Particulate carbon" = pom.stock,
#     "Mineral-associated carbon" = maom.stock
#   ) -> fig3.data
# fig3.data$Compost <-
#   recode(fig3.data$Compost,
#          "No" = "No compost",
#          "Yes" = "Compost")
# fig3.data %>%
#   select(Compost, Cover_crop_freq, `Particulate carbon`,`Mineral-associated carbon`) %>%
#   gather(-Compost, -Cover_crop_freq, key = "var", value = "value") -> fig3.data
# 
# fig3.data %>%
#   ggplot(aes(x = Compost, y = value, fill = Cover_crop_freq)) +
#   geom_boxplot(size=0.25,outlier.shape=NA,position = position_dodge2(preserve = "single")) +
#   scale_fill_manual(values = bar.colors, name = "Cover Crop Frequency") +
#   geom_point(aes(color = Cover_crop_freq, shape = Compost),
#              position=position_jitterdodge(),
#              size = 1.5) +
#   scale_color_manual(values = point.colors, name = "Cover Crop Frequency") +
#   facet_wrap(~var) + xlab("") + ylab("Carbon stock (Mg C ha-1 to 30 cm)") +
#   theme_light() +
#   guides(
#     fill = guide_legend(
#       title = "Cover Crop Frequency",
#       title.position = "top",
#       title.hjust = 0.5,
#       override.aes = list(shape = c(16, 17))
#     )
#   ) +
#   scale_shape(guide = FALSE) +
#   theme(
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 11),
#     strip.text.x = element_text(size = 12, face = "bold"),
#     legend.position = "bottom"
#   )
# #### Plot pairs for organic matter variables ####
# ### Function to generate subsetted data frame from random forest list
# pairs.rf <- function(data,resp,rf,n) {
#   data %>%
#     select(resp,total_compost:total_C_no_roots_no_exudates) -> temp
#   temp[, c(
#     1,
#     (rf[[n]]$varselect.pred %>%
#        as.numeric()
#      + 1)
#   )     ] %>%
#     return()
# }
# 
# pairs.rf(data=all_data,
#     resp="organic_matter",
#          rf=rf.res,
#          n=1) %>%
#   rename("Soil organic matter"=organic_matter,
#          "Total carbon input"=total_C_no_roots_no_exudates,
#          "Total compost input"=total_compost) %>%
#   ggscatmat() + 
#   theme_light()
# 
# pairs.rf(data=all_data,
#          resp="pom.stock",
#          rf=rf.res,
#          n=4) %>%
#   rename("Fresh organic matter input"=fresh_om,
#          "Particulate carbon"=pom.stock) %>%
#   ggscatmat() +
#   theme_light()

#### ####

coef(lm.om)

a = format(coef(lm(organic_matter ~ tom , data=all_data))[1], digits = 2)
b = format(coef(lm(organic_matter ~ tom , data=all_data))[2], digits = 1)
r2 = format(summary(lm(organic_matter ~ tom , data=all_data))$r.squared, digits = 2)

all_data %>%
  ggplot(aes(x=tom, y=organic_matter)) + 
  geom_point(aes(shape=Cover_crop_freq),size=2.5, alpha=0.75) + 
  geom_abline(intercept=1.902288611, slope=0.04871383)+
  annotate("text", x = 12, y = 3.2, label = paste0('y = ',a,' + ',b,'*x ',', R',"^",'2 = ',r2)) +
  ylab("Organic matter (%)") + xlab("Annualized organic inputs (t ha-1 y-1)") +
  theme_bw() +
  guides(shape=guide_legend(title="Cover Crop Frequency")) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 11),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = c(0.2, 0.85)
  ) -> Fig5a
Fig5a