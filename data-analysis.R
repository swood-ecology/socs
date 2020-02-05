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
library(simr)         # For power analysis
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
    sd_WHC=sd(water_holding_capacity, na.rm=T)
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
    sd_WHC=sd(water_holding_capacity, na.rm=T)
  ) %>% select(mean_om, sd_om) %>%
  write_excel_csv("tables/table1_filtered.csv")


## Fig. 1. Cumulative inputs
input.colors <- c("#e6b020","#611C16")

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
    ymin = value.mean - value.se*3.18,
    ymax = value.mean + value.se*3.18  
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
    axis.title.y = element_text(size=13),
    axis.text.y = element_text(size=11),
    axis.text.x = element_text(size=12,angle=90)
  ) -> fig1

# Get gtable object
g <- ggplotGrob(fig1)
# Add an extra top row (make some space)
g <- gtable::gtable_add_rows(x = g, heights = unit(0.65, 'cm'), pos = 2)
# First strip
g <- gtable::gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = "grey30")),
                                  textGrob(label = "Carbon", 
                                           gp = gpar(col = "white"))),
                     t = 3, l = 7, b = 3, r = 7, 
                     name = c("strip-top-1-rectg", "strip-top-1-text"))
# Second strip
g <- gtable::gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = "grey30")),
                                  textGrob(label = "Nitrogen", 
                                           gp = gpar(col = "white"))),
                     t = 3, l = 13, b = 3, r = 13, 
                     name = c("strip-top-2-rectg", "strip-top-2-text"))
# Third strip
g <- gtable::gtable_add_grob(x = g,
                             grobs = list(rectGrob(gp = gpar(col = NA, 
                                                             fill = "grey30")),
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

class(om.model) <- "lmerMod"
class(pom.model) <- "lmerMod"
class(pom.n.model) <- "lmerMod"
class(maom.model) <- "lmerMod"
class(maom.n.model) <- "lmerMod"

stargazer(om.model, maom.model, maom.n.model, pom.model, pom.n.model, out="regression.htm")

#####POWER ANALYSIS#####
fig2model %>%
  select(Replicate:organic_matter,pom.stock,maom.stock) %>%
  mutate(Rep = case_when(Replicate=='Rep 1' ~ 1,
                         Replicate=='Rep 2' ~ 2,
                         Replicate=='Rep 3' ~ 3,
                         Replicate=='Rep 4' ~ 4
  )) %>%
  fastDummies::dummy_cols() -> power.data

lmer(organic_matter ~ Cover_crop_freq_Annually + Compost_Yes + (1|Rep),
     data = power.data
) -> model1
model2 <- extend(model1,along='Rep',n=15)
pC1 <- powerCurve(model2,along='Rep')
plot(pC1)

lmer(organic_matter ~ Compost_Yes + Cover_crop_freq_Annually + (1|Rep),
     data = power.data
) -> model3
model4 <- extend(model3,along='Rep',n=15)
pC2 <- powerCurve(model4,along='Rep')
plot(pC2)

lmer(maom.stock ~ Cover_crop_freq_Annually + Compost_Yes + (1|Rep),
     data = power.data
) -> model5
model6 <- extend(model5,along='Rep',n=15)
pC3 <- powerCurve(model6,along='Rep')
plot(pC3)

lmer(maom.stock ~ Compost_Yes + Cover_crop_freq_Annually + (1|Rep),
     data = power.data
) -> model7
model8 <- extend(model7,along='Rep',n=15)
pC4 <- powerCurve(model8,along='Rep')
plot(pC4)

lmer(pom.stock ~ Cover_crop_freq_Annually + Compost_Yes + (1|Rep),
     data = power.data
) -> model9
model10 <- extend(model9,along='Rep',n=15)
pC5 <- powerCurve(model10,along='Rep')
plot(pC5)

lmer(pom.stock ~ Compost_Yes + Cover_crop_freq_Annually +  (1|Rep),
     data = power.data
) -> model11
model12 <- extend(model11,along='Rep',n=15)
pC6 <- powerCurve(model12,along='Rep')
plot(pC6)


## Fig. 4 MAOM by soil properties
# Determine best management predictors of MAOM
mgmt <- VSURF(
  y = all_data$maom.stock,
  x = all_data %>%
    select(total_compost:total_C_no_roots_no_exudates),
  parallel = TRUE
)

# Fit model with best management predictors
mgmt.model <- lm(maom.stock ~ Compost + Cover_crop_freq + total_veg_residue_shoot + annual_legume_cc_shoot, data=all_data)
mgmt.model %>% summary()

ca.model <- lm(`Ca (SP)`~ Compost, data=all_data)
k.model <- lm(`X-K...19`~ Compost, data=all_data)
mn.model <- lm(`Mn (DTPA)`~ Compost, data=all_data)

# Find best non-management predictors of model residuals
soil.prop <- VSURF(
  y = mgmt.model$residuals,
  x = all_data %>%
    select(perc_sand:perc_clay,EC:`Fe (DTPA)`),
  parallel = TRUE
)

# Create data frame of resulting variables, also including clay
resid.data <- data.frame(mgmt.model$residuals,all_data$perc_clay,all_data$`X-K...19`,all_data$`Ca (SP)`,all_data$`Mn (DTPA)`,all_data$`Na (SP)`)
names(resid.data) <- c('model_residuals','perc_clay','K','Ca','Mn','Na')

resid.model <- lm(model_residuals ~ perc_clay + K + Ca + Mn, data=resid.data)
resid.model %>% summary()

resid.data.all <- data.frame(mgmt.model$residuals,ca.model$residuals,all_data$perc_clay,k.model$residuals,mn.model$residuals)
names(resid.data.all) <- c('model_residuals','Ca','perc_clay','K','Mn')

resid.model.ca <- lm(model_residuals ~ perc_clay + Ca +K +Mn, data=resid.data.all)
resid.model.ca %>% summary()


effect_plot(model=resid.model.ca, pred = Ca, colors=all_data$System, 
            shape=all_data$Replicate, 
            x.label = "Calcium", y.label = "First-stage model residuals")
effect_plot(model=resid.model.ca, pred = K, colors=all_data$System, 
            shape=all_data$Replicate, 
            x.label = "Potassium", y.label = "First-stage model residuals")
effect_plot(model=resid.model.ca, pred = perc_clay, colors=all_data$System, 
            shape=all_data$Replicate, 
            x.label = "Percent clay", y.label = "First-stage model residuals")
effect_plot(model=resid.model.ca, pred = Mn, colors=all_data$System, 
            shape=all_data$Replicate, 
            x.label = "Manganese", y.label = "First-stage model residuals")

stargazer(mgmt.model, resid.model.ca, out="maom.htm")
Fig5a