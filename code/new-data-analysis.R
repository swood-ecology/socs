# LOAD PACKAGES
library(tidyverse)    # For data manipulation
library(readxl)       # For reading in .xls files
library(lme4)         # For mixed-effects regression
library(lmerTest)     # For p-value calculations of lmer models
library(performance)  # For R2
library(simr)         # For power analysis

# library(doBy)
# library(grid)         # For manipulating details of ggplot
# library(sjstats)
# library(stargazer)    # For generating HTML tables for export
# library(jtools)       # For interaction plots


# FUNCTIONS
process_pom <- function(data){
  massPOM <- data$`tin-dry-sample-g` - data$`tin-g`
  massMAOM_NaHMP <- (data$`nalgene-dry-sample-g`-data$`nalgene-g`)*(data$`wash-ml`/data$`subsample-ml`)
  expected_NaHMP <- (data$`nahmp-ml`/1000) * data$`nahmp-g`
  massRecovered <- massPOM + (massMAOM_NaHMP - expected_NaHMP)
  POMC <- ((massPOM*(data$`c-pom`/100))*1000)/data$`sample-g`
  MAOMC <- (((massMAOM_NaHMP-expected_NaHMP)*(data$`c-maom`/100))*1000)/data$`sample-g`
  percPOMC <- POMC*.1 
  percMAOMC <- MAOMC * .1
  sumFraction <- percPOMC + percMAOMC
  POMN <- ((massPOM*(data$`n-pom`/100))*1000)/data$`sample-g`
  MAOMN <- (((massMAOM_NaHMP-expected_NaHMP)*(data$`n-maom`/100))*1000)/data$`sample-g`
  logPOMCN <- log(POMC/POMN)
  logMAOMCN <- log(MAOMC/MAOMN)
  df <- tibble(data$`lab-id`,POMC,MAOMC,POMN,MAOMN,percPOMC,percMAOMC,sumFraction,logPOMCN,logMAOMCN,data$`d13C-maom`,data$`d13C-pom`)
  names(df) <- c('lab-id','POMC','MAOMC','POMN','MAOMN','percPOMC','percMAOMC','sumFraction','logPOMCN','logMAOMCN','d13CMAOM','d13CPOM')
  return(df)
}


# READ IN DATA
pom <- read_csv("data/new-data/analysis/pom-data.csv")        # New POM data
metadata <- read_csv("data/new-data/metadata_socs-soilc.csv") # Site metadata
total.data <- read_csv("data/new-data/analysis/total-c.csv")  # New total C data
bd <- read_excel("data/new-data/analysis/bulk-density.xlsx")  # BD data
davis <- read_csv("data/davis-soils.csv")                     # Davis total C data


# MANIPULATE DATA #
metadata_select <- metadata %>% 
  select(-country,-region,-`sub-region`,-latitute,-longitude,-depth,-`date-sampled`,-`sample-id`)

## Fraction data ##
# pom_processed <- process_pom(pom)
pom.data <- full_join(pom, metadata_select)
names(pom.data)[14] <- 'trt'
pom.data <- left_join(pom.data, bd, by=c('mgmt-1','mgmt-2','mgmt-3','mgmt-4','mgmt-5','trt','year')) %>%
  mutate(POMC_Mg_ha = (`POMC (mg g-1)`*(`Bulk Density (Mg/m3)`*3000))/1000,
         MAOMC_Mg_ha = (`MAOMC (mg g-1)`*(`Bulk Density (Mg/m3)`*3000))/1000,
         POMN_Mg_ha = (`POMN (mg g-1)`*(`Bulk Density (Mg/m3)`*3000))/1000,
         MAOMN_Mg_ha = (`MAOMN (mg g-1)`*(`Bulk Density (Mg/m3)`*3000))/1000,
         POM_CN = POMC_Mg_ha/POMN_Mg_ha,
         MAOM_CN = MAOMC_Mg_ha/MAOMN_Mg_ha)
pom.data <- pom.data %>%
  mutate(unique = paste(`mgmt-1`,`mgmt-5`))

## Total data ##
total.data <- total.data %>% select(-`socs-id`)
names(total.data)[1] <- 'lab-id'
total.c <- full_join(total.data , metadata_select, by=c('lab-id','year')) 
names(total.c)[13] <- 'trt'
total.c <- left_join(total.c, bd, by=c('mgmt-1','mgmt-2','mgmt-3','mgmt-4','mgmt-5','trt','year')) %>%
  mutate(C_Mg_ha = (`C_kg_Mg`*(`Bulk Density (Mg/m3)`*3000))/1000,
         Total_CN = `%C`/`%N`)
total.c <- total.c %>%
  mutate(unique = paste(`mgmt-1`,`mgmt-5`))

## Break out by year
pom.y8 <- pom.data %>% filter(year==8)
pom.y0 <- pom.data %>% filter(year==0)
total.y8 <- total.c %>% filter(year==8)
total.y0 <- total.c %>% filter(year==0)
all_data <- full_join(total.c, pom.data)

# Table 2 & 3
# Replace variable for all variables in the model
aggregate(C_kg_Mg ~ `mgmt-2`+`mgmt-3`+`mgmt-4`,FUN=mean,data=total.y8)

## Calculate difference
diff.data <- all_data %>% group_by(unique) %>% arrange(year, .by_group = TRUE) %>%
  mutate(POMCdiff = last(`POMC (mg g-1)`) - first(`POMC (mg g-1)`),
         MAOMCdiff = last(`MAOMC (mg g-1)`) - first(`MAOMC (mg g-1)`),
         POMCNdiff = last(POM_CN) - first(POM_CN),
         MAOMCNdiff = last(MAOM_CN) - first(MAOM_CN),
         d13CPOMdiff = last(d13CPOM) - first(d13CPOM),
         d13CMAOMdiff = last(d13CMAOM) - first(d13CMAOM),
         Cdiff = last(`C_kg_Mg`) - first(`C_kg_Mg`),
         d13Cdiff = last(d13C) - first(d13C)) %>%
  select(`mgmt-1`:`mgmt-5`,POMCdiff:d13Cdiff) %>% unique()

# ANALYSES #
## Davis Total C % vs. sum of POM & MIN %
davis <- merge(davis, metadata) %>% merge(pom.y8) %>%
  mutate(sumFraction = `POMC (mg g-1)`*.1 + `MAOMC (mg g-1)`*.1)

ggplot(davis,aes(x=sumFraction,y=perC_Davis)) + geom_point(size = 2.5, alpha=0.7) + 
  geom_abline(intercept = 0, slope = 1, linetype='longdash', alpha = 0.5) + 
  geom_smooth(method = "lm", se = FALSE) + 
  xlim(0.45,1.1) + ylim(0.45,1.1) + ylab("Total carbon (%)\n") + xlab("\nSum of POM & MAOM C (%)") + 
  theme_bw() +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  )

lm(perC_Davis ~ sumFraction,data=davis) %>% summary()


## Response of change in variables to treatment over time
### Fractions
maomdiff <- lmer(MAOMCdiff ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=diff.data)
maomdiff %>% r2()
summary(maomdiff)

pomdiff <- lmer(POMCdiff ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=diff.data)
pomdiff %>% r2()
summary(pomdiff)

maomcndiff <- lmer(MAOMCNdiff ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=diff.data)
maomcndiff %>% r2()
summary(maomcndiff)

pomcndiff <- lmer(POMCNdiff ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=diff.data)
pomcndiff %>% r2()
summary(pomcndiff)

cdiff <- lmer(Cdiff ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=diff.data)
cdiff %>% r2()
summary(cdiff)

### Isotopes
maomd13cdiff <- lmer(d13CMAOMdiff ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=diff.data)
maomd13cdiff %>% r2()
summary(maomd13cdiff)

pomd13cdiff <- lmer(d13CPOMdiff ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=diff.data)
pomd13cdiff %>% r2()
summary(pomd13cdiff)

d13cdiff <- lmer(d13Cdiff ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=diff.data)
d13cdiff %>% r2()
summary(d13cdiff)


## Response of treatment in year 8
### Fractions
maom <- lmer(`MAOMC (mg g-1)` ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=pom.y8)
maom %>% r2()
summary(maom)

pom <- lmer(`POMC (mg g-1)` ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=pom.y8)
pom %>% r2()
summary(pom)

maomcn <- lmer(MAOM_CN ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=pom.y8)
maomcn %>% r2()
summary(maomcn)

pomcn <- lmer(POM_CN ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=pom.y8)
pomcn %>% r2()
summary(pomcn)

totalc <- lmer(C_kg_Mg ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=total.y8)
totalc %>% r2()
summary(totalc)

### Isotopes
maomd13c <- lm(d13CMAOM ~ `mgmt-3` + `mgmt-4`, data=pom.y8 )
maomd13c %>% r2()
summary(maomd13c)

pomd13c <- lm(d13CPOM ~ `mgmt-3` + `mgmt-4` , data=pom.y8 )
pomd13c %>% r2()
summary(pomd13c)

d13c <- lmer(d13C ~ `mgmt-3` + `mgmt-4` + (1|`mgmt-5`), data=total.y8 )
d13c %>% r2()
summary(d13c)


## FIGURES ##
fig2data <- all_data %>% select(year,d13C:`mgmt-1`,`POMC (mg g-1)`,`MAOMC (mg g-1)`,Total_CN,C_Mg_ha,unique,d13CPOM,d13CMAOM,POMC_Mg_ha:MAOM_CN)

fig2data$year <- as.factor(fig2data$year)
fig2data$year <- recode(fig2data$year, '0' = '2003', '8' = '2011')

diff.data.numbers <- diff.data %>% select(Cdiff,MAOMCdiff,POMCdiff,POMCNdiff,MAOMCNdiff,
                                          d13Cdiff,d13CPOMdiff,d13CMAOMdiff,`mgmt-1`)

Rmisc::summarySE(diff.data.numbers, measurevar="MAOMCdiff", groupvars=c('`mgmt-1`')) %>% select(-N,-sd,-se)

### FIGURE 3 ###
## d13C
fig3aData <- Rmisc::summarySE(fig2data, measurevar="d13C", groupvars=c('year','`mgmt-1`')) %>% select(-N,-sd,-se)

d13C_plot <- ggplot(fig3aData, aes(x=`mgmt-1`, y=d13C, fill=year))+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=d13C-ci, ymax=d13C+ci), width=.2, position=position_dodge(.9)) +
  geom_point(data=fig2data,aes(y=d13C,x=`mgmt-1`), position=position_dodge(0.9), alpha=0.4) +
  coord_cartesian(ylim=c(-28,-26)) + ylab(expression(Delta*"13C")) +
  scale_fill_manual(values=c("#999999", "#E69F00"), name = "Year") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13)
  )
d13C_plot

## POM d13C
fig3bData <- Rmisc::summarySE(fig2data, measurevar="d13CPOM", groupvars=c('year','`mgmt-1`')) %>% select(-N,-sd,-se)

d13CPOM_plot <- ggplot(fig3bData, aes(x=`mgmt-1`, y=d13CPOM, fill=year))+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=d13CPOM-ci, ymax=d13CPOM+ci), width=.2,position=position_dodge(.9)) +
  geom_point(data=fig2data,aes(y=d13CPOM,x=`mgmt-1`), position=position_dodge(0.9), alpha=0.5) +
  coord_cartesian(ylim=c(-30,-22)) + ylab(expression(Delta*"13C POM")) +
  scale_fill_manual(values=c("#999999", "#E69F00"), name = "Year") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13)
  )
d13CPOM_plot

## MAOM d13C
fig3cData <- Rmisc::summarySE(fig2data, measurevar="d13CMAOM", groupvars=c('year','`mgmt-1`')) %>% select(-N,-sd,-se)

d13CMAOM_plot <- ggplot(fig3cData, aes(x=`mgmt-1`, y=d13CMAOM, fill=year))+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=d13CMAOM-ci, ymax=d13CMAOM+ci), width=.2,position=position_dodge(.9)) +
  geom_point(data=fig2data,aes(y=d13CMAOM,x=`mgmt-1`), position=position_dodge(0.9), alpha=0.5) +
  coord_cartesian(ylim=c(-32,-14)) + ylab(expression(Delta*"13C MAOM")) +
  scale_fill_manual(values=c("#999999", "#E69F00"), name = "Year") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13)
  )
d13CMAOM_plot


### FIGURE 2 ###

## Total C
fig1aData <- Rmisc::summarySE(fig2data, measurevar="C_kg_Mg", groupvars=c('year','`mgmt-1`')) %>% select(-N,-sd,-se)

TotalC_plot <- ggplot(fig1aData, aes(x=`mgmt-1`, y=C_kg_Mg, fill=year))+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=C_kg_Mg-ci, ymax=C_kg_Mg+ci), width=.2,position=position_dodge(0.9)) +
  geom_point(data=fig2data,aes(y=C_kg_Mg,x=`mgmt-1`), position=position_dodge(0.9), alpha=0.5) +
  coord_cartesian(ylim=c(0,21)) + ylab("Total carbon (kg per Mg)\n") +
  scale_fill_manual(values=c("#999999", "#E69F00"), name = "Year") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13)
  )
TotalC_plot


## POM C
fig1bData <- Rmisc::summarySE(fig2data, measurevar="POMC (mg g-1)", groupvars=c('year','`mgmt-1`')) %>% select(-N,-sd,-se)

POMC_plot <- ggplot(fig1bData, aes(x=`mgmt-1`, y=`POMC (mg g-1)`, fill=year))+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=`POMC (mg g-1)`-ci, ymax=`POMC (mg g-1)`+ci), width=.2, position=position_dodge(.9)) +
  geom_point(data=fig2data,aes(y=`POMC (mg g-1)`,x=`mgmt-1`), position=position_dodge(0.9), alpha=0.5) +
  ylab("POM C (mg per g)\n") +
  scale_fill_manual(values=c("#999999", "#E69F00"), name = "Year") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13)
  )
POMC_plot


## MAOM C
fig1cData <- Rmisc::summarySE(fig2data, measurevar="MAOMC (mg g-1)", groupvars=c('year','`mgmt-1`')) %>% select(-N,-sd,-se)

MAOMC_plot <- ggplot(fig1cData, aes(x=`mgmt-1`, y=`MAOMC (mg g-1)`, fill=year))+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=`MAOMC (mg g-1)`-ci, ymax=`MAOMC (mg g-1)`+ci), width=.2,position=position_dodge(.9)) + 
  geom_point(data=fig2data,aes(y=`MAOMC (mg g-1)`,x=`mgmt-1`), position=position_dodge(0.9), alpha=0.5) +
  ylab("MAOM C (mg per g)\n") +
  scale_fill_manual(values=c("#999999", "#E69F00"), name = "Year") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13)
  )
MAOMC_plot

## POM CN
fig1dData <- Rmisc::summarySE(fig2data, measurevar="POM_CN", groupvars=c('year','`mgmt-1`')) %>% select(-N,-sd,-se)
POMCN_plot <- ggplot(fig1dData, aes(x=`mgmt-1`, y=POM_CN, fill=year))+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=POM_CN-ci, ymax=POM_CN+ci), width=.2,position=position_dodge(.9)) +
  geom_point(data=fig2data,aes(y=POM_CN,x=`mgmt-1`), position=position_dodge(0.9), alpha=0.5) +
  ylab("POM C:N\n") +
  scale_fill_manual(values=c("#999999", "#E69F00"), name = "Year") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13)
  )
POMCN_plot

## MAOM CN
fig1eData <- Rmisc::summarySE(fig2data, measurevar="MAOM_CN", groupvars=c('year','`mgmt-1`')) %>% select(-N,-sd,-se)
MAOMCN_plot <- ggplot(fig1eData, aes(x=`mgmt-1`, y=MAOM_CN, fill=year))+
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=MAOM_CN-ci, ymax=MAOM_CN+ci), width=.2,position=position_dodge(.9)) +
  geom_point(data=fig2data,aes(y=MAOM_CN,x=`mgmt-1`), position=position_dodge(0.9), alpha=0.5) +
  ylab("MAOM C:N\n") +
  scale_fill_manual(values=c("#999999", "#E69F00"), name = "Year") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13)
  )
MAOMCN_plot


#####POWER ANALYSIS#####
simrOptions(progress=FALSE)

diff.data %>%
  mutate(Rep = case_when(`mgmt-5`=='Rep 1' ~ 1,
                         `mgmt-5`=='Rep 2' ~ 2,
                         `mgmt-5`=='Rep 3' ~ 3,
                         `mgmt-5`=='Rep 4' ~ 4
  )) %>%
  fastDummies::dummy_cols() -> power.data

lmer(POMCdiff ~ `mgmt-3` + `mgmt-4` + (1|Rep),
     data = power.data
) -> model1

compostPOM <- powerSim(model1, nsim=100, test = fcompare(y~`mgmt-3`))
compostPOM
ccPOM <- powerSim(model1, nsim=100, test = fcompare(y~`mgmt-4`))
ccPOM

mod = lm(POMCdiff ~ `mgmt-3` + `mgmt-4`, data=diff.data)  
car::Anova(mod,type="III")
etasq <- 2.1424/(2.1424+4.8964)
f2 <- etasq / (1-etasq)
pwr::pwr.f2.test(u=1,v=17,f2=f2,sig.level=0.05)

etasq <- 0.6888/(0.6888+4.8964)
f2 <- etasq / (1-etasq)
pwr::pwr.f2.test(u=1,v=17,f2=f2,sig.level=0.05)

mod = lm(MAOMCdiff ~ `mgmt-3` + `mgmt-4`, data=diff.data)  
car::Anova(mod,type="III")
etasq <- 1.3041/(1.3041+7.5853)
f2 <- etasq / (1-etasq)
pwr::pwr.f2.test(u=1,v=17,f2=f2,sig.level=0.05)

etasq <- 3.5589/(3.5589+7.5853)
f2 <- etasq / (1-etasq)
pwr::pwr.f2.test(u=1,v=17,f2=f2,sig.level=0.05)

mod = lm(Cdiff ~ `mgmt-3` + `mgmt-4`, data=diff.data)  
car::Anova(mod,type="III")
etasq <- 22.4450/(22.4450+29.0175)
f2 <- etasq / (1-etasq)
pwr::pwr.f2.test(u=1,v=17,f2=f2,sig.level=0.05)

etasq <- 5.8800/(5.8800+29.0175)
f2 <- etasq / (1-etasq)
pwr::pwr.f2.test(u=1,v=17,f2=f2,sig.level=0.05)
