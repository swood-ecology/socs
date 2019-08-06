Salinas Organic Cropping Systems Experiment
-------------------------------------------

### Read in the data files

The data are stored in separate files. usda\_soil\_data has the soil
carbon fractionation data analyzed by Wood, icp has the data analyzed by
UC Davis, and the other site-level descriptors are in the other files.

    # READ IN ALL DATA FILES
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

### Data manipulation

Now we need to do some slight manipulations to rename variables, merge
the data sets together, and select the variables that we'll want to used
for analysis.

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

    inputs$RepTreat <- paste(inputs$trt, inputs$rep)
    usda_soil_data$RepTreat <- paste(usda_soil_data$Trt_id, usda_soil_data$Replicate)
    bd$RepTreat <- paste(bd$trt, bd$rep)
    icp$RepTreat <- paste(icp$Trt_id, icp$Replicate)

    all_data <- merge(usda_soil_data, bd, by = "RepTreat") %>%
      select(RepTreat, System, Replicate, Trt_id, cover_crop, Compost:POM_CN, blkden) %>%
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
      merge(icp %>% select(IC,EC:RepTreat), by="RepTreat")

    rm(bd); rm(inputs); rm(usda_soil_data); rm(icp)

### Models of management effect on soil fractions

Now that we have the data, let's look through some models. We are going
to used hierarchical models with a random intercept for replicate block.
This is to capture random differences among plots within a block.

Let's start with a series of models, using all of the data, and looking
at the impact on particulate and mineral-associated C and microbial
biomass.

    # Particulate carbon
    lmer(pom.stock ~ Compost + Cover_crop_freq + cover_crop + (1|Replicate),
         data = all_data
    ) -> pom.model

    ## boundary (singular) fit: see ?isSingular

    summary(pom.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## pom.stock ~ Compost + Cover_crop_freq + cover_crop + (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 57.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.77387 -0.73387  0.06072  0.45096  1.96465 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.000    0.000   
    ##  Residual              1.694    1.302   
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error      df t value
    ## (Intercept)                       4.1441     1.1271 15.0000   3.677
    ## CompostYes                        1.8568     0.9203 15.0000   2.018
    ## Cover_crop_freqEvery 4th Winter  -2.0884     0.9203 15.0000  -2.269
    ## cover_cropmustard                -0.9863     0.9203 15.0000  -1.072
    ## cover_croprye                    -0.4376     0.9203 15.0000  -0.475
    ##                                 Pr(>|t|)   
    ## (Intercept)                      0.00224 **
    ## CompostYes                       0.06189 . 
    ## Cover_crop_freqEvery 4th Winter  0.03844 * 
    ## cover_cropmustard                0.30078   
    ## cover_croprye                    0.64129   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) CmpstY C__E4W cvr_crpm
    ## CompostYes  -0.816                       
    ## Cvr_crp_E4W -0.816  0.500                
    ## cvr_crpmstr -0.408  0.000  0.500         
    ## cover_crpry -0.408  0.000  0.500  0.500  
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

    r.squaredGLMM(pom.model)

    ## Warning: 'r.squaredGLMM' now calculates a revised statistic. See the help
    ## page.

    ##            R2m       R2c
    ## [1,] 0.5533291 0.5533291

    rm(pom.model)

    # Mineral-associated carbon
    lmer(maom.stock ~ Compost + Cover_crop_freq + cover_crop + (1|Replicate),
         data = all_data
    ) -> maom.model
    summary(maom.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## maom.stock ~ Compost + Cover_crop_freq + cover_crop + (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 134.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.68097 -0.35627 -0.04135  0.39491  2.01334 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 1976.0   44.45   
    ##  Residual               115.8   10.76   
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error     df t value
    ## (Intercept)                       46.275     24.101  4.025   1.920
    ## CompostYes                        10.970      7.609 12.000   1.442
    ## Cover_crop_freqEvery 4th Winter    3.399      7.609 12.000   0.447
    ## cover_cropmustard                 -4.882      7.609 12.000  -0.642
    ## cover_croprye                      1.945      7.609 12.000   0.256
    ##                                 Pr(>|t|)
    ## (Intercept)                        0.127
    ## CompostYes                         0.175
    ## Cover_crop_freqEvery 4th Winter    0.663
    ## cover_cropmustard                  0.533
    ## cover_croprye                      0.803
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) CmpstY C__E4W cvr_crpm
    ## CompostYes  -0.316                       
    ## Cvr_crp_E4W -0.316  0.500                
    ## cvr_crpmstr -0.158  0.000  0.500         
    ## cover_crpry -0.158  0.000  0.500  0.500

    r.squaredGLMM(maom.model)

    ##              R2m       R2c
    ## [1,] 0.008618742 0.9451247

    rm(maom.model)

    # Microbial biomass
    lmer(substrate_induced_respiration ~ Compost + Cover_crop_freq + cover_crop + (1|Replicate),
         data = all_data
    ) -> sir.model
    summary(sir.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## substrate_induced_respiration ~ Compost + Cover_crop_freq + cover_crop +  
    ##     (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 84.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.4699 -0.4505 -0.1223  0.5850  1.4465 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.2764   0.5258  
    ##  Residual              9.9341   3.1518  
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error      df t value
    ## (Intercept)                      23.1771     2.7422 13.6702   8.452
    ## CompostYes                        0.1661     2.2287 12.0000   0.075
    ## Cover_crop_freqEvery 4th Winter  -2.0835     2.2287 12.0000  -0.935
    ## cover_cropmustard                -3.0430     2.2287 12.0000  -1.365
    ## cover_croprye                    -2.3570     2.2287 12.0000  -1.058
    ##                                 Pr(>|t|)    
    ## (Intercept)                     8.52e-07 ***
    ## CompostYes                         0.942    
    ## Cover_crop_freqEvery 4th Winter    0.368    
    ## cover_cropmustard                  0.197    
    ## cover_croprye                      0.311    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) CmpstY C__E4W cvr_crpm
    ## CompostYes  -0.813                       
    ## Cvr_crp_E4W -0.813  0.500                
    ## cvr_crpmstr -0.406  0.000  0.500         
    ## cover_crpry -0.406  0.000  0.500  0.500

    r.squaredGLMM(sir.model)

    ##             R2m       R2c
    ## [1,] 0.09801545 0.1224342

    rm(sir.model)

Looking at the model estimates, there does not seem to be a relationship
between cover crop type and our variables of interest. We could continue
by dropping that variable from the model and re-running, but that would
take the average effect of the other variables across systems with
different cover crops. Even though cover crop type does not have an
effect in our statistical model, for the remaining comparisons we felt
that it makes more sense to only look at the variables that we know have
the same cover crop type.

Let's filter out observations that are not legume rye, and re-run the
models. We'll also select a smaller set of predictor variables that
might be of special interest.

    leg_rye <- all_data %>%
    filter(cover_crop=="legume-rye") %>%
    select(Replicate, cover_crop, Compost, Cover_crop_freq, organic_matter, 
           water_holding_capacity,substrate_induced_respiration:POM_CN,pom.stock,
           pom.n.stock, maom.stock, maom.n.stock, 
           total_C_no_roots_no_exudates, CEC, pH, `X-Ca`
    )

    # Organic matter
    lmer(organic_matter ~ Compost + Cover_crop_freq + (1|Replicate),
    data = leg_rye
    ) -> om.model
    summary(om.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: organic_matter ~ Compost + Cover_crop_freq + (1 | Replicate)
    ##    Data: leg_rye
    ## 
    ## REML criterion at convergence: 2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.12941 -0.65963 -0.08913  0.34139  1.80512 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.002765 0.05258 
    ##  Residual              0.043224 0.20790 
    ## Number of obs: 12, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error      df t value
    ## (Intercept)                       2.7110     0.1820  7.5843  14.899
    ## CompostYes                        0.6683     0.1470  6.0000   4.546
    ## Cover_crop_freqEvery 4th Winter  -0.4336     0.1470  6.0000  -2.949
    ##                                 Pr(>|t|)    
    ## (Intercept)                     6.88e-07 ***
    ## CompostYes                       0.00391 ** 
    ## Cover_crop_freqEvery 4th Winter  0.02563 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) CmpstY
    ## CompostYes  -0.808       
    ## Cvr_crp_E4W -0.808  0.500

    r.squaredGLMM(om.model)

    ##           R2m     R2c
    ## [1,] 0.829734 0.83997

    rm(om.model)

    # Particulate C
    lmer(pom.stock ~ Compost + Cover_crop_freq + (1|Replicate),
       data = leg_rye
    ) -> pom.model

    ## boundary (singular) fit: see ?isSingular

    summary(pom.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pom.stock ~ Compost + Cover_crop_freq + (1 | Replicate)
    ##    Data: leg_rye
    ## 
    ## REML criterion at convergence: 36.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5465 -0.6398  0.1154  0.3699  1.7128 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.000    0.000   
    ##  Residual              2.229    1.493   
    ## Number of obs: 12, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error     df t value
    ## (Intercept)                        4.144      1.293  9.000   3.205
    ## CompostYes                         1.857      1.056  9.000   1.759
    ## Cover_crop_freqEvery 4th Winter   -2.088      1.056  9.000  -1.978
    ##                                 Pr(>|t|)  
    ## (Intercept)                       0.0107 *
    ## CompostYes                        0.1124  
    ## Cover_crop_freqEvery 4th Winter   0.0793 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) CmpstY
    ## CompostYes  -0.816       
    ## Cvr_crp_E4W -0.816  0.500
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

    r.squaredGLMM(pom.model)

    ##            R2m       R2c
    ## [1,] 0.5597249 0.5597249

    rm(pom.model)

    # Particulate N
    lmer(pom.n.stock ~ Compost + Cover_crop_freq + (1|Replicate),
       data = leg_rye
    ) -> pom.n.model
    summary(pom.n.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pom.n.stock ~ Compost + Cover_crop_freq + (1 | Replicate)
    ##    Data: leg_rye
    ## 
    ## REML criterion at convergence: -9.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.1617 -0.4674 -0.1895  0.8013  0.9629 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.016439 0.12821 
    ##  Residual              0.006066 0.07788 
    ## Number of obs: 12, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error       df t value
    ## (Intercept)                      0.39427    0.09305  7.63242   4.237
    ## CompostYes                       0.08966    0.05507  5.99956   1.628
    ## Cover_crop_freqEvery 4th Winter -0.15178    0.05507  5.99956  -2.756
    ##                                 Pr(>|t|)   
    ## (Intercept)                      0.00317 **
    ## CompostYes                       0.15465   
    ## Cover_crop_freqEvery 4th Winter  0.03303 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) CmpstY
    ## CompostYes  -0.592       
    ## Cvr_crp_E4W -0.592  0.500

    r.squaredGLMM(pom.n.model)

    ##            R2m       R2c
    ## [1,] 0.3249311 0.8180477

    rm(pom.n.model)

    # Mineral-associated C
    lmer(maom.stock ~ Compost + Cover_crop_freq + (1|Replicate),
       data = leg_rye
    ) -> maom.model
    summary(maom.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: maom.stock ~ Compost + Cover_crop_freq + (1 | Replicate)
    ##    Data: leg_rye
    ## 
    ## REML criterion at convergence: 84.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.1281 -0.2813 -0.0933  0.2000  1.4875 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 1998.5   44.70   
    ##  Residual               113.1   10.64   
    ## Number of obs: 12, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error     df t value
    ## (Intercept)                       46.275     24.176  3.912   1.914
    ## CompostYes                        10.970      7.522  6.000   1.458
    ## Cover_crop_freqEvery 4th Winter    3.399      7.522  6.000   0.452
    ##                                 Pr(>|t|)
    ## (Intercept)                        0.130
    ## CompostYes                         0.195
    ## Cover_crop_freqEvery 4th Winter    0.667
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) CmpstY
    ## CompostYes  -0.311       
    ## Cvr_crp_E4W -0.311  0.500

    r.squaredGLMM(maom.model)

    ##           R2m       R2c
    ## [1,] 0.010745 0.9469916

    rm(maom.model)

    # Mineral-associated N
    lmer(maom.n.stock ~ Compost + Cover_crop_freq + (1|Replicate),
       data = leg_rye
    ) -> maom.n.model
    summary(maom.n.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: maom.n.stock ~ Compost + Cover_crop_freq + (1 | Replicate)
    ##    Data: leg_rye
    ## 
    ## REML criterion at convergence: 30.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.36769 -0.34958 -0.03576  0.44752  1.41898 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 9.3563   3.0588  
    ##  Residual              0.2248   0.4741  
    ## Number of obs: 12, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error     df t value
    ## (Intercept)                       4.0347     1.5836 3.3866   2.548
    ## CompostYes                        0.8487     0.3353 6.0000   2.531
    ## Cover_crop_freqEvery 4th Winter   0.5198     0.3353 6.0000   1.550
    ##                                 Pr(>|t|)  
    ## (Intercept)                       0.0746 .
    ## CompostYes                        0.0446 *
    ## Cover_crop_freqEvery 4th Winter   0.1720  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) CmpstY
    ## CompostYes  -0.212       
    ## Cvr_crp_E4W -0.212  0.500

    r.squaredGLMM(maom.n.model)

    ##             R2m       R2c
    ## [1,] 0.01370787 0.9768593

    rm(maom.n.model)

For our model of microbial biomass--measured by substrate induced
respiration--we're going to take a slightly different approach. Whereas
for the previous models we were mainly interested in how soil fractions
respond to management, now we think that microbial biomass could depend
on the stock sizes themselves, which we might want to control for. So we
can look at including stocks as predictor variables. For this we're
going to actually keep all of the data, rather than drop cover crop
type, because visual inspection suggested that microbial biomass might
be slightly higher with leguminous cover crops.

    # SIR
    lmer(substrate_induced_respiration ~ cover_crop + Compost + Cover_crop_freq + (1|Replicate),
         data = all_data
    ) -> sir.model
    summary(sir.model)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## substrate_induced_respiration ~ cover_crop + Compost + Cover_crop_freq +  
    ##     (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 84.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.4699 -0.4505 -0.1223  0.5850  1.4465 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.2764   0.5258  
    ##  Residual              9.9341   3.1518  
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error      df t value
    ## (Intercept)                      23.1771     2.7422 13.6702   8.452
    ## cover_cropmustard                -3.0430     2.2287 12.0000  -1.365
    ## cover_croprye                    -2.3570     2.2287 12.0000  -1.058
    ## CompostYes                        0.1661     2.2287 12.0000   0.075
    ## Cover_crop_freqEvery 4th Winter  -2.0835     2.2287 12.0000  -0.935
    ##                                 Pr(>|t|)    
    ## (Intercept)                     8.52e-07 ***
    ## cover_cropmustard                  0.197    
    ## cover_croprye                      0.311    
    ## CompostYes                         0.942    
    ## Cover_crop_freqEvery 4th Winter    0.368    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) cvr_crpm cvr_crpr CmpstY
    ## cvr_crpmstr -0.406                         
    ## cover_crpry -0.406  0.500                  
    ## CompostYes  -0.813  0.000    0.000         
    ## Cvr_crp_E4W -0.813  0.500    0.500    0.500

    r.squaredGLMM(sir.model)

    ##             R2m       R2c
    ## [1,] 0.09801545 0.1224342

    rm(sir.model)

    # SIR, depending on N input to microbes
    lmer(substrate_induced_respiration ~ cover_crop + Compost + Cover_crop_freq + pom.n.stock + (1|Replicate),
         data = all_data
    ) -> sir.model.cn

    ## boundary (singular) fit: see ?isSingular

    car::vif(sir.model.cn)

    ##                     GVIF Df GVIF^(1/(2*Df))
    ## cover_crop      1.828966  2        1.162924
    ## Compost         1.663863  1        1.289908
    ## Cover_crop_freq 2.674523  1        1.635397
    ## pom.n.stock     1.393720  1        1.180559

    summary(sir.model.cn)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## substrate_induced_respiration ~ cover_crop + Compost + Cover_crop_freq +  
    ##     pom.n.stock + (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 74.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.75939 -0.53354 -0.09896  0.56826  1.53206 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.000    0.000   
    ##  Residual              7.731    2.781   
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error      df t value
    ## (Intercept)                      19.0136     2.9635 14.0000   6.416
    ## cover_cropmustard                -2.4419     1.9819 14.0000  -1.232
    ## cover_croprye                    -2.0627     1.9699 14.0000  -1.047
    ## CompostYes                       -0.7807     2.0050 14.0000  -0.389
    ## Cover_crop_freqEvery 4th Winter  -0.4807     2.0755 14.0000  -0.232
    ## pom.n.stock                      10.5602     4.3812 14.0000   2.410
    ##                                 Pr(>|t|)    
    ## (Intercept)                     1.61e-05 ***
    ## cover_cropmustard                 0.2382    
    ## cover_croprye                     0.3128    
    ## CompostYes                        0.7028    
    ## Cover_crop_freqEvery 4th Winter   0.8202    
    ## pom.n.stock                       0.0303 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) cvr_crpm cvr_crpr CmpstY C__E4W
    ## cvr_crpmstr -0.402                                
    ## cover_crpry -0.367  0.503                         
    ## CompostYes  -0.536 -0.025   -0.012                
    ## Cvr_crp_E4W -0.815  0.510    0.493    0.402       
    ## pom.n.stock -0.583  0.126    0.062   -0.196  0.320
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

    r.squaredGLMM(sir.model.cn)

    ##            R2m       R2c
    ## [1,] 0.3100065 0.3100065

    rm(sir.model.cn)

### Predicting soil fractions by continuous inputs

There is a certain feedback between some of the treatments and the
drivers of soil fraction. For instance, cover crops increase soil carbon
pools partially because of root exudation. More cover crop biomass can
lead to more root inputs. On plots with compost, cover crop biomass can
be greater, which can mean different plots have different amounts of
carbon going into the soil through cover crops. In other words, cover
crop frequency is not a perfect categorical variable to capture the
amount of carbon going belowground.

We estimated the quantitative amount of nitrogen and carbon inputs going
into the system. Let's fit some models to explore if any of these
variables predict soil fractions.

First, we're going to use a machine learning algorithm called Random
Forest to see which of the quantitative input variables are the most
closely related to the fractions that we're interested in. We're going
to create a list of all of the different Random Forest models for all of
our outcomes of interest

    rf.res <- list(
      VSURF(
        y = all_data$organic_matter,
        x = all_data %>%
          select(total_compost:total_C_no_roots_no_exudates),
        parallel = TRUE
      ),
      VSURF(
        y = all_data$pom.stock,
        x = all_data %>%
          select(total_compost:total_C_no_roots_no_exudates),
        parallel = TRUE
      ),
      VSURF(
        y = all_data$maom.stock,
        x = all_data %>%
          select(total_compost:total_C_no_roots_no_exudates),
        parallel = TRUE
      ),
      VSURF(
        y = all_data$substrate_induced_respiration,
        x = all_data %>%
          select(total_compost:total_C_no_roots_no_exudates),
        parallel = TRUE
      )
    )

Now we can look at each of these model results and see which variables
are the most relevant. This Random Forest procedure generates lists of
variables most suitable for both interpretatoin and prediction. Briefly,
prediction is more strict in the number of variables. We're not using
our data for predictive modeling, so we'll look at the variables
identified for interpretation.

    for(i in 1:length(rf.res)){
      vars = rf.res[[i]]$varselect.interp
      # Print the response variable
      print(rf.res[[i]]$call$y)
      # Print the predictor variables selected
      all_data %>% 
        select(total_compost:total_C_no_roots_no_exudates) %>%
        select(vars) %>% names() %>% print()
    }

    ## all_data$organic_matter
    ## [1] "total_C_no_roots_no_exudates" "total_om"                    
    ## [3] "fresh_om"                     "total_cc_shoot"              
    ## [5] "fresh_om_perc"                "total_compost"               
    ## [7] "compost_n"                   
    ## all_data$pom.stock
    ## [1] "fresh_om"
    ## all_data$maom.stock
    ## [1] "total_veg_residue_shoot" "annual_legume_cc_shoot" 
    ## all_data$substrate_induced_respiration
    ## [1] "annual_cc_shoot"

Some of these variables are highly correlated with each other so we'll
select out the ones that we know to be somewhat independent and fit
models with these.

    # Organic matter
    lm.om <- lmer(organic_matter ~ total_C_no_roots_no_exudates + (1 | Replicate), data = all_data)

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    ## boundary (singular) fit: see ?isSingular

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    summary(lm.om)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: organic_matter ~ total_C_no_roots_no_exudates + (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 28.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.79229 -0.73045  0.04049  0.63798  1.88992 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.00000  0.0000  
    ##  Residual              0.06617  0.2572  
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                               Estimate Std. Error        df t value
    ## (Intercept)                  1.773e+00  1.931e-01 1.800e+01   9.181
    ## total_C_no_roots_no_exudates 1.738e-05  2.608e-06 1.800e+01   6.663
    ##                              Pr(>|t|)    
    ## (Intercept)                  3.27e-08 ***
    ## total_C_no_roots_no_exudates 2.98e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## ttl_C_n_r__ -0.955
    ## fit warnings:
    ## Some predictor variables are on very different scales: consider rescaling
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

    r.squaredGLMM(lm.om)

    ##            R2m       R2c
    ## [1,] 0.7002862 0.7002862

    rm(lm.om)

    # Particulate
    lm.pc <- lmer(pom.stock ~ fresh_om + (1 | Replicate), data = all_data)

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    ## boundary (singular) fit: see ?isSingular

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    summary(lm.pc)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: pom.stock ~ fresh_om + (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 87.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6891 -0.6511 -0.1953  0.6046  2.1321 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.000    0.000   
    ##  Residual              1.749    1.322   
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept) -4.528e-01  1.167e+00  1.800e+01  -0.388 0.702507    
    ## fresh_om     5.561e-05  1.265e-05  1.800e+01   4.396 0.000348 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr)
    ## fresh_om -0.967
    ## fit warnings:
    ## Some predictor variables are on very different scales: consider rescaling
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

    r.squaredGLMM(lm.pc)

    ##            R2m       R2c
    ## [1,] 0.5042523 0.5042523

    rm(lm.pc)

    # Mineral
    lm.mc <- lmer(maom.stock ~ total_veg_residue_shoot + annual_legume_cc_shoot + (1 | Replicate), data = all_data)

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    summary(lm.mc)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: maom.stock ~ total_veg_residue_shoot + annual_legume_cc_shoot +  
    ##     (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 182.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.60560 -0.27920 -0.01373  0.29544  2.01977 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 1990.2   44.61   
    ##  Residual               123.5   11.11   
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error        df t value Pr(>|t|)
    ## (Intercept)             5.023e+01  4.531e+01 1.644e+01   1.108    0.284
    ## total_veg_residue_shoot 9.887e-05  7.364e-04 1.405e+01   0.134    0.895
    ## annual_legume_cc_shoot  3.849e-04  2.397e-03 1.401e+01   0.161    0.875
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) ttl___
    ## ttl_vg_rsd_ -0.866       
    ## annl_lgm_c_ -0.467  0.472
    ## fit warnings:
    ## Some predictor variables are on very different scales: consider rescaling

    r.squaredGLMM(lm.mc)

    ##               R2m       R2c
    ## [1,] 0.0001041067 0.9415994

    rm(lm.mc)

    # SIR
    lm.sir <- lmer(substrate_induced_respiration ~ annual_cc_shoot + (1 | Replicate), data = all_data)
    summary(lm.sir)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: substrate_induced_respiration ~ annual_cc_shoot + (1 | Replicate)
    ##    Data: all_data
    ## 
    ## REML criterion at convergence: 111.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8243 -0.6817  0.1658  0.6018  1.3704 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  Replicate (Intercept) 0.2366   0.4865  
    ##  Residual              9.2802   3.0463  
    ## Number of obs: 20, groups:  Replicate, 4
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error        df t value Pr(>|t|)   
    ## (Intercept)     1.810e+01  5.520e+00 1.766e+01   3.280  0.00424 **
    ## annual_cc_shoot 4.731e-04  7.866e-04 1.743e+01   0.601  0.55527   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## annl_cc_sht -0.991

    r.squaredGLMM(lm.sir)

    ##             R2m        R2c
    ## [1,] 0.01844617 0.04285375

    rm(lm.sir)

A surprising part of these models is how poorly the fixed effects
predict mineral-associated C, while how well the random effects do. This
says that there's an important component of among-replicate variation
that we are not capturing with fixed effects. Since mineral-associated C
is known to be maintained by physical and chemical properties of the
soil, there might be some soil properties that explain this variation.
We would not expect soil properties to be as important to particulate C,
because that C is mainly in the form of partially broken down plant
material and it is not as reactive with soil surfaces.

To explore the potential for soil properties to impact
mineral-associated carbon we used a tool called two-stage least squares
regression. This approach, which is more commonly used by economists,
first models the effect of the management practices on
mineral-associated carbon (we've already done this). Then we want to
know how much of the left-over variation we can explain with other soil
properties. To do that, we extract the residuals from the first
regression and regress those against soil properties we think might be
important.

Let's start with the management model, even though we've already done
this above. First we're going to use Random Forest to identify the most
important management properties for mineral-associated C. Then we'll
used those predictors in a linear model. We're not using random effects
because we want to try to explain the variance that the random effects
were capturing with new fixed effects.

    # Determine best management predictors of MAOM
    mgmt <- VSURF(
      y = all_data$maom.stock,
      x = all_data %>%
        select(total_compost:total_C_no_roots_no_exudates),
      parallel = TRUE
    )

    # Fit model with best management predictors
    mgmt.model <- lm(maom.stock ~ Compost + Cover_crop_freq + total_veg_residue_shoot, data = all_data)
    mgmt.model %>% summary()

    ## 
    ## Call:
    ## lm(formula = maom.stock ~ Compost + Cover_crop_freq + total_veg_residue_shoot, 
    ##     data = all_data)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -49.17 -25.20 -14.44  26.13  60.60 
    ## 
    ## Coefficients:
    ##                                   Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                     442.545423 178.726290   2.476   0.0248 *
    ## CompostYes                       34.967715  29.326072   1.192   0.2505  
    ## Cover_crop_freqEvery 4th Winter -32.110372  27.566440  -1.165   0.2612  
    ## total_veg_residue_shoot          -0.007857   0.003486  -2.254   0.0386 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 38.64 on 16 degrees of freedom
    ## Multiple R-squared:  0.2468, Adjusted R-squared:  0.1056 
    ## F-statistic: 1.748 on 3 and 16 DF,  p-value: 0.1976

Now that we have our management model, we can use a Random Forest
approach to identify the variables related to the residuals of the
previous model. We can then use those in a model. We'll then plot the
effect of each of the final selected variables on the residuals. This
will show which variables are the most important to explaining the
variation in mineral-associated C that is not explained by management.
We will also look at model fit statistics (like R2) to see how well this
model does.

    # Find best non-management predictors of model residuals
    soil.prop <- VSURF(
      y = mgmt.model$residuals,
      x = all_data %>%
        select(perc_sand:perc_clay, EC:`Fe (DTPA)`),
      parallel = TRUE
    )

    # Create data frame of resulting variables, also including clay
    resid.data <- data.frame(mgmt.model$residuals, all_data$perc_clay, all_data$`X-K`, all_data$`Ca (SP)`, all_data$`Mn (DTPA)`)
    names(resid.data) <- c("model_residuals", "perc_clay", "K", "Ca", "Mn")

    resid.model <- lm(model_residuals ~ perc_clay + K + Ca + Mn, data = resid.data)
    resid.model %>% summary()

    ## 
    ## Call:
    ## lm(formula = model_residuals ~ perc_clay + K + Ca + Mn, data = resid.data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -44.518 -14.003   7.291  18.947  29.011 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  -1.2515    85.9788  -0.015  0.98858   
    ## perc_clay     0.2174     3.1264   0.070  0.94548   
    ## K            -0.3872     0.1136  -3.409  0.00389 **
    ## Ca           20.4718     5.9750   3.426  0.00375 **
    ## Mn           -1.1775     1.0419  -1.130  0.27615   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 26.95 on 15 degrees of freedom
    ## Multiple R-squared:  0.5442, Adjusted R-squared:  0.4226 
    ## F-statistic: 4.477 on 4 and 15 DF,  p-value: 0.01403

    effect_plot(model = resid.model, pred = Ca, interval = TRUE, rug = TRUE, plot.points = TRUE, partial.residuals = TRUE)

![](data-analysis_files/figure-markdown_strict/unnamed-chunk-10-1.png)

    effect_plot(model = resid.model, pred = K, interval = TRUE, rug = TRUE, plot.points = TRUE, partial.residuals = TRUE)

![](data-analysis_files/figure-markdown_strict/unnamed-chunk-10-2.png)

    effect_plot(model = resid.model, pred = perc_clay, interval = TRUE, rug = TRUE, plot.points = TRUE, partial.residuals = TRUE)

![](data-analysis_files/figure-markdown_strict/unnamed-chunk-10-3.png)

    effect_plot(model = resid.model, pred = Mn, interval = TRUE, rug = TRUE, plot.points = TRUE, partial.residuals = TRUE)

![](data-analysis_files/figure-markdown_strict/unnamed-chunk-10-4.png)

The adjusted R2 of this model is over 0.4! That's a huge improvement
over the management model, which was barely predicting 10% of the
variation.
