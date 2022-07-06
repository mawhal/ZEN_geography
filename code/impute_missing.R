###################################################################################
# IMPUTE MISSING DATA                                                             #
###################################################################################

# Note: this code runs with data assembly scripts (e.g., data_assembly.R) to umpute 
# missing data values using random forest models



# First the rationale: Model comparisons requires that alternative models use exactly the
# same source data set. Missing cells result in slightly different data sets for models
# that include, versus do not include, that variable, with different DF such that
# resulting AIC scores used to compare them will be invalid. Solving this requires
# either (1) throwing out all rows that have a missing cell, or (2) imputing (modeling)
# the missing values. The former is far the worse alternative since it will end up
# discarding a substantial part of the entire data set.

# A few other points: we cannot include the following predictors because biased by
# missing from entire site: periphyton (missing from SW.A). Also BC.A is missing 12
# of 20 samples for Zostera biomass and shoot density. I am going ahead to impute the
# missing Zostera data for BC.A but we may want to exclude this site from any analysis
# that uses Zostera mass or shoot density as either a predictor or  response.

# Note that a data set used for RF imputation model cannot include any rows that lack
# data for one of the predictor variables.

# We first impute missing values for Zostera variables (AG biomass, shoot density, %N), then
# proceed to the grazer variables that depend on these grass variables

# # Remind me which variables are missing data:
# apply(d,2,pMiss)

# create a temporary dataframe used to model (impute) the missing values
y <- d

# LEAF % NITROGEN
# # Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.Leaf.PercN)) # 14

leafN.rf = randomForest(log10.Leaf.PercN ~ Ocean + Coast + Latitude + Longitude
                        + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                        + log10.chlomean
                        + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                        + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                        + log10.periphyton.mass.per.g.zostera
                        + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                        ## + log10.Leaf.PercN
                        # + pop.density.2015 + AllelicRichness
                        ,
                        na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Leaf.PercN),
  "log10.Leaf.PercN"] = predict(leafN.rf, y[is.na(y$log10.Leaf.PercN), ])
# sum(is.na(y$log10.Leaf.PercN)) # 0


# ZOSTERA SHOOT DENSITY
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.Zostera.shoots.core)) # 15
shootdensity.rf = randomForest(log10.Zostera.shoots.core ~ Ocean + Coast + Latitude + Longitude
                               + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                               + log10.chlomean
                               # + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                               + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                               # + log10.periphyton.mass.per.g.zostera
                               # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                               + log10.Leaf.PercN
                               # + pop.density.2015
                               + AllelicRichness
                               ,
                               na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Zostera.shoots.core),
  "log10.Zostera.shoots.core"] = predict(shootdensity.rf, y[is.na(y$log10.Zostera.shoots.core), ])
# sum(is.na(y$log10.Zostera.shoots.core)) # 0


# ZOSTERA ABOVE-GROUND BIOMASS
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.Zostera.AG.mass)) # 24
ZAG.rf = randomForest(log10.Zostera.AG.mass ~ Ocean + Coast + Latitude + Longitude
                      # + sst.mean + Salinity.ppt + parmean + log10.nitrate + log10.phosphate
                      + log10.chlomean
                      # + log10.Zostera.AG.mass
                      + log10.Zostera.shoots.core
                      + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                      # + log10.periphyton.mass.per.g.zostera
                      # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                      + log10.Leaf.PercN
                      # + pop.density.2015
                      + AllelicRichness
                      ,
                      na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Zostera.AG.mass),
  "log10.Zostera.AG.mass"] = predict(ZAG.rf, y[is.na(y$log10.Zostera.AG.mass), ])
# sum(is.na(y$log10.Zostera.AG.mass)) # 0


# ZOSTERA BELOW-GROUND BIOMASS
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.Zostera.BG.mass)) # 29
ZBG.rf = randomForest(log10.Zostera.BG.mass ~ Ocean + Coast + Latitude + Longitude
                      # + sst.mean + Salinity.ppt + parmean + log10.nitrate + log10.phosphate
                      + log10.chlomean
                      # + log10.Zostera.BG.mass
                      + log10.Zostera.shoots.core
                      + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                      # + log10.periphyton.mass.per.g.zostera
                      # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                      + log10.Leaf.PercN
                      # + pop.density.2015
                      + AllelicRichness
                      ,
                      na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Zostera.BG.mass),
  "log10.Zostera.BG.mass"] = predict(ZBG.rf, y[is.na(y$log10.Zostera.BG.mass), ])
# sum(is.na(y$log10.Zostera.BG.mass)) # 0


# CRUSTACEAN mesograzer biomass
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.crustacean.mass.per.g.plant)) # 9
crust.rf = randomForest(log10.crustacean.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude
                        + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                        + log10.chlomean
                        + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                        + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                        # + log10.periphyton.mass.per.g.zostera
                        # # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                        + log10.Leaf.PercN
                        # # + pop.density.2015
                        + AllelicRichness
                        ,
                        na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.crustacean.mass.per.g.plant),
  "log10.crustacean.mass.per.g.plant"] = predict(crust.rf, y[is.na(y$log10.crustacean.mass.per.g.plant), ])
# sum(is.na(y$log10.crustacean.mass.per.g.plant)) # 0


# GASTROPOD mesograzer biomass
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.gastropod.mass.per.g.plant)) # 9
gast.rf = randomForest(log10.gastropod.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude
                       + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                       + log10.chlomean
                       + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                       + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                       # + log10.periphyton.mass.per.g.zostera
                       # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                       + log10.Leaf.PercN
                       # + pop.density.2015
                       + AllelicRichness
                       ,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.gastropod.mass.per.g.plant),
  "log10.gastropod.mass.per.g.plant"] = predict(gast.rf, y[is.na(y$log10.gastropod.mass.per.g.plant), ])
# sum(is.na(y$log10.gastropod.mass.per.g.plant)) # 0


# total MESOGRAZER biomass
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.mesograzer.mass.per.g.plant)) # 9
meso.rf = randomForest(log10.mesograzer.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude
                       + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                       + log10.chlomean
                       + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                       + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                       # + log10.periphyton.mass.per.g.zostera
                       # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                       + log10.Leaf.PercN
                       # + pop.density.2015
                       + AllelicRichness
                       ,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.mass.per.g.plant),
  "log10.mesograzer.mass.per.g.plant"] = predict(meso.rf, y[is.na(y$log10.mesograzer.mass.per.g.plant), ])
# sum(is.na(y$log10.mesograzer.mass.per.g.plant)) # 0


# total MESOGRAZER abundance
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.mesograzer.abund.per.g.plant)) # 9
meso.abund.rf = randomForest(log10.mesograzer.abund.per.g.plant ~ Ocean + Coast + Latitude + Longitude
                             + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                             + log10.chlomean
                             + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                             + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                             # + log10.periphyton.mass.per.g.zostera
                             # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                             + log10.Leaf.PercN
                             # + pop.density.2015
                             + AllelicRichness
                             ,
                             na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.abund.per.g.plant),
  "log10.mesograzer.abund.per.g.plant"] = predict(meso.abund.rf, y[is.na(y$log10.mesograzer.abund.per.g.plant), ])
# sum(is.na(y$log10.mesograzer.abund.per.g.plant)) # 0


# CRUSTACEAN mesograzer biomass PER AREA
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.crustacean.mass.per.area)) # 33
crust.area.rf = randomForest(log10.crustacean.mass.per.area ~ Ocean + Coast + Latitude + Longitude
                             + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                             + log10.chlomean
                             + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                             + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                             # + log10.periphyton.mass.per.g.zostera
                             # # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                             + log10.Leaf.PercN
                             # # + pop.density.2015
                             + AllelicRichness
                             ,
                             na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.crustacean.mass.per.area),
  "log10.crustacean.mass.per.area"] = predict(crust.area.rf, y[is.na(y$log10.crustacean.mass.per.area), ])
# sum(is.na(y$log10.crustacean.mass.per.area)) # 0


# GASTROPOD mesograzer biomass PER AREA
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.gastropod.mass.per.area)) # 33
gast.area.rf = randomForest(log10.gastropod.mass.per.area ~ Ocean + Coast + Latitude + Longitude
                            + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                            + log10.chlomean
                            + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                            + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                            # + log10.periphyton.mass.per.g.zostera
                            # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                            + log10.Leaf.PercN
                            # + pop.density.2015
                            + AllelicRichness
                            ,
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.gastropod.mass.per.area),
  "log10.gastropod.mass.per.area"] = predict(gast.area.rf, y[is.na(y$log10.gastropod.mass.per.area), ])
# sum(is.na(y$log10.gastropod.mass.per.area)) # 0


# total MESOGRAZER biomass PER AREA
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.mesograzer.mass.per.area)) # 33
meso.area.rf = randomForest(log10.mesograzer.mass.per.area ~ Ocean + Coast + Latitude + Longitude
                            + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                            + log10.chlomean
                            + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                            + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                            # + log10.periphyton.mass.per.g.zostera
                            # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                            + log10.Leaf.PercN
                            # + pop.density.2015
                            + AllelicRichness
                            ,
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.mass.per.area),
  "log10.mesograzer.mass.per.area"] = predict(meso.area.rf, y[is.na(y$log10.mesograzer.mass.per.area), ])
# sum(is.na(y$log10.mesograzer.mass.per.area)) # 0


# total MESOGRAZER abundance PER AREA
# Use random forest to impute missing values. First build predictive model:
# sum(is.na(d$log10.mesograzer.abund.per.area)) # 33
meso.abund.area.rf = randomForest(log10.mesograzer.abund.per.area ~ Ocean + Coast + Latitude + Longitude
                                  + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                                  + log10.chlomean
                                  + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                                  + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                                  # + log10.periphyton.mass.per.g.zostera
                                  # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                                  + log10.Leaf.PercN
                                  # + pop.density.2015
                                  + AllelicRichness
                                  ,
                                  na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.abund.per.area),
  "log10.mesograzer.abund.per.area"] = predict(meso.abund.area.rf, y[is.na(y$log10.mesograzer.abund.per.area), ])
# sum(is.na(y$log10.mesograzer.abund.per.area)) # 0


# Create data frame containing the imputed values and add them to the master data frame
# names(y)
imputed.values.y <- y[c("Unique.ID",  "log10.Zostera.shoots.core", "log10.Zostera.AG.mass",
                        "log10.Zostera.BG.mass", "log10.Leaf.PercN",
                        "log10.crustacean.mass.per.g.plant", "log10.crustacean.mass.per.area",
                        "log10.gastropod.mass.per.g.plant", "log10.gastropod.mass.per.area",
                        "log10.mesograzer.mass.per.g.plant", "log10.mesograzer.mass.per.area",
                        "log10.mesograzer.abund.per.g.plant", "log10.mesograzer.abund.per.area"
)]

# Rename imputed values
colnames(imputed.values.y)[2:13] <- c("log10.Zostera.shoots.core.imputed",
                                      "log10.Zostera.AG.mass.imputed",
                                      "log10.Zostera.BG.mass.imputed",
                                      "log10.Leaf.PercN.imputed",
                                      "log10.crustacean.mass.per.g.plant.imputed",
                                      "log10.crustacean.mass.per.area.imputed",
                                      "log10.gastropod.mass.per.g.plant.imputed",
                                      "log10.gastropod.mass.per.area.imputed",
                                      "log10.mesograzer.mass.per.g.plant.imputed",
                                      "log10.mesograzer.mass.per.area.imputed",
                                      "log10.mesograzer.abund.per.g.plant.imputed",
                                      "log10.mesograzer.abund.per.area.imputed"
)


# Integrate the imputed values back into master data set
d.imputed <- d

d.imputed$log10.Zostera.shoots.core.imputed <-
  imputed.values.y$log10.Zostera.shoots.core.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.Zostera.AG.mass.imputed <-
  imputed.values.y$log10.Zostera.AG.mass.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.Zostera.BG.mass.imputed <-
  imputed.values.y$log10.Zostera.BG.mass.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.Leaf.PercN.imputed <-
  imputed.values.y$log10.Leaf.PercN.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.crustacean.mass.per.g.plant.imputed <-
  imputed.values.y$log10.crustacean.mass.per.g.plant.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.crustacean.mass.per.area.imputed <-
  imputed.values.y$log10.crustacean.mass.per.area.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.gastropod.mass.per.g.plant.imputed <-
  imputed.values.y$log10.gastropod.mass.per.g.plant.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.gastropod.mass.per.area.imputed <-
  imputed.values.y$log10.gastropod.mass.per.area.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.mesograzer.mass.per.g.plant.imputed <-
  imputed.values.y$log10.mesograzer.mass.per.g.plant.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.mesograzer.mass.per.area.imputed <-
  imputed.values.y$log10.mesograzer.mass.per.area.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.mesograzer.abund.per.g.plant.imputed <-
  imputed.values.y$log10.mesograzer.abund.per.g.plant.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]

d.imputed$log10.mesograzer.abund.per.area.imputed <-
  imputed.values.y$log10.mesograzer.abund.per.area.imputed[match(d.imputed$Unique.ID, imputed.values.y$Unique.ID)]
#
# nrow(d.imputed) # 1000 - good


# PERIPHYTON

# NOTE: Here we need a reduced data set of 49 sites (i.e., excluding SW.A) for predictive model and imputation,
# because ALL plots from SW.A. had no periphyton values so it s not valid to impute values for that site.
# This requires two steps:

# First, we subset the dataframe of imputed values created above. After we derive imputed values for
# periphyton, we will paste them back inot this dataframe:
d.49_imputed <- droplevels(subset(d.imputed, Site != "SW.A"))

# Second: To rigorously estimate imputed values for a variable (in this case periphyton), we should use only
# empirical data, so we next subset the original (pre-imputation) dataframe for this purpose:
d.49 <- droplevels(subset(d, Site != "SW.A"))
x <- d.49

# PERIPHYTON mass per g Zostera
# Use random forest to impute missing values for periphyton in the reduced dataframe (49 sites).
# First build predictive model:
# sum(is.na(d.49$log10.periphyton.mass.per.g.zostera)) # 4
peri.rf = randomForest(log10.periphyton.mass.per.g.zostera ~ Ocean + Coast + Latitude + Longitude
                       + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                       + log10.chlomean
                       + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                       + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                       # # + log10.periphyton.mass.per.g.zostera
                       # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                       + log10.Leaf.PercN
                       # + pop.density.2015
                       + AllelicRichness
                       ,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = x)

# Impute missing values
x[is.na(x$log10.periphyton.mass.per.g.zostera),
  "log10.periphyton.mass.per.g.zostera"] = predict(peri.rf, x[is.na(x$log10.periphyton.mass.per.g.zostera), ])
# sum(is.na(x$log10.periphyton.mass.per.g.zostera)) # 0


# PERIPHYTON mass per AREA
# Use random forest to impute missing values for periphyton in the reduced dataframe (49 sites).
# First build predictive model:
# sum(is.na(d.49$log10.periphyton.mass.per.area)) # 28
peri.area.rf = randomForest(log10.periphyton.mass.per.area ~ Ocean + Coast + Latitude
                            + Longitude
                            + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                            + log10.chlomean
                            # + log10.Zostera.shoots.core
                            + log10.Zostera.sheath.length
                            # + log10.Zostera.longest.leaf.length
                            # # + log10.periphyton.mass.per.g.zostera
                            # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                            # + log10.Leaf.PercN
                            # + pop.density.2015
                            + AllelicRichness
                            ,
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = x)

# Impute missing values
x[is.na(x$log10.periphyton.mass.per.area),
  "log10.periphyton.mass.per.area"] = predict(peri.area.rf, x[is.na(x$log10.periphyton.mass.per.area), ])
# sum(is.na(x$log10.periphyton.mass.per.area)) # 0 (if it returns 3, then take out most predictors, and add back gradually)


# Create data frame containing the imputed value for periphyton and add to the 49-site data frame
imputed.values.x <- x[c("Unique.ID", "log10.periphyton.mass.per.g.zostera", "log10.periphyton.mass.per.area")]

# Rename imputed values
colnames(imputed.values.x)[2:3] <- c("log10.periphyton.mass.per.g.zostera.imputed", "log10.periphyton.mass.per.area.imputed" )

# Now paste the imputed values for periphyton back into the 49-site dataframe with the other
# imputed variables created above:
d.49_imputed$log10.periphyton.mass.per.g.zostera.imputed <-
  imputed.values.x$log10.periphyton.mass.per.g.zostera.imputed[match(d.49_imputed$Unique.ID, imputed.values.x$Unique.ID)]
# sum(is.na(d.49_imputed$log10.periphyton.mass.per.g.zostera.imputed)) # 0

d.49_imputed$log10.periphyton.mass.per.area.imputed <-
  imputed.values.x$log10.periphyton.mass.per.area.imputed[match(d.49_imputed$Unique.ID, imputed.values.x$Unique.ID)]
# sum(is.na(d.49_imputed$log10.periphyton.mass.per.area.imputed)) # 0
# nrow(d.imputed) # 1000 - good
# nrow(d.49_imputed) # 980 - good

# Summary: We can now build and compare models that will have same number of observations
# BUT can't include periphyton as predictor because biased by missing from entire site SW.A.


write_csv( d.imputed, "data/output/archive/ZEN_2014_imputed.csv" )
write_csv( d.49_imputed, "data/output/archive/ZEN_2014_imputed.csv" )