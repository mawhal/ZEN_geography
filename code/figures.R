###################################################################################
#                                                                                ##
# ZEN 2014: Global eelgrass ecosystem structure: Figures                         ##
# Data are current as of 2017-04-24                                              ##
# Emmett Duffy (duffye@si.edu)                                                   ##  
# Last updated 2021-05-31                                                        ##
#                                                                                ##
###################################################################################

# See this script for formatting pairs panels for publication graphs: ZEN 2014 model comparison Duffy 20191103.R

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# METADATA                                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# EXPLORE DATA DISTRIBUTIONS: GLOBAL                                              #
# FIGURE: BIVARIATE PLOTS OF EELGRASS FORM AND GENETICS AXES BY SITE              #
# FIGURE: MODEL COEFFICIENTS BY RESPONSE AND OCEAN                                #
# FIGURE: COMPARATIVE EFFECTS OF ENVIRONMENT AND EVO HISTORY (BAR CHARTS)         #
# FIGURE: ZEN SITE MAPS                                                           #
# FIGURES: BIVARIATE PLOTS (ATLANTIC)                                             #
# FIGURES: BIVARIATE PLOTS (PACIFIC)                                              #
#                                                                                 #
###################################################################################

###################################################################################
# METADATA                                                                        #
###################################################################################

# This script builds figures from ZEN 2014 global eelgrass survey. See also:
#   for preparation of the data: ZEN_2014_data_assembly_20210217.R (and subsequent)
#   for data exploration and model building: ZEN_2014_model_comparison_20210227.R (and subsequent)

# NOTE: I have deduced from exploring that partial.resid function (piecewiseSEM) returns 
# incorrect partial residual plots when the model includes interactions. Therefore, to make 
# partial residual plots correctly, I create a new variable = the interaction term, and refit 
# the model with that new interaction variable, then pliot rthe residuals. 

# For explanation of partial residuals and code for how I calculate them below, 
# see: https://rdrr.io/rforge/car/man/crPlots.html
# NOTE: NEED TO CREATE A FUNCTION FROM ALL THIS CODE SO IT CAN BE APPLIED TO ANY MODEL WITH ONE COMMAND


###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(psych) # for pairs panels
library(ggplot2)
library(nlme) # needed to run mixed models with lme
library(piecewiseSEM) # for partialResid
# library(car) # needed to extract model terms and partial residuals


###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# Read in zen2014 SITE-level data sets 
ZEN_2014_site_means <- read.csv("data/output/ZEN_2014_site_means.csv",  header = TRUE)
ZEN_2014_site_means_49 <- droplevels(subset(ZEN_2014_site_means, Site != "SW.A"))
ZEN_2014_site_means_49_Atlantic <- read.csv("data/output/ZEN_2014_site_means_49_Atlantic.csv",  header = TRUE)
ZEN_2014_site_means_Pacific <- read.csv("data/output/ZEN_2014_site_means_Pacific.csv",  header = TRUE)

# Read in zen2014 PLOT-level data sets for modeling (MINUS SW.A), with missing data imputed:
ZEN_2014_plot_49 <- read.csv("data/output/ZEN_2014_plot_49_noNA.csv", header = TRUE)
# ZEN_2014_plot_49_Atlantic <- read.csv("ZEN_2014_plot_49_Atlantic_20210227.csv", header = TRUE)
# ZEN_2014_plot_49_Pacific <- read.csv("ZEN_2014_plot_49_Pacific_20210227.csv", header = TRUE)


###################################################################################
# EXPLORE DATA DISTRIBUTIONS: GLOBAL                                              #
###################################################################################

ZEN_2014_site_means$Ocean <- as.factor(ZEN_2014_site_means$Ocean) # necessary to code pairs panels symbols by ocean ...

# Create new data frame with variable names friendlier to plotting
ZEN_2014_site_means_renamed <- ZEN_2014_site_means
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="sst.mean"] <- "SST mean"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="sst.range"] <- "SST range"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="Salinity.ppt"] <- "Salinity"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="parmean"] <- "PAR"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="cloudmean"] <- "Cloud cover"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="day.length"] <- "day length"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="sqrt.nitrate"] <- "NO3 (sqrt)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.phosphate"] <- "PO4 (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.chlomean"] <- "chl (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Leaf.PercN.site"] <- "Leaf %N (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC1.env.global"] <- "Env PCe1"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC2.env.global"] <- "Env PCe2"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC3.env.global"] <- "Env PCe3"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC1.zos.site"] <- "Eelgrass form (PCz1)"

ZEN_2014_site_means_renamed$Eelgrass.biomass <- -ZEN_2014_site_means_renamed$PC2.zos.site
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="Eelgrass.biomass"] <- "Eelgrass biomass (-PCz2)"

names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="FC1"] <- "Genetics FCA1"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="FC2"] <- "Genetics FCA2"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.AG.mass.site"] <- "AG mass (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.BG.mass.site"] <- "BG mass (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.shoots.core.site"] <- "Shoot density (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.sheath.length.site"] <- "Sheath L (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.sheath.width.site"] <- "Sheath W (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.longest.leaf.length.cm.site"] <- "Canopy ht (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.periphyton.mass.per.g.zostera.site"] <- "Periphyton mass (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.mesograzer.mass.per.g.plant.site"] <- "Mesograzer mass (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.periphyton.mass.per.area.site"] <- "Periphyton mass/area (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.mesograzer.mass.per.area.site"] <- "Mesograzer mass/area (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.mesograzer.abund.per.area.site"] <- "Mesograzer abund/area (log)"

ZEN_2014_site_means_renamed_Atlantic <- droplevels(subset(ZEN_2014_site_means_renamed, Ocean == "Atlantic"))
ZEN_2014_site_means_renamed_Pacific <- droplevels(subset(ZEN_2014_site_means_renamed, Ocean == "Pacific"))

### FIGURE S7
# Correlates of environmental/oceanographic PC axes 
png("figures/FigS7_pairs_correlates_of_envPCA.png", height = 12, width = 12, units = "in", res = 600)
pairs.panels(ZEN_2014_site_means_renamed[,c("Latitude", "Env PCe1", "Env PCe2", 
  "Env PCe3", "SST mean", "SST range", "Salinity", "PAR", "Cloud cover",
  "PO4 (log)", "Leaf %N (log)", "chl (log)")], 
  hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, 
  bg = c("blue","green")[ZEN_2014_site_means_renamed$Ocean])
dev.off()

### FIGURE S6
# Explore correlations among eelgrass characteristics by ocean 
png("figures/FigS6_pairs_correlates_of_growth_biomass.png", height = 8, width = 10, units = "in", res = 600)
pairs.panels(ZEN_2014_site_means_renamed[,c("Eelgrass form (PCz1)", "Eelgrass biomass (-PCz2)",
  "Canopy ht (log)", "Shoot density (log)", "AG mass (log)", "BG mass (log)", "Sheath L (log)", 
  "Sheath W (log)")], hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, 
  bg = c("blue","green")[ZEN_2014_site_means_renamed$Ocean])
dev.off()


# # Correlations between environment and biology
# pairs.panels(ZEN_2014_site_means_renamed[,c("Ocean", "Latitude", "Env PCe1", "Env PCe2", 
#   "Env PCe3", "Genetics FCA1", "Genetics FCA2", "Eelgrass form (PCz1)", "Eelgrass biomass (-PCz2)",
#   "Periphyton mass (log)", "Mesograzer mass (log)", "Periphyton mass/area (log)", 
#   "Mesograzer mass/area (log)")], hist.col="gray", pch = 21, 
#   smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 12, 
#   bg = c("blue","green")[ZEN_2014_site_means_renamed$Ocean])
# 
# # Atlantic: Correlations among covariates
# pairs.panels(ZEN_2014_site_means_renamed_Atlantic[,c("Latitude", "Env PCe1", "Env PCe2", "Env PCe3", 
#   "Genetics FCA1", "Genetics FCA2", "Eelgrass form (PCz1)", "Eelgrass biomass (-PCz2)",
#   "Periphyton mass (log)", "Mesograzer mass (log)")], hist.col="gray", pch = 21, 
#   smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 5, bg = c("blue")
#   )
# 
# # PACIFIC: Correlations among covariates
# pairs.panels(ZEN_2014_site_means_renamed_Pacific[,c("Latitude", "Env PCe1", "Env PCe2", "Env PCe3", 
#   "Genetics FCA1", "Genetics FCA2", "Eelgrass form (PCz1)", "Eelgrass biomass (-PCz2)",
#   "Periphyton mass (log)", "Mesograzer mass (log)")], hist.col="gray", pch = 21, 
#   smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 5, bg = c("forestgreen")
# )



###################################################################################
# FIGURE: BIVARIATE PLOTS OF EELGRASS FORM AND GENETICS AXES BY SITE              #
###################################################################################

#  Eelgrass morphology PC1 vs PC2
pc1z.pc2z <- ggplot(ZEN_2014_site_means, aes(x = PC1.zos.site, y = PC2.zos.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass morphology (PCz1)") +  
  ylab("Eelgrass morphology (PCz2)") +  
  scale_x_continuous(limits = c(-4,4.3), breaks = c(-4,-2,0,2,4)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) 
pc1z.pc2z
ggsave(pc1z.pc2z, filename = "figures/Fig1b_pc1z_pc2z.png",  width = 6.5, height = 6 ) #, bg = "transparent"


# Eelgrass genetics FC1 vs FC2
fc1.fc2 <- ggplot(ZEN_2014_site_means, aes(x = FC1, y = FC2, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass genetics (FCA1)") +  
  ylab("Eelgrass genetics (FCA2)") +  
  scale_x_continuous(limits = c(-1100,1100), breaks = c(-1000,0,1000), labels =c("-1K","0","1K")) +
  scale_y_continuous(limits=c(-1000, 2000), breaks = c(-1000,0,1000,2000), labels = c("-1","0","1", "2")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) 
fc1.fc2
ggsave(fc1.fc2, filename = "figures/Fig1c_fc1_fc2.png",  width = 6.5, height = 6 ) #, bg = "transparent"


#  Environment PC1 vs PC2
# calculate convex hull
library(tidyverse)
hull_PC12 <- ZEN_2014_site_means %>%
  group_by(Ocean) %>%
  slice(chull(PC1.env.global, PC2.env.global))
pce1.pce2 <- ggplot(ZEN_2014_site_means, aes(x = PC1.env.global, y = PC2.env.global, group = Ocean, col = Ocean, fill = Ocean)) +
  geom_polygon(data = hull_PC12, alpha = 0.1 ) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  scale_fill_manual(values = c("blue", "forestgreen")) +
  xlab("Environment (PCe1)") +
  ylab("Environment (PCe2)") +
  scale_x_continuous(limits = c(-2.75,4)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  )
pce1.pce2
ggsave(pce1.pce2, filename = "figures/FigS8a_pce1_pce2.png",  width = 6.5, height = 6 ) 


#  Environment PC2 vs PC3
hull_PC23 <- ZEN_2014_site_means %>%
  group_by(Ocean) %>%
  slice(chull(PC2.env.global, PC3.env.global))
pce2.pce3 <- ggplot(ZEN_2014_site_means, aes(x = PC2.env.global, y = PC3.env.global, group = Ocean, col = Ocean, fill = Ocean)) +
  geom_polygon(data = hull_PC23, alpha = 0.1 ) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  scale_fill_manual(values = c("blue", "forestgreen")) +
  xlab("Environment (PCe2)") +
  ylab("Environment (PCe3)") +
  scale_x_continuous(limits = c(-3.75,3.25)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  )
pce2.pce3
ggsave(pce2.pce3, filename = "figures/FigS8b_pce2_pce3.png",  width = 6.5, height = 6 ) 

# # Eelgrass Morphology PC1: Effect of FC1 (RAW) 
# pc1z.fc1 = ggplot(ZEN_2014_site_means, aes(x = FC1, y = PC1.zos.site, group = Ocean, col = Ocean)) +
#   geom_point(size = 4) +
#   # geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("blue", "forestgreen")) +
#   xlab("Eelgrass genetics (FCA1)") +  
#   ylab("Eelgrass form (PCz1)") +  
#   # scale_x_continuous(limits = c(-3,3)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
#   geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T)


# Eelgrass Morphology PC1: Effect of FC2 (PARTIAL, SITE LEVEL)
# First use best model adapted to site level 
pcz1.site <- lm(zPC1.zos.site ~ 
                      + zPC1.env.global + zPC2.env.global + zPC3.env.global + zFC2 + zFC1  
                    , data = ZEN_2014_site_means)

pc1z.fc2 = ggplot(ZEN_2014_site_means, aes(x = FC2, y = PC1.zos.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  # geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass genetics (FCA2)") +  
  ylab("Eelgrass form (PCz1)") +  
  # scale_x_continuous(limits = c(-3,3)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T)
pc1z.fc2


#  Eelgrass shoot density vs canopy height by ocean
shoots.canopy <- ggplot(ZEN_2014_site_means, aes(x = Zostera.shoots.core.site, y = Zostera.longest.leaf.length.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Shoot density") +  
  ylab("Canopy height") +  
  # scale_x_continuous(limits = c(-4,4.3), breaks = c(-4,-2,0,2,4)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
shoots.canopy
# Slope  = zero for Atlantic, negative in Pacific. Thus Pacific shows trade-off between growing 
# longer shoots vs new shoots, whereas in Atlantic more shoots means more biomass,


#  Eelgrass shoot density vs above-ground biomass by ocean
shoots.zag <- ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.shoots.core.site, y = log10.Zostera.AG.mass.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Shoot density") +  
  ylab("Eelgrass above-ground biomass") +  
  # scale_x_continuous(limits = c(-4,4.3), breaks = c(-4,-2,0,2,4)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
shoots.zag



###################################################################################
# FIGURE: MODEL COEFFICIENTS BY RESPONSE AND OCEAN                                #
###################################################################################

# ATLANTIC 

# Read in GLM model coefficients produced by modeling script (model_comparison.R)
coeffs.atl <- read.csv("data/output/zen_glm_coefficients_atlantic.csv", header = TRUE)
names(coeffs.atl)
coeffs.atl$predictor <- as.factor(coeffs.atl$predictor)

# Eelgrass growth form 1
coeffs.atl.pcz1 <- droplevels(subset(coeffs.atl, response == "zPC1.zos.global.atl"))
coeffs.atl.pcz1$predictor.order <- factor(coeffs.atl.pcz1$predictor, as.character(coeffs.atl.pcz1$predictor))
pcz1.atl.coeffs.a <- ggplot(coeffs.atl.pcz1, mapping = aes(x = predictor.order, y = coefficient, 
                                                           ymin = (coefficient - se), ymax = (coefficient + se))) +
  geom_pointrange() + coord_flip() + labs(x="", y="effect size") +
  geom_point(size = 4, shape = 21, fill = "blue") +
  scale_x_discrete(limits = rev(levels(coeffs.atl.pcz1$predictor.order)),
                   # labels = c("FC1 x PC3", "Environment PC3", "Environment PC2", "Environment PC2", "Genetic FC2", "Genetic FC1")
  ) +
  scale_y_continuous(limits=c(-1.5, 1.75)) +
  geom_hline(yintercept = 0) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme_bw() +
  ggtitle("Eelgrass form PCz1 (Atlantic)") +theme(plot.title = element_text(hjust = 0.5))
pcz1.atl.coeffs.a


# Eelgrass growth form PCz2 (ATLANTIC)
coeffs.atl.pcz2 <- droplevels(subset(coeffs.atl, response == "zPC2.zos.global.atl"))
coeffs.atl.pcz2$predictor.order <- factor(coeffs.atl.pcz2$predictor, as.character(coeffs.atl.pcz2$predictor))
pcz2.atl.coeffs.a <- ggplot(coeffs.atl.pcz2, mapping = aes(x = predictor.order, y = coefficient, 
                                                           ymin = (coefficient - se), ymax = (coefficient + se))) +
  geom_pointrange() + coord_flip() + labs(x="", y="effect size") +
  geom_point(size = 4, shape = 21, fill = "blue") +
  scale_x_discrete(limits = rev(levels(coeffs.atl.pcz2$predictor.order)),
                   # labels = c("FC1 x PC3", "Environment PC3", "Environment PC2", "Environment PC2", "Genetic FC2", "Genetic FC1")
  ) +
  scale_y_continuous(limits=c(-1.5, 1.75)) +
  geom_hline(yintercept = 0) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme_bw() +
  ggtitle("Eelgrass form PCz2 (Atlantic)") +theme(plot.title = element_text(hjust = 0.5))
pcz2.atl.coeffs.a

# Periphyton (ATLANTIC)
coeffs.atl.peri <- droplevels(subset(coeffs.atl, response == "zperi.perg.atl"))
coeffs.atl.peri$predictor.order <- factor(coeffs.atl.peri$predictor, as.character(coeffs.atl.peri$predictor))
peri.atl.coeffs.a <- ggplot(coeffs.atl.peri, mapping = aes(x = predictor.order, y = coefficient, 
                                                           ymin = (coefficient - se), ymax = (coefficient + se))) +
  geom_pointrange() + coord_flip() + labs(x="", y="effect size") +
  geom_point(size = 4, shape = 21, fill = "blue") +
  scale_x_discrete(limits = rev(levels(coeffs.atl.peri$predictor.order)),
                   # labels = c("FC1 x PC3", "Environment PC3", "Environment PC2", "Environment PC2", "Genetic FC2", "Genetic FC1")
  ) +
  scale_y_continuous(limits=c(-1.5, 1.75)) +
  geom_hline(yintercept = 0) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme_bw() +
  ggtitle("Periphyton mass per g (Atlantic)") +theme(plot.title = element_text(hjust = 0.5))
peri.atl.coeffs.a


# Mesograzers (ATLANTIC)
coeffs.atl.meso <- droplevels(subset(coeffs.atl, response == "zmeso.perg.atl"))
coeffs.atl.meso$predictor.order <- factor(coeffs.atl.meso$predictor, as.character(coeffs.atl.meso$predictor))
meso.atl.coeffs.a <- ggplot(coeffs.atl.meso, mapping = aes(x = predictor.order, y = coefficient, 
                                                           ymin = (coefficient - se), ymax = (coefficient + se))) +
  geom_pointrange() + coord_flip() + labs(x="", y="effect size") +
  geom_point(size = 4, shape = 21, fill = "blue") +
  scale_x_discrete(limits = rev(levels(coeffs.atl.meso$predictor.order)),
                   # labels = c("FC1 x PC3", "Environment PC3", "Environment PC2", "Environment PC2", "Genetic FC2", "Genetic FC1")
  ) +
  scale_y_continuous(limits=c(-1.5, 1.75)) +
  geom_hline(yintercept = 0) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme_bw() +
  ggtitle("Mesograzer mass per g (Atlantic)") +theme(plot.title = element_text(hjust = 0.5))
meso.atl.coeffs.a


# PACIFIC 

# Read in GLM model coefficients
coeffs.pac <- read.csv("data/output/zen_glm_coefficients_pacific.csv", header = TRUE)
names(coeffs.pac)
coeffs.pac$predictor <- as.factor(coeffs.pac$predictor)


# Eelgrass growth form PCz1 (PACIFIC)
coeffs.pac.pcz1 <- droplevels(subset(coeffs.pac, response == "zPC1.zos.global.pac"))
coeffs.pac.pcz1$predictor.order <- factor(coeffs.pac.pcz1$predictor, as.character(coeffs.pac.pcz1$predictor))
pcz1.pac.coeffs.p <- ggplot(coeffs.pac.pcz1, mapping = aes(x = predictor.order, y = coefficient, 
                                                           ymin = (coefficient - se), ymax = (coefficient + se))) +
  geom_pointrange() + coord_flip() + labs(x="", y="effect size") +
  geom_point(size = 4, shape = 21, fill = "forestgreen") +
  scale_x_discrete(limits = rev(levels(coeffs.pac.pcz1$predictor.order)),
                   # labels = c("FC1 x PC3", "Environment PC3", "Environment PC2", "Environment PC2", "Genetic FC2", "Genetic FC1")
  ) +
  scale_y_continuous(limits=c(-1.5, 1.75)) +
  geom_hline(yintercept = 0) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme_bw() +
  ggtitle("Eelgrass form PCz1 (Pacific)") +theme(plot.title = element_text(hjust = 0.5))
pcz1.pac.coeffs.p


# Eelgrass growth form PCz2 (PACIFIC)
coeffs.pac.pcz2 <- droplevels(subset(coeffs.pac, response == "zPC2.zos.global.pac"))
coeffs.pac.pcz2$predictor.order <- factor(coeffs.pac.pcz2$predictor, as.character(coeffs.pac.pcz2$predictor))
pcz2.pac.coeffs.p <- ggplot(coeffs.pac.pcz2, mapping = aes(x = predictor.order, y = coefficient, 
                                                           ymin = (coefficient - se), ymax = (coefficient + se))) +
  geom_pointrange() + coord_flip() + labs(x="", y="effect size") +
  geom_point(size = 4, shape = 21, fill = "forestgreen") +
  scale_x_discrete(limits = rev(levels(coeffs.pac.pcz2$predictor.order)),
                   # labels = c("FC1 x PC3", "Environment PC3", "Environment PC2", "Environment PC2", "Genetic FC2", "Genetic FC1")
  ) +
  scale_y_continuous(limits=c(-1.5, 1.75)) +
  geom_hline(yintercept = 0) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme_bw() +
  ggtitle("Eelgrass form PCz2 (Pacific)") +theme(plot.title = element_text(hjust = 0.5))
pcz2.pac.coeffs.p


# Periphyton (PACIFIC)
coeffs.pac.peri <- droplevels(subset(coeffs.pac, response == "zperi.perg.pac"))
coeffs.pac.peri$predictor.order <- factor(coeffs.pac.peri$predictor, as.character(coeffs.pac.peri$predictor))
peri.pac.coeffs.p <- ggplot(coeffs.pac.peri, mapping = aes(x = predictor.order, y = coefficient, 
                                                           ymin = (coefficient - se), ymax = (coefficient + se))) +
  geom_pointrange() + coord_flip() + labs(x="", y="effect size") +
  geom_point(size = 4, shape = 21, fill = "forestgreen") +
  scale_x_discrete(limits = rev(levels(coeffs.pac.peri$predictor.order)),
                   # labels = c("FC1 x PC3", "Environment PC3", "Environment PC2", "Environment PC2", "Genetic FC2", "Genetic FC1")
  ) +
  scale_y_continuous(limits=c(-1.5, 1.75)) +
  geom_hline(yintercept = 0) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme_bw() +
  ggtitle("Periphyton mass per g (Pacific)") +theme(plot.title = element_text(hjust = 0.5))
peri.pac.coeffs.p


# Mesograzers (PACIFIC)
coeffs.pac.meso <- droplevels(subset(coeffs.pac, response == "zmeso.perg.pac"))
coeffs.pac.meso$predictor.order <- factor(coeffs.pac.meso$predictor, as.character(coeffs.pac.meso$predictor))
meso.pac.coeffs.p <- ggplot(coeffs.pac.meso, mapping = aes(x = predictor.order, y = coefficient, 
                                                           ymin = (coefficient - se), ymax = (coefficient + se))) +
  geom_pointrange() + coord_flip() + labs(x="", y="effect size") +
  geom_point(size = 4, shape = 21, fill = "forestgreen") +
  scale_x_discrete(limits = rev(levels(coeffs.pac.meso$predictor.order)),
                   # labels = c("FC1 x PC3", "Environment PC3", "Environment PC2", "Environment PC2", "Genetic FC2", "Genetic FC1")
  ) +
  scale_y_continuous(limits=c(-1.5, 1.75)) +
  geom_hline(yintercept = 0) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme_bw() +
  ggtitle("Mesograzer mass per g (Pacific)") +theme(plot.title = element_text(hjust = 0.5))
meso.pac.coeffs.p



###################################################################################
# FIGURE: COMPARATIVE EFFECTS OF ENVIRONMENT AND EVO HISTORY (BAR CHARTS)         #
###################################################################################

# To compare influences of environment and evolutionary history,calculate the weighted 
# SUM of standardized effects of each (e.g., PCe1, PCe2, PCe3 for environment effect). 


library(dplyr)          # for data manipulation
library(tidyr)          # for data manipulation
library(magrittr)       # for easier syntax in one or two areas
library(gridExtra)      # for generating some comparison plots
library(ggplot2)        # for generating the visualizations


# Read in standardized coefficients from best models: ATLANTIC
zen_best_models_atlantic <- read.csv("data/evol_envir_path_calculations/zen_2014_model_coeffs_atlantic.csv",header = TRUE)
names(zen_best_models_atlantic)

zen_best_models_atlantic$response <- as.factor(zen_best_models_atlantic$response)
zen_best_models_atlantic$predictor <- as.factor(zen_best_models_atlantic$predictor)
zen_best_models_atlantic$effect.type <- as.factor(zen_best_models_atlantic$effect.type)

# reorder levels of variables
zen_best_models_atlantic$response <- factor(zen_best_models_atlantic$response, levels = c("eelgrass.form", "eelgrass.biomass", "periphyton", "mesograzers"))

# Create facet plot: Sums of variance-weighted effects, direct and indirect 
summed_paths_atl_high_scale <- ggplot(zen_best_models_atlantic, aes(factor(predictor), weighted.sum.abs, fill = factor(effect.type, levels = c("direct", "indirect")))) +
  geom_bar(stat = "identity", color = "grey40") +
  scale_fill_manual(values = c("blue", "#c6eafb")) +
  labs(fill = "") +
  scale_x_discrete(labels = c("Environment", "History")) +
  scale_y_continuous(limits = c(0, 3.00), breaks = c(0,1,2,3)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        legend.position = "top"
        # ,
        # panel.background = element_rect(fill = "transparent"), # bg of the panel
        # plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  ) +
  ylab("Weighted mean of standardized effects\n") + xlab("") +
  facet_grid(response ~ .) 
summed_paths_atl_high_scale
png("figures/Fig4ag_summed_paths_atl_high_scale.png", width = 3, height = 10, unit = "in", res = 600)
summed_paths_atl_high_scale
dev.off() 


# Read in standardized coefficients from best models: PACIFIC
zen_best_models_pacific <- read.csv("data/evol_envir_path_calculations/zen_2014_model_coeffs_pacific.csv",header = TRUE)
names(zen_best_models_pacific)

zen_best_models_pacific$response <- as.factor(zen_best_models_pacific$response)
zen_best_models_pacific$predictor <- as.factor(zen_best_models_pacific$predictor)

# reorder levels of variables
zen_best_models_pacific$response <- factor(zen_best_models_atlantic$response, levels = c("eelgrass.form", "eelgrass.biomass", "periphyton", "mesograzers"))

# Create facet plot: Sums of variance-weighted effects, direct and indirect 
summed_paths_pac_high_scale <- ggplot(zen_best_models_pacific, aes(factor(predictor), weighted.sum.abs, fill = factor(effect.type, levels = c("direct", "indirect")))) +
  geom_bar(stat = "identity", color = "grey40") +
  scale_fill_manual(values = c("forestgreen", "#99ff66")) +
  labs(fill = "") +
  scale_x_discrete(labels = c("Environment", "History")) +
  scale_y_continuous(limits = c(0, 3.00), breaks = c(0,1,2,3)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        legend.position = "top"
        # ,
        # panel.background = element_rect(fill = "transparent"), # bg of the panel
        # plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  ) +
  ylab("Weighted mean of standardized effects\n") + xlab("") +
  facet_grid(response ~ .) 
summed_paths_pac_high_scale
png("figures/Fig4bh_summed_paths_pac_high_scale.png", width = 3, height = 10, unit = "in", res = 600)
summed_paths_pac_high_scale
dev.off() 









###################################################################################
# FIGURES: BIVARIATE PLOTS (ATLANTIC)                                             #
###################################################################################

# ATLANTIC: EELGRASS GROWTH FORM (PCz1) vs ENVIRONMENT PCe1 (Partial, site level)
# Refit Best model with separate variable for interaction term
ZEN_2014_site_means_49_Atlantic$rPCe3.rFC1.atl <- ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl * ZEN_2014_site_means_49_Atlantic$rFC1.global.atl
pcz1.site.a.4.x <- lm(rPC1.zos.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl 
  + rFC1.global.atl + rFC2.global.atl + rPCe3.rFC1.atl
  , data = ZEN_2014_site_means_49_Atlantic)
# Extract residuals from model
pcz1.pce1.a.partial <- partialResid(rPC1.zos.atl ~ rPC1.env.global.atl, pcz1.site.a.4.x, ZEN_2014_site_means_49_Atlantic) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49_Atlantic$Ocean
Site <- ZEN_2014_site_means_49_Atlantic$Site
pcz1.pce1.a.partial <- cbind(pcz1.pce1.a.partial, Ocean, Site)
pcz1.pce1.a.partial.plot <- ggplot(pcz1.pce1.a.partial, aes(x = (-xresid), y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Latitude/Climate (-PCe1) | Z\nAtlantic") +  
  ylab("Eelgrass form (PCz1) | Z\nAtlantic") +  
  # scale_x_continuous(limits = c(30,75)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz1.pce1.a.partial.plot


# ATLANTIC: EELGRASS GROWTH FORM (PCz1) vs GENETIC FC2 (Partial, site level)
# Extract residuals from model
pcz1.fc2.a.partial <- partialResid(rPC1.zos.atl ~ rFC2.global.atl, pcz1.site.a.4.x, ZEN_2014_site_means_49_Atlantic) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49_Atlantic$Ocean
Site <- ZEN_2014_site_means_49_Atlantic$Site
pcz1.fc2.a.partial <- cbind(pcz1.fc2.a.partial, Ocean, Site)
pcz1.fc2.a.partial.plot <- ggplot(pcz1.fc2.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass genetics (FCA2) | Z\nAtlantic") +  
  ylab("Eelgrass form (PCz1) | Z\nAtlantic") +  
  # scale_x_continuous(limits = c(30,75)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz1.fc2.a.partial.plot
ggsave("figures/Fig2a_FCA2_PCz1_atlantic.png", width = 5.5, height = 5)

# NOTE: Use inverse of PCz2 so positive values indicate high biomass


# ATLANTIC: EELGRASS BIOMASS (inverse PCz2) vs ENVIRONMENT PCe1 (Partial, site level)
# Best model at site level
ZEN_2014_site_means_49_Atlantic$rPCe1.rFC2.atl <- ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl * ZEN_2014_site_means_49_Atlantic$rFC2.global.atl
pcz2.test.a.2.x <- lm(rPC2.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPCe1.rFC2.atl
                    , data = ZEN_2014_site_means_49_Atlantic)
# Extract residuals from model
pcz2.pce1.a.partial <- partialResid(rPC2.zos.atl ~ rPC1.env.global.atl, pcz2.test.a.2.x, ZEN_2014_site_means_49_Atlantic) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49_Atlantic$Ocean
Site <- ZEN_2014_site_means_49_Atlantic$Site
pcz2.pce1.a.partial <- cbind(pcz2.pce1.a.partial, Ocean, Site)
pcz2.pce1.a.partial.plot <- ggplot(pcz2.pce1.a.partial, aes(x = -xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Latitude/Climate (-PCe1) | Z\nAtlantic") +  
  ylab("Eelgrass biomass (-PCz2) | Z\nAtlantic") +  
  # scale_x_continuous(limits = c(-2.3,2.1)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2.pce1.a.partial.plot
ggsave("figures/Fig2b_FCA2_PCz2_atlantic.png", width = 5.5, height = 5)

# ATLANTIC: EELGRASS BIOMASS (inverse PCz2) vs ENVIRONMENT PCe3 (Partial, site level)
# Extract residuals from model
pcz2.pce3.a.partial <- partialResid(rPC2.zos.atl ~ rPC3.env.global.atl, pcz2.test.a.2.x, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
pcz2.pce3.a.partial <- cbind(pcz2.pce3.a.partial, Ocean, Site)
pcz2.pce3.a.partial.plot <- ggplot(pcz2.pce3.a.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Estuarine conditions (PCe3) | Z\nAtlantic") +  
  ylab("Eelgrass biomass (-PCz2) | Z\nAtlantic") +  
  # scale_x_continuous(limits = c(-2.3,2.1)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2.pce3.a.partial.plot


# ATLANTIC: EELGRASS BIOMASS (inverse PCz2) vs GENETICS FC2 (Partial, site level)
# Extract residuals from model
pcz2.fc2.a.partial <- partialResid(rPC2.zos.atl ~ rFC2.global.atl, pcz2.test.a.2.x, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
pcz2.fc2.a.partial <- cbind(pcz2.fc2.a.partial, Site)
pcz2.fc2.a.partial.plot <- ggplot(pcz2.fc2.a.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass genetics (FCA2) | Z\nAtlantic") +  
  ylab("Eelgrass biomass (-PCz2) | Z\nAtlantic") +  
  # scale_x_continuous(limits = c(-.14,.3)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2.fc2.a.partial.plot


# ATLANTIC: PERIPHYTON PER G EELGRASS vs EELGRASS FORM PCz1 (Partial, site level)
# First refit best model with separate variable for interaction 
ZEN_2014_site_means_49_Atlantic$rPCz1.rPCe2.atl <- ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl*ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl
peri.test.a.9.x <- lm(rperiphyton.perg.atl ~ 
                        + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                      + rPC1.zos.atl + rPC2.zos.atl + rPCz1.rPCe2.atl # added
                      , data = ZEN_2014_site_means_49_Atlantic)
peri.pcz1.a.partial <- partialResid(rperiphyton.perg.atl ~ rPC1.zos.atl, peri.test.a.9.x, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
peri.pcz1.a.partial <- cbind(peri.pcz1.a.partial, Site)
peri.pcz1.a.partial.plot <- ggplot(peri.pcz1.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass form (PCz1) | Z\nAtlantic") +  
  ylab("Periphyton mass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-.25,.2)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.pcz1.a.partial.plot


# ATLANTIC: PERIPHYTON PER AREA vs EELGRASS FORM PCz1 (Partial, site level)
# First refit best model with separate variable for interaction 
ZEN_2014_site_means_49_Atlantic$rPCz1.rPCe2.atl <- ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl*ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl
peri.area.test.a.9 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPCz1.rPCe2.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)
peri.pcz1.a.partial <- partialResid(rperiphyton.area.atl ~ rPC1.zos.atl, peri.area.test.a.9, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
peri.pcz1.a.partial <- cbind(peri.pcz1.a.partial, Site)
peri.pcz1.a.partial.plot <- ggplot(peri.pcz1.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass form (PCz1) | Z\nAtlantic") +  
  ylab("Periphyton mass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-.25,.2)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.pcz1.a.partial.plot
ggsave("figures/Fig3a_PCz1_periphyton_atlantic.png", width = 5.5, height = 5)


# ATLANTIC: PERIPHYTON PER AREA vs EELGRASS BIOMASS (-PCz2) (Partial, site level)
peri.pcz2.a.partial <- partialResid(rperiphyton.area.atl ~ rPC2.zos.atl, peri.area.test.a.9, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
peri.pcz2.a.partial <- cbind(peri.pcz2.a.partial, Site)
peri.pcz2.a.partial.plot <- ggplot(peri.pcz2.a.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass biomass (-PCz2) | Z\nAtlantic") +  
  ylab("Periphyton mass | Z\nAtlantic") +  
  # scale_x_continuous(limits = c(-.25,.2)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.pcz2.a.partial.plot
ggsave("figures/Fig3b_PCz2_periphyton_atlantic.png", width = 5.5, height = 5)


# ATLANTIC: PERIPHYTON PER AREA vs NUTRIENT STATUS (PCe2) (Partial, site level)
# First refit best model with separate variable for interaction 
ZEN_2014_site_means_49_Atlantic$rPCz1.rPCe2.atl <- ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl*ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl
peri.area.test.a.9 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPCz1.rPCe2.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)
peri.pce2.a.partial <- partialResid(rperiphyton.area.atl ~ rPC2.env.global.atl, peri.area.test.a.9, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
peri.pce2.a.partial <- cbind(peri.pce2.a.partial, Site)
peri.pce2.a.partial.plot <- ggplot(peri.pce2.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Nutrient status (PCe2) | Z\nAtlantic") +  
  ylab("Periphyton mass | Z\nAtlantic") +  
  # scale_x_continuous(limits = c(-.25,.2)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.pce2.a.partial.plot


# ATLANTIC: PERIPHYTON PER G EELGRASS vs EELGRASS BIOMASS PCz2 (Partial, site level)
# First refit best model with separate variable for interaction 
ZEN_2014_site_means_49_Atlantic$rPCz1.rPCe2.atl <- ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl*ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl
peri.test.a.9.x <- lm(rperiphyton.perg.atl ~ 
                        + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                      + rPC1.zos.atl + rPC2.zos.atl + rPCz1.rPCe2.atl # added
                      , data = ZEN_2014_site_means_49_Atlantic)
peri.pcz2.a.partial <- partialResid(rperiphyton.perg.atl ~ rPC2.zos.atl, peri.test.a.9.x, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
peri.pcz2.a.partial <- cbind(peri.pcz2.a.partial, Site)
peri.pcz2.a.partial.plot <- ggplot(peri.pcz2.a.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass biomass (PCz2) | Z\nAtlantic") +  
  ylab("Periphyton mass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-.25,.2)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
  # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.pcz2.a.partial.plot


# ATLANTIC: MESOGRAZER BIOMASS PER G PLANT vs ENVIRONMENT PCe3 (Partial, site level)
# Best model at site level
ZEN_2014_site_means_49_Atlantic$rPCz1.rPCe2.atl <- ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl*ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl
meso.test.a.9.x <- lm(rmesograzer.mass.perg.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl + rPCz1.rPCe2.atl
  , data = ZEN_2014_site_means_49_Atlantic)
# Extract residuals from model
meso.pce3.a.partial <- partialResid(rmesograzer.mass.perg.atl ~ rPC3.env.global.atl, meso.test.a.9.x, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
meso.pce3.a.partial <- cbind(meso.pce3.a.partial, Site)
meso.pce3.a.partial.plot <- ggplot(meso.pce3.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Estuarine conditions (PCe3) | Z\nAtlantic") +  
  ylab("Mesograzer biomass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-0.55,0.55)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.pce3.a.partial.plot


# ATLANTIC: MESOGRAZER BIOMASS PER G PLANT vs GENETICS FC2 (Partial, site level)
# Extract residuals from model
meso.fc2.a.partial <- partialResid(rmesograzer.mass.perg.atl ~ rFC2.global.atl, meso.test.a.9.x, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
meso.fc2.a.partial <- cbind(meso.fc2.a.partial, Site)
meso.fc2.a.partial.plot <- ggplot(meso.fc2.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass genetics (FCA2) | Z\nAtlantic") +  
  ylab("Mesograzer biomass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-0.2,0.25)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.fc2.a.partial.plot


# # ATLANTIC: MESOGRAZER BIOMASS PER G PLANT vs MEADOW FORM (PCz1) (Partial, site level)
# # Extract residuals from model
# meso.pcz1.a.partial <- partialResid(rmesograzer.mass.perg.atl ~ rPC1.zos.atl, meso.test.a.9.x, ZEN_2014_site_means_49_Atlantic) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_49_Atlantic$Site
# meso.pcz1.a.partial <- cbind(meso.pcz1.a.partial, Site)
# meso.pcz1.a.partial.plot <- ggplot(meso.pcz1.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("blue")) +
#   xlab("Eelgrass form (PCz1) | Z\nAtlantic") +  
#   ylab("Mesograzer biomass | Z\nAtlantic") +  
#   scale_x_continuous(limits = c(-0.25,0.2)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
#   geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# meso.pcz1.a.partial.plot
# 
# 
# # ATLANTIC: MESOGRAZER BIOMASS PER G PLANT vs EELGRASS BIOMASS (-PCz2) (Partial, site level)
# # Extract residuals from model
# meso.perg.pcz2.a.partial <- partialResid(rmesograzer.mass.perg.atl ~ rPC2.zos.atl, meso.test.a.9.x, ZEN_2014_site_means_49_Atlantic) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_49_Atlantic$Site
# meso.perg.pcz2.a.partial <- cbind(meso.perg.pcz2.a.partial, Site)
# meso.perg.pcz2.a.partial.plot <- ggplot(meso.perg.pcz2.a.partial, aes(x = xresid, y = -yresid, col = "blue")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("blue")) +
#   xlab("Eelgrass biomass (PCz2) | Z\nAtlantic") +  
#   ylab("Mesograzer biomass | Z\nAtlantic") +  
#   scale_x_continuous(limits = c(-0.35,0.38)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) # +
#   # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# meso.perg.pcz2.a.partial.plot


# ATLANTIC: MESOGRAZER BIOMASS PER AREA vs MEADOW FORM (PCz1) (Partial, site level)
# Best model at site level
meso.area.test.a.1 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl 
                         + rPC1.zos.atl + rPC2.zos.atl  
                         , data = ZEN_2014_site_means_49_Atlantic)
# Extract residuals from model
meso.pcz1.a.partial <- partialResid(rmesograzer.mass.area.atl ~ rPC1.zos.atl, meso.area.test.a.1, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
meso.pcz1.a.partial <- cbind(meso.pcz1.a.partial, Site)
meso.pcz1.a.partial.plot <- ggplot(meso.pcz1.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass form (PCz1) | Z\nAtlantic") +  
  ylab("Invertebrate biomass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-0.3,0.26)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.pcz1.a.partial.plot
ggsave("figures/Fig3c_PCz1_invertebrates_atlantic.png", width = 5.5, height = 5)


# ATLANTIC: MESOGRAZER BIOMASS PER AREA vs EELGRASS BIOMASS (-PCz2) (Partial, site level)
# Extract residuals from model
meso.area.pcz2.a.partial <- partialResid(rmesograzer.mass.area.atl ~ rPC2.zos.atl, meso.area.test.a.1, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
meso.area.pcz2.a.partial <- cbind(meso.area.pcz2.a.partial, Site)
meso.area.pcz2.a.partial.plot <- ggplot(meso.area.pcz2.a.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass biomass (-PCz2) | Z\nAtlantic") +  
  ylab("Invertebrate biomass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-0.35,0.29)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.area.pcz2.a.partial.plot
ggsave("figures/Fig3d_PCz2_invertebrates_atlantic.png", width = 5.5, height = 5)


# # ATLANTIC: MESOGRAZER BIOMASS PER AREA vs ESTUARINE CONDITIONS (PCe3) (Partial, site level)
# # Extract residuals from model
# meso.area.pce3.a.partial <- partialResid(rmesograzer.mass.area.atl ~ rPC3.env.global.atl, meso.area.test.a.1, ZEN_2014_site_means_49_Atlantic) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_49_Atlantic$Site
# meso.area.pce3.a.partial <- cbind(meso.area.pce3.a.partial, Site)
# meso.area.pce3.a.partial.plot <- ggplot(meso.area.pce3.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("blue")) +
#   xlab("Estuarine conditions (PCe3) | Z\nAtlantic") +  
#   ylab("Mesograzer biomass | Z\nAtlantic") +  
#   # scale_x_continuous(limits = c(-0.35,0.38)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
#   geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# meso.area.pce3.a.partial.plot





###################################################################################
# FIGURES: BIVARIATE PLOTS (PACIFIC)                                              #
###################################################################################


# # PACIFIC: EELGRASS GROWTH FORM (PCz1) vs LATITUDE/CLIMATE (-PCe1) (Partial, site level)
# # Refit best model 
# pcz1.test.p.1 <- lm(rPC1.zos.pac ~ 
#                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
#                     , data = ZEN_2014_site_means_Pacific)
# # Extract residuals from model
# pcz1.pce1.p.partial <- partialResid(rPC1.zos.pac ~ rPC1.env.global.pac, pcz1.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# pcz1.pce1.p.partial <- cbind(pcz1.pce1.p.partial, Site)
# pcz1.pce1.p.partial.plot <- ggplot(pcz1.pce1.p.partial, aes(x = -xresid, y = yresid, col = "blue")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Latitude/Climate (-PCe1) | Z\nPacific)") +  
#   ylab("Eelgrass form (PCz1) | Z\nPacific") +  
#   # scale_x_continuous(limits = c(-.2, .24)) +
#   # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) # +
#   # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# pcz1.pce1.p.partial.plot


# PACIFIC: EELGRASS GROWTH FORM (PCz1) vs GENETIC FC2 (Partial, site level)
# Extract residuals from model
pcz1.fc2.p.partial <- partialResid(rPC1.zos.pac ~ rFC2.global.pac, pcz1.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
pcz1.fc2.p.partial <- cbind(pcz1.fc2.p.partial, Site)
pcz1.fc2.p.partial.plot <- ggplot(pcz1.fc2.p.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass genetics (FCA2) | Z\nPacific") +  
  ylab("Eelgrass form (PCz1) | Z\nPacific") +  
  scale_x_continuous(limits = c(-0.48,0.35)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz1.fc2.p.partial.plot
ggsave("figures/Fig2c_FCA2_PCz1_pacific.png", width = 5.5, height = 5)


# # PACIFIC: EELGRASS BIOMASS (inverse PCz2) vs (inverse) ENVIRONMENT PCe1 (Partial, site level)
# # Best model at site level
# pcz2.test.p.1 <- lm(rPC2.zos.pac ~ 
#                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
#                     , data = ZEN_2014_site_means_Pacific)
# # Extract residuals from model
# pcz2.pce1.p.partial <- partialResid(rPC2.zos.pac ~ rPC1.env.global.pac, pcz2.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# pcz2.pce1.p.partial <- cbind(pcz2.pce1.p.partial, Site)
# pcz2.pce1.p.partial.plot <- ggplot(pcz2.pce1.p.partial, aes(x = -xresid, y = -yresid, col = "forestgreen")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Latitude/climate (-PCe1) | Z\nPacific") +  
#   ylab("Eelgrass biomass (-PCz2) | Z\nPacific") +  
#   # scale_x_continuous(limits = c(-2.3,2.1)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) # +
# # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# pcz2.pce1.p.partial.plot


# # PACIFIC: EELGRASS BIOMASS (inverse PCz2) vs ENVIRONMENT PCe3 (Partial, site level)
# # Best model at site level
# pcz2.test.p.1 <- lm(rPC2.zos.pac ~ 
#                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
#                     , data = ZEN_2014_site_means_Pacific)
# # Extract residuals from model
# pcz2.pce3.p.partial <- partialResid(rPC2.zos.pac ~ rPC3.env.global.pac, pcz2.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# pcz2.pce3.p.partial <- cbind(pcz2.pce3.p.partial, Site)
# pcz2.pce3.p.partial.plot <- ggplot(pcz2.pce3.p.partial, aes(x = xresid, y = -yresid, col = "forestgreen")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Estuarine conditions (PCe3) | Z\nPacific") +  
#   ylab("Eelgrass biomass (-PCz2) | Z\nPacific") +  
#   # scale_x_continuous(limits = c(-2.3,2.1)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) # +
# # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# pcz2.pce3.p.partial.plot


# PACIFIC: EELGRASS BIOMASS (inverse PCz2) vs GENETIC FC2 (Partial, site level)
# Extract residuals from model
pcz2.fc2.p.partial <- partialResid(rPC2.zos.pac ~ rFC2.global.pac, pcz2.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
pcz2.fc2.p.partial <- cbind(pcz2.fc2.p.partial, Site)
pcz2.fc2.p.partial.plot <- ggplot(pcz2.fc2.p.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass genetics (FCA2) | Z\nPacific") +  
  ylab("Eelgrass biomass (-PCz2) | Z\nPacific") +  
  # scale_x_continuous(limits = c(30,75)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2.fc2.p.partial.plot
ggsave("figures/Fig2d_FCA2_PCz2_pacific.png", width = 5.5, height = 5)


# # PACIFIC: PERIPHYTON PER G EELGRASS vs EELGRASS FORM PCz1 (Partial, site level)
# # First fit model at site level
# ZEN_2014_site_means_Pacific$rPCz1.rFC1 <- ZEN_2014_site_means_Pacific$rPC1.zos.pac * ZEN_2014_site_means_Pacific$rFC1.global.pac
# peri.test.p.11.x <- lm(rperiphyton.perg.pac ~ 
#                        + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
#                      + rPC1.zos.pac + rPC2.zos.pac
#                      + rPCz1.rFC1
#                      , data = ZEN_2014_site_means_Pacific)# Extract residuals from model
# peri.pcz1.p.partial <- partialResid(rperiphyton.perg.pac ~ rPC1.zos.pac, peri.test.p.11.x, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# peri.pcz1.p.partial <- cbind(peri.pcz1.p.partial, Site)
# peri.pcz1.p.partial.plot <- ggplot(peri.pcz1.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Eelgrass form (PCz1) | Z\nPacific") +  
#   ylab("Periphyton mass | Z\nPacific") +  
#   # scale_x_continuous(limits = c(-3,4)) +
#   # scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
#   geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# peri.pcz1.p.partial.plot


# PACIFIC: PERIPHYTON PER AREA vs EELGRASS FORM PCz1 (Partial, site level)
# First fit model at site level
peri.area.test.p.1 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                         + rPC1.zos.pac + rPC2.zos.pac
                         , data = ZEN_2014_site_means_Pacific)
peri.area.pcz1.p.partial <- partialResid(rperiphyton.area.pac ~ rPC1.zos.pac, peri.area.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
peri.area.pcz1.p.partial <- cbind(peri.area.pcz1.p.partial, Site)
peri.area.pcz1.p.partial.plot <- ggplot(peri.area.pcz1.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass form (PCz1) | Z\nPacific") +  
  ylab("Periphyton mass | Z\nPacific") +  
  # scale_x_continuous(limits = c(-3,4)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.area.pcz1.p.partial.plot
ggsave("figures/Fig3e_PCz1_periphyton_pacific.png", width = 5.5, height = 5)

# PACIFIC: PERIPHYTON PER AREA vs EELGRASS BIOMASS PCz2 (Partial, site level)
peri.area.pcz2.p.partial <- partialResid(rperiphyton.area.pac ~ rPC2.zos.pac, peri.area.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
peri.area.pcz2.p.partial <- cbind(peri.area.pcz2.p.partial, Site)
peri.area.pcz2.p.partial.plot <- ggplot(peri.area.pcz2.p.partial, aes(x = xresid, y = -yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass biomass (-PCz2) | Z\nPacific") +  
  ylab("Periphyton mass | Z\nPacific") +  
  # scale_x_continuous(limits = c(-3,4)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
# geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.area.pcz2.p.partial.plot
ggsave("figures/Fig3f_PCz2_periphyton_pacific.png", width = 5.5, height = 5)


# # PACIFIC: PERIPHYTON PER AREA vs NUTRIENT STATUS (PCe2) (Partial, site level)
# # First fit model at site level
# peri.area.test.p.1 <- lm(rperiphyton.area.pac ~ 
#                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
#                          + rPC1.zos.pac + rPC2.zos.pac
#                          , data = ZEN_2014_site_means_Pacific)
# peri.area.pce2.p.partial <- partialResid(rperiphyton.area.pac ~ rPC2.env.global.pac, peri.area.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# peri.area.pce2.p.partial <- cbind(peri.area.pce2.p.partial, Site)
# peri.area.pce2.p.partial.plot <- ggplot(peri.area.pce2.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Nutrient status (PCe2) | Z\nPacific") +  
#   ylab("Periphyton per bottom area | Z\nPacific") +  
#   # scale_x_continuous(limits = c(-3,4)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) # +
#   # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# peri.area.pce2.p.partial.plot


# 
# # PACIFIC: PERIPHYTON PER G EELGRASS vs EELGRASS BIOMASS (inverse) PCz2 (Partial, site level)
# # First fit model at site level
# peri.perg.pcz2.p.partial <- partialResid(rperiphyton.perg.pac ~ rPC2.zos.pac, peri.test.p.11.x, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# peri.perg.pcz2.p.partial <- cbind(peri.perg.pcz2.p.partial, Site)
# peri.perg.pcz2.p.partial.plot <- ggplot(peri.perg.pcz2.p.partial, aes(x = xresid, y = -yresid, col = "forestgreen")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Eelgrass biomass (PCz2) | Z\nPacific") +  
#   ylab("Periphyton mass | Z\nPacific") +  
#   # scale_x_continuous(limits = c(-3,4)) +
#   # scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) # +
#   # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# peri.perg.pcz2.p.partial.plot


# # PACIFIC: PERIPHYTON vs LATITUDE/CLIMATE PCe1 (Partial, site level)
# peri.pce1.partial <- partialResid(rperiphyton.perg.pac ~ rPC1.env.global.pac, peri.test.p.11.x, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# peri.pce1.partial <- cbind(peri.pce1.partial, Site)
# peri.pce1.partial.plot <- ggplot(peri.pce1.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Ltitude/Climate (PCe1) | Z") +  
#   ylab("Periphyton mass | Z") +  
#   # scale_x_continuous(limits = c(-3,4)) +
#   # scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
#   geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# peri.pce1.partial.plot

# 
# # PACIFIC: PERIPHYTON vs EELGRASS GENETICS FCA2 (Partial, site level)
# peri.fca2.p.partial <- partialResid(rperiphyton.perg.pac ~ rFC2.global.pac, peri.test.p.11.x, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# peri.fca2.p.partial <- cbind(peri.fca2.p.partial, Site)
# peri.fca2.p.partial.plot <- ggplot(peri.fca2.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Eelgrass genetics (FCA2) | Z\nPacific") +  
#   ylab("Periphyton mass | Z") +  
#   # scale_x_continuous(limits = c(-3,4)) +
#   # scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
#   geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# peri.fca2.p.partial.plot


# # PACIFIC: MESOGRAZER BIOMASS PER G PLANT vs ENVIRONMENT PCe3 (Partial, site level)
# # Best model at site level
# meso.test.p.1 <- lm(rmesograzer.mass.perg.pac ~ 
#                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
#                     + rPC1.zos.pac + rPC2.zos.pac  
#                     , data = ZEN_2014_site_means_Pacific)
# # Extract residuals from model
# meso.pce3.p.partial <- partialResid(rmesograzer.mass.perg.pac ~ rPC3.env.global.pac, meso.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# meso.pce3.p.partial <- cbind(meso.pce3.p.partial, Site)
# meso.pce3.p.partial.plot <- ggplot(meso.pce3.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Estuarine conditions (PCe3) | Z\nPacific") +  
#   ylab("Mesograzer biomass | Z\nPacific") +  
#   # scale_x_continuous(limits = c(-1.4,1.35)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
#   geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# meso.pce3.p.partial.plot


# # PACIFIC: MESOGRAZER BIOMASS PER G PLANT vs GENETIC FC2 (Partial, site level)
# # Extract residuals from model
# meso.fc2.p.partial <- partialResid(rmesograzer.mass.perg.pac ~ rFC2.global.pac, meso.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add Ocean and site to residuals data frame 
# Ocean <- ZEN_2014_site_means_Pacific$Ocean
# Site <- ZEN_2014_site_means_Pacific$Site
# meso.fc2.p.partial <- cbind(meso.fc2.p.partial, Ocean, Site)
# meso.fc2.p.partial.plot <- ggplot(meso.fc2.p.partial, aes(x = xresid, y = yresid, col = "blue")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Genetic FC2 | Z\nPacific") +  
#   ylab("Mesograzer biomass | Z\nPacific") +  
#   # scale_x_continuous(limits = c(30,75)) +
#   # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) # +
#   # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# meso.fc2.p.partial.plot


# # PACIFIC: MESOGRAZER BIOMASS PER G PLANT vs MEADOW FORM (PCz1) (Partial, site level)
# # Extract residuals from model
# meso.pcz1.p.partial <- partialResid(rmesograzer.mass.perg.pac ~ rPC1.zos.pac, meso.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# meso.pcz1.p.partial <- cbind(meso.pcz1.p.partial, Site)
# meso.pcz1.p.partial.plot <- ggplot(meso.pcz1.p.partial, aes(x = xresid, y = yresid, col = "blue")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Eelgrass form (PCz1) | Z\nPacific") +  
#   ylab("Invertebrate biomass | Z\nPacific") +  
#   # scale_x_continuous(limits = c(30,75)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
# geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T)
# meso.pcz1.p.partial.plot



# # PACIFIC: MESOGRAZER BIOMASS PER G PLANT vs (inverse) EELGRASS BIOMASS (PCz2) (Partial, site level)
# # Extract residuals from model
# meso.perg.pcz2.p.partial <- partialResid(rmesograzer.mass.perg.pac ~ rPC2.zos.pac, meso.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# meso.perg.pcz2.p.partial <- cbind(meso.perg.pcz2.p.partial, Site)
# meso.perg.pcz2.p.partial.plot <- ggplot(meso.perg.pcz2.p.partial, aes(x = xresid, y = -yresid, col = "blue")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Eelgrass biomass (PCz2) | Z\nPacific") +  
#   ylab("Invertebrate mass | Z\nPacific") +  
#   # scale_x_continuous(limits = c(30,75)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) # +
# # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# meso.perg.pcz2.p.partial.plot
# ggsave("figures/Fig3h_PCz2_invertebrates_pacific.png", width = 5.5, height = 5)


# PACIFIC: MESOGRAZER BIOMASS PER AREA vs MEADOW FORM (PCz1) (Partial, site level)
# Best model at site level
meso.area.test.p.1 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                         + rPC1.zos.pac + rPC2.zos.pac  
                         , data = ZEN_2014_site_means_Pacific)
# Extract residuals from model
meso.pcz1.p.partial <- partialResid(rmesograzer.mass.area.pac ~ rPC1.zos.pac, meso.area.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
meso.pcz1.p.partial <- cbind(meso.pcz1.p.partial, Site)
meso.pcz1.p.partial.plot <- ggplot(meso.pcz1.p.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass form (PCz1) | Z\nPacific") +  
  ylab("Invertebrate biomass | Z\nPacific") +  
  # scale_x_continuous(limits = c(30,75)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.pcz1.p.partial.plot
ggsave("figures/Fig3g_PCz1_invertebrates_pacific.png", width = 5.5, height = 5)


# PACIFIC: MESOGRAZER BIOMASS PER AREA vs EELGRASS BIOMASS (-PCz2) (Partial, site level)
# Extract residuals from model
meso.area.pcz2.p.partial <- partialResid(rmesograzer.mass.area.pac ~ rPC2.zos.pac, meso.area.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
meso.area.pcz2.p.partial <- cbind(meso.area.pcz2.p.partial, Site)
meso.area.pcz2.p.partial.plot <- ggplot(meso.area.pcz2.p.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass biomass (-PCz2) | Z\nPacific") +  
  ylab("Mesograzer biomass | Z\nPacific") +  
  # scale_x_continuous(limits = c(30,75)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
# geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.area.pcz2.p.partial.plot
ggsave("figures/Fig3h_PCz2_invertebrates_pacific.png", width = 5.5, height = 5)


# # PACIFIC: MESOGRAZER BIOMASS PER AREA vs ESTUARINE COINDITIONS (PCe3) (Partial, site level)
# # Best model at site level
# meso.area.test.p.1 <- lm(rmesograzer.mass.area.pac ~ 
#                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
#                          + rPC1.zos.pac + rPC2.zos.pac  
#                          , data = ZEN_2014_site_means_Pacific)
# # Extract residuals from model
# meso.area.pce3.p.partial <- partialResid(rmesograzer.mass.area.pac ~ rPC3.env.global.pac, meso.area.test.p.1, ZEN_2014_site_means_Pacific) 
# # Add site to residuals data frame 
# Site <- ZEN_2014_site_means_Pacific$Site
# meso.area.pce3.p.partial <- cbind(meso.area.pce3.p.partial, Site)
# meso.area.pce3.p.partial.plot <- ggplot(meso.area.pce3.p.partial, aes(x = xresid, y = yresid, col = "blue")) +
#   geom_point(size = 4) +
#   geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
#   scale_color_manual(values = c("forestgreen")) +
#   xlab("Estuarine conditions (PCe3) | Z\nPacific") +  
#   ylab("Mesograzer mass per area | Z\nPacific") +  
#   # scale_x_continuous(limits = c(30,75)) +
#   scale_y_continuous(position = "right") +
#   theme_bw(base_size = 12) +
#   theme(legend.position = "none", 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.5))
#   ) +
#   geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
# meso.area.pce3.p.partial.plot


