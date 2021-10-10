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
# FIGURES: BIVARIATE PLOTS (GLOBAL)                                               #
# FIGURES: BIVARIATE PLOTS (ATLANTIC)                                             #
# FIGURES: BIVARIATE PLOTS (PACIFIC)                                              #
# FIGURES: 3D PLOTS                                                               #
# BIVARIATE PLOTS (RAW DATA)                                                      #
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
ZEN_2014_site_means <- read.csv("ZEN_2014_site_means_20210315.csv",  header = TRUE)
ZEN_2014_site_means_49 <- droplevels(subset(ZEN_2014_site_means, Site != "SW.A"))
ZEN_2014_site_means_49_Atlantic <- read.csv("ZEN_2014_site_means_49_Atlantic_20210314.csv",  header = TRUE)
ZEN_2014_site_means_Pacific <- read.csv("ZEN_2014_site_means_Pacific_20210314.csv",  header = TRUE)

# Read in zen2014 PLOT-level data sets for modeling (MINUS SW.A), with missing data imputed:
ZEN_2014_plot_49 <- read.csv("ZEN_2014_plot_49_noNA_20210227.csv", header = TRUE)
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

# Correlates of environmental/oceanographic PC axes 
pairs.panels(ZEN_2014_site_means_renamed[,c("Latitude", "Env PCe1", "Env PCe2", 
  "Env PCe3", "SST mean", "SST range", "Salinity", "PAR", "Cloud cover",
  "PO4 (log)", "Leaf %N (log)", "chl (log)")], 
  hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 14, 
  bg = c("blue","green")[ZEN_2014_site_means_renamed$Ocean])

# Explore correlations among eelgrass characteristics by ocean 
pairs.panels(ZEN_2014_site_means_renamed[,c("Eelgrass form (PCz1)", "Eelgrass biomass (PCz2)",
  "Canopy ht (log)", "Shoot density (log)", "AG mass (log)", "BG mass (log)", "Sheath L (log)", 
  "Sheath W (log)")], hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 3, 
  bg = c("blue","green")[ZEN_2014_site_means_renamed$Ocean])



# Correlations between environment and biology
pairs.panels(ZEN_2014_site_means_renamed[,c("Ocean", "Latitude", "Env PCe1", "Env PCe2", 
  "Env PCe3", "Genetics FCA1", "Genetics FCA2", "Eelgrass form PCz1", "Eelgrass form PCz2",
  "Periphyton mass (log)", "Mesograzer mass (log)", "Periphyton mass/area (log)", 
  "Mesograzer mass/area (log)")], hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 12, 
  bg = c("blue","green")[ZEN_2014_site_means_renamed$Ocean])

# Atlantic: Correlations among covariates
pairs.panels(ZEN_2014_site_means_renamed_Atlantic[,c("Latitude", "Env PCe1", "Env PCe2", "Env PCe3", 
  "Genetics FCA1", "Genetics FCA2", "Eelgrass form PCz1", "Eelgrass form PCz2",
  "Periphyton mass (log)", "Mesograzer mass (log)")], hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 5, bg = c("blue")
  )

# PACIFIC: Correlations among covariates
pairs.panels(ZEN_2014_site_means_renamed_Pacific[,c("Latitude", "Env PCe1", "Env PCe2", "Env PCe3", 
  "Genetics FCA1", "Genetics FCA2", "Eelgrass form PCz1", "Eelgrass form PCz2",
  "Periphyton mass (log)", "Mesograzer mass (log)")], hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 5, bg = c("forestgreen")
)



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
# ggsave(pc1z.pc2z, filename = "./pc1z.pc2z.png",  width = 3, height = 10, bg = "transparent")
pdf("pc1z.pc2z.pdf",width = 6.5, height = 6); pc1z.pc2z;  dev.off() 


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
# ggsave(fc1.fc2, filename = "./fc1.fc2.png",  width = 6.5, height = 6, bg = "transparent")
pdf("fc1.fc2.pdf",width = 6.5, height = 6); fc1.fc2;  dev.off() 


#  Environment PC1 vs PC2
pce1.pce2 <- ggplot(ZEN_2014_site_means, aes(x = PC1.env.global, y = PC2.env.global, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
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
pdf("pce1.pce2.pdf",width = 6.5, height = 6); pce1.pce2;  dev.off() 


#  Environment PC2 vs PC3
pce2.pce3 <- ggplot(ZEN_2014_site_means, aes(x = PC2.env.global, y = PC3.env.global, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
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
pdf("pce2.pce3.pdf",width = 6.5, height = 6); pce2.pce3;  dev.off() 


# Eelgrass Morphology PC1: Effect of FC1 (RAW) 
pc1z.fc1 = ggplot(ZEN_2014_site_means, aes(x = FC1, y = PC1.zos.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  # geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass genetics (FCA1)") +  
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
pc1z.fc1
# ggsave(pc1z.fc1, filename = "./pc1z.fc1.png",  width = 3, height = 10, bg = "transparent")
pdf("pc1z.fc1.pdf",width = 4.5, height = 4); pc1z.fc1;  dev.off() 


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
# ggsave(pc1z.fc2, filename = "./pc1z.fc2",  width = 3, height = 10, bg = "transparent")
pdf("pc1z.fc2.pdf",width = 4.5, height = 4); pc1z.fc2;  dev.off() 


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
pdf("shoots.canopy.pdf",width = 6.5, height = 6); shoots.canopy;  dev.off() 
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
pdf("shoots.zag.pdf",width = 6.5, height = 6); shoots.zag;  dev.off() 



###################################################################################
# FIGURE: MODEL COEFFICIENTS BY RESPONSE AND OCEAN                                #
###################################################################################

# ATLANTIC 

# Read in GLM model coefficients produced by modeling script (ZEN_2014_model_comparison_20210227.R)
coeffs.atl <- read.csv("zen_glm_coefficients_atlantic.csv", header = TRUE)
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
pdf("pcz1.atl.coeffs.a.pdf", width = 6, height = 4); pcz1.atl.coeffs.a; dev.off() 


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
pdf("pcz2.atl.coeffs.a.pdf", width = 6, height = 4); pcz2.atl.coeffs.a; dev.off() 


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
pdf("peri.atl.coeffs.a.pdf", width = 6, height = 4); peri.atl.coeffs.a; dev.off() 


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
pdf("meso.atl.coeffs.a.pdf", width = 6, height = 4); meso.atl.coeffs.a; dev.off() 


# PACIFIC 

# Read in GLM model coefficients
coeffs.pac <- read.csv("zen_glm_coefficients_pacific.csv", header = TRUE)
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
pdf("pcz1.pac.coeffs.p.pdf", width = 6, height = 4); pcz1.pac.coeffs.p; dev.off() 


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
pdf("pcz2.pac.coeffs.p.pdf", width = 6, height = 4); pcz2.pac.coeffs.p; dev.off() 


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
pdf("peri.pac.coeffs.p.pdf", width = 6, height = 4); peri.pac.coeffs.p; dev.off() 


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
pdf("meso.pac.coeffs.p.pdf", width = 6, height = 4); meso.pac.coeffs.p; dev.off() 


###################################################################################
# FIGURE: COMPARATIVE EFFECTS OF ENVIRONMENT AND EVO HISTORY (BAR CHARTS)         #
###################################################################################

# To compare influences of environment and evolutionary history,calculate the weighted 
# SUM of standardized effects of each (e.g., PCe1, PCe2, PCe3 for environment effect). 

# NOTE: I have had residual stuff in memory foul up the plotting and throw an error. 
# Therefore it's a good idea to start with:
dev.off()

library(dplyr)          # for data manipulation
library(tidyr)          # for data manipulation
library(magrittr)       # for easier syntax in one or two areas
library(gridExtra)      # for generating some comparison plots
library(ggplot2)        # for generating the visualizations


# Read in standardized coefficients from best models: ATLANTIC
zen_best_models_atlantic <- read.csv("zen_2014_model_coeffs_atlantic_20210410.csv",header = TRUE)
names(zen_best_models_atlantic)

zen_best_models_atlantic$response <- as.factor(zen_best_models_atlantic$response)
zen_best_models_atlantic$predictor <- as.factor(zen_best_models_atlantic$predictor)
zen_best_models_atlantic$effect.type <- as.factor(zen_best_models_atlantic$effect.type)
levels(zen_best_models_atlantic$response)
levels(zen_best_models_atlantic$predictor)
levels(zen_best_models_atlantic$effect.type)

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
pdf("summed_paths_atl_high_scale.pdf", width = 3, height = 10); summed_paths_atl_high_scale; dev.off() 


# Read in standardized coefficients from best models: PACIFIC
zen_best_models_pacific <- read.csv("zen_2014_model_coeffs_pacific_20210410.csv",header = TRUE)
names(zen_best_models_pacific)

zen_best_models_pacific$response <- as.factor(zen_best_models_pacific$response)
zen_best_models_pacific$predictor <- as.factor(zen_best_models_pacific$predictor)
levels(zen_best_models_pacific$response)
levels(zen_best_models_pacific$predictor)

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
pdf("summed_paths_pac_high_scale.pdf", width = 3, height = 10); summed_paths_pac_high_scale; dev.off() 



###################################################################################
# FIGURE: ZEN SITE MAPS                                                           #
###################################################################################

# From: ZEN 2014 short MS analysis 20201125.R

# For mapping
# library(sf)
# library(sp)
# library(raster)
# library(spData)
# library(tmap)      # for static and interactive maps
# library(dplyr)     # for mutate
# library(leaflet) # for interactive maps
# library(mapview) # for interactive maps
# library(shinyjs) # for something or other ...

# POLAR PROJECTION USING GGOCEANMAPS

library(ggOceanMaps) # For base map layer
library(ggrepel) # For adding labels to points on map

zen <- read.csv("zen_2014_sites.csv")
zen.map.data <- subset(ZEN_2014_site_means, select = c(PC1.zos.site, Zostera.AG.mass.site, Zostera.BG.mass.site))
zen <- cbind(zen, zen.map.data)

# ZEN 2014 site map, colored by ocean, with site labels 
polarmap.sites <- basemap(data = zen, land.col = "gray", land.border.col = "black") +
  geom_spatial_point(data = zen, aes(x = longitude, y = latitude, group = ocean, fill = ocean), 
    size = 4, shape = 21) +
  geom_text_repel(data=transform_coord(zen), aes(x=longitude, y=latitude, label=zen$site), 
    hjust = 2, color="black", size = 3) +
  scale_fill_manual("Ocean", values = c("blue", "forestgreen")) +
  labs(x= "Longitude", y = "Latitude") +
  theme(axis.text = element_text(colour = "black"))
polarmap.sites
pdf("polarmap.sites.pdf", width = 7, height = 7); polarmap.sites; dev.off()


# ZEN 2014 site map, colored by PCz1
polarmap.pcz1 <- basemap(data = zen, land.col = "darkolivegreen3", land.border.col = "black") +
  geom_spatial_point(data = zen, aes(x = longitude, y = latitude, fill = PC1.zos.site),
    size = 3.5, shape = 21) +
  # geom_text_repel(data=transform_coord(zen), aes(x=longitude, y=latitude, label=zen$site),
  #                 hjust = 2, color="black", size = 3) +
  scale_fill_gradient("PC1.zos.site", low = "yellow", high = "red4") +
  labs(x= "Longitude", y = "Latitude") +
  theme(axis.text = element_text(colour = "black"))
polarmap.pcz1
pdf("polarmap.pcz1.pdf", width = 7, height = 7); polarmap.pcz1; dev.off()


# ZEN 2014 site map, colored by above-ground biomass
polarmap.zag <- basemap(data = zen, land.col = "darkolivegreen3", land.border.col = "black") +
  geom_spatial_point(data = zen, aes(x = longitude, y = latitude, fill = Zostera.AG.mass.site),
                     size = 3.5, shape = 21) +
  # geom_text_repel(data=transform_coord(zen), aes(x=longitude, y=latitude, label=zen$site),
  #                 hjust = 2, color="black", size = 3) +
  scale_fill_gradient("Zostera.AG.mass.site", low = "yellow", high = "red4") +
  labs(x= "Longitude", y = "Latitude") +
  theme(axis.text = element_text(colour = "black"))
polarmap.zag
pdf("polarmap.zag.pdf", width = 7, height = 7); polarmap.zag; dev.off()


###################################################################################
# FIGURES: BIVARIATE PLOTS (GLOBAL)                                               #
###################################################################################

# GLOBAL: EELGRASS GROWTH FORM (PCz1) vs ENVIRONMENT PCe1 (Partial, site level)
# Best model at site level
pcz1.site.g.4 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rPC3.env.global 
                    , data = ZEN_2014_site_means_49)
# Extract residuals from model
pcz1.pce1.partial <- partialResid(rPC1.zos.site ~ rPC1.env.global, pcz1.site.g.4, ZEN_2014_site_means_49) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49$Site
pcz1.pce1.partial <- cbind(pcz1.pce1.partial, Ocean, Site)
pcz1.pce1.partial.plot <- ggplot(pcz1.pce1.partial, aes(x = (-xresid), y = yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Latitude/Climate (-PCe1) | Z\nGlobal") +  
  ylab("Meadow form (PCz1) | Z\nGlobal") +  
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
pcz1.pce1.partial.plot
pdf("pcz1.pce1.partial.plot.pdf", width = 4.5, height = 4); pcz1.pce1.partial.plot;  dev.off() 


# GLOBAL: EELGRASS GROWTH FORM (PCz1) vs GENETIC FC2 (Partial, site level)
# Extract residuals from model
pcz1.fc2.partial <- partialResid(rPC1.zos.site ~ rFC2, pcz1.site.g.4, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
pcz1.fc2.partial <- cbind(pcz1.fc2.partial, Ocean, Site)
pcz1.fc2.partial.plot <- ggplot(pcz1.fc2.partial, aes(x = xresid, y = yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Genetic FCA2 | Z\nGlobal") +  
  ylab("Meadow form (PCz1) | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz1.fc2.partial.plot
pdf("pcz1.fc2.partial.plot.pdf",width = 4.5, height = 4); pcz1.fc2.partial.plot;  dev.off() 


# GLOBAL: EELGRASS BIOMASS (inverse PCz2) vs ENVIRONMENT PCe3 (Partial, site level)
# Best model at site level
pcz2.site.g.4 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rPC3.env.global # added
                    , data = ZEN_2014_site_means_49)
# Extract residuals from model
pcz2.pce3.partial <- partialResid(rPC2.zos.site ~ rPC3.env.global, pcz2.site.g.4, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
pcz2.pce3.partial <- cbind(pcz2.pce3.partial, Ocean, Site)
pcz2.pce3.partial.plot <- ggplot(pcz2.pce3.partial, aes(x = xresid, y = -yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Estuarine conditions (PCe3) | Z\nGlobal") +  
  ylab("Eelgrass biomass (-PCz2) | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-1.1,1.3)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2.pce3.partial.plot
pdf("pcz2.pce3.partial.plot.pdf",width = 4.5, height = 4); pcz2.pce3.partial.plot;  dev.off() 


# GLOBAL: EELGRASS BIOMASS (inverse PCz2) vs EELGRASS GENETICS FC2 (Partial, site level)
# Extract residuals from model
pcz2.fc2.partial <- partialResid(rPC2.zos.site ~ rFC2, pcz2.site.g.4, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
pcz2.fc2.partial <- cbind(pcz2.fc2.partial, Ocean, Site)
pcz2.fc2.partial.plot <- ggplot(pcz2.fc2.partial, aes(x = xresid, y = (-yresid), group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Genetic FCA2 | Z\nGlobal") +  
  ylab("Eelgrass biomass (-PCz2) | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-165,300)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2.fc2.partial.plot
pdf("pcz2.fc2.partial.plot.pdf",width = 4.5, height = 4); pcz2.fc2.partial.plot;  dev.off() 


# GLOBAL: PERIPHYTON vs EELGRASS GENETICS FC2 (Partial, site level)
# Best model at site level
peri.site.g.2 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + Ocean*rPC1.env.global # added
                    , data = ZEN_2014_site_means_49)
# Extract residuals from model
peri.fc2.partial <- partialResid(rperiphyton.perg ~ rFC2, peri.site.g.2, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
peri.fc2.partial <- cbind(peri.fc2.partial, Ocean, Site)
peri.fc2.partial.plot <- ggplot(peri.fc2.partial, aes(x = xresid, y = (-yresid), group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Genetic FCA2 | Z\nGlobal") +  
  ylab("Periphyton per g eelgrass | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.fc2.partial.plot
pdf("peri.fc2.partial.plot.pdf",width = 4.5, height = 4); peri.fc2.partial.plot;  dev.off() 


# GLOBAL: PERIPHYTON vs ENVIRONMENT PCe1 (Partial, site level)
# Extract residuals from model
peri.pce1.partial <- partialResid(rperiphyton.perg ~ rPC1.env.global, peri.site.g.2, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
peri.pce1.partial <- cbind(peri.pce1.partial, Ocean, Site)
peri.pce1.partial.plot <- ggplot(peri.pce1.partial, aes(x = xresid, y = (-yresid), group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Latitude/climate (PCe1) | Z\nGlobal") +  
  ylab("Periphyton per g eelgrass | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.pce1.partial.plot
pdf("peri.pce1.partial.plot.pdf",width = 4.5, height = 4); peri.pce1.partial.plot;  dev.off() 


# GLOBAL: PERIPHYTON PER G EELGRASS vs EELGRASS FORM PCz1 (Partial, site level)
# Best model at site level
peri.perg.site.g.2 <- lm(rperiphyton.perg ~ Ocean
  + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
  + rPC1.zos.site + rPC2.zos.site
  + Ocean*rPC1.env.global # added
  , data = ZEN_2014_site_means_49)
# Extract residuals from model
peri.perg.pcz1.partial <- partialResid(rperiphyton.perg ~ rPC1.zos.site, peri.perg.site.g.2, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
peri.perg.pcz1.partial <- cbind(peri.perg.pcz1.partial, Ocean, Site)
peri.perg.pcz1.partial.plot <- ggplot(peri.perg.pcz1.partial, aes(x = xresid, y = yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1) | Z\nGlobal") +  
  ylab("Periphyton mass | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.perg.pcz1.partial.plot
pdf("peri.perg.pcz1.partial.plot.pdf",width = 4.5, height = 4); peri.perg.pcz1.partial.plot;  dev.off() 


# GLOBAL: PERIPHYTON PER AREA vs EELGRASS FORM PCz1 (Partial, site level)
# Best model at site level
peri.area.site.g.2 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + Ocean*rPC1.env.global 
                         , data = ZEN_2014_site_means_49)
# Extract residuals from model
peri.area.pcz1.partial <- partialResid(rperiphyton ~ rPC1.zos.site, peri.area.site.g.2, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
peri.area.pcz1.partial <- cbind(peri.area.pcz1.partial, Ocean, Site)
peri.area.pcz1.partial.plot <- ggplot(peri.area.pcz1.partial, aes(x = xresid, y = yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1) | Z\nGlobal") +  
  ylab("Periphyton per area | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.area.pcz1.partial.plot
pdf("peri.area.pcz1.partial.plot.pdf",width = 4.5, height = 4); peri.area.pcz1.partial.plot;  dev.off() 


# GLOBAL: PERIPHYTON PER AREA vs EELGRASS GENETICS FC2 (Partial, site level)
# Best model at site level
# Extract residuals from model
peri.fc2.partial <- partialResid(rperiphyton ~ rFC2, peri.area.site.g.2, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
peri.fc2.partial <- cbind(peri.fc2.partial, Ocean, Site)
peri.fc2.partial.plot <- ggplot(peri.fc2.partial, aes(x = xresid, y = (yresid), group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Genetic FCA2 | Z\nGlobal") +  
  ylab("Periphyton per area | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.fc2.partial.plot
pdf("peri.fc2.partial.plot.pdf",width = 4.5, height = 4); peri.fc2.partial.plot;  dev.off() 


# GLOBAL: MESOGRAZER BIOMASS  vs ENVIRONMENT PCe3 (Partial, site level)
# Best model at site level
meso.site.g.16 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.rFC1 # added
                     , data = ZEN_2014_site_means_49)
# Extract residuals from model
meso.pce3.partial <- partialResid(rmesograzer.mass.perg ~ rPC3.env.global, meso.site.g.16, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
meso.pce3.partial <- cbind(meso.pce3.partial, Ocean, Site)
meso.pce3.partial.plot <- ggplot(meso.pce3.partial, aes(x = xresid, y = yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Estuarine conditions (PCe3) | Z\nGlobal") +  
  ylab("Mesograzer biomass | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-1.1,1.3)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.pce3.partial.plot
pdf("meso.pce3.partial.plot.pdf",width = 4.5, height = 4); meso.pce3.partial.plot;  dev.off() 


# GLOBAL: MESOGRAZER BIOMASS PER G EELGRASS vs MEADOW FORM Pz1 (Partial, site level)
# Best model at site level
meso.site.g.16 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rFC1 # added
                     , data = ZEN_2014_site_means_49)
# Extract residuals from model
meso.pcz1.partial <- partialResid(rmesograzer.mass.perg ~ rPC1.zos.site, meso.site.g.16, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
meso.pcz1.partial <- cbind(meso.pcz1.partial, Ocean, Site)
meso.pcz1.partial.plot <- ggplot(meso.pcz1.partial, aes(x = xresid, y = yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1) | Z\nGlobal") +  
  ylab("Mesograzer biomass | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-1.1,1.3)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.pcz1.partial.plot
pdf("meso.pcz1.partial.plot.pdf",width = 4.5, height = 4); meso.pcz1.partial.plot;  dev.off() 


# GLOBAL: MESOGRAZER BIOMASS PER AREA vs EELGRASS FORM PCz2 (RAW, site level)
meso.area.pcz2.plot <- ggplot(ZEN_2014_site_means_49, aes(x = -PC2.zos.site, y = log10.mesograzer.mass.per.area.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass biomass (-PCz2)\nGlobal") +  
  ylab("Mesograzer mass per area\nGlobal") +  
  # scale_x_continuous(limits = c(-1.1,1.3)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.area.pcz2.plot
pdf("meso.area.pcz2.plot.pdf",width = 4.5, height = 4); meso.area.pcz2.plot;  dev.off() 


# GLOBAL: MESOGRAZER BIOMASS PER G EELGRASS vs EELGRASS FORM PCz2 (Partial, site level)
meso.perg.pcz2.partial <- partialResid(rmesograzer.mass.perg ~ rPC2.zos.site, meso.site.g.16, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
meso.perg.pcz2.partial <- cbind(meso.perg.pcz2.partial, Ocean, Site)
meso.perg.pcz2.partial.plot <- ggplot(meso.perg.pcz2.partial, aes(x = xresid, y = -yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass biomass (PCz2) | Z\nGlobal") +  
  ylab("Mesograzer mass | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-1.1,1.3)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.perg.pcz2.partial.plot
pdf("meso.perg.pcz2.partial.plot.pdf",width = 4.5, height = 4); meso.perg.pcz2.partial.plot;  dev.off() 


# GLOBAL: MESOGRAZER BIOMASS PER AREA vs EELGRASS FORM PCz2 (Partial, site level)
# Best model at site level
meso.area.site.g.16 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rFC1 # added
                          , data = ZEN_2014_site_means_49)# Extract residuals from model
meso.area.pcz2.partial <- partialResid(rmesograzer.mass ~ rPC2.zos.site, meso.area.site.g.16, ZEN_2014_site_means_49) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_49$Ocean
Site <- ZEN_2014_site_means_49$Site
meso.area.pcz2.partial <- cbind(meso.area.pcz2.partial, Ocean, Site)
meso.area.pcz2.partial.plot <- ggplot(meso.area.pcz2.partial, aes(x = xresid, y = -yresid, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass biomass (PCz2) | Z\nGlobal") +  
  ylab("Mesograzer mass per area | Z\nGlobal") +  
  # scale_x_continuous(limits = c(-1.1,1.3)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.area.pcz2.partial.plot
pdf("meso.area.pcz2.partial.plot.pdf",width = 4.5, height = 4); meso.area.pcz2.partial.plot;  dev.off() 


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
pdf("pcz1.pce1.a.partial.plot.pdf",width = 4.5, height = 4); pcz1.pce1.a.partial.plot;  dev.off() 


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
pdf("pcz1.fc2.a.partial.plot.pdf",width = 4.5, height = 4); pcz1.fc2.a.partial.plot;  dev.off() 

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
pdf("pcz2.pce1.a.partial.plot.pdf",width = 4.5, height = 4); pcz2.pce1.a.partial.plot;  dev.off() 


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
pdf("pcz2.pce3.a.partial.plot.pdf",width = 4.5, height = 4); pcz2.pce3.a.partial.plot;  dev.off() 


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
pdf("pcz2.fc2.a.partial.plot.pdf",width = 4.5, height = 4); pcz2.fc2.a.partial.plot;  dev.off() 


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
pdf("peri.pcz1.a.partial.plot.pdf",width = 4.5, height = 4); peri.pcz1.a.partial.plot;  dev.off() 


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
pdf("peri.pcz1.a.partial.plot.pdf",width = 4.5, height = 4); peri.pcz1.a.partial.plot;  dev.off() 


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
pdf("peri.pcz2.a.partial.plot.pdf",width = 4.5, height = 4); peri.pcz2.a.partial.plot;  dev.off() 


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
pdf("peri.pce2.a.partial.plot.pdf",width = 4.5, height = 4); peri.pce2.a.partial.plot;  dev.off() 


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
pdf("peri.pcz2.a.partial.plot.pdf",width = 4.5, height = 4); peri.pcz2.a.partial.plot;  dev.off() 


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
pdf("meso.pce3.a.partial.plot.pdf",width = 4.5, height = 4); meso.pce3.a.partial.plot;  dev.off() 


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
pdf("meso.fc2.a.partial.plot.pdf",width = 4.5, height = 4); meso.fc2.a.partial.plot;  dev.off() 


# ATLANTIC: MESOGRAZER BIOMASS PER G PLANT vs MEADOW FORM (PCz1) (Partial, site level)
# Extract residuals from model
meso.pcz1.a.partial <- partialResid(rmesograzer.mass.perg.atl ~ rPC1.zos.atl, meso.test.a.9.x, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
meso.pcz1.a.partial <- cbind(meso.pcz1.a.partial, Site)
meso.pcz1.a.partial.plot <- ggplot(meso.pcz1.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass form (PCz1) | Z\nAtlantic") +  
  ylab("Mesograzer biomass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-0.25,0.2)) +
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
pdf("meso.pcz1.a.partial.plot.pdf",width = 4.5, height = 4); meso.pcz1.a.partial.plot;  dev.off() 


# ATLANTIC: MESOGRAZER BIOMASS PER G PLANT vs EELGRASS BIOMASS (-PCz2) (Partial, site level)
# Extract residuals from model
meso.perg.pcz2.a.partial <- partialResid(rmesograzer.mass.perg.atl ~ rPC2.zos.atl, meso.test.a.9.x, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
meso.perg.pcz2.a.partial <- cbind(meso.perg.pcz2.a.partial, Site)
meso.perg.pcz2.a.partial.plot <- ggplot(meso.perg.pcz2.a.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Eelgrass biomass (PCz2) | Z\nAtlantic") +  
  ylab("Mesograzer biomass | Z\nAtlantic") +  
  scale_x_continuous(limits = c(-0.35,0.38)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
  # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.perg.pcz2.a.partial.plot
pdf("meso.perg.pcz2.a.partial.plot.pdf",width = 4.5, height = 4); meso.perg.pcz2.a.partial.plot;  dev.off() 


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
  ylab("Mesograzer biomass | Z\nAtlantic") +  
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
pdf("meso.pcz1.a.partial.plot.pdf",width = 4.5, height = 4); meso.pcz1.a.partial.plot;  dev.off() 



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
  ylab("Mesograzer biomass | Z\nAtlantic") +  
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
pdf("meso.area.pcz2.a.partial.plot.pdf",width = 4.5, height = 4); meso.area.pcz2.a.partial.plot;  dev.off() 


# ATLANTIC: MESOGRAZER BIOMASS PER AREA vs ESTUARINE CONDITIONS (PCe3) (Partial, site level)
# Extract residuals from model
meso.area.pce3.a.partial <- partialResid(rmesograzer.mass.area.atl ~ rPC3.env.global.atl, meso.area.test.a.1, ZEN_2014_site_means_49_Atlantic) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_49_Atlantic$Site
meso.area.pce3.a.partial <- cbind(meso.area.pce3.a.partial, Site)
meso.area.pce3.a.partial.plot <- ggplot(meso.area.pce3.a.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Estuarine conditions (PCe3) | Z\nAtlantic") +  
  ylab("Mesograzer biomass | Z\nAtlantic") +  
  # scale_x_continuous(limits = c(-0.35,0.38)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.area.pce3.a.partial.plot
pdf("meso.area.pce3.a.partial.plot.pdf",width = 4.5, height = 4); meso.area.pce3.a.partial.plot;  dev.off() 


###################################################################################
# FIGURES: BIVARIATE PLOTS (PACIFIC)                                              #
###################################################################################


# PACIFIC: EELGRASS GROWTH FORM (PCz1) vs LATITUDE/CLIMATE (-PCe1) (Partial, site level)
# Refit best model 
pcz1.test.p.1 <- lm(rPC1.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
                    , data = ZEN_2014_site_means_Pacific)
# Extract residuals from model
pcz1.pce1.p.partial <- partialResid(rPC1.zos.pac ~ rPC1.env.global.pac, pcz1.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
pcz1.pce1.p.partial <- cbind(pcz1.pce1.p.partial, Site)
pcz1.pce1.p.partial.plot <- ggplot(pcz1.pce1.p.partial, aes(x = -xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Latitude/Climate (-PCe1) | Z\nPacific)") +  
  ylab("Eelgrass form (PCz1) | Z\nPacific") +  
  # scale_x_continuous(limits = c(-.2, .24)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
  # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz1.pce1.p.partial.plot
pdf("pcz1.pce1.p.partial.plot.pdf",width = 4.5, height = 4); pcz1.pce1.p.partial.plot;  dev.off() 


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
pdf("pcz1.fc2.p.partial.plot.pdf",width = 4.5, height = 4); pcz1.fc2.p.partial.plot;  dev.off() 


# PACIFIC: EELGRASS BIOMASS (inverse PCz2) vs (inverse) ENVIRONMENT PCe1 (Partial, site level)
# Best model at site level
pcz2.test.p.1 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
                    , data = ZEN_2014_site_means_Pacific)
# Extract residuals from model
pcz2.pce1.p.partial <- partialResid(rPC2.zos.pac ~ rPC1.env.global.pac, pcz2.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
pcz2.pce1.p.partial <- cbind(pcz2.pce1.p.partial, Site)
pcz2.pce1.p.partial.plot <- ggplot(pcz2.pce1.p.partial, aes(x = -xresid, y = -yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Latitude/climate (-PCe1) | Z\nPacific") +  
  ylab("Eelgrass biomass (-PCz2) | Z\nPacific") +  
  # scale_x_continuous(limits = c(-2.3,2.1)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
# geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2.pce1.p.partial.plot
pdf("pcz2.pce1.p.partial.plot.pdf",width = 4.5, height = 4); pcz2.pce1.p.partial.plot;  dev.off() 


# PACIFIC: EELGRASS BIOMASS (inverse PCz2) vs ENVIRONMENT PCe3 (Partial, site level)
# Best model at site level
pcz2.test.p.1 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
                    , data = ZEN_2014_site_means_Pacific)
# Extract residuals from model
pcz2.pce3.p.partial <- partialResid(rPC2.zos.pac ~ rPC3.env.global.pac, pcz2.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
pcz2.pce3.p.partial <- cbind(pcz2.pce3.p.partial, Site)
pcz2.pce3.p.partial.plot <- ggplot(pcz2.pce3.p.partial, aes(x = xresid, y = -yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Estuarine conditions (PCe3) | Z\nPacific") +  
  ylab("Eelgrass biomass (-PCz2) | Z\nPacific") +  
  # scale_x_continuous(limits = c(-2.3,2.1)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
# geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2.pce3.p.partial.plot
pdf("pcz2.pce3.p.partial.plot.pdf",width = 4.5, height = 4); pcz2.pce3.p.partial.plot;  dev.off() 


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
pdf("pcz2.fc2.p.partial.plot.pdf",width = 4.5, height = 4); pcz2.fc2.p.partial.plot;  dev.off() 


# PACIFIC: PERIPHYTON PER G EELGRASS vs EELGRASS FORM PCz1 (Partial, site level)
# First fit model at site level
ZEN_2014_site_means_Pacific$rPCz1.rFC1 <- ZEN_2014_site_means_Pacific$rPC1.zos.pac * ZEN_2014_site_means_Pacific$rFC1.global.pac
peri.test.p.11.x <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPCz1.rFC1
                     , data = ZEN_2014_site_means_Pacific)# Extract residuals from model
peri.pcz1.p.partial <- partialResid(rperiphyton.perg.pac ~ rPC1.zos.pac, peri.test.p.11.x, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
peri.pcz1.p.partial <- cbind(peri.pcz1.p.partial, Site)
peri.pcz1.p.partial.plot <- ggplot(peri.pcz1.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
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
peri.pcz1.p.partial.plot
pdf("peri.pcz1.p.partial.plot.pdf",width = 4.5, height = 4); peri.pcz1.p.partial.plot;  dev.off() 


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
pdf("peri.area.pcz1.p.partial.plot.pdf",width = 4.5, height = 4); peri.area.pcz1.p.partial.plot;  dev.off() 


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
pdf("peri.area.pcz2.p.partial.plot.pdf",width = 4.5, height = 4); peri.area.pcz2.p.partial.plot;  dev.off() 


# PACIFIC: PERIPHYTON PER AREA vs NUTRIENT STATUS (PCe2) (Partial, site level)
# First fit model at site level
peri.area.test.p.1 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                         + rPC1.zos.pac + rPC2.zos.pac
                         , data = ZEN_2014_site_means_Pacific)
peri.area.pce2.p.partial <- partialResid(rperiphyton.area.pac ~ rPC2.env.global.pac, peri.area.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
peri.area.pce2.p.partial <- cbind(peri.area.pce2.p.partial, Site)
peri.area.pce2.p.partial.plot <- ggplot(peri.area.pce2.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Nutrient status (PCe2) | Z\nPacific") +  
  ylab("Periphyton per bottom area | Z\nPacific") +  
  # scale_x_continuous(limits = c(-3,4)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
  # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.area.pce2.p.partial.plot
pdf("peri.area.pce2.p.partial.plot.pdf",width = 4.5, height = 4); peri.area.pce2.p.partial.plot;  dev.off() 



# PACIFIC: PERIPHYTON PER G EELGRASS vs EELGRASS BIOMASS (inverse) PCz2 (Partial, site level)
# First fit model at site level
peri.perg.pcz2.p.partial <- partialResid(rperiphyton.perg.pac ~ rPC2.zos.pac, peri.test.p.11.x, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
peri.perg.pcz2.p.partial <- cbind(peri.perg.pcz2.p.partial, Site)
peri.perg.pcz2.p.partial.plot <- ggplot(peri.perg.pcz2.p.partial, aes(x = xresid, y = -yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass biomass (PCz2) | Z\nPacific") +  
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
peri.perg.pcz2.p.partial.plot
pdf("peri.perg.pcz2.p.partial.plot.pdf",width = 4.5, height = 4); peri.perg.pcz2.p.partial.plot;  dev.off() 


# PACIFIC: PERIPHYTON vs LATITUDE/CLIMATE PCe1 (Partial, site level)
peri.pce1.partial <- partialResid(rperiphyton.perg.pac ~ rPC1.env.global.pac, peri.test.p.11.x, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
peri.pce1.partial <- cbind(peri.pce1.partial, Site)
peri.pce1.partial.plot <- ggplot(peri.pce1.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Ltitude/Climate (PCe1) | Z") +  
  ylab("Periphyton mass | Z") +  
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
peri.pce1.partial.plot
pdf("peri.pce1.partial.plot.pdf",width = 4.5, height = 4); peri.pce1.partial.plot;  dev.off() 


# PACIFIC: PERIPHYTON vs EELGRASS GENETICS FCA2 (Partial, site level)
peri.fca2.p.partial <- partialResid(rperiphyton.perg.pac ~ rFC2.global.pac, peri.test.p.11.x, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
peri.fca2.p.partial <- cbind(peri.fca2.p.partial, Site)
peri.fca2.p.partial.plot <- ggplot(peri.fca2.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass genetics (FCA2) | Z\nPacific") +  
  ylab("Periphyton mass | Z") +  
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
peri.fca2.p.partial.plot
pdf("peri.fca2.p.partial.plot.pdf",width = 4.5, height = 4); peri.fca2.p.partial.plot;  dev.off() 


# PACIFIC: MESOGRAZER BIOMASS PER G PLANT vs ENVIRONMENT PCe3 (Partial, site level)
# Best model at site level
meso.test.p.1 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                    + rPC1.zos.pac + rPC2.zos.pac  
                    , data = ZEN_2014_site_means_Pacific)
# Extract residuals from model
meso.pce3.p.partial <- partialResid(rmesograzer.mass.perg.pac ~ rPC3.env.global.pac, meso.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
meso.pce3.p.partial <- cbind(meso.pce3.p.partial, Site)
meso.pce3.p.partial.plot <- ggplot(meso.pce3.p.partial, aes(x = xresid, y = yresid, col = "forestgreen")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Estuarine conditions (PCe3) | Z\nPacific") +  
  ylab("Mesograzer biomass | Z\nPacific") +  
  # scale_x_continuous(limits = c(-1.4,1.35)) +
  scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.pce3.p.partial.plot
pdf("meso.pce3.p.partial.plot.pdf",width = 4.5, height = 4); meso.pce3.p.partial.plot;  dev.off() 


# PACIFIC: MESOGRAZER BIOMASS PER G PLANT vs GENETIC FC2 (Partial, site level)
# Extract residuals from model
meso.fc2.p.partial <- partialResid(rmesograzer.mass.perg.pac ~ rFC2.global.pac, meso.test.p.1, ZEN_2014_site_means_Pacific) 
# Add Ocean and site to residuals data frame 
Ocean <- ZEN_2014_site_means_Pacific$Ocean
Site <- ZEN_2014_site_means_Pacific$Site
meso.fc2.p.partial <- cbind(meso.fc2.p.partial, Ocean, Site)
meso.fc2.p.partial.plot <- ggplot(meso.fc2.p.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Genetic FC2 | Z\nPacific") +  
  ylab("Mesograzer biomass | Z\nPacific") +  
  # scale_x_continuous(limits = c(30,75)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) # +
  # geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso.fc2.p.partial.plot
pdf("meso.fc2.p.partial.plot.pdf",width = 4.5, height = 4); meso.fc2.p.partial.plot;  dev.off() 


# PACIFIC: MESOGRAZER BIOMASS PER G PLANT vs MEADOW FORM (PCz1) (Partial, site level)
# Extract residuals from model
meso.pcz1.p.partial <- partialResid(rmesograzer.mass.perg.pac ~ rPC1.zos.pac, meso.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
meso.pcz1.p.partial <- cbind(meso.pcz1.p.partial, Site)
meso.pcz1.p.partial.plot <- ggplot(meso.pcz1.p.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass form (PCz1) | Z\nPacific") +  
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
meso.pcz1.p.partial.plot
pdf("meso.pcz1.p.partial.plot.pdf",width = 4.5, height = 4); meso.pcz1.p.partial.plot;  dev.off() 


# PACIFIC: MESOGRAZER BIOMASS PER G PLANT vs (inverse) EELGRASS BIOMASS (PCz2) (Partial, site level)
# Extract residuals from model
meso.perg.pcz2.p.partial <- partialResid(rmesograzer.mass.perg.pac ~ rPC2.zos.pac, meso.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
meso.perg.pcz2.p.partial <- cbind(meso.perg.pcz2.p.partial, Site)
meso.perg.pcz2.p.partial.plot <- ggplot(meso.perg.pcz2.p.partial, aes(x = xresid, y = -yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Eelgrass biomass (PCz1) | Z\nPacific") +  
  ylab("Mesograzer mass | Z\nPacific") +  
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
meso.perg.pcz2.p.partial.plot
pdf("meso.perg.pcz2.p.partial.plot.pdf",width = 4.5, height = 4); meso.perg.pcz2.p.partial.plot;  dev.off() 


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
  ylab("Mesograzer biomass | Z\nPacific") +  
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
pdf("meso.pcz1.p.partial.plot.pdf",width = 4.5, height = 4); meso.pcz1.p.partial.plot;  dev.off() 


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
pdf("meso.area.pcz2.p.partial.plot.pdf",width = 4.5, height = 4); meso.area.pcz2.p.partial.plot;  dev.off() 


# PACIFIC: MESOGRAZER BIOMASS PER AREA vs ESTUARINE COINDITIONS (PCe3) (Partial, site level)
# Best model at site level
meso.area.test.p.1 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                         + rPC1.zos.pac + rPC2.zos.pac  
                         , data = ZEN_2014_site_means_Pacific)
# Extract residuals from model
meso.area.pce3.p.partial <- partialResid(rmesograzer.mass.area.pac ~ rPC3.env.global.pac, meso.area.test.p.1, ZEN_2014_site_means_Pacific) 
# Add site to residuals data frame 
Site <- ZEN_2014_site_means_Pacific$Site
meso.area.pce3.p.partial <- cbind(meso.area.pce3.p.partial, Site)
meso.area.pce3.p.partial.plot <- ggplot(meso.area.pce3.p.partial, aes(x = xresid, y = yresid, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("forestgreen")) +
  xlab("Estuarine conditions (PCe3) | Z\nPacific") +  
  ylab("Mesograzer mass per area | Z\nPacific") +  
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
meso.area.pce3.p.partial.plot
pdf("meso.area.pce3.p.partial.plot.pdf",width = 4.5, height = 4); meso.area.pce3.p.partial.plot;  dev.off() 


# MESOGRAZERS vs ENVIRONMENT PCe3 (RAW, site level)
meso.pce3.plot <- ggplot(ZEN_2014_site_means, aes(x = zPC3.env.global, y = log10.mesograzer.mass.per.g.plant.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Environment PCe3") +  
  ylab("Mesograzer biomass per g") +  
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
meso.pce3.plot
pdf("meso.pce3.plot.pdf",width = 4.5, height = 4); meso.pce3.plot;  dev.off() 


###################################################################################
# FIGURES: 3D PLOTS                                                               #
###################################################################################

# See: http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization

library("plot3D")
library("plot3Drgl") # for interactive viosualization of the 3D plot

# NOTES: We use the same axis ranges in all graphs so Atlantic and Pacific can be directly compared

###################################################################################
# FIGURES: 3D PLOTS OF EELGRASS FORM (PCz1)                                       #
# FIGURES: 3D PLOTS OF EELGRASS BIOMASS (PCz2)                                    #
# FIGURES: 3D PLOTS OF PERIPHYTON PER AREA                                        #
# FIGURES: 3D PLOTS OF MESOGRAZER MASS PER AREA                                   #
###################################################################################


###################################################################################
# FIGURES: 3D PLOTS OF EELGRASS FORM (PCz1)                                       #
###################################################################################

# EELGRASS FORM (PCZ1): PLOT BY CLIMATE/LATITUDE (PCE1) AND GENETIC FC2 - GLOBAL

x <- ZEN_2014_site_means_49$PC1.env.global
y <- ZEN_2014_site_means_49$FC2
z <- ZEN_2014_site_means_49$PC1.zos.site

scatter3D(x, y, z, 
          bty = "g", # back panels and grid lines are visible
          xlim = c(-2.6, 3.4), ylim = c(-900, 2000), zlim = c(-2.6, 3.3), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Climate/latitude", ylab = "Eelgrass genetic FC2", zlab = "Eelgrass form",  
          colvar = ZEN_2014_site_means_49$log10.mesograzer.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          add = FALSE, clab = c("mesograzer", "mass (log)"))

# Add regression plane to 3d figure

# Fit model with only these variables
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 16, cex = 2, col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          bty = "g",
          theta = -120, phi = 5, ticktype = "detailed",
          xlim = c(-2.6, 3.4), ylim = c(-900, 2000), zlim = c(-2.6, 3.3), 
          xlab = "Climate/latitude", ylab = "Eelgrass genetic FC2", zlab = "Eelgrass form",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() # Clear visualization of how Atlantic sites are way above plane of Pacific. 


# EELGRASS FORM (PCZ1): PLOT BY CLIMATE/LATITUDE (PCE1) AND GENETIC FC2 - ATLANTIC

x <- ZEN_2014_site_means_49_Atlantic$PC1.env.global
y <- ZEN_2014_site_means_49_Atlantic$FC2
z <- ZEN_2014_site_means_49_Atlantic$PC1.zos.site

scatter3D(x, y, z, 
          bty = "g", # back panels and grid lines are visible
          xlim = c(-2.6, 3.4), ylim = c(-900, 2000), zlim = c(-2.6, 3.3), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Climate/latitude", ylab = "Eelgrass genetic FC2", zlab = "Eelgrass form",  
          colvar = ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          add = FALSE, clab = c("mesograzer", "mass (log)"))

# Add regression plane to 3d figure

# Fit model with only these variables
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 16, cex = 2, col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          bty = "g",
          theta = -120, phi = 5, ticktype = "detailed",
          xlim = c(-2.6, 3.4), ylim = c(-900, 2000), zlim = c(-2.6, 3.3), 
          xlab = "Climate/latitude", ylab = "Eelgrass genetic FC2", zlab = "Eelgrass form",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() # This is friggin' awesome!


###################################################################################
# FIGURES: 3D PLOTS OF PERIPHYTON PER AREA                                        #
###################################################################################

# PERIPHYTON MASS PER AREA: PLOT BY EELGRASS FORM (PCZ1) AND BIOMASS (PCZ2) - ATLANTIC

x <- ZEN_2014_site_means_49_Atlantic$PC1.zos.site
y <- ZEN_2014_site_means_49_Atlantic$PC2.zos.site
z <- ZEN_2014_site_means_49_Atlantic$log10.periphyton.mass.per.area.site

scatter3D(x, y, z, 
          bty = "g", # back panels and grid lines are visible
          # xlim = c(-3.8, 4), ylim = c(-2.6, 3.3), zlim = c(-1.1, 2.6), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Eelgrass form", ylab = "Eelgrass biomass", zlab = "periphyton mass",  
          colvar = ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          add = FALSE, clab = c("mesograzer", "mass (log)"))

# Add regression plane to 3d figure

# Fit model with only these variables
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 16, cex = 2, col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          bty = "g", ticktype = "detailed",
          theta = 220, phi = 10, 
          # xlim = c(-3.8, 4), ylim = c(-2.6, 3.3), zlim = c(0, 4.7), 
          xlab = "Eelgrass form", ylab = "Eelgrass biomass", zlab = "periphyton mass",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() 



# PERIPHYTON MASS PER AREA: PLOT BY EELGRASS FORM (PCZ1) AND NUTRIENTS (PCE2) - ATLANTIC

x <- ZEN_2014_site_means_49_Atlantic$PC1.zos.site
y <- ZEN_2014_site_means_49_Atlantic$PC2.env.global
z <- ZEN_2014_site_means_49_Atlantic$log10.periphyton.mass.per.area.site

scatter3D(x, y, z, 
          bty = "g", # back panels and grid lines are visible
          xlim = c(-3.8, 3.9), ylim = c(-3.5, 2.7), zlim = c(-1, 2.6), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Eelgrass form", ylab = "Nutrient status", zlab = "periphyton",  
          colvar = ZEN_2014_site_means_49_Atlantic$log10.periphyton.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          add = FALSE, clab = c("periphyton", "mass (log)"))

# Add regression plane to 3d figure
# Fit model with only these variables
fit <- lm(z ~ x + y)
# predict values on regular xy grid
# grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid(x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane

scatter3D(x, y, z, pch = 16, cex = 2, col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          bty = "g",
          theta = 50, phi = 10, 
          ticktype = "detailed",
          xlim = c(-3.8, 3.9), ylim = c(-3.5, 2.7), zlim = c(-1, 2.6), 
          xlab = "Eelgrass form", ylab = "Nutrient status", zlab = "Periphyton mass ",   
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() 


# PERIPHYTON MASS PER AREA: PLOT BY EELGRASS FORM (PCZ1) AND NUTRIENTS (PCE2) - PACIFIC

x <- ZEN_2014_site_means_Pacific$PC1.zos.site
y <- ZEN_2014_site_means_Pacific$PC2.env.global
z <- ZEN_2014_site_means_Pacific$log10.periphyton.mass.per.area.site

scatter3D(x, y, z, 
          bty = "g", # back panels and grid lines are visible
          xlim = c(-3.8, 3.9), ylim = c(-3.5, 2.7), zlim = c(-1, 2.6), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Eelgrass form", ylab = "Nutrient status", zlab = "periphyton",  
          colvar = ZEN_2014_site_means_Pacific$log10.periphyton.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("palegreen", "limegreen", "darkgreen")),
          add = FALSE, clab = c("periphyton", "mass (log)"))

# Add regression plane to 3d figure
# Fit model with only these variables
fit <- lm(z ~ x + y)
# predict values on regular xy grid
# grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid(x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane

scatter3D(x, y, z, pch = 16, cex = 2, 
          col = ramp.col(c("palegreen", "limegreen", "darkgreen")),
          bty = "g",
          theta = 50, phi = 10, 
          ticktype = "detailed",
          xlim = c(-3.8, 3.9), ylim = c(-3.5, 2.7), zlim = c(-1, 2.6), 
          xlab = "Eelgrass form", ylab = "Nutrient status", zlab = "Periphyton mass ",   
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() 



###################################################################################
# FIGURES: 3D PLOTS OF MESOGRAZER MASS PER AREA                                   #
###################################################################################

# MESOGRAZER MASS PER AREA: PLOT BY EELGRASS FORM (PCZ1) AND BIOMASS (PCZ2) - ATLANTIC

x <- ZEN_2014_site_means_49_Atlantic$PC1.zos.site
y <- ZEN_2014_site_means_49_Atlantic$PC2.zos.site
z <- ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.area.site

scatter3D(x, y, z, 
          bty = "g", # back panels and grid lines are visible
          xlim = c(-3.8, 4), ylim = c(-2.6, 3.3), zlim = c(0, 4.7), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Eelgrass form", ylab = "Eelgrass biomass", zlab = "grazer mass ",  
          colvar = ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          add = FALSE, clab = c("mesograzer", "mass (log)"))

# Add regression plane to 3d figure

# Fit model with only these variables
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 16, cex = 2, col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          bty = "g", ticktype = "detailed",
          theta = 220, phi = 10, 
          # xlim = c(-3.8, 4), ylim = c(-2.6, 3.3), zlim = c(0, 4.7), 
          xlab = "Eelgrass form", ylab = "Eelgrass biomass", zlab = "grazer mass ",   
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() # This is friggin' awesome!



# MESOGRAZER MASS PER AREA: PLOT BY EELGRASS BIOMASS (PCZ2) AND ESTUARINE CONDITIONS (PCE3) - ATLANTIC

x2 <- -ZEN_2014_site_means_49_Atlantic$PC2.zos.site
y2 <- ZEN_2014_site_means_49_Atlantic$PC3.env.global
z2 <- ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.area.site

scatter3D(-x2, y2, z2, 
          bty = "g", # back panels and grid lines are visible
          xlim = c(-2.6, 3.3), ylim = c(-2.8, 2.6), zlim = c(0, 4.7), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Eelgrass biomass", ylab = "Estuarine conditons", zlab = "grazer mass ",  
          colvar = ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          add = FALSE, clab = c("mesograzer", "mass (log)"))

# Add regression plane to 3d figure
# Fit model with only these variables
fit2 <- lm(z2 ~ x2 + y2)
# predict values on regular xy grid
# grid.lines = 26
x2.pred <- seq(min(x2), max(x2), length.out = grid.lines)
y2.pred <- seq(min(y2), max(y2), length.out = grid.lines)
x2y2 <- expand.grid(x2 = x2.pred, y2 = y2.pred)
z2.pred <- matrix(predict(fit2, newdata = x2y2), 
                  nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit2)
# scatter plot with regression plane

scatter3D(x2, y2, z2, pch = 16, cex = 2, col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          bty = "g",
          theta =-50, phi = 10, # This looks good for this graph
          ticktype = "detailed",
          xlim = c(-3.8, 4), ylim = c(-2.6, 3.3), zlim = c(0, 4.7), 
          xlab = "Eelgrass biomass", ylab = "Estuarine conditons", zlab = "grazer mass ",   
          surf = list(x = x2.pred, y = y2.pred, z = z2.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() 


# MESOGRAZER MASS PER AREA: PLOT BY EELGRASS BIOMASS (PCZ2) AND ESTUARINE CONDITIONS (PCE3) - PACIFIC

x2 <- -ZEN_2014_site_means_Pacific$PC2.zos.site
y2 <- ZEN_2014_site_means_Pacific$PC3.env.global
z2 <- ZEN_2014_site_means_Pacific$log10.mesograzer.mass.per.area.site

scatter3D(-x2, y2, z2, 
          bty = "g", # back panels and grid lines are visible
          xlim = c(-2.6, 3.3), ylim = c(-2.8, 2.6), zlim = c(0, 4.7), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Eelgrass biomass", ylab = "Estuarine conditons", zlab = "grazer mass ",  
          colvar = ZEN_2014_site_means_Pacific$log10.mesograzer.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("palegreen", "limegreen", "darkgreen")),
          add = FALSE, clab = c("mesograzer", "mass (log)"))

# Add regression plane to 3d figure
# Fit model with only these variables
fit2 <- lm(z2 ~ x2 + y2)
# predict values on regular xy grid
# grid.lines = 26
x2.pred <- seq(min(x2), max(x2), length.out = grid.lines)
y2.pred <- seq(min(y2), max(y2), length.out = grid.lines)
x2y2 <- expand.grid(x2 = x2.pred, y2 = y2.pred)
z2.pred <- matrix(predict(fit2, newdata = x2y2), 
                  nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit2)
# scatter plot with regression plane

scatter3D(x2, y2, z2, pch = 16, cex = 2, col = ramp.col(c("palegreen", "limegreen", "darkgreen")),
          bty = "g",
          theta =-50, phi = 10, # This looks good for this graph
          ticktype = "detailed",
          xlim = c(-3.8, 4), ylim = c(-2.6, 3.3), zlim = c(0, 4.7), 
          xlab = "Eelgrass biomass", ylab = "Estuarine conditons", zlab = "grazer mass ",   
          surf = list(x = x2.pred, y = y2.pred, z = z2.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() 

# MESOGRAZER MASS PER AREA: PLOT BY EELGRASS BIOMASS (PCZ2) AND ESTUARINE CONDITIONS (PCE3) - GLOBAL

x <- ZEN_2014_site_means_49$PC2.zos.site
y <- ZEN_2014_site_means_49$PC3.env.global
z <- ZEN_2014_site_means_49$log10.mesograzer.mass.per.area.site

scatter3D(x, y, z, 
          bty = "g", # back panels and grid lines are visible
          xlim = c(-2.6, 3.3), ylim = c(-2.8, 2.6), zlim = c(0, 4.7), 
          theta = 20, phi = 25, 
          type = "h", ticktype = "detailed",
          xlab = "Eelgrass biomass", ylab = "Estuarine conditons", zlab = "grazer mass ",  
          colvar = ZEN_2014_site_means_49$log10.mesograzer.mass.per.area.site, 
          pch = 16, cex = 1.5, 
          col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          # col.var = as.factor(ZEN_2014_site_means_49$Ocean), 
          # col = c("navyblue", "green"),
          add = FALSE, clab = c("mesograzer", "mass (log)"))
# Can't get tis to color symbols by ocean

# Add regression plane to 3d figure

# Fit model with only these variables
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 16, cex = 2, 
          # col = ramp.col(c("turquoise1", "skyblue1", "navyblue")),
          col.var = as.factor(ZEN_2014_site_means_49$Ocean), 
          col = c("navyblue", "green"),
          
          bty = "g",
          theta = 30, phi = 5, ticktype = "detailed",
          xlim = c(-3.8, 4), ylim = c(-2.6, 3.3), zlim = c(0, 4.7), 
          xlab = "Eelgrass biomass", ylab = "Estuarine conditons", zlab = "grazer mass ",   
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints))

# Create an interactive plot that can be rotated in all directions
plotrgl() 


###################################################################################
# BIVARIATE PLOTS (RAW DATA)                                                      #
###################################################################################


# GLOBAL: CANOPY HEIGHT vs SHOOT DENSITY (RAW data, site level)
can.shoots <- ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.shoots.core.site, y = log10.Zostera.longest.leaf.length.cm.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Shoot density (log)") +  
  ylab("Canopy height (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
can.shoots
pdf("can.shoots.pdf",width = 4.5, height = 4); can.shoots;  dev.off() 


# GLOBAL: EELGRASS ABOVE-GROUND BIOMASS vs SHOOT DENSITY (RAW data, site level)
zag.shoots <- ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.shoots.core.site, y = log10.Zostera.AG.mass.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Shoot density (log)") +  
  ylab("Eelgrass A-G biomass (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
zag.shoots
pdf("zag.shoots.pdf",width = 4.5, height = 4); zag.shoots;  dev.off()


# GLOBAL: EELGRASS ABOVE-GROUND BIOMASS vs GROWTH FORM (PCZ1) (RAW data, site level)
zag_pcz1 <- ggplot(ZEN_2014_site_means, aes(x = PC1.zos.site, y = log10.Zostera.AG.mass.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1)") +  
  ylab("Eelgrass A-G biomass (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
zag_pcz1
pdf("zag_pcz1.pdf",width = 4.5, height = 4); zag_pcz1;  dev.off()


# GLOBAL: CANOPY HEIGHT vs GROWTH FORM (PCZ1) (RAW data, site level)
can_pcz1 <- ggplot(ZEN_2014_site_means, aes(x = PC1.zos.site, y = log10.Zostera.longest.leaf.length.cm.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1)") +  
  ylab("Canopy height (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
can_pcz1
pdf("can_pcz1.pdf",width = 4.5, height = 4); can_pcz1;  dev.off()


# GLOBAL: SHOOT DENSITY vs GROWTH FORM (PCZ1) (RAW data, site level)
shoots_pcz1 <- ggplot(ZEN_2014_site_means, aes(x = PC1.zos.site, y = log10.Zostera.shoots.core.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1)") +  
  ylab("Shoot density (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
shoots_pcz1
pdf("shoots_pcz1.pdf",width = 4.5, height = 4); shoots_pcz1;  dev.off()


# GLOBAL: EELGRASS BIOMASS (-PCz2) vs GROWTH FORM (PCZ1) (RAW data, site level)
pcz2_pcz1 <- ggplot(ZEN_2014_site_means, aes(x = PC1.zos.site, y = -PC2.zos.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1)") +  
  ylab("Eelgrass biomass (-PCz2)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz2_pcz1
pdf("pcz2_pcz1.pdf",width = 4.5, height = 4); pcz2_pcz1;  dev.off()


# GLOBAL: MESOGRAZER BIOMASS PER AREA vs GROWTH FORM (PCZ1) (RAW data, site level)
meso_pcz1 <- ggplot(ZEN_2014_site_means, aes(x = PC1.zos.site, y = log10.mesograzer.mass.per.area.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1)") +  
  ylab("Mesograzer biomass (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso_pcz1
pdf("meso_pcz1.pdf",width = 4.5, height = 4); meso_pcz1;  dev.off()


# GLOBAL: MESOGRAZER BIOMASS PER AREA vs EELGRASS ABOVE-GROUND BIOMASS (RAW data, site level)
meso_zag <- ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.AG.mass.site, y = log10.mesograzer.mass.per.area.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass AG biomass (log)") +  
  ylab("Mesograzer biomass (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso_zag
pdf("meso_zag.pdf",width = 4.5, height = 4); meso_zag;  dev.off()


# GLOBAL: MESOGRAZER BIOMASS PER AREA vs CANOIPY HEIGHT (RAW data, site level)
meso_can <- ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.longest.leaf.length.cm.site, y = log10.mesograzer.mass.per.area.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Canopy height (log)") +  
  ylab("Mesograzer biomass (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso_can
pdf("meso_can.pdf",width = 4.5, height = 4); meso_can;  dev.off()


# GLOBAL: MESOGRAZER BIOMASS PER AREA vs SHOOT DENSITY (RAW data, site level)
meso_shoots <- ggplot(ZEN_2014_site_means, aes(x = log10.Zostera.shoots.core.site, y = log10.mesograzer.mass.per.area.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Shoot density (log)") +  
  ylab("Mesograzer biomass (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
meso_shoots
pdf("meso_shoots.pdf",width = 4.5, height = 4); meso_shoots;  dev.off()



# GLOBAL: SHOOT DENSITY vs GENETIC FCA2 (RAW data, site level)
fc2.shoots <- ggplot(ZEN_2014_site_means, aes(x = FC2, y = log10.Zostera.shoots.core.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Genetic FCA2") +  
  ylab("Shoot density (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
fc2.shoots
pdf("fc2.shoots.pdf",width = 4.5, height = 4); fc2.shoots;  dev.off() 


# GLOBAL: CANOPY HEIGHT vs GENETIC FCA2 (RAW data, site level)
can.fc2 <- ggplot(ZEN_2014_site_means, aes(x = FC2, y = log10.Zostera.longest.leaf.length.cm.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Genetic FCA2") +  
  ylab("Canopy height (log)") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
can.fc2
pdf("can.fc2.pdf",width = 4.5, height = 4); can.fc2;  dev.off() 


# GLOBAL: EELGRASS ABOVE-GROUND BIOMASS vs GENETIC FCA2 (RAW data, site level)
zag.fc2 <- ggplot(ZEN_2014_site_means, aes(x = FC2, y = log10.Zostera.AG.mass.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Genetic FCA2") +  
  ylab("Eelgrass A-G biomass (log)") +  
  # scale_x_continuous(limits = c(-250,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
zag.fc2
pdf("zag.fc2.pdf",width = 4.5, height = 4); zag.fc2;  dev.off()



# GLOBAL: PERIPHYTON vs EELGRASS FORM PCz1 (RAW data, site level)
peri.pcz1 <- ggplot(ZEN_2014_site_means_49, aes(x = PC1.zos.site, y = log10.periphyton.mass.per.area.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Eelgrass form (PCz1)") +  
  ylab("Periphyton mass per area") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri.pcz1
pdf("peri.pcz1.pdf",width = 4.5, height = 4); peri.pcz1;  dev.off() 


# GLOBAL: PERIPHYTON vs NUTRIENTS (PCe2) (RAW data, site level)
peri_pce2 <- ggplot(ZEN_2014_site_means, aes(x = PC2.env.global, y = log10.periphyton.mass.per.area.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Nutrient status (PCe2)") +  
  ylab("Periphyton mass per area") +  
  # scale_x_continuous(limits = c(-165,300)) +
  # scale_y_continuous(position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
peri_pce2
pdf("peri_pce2.pdf",width = 4.5, height = 4); peri_pce2;  dev.off() 


