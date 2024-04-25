rm(list = ls())

# Check for and install required packages
for (package in c('tidyverse', 'vegan',  'betapart', 'arm', 'sjPlot', 
                  'lattice', 'ggpubr', 'viridis', 'PerformanceAnalytics')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}

theme_set(theme_bw() +  
            theme(legend.position = "bottom", panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 10), axis.text = element_text(size = 8)))

####1. Import data ####
trts <- read.csv("data/CAREER_betapart_dbrda_2018_2022.csv", row.names = 1)

trts <- trts %>% 
  mutate(year = replace(year, year == 2018, 1)) %>%
  mutate(year = replace(year, year == 2019, 2)) %>%
  mutate(year = replace(year, year == 2020, 3)) %>%
  mutate(year = replace(year, year == 2021, 4)) %>%
  mutate(year = replace(year, year == 2022, 5)) 

trts$Pool[trts$Pool == "Low diversity"] <- "Small"
trts$Pool[trts$Pool == "High diversity"] <- "Large"
trts$Pool <- factor(trts$Pool, levels = c("Small", "Large"))

# Standardize (center & scale) seed density and year for interaction. This is done
# because variables are on very different scales, make main effects interpretable
# in the presence of the interaction (main effect is effect at mean of the other
# interacting variable(s)). We standardize by 2x sd as suggest by Gelman 2008
# (to make standardized inputs comparable to unstandardized binary predictor (species pool)

trts$scaled.Seed.density <- rescale(trts$Seed.density)
trts$scaled.year <- rescale(trts$year)


#### Models ####

# Total variation in composition explained
lm.tot.mod <- lm(log(tot_dbra_r2) ~ scaled.Seed.density + Pool + scaled.year, 
                 data = trts)
par(mfrow = c(2,2))
plot(lm.tot.mod)
summary(lm.tot.mod)
car::avPlots(lm.tot.mod)
tab_model(lm.tot.mod, show.se = TRUE, show.stat = TRUE, show.df = T) 

# Balanced variation in abundance
lm.bal.mod <- lm(log(bal_dbra_r2) ~ scaled.Seed.density + Pool + scaled.year, 
                 data = trts)
par(mfrow = c(2,2))
plot(lm.bal.mod)
summary(lm.bal.mod)
car::avPlots(lm.bal.mod)
tab_model(lm.bal.mod, show.se = TRUE, show.stat = TRUE, show.df = T) 

# Abundance gradient
lm.gra.mod <- lm(log(gra_dbra_r2) ~ scaled.Seed.density + Pool + scaled.year, 
                 data = trts)
par(mfrow = c(2,2))
plot(lm.gra.mod)
summary(lm.gra.mod)
car::avPlots(lm.gra.mod)
tab_model(lm.gra.mod, show.se = TRUE, show.stat = TRUE, show.df = T) #

# Plot effects
p <- plot_models(lm.tot.mod, lm.bal.mod, lm.gra.mod, std.est = "std2", 
           axis.labels = c("Year", "Species pool\nsize", "Immigration"),
            m.labels = c("Beta-diversity", "Balanced variation in abundance", "Gradient in\nabundance"),
            legend.title = NULL)
(SS.fig <- p + ylim(-0.6, 0.6) + 
    ylab("Standardized Effect") +
    scale_color_viridis(discrete = T, direction = -1) +
    ggtitle("Composition-Environment Relationship")
)

# Figure 4
#jpeg("figures/Fig4_speciesSorting.jpg", width = 4.5, height = 5, units = "in", res = 1000)
SS.fig
#dev.off()
