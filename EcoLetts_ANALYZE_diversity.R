rm(list = ls())

#### 1. Set up environment ####

# Check for and install/load required packages
for (package in c('car', 'lme4', 'lmerTest', 'tidyverse', 'visreg', 'performance',
                  'PerformanceAnalytics', 'data.table', 'vegan', 'FD', 'mgcv', 
                  'effects', 'lattice', 'ggpubr', 'report')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}

library("emmeans")
library("viridis")
library("DHARMa")
library("sjPlot")
library("arm")


# set graphics parameters
theme_set(theme_bw() +  
            theme(legend.position = "bottom", panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 10), axis.text = element_text(size = 8)))


####2. Import & prep data ####
beta.data <- read.csv("data/CAREER_betapart_1x1_2018_2022.csv", row.names = 1)
beta.data$plot <- as.character(beta.data$plot)
beta.data$Seed.density.fac <- factor(beta.data$Seed.density,
                                 levels = c("270", "485", "700", "970"))
beta.data <- beta.data %>% 
  mutate(year = replace(year, year == 2018, 1)) %>%
  mutate(year = replace(year, year == 2019, 2)) %>%
  mutate(year = replace(year, year == 2020, 3)) %>%
  mutate(year = replace(year, year == 2021, 4)) %>%
  mutate(year = replace(year, year == 2022, 5)) 

beta.data$Pool[beta.data$Pool == "Low diversity"] <- "Small"
beta.data$Pool[beta.data$Pool == "High diversity"] <- "Large"
beta.data$Pool <- factor(beta.data$Pool, levels = c("Small", "Large"))

# Standardize (center & scale) seed density and year for interaction. This is done
# because variables are on very different scales, make main effects interpretable
# in the presence of the interaction (main effect is effect at mean of the other
# interacting variable(s)). We standardize by 2x sd as suggest by Gelman 2008
# (to make standardized inputs comparable to unstandardized binary predictor (species pool)

beta.data$scaled.Seed.density <- rescale(beta.data$Seed.density)
beta.data$scaled.Seed.density.fac <- rescale(beta.data$Seed.density.fac)
beta.data$scaled.year <- rescale(beta.data$year)


#### 3. Beta Diversity Model #### 

mod.beta <- lmer(sqrt(beta.tot + 0.1) ~ Pool * scaled.Seed.density * scaled.year 
                 + (1 + scaled.year | site), 
                 data = beta.data)

# check diagnostics: residual plots created through simulation-based appraoch
# in DHARMa package
dotplot(ranef(mod.beta, condVar = TRUE))
simulationOutput <- simulateResiduals(mod.beta, plot = F, use.u = T, n = 1000)
(diag.beta <- plot(simulationOutput))
plotResiduals(simulationOutput, form = beta.data$Pool)
plotResiduals(simulationOutput, form = beta.data$Seed.density)
plotResiduals(simulationOutput, form = beta.data$year)

# Output model statistics
tab_model(mod.beta, show.se = TRUE, show.stat = TRUE, show.df = T) 
summary(mod.beta)
anova(mod.beta)

# Plot model (Figure 2)
supp.labs <- c("Small pool: 8 species", "Large pool: 30 species")
names(supp.labs) <- c("Small", "Large")

(full.mod.beta <- emmip(mod.beta, scaled.year ~ scaled.Seed.density | Pool, mult.name = "Pool", 
                         type = "r", cov.reduce = FALSE,
                         ylab = expression(paste(beta, "-diversity")),
                         xlab = expression(paste("Immigration (seeds ", m^-2, ")"))) +
    facet_grid(~ Pool,  labeller = labeller(Pool = supp.labs))+
    scale_color_viridis(discrete=T) +
    scale_x_continuous(breaks = c(-0.6799070, -0.2630577, 0.1537915, 0.6772766), 
                       labels = c(270, 485, 700, 970)) +
    scale_color_viridis(discrete = TRUE, name = "Year", labels = c('1', '2', '3', '4', '5')) 
)  

# save model predictions
seed <- plot_model(mod.beta, type = "pred", se = F, terms = "scaled.Seed.density")
pool <- plot_model(mod.beta, type = "pred", se = F, terms = "Pool")
year <- plot_model(mod.beta, type = "pred", se = F, terms = "scaled.year")

(dispersal.plot <- seed + ylab(expression(paste(beta, "-diversity"))) +
    geom_line() +
    xlab(expression(paste("Immigration (seeds ", m^-2, ")"))) +
    scale_x_continuous(breaks = c(-0.6799070, -0.2630577, 0.1537915, 0.6772766), 
                       labels = c(270, 485, 700, 970))  +
    theme(plot.title = element_blank())
)  

(pool.plot <- pool + ylab(expression(paste(beta, "-diversity"))) +
    geom_line() +
    xlab(expression(paste("Species pool size"))) +
    theme(plot.title = element_blank())
)

(time.plot <- year + ylab(expression(paste(beta, "-diversity"))) +
    xlab(expression(paste("Year"))) +
    scale_x_continuous(breaks = c(-0.72228010, -0.36788206, -0.01348401, 0.34091403, 0.69531207), 
                       labels = c(1,2,3,4,5))+
    theme(plot.title = element_blank())
)

# Figure 2
#jpeg("figures/Fig2_betaDiversity.jpg", width = 5.5, height = 5, units = "in", res = 1000)
ggarrange(full.mod.beta, labels = c("A"),                                                
          ggarrange(dispersal.plot, pool.plot, time.plot, ncol = 1, nrow = 3,     
                    labels = c("B", "C", "D")),
          widths = c(1.35, 0.8))
#dev.off()


#### 4. Alpha Diversity Model ####

mod.alpha <- lmer(log(alpha.hill1) ~ Pool * scaled.Seed.density * poly(scaled.year, 2) 
                    + (1 + scaled.year | site), 
                    data = beta.data)

# check diagnostics: there is some hard edges in the residuals and quantile plots
dotplot(ranef(mod.alpha, condVar = TRUE))
simulationOutput <- simulateResiduals(mod.alpha, plot = F, use.u = T, n = 1000)
(plot(simulationOutput))
plotResiduals(simulationOutput, form = beta.data$Pool)
plotResiduals(simulationOutput, form = beta.data$Seed.density)
plotResiduals(simulationOutput, form = beta.data$year)

# Output model statistics
tab_model(mod.alpha, show.se = TRUE, show.stat = TRUE, show.df = TRUE)
anova(mod.alpha)
summary(mod.alpha)

# Plot model output (Figure 3)
supp.labs <- c("Small pool: 8 species", "Large pool: 30 species")
names(supp.labs) <- c("Small", "Large")

(full.mod.alpha <- emmip(mod.alpha, scaled.year ~ scaled.Seed.density | Pool, mult.name = "Pool", 
                         type = "r", cov.reduce = FALSE,
                         ylab = expression(paste(alpha, "-diversity")),
                         xlab = expression(paste("Immigration (seeds ", m^-2, ")"))) +
    facet_grid(~ Pool,  labeller = labeller(Pool = supp.labs))+
    scale_color_viridis(discrete=T) +
    scale_x_continuous(breaks = c(-0.6799070, -0.2630577, 0.1537915, 0.6772766), 
                       labels = c(270, 485, 700, 970)) +
    scale_color_viridis(discrete = TRUE, name = "Year", labels = c('1', '2', '3', '4', '5')) 
) 

seed <- plot_model(mod.alpha, type = "pred",  se = F, terms = "scaled.Seed.density")
pool <- plot_model(mod.alpha, type = "pred",  se = F, terms = "Pool")
year <- plot_model(mod.alpha, type = "pred",  se = F, terms = "scaled.year")

(dispersal.plot <- seed + ylab(expression(paste(alpha, "-diversity"))) +
    geom_line() +
    xlab(expression(paste("Immigration (seeds ", m^-2, ")"))) +
    scale_x_continuous(breaks = c(-0.6799070, -0.2630577, 0.1537915, 0.6772766), 
                       labels = c(270, 485, 700, 970))  +
    theme(plot.title = element_blank())
)  

(pool.plot <- pool + ylab(expression(paste(alpha, "-diversity"))) +
    geom_line() +
    xlab(expression(paste("Species pool size"))) +
    scale_y_continuous(breaks = c(1.50, 1.80, 2.10), 
                       labels = c("1.5", "1.8", "2.1")) +
    theme(plot.title = element_blank())
)

(time.plot <- year + ylab(expression(paste(alpha, "-diversity"))) +
    xlab(expression(paste("Year"))) +
    scale_x_continuous(breaks = c(-0.72228010, -0.36788206, -0.01348401, 0.34091403, 0.69531207), 
                       labels = c(1,2,3,4,5))+
    theme(plot.title = element_blank())
)

# Figure 3
#jpeg("figures/Fig3_alphaDiversity.jpg", width = 5.5, height = 5, units = "in", res = 1000)
ggarrange(full.mod.alpha, labels = c("A"),                                                
          ggarrange(dispersal.plot, pool.plot, time.plot, ncol = 1, nrow = 3,     
                    labels = c("B", "C", "D")),
          widths = c(1.35, 0.8))
#dev.off()


#### 5. Post-hoc treatment contrasts ####
# Test for contrast between treatments with similar average density per species.
# Small Pool & 270 Seed density VS. Large Pool & 970 Seed Density
# This lets us test whether species pool effects appear due to lower local density
# vs true diversity effect (species effect persists)

# beta-diversity response
mod.beta.fac <- lmer(sqrt(beta.tot + 0.1) ~ Pool * Seed.density.fac * scaled.year +
                       (1 + scaled.year | site), 
                     data = beta.data)

# check diagnostics: residual plots created through simulation-based approach
# in DHARMa package
dotplot(ranef(mod.beta.fac, condVar = TRUE))
simulationOutput <- simulateResiduals(mod.beta.fac, plot = F, use.u = T, n = 1000)
(diag.beta <- plot(simulationOutput))

# Compare emmeans averaged across Pool
(ls <- emmeans(mod.beta.fac, ~ Pool * Seed.density.fac))
# Get contrasts (Tukey's adjusted, 0.05)
(co <- contrast(ls, "pairwise", adjust = "none"))

#Small Seed.density.fac270 - Large Seed.density.fac970:
# estimate     SE    df   t.ratio   p.value
#-0.02759   0.0120  2566  -2.297    0.0217


# alpha-diversity response
mod.alpha.fac <- lmer(log(alpha.hill1) ~ Pool * Seed.density.fac * poly(scaled.year, 2) +
                        (1 + scaled.year | site), 
                      data = beta.data)

# check diagnostics: residual plots created through simulation-based approach
# in DHARMa package
dotplot(ranef(mod.alpha.fac, condVar = TRUE))
simulationOutput <- simulateResiduals(mod.alpha.fac, plot = F, use.u = T, n = 1000)
(diag.beta <- plot(simulationOutput))
# Compare emmeans averaged across Pool
(ls <- emmeans(mod.alpha.fac, ~ Pool * Seed.density.fac))
# Get contrasts (Tukey's adjusted, 0.05)
(co <- contrast(ls, "pairwise", adjust = "none"))
#Small Seed.density.fac270 - Large Seed.density.fac970:
# estimate   SE     df   t.ratio   p.value
#-0.5434   0.0489  2556  -11.123   <.0001

