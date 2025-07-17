# Q1: Where Are Seeds Bigger? #####

# =========================== #
## Load libraries, data #####
# =========================== #

library(dplyr)
library(lme4)
library(car)
library(ggplot2)
library(ggeffects)
library(sf)
library(sp)
library(raster)
library(tmap)

set.seed(123)

# Full Interactive Models #####

# GLM
spatialGLM_full <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                       * (Growth.Habit.short + Duration.y), 
                       data = species)
# GLmM with family
spatialGLMM_full <- lmer(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                         * (Growth.Habit.short + Duration.y) 
                         + (1 | FAMILY), 
                         data = species)

summary(spatialGLM_full)
# lat + long not significant
# lat x trees significant but no other growth form
# elevation significant and positive for all growth forms
# elevation negative for all duration
# only major differences here are tree x lat and 
# elevation x duration but not among

# Table of combinations: Family, Growth Habit, and Duration
table_combinations <- table(species$FAMILY, 
                            species$Growth.Habit.short, 
                            species$Duration.y)

table_combinations
# poor distribution among combinations, 
# full mixed model is not great, 
# break down either growth form or duraction

summary(spatialGLMM_full)
car::Anova(spatialGLMM_full)

# Calculate 95% confidence intervals
conf_interval <- confint(spatialGLMM_full, level = 0.95)
print(conf_interval)

# Growth Form Models #####

# GLM 
spatialGLM_full <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                       * (Growth.Habit.short), 
                       data = species)

# GLmM with family
spatialGLMM_full <- lmer(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                         * (Growth.Habit.short) 
                         + (1 | FAMILY), 
                         data = species)

summary(spatialGLM_full)
summary(spatialGLMM_full)

car::Anova(spatialGLM_full)
car::Anova(spatialGLMM_full)

results_df <-summary.glm(spatialGLM_full)$coefficients
results_df <- round(results_df, 3)
conf <- confint(spatialGLM_full)
conf <- round(conf, 3)

# direction is consistent among growth form 
# and only significant for shrub and tree
summary(spatialGLMM_full)
car::Anova(spatialGLMM_full)

summary(spatialGLMM_full)

summary_model <- summary(spatialGLMM_full)
coefficients <- summary_model$coefficients
coefficients <- round(coefficients, 3)

conf <- confint(spatialGLMM_full)
conf <- round(conf, 3)

# Calculate 95% confidence intervals
conf_interval <- confint(spatialGLMM_full, level = 0.95)
print(conf_interval)

# Duration Models #####

spatialGLM_full <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                       * (Duration.y), 
                       data = species)

spatialGLMM_full <- lmer(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                         * (Duration.y) 
                         + (1 | FAMILY), 
                         data = species)

summary(spatialGLM_full)
summary(spatialGLMM_full)

car::Anova(spatialGLM_full)
car::Anova(spatialGLMM_full)

results_df <-summary.glm(spatialGLM_full)$coefficients
results_df <- round(results_df, 3)
conf <- confint(spatialGLM_full)
conf <- round(conf, 3)

summary_model <- summary(spatialGLMM_full)
coefficients <- summary_model$coefficients
coefficients <- round(coefficients, 3)

conf <- confint(spatialGLMM_full)
conf <- round(conf, 3)

# Simpler Models #####

# All growth forms together
spatialGLM_full <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation), 
                       data = species)

spatialGLMM_full <- lmer(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                         + (1 | FAMILY), 
                         data = species)

summary(spatialGLM_full)
summary(spatialGLMM_full)

anova(spatialGLMM_full, spatialGLM_full)

car::Anova(spatialGLMM_full)

# Calculate 95% confidence intervals
conf_interval <- confint(spatialGLM_full, level = 0.95)
print(conf_interval)

conf_interval <- confint(spatialGLMM_full, level = 0.95)
print(conf_interval)

## Break models down by growth form #####

# Subset data by growth form
Grass <- species[species$Growth.Habit.short == 'Graminoid',]
forbs <- species[species$Growth.Habit.short == 'Forb/herb',]
shrub <- species[species$Growth.Habit.short == 'Shrub',]

# Fit models for each growth form
spatialGLM_grass <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation), 
                       data = Grass)
spatialGLM_forb <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation), 
                       data = forbs)
spatialGLM_shrub <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation), 
                       data = shrub)

## Break models down by duration #####

# Subset data by duration
annual <- species[species$Duration.y == 'Annual',]
biennial <- species[species$Duration.y == 'Biennial',]
perennial <- species[species$Duration.y == 'Perennial',]

# Fit models for each duration type
spatialGLM_annual <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation), 
                        data = annual)
spatialGLM_perennial <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation), 
                       data = perennial)
spatialGLM_biennial <- glm(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation), 
                         data = biennial)

summary(spatialGLM_annual)
summary(spatialGLM_perennial)
summary(spatialGLM_biennial)

results_df <-summary.glm(spatialGLM_annual)$coefficients
results_df <- round(results_df, 3)
conf <- confint(spatialGLM_annual)
conf <- round(conf, 3)

results_df <-summary.glm(spatialGLM_perennial)$coefficients
results_df <- round(results_df, 3)
conf <- confint(spatialGLM_perennial)
conf <- round(conf, 3)

results_df <-summary.glm(spatialGLM_biennial)$coefficients
results_df <- round(results_df, 3)
conf <- confint(spatialGLM_biennial)
conf <- round(conf, 3)

spatialGLMM_an <- lmer(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                         + (1 | FAMILY), 
                         data = annual)

spatialGLMM_bi <- lmer(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                         + (1 | FAMILY), 
                         data = biennial)

spatialGLMM_per <- lmer(log(mg.seed) ~ (coords.x1 + coords.x2 + elevation) 
                            + (1 | FAMILY), 
                            data = perennial)

summary_model <- summary(spatialGLMM_an)
coefficients <- summary_model$coefficients
coefficients <- round(coefficients, 3)
conf <- confint(spatialGLMM_an)
conf <- round(conf, 3)

summary_model <- summary(spatialGLMM_bi)
coefficients <- summary_model$coefficients
coefficients <- round(coefficients, 3)
conf <- confint(spatialGLMM_bi)
conf <- round(conf, 3)

summary_model <- summary(spatialGLMM_per)
coefficients <- summary_model$coefficients
coefficients <- round(coefficients, 3)
conf <- confint(spatialGLMM_per)
conf <- round(conf, 3)
# Breaking down by growth form, only shrubs change
# But when considering all of our data, there is no strong spatial pattern

########################### simple linear regression plots ########################### 
x_range <- range(species$coords.x2)
x_range[2]-x_range[1]
y_range <- range(log(species$mg.seed))
y_range[2]-y_range[1]
(x_range[2]-x_range[1])/(y_range[2]-y_range[1])

lat <- ggplot(data = species, aes(x = coords.x2, y = log(mg.seed), color = Growth.Habit.short)) +
  geom_point(size = 3, alpha = 0.3) +
  xlab(NULL) + ylab(NULL) + 
  theme_bw() +
  theme(legend.position = 'none', 
        plot.title = element_text(size = 26, colour = "black", face = "bold"),
        axis.text = element_text(size = 0),
        axis.title = element_text(size = 22, colour = "black", face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank()) +  # Remove tick marks on y-axis
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(name = "Growth Habit", values = habit_colors) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove buffer space around x-axis limits
  scale_y_continuous(expand = c(0, 0)) +  # Remove buffer space around y-axis limits
  coord_fixed(ratio = 1.452661)  # Set aspect ratio to 1 for a square plot

lat2 <- lat + 
        geom_smooth(method = "lm", se = TRUE, 
                    color = 'blue', fill = 'grey',
                    linewidth = 0.5,
                    inherit.aes = FALSE, aes(x = coords.x2, y = log(mg.seed)))

ggsave("Q1lat.png", lat2, width = 10, height = 10, units = "in", dpi = 300)

x_range <- range(species$Year)
x_range[2]-x_range[1]
y_range <- range(log(species$mg.seed))
y_range[2]-y_range[1]
(x_range[2]-x_range[1])/(y_range[2]-y_range[1])

year <- ggplot(data = species, aes(x = Year, y = log(mg.seed), color = Growth.Habit.short)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.3)  +
  xlab(NULL) + ylab(NULL) + 
  theme_bw() +
  theme(legend.position = 'none', 
        plot.title = element_text(size = 0, colour = "black", face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22, colour = "black", face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank()) +  # Remove tick marks on y-axis
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(name = "Growth Habit", values = habit_colors) +
  coord_fixed(ratio = 1.404915)

year2 <- year + 
  geom_smooth(method = "lm", se = TRUE, 
              color = 'red', fill = 'grey',
              linewidth = 0.5)

ggsave("Q1year.png", year2, width = 10, height = 10, units = "in", dpi = 300)

growth_form_table <- species %>%
  group_by(Growth.Habit.short) %>%
  summarise(Unique_Species_Count = n_distinct(Name2))  

########################### IDW mapping analysis ###########################

#z transform mg.seed by each species 

species <- species %>%
  group_by(Name2) %>%
  mutate(mg.seed_z = scale(mg.seed))

# Load and plot lower 48 states shapefile
states = st_read('../States_shapefile/States_shapefile.shp', )
lower_48_abb = state.abb[!(state.abb %in% c("AK", "HI"))]
lower_48 = states[states$State_Code %in% lower_48_abb,]
lower_48 <- lower_48[-c(1,2,4,5,6)]
plot(lower_48)

# Crop the lower 48 states shapefile to a new bounding box
bbox_new <- st_bbox(lower_48)
xrange <- bbox_new$xmax - bbox_new$xmin
yrange <- bbox_new$ymax - bbox_new$ymin
bbox_new[3] <- -100
bbox_new <- bbox_new %>% st_as_sfc()
lower_48 <- st_crop(lower_48, bbox_new)
plot(lower_48)

# Create a spatial points data frame
species2 <- na.omit(species)
xy <- species2[, c(4,5)]

speciesSPDF <- SpatialPointsDataFrame(coords = xy, data = species2, proj4string = CRS("+proj=longlat +datum=NAD83"))
speciesSPDF@bbox <- st_bbox(lower_48)

# Create an interpolated surface
grd <- as.data.frame(spsample(speciesSPDF, "regular", n = 13004))
names(grd) <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd) <- TRUE
fullgrid(grd) <- TRUE
proj4string(speciesSPDF) <- proj4string(speciesSPDF)
proj4string(grd) <- proj4string(speciesSPDF)
P.idw <- gstat::idw(mg.seed_z ~ 1, speciesSPDF, newdata = grd, idp = 2.0)
r <- raster(P.idw)
r.m <- mask(r, lower_48)

# Plot 
tm_shape(lower_48) + tm_polygons() +
  tm_shape(r.m) + 
  tm_raster(n = 10, palette = "RdBu", auto.palette.mapping = FALSE,
            title = "Predicted seed mass (mg, log-transformed)") +
  tm_legend(legend.outside = TRUE)

# Jackknife technique to estimate CIs

img <- gstat::idw(mg.seed_z~1, speciesSPDF, newdata=grd, idp=2.0)
n   <- length(speciesSPDF)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate, repeat for each point
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(mg.seed_z~1, speciesSPDF[-i,], newdata=grd, idp=2.0)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
  print(i)
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Texas
r <- raster(img.sig, layer="v")
r.m <- mask(r, W)

# Plot the map
tm_shape(r.m) + tm_raster(n=7,title="95% confidence interval") +
  tm_shape(speciesSPDF) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
