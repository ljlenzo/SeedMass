# Q1: Where Are Seeds Bigger? #####

# =========================== #
## Load required libraries #####
# =========================== #

# Data manipulation and analysis
library(dplyr)
library(lme4)
library(car)

# Visualization
library(ggplot2)
library(ggeffects)

# Spatial packages
library(sf)
library(sp)
library(raster)
library(tmap)

# =========================== #
## Set working directory & seed #####
# =========================== #

set.seed(123)

# =========================== #
## Load data #####
# =========================== #

# Using data prepped for Q1, full dataset
# all species and populations
species <- readRDS("Q1spatialSOS2024.rds")

# Full interactive models #####

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

# Growth form models #####

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

# Duration models #####

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

# Simpler models #####

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

summary(spatialGLM_grass)
summary(spatialGLM_forb)
summary(spatialGLM_shrub)

results_df <-summary.glm(spatialGLM_grass)$coefficients
results_df <- round(results_df, 3)
conf <- confint(spatialGLM_grass)
conf <- round(conf, 3)

results_df <-summary.glm(spatialGLM_forb)$coefficients
results_df <- round(results_df, 3)
conf <- confint(spatialGLM_forb)
conf <- round(conf, 3)

results_df <-summary.glm(spatialGLM_shrub)$coefficients
results_df <- round(results_df, 3)
conf <- confint(spatialGLM_shrub)
conf <- round(conf, 3)
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

# Fit GLMM models for each duration type
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

# Plot effects of latitude #####

# =========================== #
## Shrub data prediction ######
# =========================== #

shrubpred <- predict_response(spatialGLM_shrub, 
                              terms = "coords.x2",
                              ci_level = 0.95,
                              back.transform = FALSE) 

# Prepare prediction data
shrubpred$coords.x2 <- shrubpred$x
shrubpred$mg.seed <- shrubpred$predicted
shrubpred$Growth.Habit.short <- 'Shrub'

# Plot predictions for shrub data
ggplot(shrubpred, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)

# =========================== #
## Forb data prediction #######
# =========================== #

forbpred <- predict_response(spatialGLM_forb, 
                             terms = "coords.x2",
                             ci_level = 0.95,
                             back.transform = FALSE) 

# Prepare prediction data
forbpred$coords.x2 <- forbpred$x
forbpred$mg.seed <- forbpred$predicted
forbpred$Growth.Habit.short <- 'Forb/herb'

# Plot predictions for forb data
ggplot(forbpred, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)

# =========================== #
## Grass data prediction ######
# =========================== #

grasspred <- predict_response(spatialGLM_grass, 
                              terms = "coords.x2",
                              ci_level = 0.95,
                              back.transform = FALSE) 

# Prepare prediction data
grasspred$coords.x2 <- grasspred$x
grasspred$mg.seed <- grasspred$predicted
grasspred$Growth.Habit.short <- 'Graminoid'

# Plot predictions for grass data
ggplot(grasspred, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)

# =========================== #
## Full dataset prediction #####
# =========================== #

fullPred <- predict_response(spatialGLM_full, 
                             terms = "coords.x2",
                             ci_level = 0.95,
                             back.transform = FALSE) 

# Prepare prediction data
fullPred$coords.x2 <- fullPred$x
fullPred$mg.seed <- fullPred$predicted
fullPred$Growth.Habit.short <- 'Full'

# Plot predictions for full dataset
ggplot(fullPred, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)

# Simple linear regression plots #####

# =========================== #
## Latitude plot ######
# =========================== #

# Calculate ranges
x_range <- range(species$coords.x2)
y_range <- range(log(species$mg.seed))

# Plot Latitude vs Seed Mass
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
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(name = "Growth Habit", values = habit_colors) +
  scale_x_continuous(expand = c(0, 0)) +  # Remove buffer space around x-axis limits
  scale_y_continuous(expand = c(0, 0)) +  # Remove buffer space around y-axis limits
  coord_fixed(ratio = 1.452661)  # Set aspect ratio to 1 for a square plot

# Add regression line and confidence interval
lat2 <- lat + 
  geom_smooth(method = "lm", se = TRUE, 
              color = 'blue', fill = 'grey',
              linewidth = 0.5,
              inherit.aes = FALSE, aes(x = coords.x2, y = log(mg.seed)))

# Save plot
ggsave("Q1lat.png", lat2, width = 10, height = 10, units = "in", dpi = 300)

# =========================== #
## Year plot ######
# =========================== #

# Plot Year vs Seed Mass
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
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(name = "Growth Habit", values = habit_colors) +
  coord_fixed(ratio = 1.404915)

# Add regression line and confidence interval
year2 <- year + 
  geom_smooth(method = "lm", se = TRUE, 
              color = 'red', fill = 'grey',
              linewidth = 0.5)

# Save plot
ggsave("Q1year.png", year2, width = 10, height = 10, units = "in", dpi = 300)

# =========================== #
## Growth form summary table #####
# =========================== #

# Summarize unique species counts by growth habit
growth_form_table <- species %>%
  group_by(Growth.Habit.short) %>%
  summarise(Unique_Species_Count = n_distinct(Name2))

# IDW mapping analysis #####

# =========================== #
## Z-transform mg.seed by each species ####
# =========================== #
species <- species %>%
  group_by(Name2) %>%
  mutate(mg.seed_z = scale(mg.seed))  # Z-transform seed mass

# =========================== #
## Load and plot lower 48 states shapefile ####
# =========================== #
states <- st_read('../States_shapefile/States_shapefile.shp')
lower_48_abb <- state.abb[!(state.abb %in% c("AK", "HI"))]
lower_48 <- states[states$State_Code %in% lower_48_abb, ]
lower_48 <- lower_48[-c(1, 2, 4, 5, 6)]  # Remove unnecessary states
plot(lower_48)

# =========================== #
## Crop the shapefile ####
# =========================== #
bbox_new <- st_bbox(lower_48)
xrange <- bbox_new$xmax - bbox_new$xmin
yrange <- bbox_new$ymax - bbox_new$ymin
bbox_new[3] <- -100  # Update bounding box for the region of interest
bbox_new <- bbox_new %>% st_as_sfc()
lower_48 <- st_crop(lower_48, bbox_new)
plot(lower_48)

# =========================== #
## Create spatial points data frame ####
# =========================== #
species2 <- na.omit(species)
xy <- species2[, c(4, 5)]  # Select relevant columns for coordinates

speciesSPDF <- SpatialPointsDataFrame(coords = xy, data = species2, proj4string = CRS("+proj=longlat +datum=NAD83"))
speciesSPDF@bbox <- st_bbox(lower_48)  # Set bounding box for the spatial data

# =========================== #
## Create IDW interpolated surface ####
# =========================== #
grd <- as.data.frame(spsample(speciesSPDF, "regular", n = 13004))  # Create grid for interpolation
names(grd) <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd) <- TRUE
fullgrid(grd) <- TRUE
proj4string(speciesSPDF) <- proj4string(speciesSPDF)
proj4string(grd) <- proj4string(speciesSPDF)

# Perform IDW interpolation
P.idw <- gstat::idw(mg.seed_z ~ 1, speciesSPDF, newdata = grd, idp = 2.0)
r <- raster(P.idw)
r.m <- mask(r, lower_48)  # Mask raster to lower 48 states

# =========================== #
## Plot IDW interpolation ####
# =========================== #
tm_shape(lower_48) + tm_polygons() +
  tm_shape(r.m) + 
  tm_raster(n = 10, palette = "RdBu", auto.palette.mapping = FALSE,
            title = "Predicted seed mass (mg, log-transformed)") +
  tm_legend(legend.outside = TRUE)

# =========================== #
## Jackknife technique to estimate confidence intervals ####
# =========================== #
img <- gstat::idw(mg.seed_z ~ 1, speciesSPDF, newdata = grd, idp = 2.0)
n <- length(speciesSPDF)
Zi <- matrix(nrow = length(img$var1.pred), ncol = n)

# Jackknife procedure
st <- stack()
for (i in 1:n) {
  Z1 <- gstat::idw(mg.seed_z ~ 1, speciesSPDF[-i, ], newdata = grd, idp = 2.0)
  st <- addLayer(st, raster(Z1, layer = 1))
  Zi[, i] <- n * img$var1.pred - (n - 1) * Z1$var1.pred
  print(i)
}

# Jackknife estimator
Zj <- as.matrix(apply(Zi, 1, sum, na.rm = TRUE) / n)
c1 <- apply(Zi, 2, '-', Zj)
c1 <- apply(c1^2, 1, sum, na.rm = TRUE)
CI <- sqrt(1 / (n * (n - 1)) * c1)

# =========================== #
## Plot confidence intervals ####
# =========================== #
img.sig <- img
img.sig$v <- CI / img$var1.pred  # Calculate significance

# Clip confidence raster to region of interest (e.g., Texas or other areas)
r <- raster(img.sig, layer = "v")
r.m <- mask(r, W)  # Mask confidence raster to region

# Plot confidence intervals
tm_shape(r.m) + tm_raster(n = 7, title = "95% confidence interval") +
  tm_shape(speciesSPDF) + tm_dots(size = 0.2) +
  tm_legend(legend.outside = TRUE)

# Summarize the entire map for paper #####
# =========================== #
## Load Required Libraries #####
# =========================== #
library(prism)
library(reshape2)
library(raster)

# =========================== #
## Set Up Directories and Data #####
# =========================== #
prism_set_dl_dir("SummaryPrism")

# Convert to an sf object (points from species data)
sf_points <- st_as_sf(species, coords = c("coords.x1", "coords.x2"), crs = 4326)

# Create polygon (convex hull) around points and add buffer (10 km)
study_area_polygon <- st_buffer(st_convex_hull(st_union(sf_points)), dist = 10000)

# =========================== #
## Download Climate Data #####
# =========================== #
# Download 20 years of annual temperature and precipitation data
get_prism_annual(type = "tmean", years = 2001:2020, keepZip = FALSE)
get_prism_annual(type = "ppt", years = 2001:2020, keepZip = FALSE)

# Load climate data
AnnualTemp <- pd_stack(prism_archive_subset('tmean','annual', years = 2001:2020))
AnnualPpt <- pd_stack(prism_archive_subset('ppt','annual', years = 2001:2020))

# Transform study area polygon to the same CRS as the climate data
study_area_polygon <- st_transform(study_area_polygon, crs(AnnualTemp))

# =========================== #
## Clip and Extract Data #####
# =========================== #
# Crop and mask rasters to the polygon area of interest
temp_clip <- mask(crop(AnnualTemp, extent(study_area_polygon)), study_area_polygon)
precip_clip <- mask(crop(AnnualPpt, extent(study_area_polygon)), study_area_polygon)

# Extract mean values for each year (temperature and precipitation)
temp_values <- extract(temp_clip, study_area_polygon, fun = mean, na.rm = TRUE)
ppt_values <- extract(precip_clip, study_area_polygon, fun = mean, na.rm = TRUE)

# =========================== #
## Prepare Data for Plotting #####
# =========================== #
# Ensure 20 values (one for each year)
mean_temp_vector <- as.vector(temp_values)
mean_ppt_vector <- as.vector(ppt_values)

# Create a dataframe for plotting
years <- 2001:2020
mean_temp_df <- data.frame(
  year = years,
  mean_temp = mean_temp_vector,
  mean_ppt = mean_ppt_vector
)

# Reshape data for easier plotting
mean_data_long <- reshape2::melt(mean_temp_df, id.vars = "year")

# =========================== #
## Plot Precipitation and Temperature #####
# =========================== #
# Set axis limits and scaling factors for combined plot
precip_limits <- c(0, 650)    # Precipitation limits
temp_limits <- c(7, 12)       # Temperature limits

# Scale factor for combining precipitation and temperature on the same plot
b <- diff(precip_limits) / diff(temp_limits)
a <- precip_limits[1] - b * temp_limits[1]

# Create the plot
ggplot(mean_temp_df, aes(year, mean_ppt)) +
  geom_col() +  # Bar plot for precipitation
  geom_line(aes(y = a + mean_temp * b), color = "red") +  # Line for temperature
  scale_y_continuous(
    name = "Precipitation (mm)", 
    limits = precip_limits, 
    sec.axis = sec_axis(~(. - a) / b, name = "Temperature (°C)")  # Secondary axis for temperature
  ) +
  scale_x_continuous("Year") +  # X-axis for the years
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.line.y.right = element_line(color = "red"),
    axis.ticks.y.right = element_line(color = "red"),
    axis.text.y.right = element_text(color = "red"),
    axis.title.y.right = element_text(color = "red")
  )

# =========================== #
## Bounding Box and Area Calculation #####
# =========================== #
# Get bounding box of the study area polygon
bbox <- st_bbox(study_area_polygon)

# Print latitudinal and longitudinal extent of the polygon
lat_range <- c(bbox["ymin"], bbox["ymax"])
lon_range <- c(bbox["xmin"], bbox["xmax"])

print(paste("Latitude range:", lat_range[1], "to", lat_range[2]))
print(paste("Longitude range:", lon_range[1], "to", lon_range[2]))

# Project polygon to UTM for area calculation
polygon_proj <- st_transform(study_area_polygon, 32633)  # Replace with correct UTM zone

# Calculate area in square kilometers
area_km2 <- st_area(polygon_proj) / 10^6  
print(paste("Polygon Area:", round(area_km2, 2), "km²"))
