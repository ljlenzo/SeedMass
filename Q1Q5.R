########################### Q1 Seed Mass Variation ########################### 
library(Boruta)
library(dplyr)
library(car)
library(ape)
library(geiger)
library(nlme)
library(ggplot2)
library(RColorBrewer)
library(sf)

set.seed(123)

########################### Bring in the Tree ###########################
tree <- read.tree('../v0.1/ALLMB.tre')

capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

species$GENUS.tree <- tolower(species$GENUS.tree)
species$GENUS.tree <- capFirst(species$GENUS.tree)
species$TreeName <- paste(species$GENUS.tree,species$SPECIES.tree,sep ='_')

spnames <- c(unique(species$Name2))
spnames <- sort(spnames)

correctedSp <- unique(species$TreeName)
# Get the species names from the phylogeny tree
phylo_species <- tree$tip.label

# Identify missing species
missing_species <- setdiff(spnames, phylo_species)
Treecheck <- species[species$Name2 %in% missing_species, ]
Treecheck <- Treecheck$TreeName
TreeCorrect <- intersect(Treecheck, tree$tip.label)

rows_to_update <- species$TreeName %in% TreeCorrect
species$Name2[rows_to_update] <- species$TreeName[rows_to_update]

common_species <- intersect(species$Name2, tree$tip.label)
common_species_data <- species[species$Name2 %in% common_species, ]
subset_tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% common_species])

name.check(subset_tree, species$Name2)

common_species <- intersect(species$Name2, subset_tree$tip.label)
spnames2 <- species$Name2
phylo_species2 <- subset_tree$tip.label
missing_species <- setdiff(spnames, phylo_species)

species <- species[species$Name2 %in% common_species, ]

########################### Summarize variables at species level ########################### 
colnames(species)
colnames(species[,c(12,15:32)])

#we want the mean, min, max and range for the following:

## Year, Ydays 39, 73
## Slope, Aspect, packm, heatload, SixMonthPrecip
# 48, 49, 61, 65, 69, 83
## DeltFallMin, DeltFallMax, DeltSpringMax, DeltPrecip
#  84,85,86,88
## AWC SoilpH  clay organic
# 92,93,95,96
## AnnMinTemp
## 103

# 20,40, 74,49, 50, 62, 70, 84, 85,86,87,89, 93,94,96,97,104

#also grab 75:81 
#"Sp.MaxLat"          "Sp.MinLat"         
#[76] "Sp.MaxLong"         "Sp.MinLong"         "Sp.LatExtent"      
#[79] "Sp.LongExtent"      "ChullArea.km2."

#and we want to keep duration, growth habit short, and genus tree
#8,10,111

#72 is our response CV
#72 ~ 20,40, 74,49, 50, 62, 70, 84, 85,86,87,89, 93,94,96,97,104, 75:81 , 9,11

# Aggregate data by mean, min, max, and range
#  [1] "mg.seed"       "Year"          "Ydays"         "slope2"       
#[5] "aspect2"       "packm"         "heatload"      "sixmonh20"    
#[9] "DeltFallMin"   "DeltFallMax"   "DeltSpringMax" "DeltPrecip"   
#[13] "awc2"          "SoilpH"        "clay"          "organic"      
#[17] "AnnMinTemp" 

speciesMean <- species %>%
  dplyr::group_by(Name2) %>%
  dplyr::summarise_at(c(11,16:31), mean,na.rm = TRUE)

speciesMax <- species %>%
  dplyr::group_by(Name2) %>%
  dplyr::summarise_at(c(11,14:31), max,na.rm = TRUE)

speciesMin <- species %>%
  dplyr::group_by(Name2) %>%
  dplyr::summarise_at(c(11,14:31), min,na.rm = TRUE)

species_sf <- st_as_sf(species, coords = c("coords.x1", "coords.x2"), crs = 4269)
species_sf <- st_transform(species_sf, crs = 32610)

species_areas <- species_sf %>%
  dplyr::group_by(Name2) %>%
  dplyr::summarize(ChullArea.km = st_area(st_convex_hull(st_union(geometry))) / 1e6)  # Convert from m² to km²

species$ChullArea.km <- species_areas$ChullArea.km[match(species$Name2, species_areas$Name2)]

speciesCat <- species %>%
  dplyr::group_by(Name2) %>%
  distinct(.keep_all = TRUE) %>%
  dplyr::summarise(across(c("Sp_Unbiased_CV","Duration.y","ChullArea.km",
                     "Growth.Habit.short", "GENUS.tree",
                     "SPECIES.tree","TreeName"), first))

spRange <- (speciesMax[,c(2,5:11,16:20)]) - (speciesMin[,c(2,5:11,16:20)])

colnames(speciesMean)[2:18] <- paste(colnames(speciesMean[2:18]),'_mean',sep = '')
colnames(spRange) <- paste(colnames(spRange),'_range',sep = '')
colnames(speciesMax) <- paste(colnames(speciesMax),'_Max',sep = '')
colnames(speciesMin) <- paste(colnames(speciesMin),'_Min',sep = '')

speciesAgg <- cbind(speciesMean, spRange, speciesCat, speciesMax[,c(3,4)],speciesMin[,c(3,4)])

colnames(speciesAgg)
speciesAgg2 <- speciesAgg[ , order(names(speciesAgg))]
colnames(speciesAgg2)

speciesAgg2 <- speciesAgg2[, c(23,37,35,21,22,16,18,1:15,19,20,25:34,38:43)]

colnames(speciesAgg2)

speciesAgg3 <- na.omit(speciesAgg2)
corrs <- cor(speciesAgg3[c(3,8:40)])
CVcorrs <- corrs[1,]
CVcorrs <- sort(CVcorrs, decreasing = TRUE)
CVcorrs
#drop chull area due to high correlation & VIF later on

########################### Feature selection seed mass variation ########################### 

# Perform Boruta feature selection
boruta_output <- Boruta(log(speciesAgg3$Sp_Unbiased_CV) ~ ., data=(speciesAgg3[c(6,7,8:11,13:40)]), pValue = 0.0001, doTrace=0)  

imps <- attStats(boruta_output)
imps2 = imps[imps$decision == 'Confirmed', c('meanImp', 'decision')]
imps2

Sp_vars <- rownames(imps2)
length(Sp_vars)

speciesAgg2_scale <- speciesAgg3
speciesAgg2_scale[c(8:40)] <- lapply(speciesAgg2_scale[c(8:40)], function(x) c(scale(x)))

########################### Linear Model ########################### 

# Create a linear model formula using the scaled variables as predictors and Sp_Unbiased_CV as the response
lm1 <- paste(paste("log(Sp_Unbiased_CV) ~"), paste(paste(Sp_vars, collapse="+"), sep='+'))
lm1test <- lm(as.formula(lm1), data = speciesAgg2_scale)
vif(lm1test)
summary(lm1test)
dplyr::select(speciesAgg2_scale, all_of(Sp_vars)) %>%
  cor()

# Fit a linear model using the formula
CVglm1 <- glm(as.formula(lm1), data = speciesAgg2_scale)

summary(CVglm1)
plot(CVglm1)
##

print(as.formula(lm1))

########################### Phylogenetic GLS ########################### 
fitplain <- gls(log(Sp_Unbiased_CV) ~1, data=speciesAgg2_scale,
            correlation=corPagel(1, subset_tree, form = ~Name2))
lambda_value <- coef(fitplain$modelStruct$corStruct, unconstrained=FALSE)
print(lambda_value)

fitrange <- gls(log(mg.seed_range) ~1, data=speciesAgg2_scale,
                correlation=corPagel(1, subset_tree, form = ~Name2))
lambda_value <- coef(fitrange$modelStruct$corStruct, unconstrained=FALSE)
print(lambda_value)

#### using phylo sig
meanseed_vector <- setNames(log(speciesAgg2_scale$mg.seed_mean), speciesAgg2_scale$Name2)
seedCV_vector <- setNames(log(speciesAgg2_scale$Sp_Unbiased_CV), speciesAgg2_scale$Name2)

phylosig_resultmean <- phylosig(subset_tree, meanseed_vector, method = "lambda", test = TRUE)
print(phylosig_resultmean)

phylosig_resultcv <- phylosig(subset_tree, seedCV_vector, method = "lambda", test = TRUE)
print(phylosig_resultcv)

# Fit Generalized Least Squares (PGLS) model
fit3 <- gls(as.formula(lm1), data=speciesAgg2_scale,
            correlation=corPagel(1, subset_tree, form = ~Name2))
fitnophylo <- gls(as.formula(lm1), data=speciesAgg2_scale)


LL_no_phylo <- logLik(fitnophylo)
LL_phylo <- logLik(fit3)

pseudo_R2 <- 1 - (as.numeric(LL_no_phylo) / as.numeric(LL_phylo))
pseudo_R2

AIC(fit3)
AIC(fitnophylo)

# adding growth form and duration just to check
fit4 <- gls(as.formula(log(Sp_Unbiased_CV) ~ coords.x2_Max + Ydays_range + yearof_mintemp_range
                       +Duration.y+Growth.Habit.short), data=speciesAgg2_scale,
            correlation=corPagel(1, subset_tree, form = ~Name2))
summary(fit4)
fit5 <- gls(as.formula(Sp_Unbiased_CV ~ Duration.y+Growth.Habit.short), data=speciesAgg2_scale,
            correlation=corPagel(1, subset_tree, form = ~Name2))
summary(fit5)
#no significant differences among growth form or duration

# Extract and display model results
results <- summary(fit3)
print(results, digits = 3)
confPGLS <- confint(fit3)
confPGLS <- round(confPGLS, digits = 3)

GLScoeffs <- results$coefficients
GLScoeffs <- as.data.frame(GLScoeffs)
GLScoeffs$variable <- rownames(GLScoeffs)
GLScoeffs$stdDev <- sqrt(diag(vcov(fit3)))
GLScoeffs$Pval <- coef(summary(fit3))[,4]
GLScoeffs$Tval <- coef(summary(fit3))[,3]

print(GLScoeffs, digits = 2)
round(confint(fit3), digits = 3)
round(coef(summary(fit3))[,3], digits = 3)

GLScoeffs$significance <- ifelse(GLScoeffs$Pval < 0.0001, "***",
                                     ifelse(GLScoeffs$Pval < 0.001, "**", 
                                            ifelse(GLScoeffs$Pval < 0.05, "*", "")))

GLScoeffs$variable
GLScoeffs$Labs <- c('(Intercept)','Maximum latitude',
                    'Range of collection day',
                    'Range of annual minimum temps')

GLScoeffs2 <- results$tTable

GLScoeffs2 <- round(GLScoeffs2, digits = 3)
confGLS <- round(confint(fit3), digits = 3)

########################### Plots ########################### 

# Plot coefficients with standard errors
coul <- brewer.pal(1, "YlOrRd") 

plot <- ggplot(subset(GLScoeffs, variable != '(Intercept)'), 
       aes(y=reorder(Labs, GLScoeffs), 
           x=GLScoeffs)) + 
  geom_point(size = 3, color="#f05e20") +  # Points for coefficient estimates
  geom_errorbar(aes(xmin=GLScoeffs-stdDev, 
                    xmax=GLScoeffs+stdDev), width=.5, 
                position=position_dodge(.9),
                color="#f05e20", linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +  # Add horizontal line at y = 0
  #geom_text(aes(label = significance), vjust = 0.3, size = 20) + 
  ylab(NULL) + xlab('Coefficient (scaled)')+
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 60),
        axis.title = element_text(size = 70),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())

########################### Mean seed Feature Selection ########################### 

# Perform Boruta feature selection
boruta_output <- Boruta(log(speciesAgg3$mg.seed_mean) ~ ., data=(speciesAgg3[c(6,7,8:11,13:40)]), pValue = 0.0001, doTrace=0)  

imps <- attStats(boruta_output)
imps2 = imps[imps$decision == 'Confirmed', c('meanImp', 'decision')]
imps2

Sp_vars <- rownames(imps2)
length(Sp_vars)

speciesAgg2_scale <- speciesAgg3
speciesAgg2_scale[c(8:40)] <- lapply(speciesAgg2_scale[c(8:40)], function(x) c(scale(x)))

########################### Linear Model ########################### 

# Create a linear model formula using the scaled variables as predictors and Sp_Unbiased_CV as the response
lm1 <- paste(paste("log(mg.seed_mean) ~"), paste(paste(Sp_vars, collapse="+"), sep='+'))
lm1test <- lm(as.formula(lm1), data = speciesAgg2_scale)
vif(lm1test)
summary(lm1test)
dplyr::select(speciesAgg2_scale, all_of(Sp_vars)) %>%
  cor()

# Fit a linear model using the formula
CVglm1 <- glm(as.formula(lm1), data = speciesAgg2_scale)

summary(CVglm1)
plot(CVglm1)
##

print(as.formula(lm1))

########################### Mean seed Phylogenetic GLS ########################### 

# Fit Generalized Least Squares (PGLS) model
fit3 <- gls(as.formula(lm1), data=speciesAgg2_scale,
            correlation=corPagel(1, subset_tree, form = ~Name2))
fitnophylo <- gls(as.formula(lm1), data=speciesAgg2_scale)

R2(fit3, fitnophylo)
LL_no_phylo <- logLik(fitnophylo)
LL_phylo <- logLik(fit3)

pseudo_R2 <- 1 - (as.numeric(LL_no_phylo) / as.numeric(LL_phylo))
pseudo_R2

AIC(fit3)
AIC(fitnophylo)

# Extract and display model results
results <- summary(fit3)
print(results, digits = 3)
confPGLS <- confint(fit3)
confPGLS <- round(confPGLS, digits = 3)

GLScoeffs <- results$coefficients
GLScoeffs <- as.data.frame(GLScoeffs)
GLScoeffs$variable <- rownames(GLScoeffs)
GLScoeffs$stdDev <- sqrt(diag(vcov(fit3)))
GLScoeffs$Pval <- coef(summary(fit3))[,4]
GLScoeffs$Tval <- coef(summary(fit3))[,3]

print(GLScoeffs, digits = 2)
round(confint(fit3), digits = 3)
round(coef(summary(fit3))[,3], digits = 3)

GLScoeffs$significance <- ifelse(GLScoeffs$Pval < 0.0001, "***",
                                 ifelse(GLScoeffs$Pval < 0.001, "**", 
                                        ifelse(GLScoeffs$Pval < 0.05, "*", "")))

GLScoeffs$variable
GLScoeffs$Labs <- c('(Intercept)','Maximum latitude',
                    'Range of collection day',
                    'Range of annual minimum temps')

GLScoeffs2 <- results$tTable

GLScoeffs2 <- round(GLScoeffs2, digits = 3)
confGLS <- round(confint(fit3), digits = 3)

########################### Plots ########################### 

# Plot coefficients with standard errors
coul <- brewer.pal(1, "YlOrRd") 

plot2 <- ggplot(subset(GLScoeffs, variable != '(Intercept)'), 
               aes(y=reorder(variable, GLScoeffs), 
                   x=GLScoeffs)) + 
  geom_point(size = 3, color="#f05e20") +  # Points for coefficient estimates
  geom_errorbar(aes(xmin=GLScoeffs-stdDev, 
                    xmax=GLScoeffs+stdDev), width=.5, 
                position=position_dodge(.9),
                color="#f05e20", linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +  # Add horizontal line at y = 0
  #geom_text(aes(label = significance), vjust = 0.3, size = 20) + 
  ylab(NULL) + xlab('Coefficient (scaled)')+
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 60),
        axis.title = element_text(size = 70),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())
