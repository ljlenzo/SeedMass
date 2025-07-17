########################### Q3 Individual Species Models ########################### 
getwd()
set.seed(123)

library(Boruta)
library(dplyr)
library(tidyr)
library(stringr)
library(plyr)
library(spdep)
library(car)
library(sp)
library(ggplot2)

########################### Jitter duplicate points ########################### 

# k nearest neighbor doesn't work on duplicate points so jitter slightly
coords <- data.frame(lon = species2$coords.x1, lat = species2$coords.x2)

coords_dup <- coords[duplicated(coords[,c('lon', 'lat')]),] 
coords_distinct <- distinct(coords_dup)

fin_list <- list()
for (i in 1:nrow(coords_distinct)) {
  temp_data <- species2 %>% filter(coords.x1 == coords_distinct[i,1] & coords.x2 == coords_distinct[i,2])
  
  temp_data$coords.x1 <- temp_data$coords.x1 + seq(0, 0.0001, 0.00001)[1:nrow(temp_data)]
  temp_data$coords.x2 <- temp_data$coords.x2 + seq(0, 0.0001, 0.00001)[1:nrow(temp_data)]
  
  fin_list[[i]] <- temp_data
  print(i)
}

fin_bind <- data.table::rbindlist(fin_list)

numbers <- c(fin_bind$ACC_NUM)
not_dupe <- species2 %>% filter(!ACC_NUM %in% numbers)

modified_data <- rbind(fin_bind, not_dupe)
species2 <- modified_data

###########################  create individual dataframes ########################### 
for (i in 1:length(spnames)) {
  # first, bring in the dataframe 
  Sp <- spnames[[i]]
  SpeciesDF <- subset(species2, species2$Name2 == paste(Sp))
  assign(paste(SpeciesDF$GENUS, SpeciesDF$SPECIES, sep = '_'), SpeciesDF) # rename
  print(i)
}

data.frame(colnames(species2)) 

###########################  scale data ########################### 
# center day and year (rather than fully scale)

summary(species2[,c(10,  17:19,21, 24:27, 29:31,33:35)])

for (i in 1:length(spnames)) {
  # first, bring in the dataframe 
  Sp <- get(paste(spnames[i],sep=''))
  # create a duplicate that we can scale
  Sp_scale <- Sp
  # scale the new data
  Sp_scale[,c(10,  17:19,21, 24:27, 29:31,33:35)] <- lapply(Sp_scale[,c(10,  17:19,21, 24:27, 29:31,33:35)], 
                                                                                     function(x) c(scale(x)))
  
  Sp_scale[,c(14,23)] <- lapply(Sp_scale[,c(14,23)], 
                               function(x) c(scale(x, center = TRUE, scale = FALSE)))
  
  
  # rename these columns with 'scale'
  colnames(Sp_scale) <- paste(colnames(Sp_scale),'_scale',sep = '')
  # combine the scaled data with our original data
  temp <- cbind(Sp, Sp_scale[,c(10,  17:19,21, 24:27, 29:31,33:35,14,23)])
  # Name the new dataframe based on our sp name
  assign(paste(spnames[i],'4',sep=''),temp)
  
  print(i)
}
rm(temp, Sp, Sp_scale)

########################### Run Individual Models ########################### 
# select the column names of interest
EnvVars <- colnames(`Sporobolus_airoides4`[,c(14, 17:19,21, 23:27, 29:31,33,34,35)])

# create an empty matrix to store the results
mymatrix <- matrix(ncol=44, nrow = length(spnames))
spglm <- as.data.frame(mymatrix)

# add column names to the result data frame
colnames(spglm) <- c('Species','Sample','formula', '(Intercept)','NVar','NrSqSplm','MeanMg','MedianMg','CV','MaxMg','MinMg','Pvalue',
                     paste('coefficient_',colnames(`Sporobolus_airoides4`[,c(14, 17:19,21, 23:27, 29:31,33,34,35)]),sept=''),
                     paste('pvalues_',colnames(`Sporobolus_airoides4`[,c(14, 17:19,21, 23:27, 29:31,33,34,35)]),sep=''))

# create an empty matrix to store the coefficients and p-values for each species
CoeffsSp <- matrix(ncol=14, nrow=1000)
colnames(CoeffsSp) <- c('Species','Formula','variables','coefficient','Std.Error',
                        'zvalue','pvalues','Rsquared','NVar','FAMILY','GrowthForm',
                        'SampleSize','ConfL','ConfU')

# Begin looping through each species to determine important variables and run a spatial GLM on them
for (j in 1:length(spnames)){

SpStart <- get(paste(spnames[j], sep='')) # get the unscaled dataframe
rownames(SpStart) <- NULL
SpScale <- get(paste(spnames[j],'4',sep='')) # get scaled dataframe
rownames(SpScale) <- NULL

n = length(SpStart[,c(14, 17:19,21, 23:27, 29:31,33,34,35)]) + 11

mymatrix <- matrix(ncol=n, nrow = 1) # make an empty matrix
df <- as.data.frame(mymatrix)
colnames(df) <- c('Species','Sample','formula', '(Intercept)','NVar','NrSqSplm','MeanMg','MedianMg','CV','MaxMg','MinMg',
                  paste(colnames(SpStart[,c(14, 17:19,21, 23:27, 29:31,33,34,35)]),'_scale',sep = ''))
assign(paste(SpStart$Name2[1],'df',sep = '_'),df) # rename for sp

print(j)

    boruta_output <- Boruta(SpStart$mg.seed ~ ., data=(SpStart[,c(14, 17:19,21, 23:27, 29:31,33,34,35)]), pValue = 0.0001, doTrace=0)  
    
    imps <- attStats(boruta_output)
    imps2 = imps[imps$decision == 'Confirmed', c('meanImp', 'decision')]
    
    Sp_vars <- rownames(imps2)
    varlength <- length(Sp_vars)
    
    # If the variable length is less than 1, perform a rough fix 
    if (varlength < 1) {
      
      # Determine whether tentative variables are confirmed or not
      roughFixMod <- TentativeRoughFix(boruta_output)
      
      # Get the significant attributes 
      boruta_signif2 <- getSelectedAttributes(roughFixMod)
      
      # Get the attribute stats
      imps <- attStats(roughFixMod)
      
      # Only get confirmed variables and pull mean importance and decision
      imps2 = imps[imps$decision == 'Confirmed', c('meanImp', 'decision')]
      
      # Get the names of the significant variables
      Sp_vars <- rownames(imps2)
      
      # Count up variables
      varlength <- length(Sp_vars)
    }
    # If the variable length is greater than or equal to 1
    else{
      # Create a vector of the scaled variables
      Sp_envarsScale <- paste(Sp_vars,'_scale',sep = '')
      
      # Print the total number of variables
      print('Total Variables')
      print(varlength)
    }
    
    # If the variable length is greater than or equal to 1
    if (varlength >= 1) {
      
      # Create a vector of the scaled variables and response
      Sp_vars2 <- c('ACC_NUM','mg.seed',Sp_envarsScale)
      
      # Filter the SpScale data frame to only include the important variables
      SpScale2 <- subset(SpScale, select=Sp_vars2)
      
      # Create a linear regression formula with the scaled variable names
      lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(Sp_envarsScale, collapse="+"),  sep='+'))
      
      # Fit a linear regression model with the formula and scaled data
      Sp_fulllm <- lm(as.formula(lm1),data=SpScale)
      
      # Assign the species name to the linear regression model
      assign(paste(spnames[j],'lm',sep='_'),Sp_fulllm)
      
      # Get the latitude and longitude 
      Sp_LL <- SpScale[,c(12,13)]
      Sp_LL <- coordinates(Sp_LL)
      
      test<-knn2nb(knearneigh(Sp_LL))
      print(SpStart$Name2[1])
      test
      
      cards <- card(test)
      col.w <- nb2listw(test)
      
      lmResiduals <- rep(0, length(SpScale$mg.seed))
      resIndex <- Sp_fulllm$residuals %>% names() %>% as.integer();
      lmResiduals[resIndex] <- Sp_fulllm$residuals
      
        col.w %>%
        # Calculate the Moran's I test statistic for spatial autocorrelation of residuals
        spdep::moran.test(lmResiduals, ., zero.policy = TRUE)
      
      serrRslt <- spatialreg::errorsarlm(as.formula(lm1),
                                         data = SpScale,
                                         listw = col.w,
                                         zero.policy = TRUE, 
                                         na.action = na.omit);
      assign(paste(spnames[j],'SAR',sep='_'),serrRslt)
      
      # Plot the histogram and Q-Q plot of residuals
      hist(residuals(serrRslt))
      qqnorm(residuals(serrRslt)); qqline(residuals(serrRslt))
      
      # Generate component-residual plots
      #crPlots(Sp_fulllm)
      
      # Create a data frame with model coefficients, p-values, and other metrics
      Spdf1 <- as.data.frame(serrRslt$coefficients)
      Spdf1$Species <- SpStart$Name2[1]
      Spdf1$variables <- rownames(Spdf1)
      rownames(Spdf1) <- NULL
      colnames(Spdf1)[1] = 'coefficient'
      
      testtable <- summary(serrRslt)
      testtable2 <- testtable$Coef
      tablee <- as.data.frame(testtable2)
      
      confsp <- as.data.frame(confint(serrRslt))
      confsp <- confsp[-c(1),]

      # Add test stats to the data frame
      Spdf1$pvalues <- tablee$`Pr(>|z|)`
      Spdf1$Std.Error <- tablee$`Std. Error`
      Spdf1$zvalue <- tablee$`z value`
      Spdf1$SampleSize <- nrow(SpScale)
      
      # Add number of variables used in the model minus the intercept
      Spdf1$NVar <- nrow(Spdf1) - 1
      
      test <- summary(serrRslt, Nagelkerke = TRUE)
      
      # Add Nagelkerke's R-squared to the data frame
      Spdf1$NrSqSplm <- test$NK
      
      #cat(SpStart$Name2[1], 'R squared',Spdf1$NrSqSplm[1])
      
      # Add mean, median, coefficient of variation, maximum, and minimum of the response variable
      Spdf1$MeanMg <- mean(SpScale$mg.seed)
      Spdf1$MedianMg <- median(SpScale$mg.seed)
      Spdf1$MaxMg <- max(SpScale$mg.seed)
      Spdf1$MinMg <- min(SpScale$mg.seed)
      Spdf1$ConfL <- confsp$`2.5 %`
      Spdf1$ConfU <- confsp$`97.5 %`
      # Pivot Spdf1 from long to wide format
      pivwide <- pivot_wider(Spdf1, names_from = variables, values_from = c(coefficient, pvalues)) 
      
      # Add metadata columns to the pivoted data
      pivwide$Sample <- j
      pivwide$formula <- paste(lm1)
      pivwide$Species <- SpStart$Name2[1]
      
      # Combine the pivoted data with the existing data for the species
      tempDF <- rbind.fill(pivwide, get(paste(SpStart$Name2[1],'df',sep = '_')))
      tempDF <- subset(tempDF, !is.na(tempDF$Species))
      
      # Rename the data for use with spatial packages
      assign(paste(SpStart$Name2[1],'df',sep = '_'),tempDF)
      
      # Add the combined data to the overall spatial data frame
      spglm <- rbind.fill(tempDF, spglm)
      
      Spdf1$Formula <- paste(lm1)
      Spdf1$Rsquared <- test$NK
      
      # Get the common columns between Spdf1 and CoeffsSp, and combine them
      common_cols <- intersect(colnames(Spdf1), colnames(CoeffsSp))
      CoeffsSp <- rbind(
        subset(Spdf1, select = common_cols), 
        subset(CoeffsSp, select = common_cols)
        
    )
    
    }
    
    #if there are no important variables at all print species name and "no variables"
    else {
      print(SpStart$Name2)
      print("no vars")
      next
    }
}

########################### Plots ########################### 

##### change in fall max examples ########
  
# Fall max – Erigeron_pumilus + forb 
# Hymenoclea_salsola - shrub

termsE <- names(Erigeron_pumilus_SAR$coefficients)
termsE <- termsE[-c(1,3)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsE, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Erigeron_pumilus[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noFallMax <- spatialreg::errorsarlm(as.formula(lm1),
                                 data = Erigeron_pumilus4,
                                 listw = col.w,
                                 zero.policy = TRUE, 
                                 na.action = na.omit)

noFallMaxRes <- noFallMax$residuals
plot_data <- data.frame(DeltFallMax_scale = Erigeron_pumilus4$DeltFallMax_scale, 
                        mg.seed = noFallMaxRes)
residuals_model <- lm(noFallMaxRes ~ DeltFallMax_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$DeltFallMax_original <- predicted_seed_mass$DeltFallMax_scale * sd(Erigeron_pumilus4$DeltFallMax) + mean(Erigeron_pumilus4$DeltFallMax)
predicted_seed_mass$DeltFallMax_original <- Erigeron_pumilus4$DeltFallMax

predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = DeltFallMax_original, y = mg.seedOG)) +
  geom_point(color = '#6F8FC9',size = 10,alpha=0.5) +
  labs(x = "Change in fall max temps", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50, colour = '#808080'),
        axis.title = element_text(size = 40, colour = "#808080", face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())

## atpo

termsH <- names(Hymenoclea_salsola_SAR$coefficients)
termsH <- termsH[-c(1,4)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsH, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Hymenoclea_salsola[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noFallMax <- spatialreg::errorsarlm(as.formula(lm1),
                                    data = Hymenoclea_salsola4,
                                    listw = col.w,
                                    zero.policy = TRUE, 
                                    na.action = na.omit)

noFallMaxRes <- noFallMax$residuals
plot_data <- data.frame(DeltFallMax_scale = Hymenoclea_salsola4$DeltFallMax_scale, 
                        mg.seed = noFallMaxRes)
residuals_model <- lm(noFallMaxRes ~ DeltFallMax_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$DeltFallMax_original <- predicted_seed_mass$DeltFallMax_scale * sd(Hymenoclea_salsola4$DeltFallMax) + mean(Hymenoclea_salsola4$DeltFallMax)
predicted_seed_mass$DeltFallMax_original <- Hymenoclea_salsola4$DeltFallMax

predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = DeltFallMax_original, y = mg.seedOG)) +
  geom_point(color = '#AD964D',size = 10,alpha=0.5) +
  labs(x = "Change in fall max temps", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40, colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())

##### change in precip precip #####

# Artemisia tridentata - delt precip
# Hesperostipa comata + delt precip

## HECO

termsH <- names(Hesperostipa_comata_SAR$coefficients)
termsH <- termsH[-c(1,4)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsH, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Hesperostipa_comata[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noDeltPrecip <- spatialreg::errorsarlm(as.formula(lm1),
                                       data = Hesperostipa_comata4,
                                       listw = col.w,
                                       zero.policy = TRUE, 
                                       na.action = na.omit)

noDPrecRes <- noDeltPrecip$residuals
plot_data <- data.frame(DeltPrecip_scale = Hesperostipa_comata4$DeltPrecip_scale, 
                        mg.seed = noDPrecRes)
residuals_model <- lm(noDPrecRes ~ DeltPrecip_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$DeltPrecip_original <- predicted_seed_mass$DeltPrecip_scale * sd(Hesperostipa_comata4$DeltPrecip) + mean(Hesperostipa_comata4$DeltPrecip)
predicted_seed_mass$DeltPrecip_original <- Hesperostipa_comata4$DeltPrecip
predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = DeltPrecip_original, y = mg.seedOG)) +
  geom_point(color = '#60A75E',size = 10,alpha=0.5) +
  labs(x = "Change in precip", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40,colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())

#### ARTR

termsA <- names(Artemisia_tridentata_SAR$coefficients)
termsA <- termsA[-c(1,6)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsA, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Artemisia_tridentata[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noDeltPrecip <- spatialreg::errorsarlm(as.formula(lm1),
                                       data = Artemisia_tridentata4,
                                       listw = col.w,
                                       zero.policy = TRUE, 
                                       na.action = na.omit)

noDPrecRes <- noDeltPrecip$residuals
plot_data <- data.frame(DeltPrecip_scale = Artemisia_tridentata4$DeltPrecip_scale, 
                        mg.seed = noDPrecRes)
residuals_model <- lm(noDPrecRes ~ DeltPrecip_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$DeltPrecip_original <- predicted_seed_mass$DeltPrecip_scale * sd(Artemisia_tridentata4$DeltPrecip) + mean(Artemisia_tridentata4$DeltPrecip)
predicted_seed_mass$DeltPrecip_original <- Artemisia_tridentata4$DeltPrecip

predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = DeltPrecip_original, y = mg.seedOG)) +
  geom_point(color = '#AD964D',size = 10,alpha=0.5) +
  labs(x = "Change in precip", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40,colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


##### change in spring max #####

# Pseudoroegneria_spicata  - spring max grass
# Eriogonum_heracleoides + spring max forb

## PSSP

termsP <- names(Pseudoroegneria_spicata_SAR$coefficients)
termsP <- termsP[-c(1,2)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsP, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Pseudoroegneria_spicata[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noSpringMax <- spatialreg::errorsarlm(as.formula(lm1),
                                       data = Pseudoroegneria_spicata4,
                                       listw = col.w,
                                       zero.policy = TRUE, 
                                       na.action = na.omit)

noSpringMaxRes <- noSpringMax$residuals
plot_data <- data.frame(DeltSpringMax_scale = Pseudoroegneria_spicata4$DeltSpringMax_scale, 
                        mg.seed = noSpringMaxRes)
residuals_model <- lm(noSpringMaxRes ~ DeltSpringMax_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$SpringMax_original <- predicted_seed_mass$DeltSpringMax_scale * 
#  sd(Pseudoroegneria_spicata4$DeltSpringMax) + 
#  mean(Pseudoroegneria_spicata4$DeltSpringMax)
predicted_seed_mass$SpringMax_original <- Pseudoroegneria_spicata4$DeltSpringMax

predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = SpringMax_original, y = mg.seedOG)) +
  geom_point(color = '#60A75E',size = 10,alpha=0.5) +
  labs(x = "Change in spring max temps", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40, colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


#### ERHE

termsE <- names(Eriogonum_heracleoides_SAR$coefficients)
termsE <- termsE[-c(1,4)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsE, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Eriogonum_heracleoides[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noSpringMax <- spatialreg::errorsarlm(as.formula(lm1),
                                       data = Eriogonum_heracleoides4,
                                       listw = col.w,
                                       zero.policy = TRUE, 
                                       na.action = na.omit)

noSpringMaxRes <- noSpringMax$residuals
plot_data <- data.frame(DeltSpringMax_scale = Eriogonum_heracleoides4$DeltSpringMax_scale, 
                        mg.seed = noSpringMaxRes)
residuals_model <- lm(noSpringMaxRes ~ DeltSpringMax_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$SpringMax_original <- predicted_seed_mass$DeltSpringMax_scale * 
#  sd(Eriogonum_heracleoides4$DeltSpringMax) + 
#  mean(Eriogonum_heracleoides4$DeltSpringMax)
predicted_seed_mass$SpringMax_original <- Eriogonum_heracleoides4$DeltSpringMax
predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = SpringMax_original, y = mg.seedOG)) +
  geom_point(color = '#6F8FC9',size = 10,alpha=0.5) +
  labs(x = "Change in spring max temps", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50, colour = '#808080'),
        axis.title = element_text(size = 40, colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


##### species example plots for soil pH #######
#Soil pH: BAHO +  elel5 -

## ELEL
exp(-0.140986814*sd(Elymus_elymoides$SoilpH))
termsE <- names(Elymus_elymoides_SAR$coefficients)
termsE <- termsE[-c(1,5)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsE, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Elymus_elymoides[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noSoil <- spatialreg::errorsarlm(as.formula(lm1),
                                 data = Elymus_elymoides4,
                                 listw = col.w,
                                 zero.policy = TRUE, 
                                 na.action = na.omit)

noSoilRes <- noSoil$residuals
plot_data <- data.frame(SoilpH_scale = Elymus_elymoides4$SoilpH_scale, 
                        mg.seed = noSoilRes)
residuals_model <- lm(noSoilRes ~ SoilpH_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$SoilpH_original <- predicted_seed_mass$SoilpH_scale * sd(Elymus_elymoides4$SoilpH) + mean(Elymus_elymoides4$SoilpH)
predicted_seed_mass$SoilpH_original <- Elymus_elymoides4$SoilpH
predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = SoilpH_original, y = mg.seedOG)) +
  geom_point(color = '#60A75E',size = 10,alpha=0.5) +
  labs(x = "Soil pH", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40, colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


#BAHO
exp(0.10917082*sd(Elymus_elymoides$SoilpH))
termsB <- names(Balsamorhiza_hookeri_SAR$coefficients)
termsB <- termsB[-c(1,3)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsB, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Balsamorhiza_hookeri[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noSoil <- spatialreg::errorsarlm(as.formula(lm1),
                                 data = Balsamorhiza_hookeri4,
                                 listw = col.w,
                                 zero.policy = TRUE, 
                                 na.action = na.omit)

noSoilRes <- noSoil$residuals
plot_data <- data.frame(SoilpH_scale = Balsamorhiza_hookeri4$SoilpH_scale, 
                        mg.seed = noSoilRes)
residuals_model <- lm(noSoilRes ~ SoilpH_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$SoilpH_original <- predicted_seed_mass$SoilpH_scale * sd(Balsamorhiza_hookeri4$SoilpH) + mean(Balsamorhiza_hookeri4$SoilpH)
predicted_seed_mass$SoilpH_original <- Balsamorhiza_hookeri4$SoilpH
predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = SoilpH_original, y = mg.seedOG)) +
  geom_point(color = '#6F8FC9',size = 10,alpha=0.5) +
  labs(x = "Soil pH", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40,colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


##### ann min temp #########
#  Ann min temp:  Fallugia paradoxa or plantago –

# Asclepias_speciosa + forb
# Plantago_ovata - shrub

termsA <- names(Asclepias_speciosa_SAR$coefficients)
termsA <- termsA[-c(1,3)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsA, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Asclepias_speciosa[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noMinT <- spatialreg::errorsarlm(as.formula(lm1),
                                       data = Asclepias_speciosa4,
                                       listw = col.w,
                                       zero.policy = TRUE, 
                                       na.action = na.omit)

noMinTRes <- noMinT$residuals
plot_data <- data.frame(AnnMinTemp_scale = Asclepias_speciosa4$AnnMinTemp_scale, 
                        mg.seed = noMinTRes)
residuals_model <- lm(noMinTRes ~ AnnMinTemp_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$AnnMinTemp_original <- predicted_seed_mass$AnnMinTemp_scale * sd(Asclepias_speciosa4$AnnMinTemp) + mean(Asclepias_speciosa4$AnnMinTemp)
predicted_seed_mass$AnnMinTemp_original <- Asclepias_speciosa4$AnnMinTemp
predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = AnnMinTemp_original, y = mg.seedOG)) +
  geom_point(color = '#6F8FC9',size = 10,alpha=0.5) +
  labs(x = "Annual minimum temperature", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40, colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


### PLOV

termsP <- names(Plantago_ovata_SAR$coefficients)
termsP <- termsP[-c(1,6)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsP, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Plantago_ovata[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noMinT <- spatialreg::errorsarlm(as.formula(lm1),
                                 data = Plantago_ovata4,
                                 listw = col.w,
                                 zero.policy = TRUE, 
                                 na.action = na.omit)

noMinTRes <- noMinT$residuals
plot_data <- data.frame(AnnMinTemp_scale = Plantago_ovata4$AnnMinTemp_scale, 
                        mg.seed = noMinTRes)
residuals_model <- lm(noMinTRes ~ AnnMinTemp_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$AnnMinTemp_original <- predicted_seed_mass$AnnMinTemp_scale * sd(Plantago_ovata4$AnnMinTemp) + mean(Plantago_ovata4$AnnMinTemp)
predicted_seed_mass$AnnMinTemp_original <- Plantago_ovata4$AnnMinTemp
predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = AnnMinTemp_original, y = mg.seedOG)) +
  geom_point(color = '#AD964D',size = 10,alpha=0.5) +
  labs(x = "Annual minimum temperature", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40, colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


##### 6 month precip #########

# Heterotheca_villosa + forb
# Achnatherum_thurberianum - grass

# HEVI

termsH <- names(Heterotheca_villosa_SAR$coefficients)
termsH <- termsH[-c(1,3)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsH, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Heterotheca_villosa[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noSixMonthPrecip <- spatialreg::errorsarlm(as.formula(lm1),
                                 data = Heterotheca_villosa4,
                                 listw = col.w,
                                 zero.policy = TRUE, 
                                 na.action = na.omit)

noSixMonthPrecipRes <- noSixMonthPrecip$residuals
plot_data <- data.frame(SixMonthPrecip_scale = Heterotheca_villosa4$SixMonthPrecip_scale, 
                        mg.seed = noSixMonthPrecipRes)
residuals_model <- lm(noSixMonthPrecipRes ~ SixMonthPrecip_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$SixMonthPrecip_original <- predicted_seed_mass$SixMonthPrecip_scale * sd(Heterotheca_villosa4$SixMonthPrecip) + mean(Heterotheca_villosa4$SixMonthPrecip)
predicted_seed_mass$SixMonthPrecip_original <- Heterotheca_villosa$SixMonthPrecip

predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = SixMonthPrecip_original, y = mg.seedOG)) +
  geom_point(color = '#6F8FC9',size = 10,alpha=0.5) +
  labs(x = "Six month precip", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40, colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


### ACTH

termsA <- names(Achnatherum_thurberianum_SAR$coefficients)
termsA <- termsA[-c(1,3)]
lm1 <- paste(paste("log(mg.seed) ~"), paste(paste(termsA, collapse="+"),  sep='+'))
Sp_LL <- coordinates(Achnatherum_thurberianum[,c(12,13)])
col.w<-nb2listw(knn2nb(knearneigh(Sp_LL)))
noSixMonthPrecip <- spatialreg::errorsarlm(as.formula(lm1),
                                           data = Achnatherum_thurberianum4,
                                           listw = col.w,
                                           zero.policy = TRUE, 
                                           na.action = na.omit)

noSixMonthPrecipRes <- noSixMonthPrecip$residuals
plot_data <- data.frame(SixMonthPrecip_scale = Achnatherum_thurberianum4$SixMonthPrecip_scale, 
                        mg.seed = noSixMonthPrecipRes)
residuals_model <- lm(noSixMonthPrecipRes ~ SixMonthPrecip_scale, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predicted_seed_mass <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

#predicted_seed_mass$SixMonthPrecip_original <- predicted_seed_mass$SixMonthPrecip_scale * sd(Achnatherum_thurberianum4$SixMonthPrecip) + mean(Achnatherum_thurberianum4$SixMonthPrecip)
predicted_seed_mass$SixMonthPrecip_original <- Achnatherum_thurberianum4$SixMonthPrecip

predicted_seed_mass$mg.seedOG <- exp(predicted_seed_mass$mg.seed)

plot <- ggplot(predicted_seed_mass, aes(x = SixMonthPrecip_original, y = mg.seedOG)) +
  geom_point(color = '#60A75E',size = 10,alpha=0.5) +
  labs(x = "Six month precip", y = "Seed mass")+
  geom_line(aes(y = .fitted_back), color = "red") +  # Fitted line
  geom_ribbon(aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey") +  # Confidence interval ribbon
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50,colour = '#808080'),
        axis.title = element_text(size = 40, colour = '#808080', face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())

