########################### Min 3 dataset, environmental ###########################

set.seed(123)
library(Boruta)
library(car)
library(spdep)
library(sp)
library(dplyr)
library(randomForest)
library(RColorBrewer)
library(ggplot2)
library(moments)

# z transform seed mass by each species 

species3 <- species3 %>%
  group_by(Name2) %>%
  mutate(mg.seed_z = scale(mg.seed))

########################### Feature Selection with Boruta ###########################

# Run Boruta feature selection
boruta_output <- Boruta(species3$mg.seed_z ~ ., data=na.omit(species3[c(13:28,6,7)]), pValue = 0.0001, doTrace=0)  

# Plot the Boruta output
plot(boruta_output)

# Get the attributes selected by Boruta as significant 
boruta_signif <- getSelectedAttributes(boruta_output, withTentative = FALSE)
boruta_signif

# Get the attribute statistics from the Boruta output
imps <- attStats(boruta_output)

# Filter out the rejected attributes
imps2 <- imps[imps$decision != 'Rejected', c('meanImp', 'decision')]

# Sort the attributes by mean importance
imps2 <- imps2[order(imps2$meanImp),]

# Only keep the attributes with mean importance above 10
imps3 <- subset(imps2, imps2$meanImp >= 10)

# Get the variable names from the attribute statistics
vars <- rownames(imps2)
vars
vars2 <- rownames(imps3)
vars2

# Get scaled variables
vars_scale <- paste(vars2, '_scale', sep = '')

########################### Models ###########################

# simple linear model
lm1 <- paste(paste("mg.seed_z ~"), paste(paste(vars_scale, collapse="+"), sep='+'))

# random forest model from selected features
rf_model <- randomForest(as.formula(lm1), data = species3)

# predict scaled seed mass from random forest model
predicted_values <- predict(rf_model, species3)

# get true seed mass
actual_values <- species3$mg.seed_z

# Calculate R-squared
rsquared <- 1 - sum((actual_values - predicted_values)^2) / sum((actual_values - mean(actual_values))^2)
print(paste("R-squared value:", round(rsquared, 3)))

# Fit a linear model to check vif
lmfull <- lm(as.formula(lm1), data = species3)
summary(lmfull)
vif(lmfull)

# Pull coordinates for nearest neighbor
Sp_LL <- species3[,c(11,12)]
Sp_LL <- coordinates(Sp_LL)

test<-knn2nb(knearneigh(Sp_LL))
test

cards <- card(test)
col.w <- nb2listw(test)

lmResiduals <- rep(0, length(species3$mg.seed_z))
resIndex <- lmfull$residuals %>% names() %>% as.integer();
lmResiduals[resIndex] <- lmfull$residuals

col.w %>%
  # Calculate the Moran's I to check residuals
  spdep::moran.test(lmResiduals, ., zero.policy = TRUE)

# Run spatial simultaneous auto regressive model
serrRslt <- spatialreg::errorsarlm(as.formula(lmfull),
                                   data = species3,
                                   listw = col.w,
                                   zero.policy = TRUE, 
                                   na.action = na.omit);

vars_scale2 <- vars_scale[-c(10)]
lm2 <- paste(paste("mg.seed_z ~"), paste(paste(vars_scale2, collapse="+"), sep='+'))
serrRslt_noyear <- spatialreg::errorsarlm(as.formula(lm2),
                                   data = species3,
                                   listw = col.w,
                                   zero.policy = TRUE, 
                                   na.action = na.omit)

pred2 <- predict(serrRslt)

noYres <- serrRslt_noyear$residuals
plot_data <- data.frame(Year = species3$Year, 
                        mg.seed = noYres, 
                        Growth.Habit.short = species3$Growth.Habit.short)
residuals_model <- lm(noYres ~ Year, data = plot_data)
predictions <- augment(residuals_model, newdata = plot_data, se_fit = TRUE)
predictions <- predictions %>%
  mutate(lower = .fitted - 1.96 * .se.fit,
         upper = .fitted + 1.96 * .se.fit,
         .fitted_back = exp(.fitted),
         lower_back = exp(lower),
         upper_back = exp(upper))

ggplot(predictions, aes(x = Year, y = mg.seed)) +
  geom_point() +  # Scatter plot of residuals
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +  # Confidence interval ribbon
  geom_line(aes(y = .fitted), color = "blue") +  # Fitted line
  labs(x = "Year", y = "Residuals", title = "Semi-Partials for Year") +
  theme_minimal()

habit_colors <- c("Graminoid" = "#60A75E",
                  "Tree" = "#BB608E",
                  "Shrub" = "#AD964D",
                  "Forb/herb" = "#6F8FC9",
                  "Vine" = "#000000")

p <- ggplot(data = species3, aes(x = Year, y = log(mg.seed), color = Growth.Habit.short)) +
  geom_jitter(size = 1, alpha = 0.7,
             position = position_jitter(height = 0.5)) +
  xlab('Year') + ylab('Seed mass (mg)') + 
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 22, face = "bold"),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank()) +  # Remove tick marks on y-axis
  guides(color = guide_legend(override.aes = list(size = 3)))+
  scale_color_manual(name = "Growth Habit", values = habit_colors) +
  scale_x_continuous(expand = c(0.05,0.05)) +  # Remove buffer space around x-axis limits
  scale_y_continuous(expand = c(0, 0)) +  # Remove buffer space around y-axis limits
  coord_fixed(ratio = 1)  # Set aspect ratio to 1 for a square plot

p2 <- p + geom_line(data = predictions, 
          aes(x = Year,y = .fitted_back), color = "red") +  # Fitted line
        geom_ribbon(data = predictions,
          aes(ymin = lower_back, ymax = upper_back), alpha = 0.2, fill = "grey")+   # Confidence interval ribbon
        coord_fixed(ratio = 1)
  
p2

data_range <- range(species3$mg.seed)
log_base <- 10 
log_breaks <- seq(log(data_range[1], base = log_base), 
                  log(data_range[2], base = log_base), 
                  length.out = 5) 
breaks <- round(10^log_breaks, 2)

# check the histogram and Q-Q plot of residuals
hist(residuals(serrRslt))
qqnorm(residuals(serrRslt)); qqline(residuals(serrRslt))

# Store results
results <- summary(serrRslt)

testtable <- results$Coef
tablee <- as.data.frame(testtable)
tablee$variable <- rownames(tablee)

########################### Figure ###########################
coul <- brewer.pal(1, "GnBu") 

tablee <- subset(tablee, tablee$variable != '(Intercept)')
tablee$variable <- gsub("_scale", "", tablee$variable)

tablee$Upper <- tablee$Estimate + tablee$`Std. Error`
tablee$Lower <- tablee$Estimate - tablee$`Std. Error`

tablee$significance <- ifelse(tablee$`Pr(>|z|)` < 0.0001, "***",
                                 ifelse(tablee$`Pr(>|z|)` < 0.001, "**", 
                                        ifelse(tablee$`Pr(>|z|)` < 0.05, "*", "")))


tablee$variable
tablee$Labs <- c('Previous month precip & temp', 'Change in fall max temp',
                 'Change in spring max temp','Slope',
                 'AWC','Change in precip','Soil pH','Heatload',
                 'Day of year','Year','Organic content','Clay content',
                 'Six month cumulative precip','Annual minimum temperature')

plot <- ggplot(tablee, aes(y=reorder(Labs, Estimate), x=Estimate)) + 
  geom_point(size = 10, color="#43A2CA") +  # Points for coefficient estimates
  geom_errorbar(aes(xmin=Lower, xmax=Upper), width=.5, 
                position=position_dodge(.9),
                color="#43A2CA", linewidth = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +  # Add horizontal line at y = 0
  #geom_text(aes(label = significance), vjust = 0.3, size = 20) + 
  xlab('Coefficient (scaled)') + ylab(NULL)+
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 60),
        axis.title = element_text(size = 70),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())


########################### Clade analysis ########################### 


# z transform seed mass by each species 

species3 <- species3 %>%
  group_by(Name2) %>%
  mutate(mg.seed_z = scale(mg.seed))

##### Order, family, genus level analyses #####

# Orders with > 500 populations
OrdersFive <- c('Asterales','Poales','Caryophyllales','Lamiales','Rosales')

table(species4$FAMILY)

# select the column names of interest
EnvVars <- colnames(species4[,c(14:29,7,8)])

# create an empty matrix to store the results
Ordglm <- as.data.frame(matrix(ncol=48, nrow = 5))
Famglm <- as.data.frame(matrix(ncol=48, nrow = 25))
Genglm <- as.data.frame(matrix(ncol=48, nrow = 125))

OrdRsq <- as.data.frame(matrix(ncol=5, nrow = 5))
FamRsq <- as.data.frame(matrix(ncol=5, nrow = 25))
GenRsq <- as.data.frame(matrix(ncol=5, nrow = 125))

# add column names to the result data frame
ColNameA <- c('Species','Sample','formula', '(Intercept)','NVar','NrSqSplm','MeanMg','MedianMg','CV','MaxMg','MinMg','Pvalue',
              paste('coefficient_',colnames(species4[,c(14:29,7,8)]),sept=''),
              paste('pvalues_',colnames(species4[,c(14:29,7,8)]),sep=''))
colnames(Ordglm) <- ColNameA
colnames(Famglm) <- ColNameA
colnames(Genglm) <- ColNameA

ColNameB <- c('Taxa','model','R2lm','R2m','R2c')

colnames(OrdRsq) <- ColNameB
colnames(FamRsq) <- ColNameB
colnames(GenRsq) <- ColNameB

taxaSum <- as.data.frame(matrix(ncol=3, nrow = 1))
colnames(taxaSum) <- c('Taxa','Pops','NSp')

# create an empty matrix to store the coefficients and p-values for each species
CoeffsSp <- matrix(ncol=15, nrow=1)
colnames(CoeffsSp) <- c('Species','Pops','Formula','variables','coefficient','Std.Error',
                        'zvalue','pvalues','Rsquared','NVar','FAMILY','GrowthForm',
                        'SampleSize','ConfL','ConfU')


AICTableOrd <- as.data.frame(matrix(ncol=5, nrow = 1))
colnames(AICTableOrd) <- c('Model','Taxa','AIC','BIC','loglik')

AICTableFam <- as.data.frame(matrix(ncol=5, nrow = 1))
colnames(AICTableFam) <- c('Model','Taxa','AIC','BIC','loglik')

AICTableGen <- as.data.frame(matrix(ncol=5, nrow = 1))
colnames(AICTableGen) <- c('Model','Taxa','AIC','BIC','loglik')

n = 29
for (i in 1:length(OrdersFive)) {
  print(i)
  print(OrdersFive[i])
  # Order level analysis
  
  # Select i thru 5 orders with >500 populations
  OrderA <- subset(species4, species4$order.y == paste(OrdersFive[i]))
  
  # Prep the dataset to scale potential covariates
  Sp_scale <- OrderA
  # scale the new data
  Sp_scale[,c(16:29)] <- lapply(Sp_scale[,c(16:29)], 
                                function(x) c(scale(x)))
  
  Sp_scale[,c(14,15)] <- lapply(Sp_scale[,c(14,15)], 
                                function(x) c(scale(x, center = TRUE, scale = FALSE)))
  
  
  # rename these columns with 'scale'
  colnames(Sp_scale) <- paste(colnames(Sp_scale),'_scale',sep = '')
  # combine the scaled data with our original data
  temp <- cbind(OrderA, Sp_scale[,c(14:29,7,8)])
  
  ## analysis time!
  SpStart <- OrderA
  rownames(SpStart) <- NULL
  SpScale <- temp
  rownames(SpScale) <- NULL
  
  print('running RF')
  boruta_output <- Boruta(SpStart$mg.seed ~ ., data=(SpStart[,c(14:29,7,8)]), pValue = 0.0001, doTrace=0)  
  
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
    lm1 <- paste(paste("(mg.seed) ~ coords.x1 + coords.x2 +"),
                 paste(paste(Sp_envarsScale, collapse="+"),  sep='+'))
    # Fit a linear regression model with the formula and scaled data
    Ord_fulllm <- lm(as.formula(lm1),data=SpScale)
    
    lMm1 <- paste(lm1, '+ (1 | Name2)', sep='')
    
    MixModOrd <- lmer(as.formula(lMm1), 
                      data = SpScale)
    
    # Plot the histogram and Q-Q plot of residuals
    hist(residuals(MixModOrd))
    qqnorm(residuals(MixModOrd)); qqline(residuals(MixModOrd))
    
    # pull the top covariate & make that a random effect #####
    sigCov <- as.data.frame(car::Anova(MixModOrd))
    sigCov$var <- rownames(sigCov)
    rownames(sigCov) <- NULL
    
    coeftab <- as.data.frame(fixef(MixModOrd))
    coeftab$var <- rownames(coeftab)
    rownames(coeftab) <- NULL
    coeftab <- coeftab[order(coeftab$`fixef(MixModOrd)`, decreasing = TRUE),]
    coeftab2 <- merge(coeftab,sigCov,by='var')
    coeftab2 <- subset(coeftab2, coeftab2$`Pr(>Chisq)`<=0.05)
    topCoeff <- coeftab2[c(1),]
    topCoeff <- as.character(topCoeff$var)
    
    if (!is.na(topCoeff)) {
      
      # var1|species
      topCM <- paste('+ (',topCoeff,' | Name2)',sep='')
      lMm2 <- paste(lm1, topCM, sep='')
      
      MixModOrdEnv <- lmer(as.formula(lMm2), 
                           data = SpScale)
      
      # compare AIC for all order models #####
      
      AICBIC <- anova(MixModOrd,Ord_fulllm)
      AICBIC2 <- anova(MixModOrd,MixModOrdEnv)
      
      AICTableT <- as.data.frame(matrix(ncol=0, nrow = 2))
      AICTableT$Model <- rownames(AICBIC)
      AICTableT$Taxa <- OrdersFive[i]
      AICTableT$AIC <- AICBIC$AIC
      AICTableT$BIC <- AICBIC$BIC
      AICTableT$loglik <- AICBIC$logLik
      
      AICTableTh <- as.data.frame(matrix(ncol=0, nrow = 2))
      AICTableTh$Model <- rownames(AICBIC2)
      AICTableTh$Taxa <- OrdersFive[i]
      AICTableTh$AIC <- AICBIC2$AIC
      AICTableTh$BIC <- AICBIC2$BIC
      AICTableTh$loglik <- AICBIC2$logLik
      
      AICTableOrd <- rbind(AICTableOrd, AICTableT)
      AICTableOrd <- rbind(AICTableOrd, AICTableTh)
      
      # save the r squared for all order models #####
      r.squared <- summary(Ord_fulllm)$r.squared
      r.squared1 <- r.squaredGLMM(MixModOrd)
      r.squared2 <- r.squaredGLMM(MixModOrdEnv)
      
      Rsq1 <- as.data.frame(matrix(ncol=0, nrow = 1))
      Rsq1$Taxa <- OrdersFive[i]
      Rsq1$model <- as.character('Ord_fulllm')
      Rsq1$R2lm <- r.squared
      
      Rsq2 <- as.data.frame(matrix(ncol=0, nrow = 1))
      Rsq2$Taxa <- OrdersFive[i]
      Rsq2$model <- as.character('MixModOrd')
      Rsq2$R2m <- r.squared1[1]
      Rsq2$R2c <- r.squared1[2]
      
      Rsq3 <- as.data.frame(matrix(ncol=0, nrow = 1))
      Rsq3$Taxa <- OrdersFive[i]
      Rsq3$model <- as.character('MixModOrdEnv')
      Rsq3$R2m <- r.squared2[1]
      Rsq3$R2c <- r.squared2[2]
      
      OrdRsq <- plyr::rbind.fill(OrdRsq, Rsq1)
      OrdRsq <- plyr::rbind.fill(OrdRsq, Rsq2)
      OrdRsq <- plyr::rbind.fill(OrdRsq, Rsq3)
      
      taxacount <- as.data.frame(matrix(ncol=2, nrow = 1))
      taxacount$Taxa <-  OrdersFive[i]
      taxacount$Pops <-  nrow(SpScale)
      taxacount$NSp <-  length(unique(SpScale$Name2))
      
      taxaSum <- plyr::rbind.fill(taxaSum, taxacount)
      
      
    }
    else{
      print(SpStart$order.y[1])
      print('no sig covs')
      next
    }
    # next, identify the top 5 most collected families within the order #####
    # if less than 5, OK just take those
    FamTable <- as.data.frame(table(OrderA$FAMILY))
    
    FamTable <- FamTable[order(FamTable$Freq, decreasing = TRUE),]
    
    topFam <- FamTable[c(1:5),]
    
    topFam <- as.character(unique(na.omit(topFam$Var1)))
    
    assign(paste(OrdersFive[i], 'families',sep = '_'), topFam) # rename
    
    # Family level analysis #####
    for (j in 1:length(unique(topFam))){
      print(j)
      print(topFam[j])
      
      # Select i thru 5 families with most populations
      FamilyA <- subset(species4, species4$FAMILY == paste(topFam[j]))
      
      # Prep the dataset to scale potential covariates
      Sp_scale <- FamilyA
      # scale the new data
      Sp_scale[,c(16:29)] <- lapply(Sp_scale[,c(16:29)], 
                                    function(x) c(scale(x)))
      
      Sp_scale[,c(14,15)] <- lapply(Sp_scale[,c(14,15)], 
                                    function(x) c(scale(x, center = TRUE, scale = FALSE)))
      
      # rename these columns with 'scale'
      colnames(Sp_scale) <- paste(colnames(Sp_scale),'_scale',sep = '')
      # combine the scaled data with our original data
      temp <- cbind(FamilyA, Sp_scale[,c(14:29,7,8)])
      
      ## analysis time!
      SpStart <- FamilyA
      rownames(SpStart) <- NULL
      SpScale <- temp
      rownames(SpScale) <- NULL
      
      print('running RF')
      boruta_output <- Boruta(SpStart$mg.seed ~ ., data=(SpStart[,c(14:29,7,8)]), pValue = 0.0001, doTrace=0)  
      
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
      if ((varlength >= 1) && (length(unique(SpScale$Name2))>2)){
        # Create a vector of the scaled variables and response
        Sp_vars2 <- c('ACC_NUM','mg.seed',Sp_envarsScale)
        # Filter the SpScale data frame to only include the important variables
        SpScale2 <- subset(SpScale, select=Sp_vars2)
        # Create a linear regression formula with the scaled variable names
        lm1 <- paste(paste("(mg.seed) ~ coords.x1 + coords.x2 +"), paste(paste(Sp_envarsScale, collapse="+"),  sep='+'))
        # Fit a linear regression model with the formula and scaled data
        Fam_fulllm <- lm(as.formula(lm1),data=SpScale)
        
        lMm1 <- paste(lm1, '+ (1 | Name2)', sep='')
        
        MixModFam <- lmer(as.formula(lMm1), 
                          data = SpScale)
        
        # Plot the histogram and Q-Q plot of residuals
        hist(residuals(MixModFam))
        qqnorm(residuals(MixModFam)); qqline(residuals(MixModFam))
        
        # pull the top covariate & make that a random effect #####
        sigCov <- as.data.frame(car::Anova(MixModFam))
        sigCov$var <- rownames(sigCov)
        rownames(sigCov) <- NULL
        
        coeftab <- as.data.frame(fixef(MixModFam))
        coeftab$var <- rownames(coeftab)
        rownames(coeftab) <- NULL
        coeftab <- coeftab[order(coeftab$`fixef(MixModFam)`, decreasing = TRUE),]
        coeftab2 <- merge(coeftab,sigCov,by='var')
        coeftab2 <- subset(coeftab2, coeftab2$`Pr(>Chisq)`<=0.05)
        topCoeff <- coeftab2[c(1),]
        topCoeff <- as.character(topCoeff$var)
        
        if (!is.na(topCoeff)) {
          # var1|species
          topCM <- paste('+ (',topCoeff,' | Name2)',sep='')
          lMm2 <- paste(lm1, topCM, sep='')
          
          MixModFamEnv <- lmer(as.formula(lMm2), 
                               data = SpScale)
          
          # compare AIC for all family models #####        
          AICBIC <- anova(MixModFam,Fam_fulllm)
          AICBIC2 <- anova(MixModFam,MixModFamEnv)
          
          AICTableT <- as.data.frame(matrix(ncol=0, nrow = 2))
          AICTableT$Model <- rownames(AICBIC)
          AICTableT$Taxa <- topFam[j]
          AICTableT$AIC <- AICBIC$AIC
          AICTableT$BIC <- AICBIC$BIC
          AICTableT$loglik <- AICBIC$logLik
          
          AICTableTh <- as.data.frame(matrix(ncol=0, nrow = 2))
          AICTableTh$Model <- rownames(AICBIC2)
          AICTableTh$Taxa <- topFam[j]
          AICTableTh$AIC <- AICBIC2$AIC
          AICTableTh$BIC <- AICBIC2$BIC
          AICTableTh$loglik <- AICBIC2$logLik
          
          AICTableFam <- rbind(AICTableFam, AICTableT)
          AICTableFam <- rbind(AICTableFam, AICTableTh)
          
          # save the r squared for family models ##### 
          r.squared <- summary(Fam_fulllm)$r.squared
          r.squared1 <- r.squaredGLMM(MixModFam)
          r.squared2 <- r.squaredGLMM(MixModFamEnv)
          
          Rsq1 <- as.data.frame(matrix(ncol=0, nrow = 1))
          Rsq1$Taxa <- topFam[j]
          Rsq1$model <- as.character('Fam_fulllm')
          Rsq1$R2lm <- r.squared
          
          Rsq2 <- as.data.frame(matrix(ncol=0, nrow = 1))
          Rsq2$Taxa <- topFam[j]
          Rsq2$model <- as.character('MixModFam')
          Rsq2$R2m <- r.squared1[1]
          Rsq2$R2c <- r.squared1[2]
          
          Rsq3 <- as.data.frame(matrix(ncol=0, nrow = 1))
          Rsq3$Taxa <- topFam[j]
          Rsq3$model <- as.character('MixModFamEnv')
          Rsq3$R2m <- r.squared2[1]
          Rsq3$R2c <- r.squared2[2]
          
          FamRsq <- plyr::rbind.fill(FamRsq, Rsq1)
          FamRsq <- plyr::rbind.fill(FamRsq, Rsq2)
          FamRsq <- plyr::rbind.fill(FamRsq, Rsq3)
          
          taxacount <- as.data.frame(matrix(ncol=2, nrow = 1))
          taxacount$Taxa <-  topFam[j]
          taxacount$Pops <-  nrow(SpScale)
          taxacount$NSp <-  length(unique(SpScale$Name2))
          
          taxaSum <- plyr::rbind.fill(taxaSum, taxacount)
          
        }
        else{
          print(SpStart$FAMILY[1])
          print('no sig covs')
          next
        }
        # next, identify the top 5 most collected genera within the family #####
        GenusTable <- as.data.frame(table(FamilyA$GENUS))
        
        GenusTable <- GenusTable[order(GenusTable$Freq, decreasing = TRUE),]
        
        topGen <- GenusTable[c(1:5),]
        
        topGen <- as.character(unique(na.omit(topGen$Var1)))
        
        assign(paste(topFam[j], 'genera',sep = '_'), topGen) # rename
        
        # Genus level analysis #####
        for (k in 1:length(topGen)){
          print(k)
          print(topGen[k])
          
          # Select i thru 5 genera with most populations
          GenusA <- subset(species4, species4$GENUS == paste(topGen[k]))
          
          # Prep the dataset to scale potential covariates
          Sp_scale <- GenusA
          # scale the new data
          Sp_scale[,c(16:29)] <- lapply(Sp_scale[,c(16:29)], 
                                        function(x) c(scale(x)))
          
          Sp_scale[,c(14,15)] <- lapply(Sp_scale[,c(14,15)], 
                                        function(x) c(scale(x, center = TRUE, scale = FALSE)))
          
          # rename these columns with 'scale'
          colnames(Sp_scale) <- paste(colnames(Sp_scale),'_scale',sep = '')
          # combine the scaled data with our original data
          temp <- cbind(GenusA, Sp_scale[,c(14:29,7,8)])
          
          ## analysis time!
          SpStart <- GenusA
          rownames(SpStart) <- NULL
          SpScale <- temp
          rownames(SpScale) <- NULL
          
          print('running RF')
          boruta_output <- Boruta(SpStart$mg.seed ~ ., data=(SpStart[,c(14:29,7,8)]), pValue = 0.0001, doTrace=0)  
          
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
          if ((varlength >= 1) && (length(unique(SpScale$Name2))>2)){
            
            # Create a vector of the scaled variables and response
            Sp_vars2 <- c('ACC_NUM','mg.seed',Sp_envarsScale)
            # Filter the SpScale data frame to only include the important variables
            SpScale2 <- subset(SpScale, select=Sp_vars2)
            # Create a linear regression formula with the scaled variable names
            lm1 <- paste(paste("(mg.seed) ~ coords.x1 + coords.x2 +"), 
                         paste(paste(Sp_envarsScale, collapse="+"),  sep='+'))
            # Fit a linear regression model with the formula and scaled data
            Gen_fulllm <- lm(as.formula(lm1),data=SpScale)
            lMm1 <- paste(lm1, '+ (1 | Name2)', sep='')
            
            MixModGen <- lmer(as.formula(lMm1), 
                              data = SpScale)
            
            # Plot the histogram and Q-Q plot of residuals
            hist(residuals(MixModGen))
            qqnorm(residuals(MixModGen)); qqline(residuals(MixModGen))
            
            # pull the top covariate & make that a random effect #####
            sigCov <- as.data.frame(car::Anova(MixModGen))
            sigCov$var <- rownames(sigCov)
            rownames(sigCov) <- NULL
            
            coeftab <- as.data.frame(fixef(MixModGen))
            coeftab$var <- rownames(coeftab)
            rownames(coeftab) <- NULL
            coeftab <- coeftab[order(coeftab$`fixef(MixModGen)`, decreasing = TRUE),]
            coeftab2 <- merge(coeftab,sigCov,by='var')
            coeftab2 <- subset(coeftab2, coeftab2$`Pr(>Chisq)`<=0.05)
            topCoeff <- coeftab2[c(1),]
            topCoeff <- as.character(topCoeff$var)
            
            if (!is.na(topCoeff)) {
              # var1|species
              topCM <- paste('+ (',topCoeff,' | Name2)',sep='')
              lMm2 <- paste(lm1, topCM, sep='')
              
              MixModGenEnv <- lmer(as.formula(lMm2), 
                                   data = SpScale)
              
              # compare AIC for all genus models #####            
              AICBIC <- anova(MixModGen,Gen_fulllm)
              AICBIC2 <- anova(MixModGen,MixModGenEnv)
              
              AICTableT <- as.data.frame(matrix(ncol=0, nrow = 2))
              AICTableT$Model <- rownames(AICBIC)
              AICTableT$Taxa <- topGen[k]
              AICTableT$AIC <- AICBIC$AIC
              AICTableT$BIC <- AICBIC$BIC
              AICTableT$loglik <- AICBIC$logLik
              
              AICTableTh <- as.data.frame(matrix(ncol=0, nrow = 2))
              AICTableTh$Model <- rownames(AICBIC2)
              AICTableTh$Taxa <- topGen[k]
              AICTableTh$AIC <- AICBIC2$AIC
              AICTableTh$BIC <- AICBIC2$BIC
              AICTableTh$loglik <- AICBIC2$logLik
              
              AICTableGen <- rbind(AICTableGen, AICTableT)
              AICTableGen <- rbind(AICTableGen, AICTableTh)
              
              # save the r squared for genus models #####
              r.squared <- summary(Gen_fulllm)$r.squared
              r.squared1 <- r.squaredGLMM(MixModGen)
              r.squared2 <- r.squaredGLMM(MixModGenEnv)
              
              Rsq1 <- as.data.frame(matrix(ncol=0, nrow = 1))
              Rsq1$Taxa <- topGen[k]
              Rsq1$model <- as.character('Gen_fulllm')
              Rsq1$R2lm <- r.squared
              
              Rsq2 <- as.data.frame(matrix(ncol=0, nrow = 1))
              Rsq2$Taxa <- topGen[k]
              Rsq2$model <- as.character('MixModGen')
              Rsq2$R2m <- r.squared1[1]
              Rsq2$R2c <- r.squared1[2]
              
              Rsq3 <- as.data.frame(matrix(ncol=0, nrow = 1))
              Rsq3$Taxa <- topGen[k]
              Rsq3$model <- as.character('MixModGenEnv')
              Rsq3$R2m <- r.squared2[1]
              Rsq3$R2c <- r.squared2[2]
              
              GenRsq <- plyr::rbind.fill(GenRsq, Rsq1)
              GenRsq <- plyr::rbind.fill(GenRsq, Rsq2)
              GenRsq <- plyr::rbind.fill(GenRsq, Rsq3)
              
              taxacount <- as.data.frame(matrix(ncol=2, nrow = 1))
              taxacount$Taxa <-  topGen[k]
              taxacount$Pops <-  nrow(SpScale)
              taxacount$NSp <-  length(unique(SpScale$Name2))
              
              taxaSum <- plyr::rbind.fill(taxaSum, taxacount)
              
            }
            else{
              print(SpStart$GENUS[1])
              print('no sig covs')
              next
            }
          }
          #if there are no important variables at all print family name and "no variables"
          else {
            print(SpStart$GENUS[1])
            print("no vars")
            next
          }
        } 
      }
      else {
        print(SpStart$FAMILY[1])
        print("no vars")
        next
      }
      
    }
  }
  #if there are no important variables at all for order 
  # print order name and "no variables"
  else {
    print(SpStart$order.y[1])
    print("no vars")
    next
  } 
}

#### top fam gen #####

for (i in 1:length(OrdersFive)) {
  print(i)
  print(OrdersFive[i])
  # Order level analysis
  
  # Select i thru 5 orders with >500 populations
  OrderA <- subset(species4, species4$order.y == paste(OrdersFive[i]))
  
  FamTable <- as.data.frame(table(OrderA$FAMILY))
  
  FamTable <- FamTable[order(FamTable$Freq, decreasing = TRUE),]
  
  topFam <- FamTable[c(1:5),]
  
  topFam <- as.character(unique(na.omit(topFam$Var1)))
  
  assign(paste(OrdersFive[i], 'families',sep = '_'), topFam) # rename
  
  for (j in 1:length(unique(topFam))){
    print(j)
    print(topFam[j])
    
    # Select i thru 5 families with most populations
    FamilyA <- subset(species4, species4$FAMILY == paste(topFam[j]))
    
    # next, identify the top 5 most collected genera within the family #####
    GenusTable <- as.data.frame(table(FamilyA$GENUS))
    
    GenusTable <- GenusTable[order(GenusTable$Freq, decreasing = TRUE),]
    
    topGen <- GenusTable[c(1:5),]
    
    topGen <- as.character(unique(na.omit(topGen$Var1)))
    
    assign(paste(topFam[j], 'genera',sep = '_'), topGen) # rename
    
  }
  
}

########################### Individual species year models ########################### 

spnames = c(unique(species3$Name2))
spnames <- sort(spnames)

# Check variable names and column numbers
data.frame(colnames(species3)) # print column names

# Look at the following variables:
# - Mg.seed as response 
# - Year

########################### Jitter duplicate points ########################### 

# k nearest neighbor doesn't work on duplicate points so jitter slightly
coords <- data.frame(lon = species3$coords.x1, lat = species3$coords.x2)

coords_dup <- coords[duplicated(coords[,c('lon', 'lat')]),] 
coords_distinct <- distinct(coords_dup)

fin_list <- list()
for (i in 1:nrow(coords_distinct)) {
  temp_data <- species3 %>% filter(coords.x1 == coords_distinct[i,1] & coords.x2 == coords_distinct[i,2])
  
  temp_data$coords.x1 <- temp_data$coords.x1 + seq(0, 0.0001, 0.00001)[1:nrow(temp_data)]
  temp_data$coords.x2 <- temp_data$coords.x2 + seq(0, 0.0001, 0.00001)[1:nrow(temp_data)]
  
  fin_list[[i]] <- temp_data
  print(i)
}

fin_bind <- data.table::rbindlist(fin_list)

numbers <- c(fin_bind$ACC_NUM)
not_dupe <- species3 %>% filter(!ACC_NUM %in% numbers)

modified_data <- rbind(fin_bind, not_dupe)
species3 <- modified_data

###########################  create individual dataframes ########################### 
for (i in 1:length(spnames)) {
  # first, bring in the dataframe 
  Sp <- spnames[[i]]
  SpeciesDF <- subset(species3, species3$Name2 == paste(Sp))
  assign(paste(SpeciesDF$GENUS, SpeciesDF$SPECIES, sep = '_'), SpeciesDF) # rename
  print(i)
}

data.frame(colnames(species3)) 

##### individual year by species models ####
library(tidyr)
# create an empty matrix to store the results
mymatrix <- matrix(ncol=9, nrow = length(spnames))
spglm2 <- as.data.frame(mymatrix)

# add column names to the result data frame
colnames(spglm2) <- c('Species','Pops','variables','Estimate','Std.Error','zvalue','pvalues','upper','lower')

# Begin looping through each species to determine important variables and run a spatial GLM on them
for (j in 1:length(spnames)){
  
  print(j)
  SpStart <- get(paste(spnames[j], sep='')) # get the unscaled dataframe
  rownames(SpStart) <- NULL
  mymatrix <- matrix(ncol=8, nrow = 1) # make an empty matrix
  df <- as.data.frame(mymatrix)
  colnames(df) <- c('Species','variables','Estimate','Std.Error','zvalue','pvalues','upper','lower')
  assign(paste(SpStart$Name2[1],'df',sep = '_'),df) # rename for sp
  print(SpStart$Name2[1])
  if (length(unique(SpStart$Year)) >= 3) {
  # Create a linear regression formula with the scaled variable names
  print('running model')
    lm1 <- paste("(mg.seed_z) ~ Year_scale")

  # Get the latitude and longitude 
  col.w <- nb2listw(knn2nb(knearneigh(coordinates(SpStart[,c(11,12)]))))
  tryCatch({
    serrRslt <- spatialreg::errorsarlm(as.formula(lm1),
                                       data = SpStart,
                                       listw = col.w,
                                       zero.policy = TRUE, 
                                       na.action = na.omit);
    assign(paste(spnames[j],'SAR',sep='_'),serrRslt)
    
    # Plot the histogram and Q-Q plot of residuals
    hist(residuals(serrRslt))
    qqnorm(residuals(serrRslt)); qqline(residuals(serrRslt))
    
    # Create a data frame with model coefficients, p-values, and other metrics
    
    testtable <- summary(serrRslt)
    testtable2 <- testtable$Coef
    tablee <- as.data.frame(testtable2)
    tablee$variables <- rownames(tablee)
    rownames(tablee) <- NULL
    
    confsp <- as.data.frame(confint(serrRslt))
    tablee$lower <- confsp$`2.5 %`[c(2,3)]
    tablee$upper <- confsp$`97.5 %`[c(2,3)]
    
    # Add test stats to the data frame
    tablee$pvalues <- tablee$`Pr(>|z|)`
    tablee$Std.Error <- tablee$`Std. Error`
    tablee$zvalue <- tablee$`z value`
    
    tablee$Species <- SpStart$Name2[1]
    tablee$Pops <- nrow(SpStart)
    
    
    # Rename the data for use with spatial packages
    assign(paste(SpStart$Name2[1],'df',sep = '_'),tablee)
    
    common_cols <- intersect(colnames(tablee), colnames(spglm2))
    spglm2 <- rbind(
      subset(tablee, select = common_cols), 
      subset(spglm2, select = common_cols)
    )
  }, warning = function(w) {
    # Print "error" if a warning is encountered
    message("error")
    # Move to the next iteration
    next
                            }, 
    error = function(e) {
    # Print "error" if an error is encountered
    message("error")
    # Move to the next iteration
                          }
    )
  next
  }
  else {
    print(paste(SpStart$Name2[1],'Not enough years',sep=''))
       }
}


# check how many species have changed seed mass over time (years) 

years <- subset(spglm2, spglm2$variables == 'Year_scale')

sigsummary_table <- years %>%
  filter(pvalues <= 0.05) %>%                      # Only significant results
  mutate(direction = ifelse(Estimate > 0, "Positive", "Negative")) %>%  # Direction of effect
  count(direction)                             # Tally the counts

print(sigsummary_table)

summary_table <- years %>%
  mutate(direction = ifelse(Estimate > 0, "Positive", "Negative")) %>%  # Direction of effect
  count(direction)                             # Tally the counts

print(summary_table)

yearSp <- unique(years$Species)

species3_subset <- species3 %>%
  filter(Name2 %in% yearSp)

moments::kurtosis(years$Estimate)
anscombe.test(years$Estimate, alternative = c("two.sided"))
anscombe.test(years$Estimate, alternative = c("greater"))
anscombe.test(years$Estimate, alternative = c("less"))

plot<-ggplot(years, aes(x = Estimate)) +
  geom_histogram(fill = "lightgreen", color = "black", bins = 50) +
  geom_vline(xintercept = mean(years$Estimate), color = "red", linetype = "dashed", size = 1) +  # Add a vertical line for the mean
  #annotate("text", x = 1.5, y = 200, label = paste("t =", round(t_test_result$statistic, 2), ", p =", round(t_test_result$p.value, 4)), hjust = 1, vjust = 1, size =8) +
  labs(x = "Year Slopes", y = "Frequency") +
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 50),
        axis.title = element_text(size = 40, colour = "black", face = "bold"),
        panel.border = element_rect(linewidth = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(linewidth = 1.5, colour = "#333333", fill = "#CCCCCC"),
        axis.ticks.y = element_blank())
plot

shapiro.test(years$Estimate)
qqnorm(years$Estimate)
qqline(years$Estimate, col = "red")
text(x = -1.5, y = .5, 
     labels = paste("W =", round(shapiro_result$statistic, 5), 
                    "\nP-value =", formatC(shapiro_result$p.value, format = "e", digits = 2)),
     pos = 4, cex = 0.8, col = "blue")
text(x = -2.5, y = 1, 
     labels = paste("W =", round(ks_result$statistic, 5), 
                    "\nP-value =", formatC(round(ks_result$p.value, 5), format = "e", digits = 2)),
     pos = 4, cex = 0.8, col = "blue")


estHist <- ggplot(years, aes(x = Estimate)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightgray", color = "black") +
  stat_function(fun = dnorm,
                args = list(mean = mean(years$Estimate),
                            sd = sd(years$Estimate)),
                color = "blue", size = 1) +
  labs(x = "Year slope estimate",
       y = "Density")

prophist <- ggplot(merged_data2, aes(x = prop)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightgray", color = "black") +
  stat_function(fun = dnorm,
                args = list(mean = mean(merged_data2$prop),
                            sd = sd(merged_data2$prop)),
                color = "blue", size = 1) +
  labs(x = "Year slope estimate",
       y = "Density")

t_test_result <- t.test(years$Estimate, mu = 0)
t_test_result$statistic
