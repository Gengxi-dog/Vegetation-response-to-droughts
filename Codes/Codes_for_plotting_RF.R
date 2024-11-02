load('kNDVI_vulnerability_RF_data.RData')
library(ranger)
library(edarf)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(zoo)

response_variable <- "vulner"
predictor_variables <- c('CO2', 'pre_summer', 'tmp_summer', 
                         'rad_summer', 'VPD_summer', 'SOS', 
                         'kNDVI_spring', 'biodiv','AI')

RF_data_evergreen <- dplyr::filter(RF_data, landuse == 1 | landuse == 2)
RF_data_evergreen <- na.omit(RF_data_evergreen)
data_subset_evergreen <- RF_data_evergreen[, c(response_variable, predictor_variables)]

training_indices <- sample(nrow(data_subset_evergreen), nrow(data_subset_evergreen) * 0.7)
training_data_evergreen <- data_subset_evergreen[training_indices, ]
testing_data_evergreen <- data_subset_evergreen[-training_indices, ]

RF_evergreen_training_result <- ranger(data = training_data_evergreen, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, write.forest = T)
RF_evergreen_predict_result <- predict(RF_evergreen_training_result, testing_data_evergreen[,-1])

# evaluate accuracy, in this example, R is used. We can also use other metrics to estimate the regression accuracy.
R_evergreen <- cor(testing_data_evergreen$vulner, RF_evergreen_predict_result$predictions)

RF_evergreen_result <- ranger(data = data_subset_evergreen, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                                       write.forest = T)
pd_evergreen <- partial_dependence(fit = RF_evergreen_result, vars = c('CO2', 'pre_summer', 'tmp_summer', 
     'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring', 'biodiv','AI'), data = data_subset_evergreen[,-1])

#------------------------------------------------------------
RF_data_deciduous <- dplyr::filter(RF_data, landuse == 3 | landuse == 4)
RF_data_deciduous <- na.omit(RF_data_deciduous)
data_subset_deciduous <- RF_data_deciduous[, c(response_variable, predictor_variables)]

training_indices <- sample(nrow(data_subset_deciduous), nrow(data_subset_deciduous) * 0.7)
training_data_deciduous <- data_subset_deciduous[training_indices, ]
testing_data_deciduous <- data_subset_deciduous[-training_indices, ]

RF_deciduous_training_result <- ranger(data = training_data_deciduous, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                              write.forest = T)
RF_deciduous_predict_result <- predict(RF_deciduous_training_result, testing_data_deciduous[,-1])

# evaluate accuracy, in this example, R is used. We can also use other metrics to estimate the regression accuracy.
R_deciduous <- cor(testing_data_deciduous$vulner, RF_deciduous_predict_result$predictions)

RF_deciduous_result <- ranger(data = data_subset_deciduous, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                                       write.forest = T)
pd_deciduous <- partial_dependence(fit = RF_deciduous_result, vars = c('CO2', 'pre_summer', 'tmp_summer', 
       'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring', 'biodiv','AI'), data = data_subset_deciduous[,-1])

#------------------------------------------------------
RF_data_savannas <- dplyr::filter(RF_data, landuse == 8 | landuse == 9)
RF_data_savannas <- na.omit(RF_data_savannas)
data_subset_savannas <- RF_data_savannas[, c(response_variable, predictor_variables)]

training_indices <- sample(nrow(data_subset_savannas), nrow(data_subset_savannas) * 0.7)
training_data_savannas <- data_subset_savannas[training_indices, ]
testing_data_savannas <- data_subset_savannas[-training_indices, ]

RF_savannas_training_result <- ranger(data = training_data_savannas, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                              write.forest = T)
RF_savannas_predict_result <- predict(RF_savannas_training_result, testing_data_savannas[,-1])

# evaluate accuracy, in this example, R is used. We can also use other metrics to estimate the regression accuracy.
R_savannas <- cor(testing_data_savannas$vulner, RF_savannas_predict_result$predictions)

RF_savannas_result <- ranger(data = data_subset_savannas, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                                      write.forest = T)
pd_savannas <- partial_dependence(fit = RF_savannas_result, vars = c('CO2', 'pre_summer', 'tmp_summer', 
                'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring', 'biodiv','AI'), data = data_subset_savannas[,-1])

#------------------------------------------------------
RF_data_shrubs <- dplyr::filter(RF_data, landuse == 6 | landuse == 7)
RF_data_shrubs <- na.omit(RF_data_shrubs)
data_subset_shrubs <- RF_data_shrubs[, c(response_variable, predictor_variables)]

training_indices <- sample(nrow(data_subset_shrubs), nrow(data_subset_shrubs) * 0.7)
training_data_shrubs <- data_subset_shrubs[training_indices, ]
testing_data_shrubs <- data_subset_shrubs[-training_indices, ]

RF_shrubs_training_result <- ranger(data = training_data_shrubs, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                              write.forest = T)
RF_shrubs_predict_result <- predict(RF_shrubs_training_result, testing_data_shrubs[,-1])

# evaluate accuracy, in this example, R is used. We can also use other metrics to estimate the regression accuracy.
R_shrubs <- cor(testing_data_shrubs$vulner, RF_shrubs_predict_result$predictions)

RF_shrubs_result <- ranger(data = data_subset_shrubs, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                                    write.forest = T)
pd_shrubs <- partial_dependence(fit = RF_shrubs_result, vars = c('CO2', 'pre_summer', 'tmp_summer', 
                'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring', 'biodiv','AI'), data = data_subset_shrubs[,-1])

#------------------------------------------------------
RF_data_grass <- dplyr::filter(RF_data, landuse == 10)
RF_data_grass <- na.omit(RF_data_grass)
data_subset_grass <- RF_data_grass[, c(response_variable, predictor_variables)]

training_indices <- sample(nrow(data_subset_grass), nrow(data_subset_grass) * 0.7)
training_data_grass <- data_subset_grass[training_indices, ]
testing_data_grass <- data_subset_grass[-training_indices, ]

RF_grass_training_result <- ranger(data = training_data_grass, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                              write.forest = T) 
RF_grass_predict_result <- predict(RF_grass_training_result, testing_data_grass[,-1])

# evaluate accuracy, in this example, R is used. We can also use other metrics to estimate the regression accuracy.
R_grass <- cor(testing_data_grass$vulner, RF_grass_predict_result$predictions)

RF_grass_result <- ranger(data = data_subset_grass, formula = vulner~., importance = 'permutation', num.trees = 500, oob.error = T, 
                          write.forest = T)
pd_grass <- partial_dependence(fit = RF_grass_result, vars = c('CO2', 'pre_summer', 'tmp_summer', 
             'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring', 'biodiv','AI'), data = data_subset_grass[,-1])

# plotting importance of RF models for different vegetation types
#------------------------------------------------------------
IncMSE_vulner_EVGR <- 100*(RF_evergreen_result$variable.importance/sum(RF_evergreen_result$variable.importance))
IncMSE_vulner_DEDU <- 100*(RF_deciduous_result$variable.importance/sum(RF_deciduous_result$variable.importance))
IncMSE_vulner_SAV <- 100*(RF_savannas_result$variable.importance/sum(RF_savannas_result$variable.importance))
IncMSE_vulner_SHR <- 100*(RF_shrubs_result$variable.importance/sum(RF_shrubs_result$variable.importance))
IncMSE_vulner_GRA <- 100*(RF_grass_result$variable.importance/sum(RF_grass_result$variable.importance))

importance_data <- as.data.frame(cbind(c(IncMSE_vulner_EVGR, IncMSE_vulner_DEDU, IncMSE_vulner_SAV, 
                           IncMSE_vulner_SHR, IncMSE_vulner_GRA), c(rep('EVGR', length(IncMSE_vulner_EVGR)), 
                           rep('DEDU', length(IncMSE_vulner_DEDU)), rep('SAV', length(IncMSE_vulner_SAV)), 
                          rep('SHR', length(IncMSE_vulner_SHR)), rep('GRA', length(IncMSE_vulner_GRA))), 
                          rep(c('CO2', 'pre_summer', 'tmp_summer','rad_summer', 'VPD_summer', 
                                'SOS','kNDVI_spring', 'biodiv','AI'), 5)))
names(importance_data) <- c('importance','vegetation_type','factors')
importance_data$importance <- as.numeric(importance_data$importance)

library(ggplot2)
library(ggsci)
library(forcats)

importance_EVGR_data <- dplyr::filter(importance_data, vegetation_type == 'EVGR')
importance_DEDU_data <- dplyr::filter(importance_data, vegetation_type == 'DEDU')
importance_SAV_data <- dplyr::filter(importance_data, vegetation_type == 'SAV')
importance_SHR_data <- dplyr::filter(importance_data, vegetation_type == 'SHR')
importance_GRA_data <- dplyr::filter(importance_data, vegetation_type == 'GRA')

plot_importance <- ggplot(importance_data, aes(x = factors, y=importance, color = vegetation_type)) + 
    geom_point(size = 5) + 
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
    theme_bw() + 
    theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
    # theme(legend.position = "none") + 
    labs(x="Factors", y="Importance") +
    scale_x_discrete(limit = c('AI','biodiv','tmp_summer','pre_summer','rad_summer','VPD_summer','kNDVI_spring','SOS','CO2')) +
    scale_y_continuous(limits = c(3,20)) +
    theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(axis.ticks.length = unit(-0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
    theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

names(pd_evergreen) <- c('CO2', 'pre_summer', 'tmp_summer','rad_summer', 'VPD_summer', 
                         'SOS','kNDVI_spring', 'biodiv','AI','vulner')
names(pd_deciduous) <- c('CO2', 'pre_summer', 'tmp_summer','rad_summer', 'VPD_summer', 
                         'SOS','kNDVI_spring', 'biodiv','AI','vulner')
names(pd_savannas) <- c('CO2', 'pre_summer', 'tmp_summer','rad_summer', 'VPD_summer', 
                        'SOS','kNDVI_spring', 'biodiv','AI','vulner')
names(pd_shrubs) <- c('CO2', 'pre_summer', 'tmp_summer','rad_summer', 'VPD_summer', 
                      'SOS','kNDVI_spring', 'biodiv','AI','vulner')
names(pd_grass) <- c('CO2', 'pre_summer', 'tmp_summer','rad_summer', 'VPD_summer', 
                     'SOS','kNDVI_spring', 'biodiv','AI','vulner')

pd_data <- as.data.frame(cbind(as.data.frame(rbind(pd_evergreen, pd_deciduous, pd_savannas, pd_shrubs, pd_grass)), 
                 rep(c('EVGR','DEDU','SAV','SHR','GRA'), each = 225)))
names(pd_data) <- c('CO2', 'pre_summer', 'tmp_summer','rad_summer', 'VPD_summer', 
                    'SOS','kNDVI_spring', 'biodiv','AI','vulner','Vtype')

#----------------------------------------------
pd_plot_AI <- ggplot(pd_data, aes(x=AI, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='AI', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(2000,50000)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

#----------------------------------------------
pd_plot_biodiv <- ggplot(pd_data, aes(x=biodiv, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='biodiv', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(500,4000)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

#----------------------------------------------
pd_plot_pre <- ggplot(pd_data, aes(x=pre_summer, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='pre_summer', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(-2.5,2.5)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

#----------------------------------------------
pd_plot_tmp <- ggplot(pd_data, aes(x=tmp_summer, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='tmp_summer', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(-2.5,2.5)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

#----------------------------------------------
pd_plot_rad <- ggplot(pd_data, aes(x=rad_summer, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='rad_summer', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(-2.5,2.5)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

#----------------------------------------------
pd_plot_VPD <- ggplot(pd_data, aes(x=VPD_summer, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='VPD_summer', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(-2.5,2.5)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

#----------------------------------------------
pd_plot_kNDVI <- ggplot(pd_data, aes(x=kNDVI_spring, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='kNDVI', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(-4,4)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

#----------------------------------------------
pd_plot_SOS <- ggplot(pd_data, aes(x=SOS, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='SOS', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(-4,4)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

#----------------------------------------------
pd_plot_CO2 <- ggplot(pd_data, aes(x=CO2, y=vulner, color=Vtype)) + 
  geom_smooth(method = 'auto', se=F) +
  # geom_ribbon(alpha = 0.5) +
  scale_color_cosmic("hallmarks_light", breaks = c('EVGR','DEDU','SAV','SHR','GRA')) +
  theme_bw() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 20)) +
  labs(x='CO2', y="Susceptibility (%)") +
  theme(axis.title = element_text(color = 'black', size = 20), axis.text = element_text(color = 'black', size = 20)) + 
  scale_x_continuous(limits = c(-1.5,2)) + scale_y_continuous(limits = c(0.2,0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

plot_pd <- ggarrange(pd_plot_AI, pd_plot_pre, pd_plot_tmp, pd_plot_biodiv, pd_plot_kNDVI, 
                     pd_plot_SOS, nrow = 2, ncol = 3, common.legend = T)

ggarrange(plot_pd, plot_importance, nrow = 2, heights = c(1.5,1), common.legend = T)
