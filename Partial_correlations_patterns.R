library(ggplot2)
library(raster)
library(tidyverse)
library(ggsci)

pcor_5.5 <- brick('Data/vulner_pcor_to_7variables_5.nc')
pcor_pvalue_5.5 <- brick('Data/vulner_pcor_pvalue_to_7variables_5.nc')

pcor_5.5[pcor_pvalue_5.5 > 0.05] <- NA

pcor_5.5_arr <- as.array(pcor_5.5)
pcor_5.5_arr_abs <- abs(pcor_5.5_arr)

pcor_control <- matrix(NA, nrow = 241, ncol = 1440)
for (i in 1:241){
  for (j in 1:1440){
    if(sum(!is.na(pcor_5.5_arr_abs[i,j,])) > 0){
      pcor_control[i,j] <- which.max(pcor_5.5_arr_abs[i,j,])
    }
  }
}
pcor_control_ras <- raster(pcor_control, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

pcor_control_data <- raster::as.data.frame(pcor_control_ras, xy = T)
names(pcor_control_data) <- c('x','y','layer')
pcor_control_data$layer <- as.character(pcor_control_data$layer)

wm_ggplot <- map_data("world")
wm_ggplot <- dplyr::filter(wm_ggplot, lat > 30)
x_lines <- seq(-120,180, by = 60)

#-------------------------------------------
plot_control_spatial <- ggplot() +
    geom_polygon(data = wm_ggplot, aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_tile(data = pcor_control_data, mapping = aes(x = x,y = y,fill = layer), alpha = 0.8) +
    scale_fill_brewer(palette = "Set1") +
    # scale_fill_gradient2(low = lowcolor, high = highcolor, mid = midcolor, midpoint = midpoint1, na.value = NA) +
    # geom_point(data= data_set2, mapping = aes(x=x,y=y,color=layer), size=0.005, show.legend = F) + scale_color_continuous(low = "black", high = "black", na.value=NA) +
    
    # Convert to polar coordinates
    coord_map("ortho", orientation = c(90, 0, 0)) +
    scale_y_continuous(breaks = seq(30, 90, by = 10), labels = NULL) +
    
    # Removes Axes and labels
    scale_x_continuous(breaks = NULL) +
    xlab("") + 
    ylab("") +
    
    # Adds labels
    # geom_text(aes(x = 180, y = seq(30, 90, by = 10), hjust = -0.2, label = paste0(seq(30, 90, by = 10), "??N"))) +
    # geom_text(aes(x = x_lines, y = 30, label = c("120??W", "60??W", "0??", "60??E", "120??E", "180??W"))) +
    
    # Adds axes
    geom_hline(aes(yintercept = 30), size = 0.5)  +
    geom_segment(aes(y = 30, yend = 90, x = x_lines, xend = x_lines), linetype = "dashed", color = 'grey50') +
    
    # Change theme to remove axes and ticks
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.5, linetype = 'dashed', colour = "grey50"), axis.ticks=element_blank())+
    theme(legend.direction = 'vertical', legend.title = element_blank(), legend.spacing.x = unit(0,'cm'), legend.position = 'right',legend.key = element_rect(color="black"), legend.key.height = unit(1, 'cm'), legend.key.width = unit(0.5,'cm'), legend.text = element_text(size = 20))


percent_control <- cbind(as.data.frame(c(100*sum(pcor_control_data$layer == 1, na.rm = T)/sum(!is.na(pcor_control_data$layer)),100*sum(pcor_control_data$layer == 2, na.rm = T)/sum(!is.na(pcor_control_data$layer)),100*sum(pcor_control_data$layer == 3, na.rm = T)/sum(!is.na(pcor_control_data$layer)),100*sum(pcor_control_data$layer == 4, na.rm = T)/sum(!is.na(pcor_control_data$layer)),100*sum(pcor_control_data$layer == 5, na.rm = T)/sum(!is.na(pcor_control_data$layer)),100*sum(pcor_control_data$layer == 6, na.rm = T)/sum(!is.na(pcor_control_data$layer)),100*sum(pcor_control_data$layer == 7, na.rm = T)/sum(!is.na(pcor_control_data$layer)))), c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring'))
names(percent_control) <- c('percent', 'factor')

#-------------------------------------------
plot_control_barplot <- ggplot(percent_control) + geom_col(aes(x = factor, y = percent, fill = factor), alpha = 0.8, width = 0.7, color = 'black') +
  scale_fill_brewer(palette = "Set1")  + 
  labs(x = 'Factors', y = 'Percentage (%)') +
  theme_bw() + 
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(colour = 'black', size = 20), axis.text = element_text(colour = 'black', size = 20)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) + theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1)) + scale_y_continuous(limits = c(0,20), expand = c(0,0)) + 
  theme(plot.margin = margin(1,1,1,1,unit = 'cm'))

library(ggpubr)
ggarrange(plot_control_spatial, plot_control_barplot, ncol = 2)


#------------------------------------------------
pcor_control_data$layer <- as.numeric(pcor_control_data$layer)
pcor_data_pre <- raster::as.data.frame(pcor_5.5$X2, xy = T)
pcor_data_pre$X2[which.max(pcor_data_pre$X2)] <- 1
pcor_data_pre$X2[which.min(pcor_data_pre$X2)] <- -1
pcor_data_pre$X2[pcor_control_data$layer != 2] <- NA

pcor_data_tmp <- raster::as.data.frame(pcor_5.5$X3, xy = T)
pcor_data_tmp$X3[which.max(pcor_data_tmp$X3)] <- 1
pcor_data_tmp$X3[which.min(pcor_data_tmp$X3)] <- -1
pcor_data_tmp$X3[pcor_control_data$layer != 3] <- NA

pcor_data_rad <- raster::as.data.frame(pcor_5.5$X4, xy = T)
pcor_data_rad$X4[which.max(pcor_data_rad$X4)] <- 1
pcor_data_rad$X4[which.min(pcor_data_rad$X4)] <- -1
pcor_data_rad$X4[pcor_control_data$layer != 4] <- NA

pcor_data_VPD <- raster::as.data.frame(pcor_5.5$X5, xy = T)
pcor_data_VPD$X5[which.max(pcor_data_VPD$X5)] <- 1
pcor_data_VPD$X5[which.min(pcor_data_VPD$X5)] <- -1
pcor_data_VPD$X5[pcor_control_data$layer != 5] <- NA

pcor_data_SOS <- raster::as.data.frame(pcor_5.5$X6, xy = T)
pcor_data_SOS$X6[pcor_control_data$layer != 6] <- NA
pcor_data_SOS$X6[which.max(pcor_data_SOS$X6)] <- 1
pcor_data_SOS$X6[which.min(pcor_data_SOS$X6)] <- -1

pcor_data_NDVI <- raster::as.data.frame(pcor_5.5$X7, xy = T)
pcor_data_NDVI$X7[pcor_control_data$layer != 7] <- NA
pcor_data_NDVI$X7[which.max(pcor_data_NDVI$X7)] <- 1
pcor_data_NDVI$X7[which.min(pcor_data_NDVI$X7)] <- -1

pcor_spatial <- function(data1) {ggplot() +
  geom_polygon(data = wm_ggplot, aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  geom_tile(data = data1, mapping = aes(x = x,y = y,fill = layer), alpha = 1) +
  scale_fill_gradient2(low = '#004c99', high = '#990000', mid = 'white', midpoint = 0, na.value = NA) +
  # scale_fill_gradient2(low = lowcolor, high = highcolor, mid = midcolor, midpoint = midpoint1, na.value = NA) +
  # geom_point(data= data_set2, mapping = aes(x=x,y=y,color=layer), size=0.005, show.legend = F) + scale_color_continuous(low = "black", high = "black", na.value=NA) +
  
  # Convert to polar coordinates
  coord_map("ortho", orientation = c(90, 0, 0)) +
  scale_y_continuous(breaks = seq(30, 90, by = 10), labels = NULL) +
  
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  
  # Adds labels
  # geom_text(aes(x = 180, y = seq(30, 90, by = 10), hjust = -0.2, label = paste0(seq(30, 90, by = 10), "??N"))) +
  # geom_text(aes(x = x_lines, y = 30, label = c("120??W", "60??W", "0??", "60??E", "120??E", "180??W"))) +
  
  # Adds axes
  geom_hline(aes(yintercept = 30), size = 0.5)  +
  geom_segment(aes(y = 30, yend = 90, x = x_lines, xend = x_lines), linetype = "dashed", color = 'grey50') +
  
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'dashed', colour = "grey50"), axis.ticks=element_blank())+
  theme(legend.direction = 'vertical', legend.title = element_blank(), legend.spacing.x = unit(0,'cm'), legend.position = 'right',legend.key = element_rect(color="black"), legend.key.height = unit(1, 'cm'), legend.key.width = unit(0.5,'cm'), legend.text = element_text(size = 20))
}

names(pcor_data_pre) <- c('x','y','layer')
names(pcor_data_tmp) <- c('x','y','layer')
names(pcor_data_rad) <- c('x','y','layer')
names(pcor_data_VPD) <- c('x','y','layer')
names(pcor_data_SOS) <- c('x','y','layer')
names(pcor_data_NDVI) <- c('x','y','layer')

pcor_pre_spatial <- pcor_spatial(data1 = pcor_data_pre)
pcor_tmp_spatial <- pcor_spatial(data1 = pcor_data_tmp)
pcor_rad_spatial <- pcor_spatial(data1 = pcor_data_rad)
pcor_VPD_spatial <- pcor_spatial(data1 = pcor_data_VPD)
pcor_SOS_spatial <- pcor_spatial(data1 = pcor_data_SOS)
pcor_NDVI_spatial <- pcor_spatial(data1 = pcor_data_NDVI)

p2 <- ggarrange(plot_control_spatial, plot_control_barplot, pcor_pre_spatial, pcor_tmp_spatial, pcor_rad_spatial,  pcor_VPD_spatial, pcor_SOS_spatial, pcor_NDVI_spatial, nrow = 4, ncol = 2, legend = 'none') + theme(plot.margin = margin(2,2,2,2,'cm'))

ggsave(filename = 'E:/p5.pdf', device = 'pdf', width = 10, height = 20, units = c('in'), dpi = 300)

p3 <- ggarrange(plot_control_spatial, plot_control_barplot, pcor_SOS_spatial, pcor_NDVI_spatial, nrow = 2, ncol = 2, legend = 'none') + theme(plot.margin = margin(2,2,2,2,'cm'))

ggsave(filename = 'E:/p5.pdf', device = 'pdf', width = 10, height = 20, units = c('in'), dpi = 300)

global_ai <- raster('E:/scientif_data/global-et0-ai/global-ai_et0/ai_et0/ai_et0.tif')
global_ai_NH <- resample(global_ai, y = pcor_5.5$X1)
global_ai_NH_data <- raster::as.data.frame(global_ai_NH$ai_et0, xy = T)

biodiversity_globle <- raster('E:/papers_material/vegetation drought recovery/biotic factors/biodiversity/biodiversity_globe_025res.tif')
biodiversity_NH <- resample(biodiversity_globle, y = pcor_5.5$X1)
remove(biodiversity_globle)
biodiversity_NH_data <- raster::as.data.frame(biodiversity_NH, xy = T)

pcor_SOS_data_lat <- matrix(pcor_data_SOS$X6, nrow = 1440)
pcor_SOS_data_lat_mean <- apply(pcor_SOS_data_lat, MARGIN = 2, FUN = mean, na.rm = T)
pcor_SOS_lat_sd <- apply(pcor_SOS_data_lat, MARGIN = 2, FUN = sd, na.rm = T)
pcor_SOS_lat_data <- as.data.frame(cbind(seq(90, 30, -0.25), pcor_SOS_data_lat_mean + pcor_SOS_lat_sd/2, pcor_SOS_data_lat_mean, pcor_SOS_data_lat_mean - pcor_SOS_lat_sd/2)) 
names(pcor_SOS_lat_data) <- c("lat","max", "mean", "min")

pcor_NDVI_data_lat <- matrix(pcor_data_NDVI$X7, nrow = 1440)
pcor_NDVI_data_lat_mean <- apply(pcor_NDVI_data_lat, MARGIN = 2, FUN = mean, na.rm = T)
pcor_NDVI_lat_sd <- apply(pcor_NDVI_data_lat, MARGIN = 2, FUN = sd, na.rm = T)
pcor_NDVI_lat_data <- as.data.frame(cbind(seq(90, 30, -0.25), pcor_NDVI_data_lat_mean + pcor_NDVI_lat_sd/2, pcor_NDVI_data_lat_mean, pcor_NDVI_data_lat_mean - pcor_NDVI_lat_sd/2)) 
names(pcor_NDVI_lat_data) <- c("lat","max", "mean", "min")

pcor_SOS_lat_data1 <- dplyr::filter(pcor_SOS_lat_data, lat < 65)

pcor_SOS_line_lat <- ggplot(data = pcor_SOS_lat_data1) + geom_ribbon(aes(x = lat, ymin = min, ymax = max), alpha = 0.2) + geom_line(aes(x = lat, y = mean), linewidth = 1.5) + coord_flip() +theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(colour = "black", size = 20)) +
  scale_x_continuous(limits = c(30,75), expand = c(0.01,0.01), sec.axis = sec_axis(~ ., breaks = c(30,50,70), labels = c("30N","50N","70N"))) + 
  scale_y_continuous(limits = c(-0.3,0.3), breaks = c(-0.2,-0.1,0,0.1,0.2)) +
  xlab('Latitude') +
  ylab('r') +
  theme(axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank()) + 
  theme(panel.grid = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(color = "black", size = 20), axis.title = element_text(color = "black", size = 20), axis.text.y = element_text(angle = 90), axis.ticks.length = unit(0.3, "cm"),axis.ticks = element_line(size = 1)) + 
  theme(panel.border = element_rect(color = "black", size = 1))

pcor_NDVI_lat_data1 <- dplyr::filter(pcor_NDVI_lat_data, lat < 65)

pcor_NDVI_line_lat <- ggplot(data = pcor_NDVI_lat_data1) + geom_ribbon(aes(x = lat, ymin = min, ymax = max), alpha = 0.2) + geom_line(aes(x = lat, y = mean), linewidth = 1.5) + coord_flip() +theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(colour = "black", size = 20)) +
  scale_x_continuous(limits = c(30,75), expand = c(0.01,0.01), sec.axis = sec_axis(~ ., breaks = c(30,50,70), labels = c("30N","50N","70N"))) + 
  scale_y_continuous(limits = c(-0.3,0.3), breaks = c(-0.2,-0.1,0,0.1,0.2)) +
  xlab('Latitude') +
  ylab('r') +
  theme(axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank()) + 
  theme(panel.grid = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(color = "black", size = 20), axis.title = element_text(color = "black", size = 20), axis.text.y = element_text(angle = 90), axis.ticks.length = unit(0.3, "cm"),axis.ticks = element_line(size = 1)) + 
  theme(panel.border = element_rect(color = "black", size = 1))

#-------------------------------------------------
pcor_SOS_AI_biod <- as.data.frame(cbind(pcor_data_SOS$X6, global_ai_NH_data$ai_et0, biodiversity_NH_data$biodiversity_globe_025res))

names(pcor_SOS_AI_biod) <- c('SOS','AI','biod')

pcor_SOS_AI_biod_binmean <- binMean2D(pcor_SOS_AI_biod$AI, pcor_SOS_AI_biod$biod, pcor_SOS_AI_biod$SOS, xbreaks = seq(0,21000,1000), ybreaks = seq(0,3150,150), flatten=FALSE, fill=FALSE)

pcor_SOS_AI_biod_binmean_data <- as.data.frame(cbind(rep((pcor_SOS_AI_biod_binmean$xbreaks[1:21] + 500), 21), rep((pcor_SOS_AI_biod_binmean$ybreaks[1:21] + 75), each = 21), matrix(pcor_SOS_AI_biod_binmean$result, ncol = 1)))
names(pcor_SOS_AI_biod_binmean_data) <- c('AI','biod','SOS')
pcor_SOS_AI_biod_binmean_data$AI <- pcor_SOS_AI_biod_binmean_data$AI/10000

#-------------------------------------------
pcor_SOS_AI_biod_binmean_data$SOS[1] <- 1
pcor_SOS_AI_biod_binmean_data$SOS[2] <- -1

pcor_SOS_AI_biod_bins <- ggplot(pcor_SOS_AI_biod_binmean_data) + geom_tile(aes(x = AI, y = biod, fill = SOS), alpha = 0.8) + 
  scale_fill_gradient2(low = '#0066cc', high = '#ff8000', mid = '#ffffcc', midpoint = 0, na.value = NA) +
  scale_x_continuous(limits = c(0,2), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(0,3000), expand = c(0.01,0.01)) +
  labs(fill = "r") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_text(colour = 'black', size = 20), axis.text.y = element_text(angle = 90)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), axis.ticks=element_line(color = 'black', size = 1), axis.ticks.length = unit(0.15, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1)) +
  theme(legend.direction = 'vertical', legend.spacing.x = unit(0,'cm'), legend.position = 'right',legend.key = element_rect(color="black"), legend.key.height = unit(1, 'cm'), legend.key.width = unit(0.5,'cm'), legend.text = element_text(color = 'black', size = 20), legend.title = element_text(color = 'black', size = 20, angle = 90)) + 
  guides(color = guide_legend(title.position = "right"))

pcor_SOS_AI_biod_line <- as.data.frame(cbind(seq(0.05, 2.05, 0.1), apply(pcor_SOS_AI_biod_binmean$result, FUN = mean, MARGIN = 1, na.rm = T), seq(75, 3150, 150), apply(pcor_SOS_AI_biod_binmean$result, FUN = mean, MARGIN = 2, na.rm = T)))
names(pcor_SOS_AI_biod_line) <- c('AI', 'pcor_SOS_AI', 'biod', 'pcor_SOS_biod')

pcor_SOS_AI_trend <- ggplot(pcor_SOS_AI_biod_line, aes(x=AI, y=pcor_SOS_AI)) +
  geom_point() +
  geom_smooth(method=lm, formula = y~x) +
  scale_x_continuous(limits = c(0,2), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(-0.2,0.2), expand = c(0.01,0.01), breaks = c(-0.1,0,0.1)) +
  ylab('Coincidence Rate') +
  annotate("text", label = "***", x=0.2, y=0.25, size=8)+
  theme_bw() +
  theme(axis.title.x = element_text(colour = 'black', size = 20), axis.title.y = element_blank(), axis.text = element_text(colour = 'black', size = 20), axis.text.y = element_text(angle = 90)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), axis.ticks=element_line(color = 'black', size = 1), axis.ticks.length = unit(0.15, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1)) +
  theme(legend.direction = 'vertical', legend.title = element_blank(), legend.spacing.x = unit(0,'cm'), legend.position = 'right',legend.key = element_rect(color="black"), legend.key.height = unit(3, 'cm'), legend.key.width = unit(1,'cm'), legend.text = element_text(size = 20))

pcor_SOS_biod_trend <- ggplot(pcor_SOS_AI_biod_line, aes(x=biod, y=pcor_SOS_biod)) +
  geom_point() +
  geom_smooth(method=lm, formula = y~x) +
  annotate("text", label = "***", x=500, y=-0.1, size=8, angle = 90)+
  labs(x = 'Biodiversity', y = 'r') +
  coord_flip() +
  scale_x_continuous(limits = c(0,3000), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(-0.2,0.2), expand = c(0.01,0.01), breaks = c(-0.1,0,0.1)) +
  theme_bw() +
  theme(axis.title.y = element_text(colour = 'black', size = 20), axis.title.x = element_blank(), axis.text = element_text(colour = 'black', size = 20), axis.text.y = element_text(angle = 90)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), axis.ticks=element_line(color = 'black', size = 1), axis.ticks.length = unit(0.15, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1)) +
  theme(legend.direction = 'vertical', legend.title = element_blank(), legend.spacing.x = unit(0,'cm'), legend.position = 'right',legend.key = element_rect(color="black"), legend.key.height = unit(3, 'cm'), legend.key.width = unit(1,'cm'), legend.text = element_text(size = 20))

p_AI_biod_pcor_SOS <- ggpubr::ggarrange(pcor_SOS_biod_trend, pcor_SOS_AI_biod_bins, NA, pcor_SOS_AI_trend, common.legend = T, legend = 'right', widths = c(0.5,1,NA,1), heights = c(1,0.4))

#-----------------------------------------------
#-------------------------------------------------
pcor_NDVI_AI_biod <- as.data.frame(cbind(pcor_data_NDVI$X7, global_ai_NH_data$ai_et0, biodiversity_NH_data$biodiversity_globe_025res))

names(pcor_NDVI_AI_biod) <- c('NDVI','AI','biod')

pcor_NDVI_AI_biod_binmean <- binMean2D(pcor_NDVI_AI_biod$AI, pcor_NDVI_AI_biod$biod, pcor_NDVI_AI_biod$NDVI, xbreaks = seq(0,21000,1000), ybreaks = seq(0,3150,150), flatten=FALSE, fill=FALSE)

pcor_NDVI_AI_biod_binmean_data <- as.data.frame(cbind(rep((pcor_NDVI_AI_biod_binmean$xbreaks[1:21] + 500), 21), rep((pcor_NDVI_AI_biod_binmean$ybreaks[1:21] + 75), each = 21), matrix(pcor_NDVI_AI_biod_binmean$result, ncol = 1)))
names(pcor_NDVI_AI_biod_binmean_data) <- c('AI','biod','NDVI')
pcor_NDVI_AI_biod_binmean_data$AI <- pcor_NDVI_AI_biod_binmean_data$AI/10000

#-------------------------------------------
pcor_NDVI_AI_biod_binmean_data$NDVI[1] <- 1
pcor_NDVI_AI_biod_binmean_data$NDVI[2] <- -1

pcor_NDVI_AI_biod_bins <- ggplot(pcor_NDVI_AI_biod_binmean_data) + geom_tile(aes(x = AI, y = biod, fill = NDVI), alpha = 0.8) + 
  scale_fill_gradient2(low = '#0066cc', high = '#ff8000', mid = '#ffffcc', midpoint = 0, na.value = NA) +
  scale_x_continuous(limits = c(0,2), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(0,3000), expand = c(0.01,0.01)) +
  labs(fill = "r") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_text(colour = 'black', size = 20), axis.text.y = element_text(angle = 90)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), axis.ticks=element_line(color = 'black', size = 1), axis.ticks.length = unit(0.15, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1)) +
  theme(legend.direction = 'vertical', legend.spacing.x = unit(0,'cm'), legend.position = 'right',legend.key = element_rect(color="black"), legend.key.height = unit(1, 'cm'), legend.key.width = unit(0.5,'cm'), legend.text = element_text(color = 'black', size = 20), legend.title = element_text(color = 'black', size = 20, angle = 90)) + 
  guides(color = guide_legend(title.position = "right"))

pcor_NDVI_AI_biod_line <- as.data.frame(cbind(seq(0.05, 2.05, 0.1), apply(pcor_NDVI_AI_biod_binmean$result, FUN = mean, MARGIN = 1, na.rm = T), seq(75, 3150, 150), apply(pcor_NDVI_AI_biod_binmean$result, FUN = mean, MARGIN = 2, na.rm = T)))
names(pcor_NDVI_AI_biod_line) <- c('AI', 'pcor_NDVI_AI', 'biod', 'pcor_NDVI_biod')

pcor_NDVI_AI_trend <- ggplot(pcor_NDVI_AI_biod_line, aes(x=AI, y=pcor_NDVI_AI)) +
  geom_point() +
  geom_smooth(method=lm, formula = y~x) +
  scale_x_continuous(limits = c(0,2), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(-0.2,0.2), expand = c(0.01,0.01), breaks = c(-0.1,0,0.1)) +
  ylab('Coincidence Rate') +
  annotate("text", label = "*", x=0.2, y=0.25, size=8)+
  theme_bw() +
  theme(axis.title.x = element_text(colour = 'black', size = 20), axis.title.y = element_blank(), axis.text = element_text(colour = 'black', size = 20), axis.text.y = element_text(angle = 90)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), axis.ticks=element_line(color = 'black', size = 1), axis.ticks.length = unit(0.15, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1)) +
  theme(legend.direction = 'vertical', legend.title = element_blank(), legend.spacing.x = unit(0,'cm'), legend.position = 'right',legend.key = element_rect(color="black"), legend.key.height = unit(3, 'cm'), legend.key.width = unit(1,'cm'), legend.text = element_text(size = 20))

pcor_NDVI_biod_trend <- ggplot(pcor_NDVI_AI_biod_line, aes(x=biod, y=pcor_NDVI_biod)) +
  geom_point() +
  geom_smooth(method=lm, formula = y~x) +
  annotate("text", label = "*", x=500, y=-0.1, size=8, angle = 90)+
  labs(x = 'Biodiversity', y = 'r') +
  coord_flip() +
  scale_x_continuous(limits = c(0,3000), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(-0.2,0.2), expand = c(0.01,0.01), breaks = c(-0.1,0,0.1)) +
  theme_bw() +
  theme(axis.title.y = element_text(colour = 'black', size = 20), axis.title.x = element_blank(), axis.text = element_text(colour = 'black', size = 20), axis.text.y = element_text(angle = 90)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), axis.ticks=element_line(color = 'black', size = 1), axis.ticks.length = unit(0.15, "cm")) +
  theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1)) +
  theme(legend.direction = 'vertical', legend.title = element_blank(), legend.spacing.x = unit(0,'cm'), legend.position = 'right',legend.key = element_rect(color="black"), legend.key.height = unit(3, 'cm'), legend.key.width = unit(1,'cm'), legend.text = element_text(size = 20))

p_AI_biod_pcor_NDVI <- ggpubr::ggarrange(pcor_NDVI_biod_trend, pcor_NDVI_AI_biod_bins, NA, pcor_NDVI_AI_trend, common.legend = T, legend = 'right', widths = c(0.5,1,NA,1), heights = c(1,0.4))

#--------------------------------------------------
ggarrange(pcor_SOS_spatial, pcor_NDVI_spatial, p_AI_biod_pcor_SOS, p_AI_biod_pcor_NDVI, nrow = 2, ncol = 2)

ggsave(filename = 'E:/p10.pdf', device = 'pdf', height = 14, width = 15)


# plotting SOS and kNDVI pcor at different vegetation types
#---------------------------------------------------------

landuse_2015 <- raster('E:/scientif_data/vegetation_index/MODIS/land-cover/type1/MCD12C1.A2015001.006.2018053185652_MOD12C1.tif')
landuse_2015_025res <- resample(landuse_2015, y = pcor_5.5$X1, method = 'ngb')
landuse_2015_025res_data <- raster::as.data.frame(landuse_2015_025res, xy = T)
names(landuse_2015_025res_data) <- c('x', 'y', 'layer')

pcor_SOS_kNDVI_types <- cbind(pcor_data_SOS, pcor_data_NDVI$layer, landuse_2015_025res_data$layer)
names(pcor_SOS_kNDVI_types) <- c('x','y','SOS','kNDVI','type')

pcor_SOS_kNDVI_types$type[pcor_SOS_kNDVI_types$type == 1 | pcor_SOS_kNDVI_types$type == 2] <- 'EVGR'
pcor_SOS_kNDVI_types$type[pcor_SOS_kNDVI_types$type == 3 | pcor_SOS_kNDVI_types$type == 4] <- 'DEDU'
pcor_SOS_kNDVI_types$type[pcor_SOS_kNDVI_types$type == 6 | pcor_SOS_kNDVI_types$type == 7] <- 'shrub'
pcor_SOS_kNDVI_types$type[pcor_SOS_kNDVI_types$type == 8 | pcor_SOS_kNDVI_types$type == 9] <- 'SAV'
pcor_SOS_kNDVI_types$type[pcor_SOS_kNDVI_types$type == 10] <- 'GRA'

pcor_SOS_kNDVI_types1 <- dplyr::filter(pcor_SOS_kNDVI_types, type == 'EVGR' | type == 'shrub' | type == 'DEDU' | type == 'SAV' | type == 'GRA')

pcor_SOS_types_barplot <- ggplot(data = pcor_SOS_kNDVI_types1, aes(x=type, y=SOS, fill = type)) + geom_boxplot() +
  scale_fill_manual(breaks = c('EVGR','SAV', 'shrub', 'GRA'), values = c("#00A087B2", "#3C5488B2", "#F39B7FB2", "#DC0000B2")) +
  scale_x_discrete(limits = c('EVGR','SAV', 'shrub', 'GRA')) +
  labs(y = 'Partial correlation') +
  theme_bw() + 
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(colour = 'black', size = 20), axis.text = element_text(colour = 'black', size = 20)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) + theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))

pcor_kNDVI_types_barplot <- ggplot(data = pcor_SOS_kNDVI_types1, aes(x=type, y=kNDVI, fill = type)) + geom_boxplot() +
  scale_fill_manual(breaks = c('EVGR','SAV', 'shrub', 'GRA'), values = c("#00A087B2", "#3C5488B2", "#F39B7FB2", "#DC0000B2")) +
  scale_x_discrete(limits = c('EVGR','SAV', 'shrub', 'GRA')) +
  labs(y = 'Partial correlation') +
  theme_bw() + 
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(colour = 'black', size = 20), axis.text = element_text(colour = 'black', size = 20)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.ticks.length = unit(0.3, "cm"), axis.ticks = element_line(linewidth = 1)) + theme(panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))
