setwd('files_path')
library(raster)
library(pracma)
library(apsimx)
library(parallel)
library(CoinCalc)
library(ranger)
library(edarf)
library(randomForest)

#---------------------------------------------
# extracting phenology (start of growing season, SOS) based on PKU GIMMS NDVI

NDVI_sub <- brick('PKU_GIMMS_NDVI_1982-2022_submonthly.nc')
NDVI_sub_arr <- as.array(NDVI_sub)
dates <- apsimx::doy2date(seq(8,365,15),1982)
for (yr in 1983:2022){
  dates <- c(dates, apsimx::doy2date(seq(8,365,15),yr))
}

NH_1growing_season <- matrix(NA, ncol = 1440, nrow = 241)
NH_TRS5_SOS_doy <- NH_TRS5_EOS_doy <- NH_DER_SOS_doy <- NH_DER_EOS_doy <- array(NA, dim = c(241, 1440, 41))

for (i in 1:241){
  for (j in 1:1440){
    
    if (sum(!is.na(NDVI_sub_arr[i, j, ])) > length(NDVI_sub_arr[i, j, ])/2){
      NDVI <- cbind(as.data.frame(dates), as.data.frame(NDVI_sub_arr[i, j, ]))
      names(NDVI) <- c('date', 'NDVI')
      input = check_input(NDVI$date, NDVI$NDVI, nptperyear = 24, south = T, maxgap = 24/4)
      #plot_input(INPUT = input)
      brks_mov <- season_mov(INPUT = input, options = list(rFUN = "smooth_wWHIT", wFUN = "wTSM", lambda = 10, r_min = 0.05, ypeak_min = 0.05, verbose = T))
      
      fits <- curvefits(INPUT = input, brks = brks_mov, options = list(
        methods = c("Beck"), #,"klos", "Gu"
        wFUN = "wTSM",nextend = 2, maxExtendMonth = 2, minExtendMonth = 1, minPercValid = 0.2
      ))
      phe <- get_pheno(fits, method = c('Beck'))
      phe_data <- phe$doy$Beck
      #phe_data = cbind(NA, phe_data)
      #names(phe_data)[1] <- "growing_period"
      # for (i in 1:dim(phe_data)[1]){
      #   phe_data$growing_period[i] <- strsplit(phe_data$flag[i], split = "_")[[1]][2]
      # }
      #phe_data_period1 <- dplyr::filter(phe_data, growing_period == 1)
      if(dim(phe_data)[1]==41){
        NH_TRS5_SOS_doy[i,j,] <- phe_data$TRS5.sos
        NH_DER_SOS_doy[i,j,] <- phe_data$DER.sos
        
        NH_TRS5_EOS_doy[i,j,] <- phe_data$TRS5.eos
        NH_DER_EOS_doy[i,j,] <- phe_data$DER.eos 
        
        NH_1growing_season[i,j] <- 1
      }
    }
  }
}

NH_TRS5_SOS_doy_mean <- apply(NH_TRS5_SOS_doy, FUN = mean, MARGIN = c(1,2), na.rm = T)
NH_TRS5_EOS_doy_mean <- apply(NH_TRS5_EOS_doy, FUN = mean, MARGIN = c(1,2), na.rm = T)
NH_DER_SOS_doy_mean <- apply(NH_DER_SOS_doy, FUN = mean, MARGIN = c(1,2), na.rm = T)
NH_DER_EOS_doy_mean <- apply(NH_DER_EOS_doy, FUN = mean, MARGIN = c(1,2), na.rm = T)

MK_trend_TRS5_SOS <- MK_pvalue_TRS5_SOS <- MK_trend_DER_SOS <- MK_pvalue_DER_SOS <- matrix(NA, nrow = 241, ncol = 1440)
MK_trend_TRS5_EOS <- MK_pvalue_TRS5_EOS <- MK_trend_DER_EOS <- MK_pvalue_DER_EOS <- matrix(NA, nrow = 241, ncol = 1440)
for (i in 1:241){
  for (j in 1:1440){
    if(!is.na(NH_TRS5_SOS_doy[i,j,1])){
      mk_test_sos <- Kendall::MannKendall(NH_TRS5_SOS_doy[i,j,])
      MK_trend_TRS5_SOS[i,j] <- mk_test_sos$tau
      MK_pvalue_TRS5_SOS[i,j] <- mk_test_sos$sl
      
      mk_test_eos <- Kendall::MannKendall(NH_TRS5_EOS_doy[i,j,])
      MK_trend_TRS5_EOS[i,j] <- mk_test_eos$tau
      MK_pvalue_TRS5_EOS[i,j] <- mk_test_eos$sl
      
      mk_test1_sos <- Kendall::MannKendall(NH_DER_SOS_doy[i,j,])
      MK_trend_DER_SOS[i,j] <- mk_test1_sos$tau
      MK_pvalue_DER_SOS[i,j] <- mk_test1_sos$sl
      
      mk_test1_eos <- Kendall::MannKendall(NH_DER_EOS_doy[i,j,])
      MK_trend_DER_EOS[i,j] <- mk_test1_eos$tau
      MK_pvalue_DER_EOS[i,j] <- mk_test1_eos$sl
    }
  }
}

#---------------------------------------------
#calculating VDP using Ta, Td, and altitude
dem025 <- raster('Global_DEM_025deg.tif')
setwd('/dewtas/')
dewtas_files <- list.files(pattern = '.nc')

cl <- makeCluster(8)
dewtas <- parLapply(cl, X = dewtas_files, fun = brick)
dewtas1 <- parLapply(cl, X = dewtas, fun = rotate)
stopCluster(cl)

setwd('/ERA5 daily data/Temperature/eightday/025res/')
tas_files <- list.files(pattern = '.nc')
cl <- makeCluster(8)
tas <- parLapply(cl, X = tas_files, fun = brick)
tas1 <- tas
stopCluster(cl)

for (i in 4:44){
  DEM <- dem025$Global_DEM_025deg
  DEM1 <- DEM
  for (j in 2:dim(tas[[i]])[3]){
    DEM1 <- stack(DEM1, DEM)
  }
  
  Pmst <- 1013.25*((273.16 + as.array(tas[[i]]) - 273.15)/(273.16 + as.array(tas[[i]]) - 273.15+0.0065*as.array(DEM1)))^5.625
  VPD <- 6.112*(1+7*10^(-4)+3.46*10^(-6)*Pmst)*exp((17.67*(as.array(tas[[i]]) - 273.15))/((as.array(tas[[i]]) - 273.15)+243.5)) - 6.112*(1+7*10^(-4)+3.46*10^(-6)*Pmst)*exp((17.67*(as.array(dewtas1[[i]]) - 273.15))/((as.array(dewtas1[[i]]) - 273.15)+243.5))
  
  VPD_ras <- brick(VPD, xmn = -180.125, xmx = 179.875, ymn = -90.125, ymx = 90.125)
  writeRaster(VPD_ras, filename = 'VPD', format = 'CDF')
}

#--------------------------------------
# calculating VPD anomaly 
setwd('/VPD/')
VPD_files <- list.files()[4:44]
cl <- makeCluster(8)
VPD_raster <- parLapply(cl, X = VPD_files, fun = brick)
VPD_1982_2022 <- VPD_raster[[1]]
for (i in 2:41){
  VPD_1982_2022 <- stack(VPD_1982_2022, VPD_raster[[i]])
}
VPD_1982_2022_arr <- as.array(VPD_1982_2022)

ma <- function(x, n = 4){filter(x, rep(1 / n, n), sides = 1)}

anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

VPD_1982_2022_arr[641:721,,] <- NA
VPD_moving_anomaly <- array(NA, dim = dim(VPD_1982_2022_arr))
for (i in 240:721) {
  for (j in 1:1440) {
    if (!is.na(VPD_1982_2022_arr[i,j,1])) {
      VPD_ij_month <- ma(VPD_1982_2022_arr[i,j,])
      VPD_month <- matrix(VPD_ij_month, nrow = 46)
      VPD_anomaly <- apply(VPD_month, FUN = anomaly, MARGIN = 1)
      VPD_moving_anomaly[i,j,] <- matrix(t(VPD_anomaly), ncol = 1)
    }
  }
}

#--------------------------------------
# preprocessing soil moisture
SM_V1_file <- list.files('V1/eightday_SM_V1/')
SM_V2_file <- list.files('V2/eightday_SM_V2/')
SM_V3_file <- list.files('V3/eightday_SM_V3/')
years <- c(1979:2022)
for (i in 1:44){
  SM_V1 <- brick(paste('V1/eightday_SM_V1/', SM_V1_file[i], sep = ''))
  SM_V2 <- brick(paste('V2/eightday_SM_V2/', SM_V2_file[i], sep = ''))
  SM_V3 <- brick(paste('V3/eightday_SM_V3/', SM_V3_file[i], sep = ''))
  SM <- 0.07*SM_V1+0.21*SM_V2+0.72*SM_V3
  writeRaster(SM, filename = paste('/SM/eightday_1m/', 'SM_1m_',years[i], sep = ''), format = 'CDF')
}

setwd('/SM/eightday_1m/')
SM_file <- list.files()
SM_1982_2022 <- brick(SM_file[4])
for (i in 5:44) {
  SM_1982_2022 <- stack(SM_1982_2022, brick(SM_file[i]))
}
SM_1982_2022_arr <- as.array(SM_1982_2022)
remove(SM_1982_2022)
gc()

ma <- function(x, n = 4){filter(x, rep(1 / n, n), sides = 1)}

anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}
SM_1982_2022_arr[641:721,,] <- NA
SM_moving_anomaly <- array(NA, dim = dim(SM_1982_2022_arr))
for (i in 1:721) {
  for (j in 1:1440) {
    if (!is.na(SM_1982_2022_arr[i,j,1])) {
      SM_ij_month <- ma(SM_1982_2022_arr[i,j,])
      # SM_stl <- stlplus::stlplus(ts(SM_ij_month, frequency = 46), s.window = 'periodic', t.window = 91, l.window = 47)
      # SM_stl_anomaly <- SM_stl$data$remainder
      SM_month <- matrix(SM_ij_month, nrow = 46)
      SM_anomaly <- apply(SM_month, FUN = anomaly, MARGIN = 1)
      
      SM_moving_anomaly[i,j,] <- matrix(t(SM_anomaly), ncol = 1)
    }
  }
}

#---------------------------------------------
# SIF preprocessing 
setwd('/GOSIF/8days/025res/')
GOSIF <- brick('GOSIF_2000057_2022361.nc')
GOSIF_arr <- as.array(GOSIF)
GOSIF_anomaly_2000_2022 <- array(NA, dim = c(721, 1440, 1058))
for (i in 1:721){
  for (j in 1:1440){
    if (!is.na(GOSIF_arr[i,j,1])){
      GOSIF_ij <- c(rep(NA,7), GOSIF_arr[i,j,])
      # GOSIF_stl <- stlplus::stlplus(ts(GOSIF_ij, frequency = 46), s.window = 'periodic', t.window = 91, l.window = 47)
      model <- lm(GOSIF_ij~c(1:1058), na.action = na.exclude)
      GOSIF_detrend <- residuals(model)
      
      # GOSIF_stl_anomaly[i,j,] <- GOSIF_stl$data$remainder
      GOSIF_month <- matrix(GOSIF_detrend, nrow = 46)
      GOSIF_anomaly <- apply(GOSIF_month, FUN = anomaly, MARGIN = 1)
      
      GOSIF_anomaly_2000_2022[i,j,] <- matrix(t(GOSIF_anomaly), ncol = 1)
    }
  }
}

#-------------------------------------------
# kNDVI pre-processing
NDVI_submonthly <- brick('/PKU_GIMMS_NDVI/NDVI_submonthly/025deg/PKU_GIMMS_NDVI_1982-2022_submonthly.nc')
NDVI_submonthly_arr <- as.array(NDVI_submonthly)
date_sub1 <- date_sub <- seq(1,365,15.5)
for (i in 1:40){
  date_sub1 <- c(date_sub1, (date_sub+(365*i)))
}
date_sub1[984] <- 14961

date_8day1 <- date_8day <- seq(1,365,8)
for (i in 1:40){
  date_8day1 <- c(date_8day1, (date_8day+(365*i)))
}

NDVI_1982_2022_8day_arr <- array(NA, dim = c(721, 1440, 1886))
for (i in 148:dim(NDVI_submonthly_arr)[1]){
  for (j in 1:dim(NDVI_submonthly_arr)[2]){
    if(sum(!is.na(NDVI_submonthly_arr[i,j,])) > 492){
      NDVI_8day <- approx(x = date_sub1, y = NDVI_submonthly_arr[i,j,], xout = date_8day1)
      NDVI_1982_2022_8day_arr[i,j,] <- NDVI_8day$y
    }
  }
}

kNDVI_1982_2022_8day <- brick('/PKU_GIMMS_NDVI/NDVI_8day/NDVI_1982_2022_8day.nc')
NDVI_arr <- as.array(kNDVI_1982_2022_8day)
kNDVI_arr <- tanh(NDVI_arr^2)
kNDVI_anomaly_1982_2022 <- array(NA, dim = c(721, 1440, 1886))
for (i in 1:721){
  for (j in 1:1440){
    if (!is.na(kNDVI_arr[i,j,19])){
      kNDVI_ij <- kNDVI_arr[i,j,]
      
      model <- lm(kNDVI_ij~c(1:1886), na.action = na.exclude)
      kNDVI_detrend <- residuals(model)
      
      kNDVI_month <- matrix(kNDVI_detrend, nrow = 46)
      kNDVI_anomaly <- apply(kNDVI_month, FUN = anomaly, MARGIN = 1)
      
      kNDVI_anomaly_1982_2022[i,j,] <- matrix(t(kNDVI_anomaly), ncol = 1)
    }
  }
}

#-----------------------------------------------
# SIF vulnerability calculation under VPD and soil moisture droughts 
SM_moving_anomaly
GOSIF_anomaly_2000_2022
SIF_anomaly_NH <- as.array(GOSIF_anomaly_2000_2022)[1:241,,]
SM_mov_anomaly_NH <- as.array(SM_moving_anomaly)[1:241,,829:1886]

SM_SIF_lag_summer_max <- raster('/SM_SIF_lag_summer_max.tif')
SM_SIF_lag_summer_max <- as.matrix(SM_SIF_lag_summer_max)

SM_SIF_vulner_summer_max <- SM_SIF_lag_summer_max <- SM_SIF_pvalue_summer_max <- matrix(NA, nrow = 241, ncol = 1440)
SM_SIF_vulner_summer_lag0_11 <- SM_SIF_pvalue_summer_lag0_11 <- array(NA, dim = c(241, 1440, 12))

SM_SIF_vulner_summer_2005_2017 <- SM_SIF_lag_summer_2005_2017 <- SM_SIF_pvalue_summer_2005_2017 <- array(NA, dim = c(241,1440,13))

for (i in 1:241){
  for (j in 1:1440){
    if(sum(!is.na(SIF_anomaly_NH[i,j,])) > 10 & sum(!is.na(SM_mov_anomaly_NH[i,j,])) > 10){
      SM <- SM_mov_anomaly_NH[i, j, ]
      
      SIF <- SIF_anomaly_NH[i, j, ]
      
      # rate <- c()
      # pvalue <- c()
      # 
      # for (t in 0:11){
      #   # different lags
      #   SM_lag <- c(rep(NA, t), SM[1:(length(SM_rolling)-t)])
      #   SM_drought_summer_lag <- CC.binarize(data = SM_lag[as.vector(matrix(c(1:1058), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
      #   SIF_depre_summer_lag = CC.binarize(data = SIF[as.vector(matrix(c(1:1058), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
      #   eca <- CoinCalc::CC.eca.ts(SIF_depre_summer_lag, SM_drought_summer_lag)
      #   rate <- c(rate, eca$`trigger coincidence rate`)
      #   pvalue <- c(pvalue, eca$`p-value trigger`)
      #   
      #   VPD_lag <- c(rep(NA, t), VPD[1:(length(VPD)-t)])
      #   VPD_drought_summer_lag <- CC.binarize(data = VPD_lag[as.vector(matrix(c(1:1058), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.9, event = 'higher')
      #   SIF_depre_summer_lag = CC.binarize(data = SIF[as.vector(matrix(c(1:1058), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
      #   eca_VPD <- CoinCalc::CC.eca.ts(SIF_depre_summer_lag, VPD_drought_summer_lag)
      #   rate_VPD <- c(rate_VPD, eca_VPD$`trigger coincidence rate`)
      #   pvalue_VPD <- c(pvalue_VPD, eca_VPD$`p-value trigger`)
      # }
      # rate_max <- max(rate, na.rm = T)
      # lag_max <- which.max(rate) + 1
      # pvalue_max <- pvalue[which.max(rate)]
      # 
      # SM_SIF_vulner_summer_max[i,j] <- rate_max
      # SM_SIF_lag_summer_max[i,j] <- lag_max
      # SM_SIF_pvalue_summer_max[i,j] <- pvalue_max
      # 
      # SM_SIF_vulner_summer_lag0_11[i,j,] <- rate
      # SM_SIF_pvalue_summer_lag0_11[i,j,] <- pvalue
      
      rate_max <- c()
      pvalue_max <- c()
      lag_max <- c()
      for (r in 1:13){
        # rolling for 11 years
        SM_rolling <- SM[((r-1)*46+1):((r+10)*46)]
        SIF_rolling <- SIF[((r-1)*46+1):((r+10)*46)]

        rate <- c()
        pvalue <- c()

        for (t in 0:11){
          # different lags
        SM_lag <- c(rep(NA, t), SM_rolling[1:(length(SM_rolling)-t)])
        SM_drought_summer_lag <- CC.binarize(data = SM_lag[as.vector(matrix(c(1:506), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
         SIF_depre_summer_lag = CC.binarize(data = SIF_rolling[as.vector(matrix(c(1:506), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
         eca <- CoinCalc::CC.eca.ts(SIF_depre_summer_lag, SM_drought_summer_lag)
         rate <- c(rate, eca$`trigger coincidence rate`)
         pvalue <- c(pvalue, eca$`p-value trigger`)
        }
         rate[is.nan(rate)] <- NA 
         rate_max <- c(rate_max, max(rate, na.rm = T))
         if(length(which.max(rate))==0){
           lag_r <- NA
           pvalue_r <- NA
         }
         else {
           lag_r <- which.max(rate)
           pvalue_r <- pvalue[which.max(rate)]
         }
         lag_max <- c(lag_max, (lag_r + 1))
         pvalue_max <- c(pvalue_max, pvalue_r)
     }

      SM_SIF_vulner_summer_2005_2017[i,j,] <- rate_max
      SM_SIF_lag_summer_2005_2017[i,j,] <- lag_max
      SM_SIF_pvalue_summer_2005_2017[i,j,] <- pvalue_max
    }
  }
}

# SM_SIF_vulner_summer_2005_2017_ras <- brick(SM_SIF_vulner_summer_2005_2017, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
# SM_SIF_pvalue_summer_2005_2017_ras <- brick(SM_SIF_pvalue_summer_2005_2017, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_SIF_lag_summer_2005_2017_ras <- brick(SM_SIF_lag_summer_2005_2017, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

writeRaster(SM_SIF_vulner_summer_2005_2017_ras, filename = '/SM_SIF_vulner_summer_2005_2017', format = 'CDF')
writeRaster(SM_SIF_pvalue_summer_2005_2017_ras, filename = '/SM_SIF_pvalue_summer_2005_2017', format = 'CDF')
writeRaster(SM_SIF_lag_summer_2005_2017_ras, filename = '/SM_SIF_lag_summer_2005_2017', format = 'CDF')

rate_mktrend_value <- rate_mktrend_pvalue <- matrix(NA, nrow = 241, ncol = 1440)
lag_mktrend_value <- lag_mktrend_pvalue <- matrix(NA, nrow = 241, ncol = 1440)
for(i in 1:241){
  for(j in 1:1440){
    if(sum(!is.na(SM_SIF_vulner_summer_2005_2017[i,j,])) > 7){
      # y = SM_SIF_vulner_summer_2005_2017[i,j,]
      # y[is.infinite(y)] <- NA
      # p = Kendall(c(1:13), y)
      # rate_mktrend_value[i,j] <- p$tau[1]
      # rate_mktrend_pvalue[i,j] <- p$sl[1]
      
      z = SM_SIF_lag_summer_2005_2017[i,j,]
      z[is.infinite(z)] <- NA
      p1 = Kendall(c(1:13), z)
      lag_mktrend_value[i,j] <- p1$tau[1]
      lag_mktrend_pvalue[i,j] <- p1$sl[1]
    }
  }
}
rate_mktrend_value_ras <- raster(rate_mktrend_value, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
rate_mktrend_pvalue_ras <- raster(rate_mktrend_pvalue, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

writeRaster(rate_mktrend_value_ras, '/SM_SIF_vulner_summer_trend', format = 'GTiff')
writeRaster(rate_mktrend_pvalue_ras, '/SM_SIF_vulner_summer_trend_pvalue', format = 'GTiff')

lag_mktrend_value_ras <- raster(lag_mktrend_value, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
lag_mktrend_pvalue_ras <- raster(lag_mktrend_pvalue, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

writeRaster(lag_mktrend_value_ras, 'E:/SM_SIF_lag_summer_trend', format = 'GTiff')
writeRaster(lag_mktrend_pvalue_ras, 'E:/SM_SIF_lag_summer_trend_pvalue', format = 'GTiff')

#--------------------------------------------
# kNDVI vulnerability calculation under soil moisture drought
SM_moving_anomaly
kNDVI_anomaly_1982_2022

kNDVI_anomaly_NH <- as.array(kNDVI_anomaly_1982_2022)[1:241,,]
SM_mov_anomaly_NH <- as.array(SM_moving_anomaly)[1:241,,]

remove(SM_moving_anomaly)
remove(kNDVI_1982_2022_anomaly_8day)
gc()

SM_kNDVI_vulner_summer_max <- SM_kNDVI_lag_summer_max <- SM_kNDVI_pvalue_summer_max <- matrix(NA, nrow = 241, ncol = 1440)
SM_kNDVI_vulner_summer_lag0_11 <- SM_kNDVI_pvalue_summer_lag0_11 <- array(NA, dim = c(241, 1440, 12))

SM_kNDVI_vulner_summer_max_first <- SM_kNDVI_lag_summer_max_first <- SM_kNDVI_pvalue_summer_max_first <- matrix(NA, nrow = 241, ncol = 1440)
SM_kNDVI_vulner_summer_lag0_11_first <- SM_kNDVI_pvalue_summer_lag0_11_first <- array(NA, dim = c(241, 1440, 12))

SM_kNDVI_vulner_summer_max_second <- SM_kNDVI_lag_summer_max_second <- SM_kNDVI_pvalue_summer_max_second <- matrix(NA, nrow = 241, ncol = 1440)
SM_kNDVI_vulner_summer_lag0_11_second <- SM_kNDVI_pvalue_summer_lag0_11_second <- array(NA, dim = c(241, 1440, 12))

SM_kNDVI_vulner_summer_1987_2017 <- SM_kNDVI_lag_summer_1987_2017 <- SM_kNDVI_pvalue_summer_1987_2017 <- array(NA, dim = c(241,1440,31))

for (i in 1:241){
  for (j in 1:1440){
    if(sum(!is.na(kNDVI_anomaly_NH[i,j,])) > 10 & sum(!is.na(SM_mov_anomaly_NH[i,j,])) > 10){
      
      SM <- SM_mov_anomaly_NH[i, j, ]
      kNDVI <- kNDVI_anomaly_NH[i, j, ]
      
      #-------------------------------------------
      # calculating vulnerability for the whole period
      rate <- rate_first <- rate_second <- c()
      pvalue <- pvalue_first <- pvalue_second <- c()
      
      for (t in 0:11){
        # different lags for whole period
        SM_lag <- c(rep(NA, t), SM[1:(length(SM)-t)])

        # different lags for first period
        SM_drought_summer_lag_first <- CC.binarize(data = SM_lag[as.vector(matrix(c(1:738), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
        kNDVI_depre_summer_lag_first = CC.binarize(data = kNDVI[as.vector(matrix(c(1:738), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
        eca_first <- CoinCalc::CC.eca.ts(kNDVI_depre_summer_lag_first, SM_drought_summer_lag_first)
        rate_first <- c(rate_first, eca_first$`trigger coincidence rate`)
        pvalue_first <- c(pvalue_first, eca_first$`p-value trigger`)
        
        # different lags for second period
        SM_drought_summer_lag_second <- CC.binarize(data = SM_lag[as.vector(matrix(c(1:738), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
        kNDVI_depre_summer_lag_second = CC.binarize(data = kNDVI[as.vector(matrix(c(1:738), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
        eca_second <- CoinCalc::CC.eca.ts(kNDVI_depre_summer_lag_second, SM_drought_summer_lag_second)
        rate_second <- c(rate_second, eca_second$`trigger coincidence rate`)
        pvalue_second <- c(pvalue_second, eca_second$`p-value trigger`)
      }
      
      rate_max <- max(rate, na.rm = T)
      lag_max <- which.max(rate) - 1
      pvalue_max <- pvalue[which.max(rate)]
      # 
      SM_kNDVI_vulner_summer_max[i,j] <- rate_max
      SM_kNDVI_lag_summer_max[i,j] <- lag_max
      SM_kNDVI_pvalue_summer_max[i,j] <- pvalue_max
      # 
      SM_kNDVI_vulner_summer_lag0_11[i,j,] <- rate
      SM_kNDVI_pvalue_summer_lag0_11[i,j,] <- pvalue
      
      rate_max_first <- max(rate_first, na.rm = T)
      lag_max_first <- which.max(rate_first) - 1
      pvalue_max_first <- pvalue_first[which.max(rate_first)]
      # 
      SM_kNDVI_vulner_summer_max_first[i,j] <- rate_max_first
      SM_kNDVI_lag_summer_max_first[i,j] <- lag_max_first
      SM_kNDVI_pvalue_summer_max_first[i,j] <- pvalue_max_first
      # 
      SM_kNDVI_vulner_summer_lag0_11_first[i,j,] <- rate_first
      SM_kNDVI_pvalue_summer_lag0_11_first[i,j,] <- pvalue_first
      
      rate_max_second <- max(rate_second, na.rm = T)
      lag_max_second <- which.max(rate_second) - 1
      pvalue_max_second <- pvalue_second[which.max(rate_second)]
      # 
      SM_kNDVI_vulner_summer_max_second[i,j] <- rate_max_second
      SM_kNDVI_lag_summer_max_second[i,j] <- lag_max_second
      SM_kNDVI_pvalue_summer_max_second[i,j] <- pvalue_max_second
      # 
      SM_kNDVI_vulner_summer_lag0_11_second[i,j,] <- rate_second
      SM_kNDVI_pvalue_summer_lag0_11_second[i,j,] <- pvalue_second
      
      #----------------------------------------------
      # calculating vulnerability for rolling peiod
      rate_max <- c()
      pvalue_max <- c()
      lag_max <- c()
      for (r in 1:31){
      # rolling for 11 years
        SM_rolling <- SM[((r-1)*46+1):((r+10)*46)]
         kNDVI_rolling <- kNDVI[((r-1)*46+1):((r+10)*46)]
      #   
         rate <- c()
         pvalue <- c()
      #   
         for (t in 0:11){
           # different lags
           SM_lag <- c(rep(NA, t), SM_rolling[1:(length(SM_rolling)-t)])
           SM_drought_summer_lag <- CC.binarize(data = SM_lag[as.vector(matrix(c(1:1426), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
           kNDVI_depre_summer_lag = CC.binarize(data = kNDVI_rolling[as.vector(matrix(c(1:1426), nrow = 46)[19:30,])], ev.def = 'percentile', thres = 0.1, event = 'lower')
           eca <- CoinCalc::CC.eca.ts(kNDVI_depre_summer_lag, SM_drought_summer_lag)
           rate <- c(rate, eca$`trigger coincidence rate`)
           pvalue <- c(pvalue, eca$`p-value trigger`)
         }
         rate[is.nan(rate)] <- NA 
         rate_max <- c(rate_max, max(rate, na.rm = T))
         if(length(which.max(rate))==0){
           lag_r <- NA
           pvalue_r <- NA
         }
         else {
           lag_r <- which.max(rate)
           pvalue_r <- pvalue[which.max(rate)]
         }
         lag_max <- c(lag_max, (lag_r + 1))
         pvalue_max <- c(pvalue_max, pvalue_r)
       }
      
       SM_kNDVI_vulner_summer_1987_2017[i,j,] <- rate_max
       SM_kNDVI_lag_summer_1987_2017[i,j,] <- lag_max
       SM_kNDVI_pvalue_summer_1987_2017[i,j,] <- pvalue_max
    }
  }
}

SM_kNDVI_vulner_summer_max_ras <- raster(SM_kNDVI_vulner_summer_max, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_lag_summer_max_ras <- raster(SM_kNDVI_lag_summer_max, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_pvalue_summer_max_ras <- raster(SM_kNDVI_pvalue_summer_max, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

SM_kNDVI_vulner_summer_lag0_11_ras <- brick(SM_kNDVI_vulner_summer_lag0_11, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_pvalue_summer_lag0_11_ras <- brick(SM_kNDVI_pvalue_summer_lag0_11, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

SM_kNDVI_vulner_summer_max_first_ras <- raster(SM_kNDVI_vulner_summer_max_first, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_lag_summer_max_first_ras <- raster(SM_kNDVI_lag_summer_max_first, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_pvalue_summer_max_first_ras <- raster(SM_kNDVI_pvalue_summer_max_first, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

SM_kNDVI_vulner_summer_lag0_11_first_ras <- brick(SM_kNDVI_vulner_summer_lag0_11_first, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_pvalue_summer_lag0_11_first_ras <- brick(SM_kNDVI_pvalue_summer_lag0_11_first, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

SM_kNDVI_vulner_summer_max_second_ras <- raster(SM_kNDVI_vulner_summer_max_second, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_lag_summer_max_second_ras <- raster(SM_kNDVI_lag_summer_max_second, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_pvalue_summer_max_second_ras <- raster(SM_kNDVI_pvalue_summer_max_second, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

SM_kNDVI_vulner_summer_1987_2017_ras <- brick(SM_kNDVI_vulner_summer_1987_2017, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_pvalue_summer_1987_2017_ras <- brick(SM_kNDVI_pvalue_summer_1987_2017, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
SM_kNDVI_lag_summer_1987_2017_ras <- brick(SM_kNDVI_lag_summer_1987_2017, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

rate_mktrend_value <- rate_mktrend_pvalue <- matrix(NA, nrow = 241, ncol = 1440)
lag_mktrend_value <- lag_mktrend_pvalue <- matrix(NA, nrow = 241, ncol = 1440)
for(i in 1:241){
  for(j in 1:1440){
    if(sum(!is.na(SM_kNDVI_vulner_summer_1987_2017[i,j,])) > 7){
      y = SM_kNDVI_vulner_summer_1987_2017[i,j,]
      y[is.infinite(y)] <- NA
      p = Kendall(c(1:31), y)
      rate_mktrend_value[i,j] <- p$tau[1]
      rate_mktrend_pvalue[i,j] <- p$sl[1]
      
      z = SM_kNDVI_lag_summer_1987_2017[i,j,]
      z[is.infinite(z)] <- NA
      p1 = Kendall(c(1:31), z)
      lag_mktrend_value[i,j] <- p1$tau[1]
      lag_mktrend_pvalue[i,j] <- p1$sl[1]
    }
  }
}
rate_mktrend_value_ras <- raster(rate_mktrend_value, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
rate_mktrend_pvalue_ras <- raster(rate_mktrend_pvalue, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

lag_mktrend_value_ras <- raster(lag_mktrend_value, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
lag_mktrend_pvalue_ras <- raster(lag_mktrend_pvalue, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

#--------------------------------------------
# influencing factors analyzing based on partial correlation and RF
SM_kNDVI_vulner_summer_1987_2017_ras
SM_kNDVI_vulner_summer_max_ras
pre_1901_2022 <- brick('/cru_ts4.07.1901.2022.pre.nc')
pre_1982_2022_NH <- resample(pre_1901_2022[[973:1464]], vulner)
pre_1982_2022_NH_arr <- as.array(pre_1982_2022_NH)

pre_1982_2022_NH_spring <- pre_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(3,4,5),])]
pre_1982_2022_NH_spring_mean <- apply(pre_1982_2022_NH_spring, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
pre_1982_2022_NH_spring_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  pre_1982_2022_NH_spring_mean1[,,i] <- pre_1982_2022_NH_spring_mean[i,,]
}

pre_1987_2017_NH_spring_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  pre_1987_2017_NH_spring_mean[,,r] <- apply(pre_1982_2022_NH_spring_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
pre_1987_2017_NH_spring_mean_ras <- brick(pre_1987_2017_NH_spring_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating spring precipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

pre_1987_2017_NH_spring_zscore <- apply(pre_1987_2017_NH_spring_mean, FUN = anomaly, MARGIN = c(1,2))
pre_1987_2017_NH_spring_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  pre_1987_2017_NH_spring_zscore1[,,i] <- pre_1987_2017_NH_spring_zscore[i,,]
}
pre_1987_2017_NH_spring_zscore_ras <- brick(pre_1987_2017_NH_spring_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(pre_1987_2017_NH_spring_zscore_ras, filename = 'E:/papers_material/vulnerability_time_changes/vulnerability_based_on_ECA/Data/pre_1987_2017_NH_spring_rolling_11year_zscore', format = 'CDF')

# for summer
pre_1982_2022_NH_summer <- pre_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(6,7,8),])]
pre_1982_2022_NH_summer_mean <- apply(pre_1982_2022_NH_summer, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
pre_1982_2022_NH_summer_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  pre_1982_2022_NH_summer_mean1[,,i] <- pre_1982_2022_NH_summer_mean[i,,]
}
pre_1982_2022_NH_summer_mean_ras <- brick(pre_1982_2022_NH_summer_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

pre_1987_2017_NH_summer_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  pre_1987_2017_NH_summer_mean[,,r] <- apply(pre_1982_2022_NH_summer_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
pre_1987_2017_NH_summer_mean_ras <- brick(pre_1987_2017_NH_summer_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating summer precipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

pre_1987_2017_NH_summer_zscore <- apply(pre_1987_2017_NH_summer_mean, FUN = anomaly, MARGIN = c(1,2))
pre_1987_2017_NH_summer_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  pre_1987_2017_NH_summer_zscore1[,,i] <- pre_1987_2017_NH_summer_zscore[i,,]
}
pre_1987_2017_NH_summer_zscore_ras <- brick(pre_1987_2017_NH_summer_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(pre_1987_2017_NH_summer_zscore_ras, filename = 'E:/papers_material/vulnerability_time_changes/vulnerability_based_on_ECA/Data/pre_1987_2017_NH_summer_rolling_11year_zscore', format = 'CDF')

#-------------------------------------
# temperature
tmp_1901_2022 <- brick('F:/Drought_in_future/CRU climate data/global/1901-2021/cru_ts4.07.1901.2022.tmp.nc')
tmp_1982_2022_NH <- resample(tmp_1901_2022[[973:1464]], vulner)
tmp_1982_2022_NH_arr <- as.array(tmp_1982_2022_NH)

tmp_1982_2022_NH_spring <- tmp_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(3,4,5),])]
tmp_1982_2022_NH_spring_mean <- apply(tmp_1982_2022_NH_spring, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
tmp_1982_2022_NH_spring_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  tmp_1982_2022_NH_spring_mean1[,,i] <- tmp_1982_2022_NH_spring_mean[i,,]
}
tmp_1982_2022_NH_spring_mean_ras <- brick(tmp_1982_2022_NH_spring_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

tmp_1987_2017_NH_spring_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  tmp_1987_2017_NH_spring_mean[,,r] <- apply(tmp_1982_2022_NH_spring_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
tmp_1987_2017_NH_spring_mean_ras <- brick(tmp_1987_2017_NH_spring_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating spring tmpcipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

tmp_1987_2017_NH_spring_zscore <- apply(tmp_1987_2017_NH_spring_mean, FUN = anomaly, MARGIN = c(1,2))
tmp_1987_2017_NH_spring_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  tmp_1987_2017_NH_spring_zscore1[,,i] <- tmp_1987_2017_NH_spring_zscore[i,,]
}
tmp_1987_2017_NH_spring_zscore_ras <- brick(tmp_1987_2017_NH_spring_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(tmp_1987_2017_NH_spring_zscore_ras, filename = 'E:/papers_material/vulnerability_time_changes/vulnerability_based_on_ECA/Data/tmp_1987_2017_NH_spring_rolling_11year_zscore', format = 'CDF', overwrite = T)

# for summer
tmp_1982_2022_NH_summer <- tmp_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(6,7,8),])]
tmp_1982_2022_NH_summer_mean <- apply(tmp_1982_2022_NH_summer, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
tmp_1982_2022_NH_summer_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  tmp_1982_2022_NH_summer_mean1[,,i] <- tmp_1982_2022_NH_summer_mean[i,,]
}
tmp_1982_2022_NH_summer_mean_ras <- brick(tmp_1982_2022_NH_summer_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

tmp_1987_2017_NH_summer_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  tmp_1987_2017_NH_summer_mean[,,r] <- apply(tmp_1982_2022_NH_summer_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
tmp_1987_2017_NH_summer_mean_ras <- brick(tmp_1987_2017_NH_summer_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating summer tmpcipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

tmp_1987_2017_NH_summer_zscore <- apply(tmp_1987_2017_NH_summer_mean, FUN = anomaly, MARGIN = c(1,2))
tmp_1987_2017_NH_summer_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  tmp_1987_2017_NH_summer_zscore1[,,i] <- tmp_1987_2017_NH_summer_zscore[i,,]
}
tmp_1987_2017_NH_summer_zscore_ras <- brick(tmp_1987_2017_NH_summer_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(tmp_1987_2017_NH_summer_zscore_ras, filename = 'E:/papers_material/vulnerability_time_changes/vulnerability_based_on_ECA/Data/tmp_1987_2017_NH_summer_rolling_11year_zscore', format = 'CDF', overwrite = T)

#-------------------------------------
# radition
rad_1981_2022 <- brick('E:/scientif_data/ERA5 data/ERA5 monthly data/climate data/025resolution/surface_net_radiation_198101_202212.nc')
rad_1982_2022_NH <- resample(rad_1981_2022[[13:504]], vulner)
rad_1982_2022_NH_arr <- as.array(rad_1982_2022_NH)

rad_1982_2022_NH_spring <- rad_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(3,4,5),])]
rad_1982_2022_NH_spring_mean <- apply(rad_1982_2022_NH_spring, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
rad_1982_2022_NH_spring_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  rad_1982_2022_NH_spring_mean1[,,i] <- rad_1982_2022_NH_spring_mean[i,,]
}
rad_1982_2022_NH_spring_mean_ras <- brick(rad_1982_2022_NH_spring_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

rad_1987_2017_NH_spring_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  rad_1987_2017_NH_spring_mean[,,r] <- apply(rad_1982_2022_NH_spring_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
rad_1987_2017_NH_spring_mean_ras <- brick(rad_1987_2017_NH_spring_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating spring radcipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

rad_1987_2017_NH_spring_zscore <- apply(rad_1987_2017_NH_spring_mean, FUN = anomaly, MARGIN = c(1,2))
rad_1987_2017_NH_spring_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  rad_1987_2017_NH_spring_zscore1[,,i] <- rad_1987_2017_NH_spring_zscore[i,,]
}
rad_1987_2017_NH_spring_zscore_ras <- brick(rad_1987_2017_NH_spring_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(rad_1987_2017_NH_spring_zscore_ras, filename = 'E:/papers_material/vulnerability_time_changes/vulnerability_based_on_ECA/Data/rad_1987_2017_NH_spring_rolling_11year_zscore', format = 'CDF')

# for summer
rad_1982_2022_NH_summer <- rad_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(6,7,8),])]
rad_1982_2022_NH_summer_mean <- apply(rad_1982_2022_NH_summer, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
rad_1982_2022_NH_summer_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  rad_1982_2022_NH_summer_mean1[,,i] <- rad_1982_2022_NH_summer_mean[i,,]
}
rad_1982_2022_NH_summer_mean_ras <- brick(rad_1982_2022_NH_summer_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

rad_1987_2017_NH_summer_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  rad_1987_2017_NH_summer_mean[,,r] <- apply(rad_1982_2022_NH_summer_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
rad_1987_2017_NH_summer_mean_ras <- brick(rad_1987_2017_NH_summer_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating summer radcipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

rad_1987_2017_NH_summer_zscore <- apply(rad_1987_2017_NH_summer_mean, FUN = anomaly, MARGIN = c(1,2))
rad_1987_2017_NH_summer_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  rad_1987_2017_NH_summer_zscore1[,,i] <- rad_1987_2017_NH_summer_zscore[i,,]
}
rad_1987_2017_NH_summer_zscore_ras <- brick(rad_1987_2017_NH_summer_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(rad_1987_2017_NH_summer_zscore_ras, filename = 'E:/papers_material/vulnerability_time_changes/vulnerability_based_on_ECA/Data/rad_1987_2017_NH_summer_rolling_11year_zscore', format = 'CDF')

#-------------------------------------
# VPD
VPD_1981_2022 <- brick('E:/scientif_data/ERA5 data/ERA5 monthly data/climate data/025resolution/VPD_025res_198101-202212.nc')
VPD_1982_2022_NH <- resample(VPD_1981_2022[[13:504]], vulner)
VPD_1982_2022_NH_arr <- as.array(VPD_1982_2022_NH)

VPD_1982_2022_NH_spring <- VPD_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(3,4,5),])]
VPD_1982_2022_NH_spring_mean <- apply(VPD_1982_2022_NH_spring, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
VPD_1982_2022_NH_spring_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  VPD_1982_2022_NH_spring_mean1[,,i] <- VPD_1982_2022_NH_spring_mean[i,,]
}
VPD_1982_2022_NH_spring_mean_ras <- brick(VPD_1982_2022_NH_spring_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

VPD_1987_2017_NH_spring_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  VPD_1987_2017_NH_spring_mean[,,r] <- apply(VPD_1982_2022_NH_spring_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
VPD_1987_2017_NH_spring_mean_ras <- brick(VPD_1987_2017_NH_spring_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating spring VPDcipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

VPD_1987_2017_NH_spring_zscore <- apply(VPD_1987_2017_NH_spring_mean, FUN = anomaly, MARGIN = c(1,2))
VPD_1987_2017_NH_spring_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  VPD_1987_2017_NH_spring_zscore1[,,i] <- VPD_1987_2017_NH_spring_zscore[i,,]
}
VPD_1987_2017_NH_spring_zscore_ras <- brick(VPD_1987_2017_NH_spring_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(VPD_1987_2017_NH_spring_zscore_ras, filename = '/Data/VPD_1987_2017_NH_spring_rolling_11year_zscore', format = 'CDF')

# for summer
VPD_1982_2022_NH_summer <- VPD_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(6,7,8),])]
VPD_1982_2022_NH_summer_mean <- apply(VPD_1982_2022_NH_summer, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
VPD_1982_2022_NH_summer_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  VPD_1982_2022_NH_summer_mean1[,,i] <- VPD_1982_2022_NH_summer_mean[i,,]
}
VPD_1982_2022_NH_summer_mean_ras <- brick(VPD_1982_2022_NH_summer_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

VPD_1987_2017_NH_summer_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  VPD_1987_2017_NH_summer_mean[,,r] <- apply(VPD_1982_2022_NH_summer_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
VPD_1987_2017_NH_summer_mean_ras <- brick(VPD_1987_2017_NH_summer_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating summer VPDcipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

VPD_1987_2017_NH_summer_zscore <- apply(VPD_1987_2017_NH_summer_mean, FUN = anomaly, MARGIN = c(1,2))
VPD_1987_2017_NH_summer_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  VPD_1987_2017_NH_summer_zscore1[,,i] <- VPD_1987_2017_NH_summer_zscore[i,,]
}
VPD_1987_2017_NH_summer_zscore_ras <- brick(VPD_1987_2017_NH_summer_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(VPD_1987_2017_NH_summer_zscore_ras, filename = '/Data/VPD_1987_2017_NH_summer_rolling_11year_zscore', format = 'CDF')

# for spring SOS and productivity
#-----------------------------------------------
TRS5_SOS_1982_2022 <- brick('/NH_TRS5_SOS_1982_2022.nc')
TRS5_SOS_1982_2022_arr <- as.array(TRS5_SOS_1982_2022)
TRS5_SOS_1987_2017 <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  TRS5_SOS_1987_2017[,,r] <- apply(TRS5_SOS_1982_2022_arr[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
TRS5_SOS_1987_2017_ras <- brick(TRS5_SOS_1987_2017, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

TRS5_SOS_1987_2017_zscore <- apply(TRS5_SOS_1987_2017, MARGIN = c(1,2), FUN = anomaly)
TRS5_SOS_1987_2017_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  TRS5_SOS_1987_2017_zscore1[,,i] <- TRS5_SOS_1987_2017_zscore[i,,]
}
TRS5_SOS_1987_2017_zscore_ras <- brick(TRS5_SOS_1987_2017_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(TRS5_SOS_1987_2017_zscore_ras, filename = '/Data/TRS5_SOS_1987_2017_NH_rolling_11year_zscore', format = 'CDF')

# spring kNDVI
kNDVI_1982_2022 <- brick('E:/scientif_data/vegetation_index/PKU_GIMMS_NDVI/kNDVI_monthly/025res/kNDVI_198201_202212.nc')
kNDVI_1982_2022_NH <- resample(kNDVI_1982_2022, vulner)
kNDVI_1982_2022_NH_arr <- as.array(kNDVI_1982_2022_NH)

kNDVI_1982_2022_NH_spring <- kNDVI_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(3,4,5),])]
kNDVI_1982_2022_NH_spring_mean <- apply(kNDVI_1982_2022_NH_spring, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
kNDVI_1982_2022_NH_spring_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  kNDVI_1982_2022_NH_spring_mean1[,,i] <- kNDVI_1982_2022_NH_spring_mean[i,,]
}
kNDVI_1982_2022_NH_spring_mean_ras <- brick(kNDVI_1982_2022_NH_spring_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

kNDVI_1987_2017_NH_spring_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  kNDVI_1987_2017_NH_spring_mean[,,r] <- apply(kNDVI_1982_2022_NH_spring_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
kNDVI_1987_2017_NH_spring_mean_ras <- brick(kNDVI_1987_2017_NH_spring_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating spring kNDVIcipitation z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

kNDVI_1987_2017_NH_spring_zscore <- apply(kNDVI_1987_2017_NH_spring_mean, FUN = anomaly, MARGIN = c(1,2))
kNDVI_1987_2017_NH_spring_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  kNDVI_1987_2017_NH_spring_zscore1[,,i] <- kNDVI_1987_2017_NH_spring_zscore[i,,]
}
kNDVI_1987_2017_NH_spring_zscore_ras <- brick(kNDVI_1987_2017_NH_spring_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(kNDVI_1987_2017_NH_spring_zscore_ras, filename = '/Data/kNDVI_1987_2017_NH_spring_rolling_11year_zscore', format = 'CDF')

# for summer
kNDVI_1982_2022_NH_summer <- kNDVI_1982_2022_NH_arr[,,as.vector(matrix(1:492, nrow = 12)[c(6,7,8),])]
kNDVI_1982_2022_NH_summer_mean <- apply(kNDVI_1982_2022_NH_summer, FUN = month2year, scal = 3, method = 'mean', MARGIN = c(1,2))
kNDVI_1982_2022_NH_summer_mean1 <- array(NA, dim = c(241,1440,41))
for (i in 1:41){
  kNDVI_1982_2022_NH_summer_mean1[,,i] <- kNDVI_1982_2022_NH_summer_mean[i,,]
}
kNDVI_1982_2022_NH_summer_mean_ras <- brick(kNDVI_1982_2022_NH_summer_mean1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

kNDVI_1987_2017_NH_summer_mean <- array(NA, dim = c(241, 1440, 31))
for (r in 1:31){
  kNDVI_1987_2017_NH_summer_mean[,,r] <- apply(kNDVI_1982_2022_NH_summer_mean1[,,r:(r+10)], FUN = mean, MARGIN = c(1,2))
}
kNDVI_1987_2017_NH_summer_mean_ras <- brick(kNDVI_1987_2017_NH_summer_mean, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)

# calculating summer z-score
anomaly <- function(x) {
  y = (x-mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

kNDVI_1987_2017_NH_summer_zscore <- apply(kNDVI_1987_2017_NH_summer_mean, FUN = anomaly, MARGIN = c(1,2))
kNDVI_1987_2017_NH_summer_zscore1 <- array(NA, dim = c(241,1440,31))
for (i in 1:31){
  kNDVI_1987_2017_NH_summer_zscore1[,,i] <- kNDVI_1987_2017_NH_summer_zscore[i,,]
}
kNDVI_1987_2017_NH_summer_zscore_ras <- brick(kNDVI_1987_2017_NH_summer_zscore1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
writeRaster(kNDVI_1987_2017_NH_summer_zscore_ras, filename = '/Data/kNDVI_1987_2017_NH_summer_rolling_11year_zscore', format = 'CDF')

# CO2 concentration variations
#---------------------------------------------
CO2_annual <- read.csv('E:/scientif_data/mlo_spo_annual_mean.csv', skip = 10, header = F)

CO2_1982_2022 <- CO2_annual[24:64,]
CO2_1987_2017 <- c()
for (i in 1:31){
  CO2_1987_2017[i] <- mean(CO2_1982_2022$V2[i:(i+10)])
}
CO2_1987_2017_zscore <- anomaly(CO2_1987_2017)

SM_kNDVI_vulner_summer_1987_2017
pre_1987_2017_NH_summer_zscore1 
tmp_1987_2017_NH_summer_zscore1 
rad_1987_2017_NH_summer_zscore1 
VPD_1987_2017_NH_summer_zscore1 
TRS5_SOS_1987_2017_zscore1 
kNDVI_1987_2017_NH_spring_zscore1 

# analyzing partial correlation 
vulner_pcor_to_7variables_1.1 <- vulner_pcor_pvalue_to_7variables_1.1 <- vulner_pcor_to_7variables_3.3 <- vulner_pcor_pvalue_to_7variables_3.3 <- vulner_pcor_to_7variables_5.5 <- vulner_pcor_pvalue_to_7variables_5.5 <- vulner_pcor_to_7variables_7.7 <- vulner_pcor_pvalue_to_7variables_7.7 <- vulner_pcor_to_7variables_9.9 <- vulner_pcor_pvalue_to_7variables_9.9 <- array(NA, dim = c(241, 1440, 7))
for (i in 1:237){
  for (j in 5:1436){
    if(!is.na(sum(TRS5_SOS_1987_2017_zscore1[i,j,])) & !is.na(sum(SM_kNDVI_vulner_summer_1987_2017[i,j,])) & !is.nan(sum(kNDVI_1987_2017_NH_spring_zscore1[i,j,])) & !is.na(sum(kNDVI_1987_2017_NH_spring_zscore1[i,j,])) & !is.nan(sum(pre_1987_2017_NH_summer_zscore1[i,j,])) & !is.na(sum(pre_1987_2017_NH_summer_zscore1[i,j,])) & !is.nan(sum(VPD_1987_2017_NH_summer_zscore1[i,j,]))){
      ppcor_data_1.1 <- as.data.frame(cbind(SM_kNDVI_vulner_summer_1987_2017[i,j,], CO2_1987_2017_zscore, pre_1987_2017_NH_summer_zscore1[i,j,], tmp_1987_2017_NH_summer_zscore1[i,j,], rad_1987_2017_NH_summer_zscore1[i,j,], VPD_1987_2017_NH_summer_zscore1[i,j,], TRS5_SOS_1987_2017_zscore1[i,j,], kNDVI_1987_2017_NH_spring_zscore1[i,j,]))
      ppcor_data_1.1 <- na.omit(ppcor_data_1.1)
      ppcor_data_3.3 <- as.data.frame(cbind(SM_kNDVI_vulner_summer_1987_2017[((i-1):(i+1)),((j-1):(j+1)),], CO2_1987_2017_zscore, pre_1987_2017_NH_summer_zscore1[((i-1):(i+1)),((j-1):(j+1)),], tmp_1987_2017_NH_summer_zscore1[((i-1):(i+1)),((j-1):(j+1)),], rad_1987_2017_NH_summer_zscore1[((i-1):(i+1)),((j-1):(j+1)),], VPD_1987_2017_NH_summer_zscore1[((i-1):(i+1)),((j-1):(j+1)),], TRS5_SOS_1987_2017_zscore1[((i-1):(i+1)),((j-1):(j+1)),], kNDVI_1987_2017_NH_spring_zscore1[((i-1):(i+1)),((j-1):(j+1)),]))
      ppcor_data_3.3 <- na.omit(ppcor_data_3.3) 
      ppcor_data_5.5 <- as.data.frame(cbind(SM_kNDVI_vulner_summer_1987_2017[((i-2):(i+2)),((j-2):(j+2)),], CO2_1987_2017_zscore, pre_1987_2017_NH_summer_zscore1[((i-2):(i+2)),((j-2):(j+2)),], tmp_1987_2017_NH_summer_zscore1[((i-2):(i+2)),((j-2):(j+2)),], rad_1987_2017_NH_summer_zscore1[((i-2):(i+2)),((j-2):(j+2)),], VPD_1987_2017_NH_summer_zscore1[((i-2):(i+2)),((j-2):(j+2)),], TRS5_SOS_1987_2017_zscore1[((i-2):(i+2)),((j-2):(j+2)),], kNDVI_1987_2017_NH_spring_zscore1[((i-2):(i+2)),((j-2):(j+2)),]))
      ppcor_data_5.5 <- na.omit(ppcor_data_5.5) 
      ppcor_data_7.7 <- as.data.frame(cbind(SM_kNDVI_vulner_summer_1987_2017[((i-3):(i+3)),((j-3):(j+3)),], CO2_1987_2017_zscore, pre_1987_2017_NH_summer_zscore1[((i-3):(i+3)),((j-3):(j+3)),], tmp_1987_2017_NH_summer_zscore1[((i-3):(i+3)),((j-3):(j+3)),], rad_1987_2017_NH_summer_zscore1[((i-3):(i+3)),((j-3):(j+3)),], VPD_1987_2017_NH_summer_zscore1[((i-3):(i+3)),((j-3):(j+3)),], TRS5_SOS_1987_2017_zscore1[((i-3):(i+3)),((j-3):(j+3)),], kNDVI_1987_2017_NH_spring_zscore1[((i-3):(i+3)),((j-3):(j+3)),]))
      ppcor_data_7.7 <- na.omit(ppcor_data_7.7) 
      ppcor_data_9.9 <- as.data.frame(cbind(SM_kNDVI_vulner_summer_1987_2017[((i-4):(i+4)),((j-4):(j+4)),], CO2_1987_2017_zscore, pre_1987_2017_NH_summer_zscore1[((i-4):(i+4)),((j-4):(j+4)),], tmp_1987_2017_NH_summer_zscore1[((i-4):(i+4)),((j-4):(j+4)),], rad_1987_2017_NH_summer_zscore1[((i-4):(i+4)),((j-4):(j+4)),], VPD_1987_2017_NH_summer_zscore1[((i-4):(i+4)),((j-4):(j+4)),], TRS5_SOS_1987_2017_zscore1[((i-4):(i+4)),((j-4):(j+4)),], kNDVI_1987_2017_NH_spring_zscore1[((i-4):(i+4)),((j-4):(j+4)),]))
      ppcor_data_9.9 <- na.omit(ppcor_data_9.9) 
      names(ppcor_data_1.1) <- c('vulner', 'CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')
      names(ppcor_data_3.3) <- c('vulner', 'CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')
      names(ppcor_data_5.5) <- c('vulner', 'CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')
      names(ppcor_data_7.7) <- c('vulner', 'CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')
      names(ppcor_data_9.9) <- c('vulner', 'CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')
      
      pcor_data_1.1 <- pcor(ppcor_data_1.1)
      vulner_pcor_to_7variables_1.1[i,j,] <- pcor_data_1.1$estimate[2:8]
      vulner_pcor_pvalue_to_7variables_1.1[i,j,] <- pcor_data_1.1$p.value[2:8]
      
      pcor_data_3.3 <- pcor(ppcor_data_3.3)
      vulner_pcor_to_7variables_3.3[i,j,] <- pcor_data_3.3$estimate[2:8]
      vulner_pcor_pvalue_to_7variables_3.3[i,j,] <- pcor_data_3.3$p.value[2:8]
      
      pcor_data_5.5 <- pcor(ppcor_data_5.5)
      vulner_pcor_to_7variables_5.5[i,j,] <- pcor_data_5.5$estimate[2:8]
      vulner_pcor_pvalue_to_7variables_5.5[i,j,] <- pcor_data_5.5$p.value[2:8]
      
      pcor_data_7.7 <- pcor(ppcor_data_7.7)
      vulner_pcor_to_7variables_7.7[i,j,] <- pcor_data_7.7$estimate[2:8]
      vulner_pcor_pvalue_to_7variables_7.7[i,j,] <- pcor_data_7.7$p.value[2:8]
      
      pcor_data_9.9 <- pcor(ppcor_data_9.9)
      vulner_pcor_to_7variables_9.9[i,j,] <- pcor_data_9.9$estimate[2:8]
      vulner_pcor_pvalue_to_7variables_9.9[i,j,] <- pcor_data_9.9$p.value[2:8]
    }
  }
}
vulner_pcor_to_7variables_1.1_ras <- brick(vulner_pcor_to_7variables_1.1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_to_7variables_1.1_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

vulner_pcor_pvalue_to_7variables_1.1_ras <- brick(vulner_pcor_pvalue_to_7variables_1.1, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_pvalue_to_7variables_1.1_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

writeRaster(vulner_pcor_to_7variables_1.1_ras, filename = '/vulner_pcor_to_7variables_1.1', format = 'CDF')
writeRaster(vulner_pcor_pvalue_to_7variables_1.1_ras, filename = '/vulner_pcor_pvalue_to_7variables_1.1', format = 'CDF')

vulner_pcor_to_7variables_3.3_ras <- brick(vulner_pcor_to_7variables_3.3, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_to_7variables_3.3_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

vulner_pcor_pvalue_to_7variables_3.3_ras <- brick(vulner_pcor_pvalue_to_7variables_3.3, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_pvalue_to_7variables_3.3_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

writeRaster(vulner_pcor_to_7variables_3.3_ras, filename = '/vulner_pcor_to_7variables_3.3', format = 'CDF')
writeRaster(vulner_pcor_pvalue_to_7variables_3.3_ras, filename = '/vulner_pcor_pvalue_to_7variables_3.3', format = 'CDF')

vulner_pcor_to_7variables_5.5_ras <- brick(vulner_pcor_to_7variables_5.5, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_to_7variables_5.5_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

vulner_pcor_pvalue_to_7variables_5.5_ras <- brick(vulner_pcor_pvalue_to_7variables_5.5, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_pvalue_to_7variables_5.5_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

writeRaster(vulner_pcor_to_7variables_5.5_ras, filename = '/vulner_pcor_to_7variables_5.5', format = 'CDF')
writeRaster(vulner_pcor_pvalue_to_7variables_5.5_ras, filename = '/vulner_pcor_pvalue_to_7variables_5.5', format = 'CDF')

vulner_pcor_to_7variables_7.7_ras <- brick(vulner_pcor_to_7variables_7.7, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_to_7variables_7.7_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

vulner_pcor_pvalue_to_7variables_7.7_ras <- brick(vulner_pcor_pvalue_to_7variables_7.7, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_pvalue_to_7variables_7.7_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

writeRaster(vulner_pcor_to_7variables_7.7_ras, filename = '/vulner_pcor_to_7variables_7.7', format = 'CDF')
writeRaster(vulner_pcor_pvalue_to_7variables_7.7_ras, filename = '/vulner_pcor_pvalue_to_7variables_7.7', format = 'CDF')

vulner_pcor_to_7variables_9.9_ras <- brick(vulner_pcor_to_7variables_9.9, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_to_7variables_9.9_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

vulner_pcor_pvalue_to_7variables_9.9_ras <- brick(vulner_pcor_pvalue_to_7variables_9.9, xmn = -180.125, xmx = 179.875, ymn = 29.875, ymx = 90.125)
names(vulner_pcor_pvalue_to_7variables_9.9_ras) <- c('CO2', 'pre_summer', 'tmp_summer', 'rad_summer', 'VPD_summer', 'SOS', 'kNDVI_spring')

writeRaster(vulner_pcor_to_7variables_9.9_ras, filename = '/vulner_pcor_to_7variables_9.9', format = 'CDF')
writeRaster(vulner_pcor_pvalue_to_7variables_9.9_ras, filename = '/vulner_pcor_pvalue_to_7variables_9.9', format = 'CDF')

# random forest modeling and partial dependence
#-----------------------------------------
landuse <- raster('/MCD12C1.A2015001.006.2018053185652_MOD12C1.tif')
landuse_025res <- resample(landuse, y = vulner, method = 'ngb')
landuse_data <- raster::as.data.frame(landuse_025res, xy = T)
names(landuse_data) <- c('x', 'y', 'layer')
landuse_data <- do.call("rbind", replicate(31, landuse_data, simplify = FALSE))

vulner_data <- raster::as.data.frame(SM_kNDVI_vulner_summer_1987_2017_ras[[1]], xy = T)
names(vulner_data) <- c('x','y','layer')
for (i in 2:31){
  vulner_i <- raster::as.data.frame(SM_kNDVI_vulner_summer_1987_2017_ras[[i]], xy = T)
  names(vulner_i) <- c('x','y','layer')
  vulner_data <- rbind(vulner_data, vulner_i)
}

CO2_data <- vulner_data
CO2_data$layer <- rep(CO2_1987_2017_zscore, each = 241*1440)

pre_spring_data <- raster::as.data.frame(pre_1987_2017_NH_spring_zscore_ras[[1]], xy = T)
names(pre_spring_data) <- c('x','y','layer')
for (i in 2:31){
  pre_spring_i <- raster::as.data.frame(pre_1987_2017_NH_spring_zscore_ras[[i]], xy = T)
  names(pre_spring_i) <- c('x','y','layer')
  pre_spring_data <- rbind(pre_spring_data, pre_spring_i)
}

pre_summer_data <- raster::as.data.frame(pre_1987_2017_NH_summer_zscore_ras[[1]], xy = T)
names(pre_summer_data) <- c('x','y','layer')
for (i in 2:31){
  pre_summer_i <- raster::as.data.frame(pre_1987_2017_NH_summer_zscore_ras[[i]], xy = T)
  names(pre_summer_i) <- c('x','y','layer')
  pre_summer_data <- rbind(pre_summer_data, pre_summer_i)
}

tmp_spring_data <- raster::as.data.frame(tmp_1987_2017_NH_spring_zscore_ras[[1]], xy = T)
names(tmp_spring_data) <- c('x','y','layer')
for (i in 2:31){
  tmp_spring_i <- raster::as.data.frame(tmp_1987_2017_NH_spring_zscore_ras[[i]], xy = T)
  names(tmp_spring_i) <- c('x','y','layer')
  tmp_spring_data <- rbind(tmp_spring_data, tmp_spring_i)
}

tmp_summer_data <- raster::as.data.frame(tmp_1987_2017_NH_summer_zscore_ras[[1]], xy = T)
names(tmp_summer_data) <- c('x','y','layer')
for (i in 2:31){
  tmp_summer_i <- raster::as.data.frame(tmp_1987_2017_NH_summer_zscore_ras[[i]], xy = T)
  names(tmp_summer_i) <- c('x','y','layer')
  tmp_summer_data <- rbind(tmp_summer_data, tmp_summer_i)
}

rad_spring_data <- raster::as.data.frame(rad_1987_2017_NH_spring_zscore_ras[[1]], xy = T)
names(rad_spring_data) <- c('x','y','layer')
for (i in 2:31){
  rad_spring_i <- raster::as.data.frame(rad_1987_2017_NH_spring_zscore_ras[[i]], xy = T)
  names(rad_spring_i) <- c('x','y','layer')
  rad_spring_data <- rbind(rad_spring_data, rad_spring_i)
}

rad_summer_data <- raster::as.data.frame(rad_1987_2017_NH_summer_zscore_ras[[1]], xy = T)
names(rad_summer_data) <- c('x','y','layer')
for (i in 2:31){
  rad_summer_i <- raster::as.data.frame(rad_1987_2017_NH_summer_zscore_ras[[i]], xy = T)
  names(rad_summer_i) <- c('x','y','layer')
  rad_summer_data <- rbind(rad_summer_data, rad_summer_i)
}

VPD_spring_data <- raster::as.data.frame(VPD_1987_2017_NH_spring_zscore_ras[[1]], xy = T)
names(VPD_spring_data) <- c('x','y','layer')
for (i in 2:31){
  VPD_spring_i <- raster::as.data.frame(VPD_1987_2017_NH_spring_zscore_ras[[i]], xy = T)
  names(VPD_spring_i) <- c('x','y','layer')
  VPD_spring_data <- rbind(VPD_spring_data, VPD_spring_i)
}

VPD_summer_data <- raster::as.data.frame(VPD_1987_2017_NH_summer_zscore_ras[[1]], xy = T)
names(VPD_summer_data) <- c('x','y','layer')
for (i in 2:31){
  VPD_summer_i <- raster::as.data.frame(VPD_1987_2017_NH_summer_zscore_ras[[i]], xy = T)
  names(VPD_summer_i) <- c('x','y','layer')
  VPD_summer_data <- rbind(VPD_summer_data, VPD_summer_i)
}

SOS_data <- raster::as.data.frame(TRS5_SOS_1987_2017_zscore_ras[[1]], xy = T)
names(SOS_data) <- c('x','y','layer')
for (i in 2:31){
  SOS_i <- raster::as.data.frame(TRS5_SOS_1987_2017_zscore_ras[[i]], xy = T)
  names(SOS_i) <- c('x','y','layer')
  SOS_data <- rbind(SOS_data, SOS_i)
}

kNDVI_spring_data <- raster::as.data.frame(kNDVI_1987_2017_NH_spring_zscore_ras[[1]], xy = T)
names(kNDVI_spring_data) <- c('x','y','layer')
for (i in 2:31){
  kNDVI_spring_i <- raster::as.data.frame(kNDVI_1987_2017_NH_spring_zscore_ras[[i]], xy = T)
  names(kNDVI_spring_i) <- c('x','y','layer')
  kNDVI_spring_data <- rbind(kNDVI_spring_data, kNDVI_spring_i)
}

kNDVI_summer_data <- raster::as.data.frame(kNDVI_1987_2017_NH_summer_zscore_ras[[1]], xy = T)
names(kNDVI_summer_data) <- c('x','y','layer')
for (i in 2:31){
  kNDVI_summer_i <- raster::as.data.frame(kNDVI_1987_2017_NH_summer_zscore_ras[[i]], xy = T)
  names(kNDVI_summer_i) <- c('x','y','layer')
  kNDVI_summer_data <- rbind(kNDVI_summer_data, kNDVI_summer_i)
}

biodiversity_globle <- raster('/biodiversity_globe_025res.tif')
biodiversity_NH <- resample(biodiversity_globle, y = vulner)
biodiversity_NH_data <- raster::as.data.frame(biodiversity_NH, xy = T)
biodiversity_data <- do.call("rbind", replicate(31, biodiversity_NH_data, simplify = FALSE))

global_ai <- raster('/ai_et0.tif')
global_ai_NH <- resample(global_ai, y = vulner)
global_ai_NH_data <- raster::as.data.frame(global_ai_NH$ai_et0, xy = T)
global_ai_data <- do.call("rbind", replicate(31, global_ai_NH_data, simplify = FALSE))

RF_data <- cbind(vulner_data, CO2_data$layer, pre_spring_data$layer, pre_summer_data$layer, tmp_spring_data$layer, tmp_summer_data$layer, rad_spring_data$layer, rad_summer_data$layer, VPD_spring_data$layer, VPD_summer_data$layer, SOS_data$layer, kNDVI_spring_data$layer, kNDVI_summer_data$layer, biodiversity_data$biodiversity_globe_025res, global_ai_data$ai_et0, landuse_data$layer)
names(RF_data) <- c('long','lat','vulner', 'CO2', 'pre_spring', 'pre_summer', 'tmp_spring', 'tmp_summer', 'rad_spring', 'rad_summer', 'VPD_spring', 'VPD_summer', 'SOS', 'kNDVI_spring', 'kNDVI_summer','biodiv','AI','landuse')
save(RF_data, file = 'kNDVI_vulnerability_RF_data.RData')
