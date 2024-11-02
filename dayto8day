dayto8day <- function(daily_data, method = "mean") {
  # daily_data is daily SM in the format of brick, or precipitation, etc.
  doy_start <- c(seq(1,360,8), 361)
  doy_end1 <- doy_end2 <- doy_start+7
  doy_end1[46] <- doy_end1[46] - 3
  doy_end2[46] <- doy_end2[46] - 2
  daily_data_arr <- as.array(daily_data)
  if(method == "mean"){
    data_8day <- apply(daily_data_arr[,,1:8], FUN = mean, MARGIN = c(1,2))
    for (i in 2:46){
      if (dim(daily_data_arr)[3] == 365){
        data_8day <- abind::abind(data_8day, apply(daily_data_arr[,,doy_start[i]:doy_end1[i]], FUN = mean, MARGIN = c(1,2)), along = 3)
      }
      else{
        data_8day <- abind::abind(data_8day, apply(daily_data_arr[,,doy_start[i]:doy_end2[i]], FUN = mean, MARGIN = c(1,2)), along = 3)
      }
    }
  }
  else {
    data_8day <- apply(daily_data_arr[,,1:8], FUN = sum, MARGIN = c(1,2))
    for (i in 2:46){
      if (dim(daily_data_arr)[3] == 365){
        data_8day <- abind::abind(data_8day, apply(daily_data_arr[,,doy_start[i]:doy_end1[i]], FUN = sum, MARGIN = c(1,2)), along = 3)
      }
      else{
        data_8day <- abind::abind(data_8day, apply(daily_data_arr[,,doy_start[i]:doy_end2[i]], FUN = sum, MARGIN = c(1,2)), along = 3)
      }
    }
  }
  data_8day_brick <- brick(data_8day, xmn = daily_data@extent[1], xmx = daily_data@extent[2], ymn = daily_data@extent[3], ymx = daily_data@extent[4])
  return(data_8day_brick)
}
