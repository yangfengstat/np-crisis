source('corrected/data.R');inp<-load_inputs();p<-inp;p$data$SuddenStop_GI[p$data$Year>=2007]<-NA
z<-annual_fit(inp,2007,'SE',seed=17);zz<-annual_fit(p,2007,'SE',seed=17)
stopifnot(z$n_forecast==53,zz$n_forecast==53,identical(z$x,zz$x),identical(z$fit,zz$fit))
p<-inp;p$crisis$RGDP_Growth[p$crisis$Year>=2007]<-10000;p$global$risk_premium[p$global$Year>=2007]<-10000
stopifnot(identical(calibration_theta(inp,2007),calibration_theta(p,2007)))
stopifnot(z$observed==median(inp$data$RES_GDP[inp$data$Year==2008],na.rm=TRUE)/100)
cat('PASS future outcomes do not filter forecasts or affect fitted scores\nPASS future history does not change insurance calibration\nPASS reserve measurement uses lag metadata and stated target year\n')
