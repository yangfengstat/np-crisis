source('corrected/core.R')
load_inputs <- function() {
 file<-'data/External_database_withSuddenStopGI_19802019_public_09_18_2019_forecasting_h1_v12.dta'
 d<-foreign::read.dta(file,convert.factors=FALSE,convert.underscore=FALSE)
 d<-d[order(d$Country_Code,d$Year),]
 c<-foreign::read.dta('data/Calibration_CrisisIdentificationData_53EM_19802017_10_29_2018.dta',convert.factors=FALSE,convert.underscore=FALSE)
 g<-foreign::read.dta('data/Calibration_GlobalShocksData_19802017.dta',convert.factors=FALSE,convert.underscore=FALSE)
 features<-setdiff(names(d)[5:ncol(d)],c(RESERVE_PREDICTORS,'DB_GDP_BIS_FC'))
 list(data=d,crisis=c,global=g,features=features)
}
calibration_theta <- function(inputs,forecast_year,countries=NULL) {
 c<-inputs$crisis;g<-inputs$global
 if(!is.null(countries))c<-do.call(rbind,lapply(countries,function(k)c[c$Country_Code==k,]))
 c<-c[c$Year<=forecast_year-1,]
 cri<-c[c$Crisis_GI_Adjusted==1,]
 normal<-c[c$Crisis_GI_Adjusted==0,]
 gamma_all<-median(cri$d_RGDP_Growth,na.rm=TRUE)
 gamma_fall<-median(cri$d_RGDP_Growth[cri$d_RGDP_Growth<0],na.rm=TRUE)
 th<-list(sigma=2,r=.05,g=median(normal$RGDP_Growth,na.rm=TRUE)/100,
  lambda=-median(cri$d_capflp_pch_gdp,na.rm=TRUE)/100,
  gamma=-(gamma_all+gamma_fall)/200,
  pi.bar=mean(c$Crisis_GI_Adjusted,na.rm=TRUE),
  delta.op=median(g$risk_premium[g$Year<=forecast_year-1],na.rm=TRUE)/100)
 insurance_parts(th);th
}
annual_fit <- function(inputs,forecast_year,model='SE',seed=1,countries=NULL) {
 d<-inputs$data
 if(!is.null(countries))d<-do.call(rbind,lapply(countries,function(k)d[d$Country_Code==k,]))
 tr<-d[d$Year<=forecast_year-1 & !is.na(d$SuddenStop_GI),]
 # Forecast every country; future crisis labels never filter predictions.
 te<-d[d$Year==forecast_year,]
 obs<-d[d$Year==forecast_year+1,]
 features<-inputs$features
 f<-fit_np(tr[,features],tr$SuddenStop_GI,model,seed=seed)
 list(fit=f,x=as.matrix(te[,features]),theta=calibration_theta(inputs,forecast_year,countries),
      observed=median(obs$RES_GDP,na.rm=TRUE)/100,
      n_train=nrow(tr),n_crisis=sum(tr$SuddenStop_GI),n_forecast=nrow(te))
}
annual_result <- function(z,forecast_year,model,alphas=NULL,delta=.1) {
 if(is.null(alphas))alphas<-np_candidates(length(z$fit$cal_idx),delta)
 path<-alpha_path(z$fit,z$x,alphas,delta)
 path$reserve<-reserve_policy(path$share,z$theta)
 path$observed<-z$observed
 a<-argmin_intervals(alphas,path$reserve,z$observed)
 i<-match(a$alpha,path$alpha)
 summary<-data.frame(year=forecast_year,model=model,alpha=a$alpha,alpha_max=a$alpha_max,n_ties=a$n_ties,argmin=a$argmin,alpha_upper=a$alpha_upper,upper_included=a$upper_included,intervals=a$intervals,
  mse=a$loss,observed=z$observed,reserve=path$reserve[i],share=path$share[i],share_raw=path$share_raw[i],
  mode=path$mode[i],n_train=z$n_train,n_crisis=z$n_crisis,m=length(z$fit$cal_idx),n_forecast=z$n_forecast,
  d=length(z$fit$selected),alpha_min_finite=1-delta^(1/length(z$fit$cal_idx)))
 path$year<-forecast_year;path$model<-model
 list(summary=summary,path=path,theta=data.frame(year=forecast_year,as.data.frame(z$theta)))
}
