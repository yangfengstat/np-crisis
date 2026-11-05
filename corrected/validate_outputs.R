# Validate completed release outputs without refitting models or trusting caches.
read<-function(path)read.csv(file.path('outputs_corrected',path),stringsAsFactors=FALSE)
p<-read('empirical/annual_estimates.csv')
stopifnot(nrow(p)==50,identical(sort(unique(p$year)),1993:2017),all(table(p$model)==25),
          all(is.finite(p$alpha)),all(p$alpha>=0 & p$alpha<=1),all(p$reserve>=0),
          all(p$alpha<=p$alpha_upper),all(p$alpha_upper<=1),all(p$n_forecast==53),
          all(p$share>=.001 & p$share<=.999))
for(model in c('SE','ada')){
 B<-if(model=='SE')500L else 200L
 r<-read(paste0('empirical/bootstrap_',model,'_raw.csv'))
 s<-read(paste0('empirical/bootstrap_',model,'_summary.csv'))
 stopifnot(nrow(r)==25*B,length(unique(r$rep))==B,all(table(r$rep)==25),
           all(table(r$year)==B),anyDuplicated(r[c('rep','year')])==0,nrow(s)==25,all(s$B==B),all(r$model==model),
           all(is.finite(r$alpha)),all(is.finite(r$reserve)),all(r$reserve>=0),
           all(r$alpha>=0 & r$alpha<=1),all(r$alpha<=r$alpha_upper),all(r$alpha_upper<=1))
 for(t in 1993:2017)stopifnot(abs(s$sd[s$year==t]-sd(r$alpha[r$year==t]))<1e-12)
 cat('PASS',model,B,'complete bootstrap trajectories and summary agreement\n')
}
b<-read('backtest/metrics.csv')
stopifnot(nrow(b)==35,all(table(b$model)==7),all(b$delta==.1),all(b$m==25),
          all(b$np_auc>=0 & b$np_auc<=1),all(b$np_pauc>=0 & b$np_pauc<=.25),
          all(b$roc_auc>=0 & b$roc_auc<=1))
for(n in c(1000,2000,4000)){
 a<-read(paste0('simulations/panelA_n',n,'_raw.csv'))
 stopifnot(nrow(a)==42000,length(unique(a$rep))==1000,all(a$R0>=-1e-12 & a$R0<=1+1e-12),all(a$R1>=-1e-12 & a$R1<=1+1e-12))
}
w<-read('simulations/panelB2_raw.csv')
stopifnot(nrow(w)==42000,length(unique(w$rep))==1000,all(w$loss>=-1e-12))
r<-read('simulations/panelCD_raw.csv')
stopifnot(nrow(r)==18750,sum(r$scenario=='baseline')==7500,all(table(r$scenario)[names(table(r$scenario))!='baseline']==3750))
for(f in c('backtest/figure2_corrected.pdf','empirical/figure3_corrected.pdf','empirical/alpha_model_comparison_corrected.pdf','simulations/panelA_error_control.pdf','simulations/panelB_welfare.pdf','simulations/panelC_alpha_recovery.pdf','simulations/panelD_recovery_robustness.pdf'))stopifnot(file.exists(file.path('outputs_corrected',f)))
cat('PASS all backtest, simulation, figure, and empirical output dimensions\n')
cat('ALL RELEASE OUTPUT CHECKS PASSED\n')
