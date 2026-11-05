# Corrected simulations, base R only. Run from code/ with --vanilla.
# Rscript --vanilla corrected/simulations.R full
# Population errors below integrate the actual fixed SE score distribution.
source('corrected/core.R')

simulation_config <- function(mode='full') {
 full_counts<-mode%in%c('full','panelA_repair')
 list(mode=mode,seed=20261015L,pi0=82/1932,d=30L,s=5L,auc=.70,
      alphas=seq(0,1,.05),delta=.1,ns=c(1000L,2000L,4000L),
      B_A=if(full_counts)1000L else 3L,B_B=if(full_counts)1000L else 3L,
      B_C=if(full_counts)300L else 2L,B_D=if(full_counts)150L else 2L,
      n_train=2000L,n_eval=53L,tau=.03,cores=4L,
      output=if(full_counts)'outputs_corrected/simulations'else paste0('outputs_corrected/simulations/',mode))
}
sim_parameters <- function(cfg) {
 mu<-c(rep(sqrt(2)*qnorm(cfg$auc)/sqrt(cfg$s),cfg$s),rep(0,cfg$d-cfg$s))
 names(mu)<-paste0('x',seq_len(cfg$d))
 list(mu1=mu,d=cfg$d,s=cfg$s,pi0=cfg$pi0)
}
sim_theta <- function(cfg)list(sigma=2,r=.05,g=.034,lambda=.069,gamma=.065,
                             pi.bar=cfg$pi0,delta.op=.011)
draw_gaussian <- function(n,pars,prevalence=pars$pi0) {
 y<-rbinom(n,1,prevalence);x<-matrix(rnorm(n*pars$d),n,pars$d)
 x<-x+outer(y,pars$mu1);colnames(x)<-names(pars$mu1)
 list(x=x,y=y)
}
bernoulli_half <- function(weights,prob) {
 value<-0;mass<-1
 for(j in seq_along(weights)) {
  value<-c(value,value+weights[j]);mass<-c(mass*(1-prob[j]),mass*prob[j])
 }
 o<-order(value);list(value=value[o],mass=mass[o],id=(seq_along(value)-1L)[o])
}
weighted_score_cdf <- function(thresholds,weights,prob,tol=64*.Machine$double.eps) {
 stopifnot(length(weights)==length(prob),all(weights>=0),all(prob>=0 & prob<=1))
 full_weights<-weights;keep<-which(weights>0);weights<-weights[keep];prob<-prob[keep]
 if(!length(weights))return(as.numeric(thresholds>0))
 h<-floor(length(weights)/2);left<-seq_len(h);right<-setdiff(seq_along(weights),left)
 a<-bernoulli_half(weights[left],prob[left]);b<-bernoulli_half(weights[right],prob[right])
 cumulative<-c(0,cumsum(b$mass));eps<-tol*max(1,sum(weights))
 # Epsilon only identifies ambiguous boundary patterns. Reconstruct these
 # patterns in the full original column order and apply the same literal
 # matrix score and strict comparison as core, with no threshold tolerance.
 vapply(thresholds,function(thr){
  if(thr==-Inf)return(0);if(thr==Inf)return(1)
  lo<-findInterval(thr-a$value-eps,b$value);hi<-findInterval(thr-a$value+eps,b$value)
  result<-sum(a$mass*cumulative[lo+1L]);ii<-which(hi>lo)
  if(!length(ii))return(result)
  ai<-rep(ii,hi[ii]-lo[ii]);bi<-unlist(lapply(ii,function(i)seq.int(lo[i]+1L,hi[i])),use.names=FALSE)
  for(start in seq.int(1L,length(ai),by=10000L)) {
   jj<-seq.int(start,min(start+9999L,length(ai)))
   bits<-matrix(0,length(jj),length(full_weights))
   for(j in left)bits[,keep[j]]<-as.integer(bitwAnd(a$id[ai[jj]],bitwShiftL(1L,j-1L))!=0L)
   for(j in seq_along(right))bits[,keep[right[j]]]<-as.integer(bitwAnd(b$id[bi[jj]],bitwShiftL(1L,j-1L))!=0L)
   values<-drop(bits%*%full_weights)
   result<-result+sum(a$mass[ai[jj]]*b$mass[bi[jj]]*(values<thr))
  }
  result
 },numeric(1))
}
se_population_errors <- function(f,thresholds,pars) {
 stopifnot(f$model=='SE')
 direction<-f$base$direction;cut<-f$base$cut;mu<-pars$mu1[f$selected]
 p0<-ifelse(direction==1,pnorm(cut,lower.tail=FALSE),pnorm(cut))
 p1<-ifelse(direction==1,pnorm(cut-mu,lower.tail=FALSE),pnorm(cut-mu))
 data.frame(R0=1-weighted_score_cdf(thresholds,f$base$weight,p0),
            R1=weighted_score_cdf(thresholds,f$base$weight,p1))
}
panel_a_rep <- function(b,n,alphas,cfg) {
 seed<-cfg$seed+n*10000L+b;set.seed(seed);pars<-sim_parameters(cfg)
 dat<-draw_gaussian(n,pars);f<-fit_np(dat$x,dat$y,'SE',seed=seed+1L)
 rr<-lapply(alphas,function(a)np_rule(f$cal_scores,a,cfg$delta))
 np_thr<-vapply(rr,`[[`,0,'threshold')
 plugin<-vapply(alphas,function(a)if(a==0)-Inf else if(a==1)Inf else
                    as.numeric(quantile(f$cal_scores,a,type=1)),0)
 cuts<-c(np_thr,plugin);er<-se_population_errors(f,cuts,pars)
 data.frame(rep=b,n=n,method=rep(c('NP umbrella','plug-in quantile'),each=length(alphas)),
            alpha=rep(alphas,2),m=length(f$cal_scores),threshold=cuts,
            mode=c(vapply(rr,`[[`,'','mode'),rep('plugin',length(alphas))),
            R0=er$R0,R1=er$R1,violation=er$R1>rep(alphas,2)+1e-12)
}
summarize_a <- function(raw)do.call(rbind,lapply(split(raw,list(raw$n,raw$method,raw$alpha),drop=TRUE),function(z){
 ci<-binom.test(sum(z$violation),nrow(z))$conf.int
 data.frame(n=z$n[1],method=z$method[1],alpha=z$alpha[1],B=nrow(z),
            R0=mean(z$R0),R1=mean(z$R1),R1_sd=sd(z$R1),violation=mean(z$violation),
            violation_low=ci[1],violation_high=ci[2],
            fallback=mean(z$mode=='all_warning'),mean_m=mean(z$m))
}))
panel_b_exact <- function(cfg) {
 th<-sim_theta(cfg)
 z<-expand.grid(p=seq(.03,.30,.005),relative_error=c(.25,.5,.75))
 z$q_down<-z$p*(1-z$relative_error);z$q_up<-z$p*(1+z$relative_error)
 z$loss_down<-policy_loss(z$p,z$q_down,th);z$loss_up<-policy_loss(z$p,z$q_up,th)
 z$ratio<-ifelse(z$loss_up>1e-14,z$loss_down/z$loss_up,NA_real_)
 z$rho_truth<-reserve_policy(z$p,th);z$rho_down<-reserve_policy(z$q_down,th)
 z$rho_up<-reserve_policy(z$q_up,th);z
}
panel_b_rep <- function(b,cfg) {
 seed<-cfg$seed+50000000L+b;set.seed(seed);pars<-sim_parameters(cfg);th<-sim_theta(cfg)
 dat<-draw_gaussian(cfg$n_train,pars);f<-fit_np(dat$x,dat$y,'SE',seed=seed+1L)
 do.call(rbind,lapply(c(cfg$pi0,.20),function(p){
  ev<-draw_gaussian(cfg$n_eval,pars,p);path<-alpha_path(f,ev$x,cfg$alphas,cfg$delta)
  data.frame(rep=b,p=p,alpha=path$alpha,warning_share=path$share,
             loss=policy_loss(p,path$share,th),rho=reserve_policy(path$share,th))
 }))
}
summarize_b <- function(raw)do.call(rbind,lapply(split(raw,list(raw$p,raw$alpha),drop=TRUE),function(z){
 data.frame(p=z$p[1],alpha=z$alpha[1],B=nrow(z),mean_share=mean(z$warning_share),
            expected_loss=mean(z$loss),loss_se=sd(z$loss)/sqrt(nrow(z)),mean_rho=mean(z$rho))
}))
true_alpha <- function(periods)approx(c(1,6,10,18,25),c(.50,.90,.45,.90,.40),xout=periods)$y
in_minimizing_intervals <- function(alphas,reserves,observed,truth,tolerance=1e-10) {
 loss<-(reserves-observed)^2;best<-min(loss)
 ii<-which(abs(loss-best)<=tolerance*max(1,abs(best)))
 upper<-c(tail(alphas,-1),1)
 any((truth>=alphas[ii]-1e-12 & truth<upper[ii]-1e-12) | (alphas[ii]==1 & truth==1))
}
recovery_rep <- function(b,cfg,periods=1:25,scenarios=c('baseline','r_known','r_ignored','calibration_jitter')) {
 seed<-cfg$seed+70000000L+b*1000L;pars<-sim_parameters(cfg);th<-sim_theta(cfg)
 do.call(rbind,lapply(periods,function(t){
  # Period seeds make outputs invariant to worker count or scenario selection.
  set.seed(seed+t*10L);ev<-draw_gaussian(cfg$n_eval,pars);dp<-draw_gaussian(cfg$n_train,pars)
  da<-draw_gaussian(cfg$n_train,pars)
  fp<-fit_np(dp$x,dp$y,'SE',seed=seed+t*10L+1L)
  fa<-fit_np(da$x,da$y,'SE',seed=seed+t*10L+2L)
  pp<-alpha_path(fp,ev$x,np_candidates(length(fp$cal_scores),cfg$delta),cfg$delta)
  pa<-alpha_path(fa,ev$x,np_candidates(length(fa$cal_scores),cfg$delta),cfg$delta)
  atrue<-true_alpha(t);qpol<-alpha_path(fp,ev$x,atrue,cfg$delta)$share
  set.seed(seed+t*10L+3L);noise<-1+rnorm(1,0,cfg$tau)
  jit<-runif(4,.85,1.15);names(jit)<-c('g','lambda','gamma','delta.op')
  do.call(rbind,lapply(scenarios,function(sc){
   tg<-tr<-th
   if(sc%in%c('r_known','r_ignored'))tg$r<-.03+.04*(t-1)/24
   if(sc=='r_known')tr$r<-tg$r
   if(sc=='calibration_jitter')for(nm in names(jit))tg[[nm]]<-th[[nm]]*jit[[nm]]
   rho_true<-reserve_policy(qpol,tg);observed<-max(0,rho_true*noise)
   candidate<-reserve_policy(pa$share,tr);ch<-argmin_intervals(pa$alpha,candidate,observed)
   oracle<-argmin_intervals(pp$alpha,reserve_policy(pp$share,tg),rho_true)
   data.frame(rep=b,t=t,scenario=sc,alpha_true=atrue,alpha_hat=ch$alpha,
              alpha_max=ch$alpha_max,alpha_upper=ch$alpha_upper,n_ties=ch$n_ties,argmin=ch$argmin,
              minimizing_intervals=ch$intervals,candidate_count=nrow(pa),
              rho_observed=observed,rho_true=rho_true,share_policy=qpol,
              distinct_reserves=length(unique(round(candidate,12))),
              zero_reserve_candidates=sum(candidate==0),
              oracle_n_ties=oracle$n_ties,oracle_alpha_min=oracle$alpha,
              oracle_alpha_max=oracle$alpha_max,oracle_alpha_upper=oracle$alpha_upper,
              true_in_argmin=in_minimizing_intervals(pa$alpha,candidate,observed,atrue),
              matched_squared_error=ch$loss)
  }))
 }))
}
summarize_recovery <- function(raw)do.call(rbind,lapply(split(raw,list(raw$scenario,raw$t),drop=TRUE),function(z){
 data.frame(scenario=z$scenario[1],t=z$t[1],alpha_true=z$alpha_true[1],B=nrow(z),
            alpha_hat_mean=mean(z$alpha_hat),alpha_hat_sd=sd(z$alpha_hat),
            alpha_hat_q025=unname(quantile(z$alpha_hat,.025)),alpha_hat_q975=unname(quantile(z$alpha_hat,.975)),
            bias=mean(z$alpha_hat-z$alpha_true),rmse=sqrt(mean((z$alpha_hat-z$alpha_true)^2)),
            tie_rate=mean(z$n_ties>1),mean_ties=mean(z$n_ties),
            mean_argmin_width=mean(z$alpha_upper-z$alpha_hat),
            true_in_argmin=mean(z$true_in_argmin),mean_distinct_reserves=mean(z$distinct_reserves),
            oracle_tie_rate=mean(z$oracle_n_ties>1),oracle_mean_ties=mean(z$oracle_n_ties),
            oracle_bias=mean(z$oracle_alpha_min-z$alpha_true),
            oracle_rmse=sqrt(mean((z$oracle_alpha_min-z$alpha_true)^2)))
}))
recovery_overall <- function(raw)do.call(rbind,lapply(split(raw,raw$scenario),function(z){
 per<-aggregate(I((alpha_hat-alpha_true)^2)~rep,z,mean)
 data.frame(scenario=z$scenario[1],replications=length(unique(z$rep)),periods=length(unique(z$t)),
            bias=mean(z$alpha_hat-z$alpha_true),rmse=sqrt(mean((z$alpha_hat-z$alpha_true)^2)),
            rmse_mc_se=sd(per[[2]])/sqrt(nrow(per))/(2*sqrt(mean(per[[2]]))),
            tie_rate=mean(z$n_ties>1),mean_ties=mean(z$n_ties),true_in_argmin=mean(z$true_in_argmin),
            mean_distinct_reserves=mean(z$distinct_reserves),oracle_tie_rate=mean(z$oracle_n_ties>1),
            oracle_bias=mean(z$oracle_alpha_min-z$alpha_true),
            oracle_rmse=sqrt(mean((z$oracle_alpha_min-z$alpha_true)^2)))
}))
run_parallel <- function(ids,fun,cores) {
 zz<-parallel::mclapply(ids,fun,mc.cores=cores,mc.set.seed=FALSE)
 if(any(vapply(zz,inherits,logical(1),'try-error')))stop('A simulation worker failed')
 do.call(rbind,zz)
}
write_sim <- function(z,name,cfg)write.csv(z,file.path(cfg$output,name),row.names=FALSE)
plot_simulations <- function(a,exact,b,c,d,cfg,only_a=FALSE) {
 out<-cfg$output;cols<-c('#007C83','#B65B27')
 pdf(file.path(out,'panelA_error_control.pdf'),width=11,height=6.7,useDingbats=FALSE)
 par(mfrow=c(2,3),mar=c(4,4,2.7,.8))
 for(n in cfg$ns){
  plot(NA,xlim=c(0,1),ylim=c(0,.65),xlab=expression(alpha),ylab='Population violation frequency',main=paste('SE, n =',n));abline(h=cfg$delta,lty=3,col='grey40')
  for(i in 1:2){z<-a[a$n==n & a$method==c('NP umbrella','plug-in quantile')[i],];z<-z[order(z$alpha),];lines(z$alpha,z$violation,type='b',pch=c(16,1)[i],col=cols[i],cex=.7);arrows(z$alpha,z$violation_low,z$alpha,z$violation_high,angle=90,code=3,length=.025,col=cols[i])}
  if(n==min(cfg$ns))legend('topleft',c('NP umbrella','Plug-in quantile','Nominal delta = 0.1'),col=c(cols,'grey40'),lty=c(1,1,3),pch=c(16,1,NA),bty='n',cex=.7)
 }
 for(n in cfg$ns){z<-a[a$n==n & a$method=='NP umbrella',];z<-z[order(z$alpha),];plot(z$alpha,z$R1,type='b',pch=16,col=cols[1],ylim=c(0,1),xlab=expression(alpha),ylab='Mean population error',main=paste('NP operating curve, n =',n));lines(z$alpha,z$R0,type='b',pch=1,col=cols[2]);abline(0,1,lty=3,col='grey50');if(n==min(cfg$ns))legend('topright',c('Missed crisis','False alarm'),col=cols,pch=c(16,1),lty=1,bty='n',cex=.8)}
 dev.off()
 if(only_a)return(invisible(NULL))
 pdf(file.path(out,'panelB_welfare.pdf'),width=10,height=4.5,useDingbats=FALSE)
 par(mfrow=c(1,2),mar=c(4.3,4.4,2.4,.8));plot(NA,xlim=range(exact$p),ylim=range(c(0,exact$ratio),finite=TRUE),xlab='True crisis probability',ylab='Loss ratio: downward / upward error',main='Constrained reserve decision');abline(h=1,lty=3,col='grey50')
 for(i in 1:3){z<-exact[exact$relative_error==c(.25,.5,.75)[i],];lines(z$p,z$ratio,lty=i,col='#007C83')};legend('topright',c('25% error','50% error','75% error'),lty=1:3,col='#007C83',bty='n',cex=.8)
 plot(NA,xlim=c(0,1),ylim=range(b$expected_loss),log='y',xlab=expression(alpha),ylab='Mean utility loss (log scale)',main='Fixed SE with 53 evaluation countries')
 for(i in 1:2){z<-b[b$p==c(cfg$pi0,.2)[i],];z<-z[order(z$alpha),];lines(z$alpha,z$expected_loss,col=cols[i],type='b',pch=c(16,1)[i]);zi<-z[z$loss_se>1e-14,];arrows(zi$alpha,pmax(min(b$expected_loss)/10,zi$expected_loss-1.96*zi$loss_se),zi$alpha,zi$expected_loss+1.96*zi$loss_se,angle=90,code=3,length=.025,col=cols[i])}
 legend('topright',sprintf('True risk = %.3f',c(cfg$pi0,.2)),col=cols,pch=c(16,1),lty=1,bty='n',cex=.8);dev.off()
 pdf(file.path(out,'panelC_alpha_recovery.pdf'),width=10,height=4.5,useDingbats=FALSE)
 par(mfrow=c(1,2),mar=c(4,4.3,2.3,.8));z<-c[c$scenario=='baseline',];z<-z[order(z$t),]
 plot(z$t,z$alpha_true,type='n',ylim=c(0,1),xlab='Period',ylab=expression(alpha),main='Recovery with independent SE fits')
 polygon(c(z$t,rev(z$t)),c(z$alpha_hat_q025,rev(z$alpha_hat_q975)),col='#007C8325',border=NA)
 lines(z$t,z$alpha_true,lwd=2);lines(z$t,z$alpha_hat_mean,col=cols[1],lwd=2,lty=2)
 legend('bottomleft',c('True path','Recovered mean','95% replication interval'),col=c('black',cols[1],'#007C8340'),lty=c(1,2,1),lwd=c(2,2,8),bty='n',cex=.75)
 plot(z$t,z$tie_rate,type='b',ylim=c(0,1),pch=16,col=cols[1],xlab='Period',ylab='Fraction of replications',main='Recovery and identifiability');lines(z$t,z$true_in_argmin,type='b',pch=1,col=cols[2]);lines(z$t,z$oracle_tie_rate,lty=3)
 legend('topleft',c('Analyst argmin has ties','True alpha in analyst argmin','Oracle argmin has ties'),col=c(cols,'black'),lty=c(1,1,3),pch=c(16,1,NA),bty='n',cex=.72);dev.off()
 pdf(file.path(out,'panelD_recovery_robustness.pdf'),width=7,height=4.5,useDingbats=FALSE)
 par(mar=c(7,4,2,1));bp<-barplot(d$rmse,names.arg=gsub('_',' ',d$scenario),las=2,col='#007C83',ylab='Recovery RMSE',ylim=c(0,max(d$rmse+2*d$rmse_mc_se)*1.1));arrows(bp,d$rmse-1.96*d$rmse_mc_se,bp,d$rmse+1.96*d$rmse_mc_se,angle=90,code=3,length=.04);dev.off()
}
main_simulations <- function(mode='full') {
 cfg<-simulation_config(mode);dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
 started<-Sys.time();prefix<-if(mode=='panelA_repair')'panelA_repair_'else ''
 capture.output(sessionInfo(),file=file.path(cfg$output,paste0(prefix,'sessionInfo.txt')))
 saveRDS(list(config=cfg,pipeline=PIPELINE_VERSION,source_md5=tools::md5sum(c('corrected/core.R','corrected/simulations.R'))),file.path(cfg$output,paste0(prefix,'run_metadata.rds')))
 all_a<-list()
 for(n in cfg$ns){cat('Panel A n=',n,' started ',format(Sys.time()),'\n',sep='');flush.console()
  zz<-run_parallel(seq_len(cfg$B_A),function(b)panel_a_rep(b,n,cfg$alphas,cfg),cfg$cores)
  write_sim(zz,paste0('panelA_n',n,'_raw.csv'),cfg);all_a[[as.character(n)]]<-zz
  print(summarize_a(zz)[,c('n','method','alpha','violation','fallback')]);flush.console()
 }
 a<-summarize_a(do.call(rbind,all_a));write_sim(a,'panelA_summary.csv',cfg)
 if(mode=='panelA_repair') {
  plot_simulations(a,NULL,NULL,NULL,NULL,cfg,only_a=TRUE)
  elapsed<-as.numeric(difftime(Sys.time(),started,units='secs'))
  writeLines(sprintf('Corrected Panel A: 3000 replications, elapsed %.6f seconds.',elapsed),file.path(cfg$output,'panelA_repair_run_summary.txt'))
  cat('Corrected Panel A completed in ',round(elapsed,1),' seconds.\n',sep='')
  return(invisible(a))
 }
 exact<-panel_b_exact(cfg);write_sim(exact,'panelB1_constrained_asymmetry.csv',cfg)
 cat('Panel B pipeline started ',format(Sys.time()),'\n');flush.console()
 br<-run_parallel(seq_len(cfg$B_B),function(b)panel_b_rep(b,cfg),cfg$cores)
 write_sim(br,'panelB2_raw.csv',cfg);b<-summarize_b(br);write_sim(b,'panelB2_expected_loss.csv',cfg)
 cat('Panels C and D started ',format(Sys.time()),'\n');flush.console()
 cr<-run_parallel(seq_len(cfg$B_C),function(b)recovery_rep(b,cfg,scenarios=if(b<=cfg$B_D)c('baseline','r_known','r_ignored','calibration_jitter')else 'baseline'),cfg$cores)
 write_sim(cr,'panelCD_raw.csv',cfg);cs<-summarize_recovery(cr);write_sim(cs,'panelCD_by_period.csv',cfg)
 co<-recovery_overall(cr[cr$scenario=='baseline',]);write_sim(co,'panelC_overall.csv',cfg)
 # All robustness scenarios use the same first B_D replications and innovations.
 ds<-recovery_overall(cr[cr$rep<=cfg$B_D,]);write_sim(ds,'panelD_overall.csv',cfg)
 plot_simulations(a,exact,b,cs,ds,cfg)
 elapsed<-as.numeric(difftime(Sys.time(),started,units='secs'))
 capture.output(list(elapsed_seconds=elapsed,panelC=co,panelD=ds),file=file.path(cfg$output,'run_summary.txt'))
 cat('Completed in ',round(elapsed,1),' seconds.\n',sep='');print(co);print(ds)
 invisible(list(A=a,B1=exact,B2=b,C=cs,D=ds))
}
if(sys.nframe()==0L) {
 args<-commandArgs(trailingOnly=TRUE);main_simulations(if(length(args))args[1]else 'full')
}
