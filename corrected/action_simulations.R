# Direct reserve-action simulation with exact Gaussian population risks.
# Run from code/: Rscript --vanilla corrected/action_simulations.R
# No score is estimated. All randomness is independent calibration/validation.
source('corrected/core.R')
source('corrected/actions.R')

gaussian_action_rates<-function(cut,mu=1.5) {
  stopifnot(all(!is.na(cut)),length(mu)==1L,is.finite(mu),mu>0)
  data.frame(R0=pnorm(cut,lower.tail=FALSE),R1=pnorm(cut-mu))
}
gaussian_action_value<-function(cut,p,th,low=.05,high=.15,mu=1.5) {
  stopifnot(length(p)==1L,is.finite(p),p>0,p<1)
  e<-gaussian_action_rates(cut,mu);z<-action_components(th,low,high)
  base<-(1-p)*z$u_n_low+p*z$u_s_low
  U<-(1-p)*((1-e$R0)*z$u_n_low+e$R0*z$u_n_high)+
    p*(e$R1*z$u_s_low+(1-e$R1)*z$u_s_high)
  share<-(1-p)*e$R0+p*(1-e$R1)
  data.frame(e,utility=U,utility_base=base,gain=U-base,
    identity_gain=p*z$b*(1-e$R1)-(1-p)*z$a*e$R0,
    CE=vapply(U,ce_gain,0,U_base=base,sigma=th$sigma),
    mean_reserves=low+(high-low)*share,
    normal_excess_reserves=(high-low)*(1-p)*e$R0,
    missed_crisis_shortfall=(high-low)*p*e$R1,
    normal_utility_cost=(1-p)*z$a*e$R0,
    crisis_utility_benefit=p*z$b*(1-e$R1))
}
gaussian_welfare_cutoff<-function(p,components,mu=1.5) {
  stopifnot(p>0,p<1,components$a>0,components$b>0,mu>0)
  mu/2+log((1-p)*components$a/(p*components$b))/mu
}
empirical_action_cutoff<-function(cal,alpha) {
  stopifnot(length(cal)>0L,all(is.finite(cal)),alpha>=0,alpha<=1)
  if(alpha==1)return(Inf)
  # At a tied score, every tied observation receives high reserves. Selecting
  # the largest feasible observed value also gives the largest real cutoff.
  s<-sort(unique(cal));misses<-vapply(s,function(c)sum(cal<c),0)
  max(s[misses<=floor(alpha*length(cal)+1e-12)])
}
action_binomial_upper<-function(k,n,eta=.05) {
  stopifnot(n>0,all(k>=0 & k<=n),all(k==floor(k)),eta>0,eta<1)
  ans<-rep(1,length(k));i<-k<n
  ans[i]<-qbeta(1-eta,k[i]+1,n-k[i]);ans
}
action_welfare_lower<-function(p,alpha,U0,components) {
  p*components$b*(1-alpha)-(1-p)*components$a*U0
}
action_sim_config<-function(B=2000L)list(seed=20261015L,B=B,
  mu=1.5,prevalences=c(.05,.20),m=c(25L,100L),alphas=c(.1,.2,.25),
  delta=.1,n0=500L,eta=.05,low=.05,high=.15,
  output='outputs_action/simulations',theta_history_end=2007L)

action_simulate<-function(cfg,th) {
  rows<-list();i<-0L;co<-action_components(th,cfg$low,cfg$high)
  methods<-c('NP calibrated','Empirical calibration','NP with welfare cap',
    'Population NP','Constrained welfare oracle','Welfare oracle','Always low','Always high')
  for(m in cfg$m) {
    # Independent seeds by m preserve reproducibility if settings are extended.
    set.seed(cfg$seed+1000L*m)
    crisis<-matrix(rnorm(cfg$B*m,cfg$mu),cfg$B,m)
    normal<-matrix(rnorm(cfg$B*cfg$n0),cfg$B,cfg$n0)
    sorted<-t(apply(crisis,1,sort))
    for(alpha in cfg$alphas) {
      rule<-np_rule(seq_len(m),alpha,cfg$delta)
      np_cut<-if(rule$mode=='all_warning')rep(-Inf,cfg$B)else sorted[,rule$k]
      # Continuous Gaussian draws have no ties; the generic comparator above
      # is tested separately on tied samples.
      emp_k<-floor(alpha*m+1e-12)+1L
      emp_cut<-sorted[,emp_k]
      c_alpha<-cfg$mu+qnorm(alpha)
      for(p in cfg$prevalences) {
        cw<-gaussian_welfare_cutoff(p,co,cfg$mu)
        cuts<-list(np_cut,emp_cut,pmin(np_cut,cw),rep(c_alpha,cfg$B),
          rep(min(c_alpha,cw),cfg$B),rep(cw,cfg$B),rep(Inf,cfg$B),rep(-Inf,cfg$B))
        for(j in seq_along(methods)) {
          cut<-cuts[[j]];v<-gaussian_action_value(cut,p,th,cfg$low,cfg$high,cfg$mu)
          k0<-rowSums(normal>=cut)
          upper<-action_binomial_upper(k0,cfg$n0,cfg$eta)
          lb<-action_welfare_lower(p,alpha,upper,co)
          eligible<-j%in%c(1L,3L)
          i<-i+1L
          rows[[i]]<-data.frame(rep=seq_len(cfg$B),m=m,p=p,alpha=alpha,
            method=methods[j],threshold=cut,v,calibration_k=if(j==1L)rule$k else if(j==2L)emp_k else NA_integer_,
            normal_count=k0,R0_upper=upper,R1_violation=v$R1>alpha+1e-12,
            R0_upper_failure=v$R0>upper+1e-12,
            certificate_lower=if(eligible)lb else NA_real_,
            certificate_covered=if(eligible)v$gain>=lb-1e-12 else NA,
            certificate_positive=if(eligible)lb>0 else NA,
            certificate_false_positive=if(eligible)(lb>0 & v$gain<=0)else NA,
            joint_failure=if(eligible)(v$R1>alpha+1e-12 | v$R0>upper+1e-12)else NA,
            true_positive_gain=v$gain>0,
            population_constraint_binding=pnorm(cw-cfg$mu)>alpha+1e-12)
        }
      }
    }
  }
  do.call(rbind,rows)
}
action_sim_summary<-function(raw) {
  metrics<-c('R0','R1','utility','gain','CE','mean_reserves','normal_excess_reserves',
    'missed_crisis_shortfall','normal_utility_cost','crisis_utility_benefit',
    'certificate_lower','R1_violation','R0_upper_failure','certificate_covered',
    'certificate_positive','certificate_false_positive','joint_failure','true_positive_gain')
  binary<-c('R1_violation','R0_upper_failure','certificate_covered','certificate_positive',
    'certificate_false_positive','joint_failure','true_positive_gain')
  groups<-split(raw,list(raw$m,raw$p,raw$alpha,raw$method),drop=TRUE)
  do.call(rbind,lapply(groups,function(z){
    ans<-z[1,c('m','p','alpha','method','population_constraint_binding')];ans$B<-nrow(z)
    for(nm in metrics) {
      v<-z[[nm]];v<-v[is.finite(v)]
      meanv<-if(length(v))mean(v)else NA_real_
      se<-if(length(v)>1L)sd(v)/sqrt(length(v))else NA_real_
      ci<-if(!length(v))c(NA_real_,NA_real_)else if(nm%in%binary)
        binom.test(sum(v),length(v))$conf.int else meanv+c(-1,1)*qnorm(.975)*se
      ans[[paste0(nm,'_mean')]]<-meanv;ans[[paste0(nm,'_mcse')]]<-se
      ans[[paste0(nm,'_low95')]]<-ci[1];ans[[paste0(nm,'_high95')]]<-ci[2]
    }
    ans
  }))
}
action_sim_theory<-function(cfg,th) {
  co<-action_components(th,cfg$low,cfg$high)
  z<-expand.grid(m=cfg$m,p=cfg$prevalences,alpha=cfg$alphas)
  do.call(rbind,lapply(seq_len(nrow(z)),function(i){
    r<-z[i,];np<-np_rule(seq_len(r$m),r$alpha,cfg$delta)
    ek<-floor(r$alpha*r$m+1e-12)+1L
    cw<-gaussian_welfare_cutoff(r$p,co,cfg$mu)
    data.frame(r,NP_k=np$k,empirical_k=ek,
      NP_exact_violation=if(np$mode=='all_warning')0 else pbinom(np$k-1,r$m,r$alpha),
      empirical_exact_violation=pbinom(ek-1,r$m,r$alpha),
      welfare_cutoff=cw,welfare_R1=pnorm(cw-cfg$mu),
      NP_population_cutoff=cfg$mu+qnorm(r$alpha),
      constrained_welfare_cutoff=min(cw,cfg$mu+qnorm(r$alpha)),
      binding=pnorm(cw-cfg$mu)>r$alpha,
      joint_failure_upper=cfg$delta+cfg$eta)
  }))
}
action_sim_plot<-function(summary,theory,cfg,path) {
  pdf(path,width=10.6,height=4.5,family='Helvetica',pointsize=11)
  old<-par(mfrow=c(1,2),mar=c(4.4,4.3,2.3,1),mgp=c(2.5,.8,0))
  on.exit({par(old);dev.off()})
  cols<-c('#006D77','#C55A11','#6B4C9A','#222222')
  plot(NA,xlim=c(.075,.275),ylim=c(0,.68),xlab='Missed-crisis bound alpha',
    ylab='Probability of exceeding alpha',main='A. Population error control',xaxt='n')
  axis(1,at=cfg$alphas);abline(h=cfg$delta,col='gray45',lty=3)
  k<-0
  for(method in c('NP calibrated','Empirical calibration'))for(m in cfg$m) {
    k<-k+1;z<-summary[summary$method==method & summary$m==m & summary$p==cfg$prevalences[1],]
    z<-z[order(z$alpha),];x<-z$alpha+(k-2.5)*.002
    arrows(x,z$R1_violation_low95,x,z$R1_violation_high95,angle=90,code=3,length=.035,col=cols[k])
    lines(x,z$R1_violation_mean,type='b',pch=c(16,17,1,2)[k],col=cols[k],lty=if(k<=2)1 else 2)
  }
  legend('topleft',c('NP: m=25','NP: m=100','Empirical: m=25','Empirical: m=100'),
    col=cols,pch=c(16,17,1,2),lty=c(1,1,2,2),bty='n',cex=.8)
  text(.27,.112,'delta=0.10',adj=1,cex=.75,col='gray35')
  zz<-summary[summary$p==.05,]
  yr<-range(zz$CE_low95,zz$CE_high95,na.rm=TRUE)*10000
  plot(NA,xlim=c(.085,.265),ylim=yr,xlab='Missed-crisis bound alpha',
    ylab='Consumption-equivalent gain (basis points)',main='B. Policy value at prevalence 0.05',xaxt='n')
  axis(1,at=cfg$alphas);abline(h=0,col='gray80')
  meth<-c('NP calibrated','NP calibrated','Empirical calibration','NP with welfare cap','Constrained welfare oracle','Welfare oracle')
  ms<-c(25,100,100,100,100,100)
  cl<-c(cols[1],cols[2],cols[3],'#468B45','#222222','#666666')
  for(j in seq_along(meth)) {
    z<-zz[zz$method==meth[j] & zz$m==ms[j],];z<-z[order(z$alpha),]
    lines(z$alpha,z$CE_mean*10000,type='b',pch=c(16,17,1,15,2,NA)[j],
      col=cl[j],lty=c(1,1,2,1,2,3)[j])
    if(j<=4)arrows(z$alpha,z$CE_low95*10000,z$alpha,z$CE_high95*10000,
      angle=90,code=3,length=.035,col=cl[j])
  }
  legend('bottomright',c('NP: m=25','NP: m=100','Empirical: m=100',
    'NP + known welfare cap: m=100','Constrained welfare oracle','Welfare oracle'),
    col=cl,pch=c(16,17,1,15,2,NA),lty=c(1,1,2,1,2,3),bty='n',cex=.75)
}
run_action_simulations<-function(B=2000L) {
  source('corrected/data.R')
  cfg<-action_sim_config(B);th<-calibration_theta(load_inputs(),2008)
  dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
  start<-proc.time()[3];raw<-action_simulate(cfg,th)
  summary<-action_sim_summary(raw);theory<-action_sim_theory(cfg,th)
  stopifnot(max(abs(raw$gain-raw$identity_gain))<1e-12)
  saveRDS(list(config=cfg,theta=th,raw=raw),file.path(cfg$output,'action_simulation_replicates.rds'),compress='xz')
  raw_csv<-gzfile(file.path(cfg$output,'action_simulation_replicates.csv.gz'),'wt')
  write.csv(raw,raw_csv,row.names=FALSE);close(raw_csv)
  write.csv(summary,file.path(cfg$output,'action_simulation_summary.csv'),row.names=FALSE)
  write.csv(theory,file.path(cfg$output,'action_simulation_exact_benchmarks.csv'),row.names=FALSE)
  write.csv(as.data.frame(c(th,action_components(th,cfg$low,cfg$high))),file.path(cfg$output,'action_simulation_parameters.csv'),row.names=FALSE)
  select<-c('m','p','alpha','method','B','R0_mean','R1_mean','R1_violation_mean',
    'R1_violation_low95','R1_violation_high95','utility_mean','CE_mean','CE_mcse',
    'mean_reserves_mean','normal_excess_reserves_mean','missed_crisis_shortfall_mean')
  write.csv(summary[,select],file.path(cfg$output,'action_simulation_policy_table.csv'),row.names=FALSE)
  select<-c('m','p','alpha','method','B','certificate_lower_mean','certificate_covered_mean',
    'certificate_covered_low95','certificate_covered_high95','certificate_positive_mean',
    'certificate_positive_low95','certificate_positive_high95','joint_failure_mean',
    'joint_failure_low95','joint_failure_high95','certificate_false_positive_mean')
  write.csv(summary[summary$method%in%c('NP calibrated','NP with welfare cap'),select],
    file.path(cfg$output,'action_simulation_certificate_table.csv'),row.names=FALSE)
  action_sim_plot(summary,theory,cfg,file.path(cfg$output,'action_simulation_two_panels.pdf'))
  capture.output(sessionInfo(),file=file.path(cfg$output,'sessionInfo.txt'))
  elapsed<-proc.time()[3]-start
  writeLines(c(sprintf('Replicates per setting: %d',B),sprintf('Rows: %d',nrow(raw)),
    sprintf('Elapsed seconds: %.3f',elapsed),sprintf('Maximum identity error: %.17g',max(abs(raw$gain-raw$identity_gain)))),
    file.path(cfg$output,'run_validation.txt'))
  cat('Completed',B,'replicates per setting in',round(elapsed,2),'seconds.\n')
  invisible(list(config=cfg,theta=th,summary=summary,theory=theory))
}
if(sys.nframe()==0L)run_action_simulations()
