# Reproducible retrospective evaluation of directly implemented reserve actions.
source('corrected/data.R')
source('corrected/actions.R')
ACTION_VERSION<-'2026-10-15-direct-reserve-actions-v3.2'
choose_policy_threshold<-function(scores,y,weights,a,b,criterion=c('economic','equal_error')){
 criterion<-match.arg(criterion)
 stopifnot(length(scores)==length(y),length(y)==length(weights),all(y%in%0:1),all(weights>0),sum(y==0)>0,sum(y==1)>0)
 cand<-c(-Inf,sort(unique(scores)),Inf)
 p<-sum(weights*y)/sum(weights)
 r0<-vapply(cand,function(t)weighted.mean(scores[y==0]>=t,weights[y==0]),0)
 r1<-vapply(cand,function(t)weighted.mean(scores[y==1]<t,weights[y==1]),0)
 loss<-if(criterion=='economic')(1-p)*a*r0+p*b*r1 else r0+r1
 ii<-which(loss<=min(loss)+1e-12);i<-tail(ii,1)
 list(threshold=cand[i],loss=loss[i],R0=r0[i],R1=r1[i],n_ties=length(ii),criterion=criterion)
}
fit_policy_models<-function(tr,features,seed=20261015L){
 stopifnot(all(!is.na(tr$SuddenStop_GI)))
 y<-as.integer(tr$SuddenStop_GI);x<-as.matrix(tr[,features,drop=FALSE])
 se<-fit_np(x,y,'SE',seed=seed,d=Inf)
 # Reuse exactly the SE fitting/calibration split and fitting-only preprocessing.
 # Warm starts improve numerical convergence without changing target lambda=.01.
 logfit<-se;logfit$model<-'logistic'
 xf<-x[se$fit_idx,se$selected,drop=FALSE]
 for(j in seq_len(ncol(xf)))xf[!is.finite(xf[,j]),j]<-se$median[colnames(xf)[j]]
 ridge<-glmnet::glmnet(xf,y[se$fit_idx],family='binomial',alpha=0,lambda=c(.1,.03,.01),
        standardize=TRUE,thresh=1e-8,maxit=1000000)
 stopifnot(ridge$jerr==0,min(ridge$lambda)==.01)
 logfit$base<-list(fit=ridge,lambda=.01,converged=TRUE)
 logfit$cal_scores<-predict_score(logfit,x[se$cal_idx,,drop=FALSE])
 fits<-list(SE=se,logistic=logfit)
 stopifnot(identical(fits$SE$fit_idx,fits$logistic$fit_idx),identical(fits$SE$cal_idx,fits$logistic$cal_idx))
 idx<-fits$SE$fit_idx;ratio<-sum(y==1)/sum(y[idx]==1)
 scores<-lapply(fits,function(f)predict_score(f,x[idx,,drop=FALSE]))
 list(fits=fits,fit_scores=scores,fit_y=y[idx],weights=ifelse(y[idx]==1,ratio,1),
      p_train=mean(y),n_train=length(y),n_crisis=sum(y),odds_multiplier=ratio,seed=seed,features=features)
}
predict_policy_actions<-function(bundle,te,th,low=.05,high=.15,alphas=c(.10,.20,.25),delta=.10,prob_eps=.001){
 # Register glmnet prediction methods when reusing a saved model in a fresh session.
 stopifnot(requireNamespace('glmnet',quietly=TRUE))
 # Only named predictor columns are accessed. Evaluation labels cannot enter decisions.
 x<-as.matrix(te[,bundle$features,drop=FALSE]);n<-nrow(x)
 score<-lapply(bundle$fits,predict_score,x=x)
 qraw<-plogis(score$logistic+log(bundle$odds_multiplier));q<-clip_probability(qraw,prob_eps)
 comp<-action_components(th,low,high);ans<-list();i<-0L
 add<-function(method,rho,flag=rep(NA,n),threshold=NA_real_,score_value=rep(NA_real_,n),alpha=NA_real_,mode='none'){
  i<<-i+1L;ans[[i]]<<-data.frame(row=seq_len(n),method=method,rho=rep(rho,length.out=n),flag=rep(flag,length.out=n),
   threshold=threshold,score=score_value,alpha=alpha,mode=mode,q_ridge_raw=qraw,q_ridge=q,p_train=bundle$p_train,low=low,high=high)
 }
 add('always_low',low,FALSE);add('always_high',high,TRUE)
 add('historical_constant',reserve_policy(clip_probability(bundle$p_train),th))
 add('probability_unrestricted',reserve_policy(q,th),score_value=score$logistic)
 add('probability_bounded',pmin(high,pmax(low,reserve_policy(q,th))),score_value=score$logistic)
 for(m in names(bundle$fits)){
  f<-bundle$fits[[m]]
  for(alpha in alphas){r<-np_rule(f$cal_scores,alpha,delta);flag<-apply_rule(r,score[[m]])
   add(paste0(m,'_NP_',sprintf('%.2f',alpha)),ifelse(flag,high,low),flag,r$threshold,score[[m]],alpha,r$mode)}
  for(criterion in c('economic','equal_error')){
   r<-choose_policy_threshold(bundle$fit_scores[[m]],bundle$fit_y,bundle$weights,comp$a,comp$b,criterion)
   flag<-score[[m]]>=r$threshold
   add(paste0(m,'_',criterion),ifelse(flag,high,low),flag,r$threshold,score[[m]],mode='fitting_only')
  }
 }
 do.call(rbind,ans)
}
evaluate_policy_rows<-function(pred,te,th,seed,analysis){
 z<-cbind(te[pred$row,c('Country_Code','Country_Name','Year','SuddenStop_GI')],pred)
 z$seed<-seed;z$analysis<-analysis;z$utility<-NA_real_;labeled<-!is.na(z$SuddenStop_GI)
 z$utility[labeled]<-realized_utility(z$SuddenStop_GI[labeled],z$rho[labeled],th)
 z$u_low<-NA_real_;z$u_low[labeled]<-realized_utility(z$SuddenStop_GI[labeled],z$low[labeled],th)
 z$u_high<-NA_real_;z$u_high[labeled]<-realized_utility(z$SuddenStop_GI[labeled],z$high[labeled],th)
 z$u_normal<-realized_utility(rep(0,nrow(z)),z$rho,th)
 z$u_crisis<-realized_utility(rep(1,nrow(z)),z$rho,th)
 z$u_normal_low<-realized_utility(rep(0,nrow(z)),z$low,th)
 z$u_crisis_low<-realized_utility(rep(1,nrow(z)),z$low,th)
 z$expected_ridge<-expected_action_utility(z$q_ridge,z$rho,th)
 z$expected_prevalence<-expected_action_utility(z$p_train,z$rho,th)
 z$sigma<-th$sigma
 z
}
summarize_policy_rows<-function(z){
 lab<-z[!is.na(z$utility),];y<-lab$SuddenStop_GI;U<-mean(lab$utility);Ub<-mean(lab$u_low)
 binary<-all(!is.na(lab$flag));normal<-y==0;crisis<-y==1
 data.frame(method=z$method[1],low=z$low[1],high=z$high[1],seed=z$seed[1],analysis=z$analysis[1],
  n_forecast=nrow(z),n_labeled=nrow(lab),n_crisis=sum(crisis),n_normal=sum(normal),
  mean_reserve=mean(lab$rho),mean_reserve_all_forecasts=mean(z$rho),max_reserve=max(z$rho),
  utility=U,utility_low=Ub,ce_gain=ce_gain(U,Ub,z$sigma[1]),ce_bp=10000*ce_gain(U,Ub,z$sigma[1]),
  expected_ridge_utility=mean(lab$expected_ridge),prevalence_standardized_utility=mean(lab$expected_prevalence),
  expected_ridge_all_forecasts=mean(z$expected_ridge),prevalence_standardized_all_forecasts=mean(z$expected_prevalence),
  warnings=if(binary)sum(lab$flag)else NA_integer_,detections=if(binary)sum(lab$flag & crisis)else NA_integer_,
  false_alarms=if(binary)sum(lab$flag & normal)else NA_integer_,misses=if(binary)sum(!lab$flag & crisis)else NA_integer_,
  R0=if(binary && any(normal))mean(lab$flag[normal])else NA_real_,R1=if(binary && any(crisis))mean(!lab$flag[crisis])else NA_real_,
  normal_cost=if(any(normal))mean(lab$u_normal_low[normal]-lab$u_normal[normal])else NA_real_,
  crisis_benefit=if(any(crisis))mean(lab$u_crisis[crisis]-lab$u_crisis_low[crisis])else NA_real_,
  normal_cost_contribution=sum(lab$u_normal_low[normal]-lab$u_normal[normal])/nrow(lab),
  crisis_benefit_contribution=sum(lab$u_crisis[crisis]-lab$u_crisis_low[crisis])/nrow(lab),
  q_clipped_low=sum(z$q_ridge_raw<.001),q_clipped_high=sum(z$q_ridge_raw>.999))
}
summarize_grouped<-function(z){
 key<-interaction(z$analysis,z$seed,z$low,z$high,z$method,drop=TRUE)
 do.call(rbind,lapply(split(z,key),summarize_policy_rows))
}
split_membership<-function(bundle,tr,analysis,year=NA_integer_){
 data.frame(analysis=analysis,seed=bundle$seed,forecast_year=year,Country_Code=tr$Country_Code,Year=tr$Year,
  SuddenStop_GI=tr$SuddenStop_GI,role=ifelse(seq_len(nrow(tr))%in%bundle$fits$SE$fit_idx,'fit','crisis_calibration'))
}
bootstrap_policy<-function(rows,B=1000L,seed=20261016L){
 # Paired resampling conditional on fitted models; no refitting or causal interpretation.
 sigma<-unique(rows$sigma)
 if(length(sigma)!=1L || !is.finite(sigma) || sigma<=0)
  stop('Bootstrap CE requires one common positive risk-aversion parameter.')
 z<-rows[!is.na(rows$utility),];methods<-unique(z$method);ct<-sort(unique(z$Country_Code));years<-sort(unique(z$Year))
 units<-unique(z[,c('Country_Code','Year')]);units<-units[order(units$Country_Code,units$Year),];key<-paste(units$Country_Code,units$Year)
 U<-sapply(methods,function(m){d<-z[z$method==m,];d$utility[match(key,paste(d$Country_Code,d$Year))]})
 BL<-z[z$method==methods[1],];base<-BL$u_low[match(key,paste(BL$Country_Code,BL$Year))]
 set.seed(seed);ans<-vector('list',2L*B);ii<-0L
 for(kind in c('country','two_year_moving_block'))for(b in seq_len(B)){
  if(kind=='country'){
   draw<-sample(ct,length(ct),replace=TRUE);w<-tabulate(match(draw,ct),nbins=length(ct))[match(units$Country_Code,ct)]
  }else{
   starts<-sample(seq_len(length(years)-1L),ceiling(length(years)/2),replace=TRUE)
   draw<-unlist(lapply(starts,function(j)years[j:(j+1L)]))[seq_len(length(years))]
   w<-tabulate(match(draw,years),nbins=length(years))[match(units$Year,years)]
  }
  ub<-weighted.mean(base,w);u<-colSums(U*w)/sum(w)
  ii<-ii+1L;ans[[ii]]<-data.frame(resampling=kind,rep=b,method=methods,utility=u,sigma=sigma,
      ce_bp=vapply(u,function(v)10000*ce_gain(v,ub,sigma),0))
 }
 raw<-do.call(rbind,ans)
 summary<-do.call(rbind,lapply(split(raw,interaction(raw$resampling,raw$method,drop=TRUE)),function(d){
  data.frame(resampling=d$resampling[1],method=d$method[1],B=nrow(d),ce_bp_mean=mean(d$ce_bp),
   ce_bp_lo=unname(quantile(d$ce_bp,.025)),ce_bp_hi=unname(quantile(d$ce_bp,.975)),fraction_positive=mean(d$ce_bp>0))}))
 list(raw=raw,summary=summary)
}
paired_policy_comparisons<-function(boot_raw,point_summary){
 sigma<-unique(boot_raw$sigma)
 if(length(sigma)!=1L || !is.finite(sigma) || sigma<=0)
  stop('Paired CE requires one common positive risk-aversion parameter.')
 methods<-grep('_NP_',unique(boot_raw$method),value=TRUE)
 result<-list();i<-0L
 for(m in methods)for(comparator in c('always_high',sub('_NP_.*$','_economic',m))){
  a<-boot_raw[boot_raw$method==m,];b<-boot_raw[boot_raw$method==comparator,]
  idx<-match(paste(a$resampling,a$rep),paste(b$resampling,b$rep))
  stopifnot(!anyNA(idx))
  d<-data.frame(resampling=a$resampling,rep=a$rep,method=m,comparator=comparator,
   delta_utility=a$utility-b$utility[idx],
   relative_ce_bp=mapply(function(u,v)10000*ce_gain(u,v,sigma),a$utility,b$utility[idx]))
  i<-i+1L;result[[i]]<-d
 }
 raw<-do.call(rbind,result)
 sm<-do.call(rbind,lapply(split(raw,interaction(raw$resampling,raw$method,raw$comparator,drop=TRUE)),function(d){
  a<-point_summary$utility[point_summary$method==d$method[1]];b<-point_summary$utility[point_summary$method==d$comparator[1]]
  data.frame(resampling=d$resampling[1],method=d$method[1],comparator=d$comparator[1],B=nrow(d),
   point_delta_utility=a-b,point_relative_ce_bp=10000*ce_gain(a,b,sigma),
   delta_utility_lo=unname(quantile(d$delta_utility,.025)),delta_utility_hi=unname(quantile(d$delta_utility,.975)),
   relative_ce_bp_lo=unname(quantile(d$relative_ce_bp,.025)),relative_ce_bp_hi=unname(quantile(d$relative_ce_bp,.975)))
 }))
 list(raw=raw,summary=sm)
}
run_economic_sensitivity<-function(bundle,te,theta,out){
 scenarios<-list(baseline=theta,r_003=modifyList(theta,list(r=.03)),r_007=modifyList(theta,list(r=.07)),
                 sigma_1=modifyList(theta,list(sigma=1)),sigma_4=modifyList(theta,list(sigma=4)))
 rows<-lapply(names(scenarios),function(nm){
  th<-scenarios[[nm]];z<-evaluate_policy_rows(predict_policy_actions(bundle,te,th),te,th,bundle$seed,'economic_sensitivity')
  z$scenario<-nm;z$r<-th$r;z
 })
 summary<-do.call(rbind,lapply(rows,function(z){s<-summarize_grouped(z);s$scenario<-z$scenario[1];s$r<-z$r[1];s$sigma<-z$sigma[1];s}))
 write.csv(do.call(rbind,rows),file.path(out,'economic_sensitivity_predictions.csv'),row.names=FALSE)
 write.csv(summary,file.path(out,'economic_sensitivity_summary.csv'),row.names=FALSE)
 invisible(summary)
}
run_policy_evaluation<-function(){
 out<-'outputs_action/empirical';dir.create(out,recursive=TRUE,showWarnings=FALSE)
 inputs<-load_inputs();d<-inputs$data;tr<-d[d$Year<=2007 & !is.na(d$SuddenStop_GI),];te<-d[d$Year%in%2008:2017,]
 theta<-calibration_theta(inputs,2008);seed0<-20261015L
 tiers<-data.frame(low=c(.05,0,.10,.05),high=c(.15,.10,.20,.25))
 write.csv(data.frame(forecast_year=2008,as.data.frame(theta)),file.path(out,'frozen_theta.csv'),row.names=FALSE)
 comp<-do.call(rbind,lapply(seq_len(nrow(tiers)),function(j)cbind(tiers[j,],as.data.frame(action_components(theta,tiers$low[j],tiers$high[j])))))
 write.csv(comp,file.path(out,'tier_components.csv'),row.names=FALSE)
 all_primary<-list();all_summaries<-list();members<-list();solvers<-list();primary_bundle<-NULL
 for(b in 0:99){
  seed<-seed0+b;bundle<-fit_policy_models(tr,inputs$features,seed)
  if(b==0)primary_bundle<-bundle
  solvers[[b+1L]]<-data.frame(seed=seed,jerr=bundle$fits$logistic$base$fit$jerr,npasses=bundle$fits$logistic$base$fit$npasses,lambda=.01)
  pred<-predict_policy_actions(bundle,te,theta)
  rows<-evaluate_policy_rows(pred,te,theta,seed,'static')
  all_primary[[b+1L]]<-rows;all_summaries[[b+1L]]<-summarize_grouped(rows)
  members[[b+1L]]<-split_membership(bundle,tr,'static')
  if(b%%10==0){cat('Static split',b+1L,'of 100 completed\n');flush.console()}
 }
 write.csv(do.call(rbind,solvers),file.path(out,'ridge_solver_diagnostics.csv'),row.names=FALSE)
 rows<-do.call(rbind,all_primary);sm<-do.call(rbind,all_summaries);membership<-do.call(rbind,members)
 saveRDS(rows,file.path(out,'static_all_splits_predictions.rds'));write.csv(sm,file.path(out,'static_all_splits_summary.csv'),row.names=FALSE)
 write.csv(membership,file.path(out,'static_split_membership.csv'),row.names=FALSE)
 primary<-all_primary[[1]];pd<-primary[primary$method=='probability_unrestricted' & !is.na(primary$SuddenStop_GI),]
 write.csv(data.frame(n=nrow(pd),prevalence=mean(pd$SuddenStop_GI),mean_q_raw=mean(pd$q_ridge_raw),mean_q=mean(pd$q_ridge),brier=mean((pd$q_ridge-pd$SuddenStop_GI)^2),auc=score_auc(pd$q_ridge_raw,pd$SuddenStop_GI),odds_multiplier=primary_bundle$odds_multiplier),file.path(out,'ridge_probability_diagnostics.csv'),row.names=FALSE)
 write.csv(primary,file.path(out,'primary_predictions.csv'),row.names=FALSE)
 write.csv(all_summaries[[1]],file.path(out,'primary_summary.csv'),row.names=FALSE)
 saveRDS(primary_bundle,file.path(out,'primary_fit_bundle.rds'))
 tier_rows<-lapply(seq_len(nrow(tiers)),function(j)evaluate_policy_rows(predict_policy_actions(primary_bundle,te,theta,tiers$low[j],tiers$high[j]),te,theta,seed0,'tier'))
 write.csv(do.call(rbind,tier_rows),file.path(out,'tier_predictions.csv'),row.names=FALSE)
 write.csv(do.call(rbind,lapply(tier_rows,summarize_grouped)),file.path(out,'tier_summary.csv'),row.names=FALSE)
 clipping<-lapply(c(1e-6,.001,.01),function(e){p<-predict_policy_actions(primary_bundle,te,theta,prob_eps=e);p<-p[p$method%in%c('probability_unrestricted','probability_bounded'),];z<-summarize_grouped(evaluate_policy_rows(p,te,theta,seed0,'clipping'));z$epsilon<-e;z})
 write.csv(do.call(rbind,clipping),file.path(out,'ridge_clipping_sensitivity.csv'),row.names=FALSE)
 roll<-list();rmem<-list();rtheta<-list()
 for(t in 2008:2017){
  train<-d[d$Year<t & !is.na(d$SuddenStop_GI),];test<-d[d$Year==t,]
  th<-calibration_theta(inputs,t);bu<-fit_policy_models(train,inputs$features,seed0+t-2008L)
  roll[[as.character(t)]]<-evaluate_policy_rows(predict_policy_actions(bu,test,th),test,th,bu$seed,'rolling')
  rmem[[as.character(t)]]<-split_membership(bu,train,'rolling',t);rtheta[[as.character(t)]]<-data.frame(year=t,as.data.frame(th))
 }
 rolling<-do.call(rbind,roll);write.csv(rolling,file.path(out,'rolling_predictions.csv'),row.names=FALSE)
 write.csv(do.call(rbind,lapply(roll,function(z){s<-summarize_grouped(z);s$year<-z$Year[1];s})),file.path(out,'rolling_annual_summary.csv'),row.names=FALSE)
 pooled<-rolling;pooled$seed<-seed0;write.csv(summarize_grouped(pooled),file.path(out,'rolling_pooled_summary.csv'),row.names=FALSE)
 write.csv(do.call(rbind,rmem),file.path(out,'rolling_split_membership.csv'),row.names=FALSE)
 write.csv(do.call(rbind,rtheta),file.path(out,'rolling_theta.csv'),row.names=FALSE)
 loo<-lapply(c('none','2008-2009',as.character(2008:2017)),function(omit){
  zz<-primary;if(omit=='2008-2009')zz<-zz[!zz$Year%in%2008:2009,]else if(omit!='none')zz<-zz[zz$Year!=as.integer(omit),]
  s<-summarize_grouped(zz);s$omitted_period<-omit;s
 })
 write.csv(do.call(rbind,loo),file.path(out,'leave_period_out.csv'),row.names=FALSE)
 for(analysis in c('static','rolling')){
  z<-if(analysis=='static')primary else rolling
  boot<-bootstrap_policy(z,B=1000L)
  write.csv(boot$summary,file.path(out,paste0(analysis,'conditional_bootstrap_summary.csv')),row.names=FALSE)
  saveRDS(boot$raw,file.path(out,paste0(analysis,'conditional_bootstrap_raw.rds')))
  zz<-z;zz$seed<-seed0;paired<-paired_policy_comparisons(boot$raw,summarize_grouped(zz))
  write.csv(paired$summary,file.path(out,paste0(analysis,'_paired_policy_comparisons.csv')),row.names=FALSE)
  saveRDS(paired$raw,file.path(out,paste0(analysis,'_paired_policy_comparisons_raw.rds')))
 }
 run_economic_sensitivity(primary_bundle,te,theta,out)
 # Verify exact decomposition on every realized observation, including continuous actions.
 exact<-sm$crisis_benefit_contribution-sm$normal_cost_contribution-(sm$utility-sm$utility_low)
 stopifnot(max(abs(exact))<1e-12)
 split_summary<-do.call(rbind,lapply(split(sm,sm$method),function(z)data.frame(method=z$method[1],B=nrow(z),
  ce_bp_mean=mean(z$ce_bp),ce_bp_sd=sd(z$ce_bp),ce_bp_q025=unname(quantile(z$ce_bp,.025)),ce_bp_median=median(z$ce_bp),ce_bp_q975=unname(quantile(z$ce_bp,.975)),
  fraction_positive=mean(z$ce_bp>0),R0_mean=mean(z$R0),R1_mean=mean(z$R1),mean_reserve=mean(z$mean_reserve))))
 write.csv(split_summary,file.path(out,'split_variability_summary.csv'),row.names=FALSE)
 meta<-list(version=ACTION_VERSION,train_end=2007,eval_years=2008:2017,seeds=seed0+0:99,tiers=tiers,alphas=c(.10,.20,.25),delta=.1,
   country_bootstrap_B=1000L,moving_block_years=2L,bootstrap_conditional_on_fits=TRUE,historically_inspected=TRUE,
   n_train=nrow(tr),n_train_crisis=sum(tr$SuddenStop_GI),n_test=nrow(te),n_labeled=sum(!is.na(te$SuddenStop_GI)),n_test_crisis=sum(te$SuddenStop_GI,na.rm=TRUE),
   n_fit=length(primary_bundle$fits$SE$fit_idx),n_calibration=length(primary_bundle$fits$SE$cal_idx),
   odds_multiplier=primary_bundle$odds_multiplier,theta=theta,economic_scenarios=c('baseline','r=.03','r=.07','sigma=1','sigma=4'),
   scripts_md5=tools::md5sum(c('corrected/core.R','corrected/data.R','corrected/actions.R','corrected/policy_evaluation.R')),
   data_md5=tools::md5sum(list.files('data',pattern='\\.dta$',full.names=TRUE)),sessionInfo=sessionInfo(),max_identity_error=max(abs(exact)))
 saveRDS(meta,file.path(out,'metadata.rds'));writeLines(capture.output(str(meta)),file.path(out,'metadata.txt'))
 writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
 cat('DIRECT ACTION EMPIRICAL EVALUATION COMPLETE\n')
 invisible(list(primary=all_summaries[[1]],split=split_summary))
}
if(sys.nframe()==0L)run_policy_evaluation()
