# Corrected NP pipeline. Higher scores indicate crisis. No legacy caches are read.
PIPELINE_VERSION <- '2026-10-15-fixed-score-np-constrained-reserves-v2.2'
RESERVE_PREDICTORS <- c('RES_BMGS','RES_BM','RES_GDP','RES_CH_GDP','RES_D_S')
clip_probability <- function(x,eps=.001) pmin(pmax(x,eps),1-eps)
insurance_parts <- function(th) {
 k<-th$pi.bar+th$delta.op
 A<-1-(th$r-th$g)/(1+th$g)*th$lambda
 B<-1-th$gamma-(1+th$r)/(1+th$g)*th$lambda
 stopifnot(is.finite(k),k>0,k<1,A>0,B>0,th$sigma>0)
 list(k=k,A=A,B=B)
}
reserve_policy <- function(q,th) {
 z<-insurance_parts(th);stopifnot(all(is.finite(q)),all(q>0 & q<1))
 eta<-(q*(1-z$k)/((1-q)*z$k))^(1/th$sigma)
 pmax(0,(eta*z$A-z$B)/(1-z$k+z$k*eta))
}
policy_utility <- function(p,q,th) {
 z<-insurance_parts(th);rho<-reserve_policy(q,th)
 cn<-z$A-z$k*rho;cs<-z$B+(1-z$k)*rho
 stopifnot(all(cn>0),all(cs>0))
 u<-function(x)if(th$sigma==1)log(x)else x^(1-th$sigma)/(1-th$sigma)
 p*u(cs)+(1-p)*u(cn)
}
policy_loss <- function(p,q,th) policy_utility(p,p,th)-policy_utility(p,q,th)
warning_share <- function(flags,eps=.001) {
 stopifnot(length(flags)>0,!anyNA(flags));clip_probability(mean(flags),eps)
}
np_rule <- function(cal_scores,alpha,delta=.1) {
 stopifnot(alpha>=0,alpha<=1,delta>0,delta<1,all(is.finite(cal_scores)))
 m<-length(cal_scores)
 if(alpha==1)return(list(mode='all_normal',threshold=Inf,k=NA_integer_,m=m,alpha=alpha,delta=delta))
 k<-which(pbinom(seq_len(m)-1,m,alpha)<=delta+1e-12)
 if(!length(k))return(list(mode='all_warning',threshold=-Inf,k=NA_integer_,m=m,alpha=alpha,delta=delta))
 k<-max(k)
 list(mode='finite',threshold=sort(cal_scores)[k],k=k,m=m,alpha=alpha,delta=delta)
}
apply_rule <- function(rule,scores) {
 stopifnot(all(is.finite(scores)))
 if(rule$mode=='all_warning')return(rep(TRUE,length(scores)))
 if(rule$mode=='all_normal')return(rep(FALSE,length(scores)))
 scores>=rule$threshold
}
choose_alpha <- function(alphas,reserves,observed,tolerance=1e-10) {
 loss<-(reserves-observed)^2;best<-min(loss)
 ii<-which(abs(loss-best)<=tolerance*max(1,abs(best)))
 aa<-sort(alphas[ii])
 list(alpha=min(aa),alpha_max=max(aa),n_ties=length(aa),argmin=paste(sprintf('%.17g',aa),collapse=';'),loss=best)
}
read_result <- function(file,metadata) {
 z<-readRDS(file)
 if(!identical(z$metadata,metadata))stop('Result metadata mismatch; recompute this output.')
 z$results
}
write_result <- function(results,file,metadata)saveRDS(list(metadata=metadata,results=results),file)
# Exact two-sample KS statistic without repeated ks.test warnings from ties.
ks_strength <- function(x,y) {
 o<-order(x);x<-x[o];y<-y[o];n1<-sum(y==1);n0<-sum(y==0)
 if(!n1||!n0)return(0)
 z<-cumsum(y==1)/n1-cumsum(y==0)/n0
 max(abs(z[!duplicated(x,fromLast=TRUE)]))
}
fit_se <- function(x,y) {
 n1<-sum(y==1);n0<-sum(y==0)
 rules<-lapply(seq_len(ncol(x)),function(j){
  v<-x[,j];cand<-unique(as.numeric(quantile(v,seq(.10,.90,.05),type=1)))
  grid<-expand.grid(cut=cand,direction=c(1,-1))
  tpr<-fpr<-numeric(nrow(grid))
  for(i in seq_len(nrow(grid))) {
   flag<-if(grid$direction[i]==1)v>=grid$cut[i]else v<=grid$cut[i]
   tpr[i]<-(sum(flag & y==1)+.5)/(n1+1)
   fpr[i]<-(sum(flag & y==0)+.5)/(n0+1)
  }
  ratio<-fpr/tpr;best<-which.min(ratio)
  c(cut=grid$cut[best],direction=grid$direction[best],weight=max(0,1/ratio[best]-1))
 })
 z<-as.data.frame(do.call(rbind,rules));rownames(z)<-colnames(x)
 if(sum(z$weight)==0)z$weight[]<-1
 z$weight<-z$weight/sum(z$weight)
 z
}
fit_np <- function(x,y,model='SE',seed=1,d=Inf) {
 x<-as.matrix(x);storage.mode(x)<-'double';y<-as.integer(y)
 stopifnot(nrow(x)==length(y),all(y%in%c(0L,1L)),sum(y==1)>=2,sum(y==0)>=2)
 if(is.null(colnames(x)))colnames(x)<-paste0('x',seq_len(ncol(x)))
 set.seed(seed);i1<-which(y==1);fit1<-sample(i1,floor(length(i1)/2));cal_idx<-setdiff(i1,fit1);fit_idx<-c(which(y==0),fit1)
 # Every transformation is fitted before inspecting calibration scores.
 med<-apply(x[fit_idx,,drop=FALSE],2,function(z){z<-z[is.finite(z)];if(length(z))median(z)else NA_real_})
 usable<-names(med)[is.finite(med)];stopifnot(length(usable)>0)
 x<-x[,usable,drop=FALSE];med<-med[usable]
 for(j in seq_len(ncol(x)))x[!is.finite(x[,j]),j]<-med[j]
 strength<-vapply(seq_len(ncol(x)),function(j)ks_strength(x[fit_idx,j],y[fit_idx]),0)
 selected<-colnames(x)[order(-strength,colnames(x))][seq_len(min(d,ncol(x)))]
 xf<-x[fit_idx,selected,drop=FALSE];yf<-y[fit_idx]
 set.seed(seed+1L)
 base<-switch(model,
  SE=fit_se(xf,yf),
  logistic={m<-glmnet::glmnet(xf,yf,family='binomial',alpha=0,lambda=.01,standardize=TRUE,thresh=1e-8,maxit=100000);if(m$jerr!=0)stop('Ridge logistic failed to converge');list(fit=m,lambda=.01,converged=TRUE)},
  svm=e1071::svm(x=xf,y=as.factor(yf),kernel='radial',scale=TRUE,probability=FALSE),
  randomforest=randomForest::randomForest(x=xf,y=as.factor(yf),ntree=500),
  ada=ada::ada(x=xf,y=as.factor(yf),iter=50,nu=.1,type='discrete',loss='exponential'),
  stop('Unknown model'))
 ans<-list(model=model,base=base,median=med,selected=selected,fit_idx=fit_idx,cal_idx=cal_idx,seed=seed)
 ans$cal_scores<-predict_score(ans,x[cal_idx,,drop=FALSE]);ans
}
predict_score <- function(f,x) {
 x<-as.matrix(x)
 if(is.null(colnames(x)))colnames(x)<-paste0('x',seq_len(ncol(x)))
 x<-x[,f$selected,drop=FALSE]
 for(j in seq_len(ncol(x)))x[!is.finite(x[,j]),j]<-f$median[colnames(x)[j]]
 z<-switch(f$model,
  SE={v<-matrix(FALSE,nrow(x),ncol(x));for(j in seq_len(ncol(x)))v[,j]<-if(f$base$direction[j]==1)x[,j]>=f$base$cut[j]else x[,j]<=f$base$cut[j];drop(v%*%f$base$weight)},
  logistic=as.numeric(predict(f$base$fit,newx=x,s=f$base$lambda,type='link')),
  svm={p<-predict(f$base,x,decision.values=TRUE);dd<-drop(attr(p,'decision.values'));if(f$base$levels[f$base$labels[1]]=='1')dd else -dd},
  randomforest=as.numeric(predict(f$base,x,type='prob')[,'1']),
  ada=as.numeric(predict(f$base,as.data.frame(x),type='prob')[,2])
 )
 stopifnot(length(z)==nrow(x),all(is.finite(z)));as.numeric(z)
}
score_auc <- function(score,y) {
 n1<-sum(y==1);n0<-sum(y==0)
 (sum(rank(score)[y==1])-n1*(n1+1)/2)/(n1*n0)
}
trapezoid <- function(x,y)sum(diff(x)*(head(y,-1)+tail(y,-1))/2)
alpha_path <- function(f,x,alphas,delta=.1) {
 s<-predict_score(f,x)
 do.call(rbind,lapply(alphas,function(a){r<-np_rule(f$cal_scores,a,delta);flag<-apply_rule(r,s);data.frame(alpha=a,share_raw=mean(flag),share=warning_share(flag),mode=r$mode,threshold=r$threshold,m=r$m)}))
}

np_candidates <- function(m,delta=.1) {
 stopifnot(m>=1,delta>0,delta<1)
 c(0,qbeta(1-delta,seq_len(m),m+1-seq_len(m)),1)
}
argmin_intervals <- function(alphas,reserves,observed,tolerance=1e-10) {
 a<-choose_alpha(alphas,reserves,observed,tolerance)
 loss<-(reserves-observed)^2;best<-min(loss)
 ii<-which(abs(loss-best)<=tolerance*max(1,abs(best)))
 upper<-c(tail(alphas,-1),1)
 a$alpha_upper<-max(upper[ii]);a$upper_included<-any(alphas[ii]==1)
 a$intervals<-paste(vapply(ii,function(i)if(alphas[i]==1)'{1}'else sprintf('[%.10f,%.10f)',alphas[i],upper[i]),''),collapse=';')
 a$argmin<-paste(sprintf('%.10f',alphas[ii]),collapse=';')
 a
}
