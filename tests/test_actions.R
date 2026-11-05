# Independent checks reject sign swaps, omitted holding costs, boundary leakage,
# and incorrect CRRA consumption-equivalent transformations.
source('corrected/core.R')
if(file.exists('corrected/actions.R'))source('corrected/actions.R')
nfail<-0L
check<-function(name,expr){ok<-tryCatch({force(expr);TRUE},error=function(e){cat('FAIL',name,':',conditionMessage(e),'\n');FALSE});if(ok)cat('PASS',name,'\n')else nfail<<-nfail+1L}
th<-list(sigma=2,r=.05,g=.034,lambda=.069,gamma=.065,pi.bar=.042,delta.op=.011)
check('direct budgets and exact finite-sample two-action identity',{
 y<-c(1,0,0,1,0,0,0,1,0,0);flag<-c(1,1,0,0,0,1,0,1,0,0)==1
 low<-.05;high<-.15;k<-th$pi.bar+th$delta.op
 A<-1-(th$r-th$g)/(1+th$g)*th$lambda;B<-1-th$gamma-(1+th$r)/(1+th$g)*th$lambda
 u<-function(c)-1/c
 b<-u(B+(1-k)*high)-u(B+(1-k)*low);a<-u(A-k*low)-u(A-k*high)
 z<-action_components(th,low,high);stopifnot(abs(z$a-a)<1e-14,abs(z$b-b)<1e-14,a>0,b>0)
 direct<-u(ifelse(y==1,B+(1-k)*ifelse(flag,high,low),A-k*ifelse(flag,high,low)))
 got<-realized_utility(y,ifelse(flag,high,low),th)
 base<-u(ifelse(y==1,B+(1-k)*low,A-k*low))
 R0<-mean(flag[y==0]);R1<-mean(!flag[y==1]);p<-mean(y)
 stopifnot(max(abs(got-direct))<1e-14,abs(mean(got-base)-(p*b*(1-R1)-(1-p)*a*R0))<1e-14)
})
check('feasible reserve boundaries reject zero consumption and negative reserves',{
 z<-insurance_parts(th)
 stopifnot(is.finite(realized_utility(0,0,th)),is.finite(realized_utility(0,z$A/z$k*(1-1e-10),th)))
 for(r in c(-.001,z$A/z$k,z$A/z$k+1))stopifnot(inherits(try(realized_utility(1,r,th),silent=TRUE),'try-error'))
 stopifnot(inherits(try(action_components(th,.2,.1),silent=TRUE),'try-error'))
})
check('log utility and consumption-equivalent gain use exact scale identities',{
 for(s in c(1,2,.5,4)){
  tt<-th;tt$sigma<-s
  y<-c(0,1,0,1);rho<-c(.05,.15,.15,.05);z<-insurance_parts(tt)
  c<-ifelse(y==1,z$B+(1-z$k)*rho,z$A-z$k*rho)
  u<-function(x)if(s==1)log(x)else x^(1-s)/(1-s)
  stopifnot(max(abs(realized_utility(y,rho,tt)-u(c)))<1e-14)
  stopifnot(abs(ce_gain(mean(u(c*1.025)),mean(u(c)),s)-.025)<1e-12)
 }
})
if(file.exists('corrected/policy_evaluation.R'))source('corrected/policy_evaluation.R')
check('weighted thresholds recover target population loss with unequal sampling',{
 sc<-c(.1,.3,.5,.7,.9);yy<-c(0,0,0,1,1);ww<-c(1,1,1,4,4)
 z<-choose_policy_threshold(sc,yy,ww,a=.2,b=1,criterion='economic')
 cand<-c(-Inf,sort(unique(sc)),Inf)
 loss<-vapply(cand,function(t)sum(ww*(.2*(yy==0 & sc>=t)+1*(yy==1 & sc<t)))/sum(ww),0)
 stopifnot(abs(z$loss-min(loss))<1e-14,z$threshold==max(cand[abs(loss-min(loss))<1e-12]))
})
check('future outcome changes cannot change frozen fit predictions or actions',{
 set.seed(18);tr<-data.frame(x1=rnorm(120),x2=rnorm(120),SuddenStop_GI=rep(0:1,each=60))
 te<-data.frame(x1=c(-2,0,2),x2=c(1,-1,0),SuddenStop_GI=c(0,1,0))
 bundle<-fit_policy_models(tr,c('x1','x2'),seed=172)
 p1<-predict_policy_actions(bundle,te,th,.05,.15)
 te$SuddenStop_GI<-c(1,0,1)
 p2<-predict_policy_actions(bundle,te,th,.05,.15)
 stopifnot(identical(p1,p2),identical(bundle$fits$SE$fit_idx,bundle$fits$logistic$fit_idx),
           length(intersect(bundle$fits$SE$fit_idx,bundle$fits$SE$cal_idx))==0)
 # Changing calibration features must not change economic or equal-error thresholds.
 tr2<-tr;tr2[bundle$fits$SE$cal_idx,c('x1','x2')]<-10000
 bundle2<-fit_policy_models(tr2,c('x1','x2'),seed=172)
 q<-predict_policy_actions(bundle2,te,th,.05,.15)
 ix<-grepl('economic|equal_error',p1$method)
 stopifnot(identical(p1[ix,c('method','rho','threshold')],q[ix,c('method','rho','threshold')]))
})
check('continuous comparator respects the common action interval',{
 set.seed(10);tr<-data.frame(x1=rnorm(120),x2=rnorm(120),SuddenStop_GI=rep(0:1,each=60))
 b<-fit_policy_models(tr,c('x1','x2'),seed=18)
 te<-tr[1:8,];p<-predict_policy_actions(b,te,th,.05,.15)
 bd<-p[p$method=='probability_bounded',];un<-p[p$method=='probability_unrestricted',]
 stopifnot(nrow(bd)==nrow(te),nrow(un)==nrow(te),all(bd$rho>=.05 & bd$rho<=.15),
           max(abs(bd$rho-pmin(.15,pmax(.05,un$rho))))<1e-14)
})
check('historical alternate split ridge solver converges at fixed target penalty',{
 inp<-load_inputs();tr<-subset(inp$data,Year<=2007 & !is.na(SuddenStop_GI))
 b<-fit_policy_models(tr,inp$features,seed=20261016)
 stopifnot(b$fits$logistic$base$fit$jerr==0,b$fits$logistic$base$lambda==.01,
           all(is.finite(predict_score(b$fits$logistic,as.matrix(tr[,inp$features])))))
})
check('bootstrap and paired CE retain log utility and reject mixed risk aversion',{
 units<-expand.grid(Country_Code=c('A','B'),Year=2008:2009)
 scales<-c(always_high=1.10,SE_economic=1.05,SE_NP_0.20=1.02)
 rows<-do.call(rbind,lapply(names(scales),function(m)
  data.frame(units,method=m,utility=log(scales[m]),u_low=0,sigma=1)))
 boot<-bootstrap_policy(rows,B=4,seed=19)
 expected<-10000*(scales[boot$raw$method]-1)
 stopifnot(max(abs(boot$raw$ce_bp-expected))<1e-10)
 points<-data.frame(method=names(scales),utility=log(scales))
 paired<-paired_policy_comparisons(boot$raw,points)
 expected_pair<-10000*(1.02/scales[paired$summary$comparator]-1)
 stopifnot(max(abs(paired$summary$point_relative_ce_bp-expected_pair))<1e-10,
           max(abs(paired$summary$relative_ce_bp_lo-expected_pair))<1e-10,
           max(abs(paired$summary$relative_ce_bp_hi-expected_pair))<1e-10)
 rows$sigma[1]<-2
 stopifnot(inherits(try(bootstrap_policy(rows,B=1),silent=TRUE),'try-error'))
})
cat('Failures:',nfail,'\n');quit(status=as.integer(nfail>0))
