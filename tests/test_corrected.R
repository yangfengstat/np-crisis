# Regression tests target scientific behavior and isolation, not file formatting.
dir.create('outputs_corrected',showWarnings=FALSE)
if(file.exists('corrected/core.R')) source('corrected/core.R') else source('np-crisis-helper-funcs-r2.R')
nfail <- 0L
check <- function(name,expr) { ans<-tryCatch({force(expr);TRUE},error=function(e){cat('FAIL',name,':',conditionMessage(e),'\n');FALSE});if(ans)cat('PASS',name,'\n')else nfail<<-nfail+1L }
th <- list(sigma=2,r=.05,g=.034,lambda=.069,gamma=.065,pi.bar=.042,delta.op=.011)
check('reserve optimizer equals direct constrained utility optimization',{
  A<-1-(th$r-th$g)/(1+th$g)*th$lambda;B<-1-th$gamma-(1+th$r)/(1+th$g)*th$lambda;k<-th$pi.bar+th$delta.op
  for(p in c(.001,.02,.04,.1,.5,.95)) {
    u<-function(rho) p*(B+(1-k)*rho)^(-1)*(-1)+(1-p)*(A-k*rho)^(-1)*(-1)
    opt<-optimize(u,c(0,(A/k)*(1-1e-9)),maximum=TRUE,tol=1e-10)$maximum
    got<-if(exists('reserve_policy'))reserve_policy(p,th)else rho.star(p,th$lambda,th$pi.bar,th$gamma,th$g,th$delta.op,th$r,th$sigma)
    stopifnot(abs(got-opt)<1e-5)
  }
})
check('calibration preprocessing and base fit ignore held-out features',{
 stopifnot(exists('fit_np'))
 set.seed(8);x<-matrix(rnorm(400*5),400,5);colnames(x)<-paste0('x',1:5);y<-rep(0:1,each=200)
 f<-fit_np(x,y,'SE',seed=42,d=3);x2<-x;x2[f$cal_idx,]<-x2[f$cal_idx,]+10000
 g<-fit_np(x2,y,'SE',seed=42,d=3)
 stopifnot(length(intersect(f$fit_idx,f$cal_idx))==0,all(y[f$cal_idx]==1),identical(f$base,g$base),identical(f$median,g$median),identical(f$selected,g$selected))
})
check('infeasible alpha flags every observation and alpha one flags none',{
 stopifnot(exists('np_rule'));a<-np_rule(1:4,.05,.1);b<-np_rule(1:4,1,.1)
 stopifnot(a$mode=='all_warning',all(apply_rule(a,c(-100,0,100))),!any(apply_rule(b,c(-100,0,100))))
})
check('discrete score ties use inclusive warnings',{
 stopifnot(exists('np_rule'));a<-np_rule(rep(1,100),.1,.1);stopifnot(all(apply_rule(a,rep(1,100))))
})
check('warning share clips the final mean only',{
 stopifnot(exists('warning_share'));stopifnot(abs(warning_share(c(TRUE,rep(FALSE,52)))-1/53)<1e-14,warning_share(rep(FALSE,53))==.001)
})
check('reserve argmin reports all ties and smallest-alpha selection',{
 stopifnot(exists('choose_alpha'));a<-choose_alpha(c(.1,.2,.3),c(.4,.4,.6),.4);stopifnot(a$alpha==.1,a$n_ties==2,a$alpha_max==.2,isTRUE(all.equal(as.numeric(strsplit(a$argmin,';',fixed=TRUE)[[1]]),c(.1,.2))))
})
check('cache identity includes model and refuses an incompatible payload',{
 stopifnot(exists('read_result'));f<-tempfile(tmpdir='outputs_corrected',fileext='.rds');saveRDS(list(metadata=list(model='ada'),results=1),f)
 ans<-tryCatch(read_result(f,list(model='SE')),error=function(e)e);unlink(f);stopifnot(inherits(ans,'error'))
})
check('actual constrained welfare cost is nonnegative including boundary',{
 stopifnot(exists('policy_loss'));z<-outer(c(.005,.02,.04,.1,.5),c(.001,.02,.1,.9),Vectorize(function(p,q)policy_loss(p,q,th)));stopifnot(all(z>=-1e-10))
})
cat('Failures:',nfail,'\n');quit(status=if(nfail)1 else 0)
