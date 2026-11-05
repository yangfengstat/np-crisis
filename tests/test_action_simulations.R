# Run from code/: Rscript --vanilla tests/test_action_simulations.R
# Independent arithmetic and quadrature catch class reversal, utility signs,
# order-statistic indexing, and confusion between NP and welfare objectives.
if(!file.exists('corrected/action_simulations.R'))
  stop('FAIL: direct-action simulation implementation is absent')
source('corrected/action_simulations.R')
failures<-0L
check<-function(name,expr) {
  tryCatch({force(expr);cat('PASS:',name,'\n')},error=function(e){
    failures<<-failures+1L;cat('FAIL:',name,conditionMessage(e),'\n')
  })
}
near<-function(x,y,tol=1e-10)stopifnot(all(is.finite(x)),all(abs(x-y)<tol))
th<-list(sigma=2,r=.05,g=.03,lambda=.07,gamma=.065,pi.bar=.04,delta.op=.01)
check('sigma 2 action utility agrees with hand arithmetic',{
  # k=.05, A=1-.02*.07/1.03, B=.935-1.05*.07/1.03.
  z<-action_components(th,.05,.15)
  near(z$u_n_low,-1/(1-.02*.07/1.03-.05*.05))
  near(z$u_s_high,-1/(.935-1.05*.07/1.03+.95*.15))
  near(z$a,1/(1-.02*.07/1.03-.05*.15)-1/(1-.02*.07/1.03-.05*.05))
  near(z$b,1/(.935-1.05*.07/1.03+.95*.05)-1/(.935-1.05*.07/1.03+.95*.15))
})
check('sigma 1 and CE scaling use proportional consumption',{
  tl<-th;tl$sigma<-1
  z<-action_components(tl,.05,.15)
  near(z$a,log((1-.02*.07/1.03-.05*.05)/(1-.02*.07/1.03-.05*.15)))
  near(ce_gain(log(1.1)+.2,.2,1),.1)
  near(ce_gain(-2/1.1,-2,2),.1)
})
check('Gaussian risk has correct class orientation and limits',{
  z<-gaussian_action_rates(c(-Inf,0,1.5,Inf),1.5)
  near(z$R0,c(1,.5,.0668072012688581,0))
  near(z$R1,c(0,.0668072012688581,.5,1))
})
check('population welfare identity equals independent density integration',{
  for(sig in c(1,2))for(p in c(.05,.2))for(cut in c(-.8,.4,2.1)){
    tx<-th;tx$sigma<-sig
    z<-gaussian_action_value(cut,p,tx,.05,.15,1.5)
    A<-1-.02*.07/1.03;B<-.935-1.05*.07/1.03
    u<-function(x)if(sig==1)log(x)else -1/x
    # Integrate constant action utility separately on either side of cutoff.
    integ<-function(y,lo,hi,rho)integrate(function(x){
      cns<-if(y==0)A-.05*rho else B+.95*rho
      u(cns)*dnorm(x,1.5*y)
    },lo,hi,rel.tol=1e-12)$value
    truth<-(1-p)*(integ(0,-Inf,cut,.05)+integ(0,cut,Inf,.15))+
      p*(integ(1,-Inf,cut,.05)+integ(1,cut,Inf,.15))
    base<-(1-p)*u(A-.05*.05)+p*u(B+.95*.05)
    near(z$utility,truth);near(z$gain,truth-base)
    near(z$gain,z$identity_gain)
    near(z$normal_excess_reserves,.1*(1-p)*pnorm(cut,lower.tail=FALSE))
    near(z$missed_crisis_shortfall,.1*p*pnorm(cut-1.5))
  }
})
check('empirical cutoff maximizes threshold with allowed misses including ties',{
  near(empirical_action_cutoff(c(0,0,1,2,2),.4),1)
  near(empirical_action_cutoff(c(0,0,1,2,2),.2),0)
  near(empirical_action_cutoff(c(0,1,2,3,4),.4),2)
  stopifnot(is.infinite(empirical_action_cutoff(1:5,1)))
})
check('NP index achieves the advertised probability and fallback is safe',{
  rr<-np_rule(seq_len(25),.1,.1);stopifnot(rr$k==1L)
  near(pbinom(rr$k-1,25,.1),.0717897987691853)
  rr<-np_rule(seq_len(10),.1,.1);stopifnot(rr$threshold==-Inf)
})
check('one-sided binomial endpoint and finite-sample coverage are exact',{
  near(action_binomial_upper(0,20,.05),1-.05^(1/20))
  near(action_binomial_upper(20,20,.05),1)
  for(p in c(.01,.1,.5,.9)) {
    bounds<-action_binomial_upper(0:50,50,.05)
    stopifnot(sum(dbinom(0:50,50,p)[bounds>=p-1e-14])>=.95-1e-12)
  }
})
check('welfare oracle and constrained optimizer separate binding from slack',{
  z<-action_components(th,.05,.15)
  for(p in c(.05,.2)) {
    cw<-gaussian_welfare_cutoff(p,z,1.5)
    vals<-gaussian_action_value(c(cw-.05,cw,cw+.05),p,th,.05,.15,1.5)
    stopifnot(which.max(vals$utility)==2)
    for(a in c(.1,.2,.25)) {
      cc<-min(cw,1.5+qnorm(a))
      cand<-seq(-8,1.5+qnorm(a),length.out=10001)
      vv<-gaussian_action_value(cand,p,th,.05,.15,1.5)$utility
      stopifnot(gaussian_action_value(cc,p,th,.05,.15,1.5)$utility>=max(vv)-1e-11)
    }
  }
})
check('certificate is a welfare lower bound on its two claimed events',{
  a<-action_components(th,.05,.15)
  grid<-expand.grid(p=c(.05,.2),alpha=c(.1,.2,.25),R0=c(0,.2,.6,1))
  for(i in seq_len(nrow(grid))) {
    z<-grid[i,];R1<-z$alpha*.75;up<-min(1,z$R0+.05)
    gain<-z$p*a$b*(1-R1)-(1-z$p)*a$a*z$R0
    lb<-action_welfare_lower(z$p,z$alpha,up,a)
    stopifnot(gain>=lb-1e-14)
  }
})
cat('Action simulation test failures:',failures,'\n')
if(failures)quit(status=1L)
