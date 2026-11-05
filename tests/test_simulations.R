# Run from code/: Rscript --vanilla tests/test_simulations.R
# Catches use of <= at SE score atoms and incorrect Bernoulli weighting.
if (file.exists('corrected/simulations.R')) source('corrected/simulations.R')
stopifnot(exists('weighted_score_cdf', mode='function'))
z <- weighted_score_cdf(c(-Inf,0,.25,.5,.75,1,Inf),c(.25,.75),c(.5,.5))
stopifnot(isTRUE(all.equal(z,c(0,0,.25,.5,.5,.75,1),tolerance=1e-13)))
z <- weighted_score_cdf(c(.25,.5,.75,1),c(.25,.75),c(.75,.5))
stopifnot(isTRUE(all.equal(z,c(.125,.5,.5,.625),tolerance=1e-13)))
# A zero-weight rule must not alter the score distribution.
stopifnot(identical(weighted_score_cdf(.5,c(.25,.75,0),c(.75,.5,.9)),.5))
# Independently enumerate ten binary indicators with exactly representable
# dyadic weights, including all score atoms, against the split CDF algorithm.
bits<-as.matrix(expand.grid(rep(list(0:1),10)))
w<-c(1,3,5,9,17,31,65,129,257,513)/2048
pr<-rep(c(.25,.5,.75,.5,.25),2)
vals<-drop(bits%*%w)
mass<-apply(bits,1,function(v)prod(ifelse(v==1,pr,1-pr)))
cuts<-sort(unique(vals))
brute<-vapply(cuts,function(v)sum(mass[vals<v]),0)
stopifnot(max(abs(weighted_score_cdf(cuts,w,pr)-brute))<1e-13)
# Catches collapsing distinct floating-point score atoms, such as .1+.2 and
# .3. Expected errors use the same literal matrix score as predict_score.
bits_fp<-as.matrix(expand.grid(rep(list(0:1),5)))
w_fp<-c(.1,.2,.3,1e-17,0)
values_fp<-drop(bits_fp%*%w_fp)
cuts_fp<-sort(unique(values_fp))
expected_fp<-vapply(cuts_fp,function(v)mean(values_fp<v),0)
stopifnot(max(abs(weighted_score_cdf(cuts_fp,w_fp,rep(.5,5))-expected_fp))<1e-13)
stopifnot(identical(weighted_score_cdf(c(-1e-17,0,1e-17),c(0,0),c(.1,.8)),c(0,0,1)))
# Catches a class-mean sign reversal or ignoring rule direction.
f <- list(model='SE',selected=c('x1','x2'),base=data.frame(cut=c(0,0),direction=c(1,-1),weight=c(.25,.75)))
p <- list(mu1=c(x1=qnorm(.75),x2=0),d=2,s=1,pi0=.04)
e <- se_population_errors(f,c(.25,.5,.75,1),p)
stopifnot(isTRUE(all.equal(e$R1,c(.125,.5,.5,.625),tolerance=1e-13)))
stopifnot(isTRUE(all.equal(e$R0,c(.75,.5,.5,.25),tolerance=1e-13)))
# Catches treating a minimizing alpha interval as only its lower endpoint.
stopifnot(exists('in_minimizing_intervals',mode='function'))
stopifnot(in_minimizing_intervals(c(0,.9,1),c(1,.2,0),.2,.95),
          !in_minimizing_intervals(c(0,.9,1),c(1,.2,0),.2,1),
          in_minimizing_intervals(c(0,.9,1),c(1,.2,0),0,1))
# Exercise actual fit, endpoint fallbacks, deterministic seeds, and output coverage.
cfg <- simulation_config('smoke')
a <- panel_a_rep(1,1000,c(0,.1,1),cfg)
stopifnot(nrow(a)==6,all(a$R1[a$alpha==0]==0),all(a$R0[a$alpha==0]==1),
          all(a$R1[a$alpha==1]==1),all(a$R0[a$alpha==1]==0),
          identical(a,panel_a_rep(1,1000,c(0,.1,1),cfg)))
c <- recovery_rep(1,cfg,periods=c(1,13,25),scenarios='baseline')
stopifnot(nrow(c)==3,all(c$alpha_hat>=0 & c$alpha_hat<=1),all(c$n_ties>=1),
          all(c$rho_observed>=0),all(c$alpha_max>=c$alpha_hat))
cat('Simulation tests passed: exact discrete CDF, SE class errors, endpoints, reproducibility, recovery integration.\n')
