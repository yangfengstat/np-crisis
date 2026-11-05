source('corrected/core.R');set.seed(17);x<-matrix(rnorm(400*5),400,5);y<-rbinom(400,1,.2)
f<-fit_np(x,y,'ada',seed=17);stopifnot(length(predict_score(f,x))==400);cat('PASS AdaBoost model scores new observations\n')
