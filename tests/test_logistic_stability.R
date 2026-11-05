source('corrected/core.R');set.seed(4);y<-rep(0:1,each=30);x<-matrix(rnorm(60*65),60,65);x[,1]<-y*100
f<-fit_np(x,y,'logistic',seed=11);stopifnot(isTRUE(f$base$converged),all(is.finite(predict_score(f,x))));cat('PASS ridge logistic converges on separated high-dimensional data\n')
