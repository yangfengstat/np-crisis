# AdaBoost variant: disable unused rpart cross-validation explicitly.
# The immutable core retains SE, NP, reserve, and preprocessing definitions.
CORE_FIT_NP <- fit_np
fit_np_ada <- function(x,y,seed=1,d=Inf) {
 # Reuse exactly the common partition, imputation, and screening.
 f<-CORE_FIT_NP(x,y,'SE',seed,d)
 x<-as.matrix(x);storage.mode(x)<-'double'
 if(is.null(colnames(x)))colnames(x)<-paste0('x',seq_len(ncol(x)))
 xx<-x[,f$selected,drop=FALSE]
 for(j in seq_len(ncol(xx)))xx[!is.finite(xx[,j]),j]<-f$median[colnames(xx)[j]]
 set.seed(seed+1L)
 f$base<-ada::ada(x=xx[f$fit_idx,,drop=FALSE],y=as.factor(y[f$fit_idx]),iter=50,nu=.1,type='discrete',loss='exponential',control=rpart::rpart.control(xval=0))
 f$model<-'ada';f$cal_scores<-predict_score(f,xx[f$cal_idx,,drop=FALSE]);f
}
fit_np <- function(x,y,model='SE',seed=1,d=Inf) {
 if(model=='ada')fit_np_ada(x,y,seed,d)else CORE_FIT_NP(x,y,model,seed,d)
}
