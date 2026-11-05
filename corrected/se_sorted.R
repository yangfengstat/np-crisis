# Exact integer-count implementation of core.R fit_se, without changing the model.
# For direction +1, # >= cutoff is n - # < cutoff; for -1, use # <= cutoff.
sorted_se <- function(x,y) {
 n1<-sum(y==1);n0<-sum(y==0)
 rules<-lapply(seq_len(ncol(x)),function(j){
  v<-x[,j];cand<-unique(as.numeric(quantile(v,seq(.10,.90,.05),type=1)))
  grid<-expand.grid(cut=cand,direction=c(1,-1))
  s1<-sort(v[y==1]);s0<-sort(v[y==0])
  c1<-c(n1-findInterval(cand,s1,left.open=TRUE),findInterval(cand,s1))
  c0<-c(n0-findInterval(cand,s0,left.open=TRUE),findInterval(cand,s0))
  tpr<-(c1+.5)/(n1+1);fpr<-(c0+.5)/(n0+1)
  ratio<-fpr/tpr;best<-which.min(ratio)
  c(cut=grid$cut[best],direction=grid$direction[best],weight=max(0,1/ratio[best]-1))
 })
 z<-as.data.frame(do.call(rbind,rules));rownames(z)<-colnames(x)
 if(sum(z$weight)==0)z$weight[]<-1
 z$weight<-z$weight/sum(z$weight)
 z
}
fit_se<-sorted_se
