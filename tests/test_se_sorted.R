source('corrected/data.R')
slow_se<-fit_se
stopifnot(file.exists('corrected/se_sorted.R'))
source('corrected/se_sorted.R')
set.seed(121)
for(n in c(12,25,100,1000))for(d in c(1,4,66)){
 x<-matrix(sample(c(-1,0,0,0,1,2,pi),n*d,TRUE),n,d);colnames(x)<-paste0('x',seq_len(d));y<-rep(c(0,1),length.out=n)
 stopifnot(identical(slow_se(x,y),fit_se(x,y)))
}
cat('PASS exact sorted counts on tied, constant, one-feature, and dense matrices\n')
inp<-load_inputs()
for(b in c(0L,1L,7L,17L)){
 ct<-NULL
 if(b){set.seed(202610150L+b);cc<-sort(unique(inp$data$Country_Code));ct<-sample(cc,length(cc),TRUE)}
 for(t in 1993:2017){
  fit_se<-slow_se;aa<-annual_fit(inp,t,'SE',20261015L+t+1000L*b,ct)
  fit_se<-sorted_se;bb<-annual_fit(inp,t,'SE',20261015L+t+1000L*b,ct)
  stopifnot(identical(aa,bb),identical(annual_result(aa,t,'SE'),annual_result(bb,t,'SE')))
 }
 cat('PASS all 25 annual fits and complete paths identical, bootstrap draw',b,'\n');flush.console()
}
tr<-inp$data[inp$data$Year<=2016 & !is.na(inp$data$SuddenStop_GI),]
x<-as.matrix(tr[,inp$features]);y<-tr$SuddenStop_GI
for(j in seq_len(ncol(x)))x[!is.finite(x[,j]),j]<-median(x[,j],na.rm=TRUE)
x[!is.finite(x)]<-0
print(system.time(for(i in 1:10)slow_se(x,y)))
print(system.time(for(i in 1:10)sorted_se(x,y)))
