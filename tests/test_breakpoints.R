source('corrected/core.R');stopifnot(exists('np_candidates'))
for(m in c(1,4,25,50)){
 a<-np_candidates(m,.1);stopifnot(length(a)==m+2)
 for(k in seq_len(m)){
  p<-a[k+1];stopifnot(abs(pbinom(k-1,m,p)-.1)<1e-10)
  r<-np_rule(seq_len(m),p,.1);stopifnot(r$k==k)
  mid<-(p+a[k+2])/2;r2<-np_rule(seq_len(m),mid,.1);stopifnot(r2$k==k)
 }
}
cat('PASS exact alpha breakpoint inversion and constant policy intervals\n')
