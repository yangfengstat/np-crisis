source('corrected/data.R')
source('corrected/ada_no_cv.R')
input<-load_inputs();dd<-input$data
tr<-dd[dd$Year<=2007 & !is.na(dd$SuddenStop_GI),]
te<-dd[dd$Year>=2008 & dd$Year<=2017 & !is.na(dd$SuddenStop_GI),]
models<-c('SE','logistic','svm','randomforest','ada');sizes<-c(10,20,30,40,50,60,Inf)
metric<-curves<-roc<-list();idx<-0L
for(model in models)for(d in sizes){
 idx<-idx+1L;tt<-Sys.time();f<-fit_np(tr[,input$features],tr$SuddenStop_GI,model,seed=20261015,d=d)
 scores<-predict_score(f,te[,input$features]);m<-length(f$cal_scores)
 # Exact partition in alpha: each fitted threshold is constant on an interval.
 aa<-c(0,qbeta(.9,seq_len(m),m+1-seq_len(m)),1)
 vals<-lapply(aa,function(a){rule<-np_rule(f$cal_scores,min(1,a+1e-12),.1);flag<-apply_rule(rule,scores);c(specificity=mean(!flag[te$SuddenStop_GI==0]),missed=mean(!flag[te$SuddenStop_GI==1]))})
 vals<-do.call(rbind,vals)
 widths<-diff(aa);partial_widths<-pmax(0,pmin(tail(aa,-1),.25)-pmin(head(aa,-1),.25))
 metric[[idx]]<-data.frame(model=model,d_requested=if(is.infinite(d))'All'else as.character(d),d_actual=length(f$selected),np_auc=sum(widths*head(vals[,1],-1)),np_pauc=sum(partial_widths*head(vals[,1],-1)),roc_auc=score_auc(scores,te$SuddenStop_GI),m=m,delta=.1)
 curves[[idx]]<-data.frame(model=model,d_requested=if(is.infinite(d))'All'else as.character(d),alpha=aa,specificity=vals[,1],missed=vals[,2])
 if(is.infinite(d)){
  rr<-sort(unique(scores),decreasing=TRUE);rr<-c(Inf,rr,-Inf)
  roc[[model]]<-data.frame(model=model,fpr=vapply(rr,function(c)mean(scores[te$SuddenStop_GI==0]>=c),0),tpr=vapply(rr,function(c)mean(scores[te$SuddenStop_GI==1]>=c),0))
 }
 cat(model,d,'done',round(as.numeric(difftime(Sys.time(),tt,units='secs')),2),'s\n');flush.console()
}
out<-'outputs_corrected/backtest';dir.create(out,recursive=TRUE,showWarnings=FALSE)
metrics<-do.call(rbind,metric);curves<-do.call(rbind,curves);roc<-do.call(rbind,roc)
write.csv(metrics,file.path(out,'metrics.csv'),row.names=FALSE);write.csv(curves,file.path(out,'np_curves.csv'),row.names=FALSE);write.csv(roc,file.path(out,'roc_curves.csv'),row.names=FALSE)
for(nm in c('np_auc','np_pauc','roc_auc')) {
 tab<-matrix(metrics[[nm]],nrow=5,byrow=TRUE);colnames(tab)<-c(paste0('d',sizes[1:6]),'All');rownames(tab)<-models
 write.csv(tab,file.path(out,paste0(nm,'.csv')))
}
pdf(file.path(out,'figure2_corrected.pdf'),width=10,height=5)
par(mfrow=c(1,2),mar=c(4,4,2,1));cols<-c('#222222','#0072B2','#D55E00','#009E73','#CC79A7')
plot(0,0,type='n',xlim=c(0,1),ylim=c(0,1),xlab='False-alarm rate',ylab='Detected-crisis rate',main='Classic ROC');abline(0,1,col='grey80')
for(j in seq_along(models)){z<-roc[roc$model==models[j],];lines(z$fpr,z$tpr,col=cols[j],lty=j,lwd=1.6)}
plot(0,0,type='n',xlim=c(0,1),ylim=c(0,1),xlab=expression(alpha),ylab='Specificity',main=expression(paste('NP curves, ',delta,' = 0.1')));rect(0,0,.25,1,col=adjustcolor('grey80',.3),border=NA)
for(j in seq_along(models)){z<-curves[curves$model==models[j]&curves$d_requested=='All',];lines(z$alpha,z$specificity,type='s',col=cols[j],lty=j,lwd=1.6)}
legend('topleft',c('Signal extraction','Ridge logistic','SVM','Random forests','AdaBoost'),col=cols,lty=seq_along(models),lwd=1.6,bty='n',cex=.72)
dev.off()
writeLines(c(capture.output(print(metrics)),paste('training',nrow(tr),'crises',sum(tr$SuddenStop_GI),'test',nrow(te),'crises',sum(te$SuddenStop_GI))),file.path(out,'summary.txt'))
cat('BACKTEST COMPLETE\n');print(metrics)
