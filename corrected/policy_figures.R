# Vector publication figures, generated directly from saved numerical outputs.
out<-'outputs_action/empirical';fig<-file.path(out,'figures');dir.create(fig,showWarnings=FALSE)
s<-read.csv(file.path(out,'primary_summary.csv'));boot<-read.csv(file.path(out,'staticconditional_bootstrap_summary.csv'))
sp<-read.csv(file.path(out,'static_all_splits_summary.csv'));roll<-read.csv(file.path(out,'rolling_annual_summary.csv'))
order_methods<-c('always_low','always_high','historical_constant','probability_bounded','SE_NP_0.10','SE_NP_0.20','SE_NP_0.25','SE_economic','SE_equal_error','logistic_NP_0.10','logistic_NP_0.20','logistic_NP_0.25','logistic_economic','logistic_equal_error')
labels<-c('Always low','Always high','Historical constant','Ridge probability (bounded)','SE NP 0.10','SE NP 0.20','SE NP 0.25','SE economic','SE equal error','Ridge NP 0.10','Ridge NP 0.20','Ridge NP 0.25','Ridge economic','Ridge equal error')
colors<-ifelse(grepl('^SE',order_methods),'#167C80',ifelse(grepl('^logistic',order_methods),'#B14F32','#536273'))
pdf(file.path(fig,'policy_value_primary.pdf'),width=8.5,height=6.8,family='Helvetica',useDingbats=FALSE)
par(mar=c(4.2,11.7,2.2,1.2),las=1,mgp=c(2.7,.6,0),cex=.88)
x<-s$ce_bp[match(order_methods,s$method)];y<-rev(seq_along(x));bc<-boot[boot$resampling=='country',];bt<-boot[boot$resampling=='two_year_moving_block',]
bc<-bc[match(order_methods,bc$method),];bt<-bt[match(order_methods,bt$method),]
xlim<-range(c(x,bc$ce_bp_lo,bc$ce_bp_hi,bt$ce_bp_lo,bt$ce_bp_hi),finite=TRUE);xlim<-xlim+c(-1,1)*diff(xlim)*.04
plot(x,y,type='n',ylim=c(.5,length(x)+.5),xlim=xlim,yaxt='n',ylab='',xlab='Consumption-equivalent gain relative to always low (basis points)',main='Retrospective policy value, 2008-2017')
abline(v=0,col='gray65',lty=2);abline(h=y,col='#ECEFF1')
segments(bc$ce_bp_lo,y+.11,bc$ce_bp_hi,y+.11,col=colors,lwd=2)
segments(bt$ce_bp_lo,y-.11,bt$ce_bp_hi,y-.11,col=colors,lwd=2,lty=3)
points(x,y,pch=19,col=colors,cex=.8);axis(2,at=y,labels=labels,tick=FALSE)
legend('bottomleft',c('Country bootstrap','Two-year block bootstrap'),lty=c(1,3),lwd=2,bty='n',cex=.77)
dev.off()
pdf(file.path(fig,'policy_split_sensitivity.pdf'),width=8.5,height=6.8,family='Helvetica',useDingbats=FALSE)
par(mar=c(4.2,11.7,2.2,1.2),las=1,mgp=c(2.7,.6,0),cex=.88)
vals<-lapply(rev(order_methods),function(m)sp$ce_bp[sp$method==m])
boxplot(vals,horizontal=TRUE,names=rev(labels),las=1,col=rev(colors),border=rev(colors),outline=FALSE,
 xlab='Consumption-equivalent gain relative to always low (basis points)',main='Variation across 100 common fitting/calibration splits',boxwex=.55)
abline(v=0,col='gray55',lty=2);points(x,y,pch=23,bg='white',cex=.8)
dev.off()
pdf(file.path(fig,'policy_error_tradeoff.pdf'),width=7.4,height=5.3,family='Helvetica',useDingbats=FALSE)
par(mar=c(4.3,4.6,2.3,1.2),las=1,mgp=c(2.7,.6,0),cex=.9)
plot(NA,xlim=c(0,1),ylim=c(-.03,1),xlab='False-alarm rate R0',ylab='Missed-crisis rate R1',main='Primary split: realized conditional error rates')
abline(h=seq(0,1,.2),v=seq(0,1,.2),col='#ECEFF1')
# Score families differ by line type and filled/open symbols as well as color.
# NP labels are placed by hand so that they clear the points and both lines.
lab<-list(SE=list(dx=c(.012,.052,-.040),dy=c(.042,.036,.034)),logistic=list(dx=c(-.040,0,0),dy=c(-.036,-.036,-.036)))
for(m in c('SE','logistic')){
 se<-m=='SE';col<-if(se)'#167C80'else '#B14F32'
 z<-s[grepl(paste0('^',m,'_NP'),s$method),];z<-z[order(z$method),]
 lines(z$R0,z$R1,col=col,lwd=1.7,lty=if(se)1 else 2);points(z$R0,z$R1,col=col,pch=if(se)19 else 21,bg='white',lwd=1.4)
 text(z$R0+lab[[m]]$dx,z$R1+lab[[m]]$dy,labels=sub('.*_NP_','',z$method),col=col,cex=.77)
 ze<-s[s$method==paste0(m,'_economic'),];points(ze$R0,ze$R1,col=col,pch=if(se)17 else 24,bg='white',cex=1.1,lwd=1.4)
 zb<-s[s$method==paste0(m,'_equal_error'),];points(zb$R0,zb$R1,col=col,pch=0,cex=2,lwd=1.2)
}
legend('topright',c('SE','Ridge logistic'),col=c('#167C80','#B14F32'),lty=c(1,2),pch=c(19,21),pt.bg='white',lwd=1.7,bty='n')
legend('bottomleft',c('NP (label is missed-crisis bound)','Economic threshold','Equal-error threshold'),pch=c(19,17,0),pt.cex=c(1,1,1.4),bty='n',cex=.8)
dev.off()
pdf(file.path(fig,'policy_rolling_replay.pdf'),width=8.5,height=4.7,family='Helvetica',useDingbats=FALSE)
par(mar=c(4.2,4.5,2.3,1.2),las=1,mgp=c(2.7,.6,0),cex=.88)
mm<-c('always_high','historical_constant','probability_bounded','SE_NP_0.20','logistic_NP_0.20')
cc<-c('#929AA1','#5D4673','#4774AA','#167C80','#B14F32');yl<-range(roll$ce_bp[roll$method%in%mm])
plot(NA,xlim=c(2008,2017),ylim=yl,xlab='Evaluation outcome year',ylab='Consumption-equivalent gain (basis points)',main='Expanding-window replay with annual economic calibration',xaxt='n')
axis(1,at=2008:2017,labels=2008:2017);abline(h=0,col='gray65',lty=2)
for(i in seq_along(mm)){z<-roll[roll$method==mm[i],];z<-z[order(z$year),];lines(z$year,z$ce_bp,type='b',pch=i,lwd=1.4,col=cc[i])}
legend('topright',c('Always high','Historical constant','Probability bounded','SE NP 0.20','Ridge NP 0.20'),col=cc,lty=1,pch=1:5,bty='n',cex=.76)
dev.off()
cat('FOUR POLICY FIGURES COMPLETE\n')
