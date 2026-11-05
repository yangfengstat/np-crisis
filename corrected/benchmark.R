source('corrected/data.R');inp<-load_inputs();dir.create('outputs_corrected',showWarnings=FALSE)
cat(R.version.string,'\n');print(sapply(c('ada','randomForest','e1071','foreign'),function(p)as.character(packageVersion(p))))
for(m in c('SE','logistic','svm','randomforest','ada')) {
 tt<-system.time(z<-annual_fit(inp,2007,m,seed=20261015))
 cat(m,'elapsed',tt[['elapsed']],'s, d=',length(z$fit$selected),'m=',length(z$fit$cal_idx),'\n')
 a<-annual_result(z,2007,m);print(a$summary)
}
cat('data calendar1980to2017:',nrow(subset(inp$data,Year<=2017)),'countries',length(unique(inp$data$Country_Code)),'events',sum(subset(inp$data,Year<=2017)$SuddenStop_GI,na.rm=TRUE),'labeled',sum(!is.na(subset(inp$data,Year<=2017)$SuddenStop_GI)),'\n')
