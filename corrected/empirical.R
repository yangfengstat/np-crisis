source('corrected/data.R')
args<-commandArgs(trailingOnly=TRUE);mode<-if(length(args))args[1]else 'point'
models<-if(length(args)>=2)strsplit(args[2],',',fixed=TRUE)[[1]]else 'SE'
B<-if(length(args)>=3)as.integer(args[3])else 500L
cores<-if(length(args)>=4)as.integer(args[4])else 4L
stopifnot(identical(models,'SE'))
inp<-load_inputs();years<-1993:2017;out<-'outputs_corrected/empirical';dir.create(out,recursive=TRUE,showWarnings=FALSE)
metadata<-function(model,rep=NULL)list(version=PIPELINE_VERSION,model=model,rep=rep,years=years,seed=20261015L,delta=.1,
 core_md5=unname(tools::md5sum('corrected/core.R')),data_code_md5=unname(tools::md5sum('corrected/data.R')),
 data_md5=unname(tools::md5sum(sort(list.files('data',pattern='\\.dta$',full.names=TRUE)))),packages=sapply(c('ada','foreign'),function(p)as.character(packageVersion(p))))
if(mode=='point'){
 results<-list();ii<-0L
 for(model in models)for(t in years){
  z<-annual_fit(inp,t,model,seed=20261015L+t);a<-annual_result(z,t,model);ii<-ii+1L;results[[ii]]<-a
  cat(model,t,'alpha',a$summary$alpha,'ties',a$summary$n_ties,'observed',z$observed,'fit',a$summary$reserve,'\n');flush.console()
 }
 write.csv(do.call(rbind,lapply(results,`[[`,'summary')),file.path(out,'annual_estimates.csv'),row.names=FALSE)
 write.csv(do.call(rbind,lapply(results,`[[`,'path')),file.path(out,'reserve_matching_paths.csv'),row.names=FALSE)
 write.csv(unique(do.call(rbind,lapply(results,`[[`,'theta'))),file.path(out,'insurance_calibration.csv'),row.names=FALSE)
 write_result(results,file.path(out,'point_results_SE.rds'),lapply(models,metadata))
 d<-subset(inp$data,Year<=2017);writeLines(c(paste('calendar_rows',nrow(d)),paste('labeled_rows',sum(!is.na(d$SuddenStop_GI))),paste('crises',sum(d$SuddenStop_GI,na.rm=TRUE)),paste('countries',length(unique(d$Country_Code))),paste('prevalence_labeled',mean(d$SuddenStop_GI,na.rm=TRUE))),file.path(out,'data_counts.txt'))
}else if(mode=='bootstrap'){
 countries<-sort(unique(inp$data$Country_Code))
 for(model in models){
  dest<-file.path(out,paste0('bootstrap_',model));dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  one<-function(b){
   md<-metadata(model,b);file<-file.path(dest,sprintf('rep_%04d.rds',b))
   if(file.exists(file)){old<-tryCatch(read_result(file,md),error=function(e)NULL);if(!is.null(old))return(nrow(old))}
   set.seed(202610150L+b);ct<-sample(countries,length(countries),replace=TRUE)
   ans<-lapply(years,function(t){z<-annual_fit(inp,t,model,seed=20261015L+t+1000L*b,countries=ct);a<-annual_result(z,t,model)$summary;a$rep<-b;a})
   ans<-do.call(rbind,ans);write_result(ans,file,md);cat(model,'bootstrap',b,'done\n');flush.console();nrow(ans)
  }
  ans<-parallel::mclapply(seq_len(B),one,mc.cores=cores,mc.preschedule=FALSE)
  if(any(vapply(ans,inherits,TRUE,'try-error')))stop('Bootstrap worker failed; inspect logs.')
  raw<-do.call(rbind,lapply(seq_len(B),function(b)read_result(file.path(dest,sprintf('rep_%04d.rds',b)),metadata(model,b))))
  write.csv(raw,file.path(out,paste0('bootstrap_',model,'_raw.csv')),row.names=FALSE)
  sm<-do.call(rbind,lapply(split(raw,raw$year),function(z)data.frame(year=z$year[1],model=model,B=nrow(z),mean=mean(z$alpha),sd=sd(z$alpha),q05=quantile(z$alpha,.05),q95=quantile(z$alpha,.95),tie_fraction=mean(z$n_ties>1),boundary_fraction=mean(z$alpha%in%c(0,1)))))
  write.csv(sm,file.path(out,paste0('bootstrap_',model,'_summary.csv')),row.names=FALSE)
 }
}else stop('mode must be point or bootstrap')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'));cat('EMPIRICAL',mode,'COMPLETE\n')
