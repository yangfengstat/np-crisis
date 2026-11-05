# Complete corrected submission pipeline, run from code/ with --vanilla.
# Optional argument: bootstrap worker count (default 4).
args<-commandArgs(trailingOnly=TRUE)
cores<-if(length(args))as.integer(args[1])else 4L
stopifnot(length(cores)==1L,is.finite(cores),cores>=1)
dir.create('outputs_corrected',showWarnings=FALSE)
rs<-file.path(R.home('bin'),'Rscript')
run<-function(script,args=character(),log){
 cat('Running',script,paste(args,collapse=' '),'\n');flush.console()
 status<-system2(rs,c('--vanilla',script,args),stdout=log,stderr=log)
 if(status!=0)stop('Analysis failed: ',script,'; see ',log)
}
run('tests/run_all.R',log='outputs_corrected/tests_all_final.log')
run('corrected/backtest.R',log='outputs_corrected/backtest_run.log')
run('corrected/empirical.R',c('point','SE'),log='outputs_corrected/empirical_SE_point.log')
run('corrected/empirical_ada.R',c('point','ada'),log='outputs_corrected/empirical_ada_point.log')
run('corrected/empirical.R',c('bootstrap','SE','500',cores),log='outputs_corrected/bootstrap_SE_run.log')
run('corrected/empirical_ada.R',c('bootstrap','ada','200',cores),log='outputs_corrected/bootstrap_ada_run.log')
run('corrected/figures.R',log='outputs_corrected/empirical_figures.log')
run('corrected/simulations.R','full',log='outputs_corrected/simulations_full_run.log')
run('corrected/validate_outputs.R',log='outputs_corrected/output_validation.log')
cat('ALL CORRECTED SUBMISSION COMPUTATIONS COMPLETE\n')
