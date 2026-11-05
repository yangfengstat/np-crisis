# Run each test in a fresh R session so test-local monkey patches cannot leak.
files<-c('test_corrected.R','test_breakpoints.R','test_ada_prediction.R',
         'test_logistic_stability.R','test_ada_no_cv.R','test_data_alignment.R',
         'test_simulations.R','test_se_sorted.R','test_actions.R',
         'test_action_simulations.R')
rs<-file.path(R.home('bin'),'Rscript')
for(f in files){
 cat('\nRunning',f,'\n');flush.console()
 status<-system2(rs,c('--vanilla',file.path('tests',f)))
 if(status!=0)stop('Failed test: ',f)
}
cat('\nALL DIRECT-ACTION AND INHERITED PIPELINE TESTS PASSED\n')
