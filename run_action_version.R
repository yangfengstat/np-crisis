# Canonical entrypoint for the direct reserve-action research version.
# Run from this directory: Rscript --vanilla run_action_version.R
# Existing outputs_action results are regenerated. No network access is used.
required <- c('foreign','glmnet','ada','rpart')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install required R packages before running: ',
                         paste(missing, collapse=', '))
if(!file.exists('corrected/policy_evaluation.R'))
 stop('Run this script from the code directory.')
rs <- file.path(R.home('bin'),'Rscript')
for(script in c('tests/run_all.R','corrected/policy_evaluation.R',
                'corrected/policy_figures.R','corrected/action_simulations.R')) {
 cat('\nRunning ',script,'\n',sep=''); flush.console()
 status <- system2(rs,c('--vanilla',script))
 if(status != 0L) stop('Failed: ',script)
}
cat('\nDIRECT RESERVE-ACTION VERSION COMPLETE\n')
