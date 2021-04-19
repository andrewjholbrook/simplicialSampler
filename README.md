# simplicialSampler
Contains code, data and results associated with paper _A simple MCMC algorithm that chooses from multiple proposals at each step_.

## File list

### R/

`simplicialSampler.R` contains all project functions

### inst/

`analyze_acceptance_targets.R` takes experimental output and creates left two plots of Figure 3

`analyze_performance.R` takes experimental output and creates Figure 4

`bimodalExp.R` runs experiment portrayed in Figure 5

`bimodalTest.R` tests validity of samplers on bimodal target

`fullyIllNormalExpSequential.R` runs experiment portrayed in right plots of Figure 4

`gpAnalysis.R`  takes experimental output and creates Table 1

`gpExample.R` runs experiment leading to vanilla algorithm results in Table 1

`gpExamplePrecond.R` runs experiment leading to vanilla algorithm results in Table 1

`illNormalExp.R` runs experiment leading to middle plots of Figure 4

`illNormalScaling_scaled.R` runs experiment leading to preconditioned results of left two plots from Figure 3

`illNormalScaling_unscaled.R` runs experiment leading to vanilla results of left two plots from Figure 3







#### output/
`bimodalComparison.txt` contains results presented in Figure 5

`fullyIllGaussianSequential.txt` contains results presented in the right two plots of Figure 4

`gpClassification.txt` contains vanilla algorithm results of Table 1

`gpClassificationPrecond.txt` contains preconditioned algorithm results of Table 1

`lambdaSearch.txt` contains results for left two plots of Figure 3

`nPropsResults.txt` contains results for right plot of Figure 3

`rwmComparison.txt` contains results for left plots of Figure 4

`scaledRwmComparison.txt` contains results for middle plots of Figure 4

#### inst/data/
`electionResults.rds` contains data used in GP classification example

#### inst/savedCovs/
`cov*.rds` contain randomly generated covariance matrices of `*` dimensions used in experiments leading to the right plots of Figure 4
