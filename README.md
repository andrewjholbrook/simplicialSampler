# simplicialSampler
Contains code, data and results associated with paper _A simple MCMC algorithm that chooses from multiple proposals at each step_.

## File list

### R/

`simplicialSampler.R` contains all project functions

### inst/

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
