#  Data-driven Covariate Selection for Confounding Adjustment by Focusing on the Stability of the Effect Estimator

These files contain the R scripts to implement the method introduced in the paper above.

Preprint: https://psyarxiv.com/zkdqa/ (migrated from https://osf.io/yve6u/)

Slides from the [2022 meeting of the Working Group SEM](https://www.tilburguniversity.edu/about/schools/socialsciences/organization/departments/methodology-statistics/events/structural-equation-modeling): https://github.com/wwloh/covariate-selection-effect-stability/blob/main/slides-SEMlab-2022-03.pdf

## R scripts
selection-stability.R
- custom functions implementing the introduced covariate selection method

selection-other_methods.R
- wrapper functions implementing other methods for comparison

### Simulation study
selection-sims-data_generating.R
- function to generate a single dataset based on given parameters

selection-sim3.R
- simulation study implementing different methods and saving the results

selection-sims-results.R
- collating and presenting the results using figures and tables

## Implemented methods and results for the applied data examples

### Applied Data Example 1: Social interaction quantity on well-being
This example is in the "interaction_quantity_study1" folder.

example-interaction_quantity_study1.R
- downloads the data set from OSF and processes it for applying the method,
including renaming the variables for use in the custom functions
- carries out the steps for covariate selection
- plots the figures with the trajectory of the effect estimator
- implements the methods for comparison 

### Applied Data Example 2: New math curriculum on test scores 
This example is in the "math" folder.

example-math.R
- downloads the data set from OSF and processes it for applying the method,
including renaming the variables for use in the custom functions
- carries out the steps for covariate selection
- plots the figures with the trajectory of the effect estimator
- implements the methods for comparison 


## Step-by-step guide for the applied data examples

### Applied Data Example 1: Social interaction quantity on well-being
This example is in the "step_by_step-example_1" folder.
1) example-interaction_quantity_study1-1_dataprep.R
- downloads the data set from OSF and processes it for applying the method,
including renaming the variables for use in the custom functions
2) example-interaction_quantity_study1-2_stability.R
- carries out the steps for covariate selection

### Applied Data Example 2: New math curriculum on test scores 
This example is in the "step_by_step-example_2" folder.
1) example-math-1_dataprep.R
- downloads the data set from OSF and processes it for applying the method,
including renaming the variables for use in the custom functions
2) example-math-2_stability.R
- carries out the steps for covariate selection
