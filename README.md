# bootstrap-simulation

This repository contains the materials needed to reproduce the results from our paper **Quantifying Uncertainty: All We Need is the Bootstrap?**.

**Abstract:** Standard errors, confidence intervals, hypothesis tests, and other quantifications of uncertainty are essential to statistical practice. However, they feature a plethora of different methods, mathematical formulas, and concepts. Could we not just replace them all with the general and relatively easy-to-understand non-parametric bootstrap? We contribute to answering this question with a review of related work and a simulation study of one and two-sided confidence intervals over several different sample sizes, confidence levels, data generating processes, and functionals. Results show that double bootstrap is the overall best method and a viable alternative to typically used approaches in all but the smallest sample sizes.

## Reproducibility
If you want to reproduce our results, or are interested in adding new experiments, here are the instructions.

### Environment setup
Clone this repository, go into `bootstrap-simulation` folder and prepare the environment:
```
conda env create -f environment.yml
conda activate bootstrap
```

### Repeating the experiments
Data generating processes used are implemented in the file `generators.py`, where you can also add your custom DGP, by
extending the `DGP` class. 

To get the results of all experiments, both for non-parametric and hierarchical case, you can run the file 
`ci_comparison.py`. You can change the desired DGP's, sample sizes, functionals, methods and confidence levels in the 
main function of the file to serve your interests.

### Preprocess data
R scripts in the folder `preprocess_data` can be used to produce the aggregated data we use in the analysis from raw 
experiment data. A small subsample of the experiment data is used as an example (the entire dataset is 10GB).

Notes:
- data_raw contains a subsample of the raw experiment results.
- 0_to_rds.R should be run first, 1a_aggregate_one_sided.R and 1b_aggregate_two_sided.R can be run in any order.
- data_rds contains the processed files.

### Analyze data
The same steps for data analysis as we did can be repeated with script `analyze.R` from folder `analyze data`.

