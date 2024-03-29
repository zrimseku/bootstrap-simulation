# bootstrap-simulation

**Comparison of different methods for confidence interval calculation.**

## Table of contents
- [Abstract](#abstract)
- [Methods](#methods)
- [Results](#results)
- [Reproducibility](#reproducibility)

# Abstract
Standard errors, confidence intervals, hypothesis tests, and other
quantifications of uncertainty are essential to statistical practice. However,
they feature a plethora of different methods, mathematical formulas, and concepts. Could we not just replace them all with the general and relatively
easy-to-understand non-parametric bootstrap? We contribute to answering
this question with a review of related work and a simulation study of one and two-sided confidence intervals over several different sample sizes, confidence levels, data generating processes, and functionals. Results show that
double bootstrap is the overall best method and a viable alternative to typically used approaches in all but the smallest sample sizes.

# Methods
We compared the most commonly used bootstrap methods to traditional methods for confidence interval
calculation. 
We then compared their accuracy and correctness, to see where choosing bootstrap over the traditional ones 
would be a mistake. 
For a more detailed description of methods for confidence interval estimation and experiment setup we refer the reader 
to the article available on link TODO.

## Methods for calculation of confidence intervals
We used the bootstrap methods implemented in the `bootstrap-ci` library: *percentile*, *basic*, *standard*, *BC*, 
*BC_a*, *smoothed*, *studentized* and *double*. 
Their description can be seen in the library's
[repository](https://github.com/zrimseku/bootstrap-ci). 

## Experiment dimensions
To get the most general results possible, we compared all methods over many combinations of different DGP's, statistics, 
dataset sizes and coverage levels.
We used the following distributions:
- standard normal,
- uniform from $0$ to $1$,
- Laplace with $\mu = 0, b = 1$,
- beta with $\alpha = 2, \beta = 10$,
- exponential with $\lambda = 1$,
- log-normal with $\mu = 0, \sigma = 1$
- bi-normal with $\mu = \[1, 1\]^T$ 
            and $\Sigma = \[2, 0.5; 0.5, 1\]$.
  
We used samples of sizes $n \in \{4, 8, 16, 32, 64, 128, 256\}$ randomly generated from these distributions to estimate 
confidence intervals for the *mean*, *median*, *standard deviation*, *5<sup>th</sup>* and *95<sup>th</sup> percentile* and *correlation*.
We were interested in confidence levels $\alpha \in {0.025, 0.05, 0.25, 0.75, 0.95, 0.975}$.

## Framework
To compare the methods we used two criteria: *accuracy* and *correctness*. Accuracy is the more important one, telling 
us how close the method's achieved coverage is to the desired one. If two methods achieve the same accuracy, we compared
their correctness, which is calculated by the distance of each method's predicted confidence intervals to the *exact* 
intervals.

The study was done in three steps:
1. Choosing the best bootstrap method.
2. Comparing the best bootstrap method to all other methods (bootstrap and traditional ones), to see where another 
method gives better one-sided confidence interval estimations.
3. Repeating step 2. for two-sided confidence intervals.

# Results
More detailed results can again be found in the article, in chapter 3.
In short, we answered to the above steps:
1. The best general bootstrap method is the **double** bootstrap. Additionally we recommend to use the **standard** 
   bootstrap when estimating confidence intervals of extreme percentiles.
2. There is **no method** (bootstrap or traditional) that would have **significantly better accuracy** in most of the 
   repetitions for experiments on any DGP. Only for the correlation, **Fisher's method** is equally accurate but **more
   correct**.
3. Results for **two-sided intervals** have the **same conclusions**.

## Visualizations
Results for separate experiments and comparisons over several dimensions can be observed interactively on 
[this site](https://zrimseku.github.io/bootstrap-simulation/).

# Reproducibility
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

