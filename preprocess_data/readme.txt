These R scripts can be used to produce the aggregated data we use in the analysis from raw experiment data. A small subsample of the experiment data is used as an example (the entire dataset is 10GB).

Notes:
- data_raw contains a subsample of the raw experiment results.
- 0_to_rds.R should be run first, 1a_aggregate_one_sided.R and 1b_aggregate_two_sided.R can be run in any order.
- data_rds contains the processed files.
