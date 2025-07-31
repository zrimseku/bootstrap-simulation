rename_for_paper <- function(x) {
  methods <- c("exact", "percentile", "basic", "bca", "bc", "standard", 
               "smoothed", "double", "studentized", "ttest", "wilcoxon",
               "ci_quant_param", "ci_quant_nonparam", "maritz-jarrett", 
               "chi_sq", "ci_corr_pearson", "clopper_pearson", "agresti_coull")
  renamed <- c("exact", "PB", "BB", "BCa", "BC", "B-n", "SB", "DB", "B-t", 
               "t-test", "wilcoxon", "q-par", "q-nonpar", "m-j", "chi-sq", "fisher", "c-p", "a-c")
  names(x)[names(x) == "alpha"] <- "nominal_coverage"
  names(x)[names(x) == "statistic"] <- "functional"
  names(x)[names(x) == "coverage_n"] <- "failed"
  x$functional[x$functional == "percentile_5"] <- "Q(0.05)"
  x$functional[x$functional == "percentile_95"] <- "Q(0.95)"
  x$failed <- 10000 - x$failed
  x <- x[, !names(x) %in% c("abs_dist_to_exact_n", "value_n", "length_n")]
  x$method <- renamed[match(x$method, methods)]
  x
}