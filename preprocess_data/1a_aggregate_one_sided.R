setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("util.R")

df <- readRDS("./data_rds/intervals.rds")

methods <- names(df)[8:ncol(df)]

res <- NULL

for (method in methods) {
  print(method)
  df$value <- df[, names(df) == method]
  df$coverage <- ifelse(df$value >= df$true_value, 1, 0)
  df$abs_dist_to_exact <- abs(df$exact - df$value)
  y <- aggregate(cbind(coverage, abs_dist_to_exact, value) ~ alpha + dgp + statistic + n + B, FUN = mean, data = df, na.rm=TRUE, na.action=NULL)
  y_sd <- aggregate(cbind(coverage, abs_dist_to_exact, value) ~ alpha + dgp + statistic + n + B, FUN = sd, data = df, na.rm=TRUE, na.action=NULL)
  y_n <- aggregate(cbind(coverage, abs_dist_to_exact, value) ~ alpha + dgp + statistic + n + B, FUN = function(x){sum(!is.na(x))}, data = df, na.action=NULL)
  names(y_n) <- paste(names(y_n), "_n", sep = "")
  names(y_sd) <- paste(names(y_sd), "_sd", sep = "")  
  y <- cbind(y, y_sd[,6:8], y_n[,6:8])
  y$method <- method
  res <- rbind(res, y)
}

saveRDS(rename_for_paper(res), "./data_rds/aggregated_one_sided.rds")