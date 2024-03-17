setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("util.R")

df <- readRDS("./data_rds/intervals.rds")

methods <- names(df)[8:ncol(df)]
alphas  <- c(0.025, 0.050, 0.250, 0.750, 0.950, 0.975)

instr <- NULL
instr <- rbind(instr, data.frame(alpha = 0.95, alpha_ub = 0.975, alpha_lb = 0.025))
instr <- rbind(instr, data.frame(alpha = 0.90, alpha_ub = 0.950, alpha_lb = 0.050))

res <- NULL
for (i in 1:nrow(instr)) {
  ins <- instr[i,]
  df_lb <- df[df$alpha == ins$alpha_lb,]   
  df_ub <- df[df$alpha == ins$alpha_ub,]
  for (method in methods) {
    print(method)
    df_tmp <- df_lb

    df_tmp$value_lb <- df_lb[, names(df_lb) == method]
    df_tmp$value_ub <- df_ub[, names(df_ub) == method]
    df_tmp$length   <- df_tmp$value_ub - df_tmp$value_lb
    df_tmp$coverage <- df_tmp$value_ub >= df_tmp$true_value & df_tmp$value_lb <= df_tmp$true_value 
    df_tmp$alpha    <- ins$alpha
    y <- aggregate(cbind(coverage, length) ~ alpha + dgp + statistic + n + B, FUN = mean, data = df_tmp, na.rm=TRUE, na.action=NULL)
    y_sd <- aggregate(cbind(coverage, length) ~ alpha + dgp + statistic + n + B, FUN = sd, data = df_tmp, na.rm=TRUE, na.action=NULL)  
    y_n <- aggregate(cbind(coverage, length) ~ alpha + dgp + statistic + n + B, FUN = function(x){sum(!is.na(x))}, data = df_tmp, na.action=NULL)
    names(y_sd) <- paste(names(y_sd), "_sd", sep = "")
    names(y_n) <- paste(names(y_n), "_n", sep = "")
    y <- cbind(y, y_sd[,6:7], y_n[,6:7])
    y$method <- method
    res <- rbind(res, y)
  }
}

saveRDS(rename_for_paper(res), paste0("./data_rds/aggregated_two_sided.rds"))