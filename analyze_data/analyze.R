kl <- function(p, alpha) {
  p * log2(p / alpha) + (1 - p) * log2((1 - p) / (1 - alpha))
}

prepare_table <- function(res, methods, FUN = mean, FUN2 = min, decimals = 2) {
  names(res)[1:length(methods)] <- methods
  tbl <- NULL
  
  tbl <- rbind(tbl, apply(res[,1:length(methods)], 2, FUN, na.rm = T)) 

  rownames(tbl)[1] <- "overall"
  for (n in unique(res$n)) {
    tbl <- rbind(tbl, apply(res[res$n == n,][,1:length(methods)], 2, FUN, na.rm = T)) 
    rownames(tbl)[nrow(tbl)] <- n
  }
  
  for (functional in sort(unique(res$functional))) {
    tbl <- rbind(tbl, apply(res[res$functional == functional,]  [,1:length(methods)], 2, FUN, na.rm = T)) 
    rownames(tbl)[nrow(tbl)] <- functional  
  }
  
  tbl <- t(tbl[,order(tbl[1,], decreasing = F)])
  tbl2 <- as.data.frame(format(round(tbl, decimals), nsmall = 2))
  
  for (i in 1:ncol(tbl2)) {
    x <- tbl2[,i]
    new <- c()
    for (j in 1:length(x)) {
      if (decimals == 3) {
        print(x)
        if (x[j] == FUN2(x)) {
          new <- c(new, sprintf("$\\underline{%s}$", x[j]))
        } else {
          new <- c(new, sprintf("$%s$", x[j])) 
        }
      } else {
        if (x[j] == FUN2(x)) {
          new <- c(new, sprintf("$\\underline{%s}$", x[j]))
        } else {
          new <- c(new, sprintf("$%s$", x[j])) 
        }
      }
    }
    tbl2[,i] <- new
  }
  print(xtable(tbl2), sanitize.text.function = identity)
}


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(xtable)

df <- readRDS("aggregated_one_sided_final.rds")
df <- df[df$B == 1000,]

# remove experiments where too many replications failed ------------------------
# NOTE: It's only (DGPLaplace_0_1, percentile_95, 4) across all alpha for B-t.
df[df$failed > 0 & df$failed < 9990,c(6:11)] <- NA

# precompute stuff -------------------------------------------------------------
df$kl <- kl(df$coverage, df$nominal_coverage)

boot_methods     <- c("BB", "BC", "BCa", "B-n", "PB", "SB", "B-t", "DB")
baseline_methods <- setdiff(unique(df$method), c(boot_methods, "exact"))

# ------------------------------------------------------------------------------
# Table 2: Mean KL (overall, by n, by statistic) -------------------------------
# ------------------------------------------------------------------------------
df2 <- df[df$method %in% boot_methods,]
experiments <- unique(df2[,1:4])
res <- NULL
for (i in 1:nrow(experiments)) {
  ex <- experiments[i,]
  tmp <- df2[df2$nominal_coverage == ex$nominal_coverage &
               df2$dgp == ex$dgp &
               df2$functional == ex$functional &
               df2$n == ex$n,]
  x <- tmp$kl[match(boot_methods, tmp$method)] 
  if (sum(is.nan(x) | is.na(x)) == 0) {
    res <- rbind(res, data.frame(t(x), ex))
  } else {
  }
}

prepare_table(res, boot_methods, FUN = mean, decimals = 3, FUN2 = min)


# ------------------------------------------------------------------------------
# Table 3: Percentage of "not good enough KL" (overall, by n, by statistic) ----
# ------------------------------------------------------------------------------

df2 <- df[df$method %in% boot_methods,]
experiments <- unique(df2[,1:4])
res <- NULL
for (i in 1:nrow(experiments)) {
  ex <- experiments[i,]
  tmp <- df2[df2$nominal_coverage == ex$nominal_coverage & 
               df2$dgp == ex$dgp &
               df2$functional == ex$functional &
               df2$n == ex$n,]
  x <- tmp$kl[match(boot_methods, tmp$method)] > kl(0.945,0.95)
  if (sum(is.nan(x) | is.na(x)) == 0) {
    res <- rbind(res, data.frame(t(x), ex))
  } else {
  }
}

prepare_table(res, boot_methods, FUN = mean, decimals = 2)

# ------------------------------------------------------------------------------
# Table 7: Mean on distance from exact (overall, by n, by statistic) -----------
# NOTE: Divided by 2 times standard deviation of exact intervals ---------------
# ------------------------------------------------------------------------------
df2 <- df[df$method %in% boot_methods,]
experiments <- unique(df2[,1:4])
res <- NULL
for (i in 1:nrow(experiments)) {
  ex <- experiments[i,]
  tmp <- df2[df2$nominal_coverage == ex$nominal_coverage & 
               df2$dgp == ex$dgp &
               df2$functional == ex$functional &
               df2$n == ex$n,]
  
  ex_row <- df[df$nominal_coverage == ex$nominal_coverage & 
               df$dgp == ex$dgp &
               df$functional == ex$functional &
               df$n == ex$n & df$method == "exact",]

  x <- tmp$abs_dist_to_exact[match(boot_methods, tmp$method)]
  x <- x / (2.0 * ex_row$value_sd)
  res <- rbind(res, data.frame(t((x)), ex))
}

prepare_table(res, boot_methods, FUN = mean, decimals = 3)

# ------------------------------------------------------------------------------
# -------- FIND METHODS THAT PERFORM BETTER ------------------------------------
# ------------------------------------------------------------------------------
get_better_methods <- function(methods, competitors, df) {
  experiments <- unique(df[,1:4])
  res <- NULL
  for (i in 1:nrow(experiments)) {
    ex <- experiments[i,]
    tmp <- df[df$nominal_coverage == ex$nominal_coverage & df$dgp == ex$dgp & df$functional == ex$functional & df$n == ex$n,]
    tgt_row <- tmp[tmp$method %in% methods,]
    exa_row <- tmp[tmp$method == "exact",]
    tmp <- tmp[tmp$method != "exact" & !tmp$method %in% methods & tmp$method %in% competitors,] # don't compare with exact or themselves
    tgt_kl <- Inf
    tgt_ex <- Inf
    tgt_le <- Inf
    if (sum(!is.nan(tgt_row$kl)) > 0) {   # at least one target method has interval
      tgt_kl <- min(tgt_row$kl, na.rm = T)
      tgt_ex <- min(tgt_row$dist_abs, na.rm = T)
      tgt_le <- min(tgt_row$length, na.rm = T)
    }
    
    superior_in_cover <- ""
    best_kl <- NaN
    criterion <- (5 * tmp$kl) <= tgt_kl  & tgt_kl > kl(0.945, 0.95) * 25
    x <- tmp$method[criterion]
    x <- x[!is.na(x)]
    best_kl <- min(tmp$kl[criterion], na.rm = T)
    superior_in_cover <- paste0(superior_in_cover, x, collapse  = ";")  
    
    superior_in_exact <- ""
    best_ex <- NaN
    criterion <- tgt_ex > tmp$abs_dist_to_exact & tmp$kl <= tgt_kl * 5
    x <- tmp$method[criterion]
    x <- x[!is.na(x)]
    best_ex <- min(tmp$abs_dist_to_exact[criterion], na.rm = T)
    superior_in_exact <- paste0(superior_in_exact, x, collapse  = ";") 
    
    res <- rbind(res, data.frame(ex, 
                                 tgt_kl, 
                                 best_kl,
                                 tgt_ex,
                                 best_ex,
                                 tgt_le,
                                 superior_in_cover,
                                 superior_in_exact))
  }
  return (res)
}

process_string <- function(z) {
  c_string <- ""
  x <- sort(table(unlist(strsplit(z$superior_in_cover, ";", fixed = T))), decreasing = T)
  w <- match(names(x), boot_methods)
  x <- x[order(names(x))]
  for (i in 1:length(x)) {
    c_string <- paste0(c_string, names(x)[i], " (", x[i] ,")")
    if (i != length(x)) {
      c_string <- paste0(c_string, "; ")
    }
  }
  c_string
}

# Table 5 preprocess
resA <- get_better_methods("DB", baseline_methods, df)
resB <- get_better_methods(baseline_methods, "DB", df)

# Table 6 preprocess
#resA <- get_better_methods("B-n", baseline_methods, df)
#resB <- get_better_methods(baseline_methods, "B-n", df)

# Table 8 preprocess
#resA <- get_better_methods("DB", c(baseline_methods, "B-n"), df)
#resB <- get_better_methods(c(baseline_methods, "B-n"), "DB", df)
#resA$superior_in_cover <- resA$superior_in_exact
#resB$superior_in_cover <- resB$superior_in_exact


# Prepare table
tmpA <- resA[resA$superior_in_cover != "",]
tmpB <- resB[resB$superior_in_cover != "",]
res2 <- NULL
for (n in unique(df$n)) {
  for (s in unique(df$functional)) {
    zA <- tmpA[tmpA$n == n & tmpA$functional == s,]
    zB <- tmpB[tmpB$n == n & tmpB$functional == s,] 
    
    c_stringA <- ifelse(nrow(zA) > 0, process_string(zA), "")
    c_stringB <- ifelse(nrow(zB) > 0, process_string(zB), "")

    if (nrow(zA) > 0 | nrow(zB) > 0) {
      res2 <- rbind(res2, data.frame(n = n, s = s, A_is_better = c_stringA, B_is_better = c_stringB))
    }
  }
}

print(xtable(res2), include.rownames = FALSE)

# ------------------------------------------------------------------------------
# TWO SIDED --------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Table 4: ---------------------------------------------------------------------
# ------------------------------------------------------------------------------
df <- readRDS("aggregated_two_sided_final.rds")
df[df$failed > 0 & df$failed < 9990,c(6:11)] <- NA
df <- df[df$B == 1000,]
df$kl <- kl(df$coverage, df$nominal_coverage)

df2 <- df[df$method %in% boot_methods,]

experiments <- unique(df2[,1:4])
res <- NULL
for (i in 1:nrow(experiments)) {
  ex <- experiments[i,]
  tmp <- df2[df2$nominal_coverage == ex$nominal_coverage & 
               df2$dgp == ex$dgp &
               df2$functional == ex$functional &
               df2$n == ex$n,]
  x <- tmp$kl[match(boot_methods, tmp$method)] > kl(0.945,0.95) * 25
  if (sum(is.nan(x) | is.na(x)) == 0) {
    res <- rbind(res, data.frame(t(x), ex))
  } else {
  }
}

prepare_table(res, boot_methods, FUN = mean, decimals = 2, FUN2 = min)