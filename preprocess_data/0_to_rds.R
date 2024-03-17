setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(data.table)
df <- fread("./data_raw/intervals_sample.csv", header = TRUE)
saveRDS(as.data.frame(df), "./data_rds/intervals.rds")