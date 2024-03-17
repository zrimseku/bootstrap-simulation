setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(reshape2)

df <- readRDS("./data_rds/intervals_subsample.rds")

x <- df[df$n == 32 & 
          df$alpha == 0.95 & 
          df$dgp == "DGPNorm_0_1" & 
          df$statistic == "std" & 
          df$B == 1000,]

y <- data.frame(exact = x$exact, 
                standard = x$standard, 
                double = x$double,
                true_value = x$true_value)

# correlation
cor(y)

# distance to exact
print(mean(abs(y$exact - y$double)))
print(mean(abs(y$exact - y$standard)))

# coverage
print(mean(y$true_value < y$double))
print(mean(y$true_value < y$standard))

# plot
names(y) <- c("exact", "B-n", "DB")
y <- melt(y)
names(y) <- c("method", "value")

g1 <- ggplot(y, aes(x = value, group = method, lty = method)) + 
    geom_density(adjust = 1.75) + 
    theme_bw() + geom_vline(xintercept = 1.0) + 
    annotate("text", x = 0.95, y = 2.5, label = "(true parameter value)", angle = 90) +
    annotate("text", x = 1.55, y = 1.0, label = "DB") +
    annotate("text", x = 0.91, y = 1.4, label = "B-n") +
    annotate("text", x = 1.32, y = 3.0, label = "exact") + ylab("density") + 
  theme(legend.position = "none") + xlab("endpoint")
ggsave("distance_to_exact_example.pdf", g1, width = 6, height = 4)
