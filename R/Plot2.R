
source("R/sc2col.R")
source("R/sc2col2heap.R")
library(anticlust)
library(ggplot2)
library(dplyr)

set.seed(1)

out_csv <- "data/Plot2.csv"

M <- 2

N_values <- seq(1000, 10000, 1000)

runs <- 10

method_names <- c("2-COL-CC", "2-COL-CC HEAP")

call_method <- function(method_name, data, target_groups) {
  if (method_name == "2-COL-CC") {
    optimal_dispersion_sc2col(data, target_groups)
  } else if (method_name == "2-COL-CC HEAP") {
    optimal_dispersion_sc2col2stage(data, target_groups)
  }  else {
    stop("Unknown method: ", method_name)
  }
}

results_list <- vector("list", length(N_values) * runs * length(method_names))
idx <- 1

total_steps <- length(N_values) * runs * length(method_names)
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
step <- 0

for (N in N_values) {
  for (r in seq_len(runs)) {
    
    data <- matrix(rnorm(N * M), ncol = M)
    target_groups <- c(N/2, N/2)
    
    for (method_name in method_names) {
      
      gc()
      
      t0 <- proc.time()[["elapsed"]]
      invisible(call_method(method_name, data, target_groups))
      t1 <- proc.time()[["elapsed"]]
      
      results_list[[idx]] <- data.frame(
        N = N,
        M = M,
        method = method_name,
        run = r,
        elapsed_sec = as.numeric(t1 - t0),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
      
      step <- step + 1
      setTxtProgressBar(pb, step)
    }
  }
}

close(pb)

results <- bind_rows(results_list)

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(results, out_csv, row.names = FALSE)

## PRINT PLOT

results <- read.csv(out_csv)
summary_results <- results %>%
  group_by(N, method) %>%
  summarise(
    mean_time = mean(elapsed_sec, na.rm = TRUE), 
    sd_time = sd(elapsed_sec, na.rm = TRUE),
    .groups = "drop"
  )

summary_results$method <- factor(
  summary_results$method,
  levels = c("2-COL-CC", "2-COL-CC HEAP"),
  labels = c("2-COL-CC","2-COL-CC HEAP")
)

p <- ggplot(summary_results, aes(x = N, y = mean_time, color = method)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = pmax(mean_time - sd_time, 0),
      ymax = mean_time + sd_time,
      fill = method
    ),
    alpha = 0.15,
    color = NA
  )+
  geom_point(size = 1.5) +
  xlab("N") +
  ylab("Average Runtime (seconds)") +
  theme_minimal(base_size = 12) +
  scale_color_discrete(labels = c(
    SC2COL = "2-COL-CC",
    "2-COL-CC HEAP" = expression("2-COL-CC"["HEAP"])
  )) +
  scale_fill_discrete(labels = c(
    SC2COL = "2-COL-CC",
    "2-COL-CC HEAP" = expression("2-COL-CC"["HEAP"])
  )) +
  scale_x_continuous(breaks = N_values, labels = format(N_values, big.mark=",")) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.title=element_blank())

p


pdf("Plot2.pdf", width = 6, height = 4.5)
p
dev.off()
