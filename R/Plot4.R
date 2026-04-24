library(callr)
library(ps)
library(anticlust)
library(here)

library(dplyr)
library(ggplot2)
library(readr)

set.seed(1)

source_files <- c(
  normalizePath(here::here("R", "sc2col.R"), winslash = "/"),
  normalizePath(here::here("R", "sc2col2heap.R"), winslash = "/")
)

out_csv <- "data/Plot4.csv"

M <- 2
N_values <- seq(1000, 10000, 1000)
runs <- 2
poll_interval <- 0.01

method_names <- c("2-COL-CC", "2-COL-CC HEAP")


get_rss_mib <- function(pid) {
  h <- tryCatch(ps_handle(pid), error = function(e) NULL)
  if (is.null(h)) return(NA_real_)
  
  mi <- tryCatch(ps_memory_info(h), error = function(e) NULL)
  if (is.null(mi)) return(NA_real_)
  
  rss_bytes <- NA_real_
  if (is.numeric(mi) && !is.null(names(mi)) && "rss" %in% names(mi)) {
    rss_bytes <- as.numeric(mi[["rss"]])
  } else if (is.list(mi) && "rss" %in% names(mi)) {
    rss_bytes <- as.numeric(mi[["rss"]])
  }
  
  if (!is.finite(rss_bytes)) return(NA_real_)
  rss_bytes / 1024^2
}


run_method_and_measure_memory <- function(data, target_groups, method_name, source_files,
                                          poll_interval = 0.01) {
  
  job <- callr::r_bg(
    func = function(data, target_groups, method_name, source_files) {
      for (f in source_files) source(f)
      
      call_method <- function() {
        if (method_name == "2-COL-CC") {
          optimal_dispersion_sc2col(data, target_groups)
        } else if (method_name == "2-COL-CC HEAP") {
          optimal_dispersion_sc2col2stage(data, target_groups)
        } else {
          stop("Unknown method: ", method_name)
        }
      }
      

      invisible(call_method())
      gc()
      

      invisible(call_method())
      TRUE
    },
    args = list(
      data = data,
      target_groups = target_groups,
      method_name = method_name,
      source_files = source_files
    )
  )
  
  pid <- job$get_pid()
  peak <- 0
  
  while (job$is_alive()) {
    rss_mib <- get_rss_mib(pid)
    if (is.finite(rss_mib) && rss_mib > peak) peak <- rss_mib
    Sys.sleep(poll_interval)
  }
  
  job$get_result()
  peak
}


results <- data.frame(
  N = integer(),
  M = integer(),
  method = character(),
  run = integer(),
  peak_ram_mib = numeric(),
  stringsAsFactors = FALSE
)

total_steps <- length(N_values) * runs * length(method_names)
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
step <- 0

for (N in N_values) {
  for (r in seq_len(runs)) {
    
    set.seed(123 + N + r)
    data <- matrix(rnorm(N * M), ncol = M)
    target_groups <- c(N %/% 2, N - (N %/% 2))
    
    for (method_name in method_names) {
      peak_mem <- run_method_and_measure_memory(
        data = data,
        target_groups = target_groups,
        method_name = method_name,
        source_files = source_files,
        poll_interval = poll_interval
      )
      
      results <- rbind(results, data.frame(
        N = N,
        M = M,
        method = method_name,
        run = r,
        peak_ram_mib = peak_mem,
        stringsAsFactors = FALSE
      ))
      
      step <- step + 1
      setTxtProgressBar(pb, step)
    }
  }
}

close(pb)


dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(results, out_csv, row.names = FALSE)

cat("\nMemory benchmark completed.\nSaved as:", out_csv, "\n")
print(results)


results<-read.csv(out_csv)

summary_results <- results %>%
  group_by(N, method) %>%
  summarise(
    mean_mem = mean(peak_ram_mib, na.rm = TRUE),
    sd_mem = sd(peak_ram_mib, na.rm = TRUE),
    .groups = "drop"
  )

summary_results$method <- factor(
  summary_results$method,
  levels = c("2-COL-CC", "2-COL-CC HEAP"),
  labels = c("2-COLCC", "2-COLCC HEAP")
)


p <- ggplot(summary_results, aes(x = N, y = mean_mem, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_ribbon(
    aes(
      ymin = pmax(mean_mem - sd_mem, 0),
      ymax = mean_mem + sd_mem,
      fill = method
    ),
    alpha = 0.15,
    color = NA
  ) +
  labs(
    x = "n",
    y = "Peak Memory (MiB)",
    color = "Solver",
    fill = "Solver"
  ) +
  theme_minimal(base_size = 12) +
  scale_color_discrete(labels = c(
    SC2COL = "2-COL-CC",
    "2-COLCC HEAP" = expression("2-COLCC"["HEAP"])
  )) +
  scale_fill_discrete(labels = c(
    SC2COL = "2-COL-CC",
    "2-COLCC HEAP" = expression("2-COLCC"["HEAP"])
  ))+
  scale_x_continuous(breaks = N_values, labels = format(N_values, big.mark=",")) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.title=element_blank())+
  scale_color_manual(values = palette()[c(3,6)])+
  scale_fill_manual(values=palette()[c(3,6)])

print(p)



pdf("Plot4.pdf", width = 6, height = 4.5)
p
dev.off()
