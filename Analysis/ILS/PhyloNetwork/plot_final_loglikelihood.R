library(ggplot2)

log_files <- list.files(pattern = "net\\d+.*\\.log")

results <- data.frame(h = integer(), loglik = numeric(), filename = character())

for (filepath in log_files) {
  h <- as.integer(regmatches(filepath, regexpr("(?<=net)\\d+", filepath, perl = TRUE)))
  
  content <- paste(readLines(filepath, warn = FALSE), collapse = "\n")
  
  matches <- regmatches(content, gregexpr("-loglik\\s+[0-9.]+", content))[[1]]
  
  if (length(matches) == 0) {
    cat("警告:", filepath, "中未找到 -loglik 值\n")
    next
  }
  
  vals <- as.numeric(gsub("-loglik\\s+", "", matches))
  best <- min(vals)
  
  cat(sprintf("h=%d: -loglik = %.4f  (%s)\n", h, best, filepath))
  results <- rbind(results, data.frame(h = h, loglik = best, filename = filepath))
}

results <- results[order(results$h), ]

ggplot(results, aes(x = h, y = loglik)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  scale_x_continuous(breaks = results$h) +
  labs(
    x = "Number of hybridizations (h)",
    y = "-loglik",
    title = "PhyloNet network score across reticulations"
  ) +
  theme_bw(base_size = 13)

ggsave("network_score.png", width = 8, height = 5, dpi = 150)
cat("图已保存为 network_score.png\n")