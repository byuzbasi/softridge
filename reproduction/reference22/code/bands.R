###############################################################################
## bands.R
##
## Summary of Studies 2, 4 and 5 by design family and signal-to-noise ratio.
## For every simulation setting the mean prediction error of each of the
## thirteen methods of the paper is divided by the smallest mean prediction
## error among them in that setting. The table gives, for each method and each
## group of settings, the median and the largest of these ratios.
##
## Usage: Rscript bands.R   (in the directory that holds results/, after
## simulation.R has been run for Studies 2, 4 and 5)
###############################################################################

METH <- c("Ridge", "Lasso", "ENet", "ALasso", "SCAD", "MCP", "BAR",
          "HardRidge", "HR", "HR_SURE", "HR_OLS", "Ridge_S", "HR_S")
LAB <- c("Ridge", "Lasso", "ENet", "ALasso", "SCAD", "MCP", "BAR",
         "TRidge", "HR", "HR-SURE", "HR-OLS", "Ridge-S", "HR-S")

rd <- function(f) readRDS(file.path("results", f))
s2 <- do.call(rbind, lapply(1:8, function(s) {
  d <- rd(sprintf("study2_scenario%d.rds", s))
  data.frame(group = paste0("sparse_", d$snr), setting = paste0("2.", s), d[, c("method", "PE")])
}))
s4 <- do.call(rbind, lapply(1:8, function(s) {
  d <- rd(sprintf("study4_scenario%d.rds", s))
  g <- if (s <= 4) paste0("block_", d$snr[1]) else "departure_5"
  data.frame(group = g, setting = paste0("4.", s), d[, c("method", "PE")])
}))
s5 <- do.call(rbind, lapply(1:6, function(s) {
  d <- rd(sprintf("study5_scenario%d.rds", s))
  g <- if (d$design[1] == "block") "block_0.2" else "sparse_0.2"
  data.frame(group = g, setting = paste0("5.", s), d[, c("method", "PE")])
}))
all <- rbind(s2, s4, s5)
all <- all[all$method %in% METH, ]
m <- aggregate(PE ~ group + setting + method, data = all, FUN = mean)
best <- aggregate(PE ~ setting, data = m, FUN = min)
m$ratio <- m$PE / best$PE[match(m$setting, best$setting)]
write.csv(m, "results/pe_ratios.csv", row.names = FALSE)

groups <- c("sparse_0.2", "sparse_1", "sparse_5", "block_0.2", "block_1", "block_5", "departure_5")
cell <- function(meth, g) {
  r <- m$ratio[m$method == meth & m$group == g]
  sprintf("%.2f (%.2f)", median(r), max(r))
}
nset <- sapply(groups, function(g) length(unique(m$setting[m$group == g])))
print(nset)
rows <- sapply(seq_along(METH), function(i)
  paste0(LAB[i], " & ", paste(sapply(groups, function(g) cell(METH[i], g)), collapse = " & "), "\\\\"))
writeLines(rows, "results/table_bands.tex")
cat(rows, sep = "\n")

## Number of settings in which each method has the smallest mean error.
wins <- table(factor(m$method[m$ratio == 1], levels = METH))
print(wins)
