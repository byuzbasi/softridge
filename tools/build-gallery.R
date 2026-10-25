# Deterministic illustrations; neither is a reported simulation comparison.
library(softridge)
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
z <- seq(-3, 3, length.out = 241)
png("man/figures/threshold-geometry.png", width = 1500, height = 850, res = 160)
par(mar = c(4.5, 4.5, 3, 1.2), family = "sans", fg = "#29464d",
    col.axis = "#29464d", col.lab = "#29464d")
plot(z, z, type = "l", lty = 2, lwd = 2, col = "#b1c0bd",
     xlab = "Input coefficient", ylab = "Output coefficient",
     main = "Soft thresholding at a fixed threshold of 1")
abline(h = 0, v = 0, col = "#e3e9e6")
lines(z, soft_threshold(z, 1), col = "#167579", lwd = 4)
legend("topleft", c("Soft threshold", "Identity"),
       col = c("#167579", "#b1c0bd"), lty = c(1, 2), lwd = c(4, 2), bty = "n")
dev.off()

set.seed(26)
X <- matrix(rnorm(20 * 30), 20, 30)
beta <- c(1, -1, rep(0, 28))
threshold <- seq(0, 0.3, length.out = 31)
exact <- vapply(threshold, function(g) hr_risk(X, beta, sigma = 1,
                                             k = 2, gamma = g)$risk, numeric(1))
png("man/figures/exact-risk-example.png", width = 1500, height = 850, res = 160)
par(mar = c(4.5, 4.5, 3, 1.2), family = "sans", fg = "#29464d",
    col.axis = "#29464d", col.lab = "#29464d")
plot(threshold, exact, type = "l", lwd = 4, col = "#167579",
     xlab = "Fixed soft threshold", ylab = "Exact projection coefficient risk",
     main = "An analytic risk curve for one fixed design")
abline(h = exact[1], col = "#c18b53", lty = 2, lwd = 2)
legend("topright", c("STR", "Ridge at the same penalty"),
       col = c("#167579", "#c18b53"), lty = c(1, 2), lwd = c(4, 2), bty = "n")
dev.off()
stopifnot(all(is.finite(exact)))
cat("Two analytic demonstration figures created; no Monte Carlo study.\n")
