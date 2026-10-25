# Run explicitly from the repository root after installing both packages.
stopifnot(file.exists("DESCRIPTION"), file.exists("companion/scaleCorrection/DESCRIPTION"))
if (!requireNamespace("pkgdown", quietly = TRUE)) stop("The declared documentation dependency pkgdown is needed.")
options(downlit.local_packages = c(softridge = "https://byuzbasi.github.io/softridge/",
  scaleCorrection = "https://byuzbasi.github.io/softridge/scaleCorrection/"))
pkgdown::build_site(pkg = ".", new_process = FALSE, install = FALSE, preview = FALSE)
pkgdown::build_site(pkg = "companion/scaleCorrection", new_process = FALSE,
                    install = FALSE, preview = FALSE)
dir.create("docs/scaleCorrection", showWarnings = FALSE)
paths <- list.files("companion/scaleCorrection/docs", full.names = TRUE, all.files = TRUE, no.. = TRUE)
stopifnot(all(file.copy(paths, "docs/scaleCorrection", recursive = TRUE, overwrite = TRUE)))
file.create("docs/.nojekyll")
cat("Built the shared softridge / scaleCorrection site in docs/.\n")
