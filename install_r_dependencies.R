required_packages <- c(
  "CARBayes",
  "CARBayesdata",
  "sf",
  "spdep",
  "Rcpp",
  "RcppArmadillo",
  "coda",
  "tictoc"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

versions <- vapply(
  required_packages,
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
print(versions)
