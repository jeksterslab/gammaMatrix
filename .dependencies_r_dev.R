# R dependencies for development
lib <- ".library"
installed <- installed.packages(lib = lib)
pkg <- c(
  "covr",
  "devtools",
  "knitr",
  "lintr",
  "pkgdown",
  "rmarkdown",
  "rprojroot",
  "styler",
  "testthat"
)
for (i in seq_along(pkg)) {
  if (!(pkg[i] %in% installed)) {
    install.packages(
      pkg[i],
      lib = lib
    )
  }
}
