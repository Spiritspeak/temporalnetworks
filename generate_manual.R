setwd(this.path::this.dir())

roxygen2::roxygenize(package.dir = ".", clean = T)
print(Rcpp::compileAttributes(pkgdir = ".", verbose = T))

devtools::document()
devtools::check(args="--as-cran")
devtools::build_manual(path=".")
