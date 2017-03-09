#--- create package -------------------------------------------------------------

require(devtools)
require(Rcpp)
require(RcppEigen)

pkg.name <- "LMN"
pkg.path <- "C:/Users/Jerome/Documents/R/LMN"
build.path <- "C:/Users/Jerome/Documents/R/build"

compileAttributes(pkgdir = pkg.path)
document(pkg = pkg.path)
install(pkg = pkg.path)

build(pkg = pkg.path, path = build.path)


#--- ONLY EXECUTE THIS ONCE -----------------------------------------------------

#Rcpp.package.skeleton(name = pkg.name, code_files = "lmn-package.R",
#                      example_code = FALSE)

## add RcppEigen
#DESCRIPTION <- read.dcf(file = file.path(pkg.path,
#                          pkg.name, "DESCRIPTION"))
#DESCRIPTION[,"LinkingTo"] <- paste0(DESCRIPTION[,"LinkingTo"],
#                                    ", RcppEigen")
#write.dcf(DESCRIPTION, file = file.path(pkg.path,
#                         pkg.name, "DESCRIPTION"))
#compileAttributes(pkgdir = file.path(pkg.path, pkg.name))

## unit tests
#use_testthat()

#--- install --------------------------------------------------------------------

compileAttributes(pkgdir = pkg.path)
document(pkg = pkg.path)
install(pkg = pkg.path)

build(pkg = pkg.path)

#--- run tests ------------------------------------------------------------------

test(pkg = pkg.path)
