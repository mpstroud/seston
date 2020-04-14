

remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")


install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

#check the C++ toolchain
pkgbuild::has_build_tools(debug = TRUE)
