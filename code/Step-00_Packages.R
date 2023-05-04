Sys.setenv(LANG = "en_US.UTF-8")
source("functions/Functions.R")
packages <- read.table("required_packages.txt")
package.check(packages[, 1])
