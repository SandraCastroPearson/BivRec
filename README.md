BivRec Package - Installation requires devtools and Rtools35


Recommended commands for R installation from repository:


library(devtools)

 assignInNamespace("version_info", c(devtools:::version_info, list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")
 
find_rtools()

install_github("SandraCastroPearson/BivRec")
