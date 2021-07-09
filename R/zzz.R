#' @importFrom utils packageDescription
.onAttach = function(libname, pkgname) {
	version = packageDescription(pkgname, fields = "Version")
	
	msg = paste0("========================================
", pkgname, " version ", version, "
If you use ggmagnolia in published research, please cite:
1) Mangiola et al. ggmagnolia: an R package for modular heatmap production 
  based on tidy principles. JOSS 2020.
2) Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
  genomic data. Bioinformatics 2016.
This message can be suppressed by:
  suppressPackageStartupMessages(library(ggmagnolia))
========================================
")	
	
	packageStartupMessage(msg)
}

rv = R.Version()

if(getRversion() >= "4.0.0" && as.numeric(rv$`svn rev`) >= 77889) {
	unitType = get("unitType", envir = asNamespace("grid"))
} else {
	unitType = function(x, recurse = TRUE) attr(x, "unit")
}