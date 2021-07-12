context('tests')

# For resetting
# testthat::snapshot_review()

library(dplyr)
library(magrittr)
library(ggmagnolia)

test_that("basic plot",{

#   armet_estimates = readRDS("~/PhD/deconvolution/ARMET/dev/armet_estimate.rds")
# 	armet_test = readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/ARMET_dev/dev/armet_TCGA_May25_2021_improved_regression_relative/armet_ACC_lv_4_regression.rds")
#   
	p = 
	  magnolia_input |> 
		plot_polar()
	
	vdiffr::expect_doppelganger("base", p)
	

})
