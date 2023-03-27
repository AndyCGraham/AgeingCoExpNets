#Run all rmd files in order to replicate results in paper
rmd_files <- list("StillingPreProcess.Rmd", "ZhaoPreProcess.Rmd", "SierksmaPreProcess.Rmd", "MouseAC_CoExpNets.Rmd","GTEx Pre Process.Rmd",
                  "NugentSCDemyelinationPreProcess.Rmd", 
                  "OneilPreProcess.Rmd", "NugentBulkDemyelinationPreProcess.Rmd",
                    "Single-Cell_CoExpNets.Rmd", "Monocle.Rmd")
lapply(rmd_files, rmarkdown::render)
