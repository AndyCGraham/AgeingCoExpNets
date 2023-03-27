# AgeingCoExpNets

This repository contains the code to replicate results from the Ageing CoExpNets Paper

Software Versions:

## R version 4.2.2 (2022-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 22.04.1 LTS
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] plyr_1.8.8                  ggpubr_0.5.0                tictoc_1.1                  ensembldb_2.20.2           
##  [5] AnnotationFilter_1.20.0     GenomicFeatures_1.48.4      AnnotationDbi_1.60.0        WGCNA_1.71                 
##  [9] fastcluster_1.2.3           dynamicTreeCut_1.63-1       gprofiler2_0.2.1            AnnotationHub_3.4.0        
## [13] BiocFileCache_2.6.0         dbplyr_2.2.1                patchwork_1.1.2             SeuratObject_4.1.3         
## [17] Seurat_4.3.0                scales_1.2.1                cowplot_1.1.1               Matrix_1.5-3               
## [21] forcats_1.0.0               stringr_1.5.0               purrr_0.3.5                 tidyr_1.2.1                
## [25] ggplot2_3.4.0               tidyverse_1.3.2             DESeq2_1.36.0               SummarizedExperiment_1.28.0
## [29] Biobase_2.58.0              MatrixGenerics_1.10.0       matrixStats_0.63.0          GenomicRanges_1.50.2       
## [33] GenomeInfoDb_1.34.4         IRanges_2.32.0              S4Vectors_0.36.1            BiocGenerics_0.44.0        
## [37] readr_2.1.3                 limma_3.52.4                tibble_3.1.8                dplyr_1.0.10               
## [41] CoExpNets_0.1.0            
## 
## loaded via a namespace (and not attached):
##   [1] rappdirs_0.3.3                rtracklayer_1.58.0            scattermore_0.8               R.methodsS3_1.8.2            
##   [5] bit64_4.0.5                   knitr_1.41                    irlba_2.3.5.1                 DelayedArray_0.24.0          
##   [9] R.utils_2.12.2                data.table_1.14.8             rpart_4.1.19                  KEGGREST_1.38.0              
##  [13] RCurl_1.98-1.10               GEOquery_2.64.2               doParallel_1.0.17             generics_0.1.3               
##  [17] preprocessCore_1.58.0         RSQLite_2.3.0                 RANN_2.6.1                    future_1.31.0                
##  [21] bit_4.0.5                     tzdb_0.3.0                    spatstat.data_3.0-0           xml2_1.3.3                   
##  [25] lubridate_1.9.2               httpuv_1.6.9                  assertthat_0.2.1              gargle_1.2.1                 
##  [29] xfun_0.37                     hms_1.1.2                     jquerylib_0.1.4               evaluate_0.20                
##  [33] promises_1.2.0.1              restfulr_0.0.15               progress_1.2.2                fansi_1.0.4                  
##  [37] readxl_1.4.2                  igraph_1.4.1                  DBI_1.1.3                     geneplotter_1.74.0           
##  [41] htmlwidgets_1.5.4             spatstat.geom_3.0-6           googledrive_2.0.0             ellipsis_0.3.2               
##  [45] crosstalk_1.2.0               backports_1.4.1               annotate_1.76.0               biomaRt_2.52.0               
##  [49] deldir_1.0-6                  vctrs_0.5.1                   ROCR_1.0-11                   abind_1.4-5                  
##  [53] cachem_1.0.7                  withr_2.5.0                   progressr_0.13.0              vroom_1.6.0                  
##  [57] checkmate_2.1.0               sctransform_0.3.5             GenomicAlignments_1.34.0      prettyunits_1.1.1            
##  [61] goftest_1.2-3                 cluster_2.1.4                 lazyeval_0.2.2                crayon_1.5.2                 
##  [65] genefilter_1.78.0             spatstat.explore_3.0-5        labeling_0.4.2                pkgconfig_2.0.3              
##  [69] ProtGenerics_1.28.0           nlme_3.1-160                  nnet_7.3-18                   rlang_1.0.6                  
##  [73] globals_0.16.2                lifecycle_1.0.3               miniUI_0.1.1.1                filelock_1.0.2               
##  [77] modelr_0.1.10                 cellranger_1.1.0              polyclip_1.10-4               lmtest_0.9-40                
##  [81] carData_3.0-5                 zoo_1.8-11                    reprex_2.0.2                  base64enc_0.1-3              
##  [85] ggridges_0.5.4                googlesheets4_1.0.1           rjson_0.2.21                  png_0.1-8                    
##  [89] viridisLite_0.4.1             bitops_1.0-7                  R.oo_1.25.0                   KernSmooth_2.23-20           
##  [93] Biostrings_2.66.0             blob_1.2.3                    parallelly_1.34.0             spatstat.random_3.1-3        
##  [97] rstatix_0.7.1                 jpeg_0.1-10                   ggsignif_0.6.4                memoise_2.0.1                
## [101] magrittr_2.0.3                ica_1.0-3                     zlibbioc_1.44.0               compiler_4.2.2               
## [105] BiocIO_1.8.0                  RColorBrewer_1.1-3            fitdistrplus_1.1-8            Rsamtools_2.14.0             
## [109] cli_3.6.0                     XVector_0.38.0                listenv_0.9.0                 pbapply_1.7-0                
## [113] htmlTable_2.4.1               Formula_1.2-5                 mgcv_1.8-41                   MASS_7.3-58                  
## [117] tidyselect_1.2.0              stringi_1.7.12                highr_0.10                    yaml_2.3.7                   
## [121] locfit_1.5-9.7                latticeExtra_0.6-30           ggrepel_0.9.2                 grid_4.2.2                   
## [125] sass_0.4.4                    tools_4.2.2                   timechange_0.2.0              future.apply_1.10.0          
## [129] parallel_4.2.2                rstudioapi_0.14               foreach_1.5.2                 foreign_0.8-82               
## [133] gridExtra_2.3                 farver_2.1.1                  Rtsne_0.16                    digest_0.6.31                
## [137] BiocManager_1.30.20           shiny_1.7.3                   Rcpp_1.0.10                   car_3.1-1                    
## [141] broom_1.0.1                   BiocVersion_3.15.2            later_1.3.0                   RcppAnnoy_0.0.20             
## [145] httr_1.4.5                    colorspace_2.1-0              rvest_1.0.3                   XML_3.99-0.13                
## [149] fs_1.6.1                      tensor_1.5                    reticulate_1.26               splines_4.2.2                
## [153] uwot_0.1.14                   spatstat.utils_3.0-1          sp_1.6-0                      plotly_4.10.1                
## [157] xtable_1.8-4                  jsonlite_1.8.4                flashClust_1.01-2             R6_2.5.1                     
## [161] Hmisc_4.7-2                   pillar_1.8.1                  htmltools_0.5.4               mime_0.12                    
## [165] glue_1.6.2                    fastmap_1.1.1                 BiocParallel_1.32.4           interactiveDisplayBase_1.34.0
## [169] codetools_0.2-18              utf8_1.2.3                    lattice_0.20-45               bslib_0.4.1                  
## [173] spatstat.sparse_3.0-0         curl_5.0.0                    leiden_0.4.3                  GO.db_3.16.0                 
## [177] interp_1.1-3                  survival_3.4-0                rmarkdown_2.18                munsell_0.5.0                
## [181] GenomeInfoDbData_1.2.9        iterators_1.0.14              impute_1.70.0                 haven_2.5.1                  
## [185] reshape2_1.4.4                gtable_0.3.1
