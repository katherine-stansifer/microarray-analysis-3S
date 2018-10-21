"# microarray-analysis-3S" 

Does the analysis for Fig. S1 of
https://www.ncbi.nlm.nih.gov/pubmed/30149006
"Isolating mitotic and meiotic germ cells from male mice by developmental synchronization, staging, and sorting."

Compares 3 publicly available expression microarray datasets
*GSE54408: WIN18,446/RA synchronized timecourse
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54408
*GSE12769: unsychronized developmental timecourse https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12769 
*GSE926: unsynchronized developmental timecourse
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE926

Note that all 3 use different array sets; all 3 are from the Griswold lab

1) Setup: make subfolders with downloaded raw data
data->GSE54408_RAW->[total,IP]->[.cel.gz files for each array]
data->GSE12769_RAW->[.cel.gz files for each array]
data->GSE926_RAW->MGU75V2_[A,B,C]->[.cel.gz files for each array]

If you have a different directory setup, just change the paths in process_raw_data.R accordingly

2) Run process_raw_data.R: does basic array processing (rma, pma) and saves the workspace as basic_proc.RData

3) Run/knit secondary_process_and_merge.Rmd: loads data from basic_proc.RData, gets gene-level expression values and zscores, runs limma to get fold-changes, merges all datasets together, and saves results. Saves a bunch of intermediate .RData files in the "results/" subdirectory, of which the most important are:
"limma_merge_to0.RData" (fold-changes relative to d0 data points)
"limma_merge_toave.RData" (fold-changes relative to average expression in each dataset)
"zscore_merge_selected.RData" (expression zscores)

