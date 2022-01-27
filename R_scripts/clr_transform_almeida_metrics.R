rm(list = ls(all.names = TRUE))

source("/home/gdouglas/scripts/contrib_div/R_package/alpha_diversity.R")
source("/home/gdouglas/scripts/contrib_div/R_package/utils.R")

contrib_div <- readRDS("/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics.rds")

names(contrib_div) <- names(calc_diversity)

# Fill in all NA's with 0's in one dataset.
contrib_div_final <- contrib_div
for (div_metric in names(contrib_div)) {
  contrib_div_final[[div_metric]][is.na(contrib_div_final[[div_metric]])] <- 0
}
  
contrib_div_filled <- contrib_div

# Fill in all 0's and NAs with smallest non-zero value and then compute CLR
# Add "clr" to end of div metric names.
for (div_metric in names(contrib_div)) {
  
  smallest_value <- min(contrib_div_filled[[div_metric]][contrib_div_filled[[div_metric]] > 0], na.rm = TRUE)
  contrib_div_filled[[div_metric]][contrib_div_filled[[div_metric]] == 0] <- smallest_value
  contrib_div_filled[[div_metric]][is.na(contrib_div_filled[[div_metric]])] <- smallest_value
  
  contrib_div_final[[paste(div_metric, "clr", sep = "_")]] <- clr_transform_by_col(contrib_div_filled[[div_metric]])
}

saveRDS(object = contrib_div_final, "/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics_filled_w_clr.rds")


# Quick check by eye
for (div_metric in names(contrib_div_final)) {
  print(div_metric)
  print(contrib_div_final[[div_metric]][1:5, 1:5])
}