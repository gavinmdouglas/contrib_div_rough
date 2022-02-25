rm(list = ls(all.names = TRUE))

library(parallel)

source("/home/gdouglas/scripts/contrib_div/R_package/alpha_diversity.R")
source("/home/gdouglas/scripts/contrib_div/R_package/utils.R")
source("/home/gdouglas/scripts/contrib_div/R_package/wrappers.R")


calc_metric_for_prepped_sim <- function(rep_i, RDS_folder, RDS_prefix, RDS_suffix, out_prefix) {
  
  infile_rep_i <- paste(RDS_folder, "/", RDS_prefix, as.character(rep_i), RDS_suffix, sep = "")
  
  inobject_rep_i <- readRDS(infile_rep_i)
  
  div_metric_out <- div_metric_workflow_wide_nophylo(metrics = metrics_to_test,
                                                     func_tab = func_annot,
                                                     abun_tab = inobject_rep_i$taxa_perturb_abun)

  outfile <- paste(out_prefix, "/rep", as.character(rep_i), ".rds", sep = "")
  
  saveRDS(object = div_metric_out,
          file = outfile)
  
}

metrics_to_test <- sort(read.table("/data1/gdouglas/projects/contrib_div/data/mapfiles/metrics_to_test.txt",
                              header = FALSE, stringsAsFactors = FALSE)$V1)

# Drop Faith's PD as a metric.
metrics_to_test <- metrics_to_test[which(metrics_to_test != "faiths_pd")]

func_annot <- read.table("/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_annot/kegg_summary.tsv.gz",
                         header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

func_sim_prep <- mclapply(X = 1:1000,
                          FUN = calc_metric_for_prepped_sim,
                          RDS_folder = "/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/MAG.based_prepped_func_sim_info_sel1.5/",
                          RDS_prefix = "func_sim_info_rep",
                          RDS_suffix = ".rds",
                          out_prefix = "/data1/gdouglas/projects/contrib_div/output/sim_validation/func_sim_metrics/",
                          mc.cores = 10)

taxa_files_done <- list.files("/data1/gdouglas/projects/contrib_div/output/sim_validation/taxa_sim_metrics/")
taxa_files_done <- gsub("rep", "", taxa_files_done)
taxa_files_done <- as.numeric(gsub(".rds", "", taxa_files_done))

taxa_reps_to_run <- 1:1000
taxa_reps_to_run <- taxa_reps_to_run[which(! taxa_reps_to_run %in% taxa_files_done)]
taxa_sim_prep <- mclapply(X = taxa_reps_to_run,
                          FUN = calc_metric_for_prepped_sim,
                          RDS_folder = "/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/MAG.based_prepped_taxa_sim_info_sel1.5/",
                          RDS_prefix = "taxa_sim_info_rep",
                          RDS_suffix = ".rds",
                          out_prefix = "/data1/gdouglas/projects/contrib_div/output/sim_validation/taxa_sim_metrics/",
                          mc.cores = 10)

