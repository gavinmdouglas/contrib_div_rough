rm(list = ls(all.names = TRUE))

library(parallel)

source("/home/gdouglas/scripts/contrib_div/R_package/utils.R")

twogroup_wilcox_for_p_and_prop <- function(in_df, group1, group2, focal_func, corr_p = TRUE) {
  
  # Return NA if more than 50% of rows have only NA/NaN/Inf/-Inf
  if (length(which(in_df == NaN)) > 0) {
    in_df[in_df == NaN] <- NA
  }
  
  if (length(which(in_df == Inf)) > 0) {
    in_df[in_df == Inf] <- NA
  }
  
  if (length(which(in_df == -Inf)) > 0) {
    in_df[in_df == -Inf] <- NA
  }  

  if (any(colSums(is.na(in_df)) > 0.5 * nrow(in_df))) {
    return(list(focal_rank=NA, prop_sig=NA))
  }
  
  wilcox_p <- as.numeric()
  
  for (f in rownames(in_df)) {
    wilcox_p <- c(wilcox_p,
                  wilcox.test(as.numeric(in_df[f, group1]),
                              as.numeric(in_df[f, group2]),
                              exact = FALSE)$p.value)
  }
  
  if (corr_p) {
    wilcox_p <- p.adjust(wilcox_p, "BH")
  }
  
  prop_sig <- length(which(wilcox_p < 0.05)) / nrow(in_df)
  
  names(wilcox_p) <- rownames(in_df)
  
  if (wilcox_p[focal_func] < 0.05) {
    focal_rank <- rank(wilcox_p)[focal_func]
  } else {
    focal_rank <- NA
  }
  
  return(list(focal_rank=focal_rank,
              prop_sig=prop_sig))
  
}

wilcox_test_on_all_metrics_and_relabun_relabun_for_prepped_sim <- function(rep_i, info_rds_prefix, metrics_rds_prefix, relabun_rds_prefix) {
  
  RDS_suffix <- ".rds"
  
  info_rds_rep_i <- paste(info_rds_prefix, as.character(rep_i), RDS_suffix, sep = "")
  info_rep_i <- readRDS(info_rds_rep_i)
  
  metrics_rds_rep_i <- paste(metrics_rds_prefix, as.character(rep_i), RDS_suffix, sep = "")
  metrics_rep_i <- readRDS(metrics_rds_rep_i)
  
  relabun_rds_rep_i <- paste(relabun_rds_prefix, as.character(rep_i), RDS_suffix, sep = "")
  relabun_rep_i <- readRDS(relabun_rds_rep_i)
  
  names(metrics_rep_i$clr) <- paste("clr", names(metrics_rep_i$clr), sep = "_")
  
  unique_funcs <- rownames(metrics_rep_i$orig$gini_simpson_index)
  
  wilcox_out <- list()
  wilcox_out[["focal_rank"]] <- as.numeric()
  wilcox_out[["prop_sig"]] <- as.numeric()
  
  relabun_wilcox_out <- twogroup_wilcox_for_p_and_prop(in_df = relabun_rep_i$relabun,
                                                       group1 = info_rep_i$group1,
                                                       group2 = info_rep_i$group2,
                                                       focal_func = info_rep_i$func,
                                                       corr_p = TRUE)
  
  wilcox_out[["focal_rank"]] <- c(wilcox_out[["focal_rank"]], relabun_wilcox_out$focal_rank)
  wilcox_out[["prop_sig"]] <- c(wilcox_out[["prop_sig"]], relabun_wilcox_out$prop_sig)
  
  
  clr_wilcox_out <- twogroup_wilcox_for_p_and_prop(in_df = relabun_rep_i$clr,
                                                   group1 = info_rep_i$group1,
                                                   group2 = info_rep_i$group2,
                                                   focal_func = info_rep_i$func,
                                                   corr_p = TRUE)
  
  wilcox_out[["focal_rank"]] <- c(wilcox_out[["focal_rank"]], clr_wilcox_out$focal_rank)
  wilcox_out[["prop_sig"]] <- c(wilcox_out[["prop_sig"]], clr_wilcox_out$prop_sig)
  
    
  for (m in names(metrics_rep_i$orig)) {
    
    metric_wilcox_out <- twogroup_wilcox_for_p_and_prop(in_df = metrics_rep_i$orig[[m]],
                                                         group1 = info_rep_i$group1,
                                                         group2 = info_rep_i$group2,
                                                         focal_func = info_rep_i$func,
                                                         corr_p = TRUE)
  
    wilcox_out[["focal_rank"]] <- c(wilcox_out[["focal_rank"]], metric_wilcox_out$focal_rank)
    wilcox_out[["prop_sig"]] <- c(wilcox_out[["prop_sig"]], metric_wilcox_out$prop_sig)
    
  }
  
  for (m in names(metrics_rep_i$clr)) {
    
    metric.clr_wilcox_out <- twogroup_wilcox_for_p_and_prop(in_df = metrics_rep_i$clr[[m]],
                                                            group1 = info_rep_i$group1,
                                                            group2 = info_rep_i$group2,
                                                            focal_func = info_rep_i$func,
                                                            corr_p = TRUE)
    
    wilcox_out[["focal_rank"]] <- c(wilcox_out[["focal_rank"]], metric.clr_wilcox_out$focal_rank)
    wilcox_out[["prop_sig"]] <- c(wilcox_out[["prop_sig"]], metric.clr_wilcox_out$prop_sig)
    
  }
  
  names(wilcox_out$focal_rank) <- c("relabun", "relabun.clr", names(metrics_rep_i$orig), names(metrics_rep_i$clr))
  names(wilcox_out$prop_sig) <- c("relabun", "relabun.clr", names(metrics_rep_i$orig), names(metrics_rep_i$clr))
  
  return(wilcox_out)
}


# func_sim_wilcox <- mclapply(X = 1:1000,
#                             FUN = wilcox_test_on_all_metrics_and_relabun_relabun_for_prepped_sim,
#                             info_rds_prefix = "/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep",
#                             metrics_rds_prefix = "/data1/gdouglas/projects/contrib_div/output/sim_validation/func_sim_metrics/rep",
#                             relabun_rds_prefix = "/data1/gdouglas/projects/contrib_div/output/sim_validation/func_sim_relabun/rep",
#                             mc.cores = 50)
# saveRDS(object = func_sim_wilcox, file = "/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/func_sim_wilcox.rds")


taxa_sim_wilcox <- mclapply(X = 1:1000,
                            FUN = wilcox_test_on_all_metrics_and_relabun_relabun_for_prepped_sim,
                            info_rds_prefix = "/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep",
                            metrics_rds_prefix = "/data1/gdouglas/projects/contrib_div/output/sim_validation/taxa_sim_metrics/rep",
                            relabun_rds_prefix = "/data1/gdouglas/projects/contrib_div/output/sim_validation/taxa_sim_relabun/rep",
                            mc.cores = 50)
saveRDS(object = taxa_sim_wilcox, file = "/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/taxa_sim_wilcox.rds")

