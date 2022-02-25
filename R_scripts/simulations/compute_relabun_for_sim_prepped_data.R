rm(list = ls(all.names = TRUE))

library(parallel)

source("/home/gdouglas/scripts/contrib_div/R_package/utils.R")

calc_relabun_for_prepped_sim <- function(rep_i, RDS_folder, RDS_prefix, RDS_suffix, out_prefix) {
  
  infile_rep_i <- paste(RDS_folder, "/", RDS_prefix, as.character(rep_i), RDS_suffix, sep = "")
  
  inobject_rep_i <- readRDS(infile_rep_i)
  
  intersecting_taxa <- colnames(func_annot)[which(colnames(func_annot) %in% rownames(inobject_rep_i$taxa_perturb_abun))]
  
  relabun_out <- calc_func_abun(in_abun = inobject_rep_i$taxa_perturb_abun[intersecting_taxa, ],
                                in_func = t(func_annot[, intersecting_taxa]),
                                ncores = 1)
  
  relabun_out <- relabun_out[which(rowSums(relabun_out) > 0), which(colSums(relabun_out) > 0)]
  
  relabun_out_filled <- relabun_out
  for (s in colnames(relabun_out_filled)) { 
    relabun_out_filled[which(relabun_out_filled[, s] == 0), s] <- min(relabun_out_filled[which(relabun_out_filled[, s] > 0), s])
  }

  relabun_out <- data.frame(sweep(x = relabun_out, MARGIN = 2, STATS = colSums(relabun_out), FUN = '/')) * 100
  
  clr_out <- clr_transform_by_col(relabun_out_filled)
  
  outfile <- paste(out_prefix, "/rep", as.character(rep_i), ".rds", sep = "")
  
  saveRDS(object = list(relabun=relabun_out,
                        clr=clr_out),
          file = outfile)

}


func_annot <- read.table("/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_annot/kegg_summary.tsv.gz",
                         header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

func_sim_prep <- mclapply(X = 1:1,
                          FUN = calc_relabun_for_prepped_sim,
                          RDS_folder = "/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/MAG.based_prepped_func_sim_info_sel1.5/",
                          RDS_prefix = "func_sim_info_rep",
                          RDS_suffix = ".rds",
                          out_prefix = "/data1/gdouglas/projects/contrib_div/output/sim_validation/func_sim_relabun/",
                          mc.cores = 5)


taxa_sim_prep <- mclapply(X = 1:1000,
                          FUN = calc_relabun_for_prepped_sim,
                          RDS_folder = "/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/MAG.based_prepped_taxa_sim_info_sel1.5/",
                          RDS_prefix = "taxa_sim_info_rep",
                          RDS_suffix = ".rds",
                          out_prefix = "/data1/gdouglas/projects/contrib_div/output/sim_validation/taxa_sim_relabun/",
                          mc.cores = 5)


