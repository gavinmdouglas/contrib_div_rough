rm(list = ls(all.names = TRUE))

func_wilcox <- readRDS("/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/func_sim_wilcox.rds")

func_wilcox_categories <- names(func_wilcox[[1]]$focal_rank)

func_wilcox_rank <- data.frame(matrix(NA, nrow = 1000, ncol = length(func_wilcox_categories) + 2))
colnames(func_wilcox_rank) <- c("func", "num_encoding", func_wilcox_categories)

func_wilcox_prop.sig <- data.frame(matrix(NA, nrow = 1000, ncol = length(func_wilcox_categories) + 2))
colnames(func_wilcox_prop.sig) <- c("func", "num_encoding", func_wilcox_categories)

for (rep_i in 1:1000) {
  
  fun_sim_info_RDS <- paste("/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep",
                            as.character(rep_i),
                            ".rds",
                            sep = "")
  fun_sim_info <- readRDS(fun_sim_info_RDS)
   
  func_wilcox_rank[rep_i, "func"] <- fun_sim_info$func
  func_wilcox_rank[rep_i, "num_encoding"] <- length(fun_sim_info$orig_contributors)
  func_wilcox_rank[rep_i, func_wilcox_categories] <- func_wilcox[[rep_i]]$focal_rank
  
  func_wilcox_prop.sig[rep_i, "func"] <- fun_sim_info$func
  func_wilcox_prop.sig[rep_i, "num_encoding"] <- length(fun_sim_info$orig_contributors)
  func_wilcox_prop.sig[rep_i, func_wilcox_categories] <- func_wilcox[[rep_i]]$prop_sig
  
}


func_wilcox_rank <- func_wilcox_rank[, -which(colSums(is.na(func_wilcox_rank)) == 1000)]
func_wilcox_prop.sig <- func_wilcox_prop.sig[, -which(colSums(is.na(func_wilcox_prop.sig)) == 1000)]

write.table(x = func_wilcox_rank, file = "/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/func_sim_metrics_wilcox_rank.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(x = func_wilcox_prop.sig, file = "/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/func_sim_metrics_wilcox_prop.sig.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



taxa_wilcox <- readRDS("/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/taxa_sim_wilcox.rds")

taxa_wilcox <- taxa_wilcox[which(sapply(taxa_wilcox, length) == 2)]

taxa_wilcox_categories <- names(taxa_wilcox[[1]]$focal_rank)

taxa_wilcox_rank <- data.frame(matrix(NA, nrow = length(taxa_wilcox), ncol = length(taxa_wilcox_categories) + 2))
colnames(taxa_wilcox_rank) <- c("func", "num_encoding", taxa_wilcox_categories)

taxa_wilcox_prop.sig <- data.frame(matrix(NA, nrow = length(taxa_wilcox), ncol = length(taxa_wilcox_categories) + 2))
colnames(taxa_wilcox_prop.sig) <- c("func", "num_encoding", taxa_wilcox_categories)

for (rep_i in 1:length(taxa_wilcox)) {
  
  fun_sim_info_RDS <- paste("/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep",
                            as.character(rep_i),
                            ".rds",
                            sep = "")
  fun_sim_info <- readRDS(fun_sim_info_RDS)
  
  taxa_wilcox_prop.sig[rep_i, "func"] <- fun_sim_info$func
  taxa_wilcox_prop.sig[rep_i, "num_encoding"] <- length(fun_sim_info$orig_contributors)
  taxa_wilcox_prop.sig[rep_i, taxa_wilcox_categories] <- taxa_wilcox[[rep_i]]$prop_sig
  
}

taxa_wilcox_prop.sig <- taxa_wilcox_prop.sig[, -which(colSums(is.na(taxa_wilcox_prop.sig)) == nrow(taxa_wilcox_prop.sig))]

write.table(x = taxa_wilcox_prop.sig, file = "/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/taxa_sim_metrics_wilcox_prop.sig.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
