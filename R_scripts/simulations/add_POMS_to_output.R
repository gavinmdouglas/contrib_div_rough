rm(list = ls(all.names = TRUE))

# Add column for POMS data in output files.

# Note that these data from the POMS manuscript were based on simulations with slightly different function tables, which would need to be fixed if this
# comparison is used for a future manuscript.

func_wilcox_rank <- read.table("/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/func_sim_metrics_wilcox_rank.tsv",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(func_wilcox_rank) <- func_wilcox_rank$func

func_wilcox_prop <- read.table("/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/func_sim_metrics_wilcox_prop.sig.tsv",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(func_wilcox_prop) <- func_wilcox_prop$func

taxa_wilcox_prop <- read.table("/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/taxa_sim_metrics_wilcox_prop.sig.tsv",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(taxa_wilcox_prop) <- taxa_wilcox_prop$func




orig_func_summary <- readRDS("/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/POMS_manuscript_sim_outputs/func_rand_summary.rds")
orig_taxa_summary <- readRDS("/data1/gdouglas/projects/contrib_div/data/POMS_simulation_prepped_RDS/POMS_manuscript_sim_outputs/taxa_rand_summary.rds")
rownames(orig_taxa_summary) <- orig_taxa_summary$focal_func
  

identical(orig_func_summary$focal_func, rownames(func_wilcox_prop))

func_wilcox_rank$POMS <- orig_func_summary$POMS_rank_0.05

func_wilcox_prop$POMS <- orig_func_summary$POMS_sig_0.05

taxa_wilcox_prop$POMS <- orig_taxa_summary[rownames(taxa_wilcox_prop), "POMS_sig_0.05"]

write.table(x = func_wilcox_rank, file = "/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/func_sim_metrics_wilcox_rank_w_POMS.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(x = func_wilcox_prop, file = "/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/func_sim_metrics_wilcox_prop.sig_w_POMS.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(x = taxa_wilcox_prop, file = "/data1/gdouglas/projects/contrib_div/output/sim_validation/wilcox_summary/taxa_sim_metrics_wilcox_prop.sig_w_POMS.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

