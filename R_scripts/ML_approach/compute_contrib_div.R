rm(list = ls(all.names = TRUE))

library("parallel")

source("/home/gdouglas/scripts/contrib_div/R_package/alpha_diversity.R")
source("/home/gdouglas/scripts/contrib_div/R_package/wrappers.R")

combined_tables <- readRDS(file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/pathway_tables.rds")

div_metrics <- read.table("/data1/gdouglas/projects/contrib_div/data/mapfiles/metrics_to_test.txt", stringsAsFactors = FALSE)$V1

div_metrics <- div_metrics[-which(div_metrics == "faiths_pd")]

CRC_contrib_div_out <- mclapply(X = div_metrics,
                                FUN = run_div_metric_strat_long,
                                strat_tab = combined_tables$CRC$strat,
                                mc.cores = length(div_metrics))

saveRDS(object = CRC_contrib_div_out,
        file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/computed_metrics/CRC_contrib_div_out.rds")