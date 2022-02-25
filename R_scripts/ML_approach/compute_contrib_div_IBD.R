rm(list = ls(all.names = TRUE))

library("parallel")

source("/home/gdouglas/scripts/contrib_div/R_package/alpha_diversity.R")
source("/home/gdouglas/scripts/contrib_div/R_package/wrappers.R")

strat_in <- readRDS(file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/pathway_tables.rds")$IBD$strat

strat_in$func <- as.character(strat_in$func)
strat_in$taxon <- as.character(strat_in$taxon)
strat_in$sample <- as.character(strat_in$sample)

div_metrics <- read.table("/data1/gdouglas/projects/contrib_div/data/mapfiles/metrics_to_test.txt", stringsAsFactors = FALSE)$V1

div_metrics <- div_metrics[-which(div_metrics %in% c("faiths_pd", "margalefs_richness"))]

IBD_contrib_div_out <- mclapply(X = div_metrics,
                                FUN = run_div_metric_strat_long,
                                strat_tab = strat_in,
                                mc.cores = length(div_metrics))

names(IBD_contrib_div_out) <- div_metrics

saveRDS(object = IBD_contrib_div_out,
        file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/computed_metrics/IBD_contrib_div_out.rds")
