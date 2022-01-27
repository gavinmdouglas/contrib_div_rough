rm(list = ls(all.names = TRUE))

library("parallel")
library("phangorn")

source("/home/gdouglas/scripts/contrib_div/R_package/alpha_diversity.R")
source("/home/gdouglas/scripts/contrib_div/R_package/wrappers.R")

almeida_ko <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_annot/kegg_summary.tsv.gz",
                         sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)

almeida_abun <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_abun/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- almeida_abun[, -which(colSums(almeida_abun) == 0)]
almeida_abun <- almeida_abun[-which(rowSums(almeida_abun) == 0), ]

intersecting_taxa <- colnames(almeida_ko)[which(colnames(almeida_ko) %in% rownames(almeida_abun))]

almeida_ko <- almeida_ko[, intersecting_taxa]

almeida_abun <- almeida_abun[intersecting_taxa, ]

almeida_tree <- read.tree("/data1/gdouglas/projects/contrib_div/data/Almeida_2019/phylogenies/raxml_hgr-umgs_phylogeny.nwk")

almeida_tree <- phangorn::midpoint(almeida_tree)

top_100_samples <- names(sort(colSums(almeida_abun), decreasing = TRUE))[1:100]


# List with a table of samples vs genes for each contrib div metric
contrib_div <- mclapply(names(calc_diversity),
                        run_div_metric_wide,
                        func_tab = almeida_ko,
                        abun_tab = almeida_abun[, top_100_samples],
                        in_tree = almeida_tree,
                        mc.cores = 15)

saveRDS(object = contrib_div,
        file = "/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics.rds")
