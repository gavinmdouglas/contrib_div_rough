library("phangorn")

almeida_ko <- read.table(file = "functional_analyses/kegg_summary.csv.gz",
           sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)

almeida_abun <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_abun/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- almeida_abun[, -which(colSums(almeida_abun) == 0)]
almeida_abun <- almeida_abun[-which(rowSums(almeida_abun) == 0), ]

almeida_tree <- read.tree("/data1/gdouglas/projects/contrib_div/data/Almeida_2019/phylogenies/raxml_hgr-umgs_phylogeny.nwk")

almeida_tree <- phangorn::midpoint(almeida_tree)

faiths_pd(tips_in_sample = rownames(almeida_abun)[which(almeida_abun$DRR042264 > 0)], tree = almeida_tree)


