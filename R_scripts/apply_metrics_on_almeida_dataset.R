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


run_div_metric_long <- function(div_metric, func_tab, abun_tab, in_tree=NULL) {
  
  num_comparisons <- nrow(func_tab) * ncol(abun_tab)
  
  out <- data.frame(matrix(NA, nrow = num_comparisons, ncol = 3))
  colnames(out) <- c("func", "sample", "metric")
  
  row_i <- 1
  
  for (func in rownames(func_tab)) {
    
    taxa_encoding <- colnames(func_tab)[which(func_tab[func, ] > 0)]
    
    for (s in colnames(abun_tab)) {
      
      out[row_i, c("func", "sample")] <- c(func, s)
      
      if (div_metric == "faiths_pd") {
        taxa_encoding_in_sample <- taxa_encoding[which(abun_tab[taxa_encoding, s] > 0)]
        
        if (length(taxa_encoding_in_sample) > 0) {
          out[row_i, "metric"] <- calc_diversity[[div_metric]](tips_in_sample = taxa_encoding_in_sample, tree = in_tree)
        }
      } else {
        taxa_encoding_sample_abun <- abun_tab[taxa_encoding, s]
        
        taxa_encoding_sample_abun <- taxa_encoding_sample_abun[which(taxa_encoding_sample_abun > 0)]
        
        if (length(taxa_encoding_sample_abun) > 0) {
          out[row_i, "metric"] <- calc_diversity[[div_metric]](taxa_encoding_sample_abun)
        }
      }
    
      row_i <- row_i + 1
    }
  }
  
  return(out)
}


# List with a table of samples vs genes for each contrib div metric
contrib_div <- mclapply(names(calc_diversity),
                        run_div_metric_wide,
                        func_tab = almeida_ko,
                        abun_tab = almeida_abun[, top_100_samples],
                        in_tree = almeida_tree,
                        mc.cores = 15)


start_time <- Sys.time()
test <- run_div_metric_long(div_metric = "shannon_index",
                              func_tab = almeida_ko[1:100, ],
                              abun_tab = almeida_abun[, top_100_samples],
                              in_tree = almeida_tree)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
test <- run_div_metric_wide(div_metric = "shannon_index",
                            func_tab = almeida_ko[1:1000, ],
                            abun_tab = almeida_abun[, top_100_samples],
                            in_tree = almeida_tree)
end_time <- Sys.time()
end_time - start_time


saveRDS(object = contrib_div,
        file = "/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics.rds")
