run_div_metric_wide <- function(div_metric, func_tab, abun_tab, in_tree=NULL) {
  
  out <- data.frame(matrix(NA, nrow = nrow(func_tab), ncol = ncol(abun_tab)))
  rownames(out) <- rownames(func_tab)
  colnames(out) <- colnames(abun_tab)
  
  for (func in rownames(func_tab)) {
    
    taxa_encoding <- colnames(func_tab)[which(func_tab[func, ] > 0)]
    
    for (s in colnames(abun_tab)) {
      
      if (div_metric == "faiths_pd") {
        taxa_encoding_in_sample <- taxa_encoding[which(abun_tab[taxa_encoding, s] > 0)]
        
        if (length(taxa_encoding_in_sample) > 0) {
          out[func, s] <- calc_diversity[[div_metric]](tips_in_sample = taxa_encoding_in_sample, tree = in_tree)
        }
      } else {
        taxa_encoding_sample_abun <- abun_tab[taxa_encoding, s]
        
        taxa_encoding_sample_abun <- taxa_encoding_sample_abun[which(taxa_encoding_sample_abun > 0)]
        
        if (length(taxa_encoding_sample_abun) > 0) {
          out[func, s] <- calc_diversity[[div_metric]](taxa_encoding_sample_abun)
        } else if (div_metric == "richness") {
          out[func, s] <- 0
        }
      }
    }
  }
  
  return(out)
}