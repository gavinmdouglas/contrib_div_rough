library(Rcpp)
library(reshape2)

frame_files <- lapply(sys.frames(), function(x) x$ofile)
frame_files <- Filter(Negate(is.null), frame_files)
PATH <- dirname(frame_files[[length(frame_files)]])

Rcpp_filepath <- paste(PATH, "prep_input_list.cpp", sep = "/")
Rcpp::sourceCpp(Rcpp_filepath)

div_metric_workflow_wide_nophylo <- function(metrics, func_tab, abun_tab) {
  
  abun_tab <- abun_tab[which(rowSums(abun_tab) > 0), ]
  abun_tab <- abun_tab[, which(colSums(abun_tab) > 0)]
  
  func_tab <- func_tab[, rownames(abun_tab), drop = FALSE]
  func_tab <- func_tab[which(rowSums(func_tab) > 0), ]
  func_tab <- func_tab[, which(colSums(func_tab) > 0)]
  
  abun_tab <- abun_tab[colnames(func_tab), ]
  
  # Note that "PATH" will only be defined if this script itself was sourced.
  # Source this file each time the function is called otherwise can get odd behaviour (sometimes).
  # But you also get odd behaviour if this is run separately, so perhaps should be commented out.
  #library(Rcpp)
  #Rcpp_filepath <- paste(PATH, "prep_input_list.cpp", sep = "/")
  #Rcpp::sourceCpp(Rcpp_filepath)
  
  prepped_abun <- prep_input_list(abun_tab = as.matrix(abun_tab),
                                  func_tab = as.matrix(func_tab))
  
  func_set <- as.character()
  sample_set <- as.character()
  for (i in 1:nrow(func_tab)) {
    func_set <- c(func_set, rep(rownames(func_tab)[i], ncol(abun_tab)))
    sample_set <- c(sample_set, colnames(abun_tab))
  }
  
  div_metric_out <- list()
  
  div_metric_out_clr <- list()
  
  for (m in metrics) {
    
    div_metric_out[[m]] <- data.frame(matrix(NA, nrow = length(prepped_abun), ncol = 3))
    colnames(div_metric_out[[m]]) <- c("sample", "func", "value")
    div_metric_out[[m]]$sample <- sample_set
    div_metric_out[[m]]$func <- func_set
    
    div_metric_out[[m]]$value <- sapply(prepped_abun, calc_diversity[[m]])
    
    div_metric_out[[m]] <- reshape2::dcast(data = div_metric_out[[m]],
                                           formula = func ~ sample)
    rownames(div_metric_out[[m]]) <- div_metric_out[[m]]$func
    div_metric_out[[m]] <- div_metric_out[[m]][, -1]
    
    div_metric_out_clr[[m]] <- div_metric_out[[m]]
    
    for (s in colnames(div_metric_out_clr[[m]])) {
      
      missing_or_zero_i <- which(is.na(div_metric_out_clr[[m]][, s]) | div_metric_out_clr[[m]][, s] == 0)
      
      if (length(missing_or_zero_i) > 0) {
        div_metric_out_clr[[m]][missing_or_zero_i, s] <- min(div_metric_out_clr[[m]][-missing_or_zero_i, s])
      }
    }
    
    div_metric_out_clr[[m]] <- clr_transform_by_col(div_metric_out_clr[[m]])
    
  }
  
  return(list(orig=div_metric_out,
              clr=div_metric_out_clr))
}




run_div_metric_strat_long <- function(div_metric, strat_tab, in_tree=NULL) {
  
  unique_funcs <- unique(strat_tab$func)
  unique_samples <- unique(strat_tab$sample)
  
  out <- data.frame(matrix(NA, nrow = length(unique_funcs), ncol = length(unique_samples)))
  rownames(out) <- unique_funcs
  colnames(out) <- unique_samples
  
  for (func in unique_funcs) {
    
    strat_tab_func <- strat_tab[which(strat_tab$func == func), ]
    
    for (s in unique_samples) {
      
      strat_tab_func_sample <- strat_tab_func[which(strat_tab_func$sample == s), ]
      
      if (nrow(strat_tab_func_sample) == 0) { next }
      
      if (div_metric == "faiths_pd") {
        
        out[func, s] <- calc_diversity[[div_metric]](tips_in_sample = strat_tab_func_sample$taxon, tree = in_tree)
        
      } else {
        out[func, s] <- calc_diversity[[div_metric]](strat_tab_func_sample$relabun)
        
      }
    }
  }
  
  return(out)
}
