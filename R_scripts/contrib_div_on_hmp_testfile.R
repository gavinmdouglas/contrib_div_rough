rm(list = ls(all.names = TRUE))

# Functions to calculate a few different
# diversity metrics.
richness <- function(x) {
  return(length(which(x > 0)))
}

shannon <- function(x) { 
  x <- x / sum(x)
  return(-1 * sum(x * log(x)))
}

invsimpson <- function(x) {
  x <- x / sum(x)
  return(1 / sum(x ** 2))
}

simpson <- function(x) {
  x <- x / sum(x)
  return(1 - sum(x ** 2))
}

pilou <- function(x) {
  return(shannon(x) / log(length(x)))
}


# A list containing different functions.
# This is done for convenience in the main function below.
calc_diversity <- list()
calc_diversity[["richness"]] <- richness
calc_diversity[["shannon"]] <- shannon
calc_diversity[["invsimpson"]] <- invsimpson
calc_diversity[["simpson"]] <- simpson
calc_diversity[["pilou"]] <- pilou


### Function to calculate contributional diversity metrics for all samples / functions in a given table.
### Output should be a table with each sample and function on a different row with all specified metrics
### included as separate columns.

calc_contrib_div <- function(intab, metrics = c("richness", "shannon", "invsimpson", "simpson", "pilou")) {
  
  # Check that all specified metrics are in possible set.
  if (length(which(! metrics %in% c("richness", "shannon", "invsimpson", "simpson", "pilou"))) > 0) {
    stop("At least one specified metric is not in the possible set.") 
  }
  
  all_samples <- unique(test_input$sample)
  
  # The first thing we need to figure out is what the size of the final dataframe should be.
  # There should be one row for every func in every sample. Make this in advance is a lot faster
  # than actually adding a new row every time. Accordingly, we'll take some time to figure this
  # out beforehand.
  num_row <- 0
  for (s in all_samples) {
    num_row <- num_row + length(unique(intab[which(intab$sample == s), "func"]))
  }

  # Then need to create this dataframe.
  contrib_div <- data.frame(matrix(NA, nrow = num_row, ncol = 2 + length(metrics)))
  colnames(contrib_div) <- c("sample", "func", metrics)
  
  
  # Now will loop through samples again and actually calculate the metrics this time.
  row_i <- 1
  
  for (s in all_samples) {
    
    sample_subset <- intab[which(intab$sample == s), ]
    
    all_sample_funcs <- unique(sample_subset$func)
    
    for (f in all_sample_funcs) {
      
      contrib_div[row_i, c("sample", "func")] <- c(s, f)
      
      sample_func_taxon_relabun <- sample_subset[which(sample_subset$func == f), "taxon_rel_abun"]
      
      for (m in metrics) {
        contrib_div[row_i, m] <- calc_diversity[[m]](sample_func_taxon_relabun)
      }
      
      row_i <- row_i + 1
    }
    
  }
  
  return(contrib_div)
  
}


# Read in input file.
test_input <- read.table("/Users/Gavin/Google_Drive/postdoc/contrib_div/data/testfiles/test_contrib_file/test_contrib_file.tsv.gz",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

test_input_div <- calc_contrib_div(test_input)

saveRDS(object = test_input_div, file = "/Users/Gavin/Google_Drive/postdoc/contrib_div/data/testfiles/test_contrib_file/test_div.rds")
