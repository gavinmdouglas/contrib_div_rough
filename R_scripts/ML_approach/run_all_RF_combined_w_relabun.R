rm(list = ls(all.names = TRUE))

library(ranger)

clr_transform_by_col <- function(in_df) {
  # Perform CLR transformation on table assuming samples are columns.
  return(data.frame(apply(in_df, 2, function(x){log(x) - mean(log(x))}), check.names=FALSE))
}

# Run RF models as before, but this time also include relative abundance features in the same model to see how this changes the model accuracy.

for (disease in c("CRC", "IBD", "STH")) {
  
  metadata <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/metadata.rds")[[disease]]
  rownames(metadata) <- metadata$sample_id
  
  metrics_file <- paste("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/computed_metrics/",
                        disease, "_contrib_div_out.rds", sep = "")
  metrics <- readRDS(metrics_file)
  
  # Read in unstratified relative abundances, which for this set of models will be added to every model (either as raw relabun or as CLR-transformed)
  unstrat <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/unstrat_pathway_tables.rds")[[disease]]
  unstrat <- unstrat[rownames(metrics$simpsons_evenness), ]
  
  # Fill in all NAs with 0s.
  for (div_metric in names(metrics)) {
    if (length(which(is.na(metrics[[div_metric]]))) > 0) {
      metrics[[div_metric]][is.na(metrics[[div_metric]])] <- 0
    }
  }
  
  
  # For all metrics, also fill in all 0's and NAs with smallest non-zero value
  # and then compute CLR. Add "clr" to end of these div metric names as they will be tested in addition
  # to the original metrics.
  for (div_metric in names(metrics)) {
    
    working_table <- metrics[[div_metric]]
    
    smallest_value <- min(working_table[working_table > 0], na.rm = TRUE)
    working_table[working_table == 0] <- smallest_value
    
    metrics[[paste(div_metric, "clr", sep = "_")]] <- clr_transform_by_col(working_table)
  }
  
  unstrat_working <- unstrat
  unstrat_smallest_value <-  min(unstrat_working[unstrat_working > 0], na.rm = TRUE)
  unstrat_working[unstrat_working == 0] <- unstrat_smallest_value
  unstrat_clr <- clr_transform_by_col(unstrat_working)
  
  rownames(unstrat_clr) <- gsub("^", "relabun.clr; ", rownames(unstrat_clr))
  rownames(unstrat) <- gsub("^", "relabun; ", rownames(unstrat))
  
  # Need to transpose tables and add "disease_state" column.
  # Also add in relevant relabun features too!
  RF_input <- list()
  
  for (div_metric in names(metrics)) {
    RF_input[[div_metric]] <- data.frame(t(metrics[[div_metric]]))
    RF_input[[div_metric]] <- RF_input[[div_metric]][which(rownames(RF_input[[div_metric]]) %in% rownames(metadata)), ]
    
    if (length(grep("_clr", div_metric)) > 0) {
      relabun_tab <- unstrat_clr
    } else {
      relabun_tab <- unstrat
    }
    
    RF_input[[div_metric]] <- cbind(RF_input[[div_metric]],
                                    data.frame(t(relabun_tab)[rownames(RF_input[[div_metric]]), ]))
    RF_input[[div_metric]]$disease_state <- as.factor(metadata[rownames(RF_input[[div_metric]]), "disease"])
  }
  
  
  # Run RF on all of these input tables.
  RF_output <- list()
  for (div_metric in names(metrics)) {
    RF_output[[div_metric]] <- ranger(formula = disease_state ~ .,
                                      importance = "permutation",
                                      data = RF_input[[div_metric]],
                                      num.trees = 10000)
  }
  
  # These are the key model accuracies that we can compare:
  RF_output_OOB_accuracy <- 1 - sapply(RF_output, function(x) { x$prediction.error})
  
  VarImp <- lapply(RF_output, function(x) { x$variable.importance })
   
  names(VarImp) <- names(RF_output)
  
  RF_output_summary <- list(OOB_accuracy = RF_output_OOB_accuracy,
                            VarImp = VarImp)
  
  out_RDS <- paste("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output_combined_w_relabun/",
                   disease, "_RF_output_summary.rds", sep = "")
  saveRDS(object = RF_output_summary, file = out_RDS)
  
}
