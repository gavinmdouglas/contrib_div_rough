rm(list = ls(all.names = TRUE))

library(ranger)

clr_transform_by_col <- function(in_df) {
  # Perform CLR transformation on table assuming samples are columns.
  return(data.frame(apply(in_df, 2, function(x){log(x) - mean(log(x))}), check.names=FALSE))
}

# Build Random Forest models to classify case and control samples for each contributional diversity metric separately.
# Also do this for the standard relative abundance table.
# For relative abundance and certain metric tables also convert to CLR as separate tables to test.

CRC_metadata <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/metadata.rds")$CRC
rownames(CRC_metadata) <- CRC_metadata$sample_id

CRC_metrics <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/computed_metrics/CRC_contrib_div_out.rds")

# Add in unstratified relative abundances in as an additional "metric".
CRC_unstrat <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/unstrat_pathway_tables.rds")$CRC
CRC_metrics[["relabun"]] <- CRC_unstrat[rownames(CRC_metrics$simpsons_evenness), ]

# Fill in all NAs with 0s.
for (div_metric in names(CRC_metrics)) {
  if (length(which(is.na(CRC_metrics[[div_metric]]))) > 0) {
    CRC_metrics[[div_metric]][is.na(CRC_metrics[[div_metric]])] <- 0
  }
}


# For all metrics, also fill in all 0's and NAs with smallest non-zero value
# and then compute CLR. Add "clr" to end of these div metric names as they will be tested in addition
# to the original metrics.

for (div_metric in names(CRC_metrics)) {
  
  working_table <- CRC_metrics[[div_metric]]
  
  smallest_value <- min(working_table[working_table > 0], na.rm = TRUE)
  working_table[working_table == 0] <- smallest_value
  
  CRC_metrics[[paste(div_metric, "clr", sep = "_")]] <- clr_transform_by_col(working_table)
}



# Need to transpose tables and add "disease_state" column.
RF_input <- list()

for (div_metric in names(CRC_metrics)) {
  RF_input[[div_metric]] <- data.frame(t(CRC_metrics[[div_metric]]))
  RF_input[[div_metric]] <- RF_input[[div_metric]][which(rownames(RF_input[[div_metric]]) %in% rownames(CRC_metadata)), ]
  RF_input[[div_metric]]$disease_state <- as.factor(CRC_metadata[rownames(RF_input[[div_metric]]), "disease"])
}


# Run RF on all of these input tables.
RF_output <- list()
for (div_metric in names(CRC_metrics)) {
  RF_output[[div_metric]] <- ranger(formula = disease_state ~ ., importance = "permutation", data = RF_input[[div_metric]], num.trees = 10000)
}

# These are the key model accuracies that we can compare:
RF_output_OOB_accuracy <- 1 - sapply(RF_output, function(x) { x$prediction.error})


# To interpret them we need to know what the accuracy would be by just classifying all samples the most abundant class.
# In this case the expectation is an accuracy close to 0.5 so these are all better than random.
table(RF_input$simpsons_evenness$disease_state)
258/(246+258)


saveRDS(object = RF_output,
        file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output/CRC_RF_output.rds")
