rm(list = ls(all.names = TRUE))

# Save only RDS with variable importance scores, to save space.
IBD_RF_output <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output/IBD_RF_output.rds")

IBD_varImp <- list()

for (category in names(IBD_RF_output)) {
  IBD_varImp[[category]] <- IBD_RF_output[[category]]$variable.importance
}

saveRDS(object = IBD_varImp,
        file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output/IBD_RF_output_varImp.rds")


STH_RF_output <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output/STH_RF_output.rds")

STH_varImp <- list()

for (category in names(STH_RF_output)) {
  STH_varImp[[category]] <- STH_RF_output[[category]]$variable.importance
}

saveRDS(object = STH_varImp,
        file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output/STH_RF_output_varImp.rds")



CRC_RF_output <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output/CRC_RF_output.rds")

CRC_varImp <- list()

for (category in names(CRC_RF_output)) {
  CRC_varImp[[category]] <- CRC_RF_output[[category]]$variable.importance
}

saveRDS(object = CRC_varImp,
        file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output/CRC_RF_output_varImp.rds")

