rm(list = ls(all.names = TRUE))

library(ComplexUpset)
library(ggplot2)

# Read in RF output object with variable importance of features ranked
IBD_varImp <- readRDS("/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/output/IBD_RF_output_varImp.rds")

pathways <- names(IBD_varImp$simpsons_evenness)

# Get a table that lists all the pathways and is 1 if the pathway was in the top 100 and 0 otherwise.
IBD_varImp_top100 <- data.frame(matrix(0, nrow = length(pathways), ncol = length(IBD_varImp)))
colnames(IBD_varImp_top100) <- sort(names(IBD_varImp), decreasing = TRUE)
rownames(IBD_varImp_top100) <- pathways

for (category in names(IBD_varImp)) {

 category_varImp <- IBD_varImp[[category]][pathways]
 # Get descending rank with "rank" command (after multiplying values by -1)
 category_varImp_rank <- rank(category_varImp * -1)
 IBD_varImp_top100[which(category_varImp_rank <= 100), category] <- 1
}

IBD_varImp_top100 <- IBD_varImp_top100[which(rowSums(IBD_varImp_top100) > 0), ]

upset(IBD_varImp_top100,
      intersect = colnames(IBD_varImp_top100),
      min_size = 2,
      set_sizes = FALSE)


# Some striking observations based on this plot are that:
# (1) a lot of the top hits in the relabun/relabun_clr models are not identified by the contributional diversity metrics.\
# (2) There are 4 hits identified by almost all contributional diversity metrics, but not by the relabun/relabun_clr models
# (3) There are 3 hits identified by ALL models as in the top 100 features

# To see the hits present in the top 100 of all models:
rownames(IBD_varImp_top100)[which(rowSums(IBD_varImp_top100) == ncol(IBD_varImp_top100))]

# To see the hits present in almost all contrib div models, but not in the relabun models.
# First get all rows where the pathways are not in the top 100 of the relabun models.
IBD_varImp_top100_nonrelabun <- IBD_varImp_top100[which(IBD_varImp_top100$relabun_clr == 0 & IBD_varImp_top100$relabun == 0), ]

# Then get the subset where simpsons_evenness_clr == 0
IBD_varImp_top100_nonrelabun_nonother <- IBD_varImp_top100_nonrelabun[which(IBD_varImp_top100_nonrelabun$simpsons_evenness_clr == 0), ]

# Then get the subset of pathways that are in the top 100 of the 17 remaining models.
rownames(IBD_varImp_top100_nonrelabun)[which(rowSums(IBD_varImp_top100_nonrelabun) == 17)]

# To look at how these specific features differ between samples you would need to regenerate the object "RF_input" as in the "IBD_run_RF.R" script.
# Then can look at boxplots like this (e.g., to see differences in relabun of specific feature between case and controls)
relabun_PWY.6471_case <- RF_input$relabun$PWY.6471..peptidoglycan.biosynthesis.IV..Enterococcus.faecium[which(RF_input$relabun$disease_state == "IBD")]
relabun_PWY.6471_control <- RF_input$relabun$PWY.6471..peptidoglycan.biosynthesis.IV..Enterococcus.faecium[which(RF_input$relabun$disease_state == "healthy")]
boxplot(relabun_PWY.6471_control, relabun_PWY.6471_case,
        names = c("Healthy", "IBD"),
        ylab = "Rel. Abun.",
        main = "PWY.6471..peptidoglycan.biosynthesis.IV..Enterococcus.faecium")

# Interesting as this particular species is known to be causally involved with IBD!
