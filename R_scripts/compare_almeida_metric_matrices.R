rm(list = ls(all.names = TRUE))

library(MatrixCorrelation)
library(ComplexHeatmap)
library(circlize)

# This code is commented out as it only needed to be run once.
# contrib_div_final <- readRDS("/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics_filled_w_clr.rds")
# 
# psi_values <- data.frame(matrix(NA, nrow = length(contrib_div_final), ncol = length(contrib_div_final)))
# colnames(psi_values) <- names(contrib_div_final)
# rownames(psi_values) <- names(contrib_div_final)
# 
# for (div_metric1 in names(contrib_div_final)) {
#  
#   for (div_metric2 in names(contrib_div_final)) {
#     
#     if (div_metric1 == div_metric2) {
#       psi_values[div_metric1, div_metric2] <- 1
#       next
#     } else {
#       psi_values[div_metric1, div_metric2] <- MatrixCorrelation::PSI(as.matrix(contrib_div_final[[div_metric1]]),
#                                                                      as.matrix(contrib_div_final[[div_metric2]]))
#     }
#   }
# }
# 
# psi_values <- as.matrix(psi_values)

# Save PSI values (run once):
# saveRDS(object = psi_values, file = "/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics_PSI_matrix.rds")


# Then can read these values in when re-running these commands and skip all of the above commands.
psi_values <- readRDS(file = "/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics_PSI_matrix.rds")

cor_break_cols = colorRamp2(c(0, 0.5, 0.95, 1), c("blue3", "white", "red3", "black"))

# Heatmap
heatmap <- draw(Heatmap(psi_values, rect_gp = gpar(type = "none"), column_dend_side = "bottom", col = cor_break_cols,
                        name = "Procrustes\nSimilarity\nIndex",
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        }))

heatmap


# Dendogram:
par(mar=c(5.1, 4.1, 4.1, 12))
plot(row_dend(heatmap), horiz = TRUE)
abline(v = 0.18, lwd = 2, lty = 2)
par(mar=c(5.1, 4.1, 4.1, 2.1))


# Classical Multidimensional Scaling
# N rows (objects) x p columns (variables)
# Each row identified by a unique row name
# Euclidean distances between rows
d <- dist(psi_values)
fit <- cmdscale(d, eig=TRUE, k=2)
#k is the number of dim
#view results
fit
#plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="div_metric", ylab = "variation",
     main= "Metric MDS", type = "n")
text(x, y, labels = row.names(psi_values), cex=.7)


# To figure out the overlapping points (which could be manually noted on a powerpoint sldie for instance):
row.names(psi_values)[which(y > 0.5)]

row.names(psi_values)[which(y > 0.051 & y < 0.07 & x < 0)]


