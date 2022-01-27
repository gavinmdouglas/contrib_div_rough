rm(list = ls(all.names = TRUE))

library(MatrixCorrelation)
library(ComplexHeatmap)
library(circlize)

contrib_div_final <- readRDS("/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics_filled_w_clr.rds")

psi_values <- data.frame(matrix(NA, nrow = length(contrib_div_final), ncol = length(contrib_div_final)))
colnames(psi_values) <- names(contrib_div_final)
rownames(psi_values) <- names(contrib_div_final)

for (div_metric1 in names(contrib_div_final)) {
 
  for (div_metric2 in names(contrib_div_final)) {
    
    if (div_metric1 == div_metric2) {
      psi_values[div_metric1, div_metric2] <- 1
      next
    } else {
      psi_values[div_metric1, div_metric2] <- MatrixCorrelation::PSI(as.matrix(contrib_div_final[[div_metric1]]),
                                                                     as.matrix(contrib_div_final[[div_metric2]]))
    }
  }
}

psi_values <- as.matrix(psi_values)


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


# Save PSI values:
saveRDS(object = psi_values, file = "/data1/gdouglas/projects/contrib_div/output/Almeida_2019/2022_01_19_top100samples_KO_metrics_PSI_matrix.rds")
