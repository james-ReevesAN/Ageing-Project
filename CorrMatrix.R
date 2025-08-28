plot_corr_matrix <- function(corr_mat, 
                             title = "Correlation Matrix", 
                             fontsize = 0.7, 
                             palette = colorRampPalette(c("#FFA46B", "#F1EAE4", "#5830F7"))(200)) {
  # Validate input
  if (!is.matrix(corr_mat)) stop("Input must be a matrix")
  if (nrow(corr_mat) != ncol(corr_mat)) stop("Matrix must be square")
  if (!all(abs(corr_mat - t(corr_mat)) < 1e-10)) stop("Matrix must be symmetric")
  
  # Plot
  corrplot(corr_mat,
           method = "color",
           col = palette,
           type = "full",
           order = "hclust",
           addCoef.col = NULL,   
           tl.col = "black",
           tl.cex = fontsize,
           cl.cex = fontsize,
           title = title,
           mar = c(0,0,2,0))
}

  plot_corr_matrix(M, title="My Correlation Heatmap", fontsize=0.8)
  