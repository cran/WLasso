Sigma_Estimation <-
function(X){
  matrix_cor=cor(X)
  dist_cor <- dist(matrix_cor)
  #dist_cor <- as.dist(1-abs(matrix_cor))
  p <- ncol(matrix_cor)
  fit_clust <- hclust(dist_cor)
  groups <- cutree(fit_clust, k=2)
  group1 <- which(groups==1)
  group2 <- which(groups==2)
  M1 <- matrix_cor[group1, group1]
  M2 <- matrix_cor[group2,group2]
  M3 <- matrix_cor[-group1, -group2]
  est_alpha1 <- mean(M1 -diag(diag(M1)))
  est_alpha3 <- mean(M2 -diag(diag(M2)))
  est_alpha2 <- mean(M3)
  New_cor <- matrix(0, nrow = p, ncol = p)
  New_cor[group1, group1] <- est_alpha1
  New_cor[group2, group2] <- est_alpha3
  New_cor[-group1, -group2] = New_cor[-group2, -group1] = est_alpha2
  diag(New_cor) <- rep(1,p)
  return(list(mat=New_cor, alpha=c(est_alpha3,est_alpha2,est_alpha1), group_act=group2))
}
