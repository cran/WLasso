Whitening_Lasso <-
function(X, Y, Sigma, gamma=0.95, maxsteps=2000){
  Sigma <- round(Sigma, 6)
  file <- try(SVD_sigma <- svd(Sigma))
  if (class(file) == "try-error") {
    message("Caught an error during SVD.\n")
    Eigen_Sigma <- eigen(Sigma)
    V_sigma <- Eigen_Sigma$vectors
    lam <- Eigen_Sigma$values
    square_root_sigma <- V_sigma%*%diag(sqrt(lam))%*%solve(V_sigma)
    inv_diag <- ifelse(lam<0.000001, 0, 1/sqrt(lam))
    inv_square_root_Sigma <- V_sigma%*%diag(inv_diag)%*%solve(V_sigma)
  } else {
    U_sigma <- SVD_sigma$u
    D_sigma <- SVD_sigma$d
    square_root_sigma <- U_sigma%*%diag(sqrt(D_sigma))%*%t(U_sigma)
    inv_diag <- ifelse(D_sigma<0.000001, 0, 1/sqrt(D_sigma))
    inv_square_root_Sigma <- U_sigma%*%diag(inv_diag)%*%t(U_sigma)
    inv_square_root_Sigma <- U_sigma%*%diag(1/sqrt(D_sigma))%*%t(U_sigma)
  }

  X_tilde0 <- X%*%inv_square_root_Sigma
  suppressWarnings({
  out0 <-  genlasso::genlasso(Y, X_tilde0, inv_square_root_Sigma, maxsteps=maxsteps)
  })
  p <- ncol(X)
  top_grill <- (1:floor((p/5)))*5
  opt_top <- opt_final_top <- c()
  beta_final_df <- matrix(NA, length(out0$lambda), p)
  mse_final=c()
  for(i in 1:length(out0$lambda)){
    beta_tilde <- out0$beta[,i]
    beta_tilde_sort <- sort(abs(beta_tilde),decreasing=TRUE)
    beta_tilde_top <- sapply(top_grill, top_thresh, sorted_vect = beta_tilde_sort, x=beta_tilde)
    yhat <- as.matrix(X_tilde0%*%beta_tilde_top)
    residuals <- sweep(yhat, 1, Y)
    mse_top <- colMeans(residuals^2)
    ratio_mse <- c()
    for(k in 1:(length(top_grill)-1)){
      ratio_mse[k] <- round(mse_top[k+1]/mse_top[k], 6)
    }
    top_ratio <- min(which(ratio_mse >= gamma))
    if(is.infinite(top_ratio)){
      opt_top[i] <- p
    } else {
      opt_top[i] <- top_grill[top_ratio]
    }
    beta_tilde_opt <- top_thresh(x=beta_tilde, thresh = opt_top[i], sorted_vect = beta_tilde_sort)
    beta_final0 <- inv_square_root_Sigma%*%beta_tilde_opt
    beta_sort <- sort(abs(beta_final0), decreasing=TRUE)
    beta_interm <- sapply(top_grill, top, x=beta_sort, sorted_vect=beta_sort)
    yhat <- as.matrix(X%*%beta_interm)
    residuals <- sweep(yhat, 1, Y)
    mse_final_top <- colMeans(residuals^2)
    ratio_mse <- c()
    for(k in 1:(length(top_grill)-1)){
      ratio_mse[k] <- round(mse_final_top[k+1]/mse_final_top[k], 6)
    }
    top_ratio <- min(which(ratio_mse >= gamma))
    if(is.infinite(top_ratio)){
      opt_top[i] <- p
    } else {
      opt_final_top[i] <- top_grill[top_ratio]
    }
    beta_final_df[i, ] <- round(top(x=beta_final0, thresh = opt_final_top[i], sorted_vect = beta_sort),6)
    mse_final[i] <- mse_final_top[top_ratio]
  }
  beta_min <- beta_final_df[which.min(mse_final), ]
  return(list(lambda = out0$lambda, beta=beta_final_df, trans_mat=inv_square_root_Sigma, beta.min=beta_min, mse=mse_final))
}
