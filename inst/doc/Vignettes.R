## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)
library(WLasso)
library(tibble)
library(ggplot2)
set.seed(123456)

## ----generate Sigma-----------------------------------------------------------
p <- 300 # number of variables 
d <- 10 # number of actives
n <- 50 # number of samples
actives <- sample(1:p, d)
nonacts <- c(1:p)[-actives]
Sigma <- matrix(0, p, p)
Sigma[actives, actives] <- 0.5
Sigma[-actives, actives] <- 0.7
Sigma[actives, -actives] <- 0.7
Sigma[-actives, -actives] <- 0.9
diag(Sigma) <- rep(1,p)

## -----------------------------------------------------------------------------
sort(actives)

## ----X------------------------------------------------------------------------
X <- MASS::mvrnorm(n = n, mu=rep(0,p), Sigma, tol = 1e-6, empirical = FALSE)
beta <- rep(0,p)
beta[actives] <- 2
Y <- X%*%beta+rnorm(n,0,1)

## ----est Sigma----------------------------------------------------------------
Sigma_est <- Sigma_Estimation(X) 

## -----------------------------------------------------------------------------
round(Sigma_est$alpha, 2)

## -----------------------------------------------------------------------------
sort(Sigma_est$group_act)

## -----------------------------------------------------------------------------
Sigma_estmat <- Sigma_est$mat

## ----WLasso model, warning = FALSE--------------------------------------------
mod <- Whitening_Lasso(X = X, Y = Y, Sigma = Sigma_estmat, gamma=0.95, maxsteps = 200)

## ----variable selection,fig.width=4,fig.height=3------------------------------
beta_min <- mod$beta.min
df_beta <- data.frame(beta_est=beta_min, Status = ifelse(beta==0, "non-active", "active"))
df_plot <- df_beta[which(beta_min!=0), ]
df_plot$index <- which(beta_min!=0)
ggplot2::ggplot(data=df_plot, mapping=aes(y=beta_est, x=index, color=Status))+geom_point()+
  theme_bw()+ylab("Estimated coefficients")+xlab("Indices of selected variables")

