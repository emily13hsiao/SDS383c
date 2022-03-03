library(dplyr)
library(MASS)

df <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/greenbuildings.csv")
df$int <- 1
n <- nrow(df)
d <- 1
p <- 6
eta <- 1

y <- df$Rent * df$leasing_rate / 100
y <- as.matrix(y)
X <- df[c("green_rating", "City_Market_Rent", "age", "class_a", "class_b", "int")]
X <- as.matrix(X)
Sigma <- diag(n)
K <- diag(0.001, p)
m <- rep(1, p)


# t-dist parameters for beta
degrees <- n + d
Sigma_st <- t(X) %*% X + K
eta_st <- eta + t(y) %*% y + t(m) %*% K %*% m - 
  (t(y) %*% Sigma %*% X + t(m) %*% K) %*%
  Sigma_st %*% t(t(y) %*% Sigma %*% X + t(m) %*% K)

loc <- solve(Sigma_st) %*% (t(X) %*% y + t(K) %*% m)
scale <- solve(degrees * Sigma_st / eta_st[1,1])

# parameters for green, so beta[1]
g_df <- degrees
g_ncp <- loc[1]

# interval 

fit <- lm(y~X)
g_lm <- fit$coefficients[2]

# Gibbs sampler
h = 5
n_iter = 10000
K = diag(0.001, p)
betas <- matrix(data=NA, nrow=n_iter, ncol=6)
w <- matrix(data=NA, nrow=n_iter, ncol=6)
lambdas <- matrix(data=NA, nrow=n_iter, ncol=n)
d <- 1
eta <- 1
h <- 5

# initialize some values 
w[1,] <- rgamma(1, d/2, eta/2)
lambdas[1,] <- rgamma(n, h/2, h/2)
betas[1,] <- mvrnorm(1, m, solve(w[1,]*K))

for (iter in 2:n_iter) {
  
  # resample betas
  Sigma <- diag(lambdas[iter-1,])
  Sigma_st <- t(X) %*% Sigma %*% X + K
  a_st <- t(X) %*% Sigma %*% y + t(K) %*% m
  
  beta_mean <- solve(Sigma_st) %*% a_st
  beta_var <- solve(Sigma_st) / w[iter-1,1]
  
  new_beta <- mvrnorm(1, beta_mean, beta_var)
  
  # resample w
  w_a <- (n+d)/2
  w_b <- 0.5 * (eta + t(y) %*% Sigma %*% y + t(m) %*% K %*% m) -
    0.5 * (
      (t(y)%*%Sigma%*%X+t(m)%*%K) %*% 
        (t(X)%*%Sigma%*%X+K) %*% 
        t((t(y)%*%Sigma%*%X+t(m)%*%K))
      )
  new_w <- rgamma(1, w_a, w_b)
  
  # resample lambdas
  new_lambda <- sapply(1:n,
                       function(i) rgamma(1,(h+1)/2, 0.5*new_w*(y[i]-X[i,]%*%new_beta^2 + h/2))
                       )
  
  # keep the values
  betas[iter,] <- new_beta
  w[iter,] <- c(new_w)
  lambdas[iter,] <- new_lambda
}

