library(dplyr)
library(MASS)
library(ggplot2)
library(bayestestR)
library(mvtnorm)

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
Lambda <- diag(n)
K <- diag(0.001, p)
m <- rep(1, p)

# t-parameters
nu_st <- n + d # df
Lambda_st <- t(X) %*% Lambda %*% X + K 
mu_st <- solve(Lambda_st) %*% ( t(X)%*% Sigma %*%y + t(K) %*% m ) # location
eta_st <- eta + t(y)%*%Lambda%*%y + t(m)%*%K%*%m - 
  (t(y)%*%Lambda%*%X + 
     t(m)%*%K)%*%solve(t(X)%*%Lambda%*%X + K)%*%t(t(y)%*%Lambda%*%X + t(m)%*%K)
Sigma_st <- solve(nu_st * Lambda_st / eta_st[1,1]) # scale

betas <- rmvt(n = 1000, sigma = Sigma_st, df = nu_st, delta = mu_st)

ci(betas[,1], method = "HDI") # 95% HDI: [0.61, 2.28]

fit <- lm(y ~ 0 + X)
confint(fit) # Xgreen_rating [0.544742925, 2.29966702]

# Plot residuals
beta_est <- apply(betas, 2, mean)
res <- y - X %*% beta_est
resdf <- data.frame(res)
ggplot(resdf) + geom_histogram(aes(x=res, y=..density..))
ggsave("./Desktop/Sp22/modeling/SDS383c/ex03/greenrating_residuals.png")

