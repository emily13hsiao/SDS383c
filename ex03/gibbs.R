library(dplyr)
library(MASS)
library(ggplot2)
library(bayestestR)

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

# Gibbs sampler
h = 5
n_iter = 10000
K = diag(0.001, p)
betas <- matrix(data=NA, nrow=n_iter, ncol=6)
w <- matrix(data=NA, nrow=n_iter, ncol=1)
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
  beta_var <- solve(t(X) %*% Sigma %*% X + K) / w[iter-1, 1]
  beta_mean <- beta_var %*% (w[iter-1,1] * t(X) %*% Sigma %*% y +
                               w[iter-1,1] * K %*% m)
  
  new_beta <- mvrnorm(1, beta_mean, beta_var)
  
  # resample w
  w_a <- (n+d)/2
  w_b <- 0.5 * (eta + t(y) %*% Sigma %*% y + t(m) %*% K %*% m) -
    0.5 * (
      (t(y) %*% Sigma %*% X + t(m) %*% K) %*% 
        solve( t(X) %*% Sigma %*% X + K ) %*% 
        t( (t(y) %*% Sigma %*% X+ t(m) %*% K) )
      )
  new_w <- rgamma(1, w_a, w_b)
  
  # resample lambdas
  new_lambda <- sapply(1:n,
                       function(i) rgamma(1,
                                          (h+1)/2, 
                                          0.5*new_w*(y[i]-X[i,]%*%new_beta^2) + h/2)
                       )
  
  # keep the values
  betas[iter,] <- new_beta
  w[iter,] <- c(new_w)
  lambdas[iter,] <- new_lambda
}

# code ran too long so just took iter to 3048...rip 
# code is VERY slow :(
betas <- betas[1:3047,]
w <- w[1:3047,]
lambdas <- lambdas[1:3047,]

lambdadf = data.frame(lambdas)
write.csv(lambdadf, "./Desktop/Sp22/modeling/SDS383c/ex03/lambdas.csv")
write.csv(data.frame(betas), "./Desktop/Sp22/modeling/SDS383c/ex03/betas.csv")
write.csv(data.frame(w), "./Desktop/Sp22/modeling/SDS383c/ex03/omega.csv")

posterior_means <- apply(lambdas, 1, mean)

ci(betas[, 1], method = "HDI") # 95% HDI: [0.58, 2.38]

# Compare to lm
lm_result <- lm(y ~ 0 + X)
summary(lm_result)
# Green_rating estimate: 1.422205   std. Error: 0.447624
ci_lm <- confint(lm_result, "Xgreen_rating", level=0.95)
# 95% CI: [0.5447429, 2.299667]

# Plot residuals
beta_est <- apply(betas, 2, mean)
res <- y - X %*% beta_est
resdf <- data.frame(res)
ggplot(resdf) + geom_histogram(aes(x=res, y=..density..))
ggsave("./Desktop/Sp22/modeling/SDS383c/ex03/heterosked_residuals.png")



