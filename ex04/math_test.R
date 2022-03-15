library(ggplot2)
library(dplyr)
library(MCMCpack)

mt <- read.csv("./modeling/SDS383D-master/data/mathtest.csv")
avg_tbl <- mt %>% 
  group_by(school, .drop=FALSE) %>% 
  summarise(counts = n(), mean_score = mean(mathscore))

ggplot(avg_tbl) +
  geom_point(aes(x=counts, y=mean_score))
ggsave("./modeling/SDS383c/ex04/counts_means.png")

P = dim(avg_tbl)[1]
N = sum(avg_tbl$counts)

sample.theta.i <- function(i, mu, tau.sq, sigma.sq) {
  Ni = avg_tbl[i,]$counts
  y_bar = avg_tbl[i,]$mean_score
  v = Ni/sigma.sq + 1/(sigma.sq * tau.sq)
  m = ((mu / (tau.sq*sigma.sq)) + (Ni*y_bar) / sigma.sq) / v
  rnorm(1, mean=m, sd=sqrt(v))
}

# Gibbs sampler
n_iter = 1000
sigma.sq = rep(0, n_iter)
tau.sq = rep(0, n_iter)
mu = rep(50, n_iter)
theta =  matrix(data=NA, nrow=n_iter, ncol=dim(avg_tbl)[1])
kappa = matrix(data=NA, nrow=n_iter, ncol=dim(avg_tbl)[1])

# Initialize values
sigma.sq[1] = 1; tau.sq[1] = 1; theta[1,] = rep(0, dim(avg_tbl)[1])
kappa[1,] = sapply(1:P, function(i) 1 / (1+tau.sq[1]*avg_tbl[i,]$counts))

for (iter in 2:n_iter) {
  
  # sample mu
  mu[iter] = rnorm(1, 
                   mean=mean(theta[iter-1,]), 
                   sd=sqrt(tau.sq[iter-1]*sigma.sq[iter-1]/P)
                   )
  
  # sample sigma.sq
  a = (N + P) / 2
  mean_vec = rep(mu[iter], P)
  dev = theta[iter-1,] - mean_vec
  b = 0.5 * ( t(dev) %*% dev + sum(sapply(1:N, function(i) 
    (mt[i,]$mathscore - theta[iter-1,mt[i,]$school])^2
    )))
  sigma.sq[iter] = rinvgamma(1, a, b)
  
  # sample tau.sq
  a = (P + 1) / 2
  b = 0.5 + t(dev) %*% dev / (2*sigma.sq[iter])
  tau.sq[iter] = rinvgamma(1, a, b)
  
  # sample theta
  theta[iter,] = sapply(1:P, 
                       function(i) sample.theta.i(i, mu[iter], tau.sq[iter], sigma.sq[iter])
                       )
  
  
  # calculate kappa
  kappa[iter,] = sapply(1:P, function(i) 1 / (1+tau.sq[iter]*avg_tbl[i,]$counts))
}

kappa_means = apply(kappa, 2, mean)
avg_tbl$kappa = kappa_means
ggplot(avg_tbl) +
  geom_point(aes(x=counts, y=kappa))
ggsave("./modeling/SDS383c/ex04/counts_kappa.png")

iterations = 1:n_iter
gibbs <- data.frame(iterations, mu, sigma.sq, tau.sq)
theta_df <- data.frame(theta)
theta_df$iteration <- 1:n_iter

ggplot(gibbs) +
  geom_line(aes(x=iterations, y=tau.sq))
ggsave("./modeling/SDS383c/ex04/sigma_trace.png")

ggplot(theta_df) +
  geom_line(aes(x=iteration, y=X2))
            