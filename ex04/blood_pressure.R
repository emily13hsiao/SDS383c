library(dplyr)
library(mrs)
library(MCMCpack)
library(ggplot2)

bp <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/bloodpressure.csv")

# (A) naive t-test

treatment <- bp %>% filter(treatment == 1)
control <- bp %>% filter(treatment == 2)

test_a = t.test(treatment$systolic, control$systolic)
test_a$stderr

# (B) mean t-test

avg_bp <- bp %>% 
  group_by(subject, .drop=FALSE) %>% 
  summarise(measurements = n(), mean_bp = mean(systolic), treatment=first(treatment))

treatment_mean <- avg_bp %>% filter(treatment==1)
control_mean <- avg_bp %>% filter(treatment==2)

test_b = t.test(treatment_mean$mean_bp, control_mean$mean_bp)
test_b$stderr


# (C) hierarchical model gibbs sampler

y <- function(i, j) {
  bp$systolic[which(bp$subject == i)[j]]
}

n.iter = 10000
P = 20

theta = matrix(data=NA, nrow=n.iter, ncol=P)
sigma.sq <- tau.sq <- mu <- beta <- rep(1, n.iter)

beta.sq <- 10^9
beta.mu <- 0

x <- c(rep(0, 10), rep(1, 10))
n.i <- sapply(1:P, function(i) nrow(filter(bp,subject==i)))
N = sum(n.i)

# Initializations
theta[1,] <- sapply(1:P, function(i) mean(filter(bp, subject==i)[,2]))
mu[1] <- mean(theta[1,])

for (iter in 2:n.iter) {
  
  # Sample mu
  mu_mean = (sum(theta[iter-1,]) - beta[iter-1] * sum(x)) / P
  mu_var = tau.sq[iter-1] * sigma.sq[iter-1] / P
  mu[iter] <- rnorm(1, mu_mean, sqrt(mu_var))
  
  # Sample sigma.sq
  a <- (N + P) / 2
  b <- 0
  for (i in 1:P) {
    for (j in 1:n.i[i]) {
      b = b + (y(i, j) - theta[iter-1, i])^2 / 2
    }
  }
  b = b + sum(sapply(1:P, function(i) (theta[iter-1] - mu[iter] - beta[iter-1] * x[i])^2 / (2 * tau.sq[iter-1])))
  sigma.sq[iter] <- rinvgamma(1, a, b)
  
  # Sample tau.sq
  at = (P + 1) / 2
  bt = 1/2 + sum(sapply(1:P, 
                        function(i) (theta[iter-1, i] - mu[iter] - beta[iter-1] * x[i])^2)
                 ) / (2*sigma.sq[iter])
  tau.sq[iter] <- rinvgamma(1, at, bt)
  
  # Sample theta
  new_theta <- rep(0, P)
  for (i in 1:P) {
    variance = 1 / (n.i[i] / sigma.sq[iter] + 1 / (tau.sq[iter] * sigma.sq[iter]))
    mean = variance * (sum(sapply(1:n.i[i], function(j) y(i, j))) / sigma.sq[iter] + 
                        (mu[iter] + beta[iter-1] * x[i]) / (tau.sq[iter] * sigma.sq[iter])
                      )
    new_theta[i] <- rnorm(1, mean, sqrt(variance))
  }
  theta[iter,] = new_theta
  
  # Sample beta 
  variance = 1 / (sum(x) / (tau.sq[iter] * sigma.sq[iter]) + 1 / beta.sq)
  mean = variance * (sum(x * (new_theta - mu[iter])) / (tau.sq[iter] * sigma.sq[iter]) +
                       beta.mu / beta.sq)
  beta[iter] <- rnorm(1, mean, sqrt(variance))
}

beta_df <- data.frame(beta)

ggplot(beta_df) + geom_histogram(aes(x=beta, y=..density..))
ggsave("./Desktop/Sp22/modeling/SDS383c/ex04/beta_posterior.png")

mean(beta) # -7.313822
sd(beta) # 4.713716

# autocorrelation
for (i in 1:P) {
  # filter to just that person 
  subject_df <- data.frame(bp %>% filter(subject == i))
  ggplot(subject_df, aes(date, systolic, group=1)) + geom_line() +
    theme(axis.text.x = element_text(angle = 60, hjust=1))
  title <- paste("./Desktop/Sp22/modeling/SDS383c/ex04/acf", i, ".png", sep="")
  ggsave(title)
}




