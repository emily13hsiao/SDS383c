library(tidyverse)
library(MCMCpack)

cheese <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/cheese.csv")

# Process data to have a log of price instead
cheese$logP <- log(cheese$price)
cheese$disp <- as.factor(cheese$disp)
cheese$logQ <- log(cheese$vol)

# Scatter plot with log P on x access and volume on y, with display as color
cheese %>% 
  ggplot(aes(x=logP, y=vol, color=disp)) +
  geom_point(size=0.5, alpha=0.5)
ggsave("./Desktop/Sp22/modeling/SDS383c/ex05/cheese_eda.png")

# Lazy fix to undo the as.factor
cheese <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/cheese.csv")

# Process data to have a log of price instead
cheese$logP <- log(cheese$price)
cheese$disp <- as.factor(cheese$disp)

# Looks like with display is heteroskedastic and with display is homoskedastic
# Just going to use heteroskedastic generally then

# Parameters are beta_i, sigma_i^2, m, s^2
stores <- unique(cheese$store)
n.stores <- length(stores)

# Matrices of the individual stores, put in a list with names as key
storedf = list()
ydf = list()
for (st in stores) {
  filtered_df = cheese %>% filter(store==st)
  mat = matrix(data=NA, nrow=nrow(filtered_df), ncol=4)
  for (i in 1:nrow(filtered_df)) {
    mat[i,] = c(1, filtered_df[i,4], filtered_df[i,5], filtered_df[i,4] * filtered_df[i,5])
  }
  storedf[[st]] = mat
  
  y = as.matrix(filtered_df$logQ, ncol=1)
  ydf[[st]] = y
}

# Initialization
n.iter = 5000

# Set up lists of matrices to store parameters
beta <- list()
sigma.sq <- list()
for (store in stores) {
  beta[[store]] = matrix(data=NA, nrow=n.iter, ncol=4)
  beta[[store]][1,] = rep(1, 4)
  sigma.sq[[store]] = matrix(data=NA, nrow=n.iter, ncol=1)
  sigma.sq[[store]][1,1] = 1
}
s.sq = matrix(data=NA, nrow=n.iter, ncol=1)
s.sq[1,1] = 1
m = matrix(data=NA, nrow=n.iter, ncol=4)
m[1,] = rep(1, 4)
  
# Set up matrix to store m and s^2
m = matrix(data=NA, nrow=n.iter, ncol=4)
m[1,] = rep(0, 4)
s.sq = matrix(data=NA, nrow=n.iter, ncol=1)
s.sq[1,1] = 1

# Sampling functions as helper
sample.beta.i <- function(store) {
  vari = sigma.sq[[store]][iter-1,1]
  X = storedf[[store]]
  y = ydf[[store]]
  s.squared = s.sq[iter-1,1]
  m.iter = m[iter-1,]
  Sigma = solve(t(X) %*% X / vari + diag(1 / s.squared, nrow=4, ncol=4))
  mu.star = Sigma %*% (t(X) %*% y / vari + m.iter / s.squared)
  rmvnorm(1, mean=mu.star, sigma=solve(Sigma))
}
sample.sigma.sq <- function(store) {
  X = storedf[[store]]
  y = ydf[[store]]
  beta.i = beta[[store]][iter,]
  Ni = dim(X)[1]
  a.star = 2 + Ni/2
  residual = y - X %*% beta.i
  b.star = 1/2 + t(residual) %*% residual / 2
  rinvgamma(1, a.star, b.star)
}
sample.m <- function() {
  s.sq.iter = s.sq[iter-1,1]
  avg.store = rep(0, 4)
  for (store in stores) {
    avg.store = avg.store + beta[[store]][iter,]
  }
  rmvnorm(1, avg.store/n.stores, diag(4) / (n.stores * s.sq.iter))
}
sample.s.sq <- function() {
  a.star = (n.stores + 1) / 2
  diffs = 0 
  for (store in stores) {
    diff = beta[[store]][iter,] - m[iter,]
    diffs = diffs + t(diff) %*% diff
  }
  b.star = 1/2 + diffs / 2
  rinvgamma(1, a.star, b.star)
}

# Start doing the MCMC
for (iter in 2:n.iter) {
  # Sample beta for each store
  for (store in stores) {
    beta.store = sample.beta.i(store)
    beta[[store]][iter,] = beta.store
  }
  
  # Sample sigma.sq for each store
  for (store in stores) {
    sigma.sq.store = sample.sigma.sq(store)
    sigma.sq[[store]][iter,1] = sigma.sq.store
  }
  
  # Update m
  m[iter,] = sample.m()
  
  # Update s.sq
  s.sq[iter,1] = sample.s.sq()
}

# Now let's take a look at the betas
# Histogram of the posterior estimates of betas
beta.means = matrix(data=NA, nrow=n.stores, ncol=4)
for (i in 1:n.stores) {
  store = stores[i]
  beta.store = beta[[store]]
  burn = beta.store[1000:4000,]
  avg = apply(burn, 2, mean)
  beta.means[i,] = avg
}
beta.means = data.frame(beta.means)

beta.means %>% 
  ggplot(aes(x = X1, y = ..density..)) +
  geom_histogram(binwidth = 0.5)
ggsave("./Desktop/Sp22/modeling/SDS383c/ex05/alpha0.png")

beta.means %>% 
  ggplot(aes(x = X2, y = ..density..)) +
  geom_histogram(binwidth = 0.5)
ggsave("./Desktop/Sp22/modeling/SDS383c/ex05/alpha1.png")

beta.means %>% 
  ggplot(aes(x = X3, y = ..density..)) +
  geom_histogram(binwidth = 0.5)
ggsave("./Desktop/Sp22/modeling/SDS383c/ex05/beta0.png")

beta.means %>% 
  ggplot(aes(x = X4, y = ..density..)) +
  geom_histogram(binwidth = 0.5)
ggsave("./Desktop/Sp22/modeling/SDS383c/ex05/beta1.png")

# Histogram of posterior mean estimates of sigma.sq
var.means = matrix(data=NA, nrow=n.stores, ncol=1)
for (i in 1:n.stores) {
  store = stores[i]
  sigma.sq.store = sigma.sq[[store]]
  burn = sigma.sq.store[1000:4000,]
  avg = mean(burn)
  var.means[i,1] = avg
}
hist(var.means)

# True y against predicted y
true_y = cheese$logQ
pred_y = vector(length=5555)
for (i in 1:5555) {
  row = cheese[i,]
  store = row[1]$store
  store_index = which(stores == store)
  beta_store = as.matrix(beta.means[store_index,])
  x = as.matrix(c(1, row[4]$disp, row[5]$logP, row[4]$disp * row[5]$logP))
  prediction = beta_store %*% x
  pred_y[i] = prediction[,1]
}
df = data.frame(true_y, pred_y)
ggplot(df) +
  geom_point(aes(x=true_y, y=pred_y))
ggsave("./Desktop/Sp22/modeling/SDS383c/ex05/true_vs_pred.png")
