library(tidyverse)
library(hash)
library(mvtnorm)

polls <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/polls.csv") %>%
  drop_na()
X <- select(polls, c("state", "edu", "age", "female", "black"))
Y <- select(polls, "bush")

# Replace with integers for categorical variables
X$age <- as.integer(as.factor(X$age)) # this is the correct order
X$edu <- as.integer(factor(X$edu, 
                              levels=c("NoHS", "HS", "SomeColl", "Bacc"))
                    )

states <- unique(X$state)

# dictionary for dfs: key = state name, value = filtered df as matrix
state_dfs <- hash()
for (st in states) {
  mat <- X %>% filter(state == st)
  mat$state <- NULL
  state_dfs[[st]] <- as.matrix(mat)
}

# some parameters
n_iter = 5000
p = 4


# How I am going to store values 
# more dictionaries basically for each one
betas <- hash() 
for (state_str in states) {
  betas[[state_str]] = matrix(data=NA, nrow=n_iter, ncol=p)
  betas[[state_str]][1,] = rep(0, p)
}
Z <- hash()
for (state_str in states) {
  Z[[state_str]] = matrix(data=NA, 
                          nrow=n_iter, 
                          ncol=nrow(state_dfs[[state_str]])
                          )
  Z[[state_str]][1,] = rep(0, nrow(state_dfs[[state_str]]))
}
mu <- hash()
for (state_str in states) {
  mu[[state_str]] = matrix(data=NA, nrow=n_iter, ncol=p)
  mu[[state_str]][1,] = rep(0, p)
}
Sigma <- hash()
for (state_str in states) {
  Sigma[[state_str]] = lapply(1:n_iter, matrix, data=NA, nrow=p, ncol=p)
  Sigma[[state_str]][[1]] = diag(p)
}

# Gibbs sampler
for (it in 2:n_iter) {
  
  # Sample betas for each state
  for (st in states) {
    sigma_inv = solve(Sigma[[st]][[it-1]])
    df = state_dfs[[st]]
    V_inv = sigma_inv + t(df) %*% df
    m = V_inv %*% (sigma_inv %*% mu[[st]][it-1,] + t(df) %*% Z[[st]][it-1,])
    betas[[st]][it,] = rmvnorm(n=1, mean=m, sigma=solve(V_inv))
  }
  
  # Sample Z for each state
  
}

for (st in states) {
  
}


