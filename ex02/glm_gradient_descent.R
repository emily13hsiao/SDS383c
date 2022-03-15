
dataset <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/wdbc.csv", header=FALSE)
X <- as.matrix(dataset[,3:12])
X <- apply(X, 2, scale)
X <- cbind(rep(1, 569), X)
Y <- as.numeric(dataset$V2 == "M")
n <- length(Y)

# Canonical link for binomial
link <- function(mu) {
  log(mu/(1-mu))
}

invlink <- function(mu) {
  if (exp(mu) == Inf) {
    return(1)
  }
  return(exp(mu) / (1 + exp(mu)))
}

mu <- function(i, beta) {
  invlink(X[i,] %*% beta)
}

b <- function(theta) {
  if (exp(theta) == Inf) {
    return(theta)
  }
  log(1+exp(theta))
}

# Score - gradient of the log likelihood
score <- function(beta) {
  sapply(1:n, function(i) (Y[i]-mu(i, beta))*X[i,] ) |> rowSums()
}

# Negative of the log likelihood
loglik <- function(beta) {
  -sum(sapply(1:n, function(i) Y[i]*(X[i,]%*%beta) - b(X[i,]%*%beta)))
}

converged <- FALSE
iteration <- 1
tol <- 10e-5
beta <- rep(1, 11)
ll <- loglik(beta)
llchain <- c(ll)

while (!converged) {
  prev_beta <- beta
  ll_prev <- ll
  beta <- beta + 0.0001 * score(beta)
  if (norm(as.matrix(beta-prev_beta)) < tol) {
    converged <- TRUE
  }
  llchain <- c(llchain, -loglik(beta))
}

glm_beta <- glm(Y ~ 0 + X, family = binomial())$coefficients

# Plot comparison
png("./Desktop/Sp22/modeling/SDS383c/ex02/gradient_descent.png")
plot(llchain, type = "l", 
     ylab="log likelihood", 
     main="Log likelihood Gradient Descent",
     ylim = c(-200, -50))
abline(h = -loglik(glm_beta), col = "red")
dev.off()


