library(ggplot2)
library(mvtnorm)

# Takes in a dataset and a kernel function
# Outputs a function which is the kernel smoother
  # Parameter to this function is a vector of regressors
smoother_fitter <- function(X, Y, kernel, h) {
  smoother <- function(x) {
    # Calculate weights
    weights = sapply(X, function(val) kernel((x - val)/h) /h )
    # Normalize weights
    normalized = weights / sum(weights)
    
    (normalized %*% Y)[1,1]
    #normalized
  }
  smoother
}

# Some kernel functions
gaussian_kernel <- function(x) {
  exp(-x^2/2)/sqrt(2*pi)
}
indicator_kernel <- function(x) {
  abs(x) < 1
}

# Create some noisy data
# f(x) = x^2 sin(x)
x = runif(100, -5, 5)
y = sapply(x, function(val) rnorm(1, val^2 * sin(val), 3))
x = x - mean(x)
y = y - mean(y)

data <- data.frame(x, y)

data %>% 
  ggplot(aes(x=x, y=y)) +
  geom_point(size=0.5, alpha=0.5)

smoother <- smoother_fitter(x, y, gaussian_kernel, 0.5)
x_new = seq(-5, 5, by=0.01)

h_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 5, 10)
pred_mat <- matrix(data=NA, nrow=length(x_new), ncol=length(h_values))
for (i in 1:length(h_values)) {
  smoother <- smoother_fitter(x, y, gaussian_kernel, h_values[i])
  y_pred <- sapply(x_new, smoother)
  pred_mat[,i] <- y_pred
}

pred_df <- data.frame(pred_mat, x_new)

ggplot() +
  geom_point(data=data, aes(x=x, y=y), size=0.5, alpha=0.5) +
  geom_line(data=pred_df, aes(x=x_new, y=X1, color="h = 0.1")) +
  geom_line(data=pred_df, aes(x=x_new, y=X2, color="h = 0.5")) +
  geom_line(data=pred_df, aes(x=x_new, y=X3, color="h = 1")) +
  geom_line(data=pred_df, aes(x=x_new, y=X4, color="h = 5")) +
  geom_line(data=pred_df, aes(x=x_new, y=X5, color="h = 10")) +
  labs(color="Legend") +
  scale_color_manual(values=c("turquoise", "coral", "gold", "tomato", "orchid"))
ggsave("./Desktop/Sp22/modeling/SDS383c/ex06/different_h.png")


# Train and test should be in this format: matrix with colums X, Y
# Returns a list
  # Each item in list is a pair: function and error on test set
train_test_smoother <- function(train, 
                                test, 
                                kernel, 
                                h = c(0.1, 0.5, 1, 5, 10)) {
  fns <- vector(mode = "list", length = length(h))
  errors <- vector(length = length(h))
  for (i in 1:length(h)) {
    # Create function
    smoother <- smoother_fitter(train[,1], train[,2], kernel, h[i])
    # Calculate error on test set
    y_pred <- sapply(test[,1], smoother)
    residual <- test[,2] - y_pred
    fns[[i]] = smoother
    errors[i] = (residual%*%residual)[1,1]
  }
  list("functions" = fns, "errors" = errors)
}

# function to create dataset
create_test_train <- function(n, 
                              test_prop, 
                              fn, 
                              sd = 1, 
                              low = -5, 
                              high = 5) {
  x_train = runif(n * (1 - test_prop), low, high)
  exact_train = fn(x_train)
  noisy_train = sapply(exact_train, function(val) rnorm(1, val, sd))
  train <- cbind(x_train, noisy_train)
  
  x_test <- runif(n * test_prop, low, high)
  exact_test = fn(x_test)
  noisy_test = sapply(exact_test, function(val) rnorm(1, val, sd))
  test <- cbind(x_test, noisy_test)
  
  list("train" = train, "test" = test)
}

f <- function(x) x ^ 2 * sin(x)

dataset <- create_test_train(120, 0.2, f, 1, -5, 5)
train <- dataset$train
test <- dataset$test
results <- train_test_smoother(train, test, gaussian_kernel)

wiggly <- function(x) sin(10 * x)
smooth <- function(x) x^3 / 12

wiggly_noisy <- create_test_train(500, 0.2, wiggly, sd = 2)
smooth_noisy <- create_test_train(500, 0.2, smooth, sd = 2)
wiggly_quiet <- create_test_train(500, 0.2, wiggly, sd = 0.4)
smooth_quiet <- create_test_train(500, 0.2, smooth, sd = 0.4)

h_values = c(seq(0.1, 3, by=0.1), 4, 5, 6, 10)

# Helper functions
find_optimal_h <- function(dataset, kernel, h = c(0.1, 0.5, 1, 5, 10)) {
  results <- train_test_smoother(dataset$train,
                                 dataset$test,
                                 kernel,
                                 h = c(0.1, 0.5, 1, 5, 10))
  opt_index <- which(results$errors == min(results$errors))
  best_fn <- results$functions[[opt_index]]
  best_error <- results$errors[opt_index]
  best_h <- h[opt_index]
  list("fun" = best_fn, "error" = best_error, "h" = best_h)
}

create_line_df <- function(fun, low = -5, high = 5, by = 0.01) {
  x_grid <- seq(low, high, by = by)
  y_vals <- sapply(x_grid, fun)
  data.frame(x_grid, y_vals)
}

# Wiggly and Noisy
results <- find_optimal_h(wiggly_noisy, gaussian_kernel, h_values)
wn_test_df <- data.frame(wiggly_noisy$test)
wn_line_df <- create_line_df(results$fun)
ggplot() +
  geom_point(data=wn_test_df, aes(x=x_test, y=noisy_test), size=0.5, alpha=0.5) +
  geom_line(data=wn_line_df, aes(x=x_grid, y=y_vals))
results$h
ggsave("./Desktop/Sp22/modeling/SDS383c/ex06/wiggly_noisy.png")
# for wiggly and noisy, the best h = 0.1

# Wiggly and not-so-noisy
results2 <- find_optimal_h(wiggly_quiet, gaussian_kernel, h_values)
wq_test_df <- data.frame(wiggly_quiet$test)
wq_line_df <- create_line_df(results2$fun)
ggplot() +
  geom_point(data=wq_test_df, aes(x=x_test, y=noisy_test), size=0.5, alpha=0.5) +
  geom_line(data=wq_line_df, aes(x=x_grid, y=y_vals))
results2$h
ggsave("./Desktop/Sp22/modeling/SDS383c/ex06/wiggly_quiet.png")
# here the best h = 0.1 too

# Smooth and Noisy
results3 <- find_optimal_h(smooth_noisy, gaussian_kernel, h_values)
sn_test_df <- data.frame(smooth_noisy$test)
sn_line_df <- create_line_df(results3$fun)
ggplot() +
  geom_point(data=sn_test_df, aes(x=x_test, y=noisy_test), size=0.5, alpha=0.5) +
  geom_line(data=sn_line_df, aes(x=x_grid, y=y_vals))
results3$h
ggsave("./Desktop/Sp22/modeling/SDS383c/ex06/smooth_noisy.png")
# h = 0.1

# Smooth and not-so-noisy
results4 <- find_optimal_h(smooth_quiet, gaussian_kernel, h_values)
sq_test_df <- data.frame(smooth_quiet$test)
sq_line_df <- create_line_df(results4$fun)
ggplot() +
  geom_point(data=sq_test_df, aes(x=x_test, y=noisy_test), size=0.5, alpha=0.5) +
  geom_line(data=sq_line_df, aes(x=x_grid, y=y_vals))
results4$h
ggsave("./Desktop/Sp22/modeling/SDS383c/ex06/smooth_quiet.png")
# h = 0.1

# Using leave one out
leave_one <- function(dataset, kernel, h) {
  x = dataset[,1]
  y = dataset[,2]
  loocv <- c()
  H = x %*% solve(t(x) %*% x) %*% t(x)
  n = length(y)
  for (val in h) {
    fun <- smoother_fitter(x, y, kernel, val)
    y_pred <- sapply(x, fun)
    error = sum(sapply(1:n, function(i) ((y[i] - y_pred[i]) / (1 - H[i,i]))^2))
    loocv <- append(loocv, error)
  }
  loocv
}

# Process data
wn <- rbind(wiggly_noisy$train, wiggly_noisy$test)
wq <- rbind(wiggly_quiet$train, wiggly_quiet$test)
sn <- rbind(smooth_noisy$train, smooth_noisy$test)
sq <- rbind(smooth_quiet$train, smooth_quiet$test)

wn_errors = leave_one(wn, gaussian_kernel, h_values)
wq_errors = leave_one(wq, gaussian_kernel, h_values)
sn_errors = leave_one(sn, gaussian_kernel, h_values)
sq_errors = leave_one(sq, gaussian_kernel, h_values)

# table of errors 
errors = cbind(wn_errors, wq_errors, sn_errors, sq_errors)
rownames(errors) = h_values
colnames(errors) = c("wiggly-noisy", 
                     "wiggly-quiet", 
                     "smooth-noisy", 
                     "smooth-quiet")
write.csv(data.frame(errors), "./Desktop/Sp22/modeling/SDS383c/ex06/loocv.csv")


# Local polynomial
utilities <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/utilities.csv")
x = utilities$temp
y = utilities$gasbill / utilities$billingdays
H = x %*% solve(t(x) %*% x) %*% t(x)

sj <- function(x_new, h, j) {
  sum = 0
  for (i in 1:n) {
    int = (x_new - X[i]) / h
    sum = sum + gaussian_kernel(int) * (X[i] - x_new)^j
  }
  sum
}

w <- function(x_new, x_old, h) {
  x_diff = (x_new - x_old) / h
  s1 = sj(x_new, h, 1)
  s2 = sj(x_new, h, 2)
  gaussian_kernel(x_diff) * (s2 - (x_old - x_new) * s1)
}

# Gaussian Process

# Create covariance matrix with Matern function (squared exponential)
matern <- function(x1, x2, b, tau1.sq, tau2.sq) {
  distance = dist(rbind(x1, x2))
  ind = as.integer(all(x1==x2))
  tau1.sq*exp(-.5*(distance/b)^2) + tau2.sq*ind
}

# Create x_1,...,x_N
n <- 500
x_values <- sort(runif(n, 0, 1))

gaussian_process <- function(x, 
                             f = matern, # covariance function
                             m = function(x) 0, # mean function
                             b = 1, 
                             tau1.sq = 1, 
                             tau2.sq = 1
                             ) {
  n = length(x)
  covariance = matrix(data=NA, nrow=n, ncol=n)
  for (i1 in 1:n) {
    for (i2 in 1:n) {
      x1 = x[i1]
      x2 = x[i2]
      # fill in covariance matrix
      covariance[i1,i2] = matern(x1, x2, b, tau1.sq, tau2.sq)
    }
  }
  Y = rmvnorm(n=1, mean=sapply(x, m), sigma=covariance)
  Y
}

tau2 = 1e-6
tau1 = seq(0, 2, length.out = 10)
b = seq(0.0001, 1, length.out = 10)
iters = 3

for (b_i in b) {
  df = data.frame(x_values)
  names = c()
  for (t1 in tau1) {
    y = gaussian_process(x_values, b=b_i, tau1.sq = t1, tau2.sq=tau2)
    t1_round = round(t1, 2)
    name = paste("tau_1^2 =", t1_round)
    df[name] = as.vector(y)
    names = append(names, name)
  }
  # Create and save plot
  plot.df = df %>% gather(key = "pars", value = "value", -1)
  p = ggplot(plot.df, aes(x = x_values, y = value)) + 
    geom_line(aes(color = pars)) + 
    theme_classic() + 
    ggtitle(paste0("Gaussian Processes (b = ", b_i, ")")) + xlab("x")
  bi_title = str_replace(toString(b_i), pattern=".", replacement="_")
  plot_title = paste("b", bi_title, ".jpg", sep="")
  ggsave(filename=paste("./Desktop/Sp22/modeling/SDS383c/ex06/", plot_title,sep=""),
         plot=p)
}


matern_52 <- function(x1, x2, b, tau1.sq, tau2.sq) {
  d =  dist(rbind(x1, x2))
  ind = as.integer(all(x1==x2))
  tau1.sq*(1 + sqrt(5)*d/b + 5*d^2/(3*b^2))*exp(-sqrt(5)*d/b)+tau2.sq*ind
}

for (b_i in b) {
  df = data.frame(x_values)
  names = c()
  for (t1 in tau1) {
    y = gaussian_process(x_values, f=matern_52, b=b_i, tau1.sq = t1, tau2.sq=tau2)
    t1_round = round(t1, 2)
    name = paste("tau_1^2 =", t1_round)
    df[name] = as.vector(y)
    names = append(names, name)
  }
  # Create and save plot
  plot.df = df %>% gather(key = "pars", value = "value", -1)
  p = ggplot(plot.df, aes(x = x_values, y = value)) + 
    geom_line(aes(color = pars)) + 
    theme_classic() + 
    ggtitle(paste0("Gaussian Processes (b = ", b_i, ")")) + xlab("x")
  bi_title = str_replace(toString(b_i), pattern=".", replacement="_")
  plot_title = paste("b", bi_title, ".jpg", sep="")
  ggsave(filename=paste("./Desktop/Sp22/modeling/SDS383c/ex06/matern52", plot_title,sep=""),
         plot=p)
}
