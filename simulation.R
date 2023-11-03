# Reimplementation of the simulation study in
# https://doi.org/10.1186/s12874-022-01831-2

library(tidyverse)

expit <- function(x){
  exp(x)/(1 + exp(x))
}

rbern <- function(prob){
  sapply(prob, function(.x) rbinom(1, size = 1, prob = .x))
}

# number of individuals
n <- 100

# Settings from Table 1
theta <- c(9, -0.3, -0.1, 0.5, rep(0.5, 3), 0.1)
sigma_y2 <- 3
mu <- c(0, 0.1, 0, 0, 0.1)
omega <- c(2, 0.1, 0.1, 0.1, 0.1)
alpha <- c(0, 0.1, 0.1, 0.1, 0, 0.3, 0.3, 0.6, 0.1)

df <- tibble(
  id = 1:n
) %>%
  mutate(u1 = rnorm(n),
         u2 = rnorm(n),
         u3 = rnorm(n),
         x1 = exp(u1/2),
         x2 = exp(u1/(1 + exp(u2))),
         x3 = u1 * u3/25,
         eta = rnorm(n))

# Simulate next time step, conditional on the 
# previous values, and the simulation parameters
simulate <- function(n,
                     v_prev,
                     l_prev,
                     a_prev,
                     y_prev,
                     eta,
                     x,
                     theta,
                     sigma_y2,
                     mu,
                     omega,
                     alpha){
  y_mean <- theta[1] +
    theta[2] * a_prev +
    theta[3] * v_prev +
    theta[4] * l_prev +
    theta[5] * x[, 1] +
    theta[6] * x[, 2] +
    theta[7] * x[, 3] +
    theta[8] * eta
  y <- rnorm(n, mean = y_mean, sd = sqrt(sigma_y2))
  
  l_prob <- expit(mu[1] +
                    mu[2] * v_prev +
                    mu[3] * y +
                    mu[4] * l_prev +
                    mu[5] * a_prev)
  
  l <- rbern(prob = l_prob)
  
  v_prob <- expit(omega[1] + 
                    omega[2] * l +
                    omega[3] * x[,1] +
                    omega[4] * x[,2] +
                    omega[5] * x[,3])
  v <- rbern(prob = v_prob)
  
  a_prob <- expit(alpha[1] +
                    alpha[2] * v_prev +
                    alpha[3] * y +
                    alpha[4] * l +
                    alpha[5] * a_prev +
                    alpha[6] * x[,1] +
                    alpha[7] * x[,2] +
                    alpha[8] * x[,3] +
                    alpha[9] * eta)
  a <- rbern(prob = a_prob)
  
  data.frame(
    y = y,
    l = l,
    v = v,
    a = a
  )
}

# one list component for each time point
t <- list()

# initial values
t[[1]] <- data.frame(
  y = rep(0, n),
  l = rep(0, n),
  v = rep(0, n),
  a = rep(0, n)
)

# forward simulate
for (i in 1:10){
  t[[i + 1]] <- simulate(n = 100,
                     v_prev = t[[i]]$v,
                     l_prev = t[[i]]$l,
                     a_prev = t[[i]]$a,
                     y_prev = t[[i]]$y,
                     eta = df$eta,
                     x = cbind(df$x1, df$x2, df$x3),
                     theta = theta,
                     sigma_y2= sigma_y2,
                     mu = mu,
                     omega = omega,
                     alpha = alpha)
}

# proportion of "missed" visits in the simulation
1 - mean(do.call(rbind, t[2:11])$v)
