## Mixture of Gaussians (test example)

# Clear workspace
rm(list = ls())

# Load required packages
library(tidyverse)
library(readr)

# Sample test data
mu1 <- c(4,8)
sigma1 <-  matrix(c(2,0,0,2), ncol=2)
y1 <- rmvnorm(n=50, mean=mu1, sigma=sigma1)
plot(y1, xlim=c(-2,12), ylim=c(0,12))

mu2 <- c(8,4)
sigma2 <-  matrix(c(3,0,0,3), ncol=2)
y2 <- rmvnorm(n=50, mean=mu2, sigma=sigma2)
points(y2, pch=20)

y <- cbind(y1,y2)

