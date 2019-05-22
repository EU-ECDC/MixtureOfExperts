## Mixture of Gaussians (test example)

####################
## Initialisation ##
####################

# Clear workspace
rm(list=ls())

# Load required packages
library(tidyverse)
library(readr)
library(mvtnorm)
library(pdist)

# Sample test data
mu1 <- c(1,4)
sigma1 <-  matrix(c(1,0,0,1), ncol=2)
y1 <- rmvnorm(n=1000, mean=mu1, sigma=sigma1)
colnames(y1) <- c("x", "y")
y1 <- as_tibble(y1)
plot(y1, xlim=c(-2,12), ylim=c(0,12))

mu2 <- c(8,1)
sigma2 <-  matrix(c(1,0,0,1), ncol=2)
y2 <- rmvnorm(n=1000, mean=mu2, sigma=sigma2)
colnames(y2) <- c("x", "y")
y2 <- as_tibble(y2)
points(y2, pch=20)

mu3 <- c(8,8)
sigma3 <-  matrix(c(1,0,0,1), ncol=2)
y3 <- rmvnorm(n=1000, mean=mu3, sigma=sigma3)
colnames(y3) <- c("x", "y")
y3 <- as_tibble(y3)
points(y3, pch=20, col="blue")

sampleData <- bind_rows(y1,y2, y3)
nPoints <- nrow(sampleData)

# Initialise no.of clusters (k)
nClusters <- 3

# Initialise estimated paramaters of Gaussian 
mu.est <- runif((2*nClusters),0,10) # uniform sample
dim(mu.est) <- c(nClusters,2)
colnames(mu.est) <- c("x", "y")

sigma.est <-  rep(matrix(c(1,0,0,1), ncol=2),nClusters)
dim(sigma.est) <- c(2,2,nClusters)
 
# Initialise with uniform probabilities 
weight <- rep(1/nClusters, (nPoints * nClusters)) 
dim(weight) <- c(nPoints, nClusters)
colnames(weight) <- paste0("Cluster", 1:nClusters)
weight <- as_tibble(weight)

# Initialise cluster attribution weights (phi)
phi <- colSums(weight) / nrow(weight)

###################
## Expectation - ##
## maximisation  ##
###################

## Mixture model
tol <- 1e3 # tolerance level
names(tol) <- "tolerance"
i.iteration <- 0 # iterator

density.temp <- rep(NA, (nPoints * nClusters)) 
dim(density.temp) <- c(nPoints, nClusters)
colnames(density.temp) <- paste0("Cluster", 1:nClusters)
density.temp <- as_tibble(density.temp)

sumProduct = function(a,b){
	sumProduct = colSums(a * b) 
	return(sumProduct)
}

# Iterate E-M to convergence
while(i.iteration[length(i.iteration)] < 8){

# (Re-)initialise temporary covariance matrix
	sigma.temp <-  rep(matrix(rep(NA,4), ncol=2),nClusters)
	dim(sigma.temp) <- c(2,2,nClusters)
 
# Expectation step
# Calculate posterior probability of cluster given y_i and other parameters

	for(j in 1:nClusters){
		density.temp[,j] <- dmvnorm(sampleData, mean=mu.est[j,], sigma=sigma.est[,,j]) * phi[j]
	}
	
	# Calculate posterior probability (w_ij) by applying Bayes' rule
	post.prob <- density.temp / rowSums(density.temp)
	
# Maximisation step
	
	# Update estimate for cluster probability
	phi <- colMeans(post.prob)
	
	# Update estimate for Gaussian means
	weightedSum <- apply(post.prob, 2, sumProduct, b=sampleData)
    mu.temp <- t(weightedSum) / colSums(post.prob)
	
	# Update estimate for Gaussian variance
	for(j in 1:nClusters){
		diff1 <- abs(sampleData - mu.est[j,])
		#N.B. Check calculation of covariance 
		sigma.temp[,,j] <- diag(colSums(post.prob[,j] * diff1 * t(diff1)) / sum(post.prob[,j]))
	}
	
	
	# Tolerance is sum of distances between old and new mu and sigma matrices
	tol.temp <- sum(as.matrix(pdist(mu.est, mu.temp))) + sum(as.matrix(pdist(sigma.est[,,1], sigma.temp[,,1]))) + sum(as.matrix(pdist(sigma.est[,,2], sigma.temp[,,2]))) + sum(as.matrix(pdist(sigma.est[,,3], sigma.temp[,,3]))) 
	names(tol.temp) <- "tolerance"
	i.iteration <- c(i.iteration, (i.iteration[length(i.iteration)]+1))
	
	mu.est <- mu.temp
	sigma.est <- sigma.temp
	tol <- bind_rows(tol, tol.temp)
	
	
# Contour plot of estimated multivariate Gaussian and data
	data.grid <- expand.grid(s.1 = seq(0, 20, length.out=200), s.2 = seq(0, 20, length.out=200))
	q1.samp <- cbind(data.grid, prob = dmvnorm(data.grid, mean = mu.est[1,], sigma = sigma.est[,,1]))
	q2.samp <- cbind(data.grid, prob = dmvnorm(data.grid, mean = mu.est[2,], sigma = sigma.est[,,2]))
	q3.samp <- cbind(data.grid, prob = dmvnorm(data.grid, mean = mu.est[3,], sigma = sigma.est[,,3]))

	print(ggplot() + 
		geom_contour(data=q1.samp,aes(x=s.1,y=s.2,z=prob)) +    
		geom_contour(data=q2.samp,aes(x=s.1,y=s.2,z=prob),col="red") +
		geom_contour(data=q3.samp,aes(x=s.1,y=s.2,z=prob),col="green") +
		geom_point(data=y1, aes(x=x, y=y), col="blue") + 
		geom_point(data=y2, aes(x=x, y=y), col="red") +
		geom_point(data=y3, aes(x=x, y=y), col="green"))
	}
