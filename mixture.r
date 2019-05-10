## Mixture of experts (linear regression)

# Clear workspace
rm(list = ls())

# Load required packages
library(tidyverse)
# Load data

## Initial test dataset with socio-economic / demographic factors
detData <- read_csv("S:/HelenJohnson/Herpes Zoster/Data/testData.csv", 
			col_names = TRUE)
