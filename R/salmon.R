# Introduction ------------------------------------------------------------
# In this script I will do some dimensionality reduction with Salmon output
# and look at differential gene expression.

library(tidyverse)

# Importing the data ------------------------------------------------------
salmon_counts <- read.table(file = file.path(getwd(), 
                                             "data",
                                             "salmon_output.txt"),
                            sep = '\t',
                            header = TRUE)


# MDS calcuations and plotting --------------------------------------------


