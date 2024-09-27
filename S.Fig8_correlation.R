library(tidyverse)
library(ggplot2)

data <- read.csv("msk.csv")
data <- filter(data, AGE >0) 
p <- ggscatter(data, x = "AGE", y = "TMB",  #y = "CNV"
               add = "reg.line",  
               conf.int = TRUE,   
               cor.coef = TRUE,   
               cor.method = "spearman") +
  theme_test()
p

  
