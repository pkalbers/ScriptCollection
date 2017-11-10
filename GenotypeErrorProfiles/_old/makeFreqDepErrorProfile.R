#
# generate frequency dependent error profile
#

library(ggplot2)
library(grid)
library(data.table)


args = commandArgs(T)

pro.file = args[1]