
install.packages('Radviz')
library(Radviz)
library(ggplot2)
library(dplyr)
library(tidyr)

install.packages('bodenmiller')
library(bodenmiller)
data(refPhenoMat)
data(refFuncMat)
data(refAnnots)
ref.df <- data.frame(refAnnots,
                     refPhenoMat,
                     refFuncMat)

