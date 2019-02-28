setwd('~/FS12/IgA/')
library(tidyverse)



IgA <- read.csv('IgA_Results.csv')

IgA <- filter(IgA, ng.IgA.mg.prot < 10000)

hist(IgA$ng.IgA.mg.prot, breaks = 50)
hist(IgA$ng.IgA.mg.dry, breaks = 100)

IgA$Diet <- factor(IgA$Diet)

IgA$mg.prot.mg.dry <- IgA$mg.prot.ml.stock / 50

#IgA %>% filter()

ggplot(IgA, aes(x=Diet, y=ng.IgA.mg.prot, group = Diet)) + geom_boxplot(aes(fill = Diet)) + ggtitle('ng IgA / mg protein: Cecal Contents')
ggplot(IgA, aes(x=Diet, y=ng.IgA.mg.dry, group = Diet)) + geom_boxplot(aes(fill = Diet)) + ggtitle('ng IgA / mg dry weight: Cecal Contents')

ggplot(IgA, aes(x=Diet, y=mg.prot.mg.dry, group = Diet)) + geom_boxplot(aes(fill = Diet)) + ggtitle('mg protein / mg dry weight: Cecal Contents', subtitle = 'BCA assay')

pairwise.wilcox.test(IgA$mg.prot.mg.dry, IgA$Diet, p.adjust.method = 'none')






