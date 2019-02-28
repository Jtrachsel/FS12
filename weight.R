setwd("~/FS12")

D14_weight <- read.csv('D14_weight.csv', header = TRUE)
diets <- c(1:10)

names(diets) <- c('NC', 'PC', 'R-Potato', 'Maltodex', 'SBP', 'FA-mix', 'Phytogen', 'Zn+Cu', 'R-Corn', 'B-glu')
dietnames <- names(diets[D14_weight$Diet])
D14_weight$Diet <- dietnames



library(tidyverse)

weight.gather <- gather(D14_weight, key=measure, value=Kg, -Pen,-Diet)

weight.gather$Diet <- factor(weight.gather$Diet, levels = c('NC', 'PC', 'R-Potato', 'Maltodex', 'SBP', 'FA-mix', 'Phytogen', 'Zn+Cu', 'R-Corn', 'B-glu'))

weight.gather %>% filter(measure=='ADFI') %>% ggplot(aes(x=Diet, y=Kg, group=Diet, fill=Diet)) + geom_boxplot() +geom_point(aes(color=Diet), size=6, show.legend = FALSE)+ geom_text(aes(label=Pen)) + ggtitle('ADFI')
weight.gather %>% filter(measure=='ADG') %>% ggplot(aes(x=Diet, y=Kg, group=Diet, fill=Diet)) + geom_boxplot() +geom_point(aes(color=Diet), size=6, show.legend = FALSE)+ geom_text(aes(label=Pen)) + ggtitle('ADG')
weight.gather %>% filter(measure=='AV.Int') %>% ggplot(aes(x=Diet, y=Kg, group=Diet, fill=Diet)) + geom_boxplot() +geom_point(aes(color=Diet), size=6, show.legend = FALSE)+ geom_text(aes(label=Pen)) + ggtitle('Av.Int')
weight.gather %>% filter(measure=='GF') %>% ggplot(aes(x=Diet, y=Kg, group=Diet, fill=Diet)) + geom_boxplot() +geom_point(aes(color=Diet), size=6, show.legend = FALSE)+ geom_text(aes(label=Pen)) + ggtitle('GF')
weight.gather %>% filter(measure=='X14d.ABW') %>% ggplot(aes(x=Diet, y=Kg, group=Diet, fill=Diet)) + geom_boxplot()+geom_point(aes(color=Diet), size=6, show.legend = FALSE)+ geom_text(aes(label=Pen)) + ggtitle('14d.ABW')


############ All weight  ###########

all_weight <- read.csv('all_weight.csv', header = TRUE)
all_weight <- all_weight[!is.na(all_weight$Int.BW),]

all_weight$AvDGain <- (all_weight$X14d.BW - all_weight$Int.BW)/14

pig_weights <- all_weight[,-c(5,6,8,9,10,11,12,13)] %>% gather(key=measure, value=Kg, -Pen,-Diet, -Pig)


pig_weights$Diet <- names(diets[pig_weights$Diet])

pig_weights$Diet <- factor(pig_weights$Diet, levels = c('NC', 'PC', 'R-Potato', 'Maltodex', 'SBP', 'FA-mix', 'Phytogen', 'Zn+Cu', 'R-Corn', 'B-glu'))

pig_weights %>% filter(measure=='AvDGain') %>%
  ggplot(aes(x=Diet, y=Kg, group=Diet, fill=Diet)) +
  geom_boxplot() +geom_point(aes(color=Diet), size=6, show.legend = FALSE)+
  geom_text(aes(label=Pig)) + ggtitle('AvDGain')


pig_weights %>% filter(measure=='AvDGain') %>%
  ggplot(aes(x=Pen, y=Kg, group=Pen, fill=Diet)) +
  geom_boxplot() +geom_point(aes(color=Diet), size=4, show.legend = FALSE)+
  geom_text(aes(label=Pig), size=3) + ggtitle('AvDGain')


ggplot(all_weight, aes(x=Int.BW, y=AvDGain)) + geom_point() + geom_smooth(method = 'loess')+ ggtitle('Average Daily Gain VS Initial Birthweight', subtitle = "Pearson's correlation: coefficient = 0.018, p=0.6851")
cor.test(all_weight$Int.BW, all_weight$AvDGain)
