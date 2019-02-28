# setwd('')
library(tidyverse)


vfas <- read_csv('FS12b_vfas.csv')

vfas$treatment <- factor(vfas$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

vfas.gather <- vfas %>% gather(key = VFA, value = mM, -(pignum:hour))

vfas.gather$set <- paste(vfas.gather$hour, vfas.gather$treatment)

# initial peek #

filter(vfas.gather, hour == 0) %>% ggplot(aes(x=treatment, y=mM, group=set, fill=treatment)) +
  geom_boxplot(outlier.alpha = 0) +geom_jitter(shape=21, size=1.5, stroke=1, alpha=.75, width = .25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  facet_wrap(~VFA, scales = 'free') + ggtitle("Cecal SCFAs 21 days post challenge, no incubation")

filter(vfas.gather, hour == 24) %>% ggplot(aes(x=treatment, y=mM, group=set, fill=treatment)) +
  geom_boxplot(outlier.alpha = 0) +geom_jitter(shape=21, size=1.5, stroke=1, alpha=.75, width = .25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  facet_wrap(~VFA, scales = 'free') + ggtitle("Cecal SCFAs 21 days post challenge, no incubation")

### propionate in 24hr incubation seems to have suffered in ZN

# 0 and 24 hr together

filter(vfas.gather, VFA != 'succinate') %>% ggplot(aes(x=hour, y=mM, group=set, fill=treatment)) +
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  facet_wrap(~VFA, scales = 'free')

###

# tests for RPS vs control

buttest <- vfas.gather %>% filter(treatment %in% c('control', 'RPS') & hour ==0 & VFA == 'butyrate')
valtest <- vfas.gather %>% filter(treatment %in% c('control', 'RPS') & hour ==0 & VFA == 'valerate')
captest <- vfas.gather %>% filter(treatment %in% c('control', 'RPS') & hour ==0 & VFA == 'caproate')
tottest <- vfas.gather %>% filter(treatment %in% c('control', 'RPS') & hour ==0 & VFA == 'total')

wilcox.test(buttest$mM~buttest$treatment) # p-value = 0.00399
wilcox.test(valtest$mM~valtest$treatment) # p-value = 0.0007816
wilcox.test(captest$mM~captest$treatment) # p-value = 0.01061
wilcox.test(tottest$mM~tottest$treatment) # p-value = 0.1135


### for correlations ###

vfas0 <- filter(vfas, hour == 0) %>% select(-hour)

vfas24 <- filter(vfas, hour == 24) %>% select(-hour)

colnames(vfas24) <- paste('P', colnames(vfas24), sep = '_') # why 'P' ?

colnames(vfas24)[1] <- 'pignum'

vfas_for_cor <- merge(vfas0, vfas24, by = 'pignum')

vfas_for_cor$day <- 21
vfas_for_cor <- vfas_for_cor %>% select(pignum, treatment, day, everything())
