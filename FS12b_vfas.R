# setwd('')
library(tidyverse)


vfas <- read_csv('./data/FS12b_vfas.csv')

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

# filter(vfas.gather, VFA != 'succinate') %>% ggplot(aes(x=hour, y=mM, group=set, fill=treatment)) +
#   geom_boxplot(outlier.alpha = 0) +
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   facet_wrap(~VFA, scales = 'free')
# 
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

colnames(vfas24) <- paste('P', colnames(vfas24), sep = '_') # why 'P' ?  Because it shows the 'P'otential amount of VFAs for that bolus

colnames(vfas24)[1] <- 'pignum'

vfas_for_cor <- merge(vfas0, vfas24, by = 'pignum')

vfas_for_cor$day <- 21
vfas_for_cor$tissue <- 'cecal_contents'
vfas_for_cor <- vfas_for_cor %>% select(pignum, treatment, day, tissue, everything())


write_csv(vfas_for_cor, './data/FS12b_vfas_for_cor.csv')

############### FECAL DATA #################




library(tidyverse)


####### WARNING! I THINK THERE IS A PROBLEM WITH THE DATA FROM DAY 7 and 21!!!!!

### SUSPICIOUS AMOUNTS OF LACTATE AND SUCCINATE ARE PRESENT, USUALLY THESE ARE VERY RARE
### I THINK THERE WAS AN ISSUE DURING SAMPLE PREP
### MY GUESS IS THAT THE SAMPLES WERE ALLOWED TO THAW AT ROOM TEMP FOR TOO LONG AND 
### SOME BACTERIAL METABOLIC ACTIVITY GENERATED EXCESS LACTATE AND SUCCINATE
### USUALLY THESE WOULD BE CONVERTED TO OTHER SCFAS BY ANAEROBIC FERMENTATION BUT
### THE ANAEROBES PROBABLY DIED DUE TO OXYGEN EXPOSURE.


# write_csv(res.all, 'FS12b_fecal_vfas.csv')

res.all <- read_csv('./data/FS12b_fecal_vfas.csv')
tests <- filter(res.all, time == 0) %>% do(pwilx=pairwise.wilcox.test(.$butyrate, .$treatment, p.adjust.method = 'none'))
#pairwise.wilcox.test(p.adjust.method = )
tests[[1]]


tests <- filter(res.all, time == 0) %>% do(pwilx=pairwise.wilcox.test(.$caproate, .$treatment, p.adjust.method = 'none'))
#pairwise.wilcox.test(p.adjust.method = )
tests[[1]]


tests <- filter(res.all, time == 0) %>% do(pwilx=pairwise.wilcox.test(.$valerate, .$treatment, p.adjust.method = 'none'))
#pairwise.wilcox.test(p.adjust.method = )
tests[[1]]


res.all$treatment <- factor(res.all$treatment, levels=c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))


res.all %>% ggplot(aes(x=treatment, y=acetate, fill=treatment)) + facet_wrap(~time) + geom_boxplot() + geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=propionate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=butyrate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=valerate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=caproate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=isovalerate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=oxalate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=phenylacetate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=succinate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=fumarate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% ggplot(aes(x=treatment, y=lactate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()+ geom_text(aes(label=pignum))


res.gather <- res.all %>% gather(key = variable, value = concentration, -(pignum:set))

# res.gather %>% filter(variable == 'total')%>% ggplot(aes(x=time, y=concentration, group=pignum, color=treatment)) + geom_line()

res.sum <- res.gather %>%  group_by(time, treatment, variable) %>%
  summarise(mean=mean(concentration), sd=sd(concentration), n=n(), sterr=sd/sqrt(n))

res.sum <- res.sum %>% mutate(treat_var=paste(treatment, variable, sep='_'))

res.sum$treatment <- factor(res.sum$treatment, levels= c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

res.sum %>% ggplot(aes(x=time, y=mean, group=treat_var, color=treatment)) +
  geom_line() + facet_wrap(~variable, scales='free') + scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Mean concentration of fecal SCFAs over time')
 



