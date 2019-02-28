# setwd("~/FS12")
# setwd('~/Documents/FS12/')

require(pracma)
library(tidyverse)
library(Hmisc)
library(reshape2)
library(forcats) # dont think i use this
library(tidyverse)
library(cowplot)

# Read in data

sal_data <- read_csv('./data/FS12b_salmonella_data.csv')

# orders the treatment factors and makes a timepoint factor for plotting ease

sal_data$treatment <- factor(sal_data$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

sal_data$treatXtime <- factor(sal_data$treatXtime, levels = c('control_0', 'RPS_0', 'Acid_0', 'Zn+Cu_0', 'RCS_0', 'Bglu_0',
                                                              'control_2', 'RPS_2', 'Acid_2', 'Zn+Cu_2', 'RCS_2', 'Bglu_2',
                                                              'control_7', 'RPS_7', 'Acid_7', 'Zn+Cu_7', 'RCS_7', 'Bglu_7',
                                                              'control_14', 'RPS_14', 'Acid_14', 'Zn+Cu_14', 'RCS_14', 'Bglu_14',
                                                              'control_21', 'RPS_21', 'Acid_21', 'Zn+Cu_21', 'RCS_21', 'Bglu_21'))

sal_data$time_point_fact <- factor(sal_data$time_point)
# sal_data$time_point
#

sal_data %>% ggplot(aes(log_sal)) + geom_histogram()

## this is a nice one i think.....

sal_data %>% filter(time_point != 0 & pignum != 101) %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatXtime, fill=treatment)) +
  geom_boxplot(outlier.alpha = 0, position = position_dodge2(preserve = 'total')) + geom_jitter(aes(fill=treatment), width=.2,shape=21, size=2, stroke=1.1)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  facet_wrap(~time_point_fact) + ggtitle('Daily Shedding')


get_pairs <- function(df){
  pp <- pairwise.wilcox.test(df$log_sal, df$treatment, p.adjust.method = 'none')
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <- paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps
  
}

# DAILY SHEDDING WILCOX

daily_tests <- sal_data %>% filter(pignum != 101) %>% group_by(time_point) %>% nest() %>% mutate(pps = map(data, get_pairs)) %>% select(time_point, pps) %>% unnest()

daily_tests <- daily_tests %>% select(time_point, starts_with('control'))
# write.csv(daily_tests, file = 'Daily_shedding_wilcox.csv')

# PW_wilc_per_gene <- FS1.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data,  get_pairs)) %>% select(gene, pps) %>% unnest()



# sal_data %>% filter(pignum != 101) %>% ggplot(aes(x=time_point, y=log_sal)) + geom_line(aes(group=pignum, color=treatment), size=1)+scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))

sal_data %>% filter(pignum == 101 )


sal_data %>% filter(pignum != 101 ) %>% ggplot(aes(x=time_point, y=log_sal)) +
  geom_line(aes(group=pignum, color=treatment), size=1) + facet_wrap(~treatment) + geom_point()+ #geom_text(aes(label=pignum))+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))


# sal_data %>% filter(pignum != 101 & treatment %in% c('control', 'Zn+Cu')) %>% ggplot(aes(x=time_point, y=log_sal)) +
#   geom_line(aes(group=pignum, color=treatment), size=1) + facet_wrap(~treatment) + geom_point()+ #geom_text(aes(label=pignum))+
#   scale_color_manual(values=c('#33CC33', 'red', 'grey', 'purple'))


all_daily <- sal_data %>% filter(pignum != 101) %>% group_by(time_point, treatment) %>%
  summarise(mean_sal=mean(log_sal),
            sd_sal=sd(log_sal),
            num=n(),
            se_sal=sd_sal/sqrt(num))



ZN_CONT <- sal_data %>% filter(pignum != 101 & treatment %in% c('control', 'Zn+Cu')) %>%
  group_by(time_point, treatment) %>% 
  summarise(mean_sal=mean(log_sal), sd_sal=sd(log_sal), se_sal=sd_sal/sqrt(n()))


ZN_CONT %>% ggplot(aes(x=time_point, y=mean_sal, fill=treatment, group=treatment)) +geom_line(aes(color=treatment))+ geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) + 
  geom_point(shape=21, size=4) +scale_color_manual(values=c('#33CC33', 'red', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values=c('#33CC33', 'red', 'orange', 'red', 'grey', 'purple')) + ylab('log Salmonella') + annotate(x=21, y=3, geom='text', label='p=0.06')+ xlab('Day post-challenge')





all_daily %>% ggplot(aes(x=time_point, y=mean_sal, fill=treatment, group=treatment)) +geom_line(aes(color=treatment), size=1.5)+ geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) +   geom_point(shape=21, size=4) +scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ylab('log Salmonella') + xlab('Day post-challenge') + ggtitle('Daily shedding, group summary statistics')



all_daily %>% filter(treatment %in% c('control', 'RPS')) %>% ggplot(aes(x=time_point, y=mean_sal, fill=treatment, group=treatment)) +geom_line(aes(color=treatment))+ geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) + 
  geom_point(shape=21, size=4) +scale_color_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ylab('log Salmonella') + xlab('Day post-challenge') +
  annotate(x=21, y=3, geom='text', label='p=0.02')+
  annotate(x=7, y=3.25, geom='text', label='p=0.03')+
  annotate(x=2, y=3.75, geom='text', label='p=0.06')+xlab('Day post-challenge')


ZN_CONT_stats <- sal_data %>% filter(pignum != 101 & treatment %in% c('control', 'Zn+Cu'))


ZN_CONT_stats

pairwise.t.test(ZN_CONT_stats$log_sal,ZN_CONT_stats$treatXtime, p.adjust.method = 'none', var.equal=FALSE)
#### maybe some statistical test at each timepoint?


sum_sal <- sal_data %>% group_by(pignum) %>%
  summarise(AULC=trapz(time_point, log_sal),
            sum=sum(log_sal))

# this mess is to calculate AULC at each timepoint for each pig

AULCs <- list(sal_data %>% filter(time_point <=2 & pignum != 101)  %>% group_by(pignum) %>% summarise(AULC=trapz(time_point, log_sal)) %>% mutate(AULC=0, day=0),
     sal_data %>% filter(time_point <=2 & pignum != 101)  %>% group_by(pignum) %>% summarise(AULC=trapz(time_point, log_sal)) %>% mutate(day=2),
     sal_data %>% filter(time_point <=7 & pignum != 101) %>% group_by(pignum) %>% summarise(AULC=trapz(time_point, log_sal)) %>% mutate(day=7),
     sal_data %>% filter(time_point <=14 & pignum != 101)  %>% group_by(pignum) %>% summarise(AULC=trapz(time_point, log_sal)) %>% mutate(day=14),
     sal_data %>% filter(time_point <=21 & pignum != 101)  %>% group_by(pignum) %>% summarise(AULC=trapz(time_point, log_sal)) %>% mutate(day=21)) %>%
  bind_rows()

treats <- sal_data %>% select(pignum, treatment) %>% unique()
AULCs <- merge(AULCs, treats, by = 'pignum', all = TRUE)

AULCs$treatment <- factor(AULCs$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))


AULCs_av <- AULCs %>% group_by(day, treatment) %>% summarise(mean_AULC=mean(AULC), sd=sd(AULC), n=n(), sterr=sd/sqrt(n))

# AULCs_av$sd[is.na(AULCs_av$sd)] <- 0
# AULCs_av$sterr[is.na(AULCs_av$sterr)] <- 0

AULCs_av$treatment <- factor(AULCs_av$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

AULCs_av %>% filter(treatment %in% c('control', 'RPS')) %>% ggplot(aes(x=day, y=mean_AULC, group=treatment, color=treatment, fill = treatment)) +
  geom_point() +
  geom_line() + geom_ribbon(aes(ymin=mean_AULC-sterr, ymax=mean_AULC+sterr), alpha=.5) + 
  scale_color_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Mean AULC over time') + 
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) 


AULCs_av %>% filter(treatment %in% c('control', 'RPS')) %>% ggplot(aes(x=day, y=mean_AULC, group=treatment, color=treatment, fill = treatment)) +
  geom_point(size=3) +
  geom_line(size=1.5) + geom_errorbar(aes(ymin=mean_AULC-sterr, ymax=mean_AULC+sterr), size=1,width=1) + 
  scale_color_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('') + 
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ylab('AULC') 
#### for brad ####
# 
# AULCs_av %>% filter(treatment %in% c('control', 'Zn+Cu')) %>% ggplot(aes(x=day, y=mean_AULC, group=treatment, color=treatment, fill = treatment)) +
#   geom_line(size=1.5) + 
#   geom_errorbar(aes(ymin=mean_AULC-sterr, ymax=mean_AULC+sterr), color='black', size=.5,width=.5) + 
#   geom_point(size=3, shape=21, color='black') +
#   scale_color_manual(values=c('grey50', 'blue'), labels=c('Ctrl', 'Zn+Cu')) + ggtitle('') + 
#   scale_fill_manual(values=c('grey50', 'blue'), labels=c('Ctrl', 'Zn+Cu')) + labs(y=expression('AULC - cumulative CFU/g feces ('*log[10]*')'), 
#                                                                                           x='Days post-innoculation') + 
#   theme(axis.title.y=element_text(size=11.5, family = 'bold'), 
#         axis.title.x = element_text(size = 11.5, family= 'bold'), 
#         legend.title = element_blank())
#   
# 
# ZN_CONT %>% ggplot(aes(x=time_point, y=mean_sal, fill=treatment, group=treatment)) +geom_line(aes(color=treatment))+ geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.5) + 
#   geom_point(shape=21, size=4) +scale_color_manual(values=c('grey50', 'blue'), labels=c('Ctrl', 'Zn+Cu')) + ggtitle('') + 
#   scale_fill_manual(values=c('grey50', 'blue'), labels=c('Ctrl', 'Zn+Cu')) + labs(y=expression('CFU/g feces ('*log[10]*')'), 
#                                                                                   x='Days post-innoculation') + 
#   theme(axis.title.y=element_text(size=11.5, family = 'bold'), 
#         axis.title.x = element_text(size = 11.5, family= 'bold'), 
#         legend.title = element_blank()) + annotate(x=21, y=3, geom='text', label='p=0.06')
# 
# ZN_AULCs_av <- AULCs_av %>% filter(treatment %in% c('control', 'Zn+Cu'))
# ZN_AULCs <- AULCs %>% filter(treatment %in% c('control', 'Zn+Cu'))
# 
# write_tsv(ZN_AULCs, 'brad_AULCs.tsv')
# write_tsv(ZN_AULCs_av, 'brad_AULCs_sum_stats.tsv')
#########temp ##########


sum_sal <- sum_sal[match(filter(sal_data, time_point==2)$pignum,sum_sal$pignum),]



#
filter(sal_data, time_point==2)$pignum == filter(sal_data, time_point==7)$pignum
filter(sal_data, time_point==14)$pignum == filter(sal_data, time_point==21)$pignum
filter(sal_data, time_point==2)$pignum == filter(sal_data, time_point==21)$pignum
sum_sal$pignum == filter(sal_data, time_point==2)$pignum

#
sum_sal$d2_shed <- filter(sal_data, time_point==2)$log_sal
sum_sal$d7_shed <- filter(sal_data, time_point==7)$log_sal
sum_sal$d14_shed <- filter(sal_data, time_point==14)$log_sal
sum_sal$d21_shed <- filter(sal_data, time_point==21)$log_sal
sum_sal$d0_temp <- filter(sal_data, time_point==0)$temp
sum_sal$d2_temp <- filter(sal_data, time_point==2)$temp
sum_sal$d7_temp <- filter(sal_data, time_point==7)$temp
sum_sal$d2_Dtemp <- filter(sal_data, time_point==2)$temp - filter(sal_data, time_point==0)$temp
sum_sal$d7_Dtemp <- filter(sal_data, time_point==7)$temp - filter(sal_data, time_point==0)$temp
gain <- filter(sal_data, time_point == 21)$pig_weight - filter(sal_data, time_point == 0)$pig_weight
ADG <- gain/21
sum_sal$ADG <- ADG
sum_sal$init_weight <- filter(sal_data, time_point == 0)$pig_weight
sum_sal$final_weight <- filter(sal_data, time_point == 21)$pig_weight
sum_sal$P_init_gain <- (gain/sum_sal$init_weight)*100
sum_sal$P_dry_matter0 <- filter(sal_data, time_point == 0)$percent_drymatter
sum_sal$P_dry_matter2 <- filter(sal_data, time_point == 2)$percent_drymatter
sum_sal$P_dry_matter7 <- filter(sal_data, time_point == 7)$percent_drymatter


sum_sal <- merge(sum_sal, treats, by='pignum', all=TRUE)
sum_sal$treatment <- factor(sum_sal$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

hist(sum_sal$AULC, breaks = 20)


filter(sum_sal, pignum !=101) %>% ggplot(aes(x=treatment, y=AULC, fill=treatment))+
  geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(fill=treatment), shape=21, size=2, stroke=1.25, width = .12) + #geom_text(aes(label=pignum)) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Cumulative Salmonella shedding (AULC)', subtitle = 'Wilcoxon vs control: RPS p=0.013, Acid p=0.10, Bglu p=0.043')



filter(sum_sal, pignum !=101) %>% ggplot(aes(x=treatment, y=AULC, fill=treatment))+
  geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(fill=treatment), shape=21, size=2, stroke=1.25, width = .12) + #geom_text(aes(label=pignum)) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  theme_bw(base_size = 16)+ ggtitle('')

sum_sal.filter <- filter(sum_sal, pignum !=101)

####### FOr keystone poster #######
sum_sal.poster <- filter(sum_sal, pignum !=101 & treatment %in% c('control', 'RPS'))

sum_sal.poster %>% ggplot(aes(x=treatment, y=AULC, fill=treatment)) +
  geom_boxplot() + geom_jitter(width = .2, shape=21, size=2, stroke=1) +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Cumulative Salmonella shedding over 21 days', subtitle = 'Wilcox pvalue=0.013')




sum_sal_sum <- sum_sal.filter %>% group_by(treatment) %>%
  summarise(mean_AULC=mean(AULC), sd=sd(AULC), num_obs=n(), se=sd/sqrt(n()))


ggplot(sum_sal_sum, aes(x=treatment, y=mean_AULC)) + geom_col(aes(fill=treatment)) +
  geom_errorbar(aes(ymin=mean_AULC - se, ymax=mean_AULC + se), size=.5, width=.5) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Mean AULC', subtitle = 'Error bars represent +/- 1 standard error from the mean')


wilcox.test(sum_sal$AULC[sum_sal$treatment == 'control'], sum_sal$AULC[sum_sal$treatment == 'RPS'])
wilcox.test(sum_sal$AULC[sum_sal$treatment == 'control'], sum_sal$AULC[sum_sal$treatment == 'Acid'])





sal_data %>% filter(time_point < 14) %>% ggplot(aes(x=treatment, y=temp, fill=treatment, group=treatXtime)) +
  geom_boxplot() + geom_jitter(shape=21, size=2)+ geom_text(aes(label=pignum)) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  ggtitle('Rectal Temperature') + facet_wrap(~time_point)

sal_data %>% filter(time_point < 14) %>% ggplot(aes(x=time_point_fact, y=temp, fill=treatment, group=treatXtime)) +
  geom_boxplot() + geom_jitter(shape=21, size=2)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  ggtitle('Rectal temperatures') + facet_wrap(~treatment)#+ geom_text(aes(label=pignum))

#### oooh.... maybe can do an AULC for temperature


### weight figs

# this shows pig 101 not gaining as the other pigs are #
# use this at beginning
sal_data %>% filter(time_point %in% c(0,14)) %>%
  ggplot(aes(x=time_point, y=pig_weight, group=pignum, color=treatment)) + geom_line() +
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  geom_line(data = filter(sal_data, pignum ==101 & time_point %in% c(0,14)), color='black', size=1.25) +
  ggtitle('weight gain over the 1st two weeks', subtitle = 'Pig 101 shown in black') + ylim(0,32)


sal_data %>% filter(pignum !=101) %>% ggplot(aes(x=time_point_fact, y=pig_weight, fill=treatment, group=treatXtime)) +
  geom_boxplot() +scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  ggtitle('Weight') #+ geom_text(aes(label=pignum))


sal_data %>% filter(pignum !=101 & treatment %in% c('control', 'Bglu')) %>% ggplot(aes(x=time_point_fact, y=pig_weight, fill=treatment, group=treatXtime)) +
  geom_boxplot() +scale_fill_manual(values=c('#33CC33', 'purple', 'orange', 'red', 'grey', 'purple'))+
  ggtitle('Weight') #+ geom_text(aes(label=pignum))

# some sum_sal figs #

sum_sal %>% ggplot(aes(x=treatment, y=ADG, group=treatment, fill=treatment)) +
  geom_boxplot() +scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))


sum_sal %>% ggplot(aes(x=treatment, y=P_init_gain, group=treatment, fill=treatment)) +
  geom_boxplot() +scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Percent Initial gain')

################################
# setwd("~/FS12")

tis <- read.csv('./data/D21_tissues.csv')
tis$log_sal <- log10(tis$Salmonella)
tis$log_sal <- as.numeric(sub(-Inf, 0, tis$log_sal))
tis <- na.exclude(tis)

tis$treatment <- factor(tis$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))


pairwise.wilcox.test(tis$log_sal[tis$tissue == 'cecal_cont'], tis$treatment[tis$tissue == 'cecal_cont'], p.adjust.method = 'none')

pairwise.t.test(tis$log_sal[tis$tissue == 'cecal_cont'], tis$treatment[tis$tissue == 'cecal_cont'], p.adjust.method = 'none')


tis %>% filter(tissue=='cecal_cont') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0)+ geom_jitter(shape=21,width = .1, size=2.25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Cecal Contents')


tis %>% filter(pignum ==97)


###### For Shawn FSIS  ##########

ytitle <- expression(paste('Log ', italic("Salmonella")))

tis %>% filter(tissue=='cecal_cont') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0)+ geom_jitter(shape=21,width = .1, size=2.25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ylab(ytitle)


tis %>% filter(tissue=='cecal_cont') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0)+ geom_jitter(shape=21,width = .1, size=2.25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  theme_bw(base_size = 16) + ylab(ytitle)



tis %>% filter(tissue=='Cecum') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0)+ geom_jitter(shape=21,width = .2, size=2.25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Cecum')


tis %>% filter(tissue=='ICLN') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0)+ geom_jitter(shape=21,width = .1, size=2.25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('ICLN')

tis %>% filter(tissue=='IPP') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0)+ geom_jitter(shape=21,width = .1, size=2.25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('IPP')

tis %>% filter(tissue=='Tonsil') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(shape=21,width = .1, size=2.25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Tonsil')

tis$tisXtreat <- paste(tis$tissue, tis$treatment, sep = '_')

tis.sum <- tis %>% group_by(tisXtreat) %>%
  summarise(mean=mean(log_sal), sd=sd(log_sal), num_obs=n(), se=sd/sqrt(n()))

tis.sum$tisXtreat

tis.sum$tissue <- sub('(.*)_(.*)', '\\1', tis.sum$tisXtreat)
tis.sum$treatment <- sub('(.*)_(.*)', '\\2', tis.sum$tisXtreat)

tis.sum$treatment <- factor(tis.sum$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

ggplot(tis.sum, aes(x=treatment, y=mean)) + geom_col(aes(fill=treatment)) +
  geom_errorbar(aes(ymin=mean- se, ymax=mean + se), size=.25, width=.35, color = 'black', alpha=.85) +
  facet_wrap(~tissue) +  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))

###### Shawn FSIS Fig #####

ytitle <- expression(paste('Log ', italic("Salmonella")))

tis.sum %>% filter(tissue == 'cecal_cont') %>% ggplot(aes(x=treatment, y=mean)) + geom_col(aes(fill=treatment)) +
  geom_errorbar(aes(ymin=mean- se, ymax=mean + se), size=.5, width=.35, color = 'black', alpha=.85) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  ylab(ytitle <- expression(paste('Log ', italic("Salmonella")))) + 
  theme_bw(base_size = 16) + ggtitle('')

tis.sum %>% filter(tissue == 'cecal_cont') %>% ggplot(aes(x=treatment, y=mean)) + geom_col(aes(fill=treatment)) +
  geom_errorbar(aes(ymin=mean- se, ymax=mean + se), size=.5, width=.35, color = 'black', alpha=.85) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  ylab(ytitle <- expression(paste('Log ', italic("Salmonella")))) + 
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('')




tis.sp <- tis %>% select(-tisXtreat, -Salmonella) %>% spread(key = tissue, value = log_sal)

sum_sal <- sum_sal[match(tis.sp$pignum, sum_sal$pignum),]

sum_sal$pignum == tis.sp$pignum
tis.sp <- tis.sp[,-2]
sal_for_cor <- merge(sum_sal, tis.sp, by = 'pignum')

### NEED CECAL VFA DATA FOR THIS NEXT SECTION ###




sal_for_cor <- merge(sal_for_cor, vfas_for_cor, by = 'pignum')
rownames(sal_for_cor) <- sal_for_cor$pignum



ggplot(sal_for_cor, aes(x=butyrate, y=AULC)) +
  geom_point(aes(color=treatment)) + geom_smooth(method = 'lm') + scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Correlation between cecal butyrate concentration and total Salmonella shedding', subtitle = 'Spearman: -0.4, p=0.003')


ggplot(sal_for_cor, aes(x=caproate, y=AULC)) +
  geom_point(aes(color=treatment)) + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Correlation between cecal caproate concentration and total Salmonella shedding', subtitle = 'Spearman: -0.53, p=0.0005')
ggplot(sal_for_cor, aes(x=valerate, y=AULC)) +
  geom_point(aes(color=treatment)) + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Correlation between cecal valerate concentration and total Salmonella shedding', subtitle = 'Spearman: -0.38, p=0.0005')
ggplot(sal_for_cor, aes(x=total, y=AULC)) +
  geom_point(aes(color=treatment)) + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Correlation between cecal total SCFA concentration and total Salmonella shedding', subtitle = 'Spearman: -0.4, p=0.002')

##### Shawn FSIS #######

ggplot(sal_for_cor, aes(x=butyrate, y=AULC)) +
  geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme_bw(base_size = 16) + xlab('butyrate (mM)') + ggtitle('AULC vs Cecal butyrate (D21)')

ggplot(sal_for_cor, aes(x=total, y=AULC)) +
  geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme_bw(base_size = 16) + xlab('Total SCFAs (mM)') + ggtitle('')



ggplot(sal_for_cor, aes(x=butyrate, y=AULC, color=treatment)) +
  geom_smooth(method = 'lm', se = FALSE) + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
    theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('') + xlab('butyrate (mM)')

ggplot(sal_for_cor, aes(x=valerate, y=AULC, color=treatment)) +
  geom_smooth(method = 'lm', se = FALSE) + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('') + xlab('valerate (mM)')

ggplot(sal_for_cor, aes(x=caproate, y=AULC, color=treatment)) +
  geom_smooth(method = 'lm', se = FALSE) + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('') + xlab('caproate (mM)')

ggplot(sal_for_cor, aes(x=total, y=AULC, color=treatment)) +
  geom_smooth(method = 'lm', se = FALSE) + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('') + xlab('total (mM)')


############ TESTESTESTEST #################



sal_for_cor %>% filter(treatment == 'RPS') %>% ggplot(aes(x=butyrate, y=AULC)) +
  geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('') + xlab('butyrate (mM)')


sal_for_cor %>% filter(treatment == 'Bglu') %>% ggplot(aes(x=butyrate, y=AULC)) +
  geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('') + xlab('butyrate (mM)')

bglu <- sal_for_cor %>% filter(treatment == 'Bglu')
RPS <- sal_for_cor %>% filter(treatment == 'RPS')
cont <- sal_for_cor %>% filter(treatment == 'control')

test
fads <- cor.test(x = bglu$AULC, y = bglu$butyrate)
fads$statistic
fads$p.value

cor.test(sal_for_cor$AULC, sal_for_cor$butyrate, method = 'spear')

sal_for_cor %>% group_by(treatment) %>% summarise(butP = cor.test(AULC, butyrate)$p.value, 
                                                  butT = cor.test(AULC, butyrate)$statistic,
                                                  propP = cor.test(x=AULC, y=propionate)$p.value,
                                                  propT = cor.test(x=AULC, y=propionate)$statistic,
                                                  capP = cor.test(x=AULC, y=caproate)$p.value,
                                                  capT = cor.test(x=AULC, y=caproate)$statistic,
                                                  valP = cor.test(x=AULC, y=valerate)$p.value,
                                                  valT = cor.test(x=AULC, y=valerate)$statistic,
                                                  totP = cor.test(x=AULC, y=total)$p.value,
                                                  totT = cor.test(x=AULC, y=total)$statistic,
                                                  Pbutp = cor.test(x=AULC, y=P_butyrate)$p.value,
                                                  PbutT = cor.test(x=AULC, y=P_butyrate)$statistic)



ggplot(sal_for_cor, aes(x=total, y=AULC)) +
  geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('') + xlab('Total SCFAs (mM)')

filter(sal_for_cor, treatment %in% c('control', 'RPS'))
filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')


filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')



filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')


filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')


cor.test(filter(sal_for_cor, treatment =='control')$AULC,
         filter(sal_for_cor, treatment =='control')$butyrate)

cor.test(filter(sal_for_cor, treatment =='RPS')$AULC,
         filter(sal_for_cor, treatment =='RPS')$butyrate)


cor.test(filter(sal_for_cor, treatment =='control')$AULC,
         filter(sal_for_cor, treatment =='control')$valerate)

cor.test(filter(sal_for_cor, treatment =='RPS')$AULC,
         filter(sal_for_cor, treatment =='RPS')$valerate)

cor.test(filter(sal_for_cor, treatment =='control')$AULC,
         filter(sal_for_cor, treatment =='control')$caproate)

cor.test(filter(sal_for_cor, treatment =='RPS')$AULC,
         filter(sal_for_cor, treatment =='RPS')$caproate)

cor.test(filter(sal_for_cor, treatment =='control')$AULC,
         filter(sal_for_cor, treatment =='control')$total)

cor.test(filter(sal_for_cor, treatment =='RPS')$AULC,
         filter(sal_for_cor, treatment =='RPS')$total)



cor.test(filter(sal_for_cor, treatment%in%c('RPS', 'control'))$AULC,
         filter(sal_for_cor, treatment%in%c('RPS', 'control'))$total)



cor.test(filter(sal_for_cor, treatment%in%c('RPS', 'control'))$AULC,
         filter(sal_for_cor, treatment%in%c('RPS', 'control'))$butyrate)



cor.test(filter(sal_for_cor, treatment%in%c('RPS', 'control'))$AULC,
         filter(sal_for_cor, treatment%in%c('RPS', 'control'))$valerate)



cor.test(filter(sal_for_cor, treatment%in%c('RPS', 'control'))$AULC,
         filter(sal_for_cor, treatment%in%c('RPS', 'control'))$caproate)






filter(sal_for_cor, treatment =='RPS') %>% cor.test()

poster <- sal_for_cor %>% select(-starts_with('P_')) %>% 
  gather(variable, response, -pignum, -treatment) %>% 
  filter(treatment %in% c('control', 'RPS')) %>% 
  filter(variable %in% c('AULC', 'butyrate', 'valerate', 'caproate', 'total')) %>% spread(variable, response) %>% 
  gather(variable, response, -pignum, -treatment, -AULC)
poster$response <- poster$response *3


poster$variable <- factor(poster$variable, levels = c('butyrate', 'valerate', 'caproate', 'total'))
ggplot(poster, aes(x=response, y=AULC, fill=treatment)) +
  geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + facet_wrap(~variable, scales = 'free') +
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  geom_point(shape=21, stroke=1, size=3) + ggtitle('Cecal SCFA concentrations VS cumulative Salmonella shedding') + xlab('Concentration (mM)')


poster %>% filter(variable == 'butyrate') %>% ggplot(aes(x=response, y=AULC, fill=treatment)) +
  geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + #facet_wrap(~variable, scales = 'free') +
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  geom_point(shape=21, stroke=1, size=3) + ggtitle('') + xlab('Concentration (mM)')
#3399FF


poster %>% filter(variable == 'caproate') %>% ggplot(aes(x=response, y=AULC, fill=treatment)) +
  geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + #facet_wrap(~variable, scales = 'free') +
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  geom_point(shape=21, stroke=1, size=3) + ggtitle('') + xlab('Concentration (mM)')



poster %>% filter(variable == 'valerate') %>% ggplot(aes(x=response, y=AULC, fill=treatment)) +
  geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + #facet_wrap(~variable, scales = 'free') +
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  geom_point(shape=21, stroke=1, size=3) + ggtitle('') + xlab('Concentration (mM)')



poster %>% filter(variable == 'total') %>% ggplot(aes(x=response, y=AULC, fill=treatment)) +
  geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + #facet_wrap(~variable, scales = 'free') +
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  geom_point(shape=21, stroke=1, size=3) + ggtitle('') + xlab('Concentration (mM)')





cecal_cont
#########################


sal_for_cor <- sal_for_cor[,-1]

# I think I could have just merged on pignum here instead of all this bullshit
#sal_for_cor <- tis.sp %>% select(-pignum, -treatment) %>% bind_cols(sum_sal,.) %>% select(-pignum, -treatment)


# pignumblers <- sum_sal$pignum
#
# sum_sal <- sum_sal[,-c(1,17)]
#
sal_for_cor.m <- as.matrix(sal_for_cor)
# rownames(sum_sal.m) <- pignumblers

salcor <- rcorr(sal_for_cor.m)
salcor$P
apply(sal_for_cor, MARGIN = 2, str_trim) %>% apply(MARGIN = 2, as.numeric) %>% as.matrix()
str_trim()

salcor.sigs <- rcorr_to_ggnet(salcor, spearcut = -1, pcut = 0.05)

library(funfuns)
library(geomnet)
################

nodes <- rbind(gather_nodes(sal_for_cor.m, 'Sal'))

nodes$node <- as.character(nodes$node)
#nodes$node[770] <- 'y'
alll <- na.exclude(fortify(as.edgedf(salcor.sigs), nodes))
na.exclude(alll)
alll$spear
rep(1, length(alll$spear))
test <- data.frame(rep(1, length(alll$spear)),alll$spear)
colnames(test) <- c('x', 'spear')

filts_pos <- filter(alll, spear > 0)
filts_neg <- filter(alll, spear < 0)

# butyrate caproate valerate neg cor w/ AULC





ggplot()

tetetete <- ggplot(test, aes(x=x, y=spear, color=spear)) + geom_point() + scale_color_gradient2()

hmm <- ggplot_build(tetetete)

hmm$data[[1]]$colour



salnet <- ggplot(alll, aes(from_id = from_id, to_id = to_id, label=from, color=type, ecolour=spear)) +
  geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)),
           aes(color = type, label = from_id),
           linewidth = 2, size = 5, vjust = 0, alpha = 0.3,
           repel = TRUE, fontsize=3, singletons = FALSE,labelcolour="black",
           labelgeom = 'text', ecolour = hmm$data[[1]]$colour) +
  theme_net()

salnet

##################

# Compare the ADG in different intervals back to Salmonella colonization in different tissues
# prelavence of Salmonella in farms tend to have lower ADG.  Can we see some of this in our pigs here?

# can we see some

# Do larger pigs tend to have greater amounts of salmonella





##### some correlation stuff here ####



# prob need to replace some of this stuff #


### change in rectal temp? ###

filter(sal_data, time_point ==0)$temp


delt_temp_2 <- filter(sal_data, time_point ==2)$temp - filter(sal_data, time_point ==0)$temp
delt_temp_7 <- filter(sal_data, time_point ==7)$temp - filter(sal_data, time_point ==0)$temp

sal_data$delta_temp <- NA
sal_data[sal_data$time_point ==2,]$delta_temp <- delt_temp_2
sal_data[sal_data$time_point ==7,]$delta_temp <- delt_temp_7



tis$PA <- ifelse(tis$Salmonella > 0, 1, 0)


tis




fisher_map <- function(data, row_vector, col_vector){
  fisher.test(data[row_vector, col_vector])$p.value
}


PA_tab <- tis %>% select(treatment, tissue, PA) %>%
  group_by(tissue, treatment) %>% 
  summarise(pos = sum(PA == 1), neg = sum(PA == 0)) %>% group_by(tissue) %>% nest() %>% 
  mutate(RPSp = map(data, fisher_map, c(1,2), c(2,3)), 
         Acidp = map(data, fisher_map, c(1,3), c(2,3)),
         Znp = map(data, fisher_map, c(1,4), c(2,3)), 
         RCSp = map(data, fisher_map, c(1,5), c(2,3)), 
         Bglup = map(data, fisher_map, c(1,6), c(2,3))) %>% select(-data) %>% unnest()

ftes <- fisher.test(PA_tab[[2]][[1]][c(1,3),-1])

ftes$data.name

row.
cecconttest <- PA_tab[[2]][[1]]


muc <- data.frame(SB1434 = c(1050, 0), 
           SX240 = c(1832,5), row.names = c('normal', 'mucoid'))


intra <- data.frame(norm=c(1175,0), 
                    mono=c(1832, 5), 
                    row.names = c('normal', 'mucoid'))

t(muc)
fisher.test(muc)
fisher.test(t(muc))

t(intra)
fisher.test(intra)
fisher.test(t(intra))


