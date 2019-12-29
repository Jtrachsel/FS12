require(pracma)
library(tidyverse)
library(Hmisc)
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

sal_data <- sal_data %>% filter(!(treatment %in% c('Zn+Cu', 'Bglu')) & pignum != 101)


# sal_data %>% filter(time_point != 0 & !(treatment %in% c('Zn+Cu', 'Bglu'))) %>% ggplot(aes(log_sal, fill=treatment)) + geom_histogram() + facet_wrap(~time_point, ncol=1)
# sal_data %>% filter(time_point != 0 & !(treatment %in% c('Zn+Cu', 'Bglu'))) %>% ggplot(aes(log_sal, fill=treatment)) + geom_histogram()

## this is a nice one i think.....

sal_data %>% filter(time_point != 0 & pignum != 101) %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatXtime, fill=treatment)) +
  geom_boxplot(outlier.alpha = 0, position = position_dodge2(preserve = 'total')) + geom_jitter(aes(fill=treatment), width=.2,shape=21, size=2, stroke=1.1)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  facet_wrap(~time_point_fact) + ggtitle('Daily Shedding') + theme_bw() + ylab('Log10 Salmonella CFUs / g feces')

sal_data %>% #filter(time_point != 0 & pignum != 101) %>%
  ggplot(aes(x=time_point_fact, y=log_sal, group=treatXtime, fill=treatment)) +
  geom_boxplot(outlier.alpha = 0, position = position_dodge2(preserve = 'total')) + geom_jitter(aes(fill=treatment), width=.2,shape=21, size=2, stroke=1.1)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  facet_wrap(~treatment) + ggtitle('Daily Shedding')+ theme_bw()

# sal_data %>% #filter(time_point != 0 & pignum != 101) %>%
#   ggplot(aes(x=time_point_fact, y=log_sal, group=treatXtime, fill=treatment)) +
#   geom_boxplot(outlier.alpha = 0, position = position_dodge2(preserve = 'total')) + geom_jitter(aes(fill=treatment), width=.2,shape=21, size=2, stroke=1.1)+
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   ggtitle('Daily Shedding')+ theme_bw()



get_pairs <- function(df){
  pp <- pairwise.wilcox.test(df$log_sal, df$treatment, p.adjust.method = 'none')
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <- paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps
  
}

# DAILY SHEDDING WILCOX

daily_tests <- sal_data %>% filter(pignum != 101) %>%
  group_by(time_point) %>% nest() %>% mutate(pps = map(data, get_pairs)) %>%
  select(time_point, pps) %>% unnest()

daily_tests <- daily_tests %>% select(time_point, starts_with('control'))
# write.csv(daily_tests, file = 'Daily_shedding_wilcox.csv')

# PW_wilc_per_gene <- FS1.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data,  get_pairs)) %>% select(gene, pps) %>% unnest()


### Just Bglu stuff here ###

# test <- sal_data %>% filter(time_point == 21 & treatment %in% c('control', 'Bglu')) %>% 
#   write_tsv('Sal_bglu_control_D21.tsv')
# 
# test <- read_tsv('Sal_bglu_control_D21.tsv')
# AOV <- aov(data = test, formula = log_sal ~ treatment)
# summary(AOV)
# TukeyHSD(AOV)
# 
# pairwise.t.test(test$log_sal, test$treatment)
# 
# 
# 
# test2 <- sal_data %>% filter(treatment %in% c('control', 'Bglu') & time_point != 0) %>% 
#   write_tsv('Sal_bglu_control.tsv')
# 
# test <- read_tsv('Sal_bglu_control_D21.tsv')
# AOV <- aov(data = test, formula = log_sal ~ treatment)
# summary(AOV)
# TukeyHSD(AOV)
# 
# pairwise.t.test(test$log_sal, test$treatment)
# 
# 
# 
# 
# summary(aov(log_sal ~ treatment + Error(pignum/time_point_fact), data=test2))
# 
# 
#

library(lmerTest)
library(psycho)
# test2$treatment <- factor(test2$treatment, levels = c('control', 'Bglu'))
# # install.packages('psycho')
# fit <- lmer(log_sal ~ treatment + (1|pignum), data=test2)
# anova(fit)
# 
# results <- analyze(fit)
# print(results)
# 
# summary(fit)


##############


# sal_data %>% filter(pignum != 101) %>% ggplot(aes(x=time_point, y=log_sal)) + geom_line(aes(group=pignum, color=treatment), size=1)+scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))

# sal_data %>% filter(pignum == 101 )


sal_data %>% filter(pignum != 101 ) %>%
  ggplot(aes(x=time_point, y=log_sal)) +
  geom_line(aes(group=pignum, color=treatment), size=1) +
  facet_wrap(~treatment) + geom_point()+ #geom_text(aes(label=pignum))+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme_bw()




# sal_data %>% filter(pignum != 101 & treatment %in% c('control', 'Zn+Cu')) %>% ggplot(aes(x=time_point, y=log_sal)) +
#   geom_line(aes(group=pignum, color=treatment), size=1) + facet_wrap(~treatment) + geom_point()+ #geom_text(aes(label=pignum))+
#   scale_color_manual(values=c('#33CC33', 'red', 'grey', 'purple'))


all_daily <- sal_data %>% filter(pignum != 101) %>% group_by(time_point, treatment) %>%
  summarise(mean_sal=mean(log_sal),
            sd_sal=sd(log_sal),
            num=n(),
            se_sal=sd_sal/sqrt(num))



# Just ZnCu vs Control 
# ZN_CONT <- sal_data %>% filter(pignum != 101 & treatment %in% c('control', 'Zn+Cu')) %>%
#   group_by(time_point, treatment) %>% 
#   summarise(mean_sal=mean(log_sal), sd_sal=sd(log_sal), se_sal=sd_sal/sqrt(n()))
# 
# 
# ZN_CONT %>% ggplot(aes(x=time_point, y=mean_sal, fill=treatment, group=treatment)) +geom_line(aes(color=treatment))+ geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) + 
#   geom_point(shape=21, size=4) +scale_color_manual(values=c('#33CC33', 'red', 'orange', 'red', 'grey', 'purple')) + 
#   scale_fill_manual(values=c('#33CC33', 'red', 'orange', 'red', 'grey', 'purple')) + ylab('log Salmonella') + annotate(x=21, y=3, geom='text', label='p=0.06')+ xlab('Day post-challenge')




### option 1
all_daily %>% filter(!(treatment %in% c('Zn+Cu', 'Bglu'))) %>%  ggplot(aes(x=time_point, y=mean_sal, fill=treatment, group=treatment)) +
  geom_jitter(aes(x=time_point, y=log_sal, color=treatment), data=filter(sal_data, !(treatment %in% c('Zn+Cu', 'Bglu'))& pignum !=101), alpha=.5) + 
  geom_line(aes(color=treatment), size=1.5)+
  geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) +   geom_point(shape=21, size=4) +
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ylab('log Salmonella') +
  xlab('Day post-challenge') +
  ggtitle('Daily shedding, group summary statistics')+ theme_bw()


# option 2
all_daily %>% filter(!(treatment %in% c('Zn+Cu', 'Bglu'))) %>%  ggplot(aes(x=time_point, y=mean_sal, fill=treatment, group=treatment)) +
  geom_line(aes(color=treatment), size=1.5)+
  geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) +   geom_point(shape=21, size=4) +
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ylab('log Salmonella') +
  xlab('Day post-challenge') +
  ggtitle('Daily shedding, group summary statistics')+ theme_bw()



all_daily %>% filter(treatment %in% c('control', 'RPS')) %>% ggplot(aes(x=time_point, y=mean_sal, fill=treatment, group=treatment)) +geom_line(aes(color=treatment))+ geom_errorbar(aes(ymin=mean_sal-se_sal,ymax=mean_sal+se_sal), width=.2) + 
  geom_point(shape=21, size=4) +scale_color_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ylab('log Salmonella') + xlab('Day post-challenge') +
  annotate(x=21, y=3, geom='text', label='p=0.02')+
  annotate(x=7, y=3.25, geom='text', label='p=0.03')+
  annotate(x=2, y=3.75, geom='text', label='p=0.06')+xlab('Day post-challenge')+ theme_bw()


#ZN_CONT_stats <- sal_data %>% filter(pignum != 101 & treatment %in% c('control', 'Zn+Cu'))


# ZN_CONT_stats

#pairwise.t.test(ZN_CONT_stats$log_sal,ZN_CONT_stats$treatXtime, p.adjust.method = 'none', var.equal=FALSE)
#### maybe some statistical test at each timepoint?


sum_sal <- sal_data %>% group_by(pignum) %>%
  summarise(AULC=trapz(time_point, log_sal),
            sum=sum(log_sal), 
            maxshed=log_sal[which.max(log_sal)], 
            day_max=time_point[which.max(log_sal)], 
            pos_samples=sum(log_sal > 0), 
            treatment=unique(treatment))

sum_sal %>% ggplot(aes(x=treatment, y=pos_samples, fill=treatment)) + geom_boxplot() + geom_jitter(shape=21, width=.2, height = .05)
sum_sal %>% ggplot(aes(x=treatment, y=day_max, fill=treatment)) + geom_boxplot() + geom_jitter(shape=21, width=.2, height = .05)
sum_sal %>% ggplot(aes(x=treatment, y=maxshed, fill=treatment)) + geom_boxplot() + geom_jitter(shape=21, width=.2, height = .05)

sum_sal %>% ggplot(aes(x=treatment, y=pos_samples, fill=treatment)) + geom_boxplot() + geom_jitter(shape=21, width=.2, height = .05)


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


# AULCs_av <- AULCs %>% group_by(day, treatment) %>% summarise(mean_AULC=mean(AULC), sd=sd(AULC), n=n(), sterr=sd/sqrt(n))
# 
# # AULCs_av$sd[is.na(AULCs_av$sd)] <- 0
# # AULCs_av$sterr[is.na(AULCs_av$sterr)] <- 0
# 
# AULCs_av$treatment <- factor(AULCs_av$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))
# 
# AULCs_av %>% filter(treatment %in% c('control', 'RPS')) %>% ggplot(aes(x=day, y=mean_AULC, group=treatment, color=treatment, fill = treatment)) +
#   geom_point() +
#   geom_line() + geom_ribbon(aes(ymin=mean_AULC-sterr, ymax=mean_AULC+sterr), alpha=.5) + 
#   scale_color_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Mean AULC over time') + 
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) 
# 
# 
# AULCs_av %>% filter(treatment %in% c('control', 'RPS')) %>% ggplot(aes(x=day, y=mean_AULC, group=treatment, color=treatment, fill = treatment)) +
#   geom_point(size=3) +
#   geom_line(size=1.5) + geom_errorbar(aes(ymin=mean_AULC-sterr, ymax=mean_AULC+sterr), size=1,width=1) + 
#   scale_color_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('') + 
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ylab('AULC') 
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

# what's this?
#########temp ##########

# re orders sum_sal
sum_sal <- sum_sal[match(filter(sal_data, time_point==2)$pignum,sum_sal$pignum),]



#
filter(sal_data, time_point==2)$pignum == filter(sal_data, time_point==7)$pignum
filter(sal_data, time_point==14)$pignum == filter(sal_data, time_point==21)$pignum
filter(sal_data, time_point==2)$pignum == filter(sal_data, time_point==21)$pignum
sum_sal$pignum == filter(sal_data, time_point==2)$pignum

# There is probably a better way to do this....
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

sum_sal$ADG <- gain/21
sum_sal$init_weight <- filter(sal_data, time_point == 0)$pig_weight
sum_sal$final_weight <- filter(sal_data, time_point == 21)$pig_weight
sum_sal$P_init_gain <- (gain/sum_sal$init_weight)*100
sum_sal$P_dry_matter0 <- filter(sal_data, time_point == 0)$percent_drymatter
sum_sal$P_dry_matter2 <- filter(sal_data, time_point == 2)$percent_drymatter
sum_sal$P_dry_matter7 <- filter(sal_data, time_point == 7)$percent_drymatter


# sum_sal <- merge(sum_sal, treats, by='pignum', all=TRUE)
sum_sal$treatment <- factor(sum_sal$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

hist(sum_sal$AULC, breaks = 50)



filter(sum_sal, pignum !=101) %>% ggplot(aes(x=treatment, y=AULC, fill=treatment))+
  geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(fill=treatment), shape=21, size=2, stroke=1.25, width = .12) + #geom_text(aes(label=pignum)) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Cumulative Salmonella shedding (AULC)', subtitle = 'Wilcoxon vs control: RPS p=0.013, Acid p=0.10')+ theme_bw()


filter(sum_sal, pignum !=101 & treatment %in% c('control', 'RPS')) %>% ggplot(aes(x=treatment, y=AULC, fill=treatment))+
  geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(fill=treatment), shape=21, size=2, stroke=1.25, width = .12) + #geom_text(aes(label=pignum)) +
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Cumulative Salmonella shedding (AULC)', subtitle = 'Wilcoxon p=0.013')+ theme_bw()





# filter(sum_sal, pignum !=101) %>% ggplot(aes(x=treatment, y=AULC, fill=treatment))+
#   geom_text(aes(label=pignum)) + #geom_jitter(aes(fill=treatment), shape=21, size=2, stroke=1.25, width = .12) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   ggtitle('Cumulative Salmonella shedding (AULC)', subtitle = 'Wilcoxon vs control: RPS p=0.013, Acid p=0.10, Bglu p=0.043')



# filter(sum_sal, pignum !=101) %>% select(pignum, AULC, treatment)

low_conts <- c(345,77,191,458,87,160)
high_conts <- c(337,122,224,419)

# sum_sal.filter <- filter(sum_sal, pignum !=101)

####### FOr keystone poster #######
# sum_sal.poster <- filter(sum_sal, pignum !=101 & treatment %in% c('control', 'RPS'))
# 
# sum_sal.poster %>% ggplot(aes(x=treatment, y=AULC, fill=treatment)) +
#   geom_boxplot() + geom_jitter(width = .2, shape=21, size=2, stroke=1) +
#   scale_fill_brewer(palette = 'Dark2') + ggtitle('Cumulative Salmonella shedding over 21 days', subtitle = 'Wilcox pvalue=0.013')




sum_sal_sum <- sum_sal %>% group_by(treatment) %>%
  summarise(mean_AULC=mean(AULC), sd=sd(AULC), num_obs=n(), se=sd/sqrt(n()))

# not this one.  Need boxplots because highlight low vs high shedders later

# ggplot(sum_sal_sum, aes(x=treatment, y=mean_AULC)) + geom_col(aes(fill=treatment)) +
#   geom_errorbar(aes(ymin=mean_AULC - se, ymax=mean_AULC + se), size=.5, width=.5) +
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   ggtitle('Mean AULC', subtitle = 'Error bars represent +/- 1 standard error from the mean')


wilcox.test(sum_sal$AULC[sum_sal$treatment == 'control'], sum_sal$AULC[sum_sal$treatment == 'RPS'])
t.test(sum_sal$AULC[sum_sal$treatment == 'control'], sum_sal$AULC[sum_sal$treatment == 'RPS'])

wilcox.test(sum_sal$AULC[sum_sal$treatment == 'control'], sum_sal$AULC[sum_sal$treatment == 'Acid'])
t.test(sum_sal$AULC[sum_sal$treatment == 'control'], sum_sal$AULC[sum_sal$treatment == 'Acid'])

# wilcox.test(sum_sal$AULC[sum_sal$treatment == 'control'], sum_sal$AULC[sum_sal$treatment == 'Bglu'])
# t.test(sum_sal$AULC[sum_sal$treatment == 'control'], sum_sal$AULC[sum_sal$treatment == 'Bglu'])

#### Trying a little mixed model stuff....

# library(lmerTest)
library(lme4)
# lmer(y ~ time * tx + (1 | subjects), data=data)
# lmer(y ~ time * tx + (time | subjects), data=data)

# filtering T0 because no pigs were shedding salmonella at t0
# sal_data2 <- sal_data %>% filter(pignum != 101 & time_point != 0 & treatment %in% c('control', 'Zn+Cu'))
sal_data2 <- sal_data %>% filter(time_point != 0)

sal_data2$pignum <- factor(sal_data2$pignum)

sal_data2 %>% ggplot(aes(x=time_point, y=log_sal, color=treatment)) + geom_point(alpha=.2) + geom_smooth(method = 'lm', fill=NA)

sal_data2 %>% ggplot(aes(x=time_point, y=log_sal, color=treatment)) +  geom_smooth(method = 'lm', fill=NA)

###

# LRT test?
# Are these two models equivalent?
shedding.null <- lme4::lmer(log_sal ~ time_point + (1|pignum) , data=sal_data2, REML = FALSE)
shedding.model <- lme4::lmer(log_sal ~ time_point * treatment + (1|pignum) , data=sal_data2, REML = FALSE)

anova(shedding.model, shedding.null)
anova(shedding.null, shedding.model)
# No, two models are not equivalent, adding treatment factor helps explain more variance in the data
# so now extract coeficcients for treatments? 
# what do these mean?  Are they the difference in intercept between reference group (control) and other groups?
# intercepts would be shedding at time_point 0?

###
sal_data2$time_point_fact
# fit <- lmer(log_sal ~ time_point_fact * treatment + (1|pignum) , data=sal_data2) # time is factor
fit <- lmer(log_sal ~ time_point * treatment + (1|pignum) , data=sal_data2)      # time is continuous



confints <- lme4::confint.merMod(fit)
fixefs <- lme4::fixef(fit)
colnames(confints) <- c('lower', 'upper')
confints <- confints[-(1:2),]
confints <- data.frame(confints)


confints$fixef <- fixefs
confints$coef <- rownames(confints)

confints %>% ggplot(aes(x=coef, y=fixef)) +geom_point()


whatsTHIS <- lme4::fortify.merMod(fit)
#hmmm...... interesteting
whatsTHIS %>% group_by(treatment, time_point) %>% summarise(avg=mean(log_sal))
whatsTHIS %>% ggplot(aes(x=time_point, y=.fitted)) +
  geom_point(aes(x=time_point, y=log_sal), color='red', alpha=.2) + geom_line(aes(x=time_point, y=log_sal, group=pignum), color='red', alpha=.2)+
  geom_point() +
  geom_line(aes(group=pignum)) +
  facet_wrap(~treatment) + ggtitle('fitted values in black, real values in red') + xlim(0,25)

# each pig has different intercetp, but slopes are the same within treatments


whatsTHIS %>% ggplot(aes(x=time_point, y=.fitted)) +
  geom_point(aes(x=time_point, y=log_sal), color='red', alpha=.2) + geom_line(aes(x=time_point, y=log_sal, group=pignum), color='red', alpha=.2)+
  geom_point() +
  geom_line(aes(group=pignum)) +
  facet_wrap(~treatment) + ggtitle('fitted values in black, real values in red') + xlim(0,25) +
  geom_hline(yintercept = 3.26) + # Intercept (y interept of control group at T=0?)
  geom_hline(yintercept = 3.26-1.25, color='blue') + geom_vline(xintercept = 0)

3.26-1.25




lme4::ranef(fit)
lme4::show(fit)
lme4::VarCorr(fit)
plot(fit)
13*7

summary(fit)
# summary(fit2)
# with lmerTest, fixed effects get pvalues, without, no pvalues.
# what do pvalues for fixed effects mean?

# how do you evaluate the quality of the model?
# is there an R2 thingy like for lm()?

library(lmerTest)
summ <- summary(fit)



#####
# use this #

D2_tuk <- sal_data %>% filter(time_point == 2 & !(treatment %in% c('Bglu', 'Zn+Cu')))
sal_data$treatment
D2_tuk$treatment

D2_aov <- aov(data=D2_tuk, log_sal ~ treatment)

summary(D2_aov)
TukeyHSD(D2_aov)


aov_AULC <- aov(data=sum_sal, AULC~treatment)
summary(aov_AULC)

TukeyHSD(aov_AULC)


# rpart::rpart()


cut(sum_sal$AULC, breaks = 3,include.lowest = TRUE, labels = c('low','mod', 'high'))
table(cut(sum_sal$AULC, breaks = 2,include.lowest = TRUE,  labels = c('low', 'high')))



#####
fit3 <- lm(data=sum_sal, AULC~treatment)
fit3.null <- lm(data=sum_sal, AULC~1)
# plot(fit3)

summary(fit3)
summary(fit3.null)
anova(fit3, fit3.null)

sum_sal %>% filter(pignum != 101) %>%  group_by(treatment) %>% summarise(tmean=mean(AULC)) %>% 
  mutate(d_cont=tmean-tmean[1])


conf_ints <- confint(fit3)
fit3
residuals.lm(fit3)

summary(fit3)
fit3_anov <- aov(fit3)
aov(fit3)


summary(fit3_anov)
TukeyHSD(fit3_anov, which='treatment')

# I'm going to need some help interpreting this...
# summ
# 
# conf_ints <- confint(fit)
# coefs <- summ$coefficients
# mains <- data.frame(cbind(coefs[3:7,], conf_ints[5:9,]))
# interacts <- data.frame(cbind(coefs[8:12,], conf_ints[10:14,]))
# 
# 
# colnames(mains) <- c('Estimate', 'std_err', 'df', 't_val', 'P_val', 'CI_L', 'CI_U')
# colnames(interacts) <- c('Estimate', 'std_err', 'df', 't_val', 'P_val', 'CI_L', 'CI_U')
# 
# 
# ### THESE ARE DUM ###
# mains %>%
#   rownames_to_column(var = 'treat') %>%
#   mutate(treatment=sub('treatment', '', treat)) %>% 
#   select(treatment, everything(), -treat) %>% 
#   gather(key = param, value = val, -treatment) %>% 
#   filter(param %in% c('CI_L', 'CI_U', 'Estimate', 'std_err')) %>%
#   spread(key = param, value = val) %>% 
#   ggplot(aes(x=treatment, y=Estimate, ymin=CI_L, ymax=CI_U)) + geom_point() + geom_errorbar() + 
#   geom_errorbar(aes(ymin=Estimate-std_err, ymax=Estimate+std_err), color='green')
# 
# interacts %>%
#   rownames_to_column(var = 'treat') %>%
#   mutate(treatment=sub('treatment', '', treat)) %>% 
#   select(treatment, everything(), -treat) %>% 
#   gather(key = param, value = val, -treatment) %>% 
#   filter(param %in% c('CI_L', 'CI_U', 'Estimate', 'std_err')) %>%
#   spread(key = param, value = val) %>% 
#   ggplot(aes(x=treatment, y=Estimate, ymin=CI_L, ymax=CI_U)) + geom_point() + geom_errorbar() + 
#   geom_errorbar(aes(ymin=Estimate-std_err, ymax=Estimate+std_err), color='green')
# 
# ### END DUM ###
# 
# anova(fit)
# rePCA(fit)
# 
# 
# sal_data2 %>% ggplot(aes(x=time_point, y=log_sal, group=treatment, color=treatment)) + geom_point() + geom_smooth(method = 'lm')



# really dont know if this makes any sense....
### early shed ###
# 
# sal_data2 <- sal_data %>% filter(pignum != 101 & time_point %in% c(2, 7))
# sal_data2$pignum <- factor(sal_data2$pignum)
# fit <- lmer(log_sal ~ time_point * treatment + (1|pignum), data=sal_data2)
# plot(fit)
# summ <- summary(fit)
# summ
# 
# ### late shed ###
# sal_data2 <- sal_data %>% filter(pignum != 101 & time_point %in% c(14, 21))
# sal_data2$pignum <- factor(sal_data2$pignum)
# fit <- lmer(log_sal ~ time_point * treatment + (1|pignum), data=sal_data2)
# plot(fit)
# summ <- summary(fit)
# summ




######

# sal_data %>% filter(time_point < 14) %>% ggplot(aes(x=treatment, y=temp, fill=treatment, group=treatXtime)) +
#   geom_boxplot() + geom_jitter(shape=21, size=2)+ geom_text(aes(label=pignum)) +
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   ggtitle('Rectal Temperature') + facet_wrap(~time_point)
# 
# sal_data %>% filter(time_point < 14) %>% ggplot(aes(x=time_point_fact, y=temp, fill=treatment, group=treatXtime)) +
#   geom_boxplot() + geom_jitter(shape=21, size=2)+
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
#   ggtitle('Rectal temperatures') + facet_wrap(~treatment)#+ geom_text(aes(label=pignum))
# 
#### oooh.... maybe can do an AULC for temperature


### weight figs

# this shows pig 101 not gaining as the other pigs are #
# use this at beginning
# sal_data %>% filter(time_point %in% c(0,14)) %>%
#   ggplot(aes(x=time_point, y=pig_weight, group=pignum, color=treatment)) + geom_line() +
#   scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
#   geom_line(data = filter(sal_data, pignum ==101 & time_point %in% c(0,14)), color='black', size=1.25) +
#   ggtitle('weight gain over the 1st two weeks', subtitle = 'Pig 101 shown in black') + ylim(0,32)
# 

sal_data %>% filter(pignum !=101) %>% ggplot(aes(x=time_point_fact, y=pig_weight, fill=treatment, group=treatXtime)) +
  geom_boxplot() +scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  ggtitle('Weight') #+ geom_text(aes(label=pignum))


# sal_data %>% filter(pignum !=101 & treatment %in% c('control', 'Bglu')) %>% ggplot(aes(x=time_point_fact, y=pig_weight, fill=treatment, group=treatXtime)) +
#   geom_boxplot() +scale_fill_manual(values=c('#33CC33', 'purple', 'orange', 'red', 'grey', 'purple'))+
#   ggtitle('Weight') #+ geom_text(aes(label=pignum))
# 
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

tis <- tis %>% filter(!treatment %in% c('Zn+Cu', 'Bglu'))

tis$treatment <- factor(tis$treatment, levels = c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))


pairwise.wilcox.test(tis$log_sal[tis$tissue == 'cecal_cont'], tis$treatment[tis$tissue == 'cecal_cont'], p.adjust.method = 'none')

pairwise.t.test(tis$log_sal[tis$tissue == 'cecal_cont'], tis$treatment[tis$tissue == 'cecal_cont'], p.adjust.method = 'none')
pairwise.t.test(tis$log_sal[tis$tissue == 'Cecum'], tis$treatment[tis$tissue == 'Cecum'], p.adjust.method = 'none')




tis %>% filter(tissue=='cecal_cont') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) + geom_boxplot(outlier.alpha = 0)+ geom_jitter(shape=21,width = .1, size=2.25)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Cecal Contents')


tis %>% #filter(tissue=='cecal_cont') %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape=21,width = .2, size=2.25) +
  facet_wrap(~tissue) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Salmonella colonization at D21, tissues')

tis %>% filter(treatment %in% c('control', 'RPS')) %>%
  ggplot(aes(x=treatment, y=log_sal, group=treatment, fill=treatment)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape=21,width = .2, size=2.25) +
  facet_wrap(~tissue, nrow = 1) +
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Salmonella colonization at D21, tissues') + theme_bw()






tis_RPS <- tis %>% filter(treatment =='RPS')
sum_sal_RPS <- sum_sal %>% filter(treatment == 'RPS')
tis_RPS$shed <- ifelse(tis_RPS$pignum %in% c(373,321,181,392,97), 'low', 'high')
sum_sal_RPS$shed <- ifelse(sum_sal_RPS$pignum %in% c(373,321,181,392,97), 'low', 'high')

tis_RPS %>% ggplot(aes(x=shed, y=log_sal, group=shed, fill=shed)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape=21,width = .1, size=2.25) +
  facet_wrap(~tissue) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('Cecal Contents')

sum_sal_RPS %>% ggplot(aes(x=shed, y=AULC, group=shed, fill=shed)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape=21,width = .1, size=2.25) +
  #facet_wrap(~tissue) +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) +
  ggtitle('AULC')




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
# tis.sp <- tis.sp[,2]
tis.sp <- tis.sp[,-2]
sal_for_cor <- merge(sum_sal, tis.sp, by = 'pignum')

### NEED CECAL VFA DATA FOR THIS NEXT SECTION ###

vfas_for_cor <- read_csv('data/FS12b_vfas_for_cor.csv')

### NEED TO REMOVE ONE TREATMENT BEFORE MERGE
vfas_for_cor <- vfas_for_cor %>% select(-treatment)
sal_for_cor <- merge(sal_for_cor, vfas_for_cor, by = 'pignum')
rownames(sal_for_cor) <- sal_for_cor$pignum


##### fecal corrs #####
fec_cor <- res.all %>% select(-treatment) %>% full_join(sum_sal, by = 'pignum')
fec_cor <- fec_cor[-288,]

sum_sal

fec_cor %>% filter(time == 0) %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))
fec_cor %>% filter(time == 0) %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))
fec_cor %>% filter(time == 0) %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))


fec_cor %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + facet_wrap(~time)
fec_cor %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + facet_wrap(~time)



fec_cor %>% filter(time == 0) %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')+ scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))


### is this meta from the 16S stuff?
meta %>% filter(experiment == 'X12b') %>% ggplot(aes(x=caproate, y=log_sal)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~day)

meta_for_corr <- meta %>% mutate(day_fact=factor(day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))) %>% filter(experiment == 'X12b')

meta_for_corr %>% filter(tissue =='F' & day != 'D0') %>% ggplot(aes(x=caproate, y=log_sal)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~day_fact)
meta_for_corr %>% filter(tissue =='F' & day != 'D0') %>% ggplot(aes(x=valerate, y=log_sal)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~day_fact)
meta_for_corr %>% filter(tissue =='F' & day != 'D0') %>% ggplot(aes(x=butyrate, y=log_sal)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~day_fact)

# meta_for_corr$treatment
sal_vfa_cor <- meta_for_corr %>% select(pignum,day, tissue, treatment, log_sal, AULC, ends_with('ate')) %>%
  filter(pignum!=101 & tissue != 'Q' & treatment %in% c('Control', 'RPS', 'Acid', 'RCS')) %>% na.omit()

sal_vfa_cor <- sal_vfa_cor %>% mutate(total=rowSums(.[grep("ate", names(.))]))


globa_log_sal_VFA_corrs <- sal_vfa_cor %>% group_by(day, tissue) %>% 
                                          summarise(ace_corP=cor.test(log_sal, acetate)$p.value, 
                                                    pro_corP=cor.test(log_sal, propionate)$p.value, 
                                                    but_corP=cor.test(log_sal, butyrate)$p.value, 
                                                    val_corP=cor.test(log_sal, valerate)$p.value, 
                                                    cap_corP=cor.test(log_sal, caproate)$p.value, 
                                                    isob_corP=cor.test(log_sal, isobutyrate)$p.value, 
                                                    isov_corP=cor.test(log_sal, isovalerate)$p.value, 
                                                    tot_corP=cor.test(log_sal, total)$p.value)%>% na.omit() %>%
  gather(-(day:tissue), key = 'vfa', value = 'pval') %>% filter(pval < 0.05)



globa_AULC_VFA_corrs <- sal_vfa_cor %>% group_by(day, tissue) %>% 
  summarise(ace_corP=cor.test(AULC, acetate, method = 'spearman')$p.value, 
            pro_corP=cor.test(AULC, propionate, method = 'spearman')$p.value, 
            but_corP=cor.test(AULC, butyrate, method = 'spearman')$p.value, 
            val_corP=cor.test(AULC, valerate, method = 'spearman')$p.value, 
            cap_corP=cor.test(AULC, caproate, method = 'spearman')$p.value, 
            isob_corP=cor.test(AULC, isobutyrate, method = 'spearman')$p.value, 
            isov_corP=cor.test(AULC, isovalerate, method = 'spearman')$p.value, 
            tot_corP=cor.test(AULC, total)$p.value, method = 'spearman')%>% na.omit() %>%
  gather(-(day:tissue), key = 'vfa', value = 'pval') %>% filter(pval < 0.05)

 c('acetate', 'propionate', 'valerate', 'caproate', 'total')
 
globa_AULC_VFA_corrs %>% filter(day =='D0')
globa_AULC_VFA_corrs %>% filter(day =='D2')


treat_log_sal_VFA_corrs <- sal_vfa_cor %>% group_by(day, tissue, treatment) %>% 
                                          summarise(ace_corP=cor.test(log_sal, acetate, method = 'spearman')$p.value, 
                                                    pro_corP=cor.test(log_sal, propionate, method = 'spearman')$p.value, 
                                                    but_corP=cor.test(log_sal, butyrate, method = 'spearman')$p.value, 
                                                    val_corP=cor.test(log_sal, valerate, method = 'spearman')$p.value, 
                                                    cap_corP=cor.test(log_sal, caproate, method = 'spearman')$p.value, 
                                                    isob_corP=cor.test(log_sal, isobutyrate, method = 'spearman')$p.value, 
                                                    isov_corP=cor.test(log_sal, isovalerate, method = 'spearman')$p.value, 
                                                    tot_corP=cor.test(log_sal, total, method = 'spearman')$p.value) %>% na.omit() %>%
  gather(-(day:treatment), key = 'vfa', value = 'pval') %>% filter(pval < 0.15)


treat_AULC_VFA_corrs <- sal_vfa_cor %>% group_by(day, tissue, treatment) %>% 
  summarise(ace_corP=cor.test(AULC, acetate, method = 'spearman')$p.value, 
            pro_corP=cor.test(AULC, propionate, method = 'spearman')$p.value, 
            but_corP=cor.test(AULC, butyrate, method = 'spearman')$p.value, 
            val_corP=cor.test(AULC, valerate, method = 'spearman')$p.value, 
            cap_corP=cor.test(AULC, caproate, method = 'spearman')$p.value, 
            isob_corP=cor.test(AULC, isobutyrate, method = 'spearman')$p.value, 
            isov_corP=cor.test(AULC, isovalerate, method = 'spearman')$p.value, 
            tot_corP=cor.test(AULC, total, method = 'spearman')$p.value) %>% na.omit() %>%
  gather(-(day:treatment), key = 'vfa', value = 'pval') %>% filter(pval < 0.15)

##### schemezone #####
globa_AULC_VFA_corrs %>% filter(pval < 0.05 & day == 'D0')
treat_AULC_VFA_corrs %>% filter(pval < 0.15 & day == 'D0')

globa_AULC_VFA_corrs %>% filter(pval < 0.05 & day=='D21') %>% print(n=40)
treat_AULC_VFA_corrs %>% filter(pval < 0.15 & day=='D21') %>% print(n=40)

treat_log_sal_VFA_corrs %>%filter(pval <0.1) %>%  print(n = 40)
globa_log_sal_VFA_corrs

#### AULC CORRELATIONS PLOTS #####
sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D0') %>% 
  ggplot(aes(x=mM, y=AULC)) +
  geom_smooth(color='black', method = 'lm', fill=NA) +geom_point(aes(color=treatment), alpha=.5)+
  geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
  facet_wrap(~VFA, scales = 'free') + ggtitle('D0 fecal VFAs correlate with final AULC')

sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D21' & tissue =='C') %>% 
  ggplot(aes(x=mM, y=AULC)) +
  geom_smooth(color='black', method = 'lm', fill=NA) +geom_point(aes(color=treatment), alpha=.5)+
  geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
  facet_wrap(~VFA, scales = 'free') + ggtitle('D21 cecal VFAs correlate with final AULC')

#### LOG SAL CORRELATIONS PLOTS #########
sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D2') %>% 
  ggplot(aes(x=mM, y=log_sal)) +
  geom_smooth(color='black', method = 'lm', fill=NA) +
  geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
  facet_wrap(~VFA, scales = 'free') + ggtitle('D2 fecal VFAs correlate with log_sal')

sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D7') %>% 
  ggplot(aes(x=mM, y=log_sal)) +
  geom_smooth(color='black', method = 'lm', fill=NA) +
  geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
  facet_wrap(~VFA, scales = 'free') + ggtitle('D7 fecal VFAs correlate with log_sal')

sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D14') %>% 
  ggplot(aes(x=mM, y=log_sal)) +
  geom_smooth(color='black', method = 'lm', fill=NA) +
  geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
  facet_wrap(~VFA, scales = 'free') + ggtitle('D14 fecal VFAs correlate with log_sal')

sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D21' & tissue =='F') %>% 
  ggplot(aes(x=mM, y=log_sal)) +
  geom_smooth(color='black', method = 'lm', fill=NA) +
  geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
  facet_wrap(~VFA, scales = 'free') + ggtitle('D21 fecal VFAs correlate with log_sal')

sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D21' & tissue =='X') %>% 
  ggplot(aes(x=mM, y=log_sal)) +
  geom_smooth(color='black', method = 'lm', fill=NA) +
  geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
  facet_wrap(~VFA, scales = 'free') + ggtitle('D21 cecal VFAs correlate with cecal_tissue log_sal')

sal_vfa_cor %>% gather(-(pignum:AULC), key=VFA, value=mM) %>% filter(day == 'D21' & tissue =='C') %>% 
  ggplot(aes(x=mM, y=log_sal)) +
  geom_smooth(color='black', method = 'lm', fill=NA) +
  geom_smooth(aes(color=treatment), method = 'lm', fill=NA)+
  facet_wrap(~VFA, scales = 'free') + ggtitle('D21 cecal VFAs correlate with cecal log_sal')



########
treat_log_sal_VFA_corrs %>% na.omit()

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



sal_for_cor %>% filter(treatment %in% c('control','RPS')) %>% ggplot(aes(x=butyrate, y=AULC)) +
  geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme_bw()+ ggtitle('Correlation between AULC and cecal butyrate at D21') + xlab('butyrate (mM)') 


sal_for_cor %>% filter(treatment == 'Bglu') %>% ggplot(aes(x=butyrate, y=AULC)) +
  geom_smooth(method = 'lm') + geom_point(aes(fill=treatment), shape = 21, size=2) + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + ggtitle('') + xlab('butyrate (mM)')


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

# filter(sal_for_cor, treatment %in% c('control', 'RPS'))
# filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')


# filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=butyrate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# 


# filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=valerate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# 

# filter(sal_for_cor, treatment =='control') %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')
# filter(sal_for_cor, treatment =='RPS') %>% ggplot(aes(x=caproate, y=AULC)) + geom_point() + geom_smooth(method = 'lm')


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




### TO DO !!!!!!!! ########

###CORRELATIONS WITH FECAL SCFAS ###



# filter(sal_for_cor, treatment =='RPS') %>% cor.test()

# poster <- sal_for_cor %>% select(-starts_with('P_')) %>% 
#   gather(variable, response, -pignum, -treatment) %>% 
#   filter(treatment %in% c('control', 'RPS')) %>% 
#   filter(variable %in% c('AULC', 'butyrate', 'valerate', 'caproate', 'total')) %>% spread(variable, response) %>% 
#   gather(variable, response, -pignum, -treatment, -AULC)
# poster$response <- poster$response *3
# 
# 
# poster$variable <- factor(poster$variable, levels = c('butyrate', 'valerate', 'caproate', 'total'))
# ggplot(poster, aes(x=response, y=AULC, fill=treatment)) +
#   geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + facet_wrap(~variable, scales = 'free') +
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   geom_point(shape=21, stroke=1, size=3) + ggtitle('Cecal SCFA concentrations VS cumulative Salmonella shedding') + xlab('Concentration (mM)')
# 
# 
# poster %>% filter(variable == 'butyrate') %>% ggplot(aes(x=response, y=AULC, fill=treatment)) +
#   geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + #facet_wrap(~variable, scales = 'free') +
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   geom_point(shape=21, stroke=1, size=3) + ggtitle('') + xlab('Concentration (mM)')
# #3399FF
# 
# 
# poster %>% filter(variable == 'caproate') %>% ggplot(aes(x=response, y=AULC, fill=treatment)) +
#   geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + #facet_wrap(~variable, scales = 'free') +
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   geom_point(shape=21, stroke=1, size=3) + ggtitle('') + xlab('Concentration (mM)')
# 
# 
# 
# poster %>% filter(variable == 'valerate') %>% ggplot(aes(x=response, y=AULC, fill=treatment)) +
#   geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + #facet_wrap(~variable, scales = 'free') +
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   geom_point(shape=21, stroke=1, size=3) + ggtitle('') + xlab('Concentration (mM)')
# 
# 
# 
# poster %>% filter(variable == 'total') %>% ggplot(aes(x=response, y=AULC, fill=treatment)) +
#   geom_smooth(aes(x=response, y=AULC), inherit.aes = FALSE,method = 'lm', fullrange = TRUE, color='grey27') + #facet_wrap(~variable, scales = 'free') +
#   scale_fill_manual(values=c('red', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
#   geom_point(shape=21, stroke=1, size=3) + ggtitle('') + xlab('Concentration (mM)')
# 




# cecal_cont
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
# 
# salcor <- rcorr(sal_for_cor.m)
# salcor$P
# apply(sal_for_cor, MARGIN = 2, str_trim) %>% apply(MARGIN = 2, as.numeric) %>% as.matrix()
# str_trim()
# 
# salcor.sigs <- rcorr_to_ggnet(salcor, spearcut = -1, pcut = 0.05)
# 
# library(funfuns)
# library(geomnet)
# ################
# 
# nodes <- rbind(gather_nodes(sal_for_cor.m, 'Sal'))
# 
# nodes$node <- as.character(nodes$node)
# #nodes$node[770] <- 'y'
# alll <- na.exclude(fortify(as.edgedf(salcor.sigs), nodes))
# na.exclude(alll)
# alll$spear
# rep(1, length(alll$spear))
# test <- data.frame(rep(1, length(alll$spear)),alll$spear)
# colnames(test) <- c('x', 'spear')
# 
# filts_pos <- filter(alll, spear > 0)
# filts_neg <- filter(alll, spear < 0)
# 
# # butyrate caproate valerate neg cor w/ AULC
# 
# 
# 
# 
# 
# ggplot()
# 
# tetetete <- ggplot(test, aes(x=x, y=spear, color=spear)) + geom_point() + scale_color_gradient2()
# 
# hmm <- ggplot_build(tetetete)
# 
# hmm$data[[1]]$colour
# 
# 
# 
# salnet <- ggplot(alll, aes(from_id = from_id, to_id = to_id, label=from, color=type, ecolour=spear)) +
#   geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)),
#            aes(color = type, label = from_id),
#            linewidth = 2, size = 5, vjust = 0, alpha = 0.3,
#            repel = TRUE, fontsize=3, singletons = FALSE,labelcolour="black",
#            labelgeom = 'text', ecolour = hmm$data[[1]]$colour) +
#   theme_net()
# 
# salnet
# 
# ##################
# 
# # Compare the ADG in different intervals back to Salmonella colonization in different tissues
# # prelavence of Salmonella in farms tend to have lower ADG.  Can we see some of this in our pigs here?
# 
# # can we see some
# 
# # Do larger pigs tend to have greater amounts of salmonella
# 
# 
# 
# 
# 
# ##### some correlation stuff here ####
# 
# 
# 
# # prob need to replace some of this stuff #
# 
# 
# ### change in rectal temp? ###
# 
# filter(sal_data, time_point ==0)$temp
# 
# 
# delt_temp_2 <- filter(sal_data, time_point ==2)$temp - filter(sal_data, time_point ==0)$temp
# delt_temp_7 <- filter(sal_data, time_point ==7)$temp - filter(sal_data, time_point ==0)$temp
# 
# sal_data$delta_temp <- NA
# sal_data[sal_data$time_point ==2,]$delta_temp <- delt_temp_2
# sal_data[sal_data$time_point ==7,]$delta_temp <- delt_temp_7
# 
# 
# 
# tis$PA <- ifelse(tis$Salmonella > 0, 1, 0)
# 
# 
# tis
# 
# 
# 
# 
# fisher_map <- function(data, row_vector, col_vector){
#   fisher.test(data[row_vector, col_vector])$p.value
# }
# 
# 
# PA_tab <- tis %>% select(treatment, tissue, PA) %>%
#   group_by(tissue, treatment) %>% 
#   summarise(pos = sum(PA == 1), neg = sum(PA == 0)) %>% group_by(tissue) %>% nest() %>% 
#   mutate(RPSp = map(data, fisher_map, c(1,2), c(2,3)), 
#          Acidp = map(data, fisher_map, c(1,3), c(2,3)),
#          Znp = map(data, fisher_map, c(1,4), c(2,3)), 
#          RCSp = map(data, fisher_map, c(1,5), c(2,3)), 
#          Bglup = map(data, fisher_map, c(1,6), c(2,3))) %>% select(-data) %>% unnest()
# 
# ftes <- fisher.test(PA_tab[[2]][[1]][c(1,3),-1])
# 
# ftes$data.name
# 
# row.
# cecconttest <- PA_tab[[2]][[1]]
# 
# 
# muc <- data.frame(SB1434 = c(1050, 0), 
#            SX240 = c(1832,5), row.names = c('normal', 'mucoid'))
# 
# 
# intra <- data.frame(norm=c(1175,0), 
#                     mono=c(1832, 5), 
#                     row.names = c('normal', 'mucoid'))
# 
# t(muc)
# fisher.test(muc)
# fisher.test(t(muc))
# 
# t(intra)
# fisher.test(intra)
# fisher.test(t(intra))
# 
# 
