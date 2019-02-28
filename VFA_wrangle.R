getwd()
# setwd('./FS12/FS12b/')

library(tidyverse)

list.files()

fec1 <- read.csv('./data/FS12_fec1_results.csv')
fec2 <- read.csv('./data/FS12_fec2_results.csv')
fec3 <- read.csv('./data/FS12_fec3_results.csv')
key1 <- read.csv('./data/GCKey_fec_1.csv')
key2 <- read.csv('./data/GCKey_fec_2.csv')
key3 <- read.csv('./data/GCKey_fec_3.csv')



colnames(fec1)[2] <- 'position'
colnames(fec2)[2] <- 'position'
colnames(fec3)[2] <- 'position'



res1 <- merge(key1, fec1, by = 'position', all = TRUE)
res2 <- merge(key2, fec2, by = 'position', all = TRUE)
res3 <- merge(key3, fec3, by = 'position', all = TRUE)
colnames(res1) == colnames(res2)



res.all <- rbind(res1, res2, res3)
res.all[is.na(res.all)] <- 0
res.all[,-c(1:4)] <- res.all[,-c(1:4)] *3 # corrects for sample dillution



res.all[,-c(1:4)] - colMeans(res.all[grep('BLANK', res.all$sample),-c(1:4)])

res.all[,-c(1:4)] <- scale(x = res.all[,-c(1:4)],
                           center = colMeans(res.all[grep('BLANK', res.all$sample),-c(1:4)]),
                           scale = FALSE)

res.all[,-c(1:4)][res.all[,-c(1:4)] < 0] <- 0

res.all <- res.all[-grep('BLANK', res.all$sample),]

res.all <- res.all[,-c(1,4,15)]

colnames(res.all)[1] <- 'pignum'

res.all[res.all$pignum ==329,]$pignum <- 392
res.all[res.all$pignum == 129,]$pignum <- 126


res.all$lactate <- res.all$lactate1 + res.all$lactate2

rooms <- read.csv('./data/Rooms.csv')

res.all$treatment <- ifelse(res.all$pignum %in% rooms$X6, 'control',
                             ifelse(res.all$pignum %in% rooms$X7, 'RPS', 
                                    ifelse(res.all$pignum %in% rooms$X8, 'Acid', 
                                           ifelse(res.all$pignum %in% rooms$X9, 'Zn+Cu',
                                                  ifelse(res.all$pignum %in% rooms$X10, 'RCS',
                                                         ifelse(res.all$pignum %in% rooms$X11, 'Bglu', 'asdfsa'))))))

res.all$set <- paste(res.all$time, res.all$treatment, sep = '_')
res.all <- res.all %>% select(-lactate1, -lactate2) %>% select(pignum, time, treatment, set, everything()) %>% filter(pignum !=101)

res.all$total <- rowSums(select(.data = res.all, -pignum, -time, -treatment, -set))

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
res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=butyrate, fill=treatment)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + geom_jitter(shape=21, stroke=1.2, width = .2, size=2) + ggtitle("D0 Fecal Butyrate")
res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=caproate, fill=treatment)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+ geom_jitter(shape=21, stroke=1.2, width = .2, size=2) + ggtitle("D0 Fecal Caproate")
res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=valerate, fill=treatment)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+ geom_jitter(shape=21, stroke=1.2, width = .2, size=2) + ggtitle("D0 Fecal Valerate")


res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=butyrate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=butyrate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=butyrate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=butyrate)) + geom_boxplot()+ geom_text(aes(label=pignum))


res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=caproate)) + geom_boxplot()+ geom_text(aes(label=pignum))

res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()+ geom_text(aes(label=pignum))

res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=acetate)) + geom_boxplot()


res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()

res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=propionate)) + geom_boxplot()
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=valerate)) + geom_boxplot()


res.all %>% filter(time == 0) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot() + geom_text(aes(label=pignum))
res.all %>% filter(time == 2) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 7) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 14) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot()+ geom_text(aes(label=pignum))
res.all %>% filter(time == 21) %>% ggplot(aes(x=treatment, y=total)) + geom_boxplot()+ geom_text(aes(label=pignum))




res.all[res.all$treatment == 'asdfsa',]

res.all %>% filter(time == 2)
126 %in% res.all[res.all$time == 0, ]$pignum
126 %in% res.all[res.all$time == 2, ]$pignum
126 %in% res.all[res.all$time == 7, ]$pignum
126 %in% res.all[res.all$time == 14, ]$pignum
126 %in% res.all[res.all$time == 21, ]$pignum


326 %in% res.all[res.all$time == 0, ]$pignum
326 %in% res.all[res.all$time == 2, ]$pignum
326 %in% res.all[res.all$time == 7, ]$pignum
326 %in% res.all[res.all$time == 14, ]$pignum
326 %in% res.all[res.all$time == 21, ]$pignum

219 %in% res.all[res.all$time == 0, ]$pignum
219 %in% res.all[res.all$time == 2, ]$pignum
219 %in% res.all[res.all$time == 7, ]$pignum
219 %in% res.all[res.all$time == 14, ]$pignum
219 %in% res.all[res.all$time == 21, ]$pignum


392 %in% res.all[res.all$time == 0, ]$pignum
392 %in% res.all[res.all$time == 2, ]$pignum
392 %in% res.all[res.all$time == 7, ]$pignum
392 %in% res.all[res.all$time == 14, ]$pignum
392 %in% res.all[res.all$time == 21, ]$pignum

library(vegan)


rownames(res.all) <- paste(res.all$pignum, res.all$time, res.all$treatment, sep = '_')

vegdist(res.all[,-c(1:4)])

# library(devtools)
# install_github('Jtrachsel/funfuns')
library(funfuns)

res.prop <- res.all
res.prop[,-c(1:4,18)] <- res.all[,-c(1:4,18)]/rowSums(res.all[,-c(1:4,18)])

res.prop <- res.prop %>% select(-total)

res.prop %>% ggplot(aes(x=treatment, y=acetate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=propionate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=butyrate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=valerate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=caproate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=isovalerate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=oxalate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=phenylacetate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=succinate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=fumarate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()
res.prop %>% ggplot(aes(x=treatment, y=lactate, fill=treatment)) + facet_wrap(~time) + geom_boxplot()

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

res.gather %>% filter(variable == 'total')%>% ggplot(aes(x=time, y=concentration, group=pignum, color=treatment)) + geom_line()

res.sum <- res.gather %>%  group_by(time, treatment, variable) %>%
  summarise(mean=mean(concentration), sd=sd(concentration), n=n(), sterr=sd/sqrt(n))

res.sum <- res.sum %>% mutate(treat_var=paste(treatment, variable, sep='_'))

res.sum$treatment <- factor(res.sum$treatment, levels= c('control', 'RPS', 'Acid', 'Zn+Cu', 'RCS', 'Bglu'))

res.sum %>% ggplot(aes(x=time, y=mean, group=treat_var, color=treatment)) +
  geom_line() + facet_wrap(~variable, scales='free') + scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Mean concentration of fecal SCFAs over time')

SCFA.j <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'jaccard')
SCFA.j <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'time', distance_method = 'jaccard')

SCFA.b <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'bray')
SCFA.e <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'euclidean')
SCFA.man <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'manhattan')
SCFA.g <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'gower')
SCFA.can <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'canberra')
SCFA.kul <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'kulczynski')
#SCFA.mor <- NMDS_ellipse(OTU_table = res.all[,-c(1:4)], metadata = res.all[,c(1:4)], grouping_set = 'set', distance_method = 'morisita')

SCFA.j[[1]]$time <- factor(SCFA.j[[1]]$time)


ggplot(SCFA.j[[1]], aes(x=MDS1, y=MDS2, color=treatment)) + geom_point()
ggplot(SCFA.j[[1]], aes(x=MDS1, y=MDS2, color=time)) + geom_point() + geom_text(aes(label=pignum))


SCFA.j[[2]]$treatment <- gsub('([0-9]+)_([A-Za-z]+)','\\2',SCFA.j[[2]]$group)
SCFA.j[[2]]$time <- gsub('([0-9]+)_(.*)','\\1',SCFA.j[[2]]$group)
ggplot(SCFA.j[[1]], aes(x=MDS1, y=MDS2, color=time)) + geom_point() +
  geom_path(data = SCFA.j[[2]], aes(x=NMDS1, y=NMDS2, group=group))




