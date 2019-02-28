# setwd('~/FS12/vfas/')
library(tidyverse)
# 
# vfas1 <- read.csv('FS12_nec1_res.csv')
# vfas2 <- read.csv('FS12_nec2_results.csv')


vfas1 <- read.csv('data/FS12a_nec1_res.csv')
vfas2 <- read.csv('data/FS12a_nec2_res.csv')

# GCkey1 <- read.csv('GCKey1.csv')
# GCkey2 <- read.csv('GCKey2.csv')

GCkey1 <- read.csv('data/FS12a_GCKey1.csv')
GCkey2 <- read.csv('data/FS12a_GCKey2.csv')
#colnames(GCkey)[3] <- 'line'

vfas1 <- merge(GCkey1, vfas1, by='line')
vfas2 <- merge(GCkey2, vfas2, by = 'line')

vfas <- rbind(vfas1, vfas2)


vfas[is.na(vfas)] <- 0
vfas[,-c(1:5)] <- vfas[,-c(1:5)] * 3 # correcting for sample dillution during VFA prep
vfas$total <- rowSums(vfas[,-c(1:5)])
vfas$lactate <- vfas$lactate1 + vfas$lactate2

#vfas$butprop <- vfas$butyrate/vfas$propionate
#vfas <- merge(GCkey, vfas, by = 'line')

vfas <- vfas[,-c(1,5,11,12,16)]
vfas$group <- paste(vfas$pen, vfas$timepoint, sep='_')
rownames(vfas) <- vfas$group

colnames(vfas)

meta <- vfas[,c(1,2,3,18)]
vfas.matrix <- vfas[,-c(1,2,3,16,18)]
#vfas <- vfas[,-17]
vfas <- vfas %>% select(pen:treatment, group, everything(), -formate, -fumarate)

# write_csv(vfas, 'FS12a_vfas.csv')
vfas <- read_csv('FS12a_vfas.csv')


vfas.gather <- gather(vfas, key = vfa, value = concentration, -(pen:treatment), -group)
#vfas.gather$treatment <- factor(vfas.gather$treatment, levels = )
vfas.gather$treatment <- factor(vfas.gather$treatment)
vfas.gather$concentration <- as.numeric(vfas.gather$concentration)
vfas.gather$vfa <- factor(vfas.gather$vfa, levels = c('acetate', 'propionate', 'butyrate','valerate', 'caproate', 'isobutyrate', 'isovalerate', 'total', 'formate', 'fumarate', 'lactate', 'oxalate', 'phenylacetate', 'succinate'))
levels(vfas.gather$vfa)
vfas.gather$treatment_timepoint <- paste(vfas.gather$treatment, vfas.gather$timepoint, sep = '_')
vfas.gather$treatment_timepoint <- factor(vfas.gather$treatment_timepoint, 
                                          levels = c('1_0','2_0','3_0','6_0','8_0','9_0','10_0',
                                                     '1_24','2_24','3_24','6_24','8_24','9_24','10_24'))


vfas.gather$treatment_timepoint
vfas.gather$timepoint <- factor(vfas.gather$timepoint)

vfas.gather %>% filter(timepoint ==0) %>% filter(!(vfa %in% c('formate', 'fumarate', 'lactate', 'oxalate', 'phenylacetate', 'succinate'))) %>% 
  ggplot(aes(x=treatment, y=concentration, group=treatment, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~vfa, scales = 'free') + ggtitle('Cecal VFAs, 4 weeks post weaning', subtitle = )

vfas.gather %>% filter(timepoint ==24) %>% filter(!(vfa %in% c('formate', 'fumarate', 'lactate', 'oxalate', 'phenylacetate', 'succinate'))) %>% 
  ggplot(aes(x=treatment, y=concentration, group=treatment, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~vfa, scales = 'free') + ggtitle('Cecal VFAs, 4 weeks post-weaning.  24 hour incubation', subtitle = 'passed into anaerobic chamber and allowed to ferment for 24 hours')

vfas.gather %>% filter(!(vfa %in% c('formate', 'fumarate', 'lactate', 'oxalate', 'phenylacetate', 'succinate'))) %>% 
  ggplot(aes(x=timepoint, y=concentration, fill=treatment, group=treatment_timepoint)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~vfa, scales = 'free') + ggtitle('Cecal VFAs, 4 weeks post-weaning, with and without incubation')

#### succinate and total stuff ####

vfas.gather %>% filter(timepoint ==0) %>% filter((vfa %in% c('succinate', 'total') & treatment %in% c(1,8))) %>% 
  ggplot(aes(x=treatment, y=concentration, group=treatment, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~vfa, scales = 'free') + ggtitle('Cecal VFAs, 4 weeks post weaning', subtitle = )

vfas.gather %>% filter(timepoint ==24) %>% filter((vfa %in% c('succinate', 'total') & treatment %in% c(1,8))) %>% 
  ggplot(aes(x=treatment, y=concentration, group=treatment, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~vfa, scales = 'free') + ggtitle('Cecal VFAs, 4 weeks post-weaning.  24 hour incubation', subtitle = 'passed into anaerobic chamber and allowed to ferment for 24 hours')



vfas0 <- vfas[vfas$timepoint == 0,]


but0 <- pairwise.wilcox.test(vfas0$butyrate, vfas0$treatment, p.adjust.method = 'none')
tot0 <- pairwise.wilcox.test(vfas0$total, vfas0$treatment, p.adjust.method = 'none')

write.csv(but0$p.value, file = 'but0.csv')
write.csv(tot0$p.value, file = 'tot0.csv')

vfas24 <- vfas[vfas$timepoint == 24,]

but24 <- pairwise.wilcox.test(vfas24$butyrate, vfas24$treatment, p.adjust.method = 'none')
tot24 <- pairwise.wilcox.test(vfas24$total, vfas24$treatment, p.adjust.method = 'none')
val24 <- pairwise.wilcox.test(vfas24$valerate, vfas24$treatment, p.adjust.method = 'none')

write.csv(but24$p.value, file = 'but24.csv')
write.csv(tot24$p.value, file = 'tot24.csv')
write.csv(val24$p.value, file = 'val24.csv')



pairwise.wilcox.test(vfas$total, vfas$treatment, p.adjust.method = 'none')


library(vegan)

rownames(vfas.matrix)
#meta <- meta %>% filter(timepoint !=4)
#rownames(meta) <- meta$group
vfas.matrix <- vfas.matrix[rownames(vfas.matrix)%in%meta$group,]
meta$treatment_timepoint <- paste(meta$treatment, meta$timepoint, sep = '_')
vfas.treat.nmds <- NMDS_ellipse(metadata = meta, OTU_table = vfas.matrix, grouping_set = 'treatment_timepoint')

library(funfuns)

####################

flow_points <- vfas.treat.nmds[[1]]
flow_ell <- vfas.treat.nmds[[2]]
flow_ell$group

flow_points %>% filter(treatment %in% c(1, 3)) %>% ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=treatment), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) +
  geom_path(data=filter(flow_ell, group %in% c(1,3)), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Day 0', subtitle = 'NMDS projection of Bray-Curtis community structure similarities') + geom_text(aes(label=group))



flow_points %>% filter(treatment %in% c(1:10)) %>% ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=treatment), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) +
  geom_path(data=filter(flow_ell, group %in% c(1:10)), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Day 0', subtitle = 'NMDS projection of Bray-Curtis community structure similarities') #+ geom_text(aes(label=group))



