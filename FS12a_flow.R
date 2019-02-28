setwd('~/FS12/FS12_flow/')


Oct12 <- read.csv('FS12a_Flow_Oct12.csv', as.is = TRUE, check.names = FALSE)


Oct13 <- read.csv('FS12a_Flow_Oct13.csv', check.names = FALSE, as.is = TRUE)


colnames(Oct12)[-1]
colnames(Oct13)[-1]
Oct12$necropsy <- 'Oct12'
Oct13$necropsy <- 'Oct13'

Oct13[,match(colnames(Oct12), colnames(Oct13))]

Oct13 <- Oct13[,match(colnames(Oct12), colnames(Oct13))]

colnames(Oct12) == colnames(Oct13)



#flow <- rbind(Oct12)
#flow <- rbind(Oct13)
flow <- rbind(Oct12, Oct13)
flow$ID

flow$pignum <- sub('([0-9]+)\\-([0-9]+)','\\1', flow$ID)

flow$treatment <- sub('([0-9]+)\\-([0-9]+)','\\2', flow$ID)

flow$ID
#flow <- flow[-(grep('LPC', flow$ID)),]


meta <- flow[,c(1,34,35, 36)]
rownames(flow) <- flow$pignum
flow <- flow[,-c(1,34,35, 36)]

library(vegan)
library(tidyverse)


sort(rowSums(flow))
sort(colSums(flow))

flow <- flow[,colSums(flow)>2000]
sort(colSums(flow))
flow <- flow[rowSums(flow)>14120,]

flow.r <- rrarefy(flow, 14120)

flow.r_CD3 <- flow.r[,grep('CD3\\+', colnames(flow.r))]

#####
meta$treatment <- factor(as.numeric(meta$treatment))
meta <- na.exclude(meta)
rownames(meta) <- meta$pignum
# cd4 POS
CD4POS <- flow.r_CD3[,grep('CD4\\+CD8a-', colnames(flow.r_CD3))]
CD4POS <- CD4POS/rowSums(CD4POS)*100
CD4POS <- CD4POS[rownames(CD4POS) %in% rownames(meta),]

meta1 <- meta[rownames(meta) %in% rownames(CD4POS),]
meta <- meta[rownames(meta) %in% rownames(CD4POS),]

rownames(meta1) == rownames(CD4POS)

CD4POS.meta <- cbind(meta1, CD4POS)



# cd8 POS
CD8POS <- flow.r_CD3[,grep('CD4-CD8a\\+', colnames(flow.r_CD3))]
CD8POS <- CD8POS/rowSums(CD8POS)*100
CD8POS <- CD8POS[rownames(CD8POS) %in% rownames(meta),]
CD8POS.meta <- cbind(meta1, CD8POS)
# dn
DN <- flow.r_CD3[,grep('CD4-CD8a-', colnames(flow.r_CD3))]
DN <- DN/rowSums(DN)*100
DN <- DN[rownames(DN) %in% rownames(meta),]
DN.meta <- cbind(meta1, DN)
# dp
DP <- flow.r_CD3[,grep('CD4\\+CD8a\\+', colnames(flow.r_CD3))]
DP <- DP/rowSums(DP)*100
DP <- DP[rownames(DP) %in% rownames(meta),]
DP.meta <- cbind(meta1, DP)



flow.r_CD3

flow.r_CD3 <- (flow.r_CD3/rowSums(flow.r_CD3)) *100

flow.r <- (flow.r/rowSums(flow.r)) *100

meta <- meta[meta$pignum %in% rownames(flow.r),]
flow.r <- flow.r[rownames(flow.r) %in% meta$pignum,]
rownames(flow.r) == meta$pignum
rownames(meta) <- meta$pignum

flow.meta <- cbind(meta, flow.r)

flow.r_CD3 <- flow.r_CD3[rownames(flow.r_CD3) %in% meta$pignum,]
flowCD3.meta <- cbind(meta, flow.r_CD3)
#meta <- meta[-grep('LPC', meta$treatment),]
flow.r <- flow.r[rownames(flow.r) %in% rownames(meta),]

rownames(flow.r) ==  rownames(meta)

#meta$treatment <- factor(as.numeric(meta$treatment))

flow_by_treat <- NMDS_ellipse(metadata = meta, OTU_table = flow.r, grouping_set = 'treatment')
meta$treatment
flow_points <- flow_by_treat[[1]]
flow_ell <- flow_by_treat[[2]]
flow_ell$group

flow_points %>% filter(treatment %in% c(1, 3)) %>% ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=treatment), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) +
  geom_path(data=filter(flow_ell, group %in% c(1,3)), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Day 0', subtitle = 'NMDS projection of Bray-Curtis community structure similarities') #+ geom_text(aes(label=pignum))



flow_points %>% filter(treatment %in% c(1:10)) %>% ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=treatment), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) +
  geom_path(data=filter(flow_ell, group %in% c(1:10)), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Day 0', subtitle = 'NMDS projection of Bray-Curtis community structure similarities') #+ geom_text(aes(label=pignum))







#############

























flow.gather <- gather(flow.meta, key = cell_type, value = percent_CD45, -(ID:treatment))
flowCD3.gather <- gather(flowCD3.meta, key = cell_type, value = percent_CD3, -(ID:treatment))


#flow.gather$treatment <- factor(as.numeric(flow.gather$treatment))
#flowCD3.gather$treatment <- factor(as.numeric(flowCD3.gather$treatment))
CD4POS.gather <- gather(CD4POS.meta, key = cell_type, value = percent_CD45_CD3_CD4, -(ID:treatment))
CD8POS.gather <- gather(CD8POS.meta, key = cell_type, value = percent_CD45_CD3_CD8, -(ID:treatment))
DN.gather <- gather(DN.meta, key = cell_type, value = percent_CD45_CD3_CD4n_CD8n, -(ID:treatment))
DP.gather <- gather(DP.meta, key = cell_type, value = percent_CD45_CD3_CD4_CD8, -(ID:treatment))



summar <- flow.gather %>% group_by(necropsy, cell_type) %>% summarise(mean=mean(percent_CD45), n=n())# %>% ggplot(aes(x=cell_type, y=mean, color=necropsy)) + geom_point()

summar <- flow.gather %>% group_by(necropsy, treatment, cell_type) %>% summarise(mean=mean(percent_CD45), n=n())# %>% ggplot(aes(x=cell_type, y=mean, color=necropsy)) + geom_point()


ggplot(summar, aes(x=necropsy, y=mean, color=treatment)) + geom_point() + facet_wrap(~cell_type, scales = 'free')

filter(flow.gather, treatment %in% c(1:10)) %>% ggplot(aes(x=treatment, y=percent_CD45, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + geom_point(aes(color=necropsy))+ facet_wrap(~cell_type, scales = 'free')

# all treatments, both necropsies

filter(flow.gather, treatment %in% c(1:10)) %>%
  ggplot(aes(x=treatment, y=percent_CD45, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~cell_type, scales = 'free') + 
  ggtitle("Immune cell phenotypes, Cecal tissue", subtitle = 'Relative to total CD45+ events')

# all treatments Oct 12
filter(flow.gather, treatment %in% c(1:10) & necropsy == 'Oct12') %>%
  ggplot(aes(x=treatment, y=percent_CD45, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~cell_type, scales = 'free') + ggtitle('October 12 Necropsy')

# all treatments Oct 13
filter(flow.gather, treatment %in% c(1:10) & necropsy == 'Oct13') %>%
  ggplot(aes(x=treatment, y=percent_CD45, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~cell_type, scales = 'free') + ggtitle('October 13 Necropsy')


#Con and RPS both nec
filter(flow.gather, treatment %in% c(1,3)) %>% ggplot(aes(x=treatment, y=percent_CD45, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~cell_type, scales = 'free')

# Con and RPS Oct 12
filter(flow.gather, treatment %in% c(1,3) & necropsy == 'Oct12') %>%
  ggplot(aes(x=treatment, y=percent_CD45, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~cell_type, scales = 'free') + ggtitle('October 12 Necropsy')

# Con and RPS Oct 13
filter(flow.gather, treatment %in% c(1,3) & necropsy == 'Oct13') %>%
  ggplot(aes(x=treatment, y=percent_CD45, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~cell_type, scales = 'free') + ggtitle('October 13 Necropsy')

######## CD3+ only  ##########

# all treatments, both necropsies
filter(flowCD3.gather, treatment %in% c(1:10)) %>% ggplot(aes(x=treatment, y=percent_CD3, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~cell_type, scales = 'free') + 
  ggtitle('Immune cell phenotypes, Cecal Tissue, 4 weeks post weaning', subtitle = 'relative to total CD3+ events')

# all treatments Oct 12
filter(flowCD3.gather, treatment %in% c(1:10) & necropsy == 'Oct12') %>%
  ggplot(aes(x=treatment, y=percent_CD3, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~cell_type, scales = 'free') + ggtitle('October 12 Necropsy')

# all treatments Oct 13
filter(flowCD3.gather, treatment %in% c(1:10) & necropsy == 'Oct13') %>%
  ggplot(aes(x=treatment, y=percent_CD3, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~cell_type, scales = 'free') + ggtitle('October 13 Necropsy')


#Con and RPS both nec
filter(flowCD3.gather, treatment %in% c(1,3)) %>% ggplot(aes(x=treatment, y=percent_CD3, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~cell_type, scales = 'free')

# Con and RPS Oct 12
filter(flowCD3.gather, treatment %in% c(1,3) & necropsy == 'Oct12') %>%
  ggplot(aes(x=treatment, y=percent_CD3, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~cell_type, scales = 'free') + ggtitle('October 12 Necropsy')

# Con and RPS Oct 13
filter(flowCD3.gather, treatment %in% c(1,3) & necropsy == 'Oct13') %>%
  ggplot(aes(x=treatment, y=percent_CD3, fill = treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~cell_type, scales = 'free') + ggtitle('October 13 Necropsy')


ggplot(flowCD3.gather, aes(x=treatment, y=percent_CD3, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +geom_point()+ facet_wrap(~cell_type, scales = 'free')

### CD4 pos ###
CD4POS.gather %>% filter(treatment %in% c(1:10)) %>%  ggplot(aes(x=treatment, y=percent_CD45_CD3_CD4, group=treatment, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~cell_type, scales = 'free') + ggtitle('Immune Cell Phenotypes, Cecal Tissue, 4 Weeks post-weaning', subtitle = 'relative to CD45+CD3+CD4+CD8-')
CD8POS.gather %>% filter(treatment %in% c(1:10)) %>%  ggplot(aes(x=treatment, y=percent_CD45_CD3_CD8, group=treatment, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~cell_type, scales = 'free')+ ggtitle('Immune Cell Phenotypes, Cecal Tissue, 4 Weeks post-weaning', subtitle = 'relative to CD45+CD3+CD4-CD8+')
DN.gather %>% filter(treatment %in% c(1:10)) %>%  ggplot(aes(x=treatment, y=percent_CD45_CD3_CD4n_CD8n, group=treatment, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~cell_type, scales = 'free')+ ggtitle('Immune Cell Phenotypes, Cecal Tissue, 4 Weeks post-weaning', subtitle = 'relative to CD45+CD3+CD4-CD8-')
DP.gather %>% filter(treatment %in% c(1:10)) %>%  ggplot(aes(x=treatment, y=percent_CD45_CD3_CD4_CD8, group=treatment, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~cell_type, scales = 'free')+ ggtitle('Immune Cell Phenotypes, Cecal Tissue, 4 Weeks post-weaning', subtitle = 'relative to CD45+CD3+CD4+CD8+')

