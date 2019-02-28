setwd('~/Documents/FS12/final/')
getwd()
library(phyloseq)
library(tidyverse)


#library(ggscinames)


library(funfuns)
list.files()


meta <- read.csv('FS12_final_meta.csv', header = TRUE, stringsAsFactors = FALSE)
shared <- read_delim('FS12.shared', delim = '\t') %>% as.data.frame()





# shared <- as.matrix(read.csv('ESV_table.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE))
# shared <- read.table('FS12.mothur.shared', header = TRUE)
# taxa <- read.csv('dada_taxa.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
taxa <- extract_mothur_tax('FS12.taxonomy')
taxa[1:10,1:8]

rownames(shared) <- shared$Group
shared <- shared[,-c(1,2,3)]

hist(rowSums(shared), breaks = 50)

rowSums(shared)

shared[1:10,1:10]

sum(rowSums(shared) < 2000)
sum(rowSums(shared) < 2500)
sum(rowSums(shared) < 3000)
sum(rowSums(shared) < 2000)


# colnames(shared) <- gsub('ESV', 'shared', colnames(shared))
# colnames(taxa)[2] <- 'shared'

# taxa$shared <- gsub('ESV', 'shared', taxa$shared)

#rownames(shared) == meta$sample_ID


shared <- shared[rowSums(shared) > 1250,] # remove samples with less than 1500 reads
shared <- shared[,colSums(shared > 0) > 3] # removes otus that are detected in fewer than 4 samples globally (including timecourse data)
shared <- shared[,colSums(shared) > 5] # at least 10 observations globally

length(colnames(shared)) # 2993 total OTUs detected


rownames(shared) %in% meta$sample_ID
meta <- meta[meta$sample_ID %in% rownames(shared),]

shared <- shared[rownames(shared) %in% meta$sample_ID,]

rownames(shared) == meta$sample_ID

rownames(meta) <- meta$sample_ID

meta <- meta[order(meta$sample_ID),]
shared <- shared[order(rownames(shared)),]

rownames(shared) == meta$sample_ID

# taxa <- taxa[taxa$shared %in% colnames(shared),]
# taxa <- taxa %>% select(-seq,-shared, everything(), seq, shared) 

# meta$pig_pen
# 
# meta$pignum <- as.numeric(gsub('P([0-9]+)', '\\1', meta$pig_pen))
# meta$pen <- as.numeric(gsub('N([0-9]+)', '\\1', meta$pig_pen))
# rooms <- read.csv('Rooms.csv')
# 
# 
# meta$replicate <- gsub('([0-9][0-9][0-9])([0-9])', '\\2',meta$pignum)
# 
# grep('[0-9][0-9][0-9][0-9]', meta$pignum)
# 
# meta <- meta[-c(31,43,94),]
# meta$pignum[30] <- 126
# 
# meta$pignum[grep(273, meta$pignum)] <- 237
#  
# 
# shared <- shared[rownames(shared) %in% rownames(meta),]
# 
# 
# meta$treatment <- ifelse(meta$pignum %in% rooms$X6, 'control',
#                              ifelse(meta$pignum %in% rooms$X7, 'RPS',
#                                     ifelse(meta$pignum %in% rooms$X8, 'Acid',
#                                            ifelse(meta$pignum %in% rooms$X9, 'ZnCu',
#                                                   ifelse(meta$pignum %in% rooms$X10, 'RCS',
#                                                          ifelse(meta$pignum %in% rooms$X11, 'Bglu', NA))))))
# 
meta$set <- paste(meta$experiment, meta$day, meta$tissue, meta$treatment, sep = '_')

nnnn <- meta %>% group_by(set) %>% summarise(N=n())


rownames(meta) == rownames(shared)
meta <- meta[rownames(meta) %in% rownames(shared),]
rownames(shared) %in% rownames(meta)
shared <- shared[match(rownames(meta), rownames(shared)),]

rownames(meta) == rownames(shared)


# rownames(taxa) <- taxa$shared

taxa <- taxa[taxa$OTU %in% colnames(shared),]

#match(colnames(shared), rownames(taxa))
# taxa <- taxa[match(colnames(shared), rownames(taxa)),]

# all(colnames(shared) == rownames(taxa))

# is.na(meta$treatment)
########## Tree?  ###########


library(phyloseq)
p_meta <- sample_data(meta) 
p_taxa <- import_mothur(mothur_constaxonomy_file = 'FS12.taxonomy')

# p_taxa <- tax_table(taxa)
colnames(p_taxa) <- colnames(taxa)[-c(1,2,3)]
# taxa_names(p_taxa) <- rownames(taxa)
rank_names(p_taxa)

meta$treatment
FS12 <- phyloseq(p_meta, p_taxa, otu_table(shared, taxa_are_rows = FALSE))

#NMDS_ellipse(metadata = meta, OTU_table = shared, grouping_set = 'set')

FS12@sam_data$experiment

# FS12a <- subset_samples(FS12, experiment == 'X12a')
# 
# 
# 
# FS12a.otu <- FS12a@otu_table
# 
# rowSums(FS12a.otu)
# FS12a@sam_data$sample_ID == sample_names(FS12a)
# 
# 
# rownames(FS12a@otu_table) ==FS12a@sam_data$sample_ID
# 
# # FS12a_pens <- read.csv('FS12a_pentreat.csv') %>% na.exclude()
# 
# # FS12a@sam_data$pig_pen
# # 
# # 
# # 
# # 
# # FS12a@sam_data$pen
# # 
# # pentreatvec <- FS12a_pens$Diet
# # names(pentreatvec) <- FS12a_pens$Pen
# # 
# # 
# # # 1cont, 2abx, 3rps, 6acid, 8zn, 9rcs, 10bglu
# # 
# # treatvec <- c('control', 'Abx', 'RPS', '4', '5', 'Acid', '7', 'ZnCu', 'RCS', 'Bglu')
# # names(treatvec) <- c(1,2,3,4,5,6,7,8,9,10)
# # 
# # 
# # treats <- pentreatvec[FS12a@sam_data$pen]
# # 
# # FS12a@sam_data$treatment <- treatvec[treats]
# # 
# # sort(rowSums(FS12a@otu_table))
# # 
# # 
# # FS12a@sam_data$day <- sub('W','',FS12a@sam_data$day)
# # FS12a@sam_data$day <- sub('[Dd]+','',FS12a@sam_data$day)
# # 
# # 
# # FS12a@sam_data$day <- as.numeric(FS12a@sam_data$day)
# 
# library(vegan)
# 
# # ERROR day not found
# #prune_samples(day %in% c(0,23), FS12a)
# 
# sample_data(FS12a)
# 
# # FS12a <- subset_samples(FS12a, day %in% c(0,23))
# 
# FS12a <- FS12a %>% subset_samples(day %in% c('D0', 'D23', 'D30'))
# FS12a <- FS12a %>% subset_samples(!(day =='D30' & treatment =='Phyto'))
# 
# FS12a.otu.rare <- rrarefy(FS12a@otu_table, sample = min(rowSums(FS12a@otu_table)))
# 
# FS12a@sam_data$shan <- diversity(FS12a.otu.rare)
# 
# FS12a@sam_data %>% ggplot(aes(x=treatment, y=shan)) + geom_boxplot()+ facet_wrap(~day)
# 
# FS12a@sam_data$set <- paste(FS12a@sam_data$treatment, FS12a@sam_data$day)
# # Only have data from D0 and D23 
# 
# #FS12a.o
# 
# FS12a.ord <- NMDS_ellipse(metadata = FS12a@sam_data, OTU_table = FS12a.otu.rare, grouping_set = 'set')
# 
# FS12a.meta <- FS12a.ord[[1]]
# FS12a.ell <- FS12a.ord[[2]]
# FS12a.ell$treatment <- gsub('(.*) (.*)','\\1',FS12a.ell$group)
# FS12a.ell$day <- gsub('(.*) (.*)','\\2',FS12a.ell$group)
# 
# 
# 
# FS12a.meta %>% filter(day == 'D0') %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   geom_point(alpha=0.5, size=2) +  
#   geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
#   geom_path(data = filter(FS12a.ell, day == 'D0'), aes(x=NMDS1,y=NMDS2), size=1.05)
#   # scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
#   # annotate(geom='text', label='Day 0', x=-1.25, y=.6, size = 7)
# 
# FS12a.meta %>% filter(day == 'D23' & treatment %in% c('Control', 'Malto')) %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   geom_point(alpha=0.5, size=2) +  
#   geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
#   geom_path(data = filter(FS12a.ell, day == 'D23' & treatment %in% c('Control', 'Malto')), aes(x=NMDS1,y=NMDS2), size=1.05) + 
#   scale_color_manual(values = c('blue1', 'green4'), labels = c('Maltodextrin', 'Control')) + ggtitle('Community similarity at day 23 post-weaning (feces)',
#                                                               subtitle = 'uncorrected PERMANOVA pvalue = 0.04, pseudo F = 1.9231317') + 
# 
# FS12a.foadon <- subset_samples(FS12a, day == 23 & !is.na(treatment))
# 
# FS12a.adonres <- pairwise.adonis(FS12a.foadon@otu_table, factors = FS12a.foadon@sam_data$treatment, permutations = 999)
# FS12a.adonres[grep(4, FS12a.adonres$pairs),]
# 
# # Deseq
# 
# 
# Deseq.quickplot()
# FS12a.23@sam_data$set
# library(DESeq2)
# FS12a.23 <- subset_samples(FS12a, day > 0)
# FS12a.23.glom <- tax_glom(FS12a.23,taxrank = 'Genus')
# 
# rank_names(FS12a.23)
# FS12.de <- phyloseq_to_deseq2(FS12a.23.glom, ~set)
# 
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# 
# FS12a.malto <- Deseq.quickplot(DeSeq.object = FS12.de,
#                              phyloseq.object = FS12a.23.glom, pvalue = .05,
#                              contrast.vector = c('set', 'control 23', '4 23'),
#                              taxlabel = 'Genus', colors = c('green4', 'blue1')) 
# 
# FS12a.malto[[1]] + xlab('Genus') + ggtitle('Differentially abundant genera 23 days post-weaning (feces)',
#                                            subtitle = 'As determined by DeSeq2, p < 0.05') 

############ FS12b from here down ! ##############


FS12b <- subset_samples(FS12, experiment == 'X12b' & pignum != 101)

FS12b <- subset_samples(FS12b, treatment %in% c('Control', 'RPS', 'Acid', 'RCS'))

# making sure factors are set correctly

FS12b@sam_data$treatment <- factor(FS12b@sam_data$treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))

FS12b@sam_data$set

library(vegan)
library(funfuns)
library(tidyverse)
# FS12b@otu_table

# shared[shared[,'shared2424'] > 0,'shared2424']

FS12b <- prune_taxa(taxa_sums(FS12b) > 2, FS12b) # removes singletons
FS12b@sam_data
# prune_samples(!(is.na(treatment)), FS12b)
# subset_samples(!(is.na(treatment), FS12b))
# which(is.na(FS12b@sam_data$treatment))

FS12b <- prune_samples(!is.na(FS12b@sam_data$treatment), FS12b)

min(sample_sums(FS12b))

taxa_sums(FS12b) > 2

grep('Clostridium_sensu_stricto_1', FS12b@tax_table[,6])

FS12b_metanmds <- NMDS_ellipse(metadata = FS12b@sam_data,
                               OTU_table = rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table))),
                               grouping_set = 'set',distance_method = 'jaccard')

FS12b_metanmds

nums <- FS12b_metanmds[[1]] %>% group_by(set) %>% summarise(N=n())

###


FS12b_jac <- vegdist(FS12b@otu_table, method = 'jaccard')
FS12b_jac

attr(FS12b_jac, which = 'Labels') == FS12b@sam_data$sample_ID
dispers <- betadisper(FS12b_jac, group = FS12b@sam_data$set)
pdispers <- permutest(dispers, pairwise = TRUE)
pdispers$pairwise$observed
dispersdf <- data.frame(dispers$distances)
dispersdf$group <- rownames(dispersdf)
FS12b@sam_data$sample_ID == dispersdf$group


meta$sample_ID %in% dispersdf$group

colnames(dispersdf)[2] <- 'sample_ID'
FS12b_meta <- FS12b_metanmds[[1]]

FS12b_meta <- merge(FS12b_meta, dispersdf, by='sample_ID')
FS12b_meta$day <- as.numeric(gsub('D', '', FS12b_meta$day))
FS12b_meta$dayfact <- factor(FS12b_meta$day)


FS12b_meta$treatment <- factor(FS12b_meta$treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))

FS12b_meta$shan <- diversity(rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table))))

FS12b_meta$rich <- specnumber(rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table))))
FS12b_meta$even <- FS12b_meta$shan/log(FS12b_meta$rich)

# # fecal shannon
# FS12b_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=dayfact, y=shan, group=set, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Fecal Shannon Diversity (alpha) over time') #+ geom_text(aes(label=pignum))
# 
# #fecal dispersion
# FS12b_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=dayfact, y=dispers.distances, group=set, fill = treatment)) + geom_boxplot() + 
#   scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + ggtitle('Fecal community dispersion over time')#+ geom_text(aes(label=pignum))


#fecal shannon
FS12b_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~day)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Fecal Shannon Diversity (alpha) over time')  + geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)#+ geom_text(aes(label=pignum))


#fecal even
FS12b_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=even, group=set, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~day)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Fecal evenness over time')+ geom_jitter(shape = 21, stroke=1.2, size=2, width = .2) #+ geom_text(aes(label=pignum))

#fecal rich
FS12b_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=rich, group=set, fill = treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) + facet_wrap(~day)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Fecal richness (num OTUs) over time')+ geom_jitter(shape = 21, stroke=1.2, size=2, width = .2) #+ geom_text(aes(label=pignum))


#fecal dispersion
FS12b_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=dispers.distances, group=set, fill = treatment)) + geom_boxplot(position=position_dodge2(preserve = 'total')) + facet_wrap(~day)+
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + ggtitle('Fecal community dispersion over time')+ geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)#+ geom_text(aes(label=pignum))

FS12b_meta %>% filter(tissue == 'F') %>%
  ggplot(aes(x=day, y=dispers.distances, group=treatment, fill = treatment, color=treatment)) +
  #scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  +
  ggtitle('Fecal community dispersion over time') + geom_smooth()



get_pairs <- function(df){
  pp <- pairwise.wilcox.test(df$dispers.distances, df$treatment, p.adjust.method = 'none')
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <- paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps
  
}

disper_fecal_tests <- FS12b_meta %>% filter(tissue =='F') %>% group_by(day) %>% 
  nest() %>% mutate(pps = map(data, get_pairs)) %>%
  select(day, pps) %>% unnest() %>% select(day, starts_with('control'))


get_pairs <- function(df){
  pp <- pairwise.wilcox.test(df$shan, df$treatment, p.adjust.method = 'none')
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <- paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps

}

shan_fecal_tests <- FS12b_meta %>% filter(tissue =='F') %>% group_by(day) %>% 
  nest() %>% mutate(pps = map(data, get_pairs)) %>%
  select(day, pps) %>% unnest() %>% select(day, starts_with('control'))



# dispersion tissues
FS12b_meta %>% filter(tissue == 'X') %>% ggplot(aes(x=day, y=dispers.distances, group=set, fill = treatment)) + geom_boxplot() + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))

FS12b_meta %>% filter(tissue == 'C') %>% ggplot(aes(x=day, y=dispers.distances, group=set, fill = treatment)) + geom_boxplot() + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))

FS12b_meta %>% filter(tissue == 'I') %>% ggplot(aes(x=treatment, y=dispers.distances, group=set, fill = treatment)) + geom_boxplot() + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))

# shannon tissues
FS12b_meta %>% filter(tissue == 'X') %>% ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) + geom_boxplot() + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))

FS12b_meta %>% filter(tissue == 'C') %>% ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) + geom_boxplot() + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))

FS12b_meta %>% filter(tissue == 'I') %>% ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) + geom_boxplot() + 
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) #+ geom_text(aes(label=pignum))




disps <- FS12b_meta %>% group_by(day, treatment, tissue) %>% summarise(mean_disp = mean(dispers.distances),
                                                              sd_disp = sd(dispers.distances), 
                                                              se_disp = sd_disp/sqrt(n()))


FS12b_meta[which(is.na(FS12b_meta$treatment)),]

disps %>% filter(tissue == 'F') %>% 
  ggplot(aes(x=day, y=mean_disp, color=treatment, group=treatment)) +
  geom_line() + geom_errorbar(aes(ymin=mean_disp-se_disp, ymax=mean_disp+se_disp), width=.2)



#vegdist()
#FS12b_meta <- FS12b_metanmds[[1]] 
df_ell <- FS12b_metanmds[[2]]

df_ell$experiment <- gsub('(.*)_(.*)_(.*)_(.*)','\\1', df_ell$group)
df_ell$day <- gsub('(.*)_(.*)_(.*)_(.*)','\\2', df_ell$group)
df_ell$day <- gsub('D', '', df_ell$day)
df_ell$day <- factor(df_ell$day, levels = c(0, 2, 7, 14, 21))

df_ell$tissue <- gsub('(.*)_(.*)_(.*)_(.*)','\\3', df_ell$group)
df_ell$treatment <- gsub('(.*)_(.*)_(.*)_(.*)','\\4', df_ell$group)



FS12b_meta$day <- factor(FS12b_meta$day, levels = c(0, 2, 7, 14, 21))
FS12b_meta$treatment <- factor(FS12b_meta$treatment, levels = c('Control', 'RPS', 'Acid', 'RCS'))
df_ell$treatment <- factor(df_ell$treatment, levels = c('Control', 'RPS', 'Acid','RCS'))


greys <- FS12b_meta


    


FS12b_meta %>% filter(tissue == 'F' & day == 0) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +  
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
  geom_path(data = filter(df_ell, tissue == 'F' & day == 0), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
  annotate(geom='text', label='Day 0', x=-1.25, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'F' & day == 2) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +  
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
  geom_path(data = filter(df_ell, tissue == 'F' & day == 2), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
  annotate(geom='text', label='Day 2', x=-1.25, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'F' & day == 7) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +  
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
  geom_path(data = filter(df_ell, tissue == 'F' & day == 7), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
  annotate(geom='text', label='Day 7', x=-1.25, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'F' & day == 14) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +  
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
  geom_path(data = filter(df_ell, tissue == 'F' & day == 14), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
  annotate(geom='text', label='Day 14', x=-1.25, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'F' & day == 21) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +  
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
  geom_path(data = filter(df_ell, tissue == 'F' & day == 21), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
  annotate(geom='text', label='Day 21 feces', x=-.5, y=.6, size = 7)



FS12b_meta %>% filter(tissue == 'C' & day == 21) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  # geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +  
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
  geom_path(data = filter(df_ell, tissue == 'C' & day == 21), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
  annotate(geom='text', label='Day 21\nCecal Contents', x=-0, y=.0, size = 7)



FS12b_meta %>% filter(tissue == 'X' & day == 21) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  # geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +  
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
  geom_path(data = filter(df_ell, tissue == 'X' & day == 21), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
  annotate(geom='text', label='Day 21\nCecal Mucosa', x=-0, y=.6, size = 7)

FS12b_meta %>% filter(tissue == 'I' & day == 21) %>% 
  ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
  # geom_point(data=greys, inherit.aes = FALSE, color='grey', aes(x=MDS1, y=MDS2), size=2) + geom_point(alpha=0.5, size=2) +  
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + 
  geom_path(data = filter(df_ell, tissue == 'I' & day == 21), aes(x=NMDS1,y=NMDS2), size=1.05)+
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) + 
  annotate(geom='text', label='Day 21\nIleal Mucosa', x=-0, y=.6, size = 5)




linestuff <- FS12b_meta %>% group_by(set) %>%
  summarise(centx=mean(MDS1), centy=mean(MDS2),
            lineg=unique(paste(treatment, tissue)),
            day=unique(sub('D', '', day)))

linestuff$day <- as.numeric(linestuff$day)
linestuff <- linestuff[order(linestuff$day),]

#linestuff$lineg <- factor(lineg, levels = c())
linestuff$tissue <- gsub('[A-Za-z]+ ([FCX])', '\\1', linestuff$lineg)
linestuff$treatment <- gsub('([A-Za-z]+) ([FCX])', '\\1', linestuff$lineg)
linestuff$treatment <- factor(linestuff$treatment, levels = c('Control', 'RPS', 'Acid','ZnCu', 'RCS', 'Bglu'))


linestuff %>%filter(tissue == 'F') %>% 
  ggplot(aes(x=centx, y=centy, group=lineg, color=treatment)) +
  geom_path() + geom_point(size=5)+ geom_text(aes(label=day), color='black') +
  scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank()) +ylab('MDS2') + xlab('MDS1')



# D0 vs D2 : Early Salmonella infection changes
# D0 vs D7 : 


# scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + theme(panel.background = element_blank())
########## PW ADON HERE ########
#should split fecal and tissue 1st to reduce # of permutations...


PW.ad <- pairwise.adonis(x=rrarefy(FS12b@otu_table, min(rowSums(FS12b@otu_table))), factors = FS12b@sam_data$set, sim.method = 'bray', p.adjust.m = 'none', permutations = 9999)

PW.ad$pairs


goods <- PW.ad[grep('(.*)_(.*)_(.*)_(.*) vs (.*)_\\2_\\3_(.*)', PW.ad$pairs),]

times <- PW.ad[grep('(.*)_(.*)_(.*)_(.*) vs (.*)_(.*)_\\3_\\4', PW.ad$pairs),]

length(goods[,1])

# goods$p.adjusted <- p.adjust(p=goods$p.value,method = 'holm')

D0 <- goods[grep('D0', goods$pairs),]
D0$day <- 0

D2 <- goods[grep('D2_', goods$pairs),]
D2$day <- 2
D7 <- goods[grep('D7', goods$pairs),]
D7$day <- 7
D14 <- goods[grep('D14', goods$pairs),]
D14$day <- 14
D21 <- goods[grep('D21', goods$pairs),]
D21$day <- 21

fin <- rbind(D0, D2, D7, D14, D21)

fin$pairs <- gsub('X12b_', '', fin$pairs)
fin$pairs <- gsub('_F_', ' feces ', fin$pairs)
fin$pairs <- gsub('_C_', ' cec_cont ', fin$pairs)
fin$pairs <- gsub('_X_', ' cec_muc ', fin$pairs)
fin$pairs <- gsub('_I_', ' il_muc ', fin$pairs)
fin$pairs <- gsub('_Q_', ' tet ', fin$pairs)

write.csv(fin, 'mothur_PERMANOVA_results.csv')


to_conts <- fin[grep('Control', fin$pairs),]

to_conts$tissue <- gsub('D[0-9]+ ([A-Za-z_]+) [A-Za-z]+ vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\1', to_conts$pairs)

to_conts$treatment <- gsub('D[0-9]+ ([A-Za-z_]+) ([A-Za-z]+) vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\2', to_conts$pairs)

to_conts$treatment[to_conts$treatment == 'Control'] <- gsub('D[0-9]+ ([A-Za-z_]+) ([A-Za-z]+) vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\3', to_conts[to_conts$treatment == 'Control',]$pairs)



to_conts$p.hoch <- p.adjust(to_conts$p.value, method = 'hoch')
to_conts$p.holm <- p.adjust(to_conts$p.value, method = 'holm')
to_conts$p.fdr <- p.adjust(to_conts$p.value, method = 'fdr')
to_conts$p.fdr <- round(to_conts$p.fdr, digits = 3)
to_conts$p.fdr.lab <- ifelse(to_conts$p.fdr < 0.05, to_conts$p.fdr, NA)

to_conts$treatment <- factor(to_conts$treatment, levels=c('RPS', 'Acid', 'RCS'))

to_conts %>% filter(tissue == 'feces') %>% ggplot(aes(x=day, y=F.Model, group=treatment, fill=treatment, color=treatment, label=p.fdr.lab)) +
  geom_line(size=1.52) + geom_point(shape=21) + scale_color_manual(values=c('#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  geom_label(color='black') +
  scale_fill_manual(values=c('#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  ggtitle('Community differences compared to control group over time', subtitle = )




########### HIGH LOW TO CONT

##SHOULD ORDINATE THIS TOO##

FS12b_HL <- FS12b %>% subset_samples(treatment %in% c('Control', 'RPS') & tissue =='F')
FS12b_HL <- FS12b %>% subset_samples(treatment %in% c('Control', 'RPS'))

FS12b_HL %>% subset_samples(treatment == 'RPS') %>% sample_data() %>% select(pignum)


FS12b_HL@sam_data$shed <- ifelse(FS12b_HL@sam_data$pignum %in% c(373,321,181,392,97), 'low', 
                                 ifelse(FS12b_HL@sam_data$pignum %in% c(50, 93,355, 244), 'high', 'Control'))

FS12b_HL@sam_data$set <- paste(FS12b_HL@sam_data$day, FS12b_HL@sam_data$tissue, FS12b_HL@sam_data$shed, sep = '_')

FS12b_HL@sam_data


PW.ad <- pairwise.adonis(x=rrarefy(FS12b_HL@otu_table, min(rowSums(FS12b_HL@otu_table))), factors = FS12b_HL@sam_data$set, sim.method = 'bray', p.adjust.m = 'none', permutations = 9999)

PW.ad$pairs


goods <- PW.ad[grep('(.*)_(.*) vs \\1_.*', PW.ad$pairs),]

# times <- PW.ad[grep('(.*)_(.*)_(.*)_(.*) vs (.*)_(.*)_\\3_\\4', PW.ad$pairs),]

length(goods[,1])

# goods$p.adjusted <- p.adjust(p=goods$p.value,method = 'holm')

goods$pairs


D0 <- goods[grep('D0', goods$pairs),]
D0$day <- 0

D2 <- goods[grep('D2_', goods$pairs),]
D2$day <- 2
D7 <- goods[grep('D7', goods$pairs),]
D7$day <- 7
D14 <- goods[grep('D14', goods$pairs),]
D14$day <- 14
D21 <- goods[grep('D21', goods$pairs),]
D21$day <- 21

# I think fin and goods are the same thing right now.... why did I do this again?
fin <- rbind(D0, D2, D7, D14, D21)

fin$pairs <- gsub('X12b_', '', fin$pairs)
fin$pairs <- gsub('_F_', ' feces ', fin$pairs)
fin$pairs <- gsub('_C_', ' cec_cont ', fin$pairs)
fin$pairs <- gsub('_X_', ' cec_muc ', fin$pairs)
fin$pairs <- gsub('_I_', ' il_muc ', fin$pairs)
fin$pairs <- gsub('_Q_', ' tet ', fin$pairs)


# write.csv(fin, 'mothur_PERMANOVA_results.csv')

#within tissues
fin <- fin[grep('.* (.*) .* vs .* \\1 .*', fin$pairs),]


to_conts <- fin[grep('Control', fin$pairs),]

not_conts <- fin[-grep('Control', fin$pairs),]

to_conts$tissue <- gsub('D[0-9]+ (.*) ([A-Za-z_]+) vs D[0-9]+ .* ([A-Za-z]+)', '\\1', to_conts$pairs)
to_conts$treatment <- gsub('D[0-9]+ .* ([A-Za-z_]+) vs D[0-9]+ .* ([A-Za-z]+)', '\\2', to_conts$pairs)

# to_conts$treatment[to_conts$treatment == 'control'] <- gsub('D[0-9]+ ([A-Za-z_]+) ([A-Za-z]+) vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\3', to_conts[to_conts$treatment == 'control',]$pairs)



to_conts$p.hoch <- p.adjust(to_conts$p.value, method = 'hoch')
to_conts$p.holm <- p.adjust(to_conts$p.value, method = 'holm')
to_conts$p.fdr <- p.adjust(to_conts$p.value, method = 'fdr')
to_conts$p.fdr <- round(to_conts$p.fdr, digits = 3)
to_conts$p.fdr.lab <- ifelse(to_conts$p.fdr < 0.05, to_conts$p.fdr, NA)

# to_conts$treatment <- factor(to_conts$treatment, levels=c('RPS', 'Acid', 'ZnCu', 'RCS', 'Bglu'))

to_conts %>% filter(tissue =='feces') %>%  ggplot(aes(x=day, y=F.Model, group=treatment, fill=treatment, color=treatment, label=p.fdr.lab)) +
  geom_line(size=1.52) + geom_point(shape=21) + scale_color_brewer(palette = 'Set1') + 
  geom_label(color='black') +
  scale_fill_brewer(palette = 'Set1') + 
  ggtitle('Community differences compared to control group over time', subtitle = 'RPS only') + labs(fill='Shedding', 
                                                                                           color='Shedding')
to_conts %>% filter(tissue !='feces') %>% 
  ggplot(aes(x=tissue, y=F.Model, fill=treatment)) +
  geom_col(position = 'dodge', color='black') + geom_text(aes(label=p.fdr.lab), position=position_dodge(width = 1), vjust=1.5) + 
  ggtitle('PERMANOVA F.stat. : Difference compared to controls across tissues',
          subtitle = 'Higher values represent a greater difference compared to control')  + scale_fill_brewer(palette = 'Set1')
  


###

#HIGH LOW ORDINATE#

HIGH_LOW_NMDS <- NMDS_ellipse(OTU_table=rrarefy(FS12b_HL@otu_table, min(rowSums(FS12b_HL@otu_table))), metadata = FS12b_HL@sam_data, grouping_set = 'set')

HIGH_LOW_NMDS[[1]]$shed <- factor(HIGH_LOW_NMDS[[1]]$shed, levels = c('high', 'low', 'Control'))

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D0') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 0, RPS high/low & control')


HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D2') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5) + scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 2, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D7') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 7, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D14') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 14, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'F' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'X' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Cecal mucosa')

HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'C' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Cecal contents')

p <- HIGH_LOW_NMDS[[1]] %>% filter(tissue == 'I' & day =='D21') %>%
  ggplot(aes(x=MDS1, y=MDS2, group=set, color=shed)) + geom_point(size=3)+
  geom_segment(aes(xend =centroidX , yend = centroidY), alpha=0.5)+ scale_color_brewer(palette = 'Set1') + 
  ggtitle('Feces Day 21, RPS high/low & control, Ileal mucosa')



ggplot2::ggplot_build(p)
########### groups compared to their D0  #######

T0s <- times[grep('D0', times$pairs),]

T0s$pairs <- gsub('X12b_', '', T0s$pairs)
T0s$pairs <- gsub('_F_', ' feces ', T0s$pairs)

#4DAF4A, #377EB8)

#377EB8

#E41A1C



T0s$tissue <- gsub('D[0-9]+ ([A-Za-z_]+) [A-Za-z]+ vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\1', T0s$pairs)
T0s$treatment <- gsub('D[0-9]+ ([A-Za-z_]+) [A-Za-z]+ vs D[0-9]+ [A-Za-z_]+ ([A-Za-z]+)', '\\2', T0s$pairs)

T0s$day <- gsub('D[0-9]+ ([A-Za-z_]+) [A-Za-z]+ vs (D[0-9]+) [A-Za-z_]+ ([A-Za-z]+)', '\\2', T0s$pairs)

# what's this for??
T0s[T0s$day == "D0",]$day <- gsub('(D[0-9]+) ([A-Za-z_]+) [A-Za-z]+ vs (D[0-9]+) [A-Za-z_]+ ([A-Za-z]+)', '\\1', T0s[T0s$day == "D0",]$pairs)

T0s$day <- factor(gsub('D','',T0s$day), levels = c(2,7,14,21))
T0s$pairs

T0s$p.fdr <- round(p.adjust(T0s$p.value, 'fdr'),3)
T0s$p.fdr.lab <- ifelse(T0s$p.fdr <0.05, T0s$p.fdr, NA)
T0s$treatment <- factor(T0s$treatment, levels = c('Control', 'RPS', 'Acid', 'ZnCu', 'RCS', 'Bglu'))



T0s %>% filter(tissue == 'feces') %>% ggplot(aes(x=day, y=F.Model, group=treatment, fill=treatment, color=treatment, label=p.fdr.lab)) +
  geom_line(size=1.52) + geom_point(shape=21) + scale_color_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  geom_label(color='black') +
  scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))+ 
  ggtitle("Community differences compared to each group's Day 0 conformation", 
          subtitle = 'FDR corrected pvalues shown in boxes') + xlab('Day (vs Day 0)')



# T0s[T0s$day == "D0",]

# gsub('(D[0-9]+) ([A-Za-z_]+) [A-Za-z]+ vs (D[0-9]+) [A-Za-z_]+ ([A-Za-z]+)', '\\1', T0s[T0s$day == "D0",]$pairs)

# T0s[T0s$treatment == 'control',]

# FS12b_meta %>% filter(treatment %in% c('control', 'RCS')) %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   facet_wrap(~day) + geom_text(aes(label=pignum)) +
#   geom_segment(aes(xend =centroidX , yend = centroidY)) + 
#   geom_path(data = filter(df_ell, treatment %in% c('RCS', 'control') & tissue == 'F'), aes(x=NMDS1,y=NMDS2))


# FS12b_meta %>% filter(treatment %in% c('control', 'Acid')) %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   facet_wrap(~day) + geom_text(aes(label=pignum)) +
#   geom_segment(aes(xend =centroidX , yend = centroidY)) + 
#   geom_path(data = filter(df_ell, treatment %in% c('Acid', 'control') & tissue == 'F'), aes(x=NMDS1,y=NMDS2))


# FS12b_meta %>% filter(treatment %in% c('control', 'ZnCu')) %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   facet_wrap(~day) + geom_text(aes(label=pignum)) +
#   geom_segment(aes(xend =centroidX , yend = centroidY)) + 
#   geom_path(data = filter(df_ell, treatment %in% c('ZnCu', 'control') & tissue == 'F'), aes(x=NMDS1,y=NMDS2))
# 
# FS12b_meta %>% filter(treatment %in% c('control', 'Bglu') & tissue == 'F')  %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   facet_wrap(~day) + geom_text(aes(label=pignum)) +
#   geom_segment(aes(xend =centroidX , yend = centroidY)) + 
#   geom_path(data = filter(df_ell, treatment %in% c('Bglu', 'control') & tissue == 'F'), aes(x=NMDS1,y=NMDS2))
# 
# 
# 
# FS12b_meta %>% filter(tissue == 'F' & treatment %in% c('control', 'RPS', 'Bglu')) %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   facet_wrap(~day) + geom_text(aes(label=pignum)) +
#   geom_segment(aes(xend =centroidX , yend = centroidY)) + 
#   geom_path(data = filter(df_ell, tissue == 'F' & treatment %in% c('control', 'RPS', 'Bglu')), aes(x=NMDS1,y=NMDS2)) + ggtitle('Similarity of fecal communities over time')
# 
# FS12b_meta %>% filter(tissue == 'C' & treatment %in% c('control', 'RPS', 'Bglu')) %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   facet_wrap(~day, scales = 'free') + geom_text(aes(label=pignum)) +
#   geom_segment(aes(xend =centroidX , yend = centroidY)) + 
#   geom_path(data = filter(df_ell, tissue == 'C' & treatment %in% c('control', 'RPS', 'Bglu')), aes(x=NMDS1,y=NMDS2))
# 
# FS12b_meta %>% filter(tissue == 'X' & treatment %in% c('control', 'RPS', 'Bglu')) %>% 
#   ggplot(aes(x=MDS1, y=MDS2, fill=treatment, color=treatment)) +
#   facet_wrap(~day, scales = 'free') + geom_text(aes(label=pignum)) +
#   geom_segment(aes(xend =centroidX , yend = centroidY)) + 
#   geom_path(data = filter(df_ell, tissue == 'X' & treatment %in% c('control', 'RPS', 'Bglu')), aes(x=NMDS1,y=NMDS2))
# 



###############################################

library(DESeq2)

# up first is the high vs control and low vs control stuff

# FS12_RPS <- subset_samples(FS12b, treatment %in% c('Control', 'RPS'))
# FS12_RPS@sam_data$day <- factor(FS12_RPS@sam_data$day, levels = c('D0', 'D2','D7', 'D14', 'D21'))
# 
# FS12_RPS@sam_data$shed <- ifelse(FS12_RPS@sam_data$pignum %in% c(373,321,181,392,97), 'low', 
#                                  ifelse(FS12_RPS@sam_data$pignum %in% c(50, 244, 355, 373), 'high', 'Control'))
# FS12_RPS@sam_data$shed <- factor(FS12_RPS@sam_data$shed, levels = c('Control', 'high', 'low'))
# 
# # dont think i need this
# FS12_RPS@sam_data$set <- paste(FS12_RPS@sam_data$set, FS12_RPS@sam_data$shed, sep = '_')
# 
# # Keeping samples separate by day #
# 
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D0' & tissue == 'F')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# 
# 
# D0_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
#                               phyloseq.object = FS12b.glom, pvalue = .05,
#                               contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# 
# D0_conthigh[[1]] + ggtitle('Control vs RPS_high Day 0') 
# 
# D0_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                               phyloseq.object = FS12b.glom, pvalue = .05,
#                               contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# 
# D0_contlow[[1]]+ ggtitle('Control vs RPS_low Day 0') 
# 
# ## DAY 2
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D2' & tissue == 'F')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# 
# 
# D2_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                phyloseq.object = FS12b.glom, pvalue = .05,
#                                contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# 
# D2_conthigh[[1]] + ggtitle('Control vs RPS_high D2')
# 
# D2_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                               phyloseq.object = FS12b.glom, pvalue = .05,
#                               contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# 
# D2_contlow[[1]] + ggtitle('Control vs RPS_low D2')
# 
# ##
# # FS12b.glom <- FS12_RPS
# # FS12b.glom <- subset_samples(FS12b.glom, day =='D2' & tissue == 'F')
# # 
# # FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# # 
# # FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# # library(DESeq2)
# # FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# # resultsNames(FS12.de)
# # 
# # 
# # 
# # D2_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
# #                                phyloseq.object = FS12b.glom, pvalue = .05,
# #                                contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# # 
# # D2_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
# #                               phyloseq.object = FS12b.glom, pvalue = .05,
# #                               contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# # 
# ##
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D7' & tissue == 'F')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# D7_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                phyloseq.object = FS12b.glom, pvalue = .05,
#                                contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# 
# D7_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                               phyloseq.object = FS12b.glom, pvalue = .05,
#                               contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# 
# D7_contlow[[1]] + ggtitle('Control vs RPS_low Day 7')
# 
# ##
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D14' & tissue == 'F')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# 
# 
# D14_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                phyloseq.object = FS12b.glom, pvalue = .05,
#                                contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# 
# D14_conthigh[[1]] + ggtitle('Control vs RPS_high Day 14')
# 
# D14_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                               phyloseq.object = FS12b.glom, pvalue = .05,
#                               contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# D14_contlow[[1]] + ggtitle('Control vs RPS_low Day 14')
# 
# ##
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'F')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# 
# 
# D21_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                phyloseq.object = FS12b.glom, pvalue = .05,
#                                contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# D21_conthigh[[1]] + ggtitle('Control vs RPS_high Day 21 feces')
# 
# 
# D21_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                               phyloseq.object = FS12b.glom, pvalue = .05,
#                               contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# 
# D21_contlow[[1]] + ggtitle('Control vs RPS_low Day 21 feces')
# 
# #### tissue C ###
# 
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'C')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# 
# 
# D21_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                 phyloseq.object = FS12b.glom, pvalue = .05,
#                                 contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# D21_conthigh[[1]] + ggtitle('Control vs RPS_high Day 21 cec_cont')
# 
# 
# D21_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                phyloseq.object = FS12b.glom, pvalue = .05,
#                                contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# 
# D21_contlow[[1]] + ggtitle('Control vs RPS_low Day 21 cec_cont')
# 
# ## tissue X
# 
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'X')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# 
# 
# D21_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                 phyloseq.object = FS12b.glom, pvalue = .05,
#                                 contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# D21_conthigh[[1]] + ggtitle('Control vs RPS_high Day 21 cec_muc')
# 
# 
# D21_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                phyloseq.object = FS12b.glom, pvalue = .05,
#                                contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# 
# D21_contlow[[1]] + ggtitle('Control vs RPS_low Day 21 cec_muc')
# 
# ## tissue I
# 
# FS12b.glom <- FS12_RPS
# FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'I')
# 
# FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
# 
# FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
# library(DESeq2)
# FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
# resultsNames(FS12.de)
# 
# 
# 
# D21_conthigh <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                 phyloseq.object = FS12b.glom, pvalue = .05,
#                                 contrast.vector = c('shed', 'high', 'Control') ,taxlabel = 'Genus')
# D21_conthigh[[1]] + ggtitle('Control vs RPS_high Day 21 cec_muc')
# 
# 
# D21_contlow <- Deseq.quickplot(DeSeq.object = FS12.de,
#                                phyloseq.object = FS12b.glom, pvalue = .05,
#                                contrast.vector = c('shed', 'low', 'Control') ,taxlabel = 'Genus')
# 
# D21_contlow[[1]] + ggtitle('Control vs RPS_low Day 21 il_muc')
# 





##

################################ RPS SPLIT #########################
FS12b@sam_data$pignum
FS12_RPS <- subset_samples(FS12b, treatment == 'RPS')
FS12_RPS@sam_data$day <- factor(FS12_RPS@sam_data$day, levels = c('D0', 'D2','D7', 'D14', 'D21'))

FS12_RPS@sam_data$shed <- ifelse(FS12_RPS@sam_data$pignum %in% c(373,321,181,392,97), 'low', 'high')
FS12_RPS@sam_data$shed <- factor(FS12_RPS@sam_data$shed, levels = c('high', 'low'))


FS12_RPS@sam_data$set <- paste(FS12_RPS@sam_data$set, FS12_RPS@sam_data$shed, sep = '_')
#
FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~day*shed)
#
# Keeping samples separate by day #

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D0' & tissue == 'F')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
library(DESeq2)
FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)

tmpres <- results(FS12.de, name = 'shed_low_vs_high', cooksCutoff = FALSE)
tmpres <- lfcShrink(FS12.de, res=tmpres, coef = 'shed_low_vs_high', type = 'apeglm')
tmpres[tmpres$padj < 0.1,]

D0_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'apeglm', cookscut = FALSE)

### D0 Q

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D0' & tissue == 'Q')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)
library(DESeq2)
FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)

D0_Q_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                                phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                                name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'apeglm', cookscut = FALSE)


##### D2
FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D2')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D2_highlow <-Deseq.quickplot(DeSeq.object = FS12.de,
                             phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                             name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'apeglm', cookscut = FALSE)

#### D7
FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D7')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D7_highlow <-Deseq.quickplot(DeSeq.object = FS12.de,
                             phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                             name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'apeglm', cookscut = FALSE)


# D14 #

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D14')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D14_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                               phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                               name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'apeglm', cookscut = FALSE)

##### D21 F

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'F')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D21F_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                                phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                                name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'apeglm', cookscut = FALSE)

#### Tissue X


FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'X')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D21X_highlow <-Deseq.quickplot(DeSeq.object = FS12.de,
                               phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                               name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'apeglm', cookscut = FALSE)


##### tissue C



FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'C')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D21C_highlow <-Deseq.quickplot(DeSeq.object = FS12.de,
                               phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                               name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)

##### tissue I

FS12b.glom <- FS12_RPS
FS12b.glom <- subset_samples(FS12b.glom, day =='D21' & tissue == 'I')

FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)

FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~shed)

FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
resultsNames(FS12.de)


D21I_highlow <- Deseq.quickplot(DeSeq.object = FS12.de,
                                phyloseq.object = FS12b.glom, pvalue = .05, alpha = 0.05,
                                name = 'shed_low_vs_high' ,taxlabel = 'Genus', shrink_type = 'normal', cookscut = FALSE)



apeglm::apeglm()

#
D0_highlow[[2]]$set <- 'D0_feces'
D2_highlow[[2]]$set <- 'D2_feces'
D7_highlow[[2]]$set <- 'D7_feces'
D14_highlow[[2]]$set <- 'D14_feces'
D21F_highlow[[2]]$set <- 'D21_feces'
D21C_highlow[[2]]$set <- 'D21_cecal_content'
D21X_highlow[[2]]$set <- 'D21_cecal_mucosa'
D21I_highlow[[2]]$set <- 'D21_ileal_mucosa'
#

class(D0_highlow[[2]])


RPS_split_master <- bind_rows(list(D0_highlow[[2]],
                                   D2_highlow[[2]],
                                   D7_highlow[[2]],
                                   D14_highlow[[2]],
                                   D21F_highlow[[2]],
                                   D21C_highlow[[2]],
                                   D21X_highlow[[2]], 
                                   D21I_highlow[[2]]))



RPS_split_master$imp <- ifelse(RPS_split_master$padj <= 0.05, TRUE, FALSE)

RPS_split_master$set <- factor(RPS_split_master$set, levels = c('D0_feces','D2_feces' ,'D7_feces', 'D14_feces', 'D21_feces', 'D21_cecal_content', 'D21_cecal_mucosa', 'D21_ileal_mucosa'))
RPS_split_master <- RPS_split_master %>% mutate(newp2=paste0('p=', newp))
library(ggscinames)
library(grid)

RPS_split_master %>% filter(set %in% c('D0_feces' ,'D7_feces', 'D14_feces', 'D21_feces')) %>%
  ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
  geom_bar(stat='identity') + 
  geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
  facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1')

RPS_split_master %>% filter(set %in% c('D21_feces', 'D21_cecal_content', 'D21_cecal_mucosa')) %>%
  ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
  geom_bar(stat='identity') + 
  geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
  facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1') + xlab('') + labs(fill='Shedding')

RPS_split_master %>% filter(set %in% c('D21_ileal_mucosa')) %>%
  ggplot(aes(x=reorder(OTU, log2FoldChange), y=log2FoldChange, fill=Treatment)) +
  geom_bar(stat='identity') + 
  geom_text_sciname(aes(x=OTU, y=0, sci = Genus, nonsci=newp2, important=imp), size=3) + coord_flip() +
  facet_wrap(~set, ncol = 1, scales = 'free_y') + scale_fill_brewer(palette = 'Set1') + xlab('') + labs(fill='Shedding')




library(cowplot)

p <- RPS_split_master %>%
  group_by(OTU, Treatment) %>%
  filter(padj <= 0.05) %>%  tally() %>%
  ggplot(aes(x=OTU, y=n, fill=Treatment)) + geom_col() +
  scale_fill_brewer(palette = 'Pastel1') + ylab('occurences') + ggtitle('Number of times OTUs are significantly enriched (p<0.05)\n in either shedding phenotype') + 
  theme(axis.text.x = element_text(angle = 90, vjust = .4)) + xlab('')

 
ggplot2::ggplot_build(p)



c("#E41A1C", "#FBB4AE", "#377EB8", "#B3CDE3")



# # MERGE THIS WITH TAX AND PRINT TABLE
# RPS_split_master %>%
#   group_by(OTU, Treatment) %>%
#   filter(padj <= 0.1) %>%  tally()


RPS_split_master <- RPS_split_master %>% mutate(p2=ifelse(padj <= 0.05, 'p < 0.05',
                                      ifelse(padj <= 0.1, 'p < 0.1', NA)))
RPS_split_master <- RPS_split_master %>% mutate(group=factor(paste(p2, Treatment), levels=c("p < 0.05 high", 
                                                                                            "p < 0.1 high", 
                                                                                            "p < 0.05 low", 
                                                                                            "p < 0.1 low")))


RPS_split_master %>% group_by(OTU, group) %>% tally() %>% 
  ggplot(aes(x=OTU, y=n, fill=group)) + geom_col(color='black') +
  scale_fill_manual(values = c("#E41A1C", "#FBB4AE", "#377EB8", "#B3CDE3")) +
  ylab('occurences') +
  ggtitle('Number of times OTUs are enriched \n in either RPS shedding phenotype') + 
  theme(axis.text.x = element_text(angle = 90, vjust = .4)) + xlab('')



int_OTUs <- RPS_split_master %>% group_by(OTU, group) %>% tally() %>% filter(n>1) %>% select(OTU) %>% unlist(use.names = FALSE)

write_csv(RPS_split_ints, 'RPS_split_int_OTUs.csv')
RPS_split_ints <- RPS_split_master %>% filter(OTU %in% int_OTUs) %>% 
  select(OTU, Treatment, Genus) %>% unique()

tax <- as.data.frame(FS12b.glom@tax_table)
tax$OTU <- rownames(tax)



#####


#############################################

# should make a function for this....
# it would take timepoint, tissue, and return the sig diff OTUs in a dataframe
# need to add tissue and timepoint to dataframe before return


unique(FS12b@sam_data$pignum)

FS12b@sam_data$treatment


DESeq_difabund <- function(phyloseq, day, tissue, scientific = TRUE, shrink_type='normal', 
                           alpha=0.1, cooks_cut=FALSE, pAdjustMethod='BH'){
  
  # FS12b.glom <- tax_glom(FS12b, taxrank = 'Genus')
  FS12b.glom <- prune_samples(x = phyloseq, samples = phyloseq@sam_data$day == day & phyloseq@sam_data$tissue == tissue)
  FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
  FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~treatment)
  FS12.de <- DESeq(FS12.de, test = 'Wald', fitType = 'parametric')
  
  finres <- list()
  resind <- 1
  for (i in 2:length(resultsNames(FS12.de))){
    print(resultsNames(FS12.de)[i])
    treat <- sub('treatment_(.*)_vs_Control','\\1',resultsNames(FS12.de)[i])
    comp <- sub('treatment_', '', resultsNames(FS12.de)[i])
    res <- results(object = FS12.de, name = resultsNames(FS12.de)[i], alpha=alpha, cooksCutoff = cooks_cut, pAdjustMethod = pAdjustMethod)
    res <- lfcShrink(FS12.de, coef = resultsNames(FS12.de)[i], type = shrink_type)
    sigtab = res[which(res$padj < alpha), ]
    
    if (nrow(sigtab) != 0){
      # browser()
      sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(FS12b.glom)[rownames(sigtab), ], "matrix"))
      sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = scientific)
      sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, treat, paste('down',treat, sep = '_'))
      sigtab$OTU <- rownames(sigtab)
      sigtab$tissue <- tissue
      sigtab$day <- day
      sigtab$comp <- comp
      finres[[resind]] <- sigtab
      
      resind <- resind + 1
    }
    
    
    
  }
  
  finres <- bind_rows(finres)
  return(finres)
  
}



tocont <- list(DESeq_difabund(phyloseq = FS12b, day = 'D0', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D0', tissue = 'Q', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D2', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D7', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D14', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D21', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D21', tissue = 'C', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D21', tissue = 'X', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'),
               DESeq_difabund(phyloseq = FS12b, day = 'D21', tissue = 'I', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = TRUE, pAdjustMethod = 'BH'))


tocont <- bind_rows(tocont)

tocontf <- tocont[abs(tocont$log2FoldChange) > .75,]

tocontf %>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + geom_hline(yintercept = 20, color='red', size=3)


biguns <- tocontf %>% group_by(OTU) %>% summarise(tot=sum(log2FoldChange)) %>% filter(tot >20) %>% select(OTU) %>% unlist()

tocont %>% filter(OTU %in% biguns) %>% select(OTU,Genus) %>% unique()

tocontf %>% group_by(OTU, Treatment) %>% tally() %>% filter(n>3) %>% as.data.frame()
#### Ok that wasnt so bad.

# Now, which OTUs changed at Salmonella infection?
#ALL FECES
# D0 vs D2 within treatments
# D0 vs D7 within treatments
# D0 vs D14 within treatments
# D0 vs D21 within treatments
unique(FS12b@sam_data$pignum)

FS12b@sam_data$day <- factor(FS12b@sam_data$day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))
FS12b@sam_data$pignum <- factor(FS12b@sam_data$pignum)


FS12b.glom <- prune_samples(x = FS12b, samples = FS12b@sam_data$treatment == 'Control' & FS12b@sam_data$tissue == 'F')
FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
FS12.de <- phyloseq_to_deseq2(FS12b.glom, ~pignum + day)
FS12.de <- DESeq(FS12.de, test = 'LRT', reduced = ~ pignum)

resultsNames(FS12.de)

test2 <- results(object = FS12.de, name = 'day_D2_vs_D0')
test7 <- results(object = FS12.de, name = 'day_D7_vs_D0')

sigtab2 <- test2[which(test2$padj < 0.1),]
sigtab7 <- test7[which(test7$padj < 0.1),]

all(rownames(sigtab2) == rownames(sigtab7))


sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(FS12b.glom)[rownames(sigtab), ], "matrix"))
sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = scientific)
sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, treat, 'Control')
sigtab$OTU <- rownames(sigtab)
sigtab$tissue <- tissue
sigtab$day <- day
sigtab$comp <- comp












# FS12b@sam_data$tissue



C_ss_1 <- taxa[grep('Clostridium_sensu_stricto_1', taxa$Genus),]
write.csv(C_ss_1, 'C_ss_1.csv')
################### pig trips ##############


min(rowSums(FS12b@otu_table))


test <- data.frame(FS12b@otu_table)
rownames(test)
rowSums(test)

FS12.otu.rare <- rrarefy(test, min(rowSums(test)))



bray.dist <- vegdist(FS12.otu.rare, method = 'jaccard', binary = FALSE)

FS12b_meta$sample_ID == rownames(FS12.otu.rare)

#FS12b_meta <- data.frame(FS12b@sam_data)



#FS12b_meta$

FS12b_meta$shan2 <- diversity(FS12b@otu_table, base = 2)
# FS12b_meta$shan <- diversity(FS12b@otu_table)
# FS12b_meta$invsimp <- diversity(FS12b@otu_table, index = 'invsimpson')

#FS12b_meta$day <- factor(FS12b_meta$day, levels = c('D0', 'D2', 'D7', 'D14', 'D21'))
#FS12b_meta$treatment <- factor(FS12b_meta$treatment, levels = c('control', 'RPS', 'Acid', 'ZnCu', 'RCS', 'Bglu'))

# FS12b_meta %>% filter(tissue == 'F') %>% 
#   ggplot(aes(x=day, y=shan2, group=set, fill=treatment)) +
#   geom_boxplot() + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + geom_text(aes(label=pignum))

# FS12b_meta %>% filter(tissue == 'F') %>% 
#   ggplot(aes(x=day, y=invsimp, group=set, fill=treatment)) +
#   geom_boxplot() + scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + geom_text(aes(label=pignum))
# 



# FS12b_meta %>% filter(day == 'D21') %>% ggplot(aes(x=tissue, y=shan2, group=set, fill=treatment)) + geom_boxplot() + geom_text(aes(label=pignum))


# FS12b_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=day, y=shan, group=set, fill=treatment)) + geom_boxplot()
# FS12b_meta %>% filter(day == 'D21') %>% ggplot(aes(x=tissue, y=shan, group=set, fill=treatment)) + geom_boxplot()




#ggplot(FS12b_meta, aes(x=treatment, y=shan, group=set)) + geom_boxplot()

# min(rowSums(shareds_test))
# hist(rowSums(shareds_test))
# sort(rowSums(shareds_test))

dist.data <- as.data.frame(as.matrix(bray.dist))
dist.data$from <- rownames(dist.data)

dist.gather <- gather(data = dist.data, key = 'to', value = 'distance', -from)

#


dist.gather$fromPig <- gsub('([X12]+[ab]?)([NP]+)([0-9]+)([dWXDi]+)([0-9]+)([A-Z]?)', '\\3', dist.gather$from)

#

dist.gather$FT <- paste(dist.gather$from, dist.gather$to, sep = ' ')



#dist.gather$TF <- paste(dist.gather$to, dist.gather$from, sep = ' ')

######
# all pig pairwise #

total_ground_covered <- dist.gather[grep('X12bP([0-9]+)D[0-9]+F X12bP\\1D[0-9]+F', dist.gather$FT),] %>% group_by(fromPig) %>% summarise(allpw=sum(distance),
                                                                                                                                         num=n())

rooms <- read.csv('Rooms.csv')
total_ground_covered$treatment <- ifelse(total_ground_covered$fromPig %in% rooms$X6, 'control',
                                         ifelse(total_ground_covered$fromPig %in% rooms$X7, 'RPS', 
                                                ifelse(total_ground_covered$fromPig %in% rooms$X8, 'Acid', 
                                                       ifelse(total_ground_covered$fromPig %in% rooms$X9, 'Zn+Cu',
                                                              ifelse(total_ground_covered$fromPig %in% rooms$X10, 'RCS',
                                                                     ifelse(total_ground_covered$fromPig %in% rooms$X11, 'Bglu', 'asdfsa'))))))




sum_sal$fromPig <- sum_sal$pignum
total_ground_covered$fromPig

total_ground_covered <- total_ground_covered %>% filter(num == 25)

boxplot(total_ground_covered$allpw~total_ground_covered$treatment)

sum_sal
total_ground_covered <- merge(total_ground_covered, sum_sal, by = 'fromPig')

############### NEED TO READ IN SUM_SAL ################
########################################################

total_ground_covered$treatment.y == total_ground_covered$treatment.x
total_ground_covered <- total_ground_covered %>% mutate(treatment=treatment.x) %>% select(-treatment.x, -treatment.y)

#cor.test(total_ground_covered$allpw, total_ground_covered$AULC)


total_ground_covered %>% group_by(treatment) %>% summarise(AULCvTRIP_P=cor.test(AULC, allpw, method = 'pearson')$p.value,
                                                           AULCvTRIP_T=cor.test(AULC, allpw, method = 'pearson')$statistic)



total_ground_covered %>% 
  ggplot(aes(x=allpw, y=AULC, fill=treatment, color=treatment)) +
  geom_point(size=2, shape=21) + geom_smooth(method = 'lm', se=FALSE) +
  ggtitle('Correlation between cumulative community membership change and cumulative shedding', 
          subtitle = 'correlation stats: RPS pval = 0.02, control pval = 0.31, Bglu pval = 0.42') +
  xlab('Cumulative Bray-Curtis distance (presence/abscence)')





#
######
D0_2 <- dist.gather[grep('X12bP([0-9]+)D0F X12bP\\1D2F', dist.gather$FT),]
#colnames(D0_2)[1] <- 'sample_ID'
colnames(D0_2)[3] <- 'D0_2'
D0_2 <- D0_2[,c(3,4)]

D2_7 <- dist.gather[grep('X12bP([0-9]+)D2F X12bP\\1D7F', dist.gather$FT),]
#colnames(D2_7)[1] <- 'sample_ID'
colnames(D2_7)[3] <- 'D2_7'
D2_7 <- D2_7[,c(3,4)]

D7_14 <- dist.gather[grep('X12bP([0-9]+)D7F X12bP\\1D14F', dist.gather$FT),]
#colnames(D7_14)[1] <- 'sample_ID'
colnames(D7_14)[3] <- 'D7_14'
D7_14 <- D7_14[,c(3,4)]

D14_21 <- dist.gather[grep('X12bP([0-9]+)D14F X12bP\\1D21F', dist.gather$FT),]
#colnames(D14_21)[1] <- 'sample_ID'
colnames(D14_21)[3] <- 'D14_21'
D14_21 <- D14_21[,c(3,4)]

D0_21 <- dist.gather[grep('X12bP([0-9]+)D0F X12bP\\1D21F', dist.gather$FT),]
#colnames(D14_21)[1] <- 'sample_ID'
colnames(D0_21)[3] <- 'D0_21'
D0_21 <- D0_21[,c(3,4)]


#full_join(D0_2, D2_7)
pig_trips <- merge(D0_2, D2_7, all = TRUE, by = 'fromPig')
pig_trips <- merge(pig_trips, D7_14, all = TRUE, by = 'fromPig')
pig_trips <- merge(pig_trips, D14_21, all = TRUE, by = 'fromPig')
pig_trips <- merge(pig_trips, D0_21, all = TRUE, by = 'fromPig')

#rowSums(pig_trips)
#pig_trips <- na.omit(pig_trips)
colnames(pig_trips[,c(2:5)])
pig_trips$trip <- rowSums(pig_trips[,c(2:5)])
hist(pig_trips$trip, breaks = 10)



rooms <- read.csv('../FS12/Rooms.csv')

# add treatment data.  This probably isn't the best way to do this...

#colnames(sum_sal)[1] <- 'fromPig'

library(funfuns)

#NMDS_ellipse(metadata = meta_test, OTU_table = shareds_test, grouping_set = 'pig_pen')
###############################

# pig_trips$treatment <- ifelse(pig_trips$fromPig %in% rooms$X6, 'control',
#                               ifelse(pig_trips$fromPig %in% rooms$X7, 'RPS', 
#                                      ifelse(pig_trips$fromPig %in% rooms$X8, 'Acid', 
#                                             ifelse(pig_trips$fromPig %in% rooms$X9, 'Zn+Cu',
#                                                    ifelse(pig_trips$fromPig %in% rooms$X10, 'RCS',
#                                                           ifelse(pig_trips$fromPig %in% rooms$X11, 'Bglu', 'asdfsa'))))))
# 
# 

pig_trips <- merge(pig_trips, sum_sal, by = 'fromPig')

boxplot(pig_trips$trip~pig_trips$treatment)
boxplot(pig_trips$D0_21~pig_trips$treatment)

pairwise.wilcox.test(x=pig_trips$trip, g=pig_trips$treatment, p.adjust.method = 'none')

#colnames(sum_sal)[1] <- 'fromPig'
pig_trips %>% filter(treatment == "RPS") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

pig_trips %>% filter(treatment == "control") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

pig_trips %>% 
  ggplot(aes(x=treatment, y=trip, fill=treatment)) + 
  geom_boxplot() +
  geom_jitter(shape=21, color='black', stroke=1.2, size=2, width = .2) +
  scale_fill_manual(values = c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  ggtitle('Cumulative change in each individual pigs community stucture over 21 days') + ylab("Cumulative Jaccard distance")




pig_trips_cor <- pig_trips %>% group_by(treatment) %>% summarise(AULCvTRIP_P=cor.test(AULC, trip)$p.value,
                                                AULCvTRIP_T=cor.test(AULC, trip)$statistic,
                                                AULCv02_P=cor.test(AULC, D0_2)$p.value,
                                                AULCv02_T=cor.test(AULC, D0_2)$statistic,
                                                AULCv27_P=cor.test(AULC, D2_7)$p.value,
                                                AULCv27_T=cor.test(AULC, D2_7)$statistic,
                                                AULCv714_P=cor.test(AULC, D7_14)$p.value,
                                                AULCv714_T=cor.test(AULC, D7_14)$statistic,
                                                AULCv1421_P=cor.test(AULC, D14_21)$p.value,
                                                AULCv1421_T=cor.test(AULC, D14_21)$statistic)



# pig_trips %>% group_by(treatment) %>% summarise(AULCvTRIP_P=cor.test(AULC, D0_21)$p.value,
#                                                 AULCvTRIP_T=cor.test(AULC, D0_21)$statistic)


pig_trips %>% filter(treatment == "control") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
pig_trips %>% filter(treatment == "RPS") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
pig_trips %>% filter(treatment == "Bglu") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
pig_trips %>% filter(treatment == "Zn+Cu") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

pig_trips %>% 
  ggplot(aes(x=trip, y=AULC, fill=treatment, color=treatment)) +
  geom_point(size=3, shape=21, color='black') + geom_smooth(method = 'lm', se=FALSE, size=2) +
  ggtitle('Correlation between cumulative community change and cumulative shedding', 
          subtitle = 'correlation stats: RPS pval = 0.037, control pval = 0.09, Bglu pval = 0.08') +
  xlab('Cumulative Jaccard distance') + scale_color_manual(values = c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + 
  scale_fill_manual(values = c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))


c('#3399FF', 'orange', 'red', 'grey', 'purple')
c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')

#testse <- cor.test(pig_trips$trip, pig_trips$AULC)
#testse$p.value
#testse$statistic


ggplot(pig_trips, aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
ggplot(pig_trips, aes(x=D0_2, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
ggplot(pig_trips, aes(x=D2_7, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
ggplot(pig_trips, aes(x=D7_14, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
ggplot(pig_trips, aes(x=D14_21, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

ggplot(pig_trips_test, aes(x=mean_trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm', fill = NA) + geom_text(aes(label=pignum))

#ggplot(pig_trips, aes(x=sum, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

pig_trips %>% group_by(treatment) %>% summarise(num=n())

apply(X = pig_trips, MARGIN = 2, FUN = mean, na.rm=TRUE)

mean(pig_trips$D0_2, na.rm=TRUE)
mean(pig_trips$D2_7, na.rm=TRUE)
mean(pig_trips$D7_14, na.rm=TRUE)
mean(pig_trips$D14_21, na.rm=TRUE)

median(pig_trips$D0_2, na.rm=TRUE)
median(pig_trips$D2_7, na.rm=TRUE)
median(pig_trips$D7_14, na.rm=TRUE)
median(pig_trips$D14_21, na.rm=TRUE)


# looking for missing samples


sum(shared_table[grep('P50D0', rownames(shared_table)),])
sum(shared_table[grep('P181D7F', rownames(shared_table)),])



ggplot(pig_trips, aes(x=treatment, y=trip, fill=treatment)) +
  geom_boxplot() + ylab('Cumulative bray-curtis dissimilarity (each pig)') + geom_jitter(size=2.5,width = 0.2, shape=21)+
  ggtitle('Cumulative change in community structure through Salmonella infection')


########### cor stuff  ###########

# D0 correlations
fec_VFAs <- res.all

fec_VFAs_0 <- fec_VFAs %>% filter(time == 0) %>% mutate(day=time) %>% select(day, everything(),-time)

ttttt <- FS12b_meta %>% group_by(day) %>% nest()

FS12b_meta %>% group_by(day) %>% nest()

colnames(ttttt$data[[1]])





col_nams_map <- function(df){
  colnames(df) <- paste(day)
}
map()


get_pairs <- function(df){
  pp <- pairwise.wilcox.test(df$dispers.distances, df$treatment, p.adjust.method = 'none')
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <- paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps
  
}


shan_fecal_tests <- FS12b_meta %>% filter(tissue =='F') %>% group_by(day) %>% 
  nest() %>% mutate(pps = map(data, get_pairs)) %>%
  select(day, pps) %>% unnest() %>% select(day, starts_with('control'))


########## MISSING DATA??  #########



tttt <- FS12b_meta %>%filter(tissue =='F') %>%  group_by(pignum, day) %>% tally() %>% spread(key = day, value = n)

tttt <- FS12b_meta %>% select(pignum, treatment) %>% unique() %>% left_join(tttt, by = 'pignum')

pig_trips %>% ggplot(aes(x=D0_2, y = D2_7)) + geom_point(aes(color = treatment),size=3) + geom_point()
pig_trips %>% ggplot(aes(x=D0_2, y = D7_14)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D0_2, y = D14_21)) + geom_point(aes(color = treatment),size=3)

pig_trips %>% ggplot(aes(x=D2_7, y = D0_2)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D2_7, y = D7_14)) + geom_point(aes(color = treatment),size=3) + geom_smooth(method = 'lm')
pig_trips %>% ggplot(aes(x=D2_7, y = D14_21)) + geom_point(aes(color = treatment),size=3)

pig_trips %>% ggplot(aes(x=D7_14, y = D0_2)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D7_14, y = D2_7)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D7_14, y = D14_21)) + geom_point(aes(color = treatment),size=3)

pig_trips %>% ggplot(aes(x=D14_21, y = D0_2)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D14_21, y = D2_7)) + geom_point(aes(color = treatment),size=3)
pig_trips %>% ggplot(aes(x=D14_21, y = D7_14)) + geom_point(aes(color = treatment),size=3)


pig_trips$missing <- ifelse(pig_trips$pignum %in% c(50,181,211,240,253,469), TRUE, FALSE)
pig_trips_test <- pig_trips[,c(1:5, 7)]
PTgath <- pig_trips_test %>% gather(key = interval, value = distance, -fromPig)

avtrp <- PTgath %>% group_by(fromPig) %>% summarise(mean_trip = mean(distance, na.rm = TRUE))

pig_trips_test <- merge(pig_trips, avtrp, by = 'fromPig')

pig_trips_test %>% ggplot(aes(x=trip, y=mean_trip)) + geom_point()

#########

phyloseq::transform_sample_counts()


phyloseq::transform_sample_counts()

FS12_RPS <- subset_taxa(FS12_RPS, taxa_sums(FS12_RPS) > 1)

# plot_bar(FS12_RPS, x='shed')


FS12_RPS_sam <- as_data_frame(FS12_RPS@sam_data)

wht <- FS12_RPS_sam %>% group_by(pignum, tissue) %>% tally()
# missing 50 and 181 fecals

FS12_RPS_otu <- as.data.frame(FS12_RPS@otu_table)
FS12_RPS_otu <- FS12_RPS_otu/rowSums(FS12_RPS_otu) # transforms to relative abundance
FS12_RPS_tax <- as.data.frame(FS12_RPS@tax_table)
FS12_RPS_tax$OTU <- rownames(FS12_RPS_tax)
#
colSums(FS12_RPS_otu)
colsums97 <- colSums(FS12_RPS_otu[grep(97, rownames(FS12_RPS_otu)),])/nrow(FS12_RPS_otu[grep(97, rownames(FS12_RPS_otu)),])
colsums_others <- colSums(FS12_RPS_otu[grep(97, rownames(FS12_RPS_otu), invert=TRUE),])/nrow(FS12_RPS_otu[grep(97, rownames(FS12_RPS_otu), invert=TRUE),])
lowerin97 <- (colsums97 - colsums_others) < 0
higherin97 <- (colsums97 - colsums_others) > 0
#

FS12_RPS_otu$sample_ID <- rownames(FS12_RPS_otu)



FS12_RPS_all <- merge(FS12_RPS_sam, FS12_RPS_otu, by='sample_ID')

FS12_RPS_all[1:10, 1:10]

FS12_gath <- FS12_RPS_all %>% gather(key=OTU, value=relabund, -(sample_ID:shed))
FS12_RPS_tax


FS12_gath %>% ggplot(aes(x=pignum, y=relabund)) + geom_col()
