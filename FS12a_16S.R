getwd()

library(phyloseq)
library(tidyverse)
library(funfuns)
library(pairwiseAdonis)

meta <- read.csv('./data/FS12_final_meta.csv', header = TRUE, stringsAsFactors = FALSE)
shared <- read_delim('./data/FS12.shared', delim = '\t') %>% as.data.frame()

taxa <- extract_mothur_tax('./data/FS12.taxonomy')
# taxa[1:10,1:8]

rownames(shared) <- shared$Group
shared <- shared[,-c(1,2,3)]

hist(rowSums(shared), breaks = 50)


# length(colnames(shared)) # 10648 total OTUs detected (no removal mocks and NTCs and all kinds of garbage here)


# rownames(shared) %in% meta$sample_ID
meta <- meta[meta$sample_ID %in% rownames(shared),]
mocks <- shared[(!rownames(shared) %in% meta$sample_ID),]
shared <- shared[rownames(shared) %in% meta$sample_ID,]

rownames(shared) == meta$sample_ID

rownames(meta) <- meta$sample_ID

shared <- shared[rowSums(shared) > 1250,] # remove samples with less than 1250 reads
shared <- shared[,colSums(shared > 0) > 3] # removes otus that are detected in fewer than 4 samples globally (including timecourse data)
shared <- shared[,colSums(shared) > 5] # at least 6 observations globally

length(colnames(shared))
hist(rowSums(shared), breaks = 50)


meta <- meta[order(meta$sample_ID),]
shared <- shared[order(rownames(shared)),]

rownames(shared) == meta$sample_ID

meta$set <- paste(meta$experiment, meta$day, meta$tissue, meta$treatment, sep = '_')

nnnn <- meta %>% group_by(set) %>% summarise(N=n())


rownames(meta) == rownames(shared)
meta <- meta[rownames(meta) %in% rownames(shared),]
rownames(shared) %in% rownames(meta)
shared <- shared[match(rownames(meta), rownames(shared)),]

rownames(meta) == rownames(shared)

# taxa <- taxa[taxa$OTU %in% colnames(shared),]

p_meta <- sample_data(meta) 
p_taxa <- import_mothur(mothur_constaxonomy_file = './data/FS12.taxonomy')


colnames(p_taxa) <- colnames(taxa)[-c(1,2,3)]

# meta$treatment  
FS12 <- phyloseq(p_meta, p_taxa, otu_table(shared, taxa_are_rows = FALSE))  # builds the phyloseq obj

### mock stuff ###

rownames(mocks)
mocks[1:5,1:5]
mocks <- mocks[rowSums(mocks) >1000,]
rownames(mocks)
hist(rowSums(mocks), breaks = 50)

mocks <- mocks[grep('NTC', rownames(mocks)),]

mocks <- mocks[,colSums(mocks)>2]
colSums(mocks)
rownames(mocks)
mock_tax <- taxa[taxa$OTU %in% colnames(mocks),]






FS12a <- FS12 %>% subset_samples(experiment == 'X12a')

FS12a@sam_data$treatment

# 
# FS12a@sam_data$treatment <- factor(FS12a@sam_data$treatment, levels = c('Control', 
#                                                                         'RPS',
#                                                                         'Acid', 
#                                                                         'ZnCu', 
#                                                                         'RCS', 
#                                                                         'Bglu', 
#                                                                         'Pos_control', 
#                                                                         'Phyto', 
#                                                                         'Malto', 
#                                                                         'SBP'))

hist(sample_sums(FS12a), breaks = 50)


FS12a@sam_data %>% filter(day == 'D30' & treatment == 'Phyto')

# FS12a <- FS12a %>% subset_samples(day %in% c('D0', 'D23', 'D30'))
FS12a <- FS12a %>% subset_samples(day %in% c('D0', 'D23'))


rems <- !(FS12a@sam_data$sample_ID %in% c('X12N64WD0dup', 'X12N77WD0dup', 'X12N70XD30'))
FS12a <- prune_samples(x = FS12a, rems)

FS12a@sam_data %>% group_by(day, treatment) %>% tally()




# FS12a@sam_data %>% filter(day == 'D23' & treatment == 'Control')


keeps_at_d23 <- FS12a@sam_data %>% filter(day == 'D0' & treatment == 'Control') %>% select(pignum) %>% unique() %>% unlist()
REMOVE_THESE <- FS12a@sam_data %>% filter(day == 'D23' & treatment == 'Control') %>% filter(!(pignum %in% keeps_at_d23)) %>% select(sample_ID) %>% unlist()
ANDTHESE <- c('X12N16WD23', 'X12N18WD23', 'X12N32WD23', 'X12N39WD23', 'X12N45WD23', 'X12N68WD23')
REMOVE_THESE <- c(REMOVE_THESE, ANDTHESE)

FS12a <- FS12a %>% subset_samples(!(sample_ID %in% REMOVE_THESE))
FS12a@sam_data %>% group_by(day, treatment) %>% tally()

#### PERMANOVAS VS CONTROL

FS12a@sam_data %>% group_by(day, tissue) %>% mutate(minseqs=min(num_seqs)) %>% select(day, tissue, minseqs) %>% unique()


FS12a@sam_data$set <- paste(FS12a@sam_data$day, FS12a@sam_data$tissue, FS12a@sam_data$treatment, sep='_')



FS12a_rare <- phyloseq::rarefy_even_depth(FS12a,sample.size =  min(phyloseq::sample_sums(FS12a)), rngseed = 7)


FS12a_NMDS <- NMDS_ellipse(FS12a_rare@sam_data, OTU_table = FS12a_rare@otu_table, grouping_set = 'set', distance_method = 'bray', MDS_trymax = 150)

labbies <- FS12a_NMDS[[1]] %>% select(tissue, treatment, day, set, centroidX, centroidY) %>% unique()
# library(ggrepel)

#####
library(RColorBrewer)
colourCount = length(unique(FS12a@sam_data$treatment))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

fills <- c('black', getPalette(9))

unique(FS12a_NMDS[[1]]$treatment)


FS12a_NMDS[[1]]$treatment <- factor(FS12a_NMDS[[1]]$treatment, levels = c('Control', 'Acid', 'Bglu', 'Malto', 'Phyto', 'Pos_control', 'RCS', 'RPS', 'SBP', 'ZnCu'))
# getPal

labbies$treatment <- factor(labbies$treatment, levels =c('Control', 'Acid', 'Bglu', 'Malto', 'Phyto', 'Pos_control', 'RCS', 'RPS', 'SBP', 'ZnCu'))
### THIS IS A GOOD ONE. SHOW WITH GLOBAL ADONIS TEST

FS12a_NMDS[[1]] %>% ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=day), size=2)+
  geom_segment(aes(xend=centroidX, yend=centroidY, color=day)) +
  geom_point(data=labbies,aes(x=centroidX, y=centroidY, fill=treatment), size=3, inherit.aes = FALSE, shape=21)+# scale_shape_manual(values = c(21, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + scale_fill_manual(values = fills)
# geom_text_repel(data=labbies, aes(label=treatment, x=centroidX, y=centroidY))#+ geom_path(data = FS12a_NMDS[[2]], aes(x=NMDS1, y=NMDS2, group=group))

set.seed(7)
all_pwad <- pairwise.adonis(FS12a_rare@otu_table, FS12a_rare@sam_data$set, perm = 9999, sim.method = 'bray')

pwad_to_cont <- all_pwad[grep('Control', all_pwad$pairs),]
# same_day <- pwad_to_cont[grep('.*_(.*)_.*_.* vs .*_\\1_.*_.*', pwad_to_cont$pairs),]
same_day_tissue <- pwad_to_cont[grep('(.*)_(.*)_.* vs \\1_\\2_.*', pwad_to_cont$pairs),]
same_day_tissue$treatment <- sub('D[0-9]+_[FX]_([A-Za-z_]+) vs .*_.*_.*', '\\1',same_day_tissue$pairs)

same_day_tissue[same_day_tissue$treatment == 'Control',]$treatment <- sub('.*_.*_.* vs .*_.*_(.*)','\\1', same_day_tissue[same_day_tissue$treatment == 'Control',]$pairs)
# sub('(D[0-9]+)_([FX])_([A-Za-z_]+) vs .*_.*_.*', '\\2',same_day_tissue$pairs)
same_day_tissue$tissue <- sub('(D[0-9]+)_([FX])_([A-Za-z_]+) vs .*_.*_.*', '\\2',same_day_tissue$pairs)
same_day_tissue$day <- sub('(D[0-9]+)_([FX])_([A-Za-z_]+) vs .*_.*_.*', '\\1',same_day_tissue$pairs)

same_day_tissue$tissue <- ifelse(same_day_tissue$tissue == 'F', 'feces', 'cecal_mucosa')

# same_day_tissue$p.adjusted <- p.adjust(same_day_tissue$p.value, method = 'fdr')

same_day_tissue <- same_day_tissue %>% mutate(set=paste(tissue, day, sep = '_'))
same_day_tissue <- same_day_tissue %>% group_by(set) %>% mutate(p.adjusted = p.adjust(p.value, method = 'fdr'))
same_day_tissue$p.plot <- ifelse(same_day_tissue$p.adjusted <= 0.05, paste('p=', round(same_day_tissue$p.adjusted, 3), sep = ''), NA)

same_day_tissue$set <- factor(same_day_tissue$set, levels = c('feces_D0', 'feces_D23', 'feces_D30', 'cecal_mucosa_D30'))

#### NEED TO COME UP WITH COLOR TO USE
# 


same_day_tissue %>%
  ggplot(aes(x=treatment, y=F.Model, fill=treatment)) + 
  geom_col(color='black') + facet_wrap(~set) + geom_text(aes(label=p.plot), nudge_y = .2) + 
  ggtitle('Community differences compared to controls', subtitle = 'larger F = greater difference.  pvalues shown when <0.05') + 
  scale_fill_brewer(palette = 'Set1') + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))



# global adonis



FS12a_D0 <- FS12a_rare %>% subset_samples(day == 'D0')
FS12a_D23 <- FS12a_rare %>% subset_samples(day == 'D23')


global_adon <- adonis(FS12a_rare@otu_table~day + treatment,
                      data = data.frame(FS12a_rare@sam_data), method = 'bray')

global_adon_D0 <- adonis(FS12a_D0@otu_table~treatment,
                         data = data.frame(FS12a_D0@sam_data), method = 'bray')

global_adon_D23 <- adonis(FS12a_D23@otu_table~treatment,
                          data = data.frame(FS12a_D23@sam_data), method = 'bray')









write <- global_adon$aov.tab
write$variable <- rownames(write)

write <- write %>% select(variable, everything())
# write_csv(write, 'GLOB_ADON')
#### Deseq comparisons
# D23 each diet vs control
# D23 each diet vs pos_control

# D30 each diet vs control
# D30 each diet vs pos_control

# this is the same function i wrote for doing all the FS12b treatment comps for each day and tissue
# 
# FS12a@sam_data$treatment <- factor(FS12a@sam_data$treatment, levels = c('Control', 
#                                                                         'RPS',
#                                                                         'Acid', 
#                                                                         'ZnCu', 
#                                                                         'RCS', 
#                                                                         'Bglu', 
#                                                                         'Pos_control', 
#                                                                         'Phyto', 
#                                                                         'Malto', 
#                                                                         'SBP'))
FS12a@sam_data$treatment <- factor(FS12a@sam_data$treatment, levels = c('Control', 'Acid', 'Bglu', 'Malto', 'Phyto', 'Pos_control', 'RCS', 'RPS', 'SBP', 'ZnCu'))

# 
# unique(FS12a@sam_data$treatment)
# 
# FS12a.glom <- FS12a %>% tax_glom(taxrank = 'Genus')

# this can stay right?
# something weird is happening here....  treatment factor levels being replaced by integers
# it happened after I loaded te Rtsne package.... went away when I restarted R

library(DESeq2)
tmpres <- list()
for (lev in c('Genus', 'Family', 'Order')){
  FS12a.glom <- FS12a %>% tax_glom(taxrank = lev)
  tocont <- list(DESeq_difabund(phyloseq = FS12a.glom, day = 'D0', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
                 DESeq_difabund(phyloseq = FS12a.glom, day = 'D23', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'))
  tmp <- bind_rows(tocont)
  tmp$level <- lev
  tmpres[[lev]] <- tmp
}

tmpres$Order

# should write out each df individually instead of combining....
bind_rows(tmpres)

# FS12a_D23@sam_data$treatment <- factor(FS12a_D23@sam_data$treatment, levels = c('Control', 'Acid', 'Bglu', 'Malto', 'Phyto', 'Pos_control', 'RCS', 'RPS', 'SBP', 'ZnCu'))


# checkme <- phyloseq_to_deseq2(FS12a_D23, ~ treatment)
# checkme$treatment
# checkme@colData$treatment
# test <- DESeq(checkme, test = 'Wald', fitType = 'parametric')
# resultsNames(test)
# DESeq_difabund(phyloseq = FS12a, day = 'D23', tissue = 'F')

# 
# tocont.genus <- list(DESeq_difabund(phyloseq = FS12a.glom, day = 'D0', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                DESeq_difabund(phyloseq = FS12a.glom, day = 'D23', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                DESeq_difabund(phyloseq = FS12a.glom, day = 'D30', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                DESeq_difabund(phyloseq = FS12a.glom, day = 'D30', tissue = 'X', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'))


tocont.otu <- list(DESeq_difabund(phyloseq = FS12a, day = 'D0', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
                   DESeq_difabund(phyloseq = FS12a, day = 'D23', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'))

tmpres$OTU <- bind_rows(tocont.otu)
###
# 
# tocont.genus.ape <- list(DESeq_difabund(phyloseq = FS12a.glom, day = 'D0', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                      DESeq_difabund(phyloseq = FS12a.glom, day = 'D23', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                      DESeq_difabund(phyloseq = FS12a.glom, day = 'D30', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                      DESeq_difabund(phyloseq = FS12a.glom, day = 'D30', tissue = 'X', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'))
# 
# 
# tocont.otu.ape <- list(DESeq_difabund(phyloseq = FS12a, day = 'D0', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                    DESeq_difabund(phyloseq = FS12a, day = 'D23', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                    DESeq_difabund(phyloseq = FS12a, day = 'D30', tissue = 'F', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
#                    DESeq_difabund(phyloseq = FS12a, day = 'D30', tissue = 'X', scientific = TRUE, shrink_type = 'apeglm',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'))
# 

# 
# 
# 
# tocont.otu <- bind_rows(tocont.otu)
# tocont.genus <- bind_rows(tocont.genus)


# tocont.otu.ape <- bind_rows(tocont.otu.ape)
# tocont.genus.ape <- bind_rows(tocont.genus.ape)



tocontf.otu <- tocont.otu[abs(tocont.otu$log2FoldChange) > .5,]
tocontf.genus <- tocont.genus[abs(tocont.genus$log2FoldChange) > .5,]

# 
# tocontf.otu.ape <- tocont.otu.ape[abs(tocont.otu.ape$log2FoldChange) > .5,]
# tocontf.genus.ape <- tocont.genus.ape[abs(tocont.genus.ape$log2FoldChange) > .5,]

# 
# tocontf %>% group_by(OTU, Treatment) %>% tally() %>% filter(n>3) %>% as.data.frame()
# 
# taxa[taxa$OTU=='Otu00141',]
# taxa[taxa$OTU=='Otu00154',]
# taxa[taxa$OTU=='Otu00168',]


# Meh, not happy with this really...

# tocontf$Treatment <- factor(tocontf$Treatment, levels = c('RPS',
#                                                           'Acid', 
#                                                           'ZnCu', 
#                                                           'RCS', 
#                                                           'Bglu', 
#                                                           'Pos_control', 
#                                                           'Phyto', 
#                                                           'Malto', 
#                                                           'SBP', 
#                                                           'down_RPS',
#                                                           'down_Acid', 
#                                                           'down_ZnCu', 
#                                                           'down_RCS', 
#                                                           'down_Bglu', 
#                                                           'down_Pos_control', 
#                                                           'down_Phyto', 
#                                                           'down_Malto', 
#                                                           'down_SBP'))
# 
# 
# col_set1 <- RColorBrewer::brewer.pal(9, 'Set1')
# col_set2 <- sapply(col_set1, lighten)
# 
# fin_colset <- c(col_set1, col_set2)

levs <- c('RPS','Acid', 'ZnCu', 'RCS', 'Bglu', 'Pos_control', 'Phyto', 
          'Malto', 'SBP', 'down_RPS','down_Acid', 'down_ZnCu', 'down_RCS', 
          'down_Bglu', 'down_Pos_control', 'down_Phyto', 'down_Malto', 'down_SBP')

pal <- getPalette(9)

pal_light <- lapply(X = pal, FUN = lighten) %>% unlist()

my_pal <- c(pal, pal_light)
names(my_pal) <- levs

# lighten(pal)
# names(fin_colset) <- levs

# D23 and D30 all tissues all treats
# tocontf.genus %>% filter(day !='D0') %>%  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
#   geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) + 
#   geom_text(aes(label=Genus, y=0))
# 
# tocontf.otu %>% filter(day !='D0') %>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
#   geom_text(aes(label=Genus, y=0))
# 
# 


tmpres$OTU %>% filter(day !='D0') %>%  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) +
  geom_text(aes(label=Genus, y=0))

tmpres$Genus %>% filter(day !='D0') %>%  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) +
  geom_text(aes(label=Genus, y=0))


tmpres$Family %>% filter(day !='D0') %>%  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) +
  geom_text(aes(label=Family, y=0))


tmpres$Order %>% filter(day !='D0') %>%  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) +
  geom_text(aes(label=Order, y=0))
# 
# tocontf.otu %>% filter(day !='D0') %>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
#   geom_text(aes(label=Genus, y=0))
# 
# 


#### genera for each treatment at D23 feces

# tocontf.genus %>% filter(day =='D23' & grepl('Acid', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0))
# 
# tocontf.genus %>% filter(day =='D23' & grepl('RCS', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0))
# 
# tocontf.genus %>% filter(day =='D23' & grepl('Bglu', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0))

tocontf.genus %>% filter(day =='D23' & grepl('Malto', Treatment))%>% ggplot(aes(x=Genus, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0)) + theme(axis.text.y = element_blank())

tocontf.genus %>% filter(day =='D23' & grepl('RPS', Treatment))%>% ggplot(aes(x=Genus, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0)) + theme(axis.text.y = element_blank())

tocontf.genus %>% filter(day =='D23' & grepl('Pos_control', Treatment))%>% ggplot(aes(x=Genus, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))+ theme(axis.text.y = element_blank()) + ggtitle('D23 Pos_control')

# tocontf.genus %>% filter(day =='D23' & grepl('ZnCu', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0))

tocontf.genus %>% filter(day =='D23' & grepl('Phyto', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0)) + ylim(-3,3)

tocontf.genus %>% filter(day =='D23' & grepl('SBP', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
  geom_text(aes(label=Genus, y=0))

####### Genera D30

tocontf.genus %>% filter(day =='D30' & grepl('Acid', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))

tocontf.genus %>% filter(day =='D30' & grepl('RCS', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))

# tocontf.genus %>% filter(day =='D30' & grepl('Bglu', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0))
# 
# tocontf.genus %>% filter(day =='D30' & grepl('Malto', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0))

tocontf.genus %>% filter(day =='D30' & grepl('RPS', Treatment))%>% ggplot(aes(x=Genus, y=log2FoldChange, fill=tissue)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0)) + ylim(-5,15)+ theme(axis.text.y = element_blank()) + ggtitle('D30 RPS')

tocontf.genus %>% filter(day =='D30' & grepl('Pos_control', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))+ theme(axis.text.y = element_blank()) + ggtitle('D30 Pos_control')

tocontf.genus %>% filter(day =='D30' & grepl('ZnCu', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))

# tocontf.genus %>% filter(day =='D30' & grepl('Phyto', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0)) + ylim(-3,3)

# tocontf.genus %>% filter(day =='D30' & grepl('SBP', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
#   geom_text(aes(label=Genus, y=0))


########## Now OTUs ############

#### otus for each treatment at D23 feces

tocontf.otu %>% filter(day =='D23' & grepl('Acid', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
  geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D23' & grepl('RCS', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
  geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D23' & grepl('Bglu', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
  geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D23' & grepl('Malto', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D23' & grepl('RPS', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0)) + ylim(-6,7) + ggtitle('D23 RPS')

tocontf.otu %>% filter(day =='D23' & grepl('Pos_control', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))  + ggtitle('D23 Pos_control')

tocontf.otu %>% filter(day =='D23' & grepl('ZnCu', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
  geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D23' & grepl('Phyto', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D23' & grepl('SBP', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
  geom_text(aes(label=Genus, y=0))

####### OTUs D30

# tocontf.otu %>% filter(day =='D30' & grepl('Acid', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D30' & grepl('RCS', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D30' & grepl('Bglu', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
  geom_text(aes(label=Genus, y=0))
# 
# tocontf.otu %>% filter(day =='D30' & grepl('Malto', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) +
#   geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D30' & grepl('RPS', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=tissue)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0)) + ylim(-10,30) + ggtitle('D30 RPS')

tocontf.otu %>% filter(day =='D30' & grepl('Pos_control', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))

tocontf.otu %>% filter(day =='D30' & grepl('ZnCu', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
  geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
  geom_text(aes(label=Genus, y=0))

# tocontf.otu %>% filter(day =='D30' & grepl('Phyto', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
#   geom_text(aes(label=Genus, y=0)) + ylim(-3,3)

# tocontf.otu %>% filter(day =='D30' & grepl('SBP', Treatment))%>% ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) + 
#   geom_col(color='black') + coord_flip() + #scale_fill_manual(values = fin_colset) + 
#   geom_text(aes(label=Genus, y=0))

tocontf.genus <- tocontf.genus%>% mutate(comp2=sub('(.*)_vs_Control','\\1',comp))
tocontf.otu <- tocontf.otu %>%  mutate(comp2=sub('(.*)_vs_Control','\\1',comp))

tocontf.genus$comp2 <- factor(tocontf.genus$comp2, levels = c('Acid', 'Bglu', 'Malto', 'Phyto', 'Pos_control', 'RCS', 'RPS', 'SBP', 'ZnCu'))
tocontf.otu$comp2 <- factor(tocontf.otu$comp2, levels = c('Acid', 'Bglu', 'Malto', 'Phyto', 'Pos_control', 'RCS', 'RPS', 'SBP', 'ZnCu'))


tocontf.genus %>% filter(day != 'D0') %>% group_by(comp2, day, tissue) %>% summarise(num_diff_genera=n()) %>% 
  ggplot(aes(x=comp2, y=num_diff_genera)) + geom_col(aes(fill=comp2),color='black') + facet_wrap(~day) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) + scale_fill_brewer(palette = 'Set1', drop=FALSE) + scale_x_discrete(drop=FALSE)



tocontf.otu %>% filter(day != 'D0') %>% group_by(comp2, day, tissue) %>% summarise(num_diff_otus=n()) %>% 
  ggplot(aes(x=comp2, y=num_diff_otus)) + geom_col(aes(fill=comp2),color='black') + facet_wrap(~day) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) + scale_fill_brewer(palette = 'Set1', drop=FALSE) + scale_x_discrete(drop=FALSE)


tocontf.genus %>% write_csv('Diff_abund_genera_all.csv')
tocontf.otu %>% write_csv('Diff_abund_OTUs_all.csv')

FS12a@sam_data %>% group_by(day, tissue, treatment) %>% tally()


taxa$Genus
###### BLAST RESULTS FOR ALL THESE OTUS

# tocontf %>% select(OTU) %>% unlist(use.names = FALSE) %>% write_lines('./data/FS12a_int_OTUs.txt')

#FS12a diversity

### NEEDS ATTENTION ###

FS12a_jac <- vegdist(FS12a_rare@otu_table, method = 'bray')
# FS12a_jac

attr(FS12a_jac, which = 'Labels') == FS12a_rare@sam_data$sample_ID
FS12a_dispers <- betadisper(FS12a_jac, group = FS12a_rare@sam_data$set)
FS12a_pdispers <- permutest(FS12a_dispers, pairwise = TRUE)
FS12a_pdispers$pairwise$observed
FS12a_dispersdf <- data.frame(FS12a_dispers$distances)
FS12a_dispersdf$group <- rownames(FS12a_dispersdf)
FS12a@sam_data$sample_ID == FS12a_dispersdf$group


meta$sample_ID %in% FS12a_dispersdf$group

colnames(FS12a_dispersdf)[2] <- 'sample_ID'
FS12a_meta <- FS12a_NMDS[[1]]

FS12a_meta <- merge(FS12a_meta, FS12a_dispersdf, by='sample_ID')
FS12a_meta$day <- as.numeric(gsub('D', '', FS12a_meta$day))
FS12a_meta$dayfact <- factor(FS12a_meta$day)


# FS12a_meta$treatment <- factor(FS12a_meta$treatment, levels = c('Control', 'RPS', 'Acid','ZnCu', 'RCS', 'Bglu'))
# names(diversity(FS12a_rare@otu_table)) == FS12a_meta$sample_ID

FS12a_meta$shan <- diversity(FS12a_rare@otu_table)

FS12a_meta$rich <- specnumber(FS12a_rare@otu_table)
FS12a_meta$even <- FS12a_meta$shan/log(FS12a_meta$rich)


#fecal shannon
FS12a_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=shan, group=set, fill = treatment)) +
  geom_boxplot() + facet_wrap(~day) + scale_fill_manual(values = fills)+ theme(axis.text.x = element_text(angle = 90))
# scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Fecal Shannon Diversity (alpha) over time')  + geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)#+ geom_text(aes(label=pignum))


#fecal even
FS12a_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=even, group=set, fill = treatment)) +
  geom_boxplot() + facet_wrap(~day) + scale_fill_manual(values = fills) + theme(axis.text.x = element_text(angle = 90))
# scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Fecal evenness over time')+ geom_jitter(shape = 21, stroke=1.2, size=2, width = .2) #+ geom_text(aes(label=pignum))

#fecal rich
FS12a_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=rich, group=set, fill = treatment)) +
  geom_boxplot() + facet_wrap(~day) + scale_fill_manual(values = fills) + theme(axis.text.x = element_text(angle = 90))
# scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple')) + ggtitle('Fecal richness (num OTUs) over time')+ geom_jitter(shape = 21, stroke=1.2, size=2, width = .2) #+ geom_text(aes(label=pignum))


#fecal dispersion
FS12a_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=FS12a_dispers.distances, group=set, fill = treatment)) +
  geom_boxplot() + facet_wrap(~day) + scale_fill_manual(values = fills) + theme(axis.text.x = element_text(angle = 90))
# scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + ggtitle('Fecal community dispersion over time')+ geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)#+ geom_text(aes(label=pignum))


#### FS12a weight here ####
#### should do D30 here too....

FS12a_prune <- prune_samples(!(is.na(FS12a@sam_data$nurse_gain)), FS12a)

glob_weight_assoc_23 <- blarg(phyloseq_obj = FS12a_prune, day = 'D23', tissue = 'F', covariate = 'nurse_gain')


glob_weight_assoc_0 <- blarg(phyloseq_obj = FS12a_prune, day = 'D0', tissue = 'F', covariate = 'nurse_gain')
# does it make sense to look at linear relationships between OTU abundance prior to weaning and nursery preformance?

##############
# Control Pen effect stuff


# NO DONT DO IT





Pentrols <- prune_samples(FS12@sam_data$day == 'D23' & FS12@sam_data$treatment == 'Control', x = FS12)


Pen_rare <- rarefy_even_depth(Pentrols)
Pen_rare@sam_data$pen

penmds <- NMDS_ellipse(Pen_rare@otu_table, metadata = Pen_rare@sam_data, grouping_set = 'pen')



penmds[[1]] %>% ggplot(aes(x=MDS1, y=MDS2, xend=centroidX, yend=centroidY, color=pen)) + geom_point() + geom_segment()


global_adon_D23 <- adonis(FS12a_D23@otu_table~treatment,
                          data = data.frame(FS12a_D23@sam_data), method = 'bray')


adonis(Pen_rare@otu_table~pen, data = data.frame(Pen_rare@sam_data))

pairwise.adonis(x = Pen_rare@otu_table, Pen_rare@sam_data$pen, binary = FALSE, p.adjust.m = 'fdr')

# no strong detectable pen effect in terms of total community structure
# just because penmates, doesn't suggest they will have similar fecal bacterial community structures

