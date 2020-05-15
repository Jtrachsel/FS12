getwd()

library(phyloseq)
library(tidyverse)
library(funfuns)
# library(pairwiseAdonis)

theme_set(cowplot::theme_cowplot())  


new_names <- c(NC='Control',
               CT='Pos_control',
               RSpo='RPS',
               SCF='Malto',
               SBP='SBP',
               FAM='Acid',
               PHY='Phyto',
               CuZn='ZnCu',
               RScn='RCS',
               BG='Bglu')

blarg <- function(phyloseq_obj, day, tissue, covariate, shrink_type='apeglm'){
  form <- formula(paste('~', covariate))
  # print(form)
  FS12b.glom <- phyloseq_obj %>% prune_samples(samples = phyloseq_obj@sam_data$day %in% c(day) & phyloseq_obj@sam_data$tissue == tissue)
  FS12b.glom <- prune_taxa(taxa_sums(FS12b.glom) > 1, FS12b.glom)
  
  # FS12b.glom@sam_data$log_sal
  
  FS12b.de <- phyloseq_to_deseq2(FS12b.glom, form)
  FS12b.de <- DESeq(FS12b.de, test = 'Wald', fitType = 'parametric')
  
  # these are not both possible.  Right now only lfcshrink is doing anytihng
  res <- results(FS12b.de, cooksCutoff = FALSE, name = covariate)
  res <- lfcShrink(FS12b.de, coef = covariate, type = shrink_type)
  
  # resultsNames(FS12b.de)
  
  res <- res[!is.na(res$padj),]
  res <- res[res$padj < 0.1,]
  sigtab <- res[abs(res$log2FoldChange) > .1 ,]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
  sigtab$newp <- format(round(sigtab$padj, digits = 3), scientific = TRUE)
  # sigtab$Treatment <- ifelse(sigtab$log2FoldChange >=0, treat, paste('down',treat, sep = '_'))
  sigtab$OTU <- rownames(sigtab)
  sigtab[[covariate]] <- ifelse(sigtab$log2FoldChange >0 , 'increased', 'decreased')
  # sigtab$salm <- ifelse(sigtab$log2FoldChange >0 , 'increased', 'decreased')
  sigtab <- sigtab[order(sigtab$log2FoldChange),]
  sigtab$OTU <- factor(sigtab$OTU, levels = sigtab$OTU)
  sigtab$day <- day
  sigtab$tissue <- tissue
  
  
  p <- sigtab %>% ggplot(aes_string(x='OTU', y='log2FoldChange', fill=covariate)) +
    geom_col(color='black') + coord_flip() + geom_text(aes(label=Genus, y=0))
  
  return(list(p, sigtab))
  
  
}
########

#### data_wrangling ####
# meta <- read.csv('./data/FS12_final_meta.csv', header = TRUE, stringsAsFactors = FALSE)
# shared <- read_delim('./data/FS12.shared', delim = '\t') %>% as.data.frame()
# 
# taxa <- extract_mothur_tax('./data/FS12.taxonomy')
# # taxa[1:10,1:8]
# 
# rownames(shared) <- shared$Group
# shared <- shared[,-c(1,2,3)]
# 
# # hist(rowSums(shared), breaks = 50)
# 
# 
# # length(colnames(shared)) # 10648 total OTUs detected (no removal mocks and NTCs and all kinds of garbage here)
# 
# 
# # rownames(shared) %in% meta$sample_ID
# meta <- meta[meta$sample_ID %in% rownames(shared),]
# mocks <- shared[(!rownames(shared) %in% meta$sample_ID),]
# shared <- shared[rownames(shared) %in% meta$sample_ID,]
# 
# rownames(shared) == meta$sample_ID
# 
# rownames(meta) <- meta$sample_ID
# 
# shared <- shared[rowSums(shared) > 1250,] # remove samples with less than 1250 reads
# shared <- shared[,colSums(shared > 0) > 3] # removes otus that are detected in fewer than 4 samples globally (including timecourse data)
# shared <- shared[,colSums(shared) > 5] # at least 6 observations globally
# 
# length(colnames(shared))
# hist(rowSums(shared), breaks = 50)
# 
# 
# meta <- meta[order(meta$sample_ID),]
# shared <- shared[order(rownames(shared)),]
# 
# rownames(shared) == meta$sample_ID
# 
# meta$set <- paste(meta$experiment, meta$day, meta$tissue, meta$treatment, sep = '_')
# 
# nnnn <- meta %>% group_by(set) %>% summarise(N=n())
# 
# 
# rownames(meta) == rownames(shared)
# meta <- meta[rownames(meta) %in% rownames(shared),]
# rownames(shared) %in% rownames(meta)
# shared <- shared[match(rownames(meta), rownames(shared)),]
# 
# rownames(meta) == rownames(shared)
# 
# # taxa <- taxa[taxa$OTU %in% colnames(shared),]
# 
# # Fix metadata group labels here:
# meta$treatment <- factor(meta$treatment, levels = c('Control',
#                                                     'Pos_control',
#                                                     'RPS',
#                                                     'Malto',
#                                                     'SBP',
#                                                     'Acid',
#                                                     'Phyto',
#                                                     'ZnCu',
#                                                     'RCS',
#                                                     'Bglu'))
# 
# new_names <- c(NC='Control',
#                CT='Pos_control',
#                RSpo='RPS',
#                SCF='Malto',
#                SBP='SBP',
#                FAM='Acid',
#                PHY='Phyto',
#                CuZn='ZnCu',
#                RScn='RCS',
#                BG='Bglu')
# 
# meta <- meta %>% rownames_to_column() %>%
#   mutate(treatment2 = fct_recode(treatment, !!!new_names)) %>%
#   column_to_rownames()
# 
# meta$treatment2
# #
# 
# 
# 
# p_meta <- sample_data(meta)
# p_taxa <- import_mothur(mothur_constaxonomy_file = './data/FS12.taxonomy')
# 
# 
# colnames(p_taxa) <- colnames(taxa)[-c(1,2,3)]
# 
# # meta$treatment
# FS12 <- phyloseq(p_meta, p_taxa, otu_table(shared, taxa_are_rows = FALSE))  # builds the phyloseq obj
# 
# #
# # saveRDS(FS12, "FS12_phyloseq.rds")
# # readRDS('FS12_phyloseq.rds')# Close and re-open R
# 
# 
# ### mock stuff ###
# 
# # rownames(mocks)
# # mocks[1:5,1:5]
# # mocks <- mocks[rowSums(mocks) >1000,]
# # rownames(mocks)
# # hist(rowSums(mocks), breaks = 50)
# #
# # mocks <- mocks[grep('NTC', rownames(mocks)),]
# #
# # mocks <- mocks[,colSums(mocks)>2]
# # colSums(mocks)
# # rownames(mocks)
# # mock_tax <- taxa[taxa$OTU %in% colnames(mocks),]
# #
# #
# 
#
# FS12a <- FS12 %>% subset_samples(experiment == 'X12a')
# 
# 
# 
# 
# FS12a@sam_data$treatment
# 
# #
# # FS12a@sam_data$treatment <- factor(FS12a@sam_data$treatment, levels = c('Control',
# #                                                                         'RPS',
# #                                                                         'Acid',
# #                                                                         'ZnCu',
# #                                                                         'RCS',
# #                                                                         'Bglu',
# #                                                                         'Pos_control',
# #                                                                         'Phyto',
# #                                                                         'Malto',
# #                                                                         'SBP'))
# 
# hist(sample_sums(FS12a), breaks = 50)
# 
# 
# # FS12a@sam_data %>% filter(day == 'D30' & treatment == 'Phyto')
# 
# # FS12a <- FS12a %>% subset_samples(day %in% c('D0', 'D23', 'D30'))
# FS12a <- FS12a %>% subset_samples(day %in% c('D0', 'D23'))
# 
# 
# rems <- !(FS12a@sam_data$sample_ID %in% c('X12N64WD0dup', 'X12N77WD0dup', 'X12N70XD30'))
# FS12a <- prune_samples(x = FS12a, rems)
# 
# FS12a@sam_data %>% group_by(day, treatment) %>% tally()
# 
# 
# 
# 
# # FS12a@sam_data %>% filter(day == 'D23' & treatment == 'Control')
# 
# 
# keeps_at_d23 <- FS12a@sam_data %>% filter(day == 'D0' & treatment == 'Control') %>% select(pignum) %>% unique() %>% unlist()
# REMOVE_THESE <- FS12a@sam_data %>% filter(day == 'D23' & treatment == 'Control') %>% filter(!(pignum %in% keeps_at_d23)) %>% select(sample_ID) %>% unlist()
# ANDTHESE <- c('X12N16WD23', 'X12N18WD23', 'X12N32WD23', 'X12N39WD23', 'X12N45WD23', 'X12N68WD23')
# REMOVE_THESE <- c(REMOVE_THESE, ANDTHESE)
# 
# FS12a <- FS12a %>% subset_samples(!(sample_ID %in% REMOVE_THESE))
# FS12a@sam_data %>% group_by(day, treatment) %>% tally()
# 
# # MAKE FINAL FS12a OBJECT HERE AND SAVE SO CAN READ IN LATER
# 
# 
# 
# saveRDS(FS12a, "FS12a_phyloseq.rds")
#### read in phyloseq ####

FS12a <- readRDS('FS12a_phyloseq.rds')# Close and re-open R

# sample_data(FS12a)


# PERMANOVAS VS CONTROL

# FS12a@sam_data %>% group_by(day, tissue) %>% mutate(minseqs=min(num_seqs)) %>% select(day, tissue, minseqs) %>% unique()


FS12a@sam_data$set <- paste(FS12a@sam_data$day, FS12a@sam_data$tissue, FS12a@sam_data$treatment, sep='_')

FS12a_rare <- phyloseq::rarefy_even_depth(FS12a,sample.size =  min(phyloseq::sample_sums(FS12a)),
                                          rngseed = 10)


FS12a_NMDS <- NMDS_ellipse(FS12a_rare@sam_data, OTU_table = FS12a_rare@otu_table, grouping_set = 'set', distance_method = 'bray', MDS_trymax = 50)

labbies <- FS12a_NMDS[[1]] %>% select(tissue, treatment, day, set, centroidX, centroidY, treatment2) %>% unique()
# library(ggrepel)

#



library(RColorBrewer)
colourCount = length(unique(FS12a@sam_data$treatment))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

fills <- c('black', getPalette(9))

unique(FS12a_NMDS[[1]]$treatment)


#### THIS NEEDS TO BE CHANGED ####
FS12a_NMDS[[1]]$treatment <- factor(FS12a_NMDS[[1]]$treatment, levels = c('Control', 'Acid', 'Bglu', 'Malto', 'Phyto', 'Pos_control', 'RCS', 'RPS', 'SBP', 'ZnCu'))
# getPal

labbies$treatment <- factor(labbies$treatment, levels =c('Control', 'Acid', 'Bglu', 'Malto', 'Phyto', 'Pos_control', 'RCS', 'RPS', 'SBP', 'ZnCu'))
### THIS IS A GOOD ONE. SHOW WITH GLOBAL ADONIS TEST
### FIGURE 1 ###
FS12a_NMDS[[1]] %>% ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=day), size=2)+
  geom_segment(aes(xend=centroidX, yend=centroidY, color=day)) +
  geom_point(data=labbies,aes(x=centroidX, y=centroidY, fill=treatment2), size=3, inherit.aes = FALSE, shape=21)+# scale_shape_manual(values = c(21, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + scale_fill_manual(values = fills, name='Treatment')
# geom_text_repel(data=labbies, aes(label=treatment, x=centroidX, y=centroidY))#+ geom_path(data = FS12a_NMDS[[2]], aes(x=NMDS1, y=NMDS2, group=group))



### FIGURE 2 ###
set.seed(1)
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

same_day_tissue$treatment <- factor(same_day_tissue$treatment,
                                      levels = c('Pos_control',
                                                 'RPS',
                                                 'Malto',
                                                 'SBP',
                                                 'Acid',
                                                 'Phyto',
                                                 'ZnCu',
                                                 'RCS',
                                                 'Bglu'))

  
  
same_day_tissue <- same_day_tissue %>% mutate(treatment2=fct_recode(treatment, !!!new_names))
# fig2 #
p <- same_day_tissue %>% filter(set == 'feces_D23') %>% 
      ggplot(aes(x=treatment2, y=F.Model, fill=treatment2)) + 
      geom_col(color='black') + #facet_wrap(~set) +
      geom_text(aes(label=p.plot), nudge_y = .2) + 
      # ggtitle('Community differences compared to controls', subtitle = 'larger F = greater difference.  pvalues shown when <0.05') + 
      scale_fill_brewer(palette = 'Set1', name='Treatment') + xlab('Treatment')
p  

# global adonis

# end figure 2 #


FS12a_D0 <- FS12a_rare %>% subset_samples(day == 'D0')
FS12a_D23 <- FS12a_rare %>% subset_samples(day == 'D23')


global_adon <- adonis(FS12a_rare@otu_table~day + treatment,
                      data = data.frame(FS12a_rare@sam_data), method = 'bray')

global_adon_D0 <- adonis(FS12a_D0@otu_table~treatment,
                         data = data.frame(FS12a_D0@sam_data), method = 'bray')

global_adon_D23 <- adonis(FS12a_D23@otu_table~treatment,
                          data = data.frame(FS12a_D23@sam_data), method = 'bray')




global_adon_D23
#day F.model=54.605 pval=0.001
#treat day0 F.model=1.1205, pval=0.109
#treat day 24 F.model=1.4638 pval=0.001


# write <- global_adon$aov.tab
# write$variable <- rownames(write)
# 
# write <- write %>% select(variable, everything())
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
FS12a@sam_data$treatment <- fct_recode(FS12a@sam_data$treatment, !!!new_names)
# 
# unique(FS12a@sam_data$treatment)
# 
# FS12a.glom <- FS12a %>% tax_glom(taxrank = 'Genus')

# this can stay right?
# something weird is happening here....  treatment factor levels being replaced by integers
# it happened after I loaded te Rtsne package.... went away when I restarted R

# begin figure 3 #

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

tocont.otu <- list(DESeq_difabund(phyloseq = FS12a, day = 'D0', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'),
                   DESeq_difabund(phyloseq = FS12a, day = 'D23', tissue = 'F', scientific = TRUE, shrink_type = 'normal',alpha = 0.05, cooks_cut = FALSE, pAdjustMethod = 'BH'))

tmpres$OTU <- bind_rows(tocont.otu)

levs <- c('RPS','Acid', 'ZnCu', 'RCS', 'Bglu', 'Pos_control', 'Phyto', 
          'Malto', 'SBP', 'down_RPS','down_Acid', 'down_ZnCu', 'down_RCS', 
          'down_Bglu', 'down_Pos_control', 'down_Phyto', 'down_Malto', 'down_SBP')

pal <- getPalette(9)

pal_light <- lapply(X = pal, FUN = lighten) %>% unlist()

my_pal <- c(pal, pal_light)
names(my_pal) <- levs


###NEW NAMES ###

NEW_NAMES2 <- c(CT='Pos_control',
                RSpo='RPS',
                SCF='Malto',
                SBP='SBP',
                FAM='Acid',
                PHY='Phyto',
                CuZn='ZnCu',
                RScn='RCS',
                BG='Bglu',
                down_CT='down_Pos_control',
                down_RSpo='down_RPS',
                down_SCF='down_Malto',
                down_SBP='down_SBP',
                down_FAM='down_Acid',
                down_PHY='down_Phyto',
                down_CuZn='down_ZnCu',
                down_RScn='down_RCS',
                down_BG='down_Bglu')


# 
### I THINK THIS ONE IS GOOD!
nrow(tmpres$OTU)
# wrong
# tmpres$OTU %>% mutate(comp2=fct_recode(comp, !!!NEW_NAMES2))

tmpres$OTU$comp <- factor(tmpres$OTU$comp ,
                          levels = c("CT_vs_NC",
                                     "RSpo_vs_NC",
                                     "SCF_vs_NC",
                                     "SBP_vs_NC",
                                     "FAM_vs_NC",
                                     "PHY_vs_NC",
                                     "CuZn_vs_NC",
                                     "RScn_vs_NC",
                                     "BG_vs_NC"))


tmpres$OTU %>% group_by(comp) %>%
  tally() %>% 
  ggplot(aes(x=comp, y=n, fill=comp)) +
  geom_col(color='black') + 
  scale_fill_brewer(palette = 'Set1') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle('Number of significantly differentially abundant OTUs compared to the Control group')

theseones <- c('Malto_vs_Control', 'Pos_control_vs_Control', 'RPS_vs_Control')

tmpres$OTU %>% filter(day !='D0') %>%#filter(comp %in% theseones) %>%
  ggplot(aes(x=Family, y=log2FoldChange, fill=comp)) +
  geom_hline(yintercept = 0, color='grey50', size=2) +
  geom_point(color='black', shape=21, size=3) + coord_flip() + 
  scale_fill_brewer(palette = 'Set1', name='Comparison')+
  theme(axis.text.y = element_text(size=10),panel.grid.major.x = element_line(color='grey')) +
  ggtitle('')




unique(tmpres$OTU$Family)
gram_negs <- c('Enterobacteriaceae', 'Desulfovibrionaceae', 'Proteobacteria_unclassified', 'Campylobacteraceae','Deltaproteobacteria_unclassified', 'Gammaproteobacteria_unclassified')


# CHANGE LEGEND TITLE #
tmpres$OTU %>% filter(day !='D0') %>%
  filter(Family %in% gram_negs) %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=comp)) +
  geom_hline(yintercept = 0, color='grey50', size=2) +
  geom_point(color='black', shape=21, size=4) + coord_flip() + 
  # scale_fill_discrete(drop=FALSE)+
  scale_fill_brewer(palette = 'Set1',drop=FALSE, name='Comparison')+
  theme(panel.grid.major.x = element_line(color='grey'), 
        legend.title = element_text()) +
  ylim(-4.5,4.5)




# tmpres$Genus %>% filter(day !='D0') %>%  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) +
#   geom_text(aes(label=Genus, y=0))+ facet_wrap(~comp, scales = 'free')
# 
# 
# tmpres$Family %>% filter(day !='D0') %>%  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) +
#   geom_text(aes(label=Family, y=0))
# 
# 
# tmpres$Order %>% filter(day !='D0') %>%  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
#   geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) +
#   geom_text(aes(label=Order, y=0))
# 
# FS12a@sam_data %>% group_by(day, tissue, treatment) %>% tally()
# 
# 
# taxa$Genus
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


# meta$sample_ID %in% FS12a_dispersdf$group

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

# NEED TESTS
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


FS12a_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=nurse_gain, group=set, fill = treatment)) +
  geom_boxplot() + facet_wrap(~day) + scale_fill_manual(values = fills) + theme(axis.text.x = element_text(angle = 90))

FS12a_meta %>% filter(tissue == 'F') %>% ggplot(aes(x=treatment, y=weight_10_9, group=set, fill = treatment)) +
  geom_boxplot() + facet_wrap(~day) + scale_fill_manual(values = fills) + theme(axis.text.x = element_text(angle = 90))



# scale_fill_manual(values=c('#33CC33', '#3399FF', 'orange', 'red', 'grey', 'purple'))  + ggtitle('Fecal community dispersion over time')+ geom_jitter(shape = 21, stroke=1.2, size=2, width = .2)#+ geom_text(aes(label=pignum))


######### tests  ##########


FS12a_meta$day
tests <- FS12a_meta %>% filter(day ==23) 

res.aov <- aov(shan ~ treatment, data = tests)
summary.aov(res.aov)

res.aov <- aov(rich ~ treatment, data = tests)
summary.aov(res.aov)




plant.lm <- lm(rich ~ treatment, data = tests)
plant.av <- aov(plant.lm)
summary(plant.av)

tukey.test <- TukeyHSD(plant.av)
tukey.test
tukey.test$treatment[tukey.test$treatment[,4] < 0.05,]
#### FS12a weight here ####
#### should do D30 here too....

FS12a_prune <- prune_samples(!(is.na(FS12a@sam_data$nurse_gain)), FS12a)
FS12a_prune@sam_data %>% group_by(day, treatment) %>% tally()

### NEED BLARG FUNCTION ###

glob_weight_assoc_23 <- blarg(phyloseq_obj = FS12a_prune, day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
nrow(glob_weight_assoc_23[[2]])

glob_weight <- glob_weight_assoc_23[[2]]

p <- glob_weight_assoc_23[[1]] #+ xlab("Strength of ")

ggsave(filename = './side_proj/Kerr_figs/NURSE_GAIN_OTUs.jpeg', plot = p, 
        device = 'jpeg', width = 170, height = 185, units = 'mm')

##### ASIDE ####

# THESE ARE OTUS that have a linear relationship with nursery preformance and are also
# enriched in treatments,
tmpres$OTU %>% filter(day !='D0') %>% filter(OTU %in% glob_weight$OTU) %>% 
  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_col(color='black') + coord_flip() + scale_fill_manual(values = my_pal) +
  geom_text(aes(label=Genus, y=0)) + facet_wrap(~comp, scales = 'free')



glob_weight$log2FoldChange > 0

tmpres$OTU %>% filter(day !='D0') %>% filter(OTU %in% glob_weight$OTU) %>% 
  ggplot(aes(x=OTU, y=log2FoldChange, fill=Treatment)) +
  geom_point(color='black', shape=21, size=3) + coord_flip() + scale_fill_manual(values = my_pal) +
  geom_text(aes(label=Genus, y=0))# + facet_wrap(~comp, scales = 'free')


tst <- tmpres$OTU
tryme <- glob_weight %>% transmute(OTU=OTU, 
                          weight_L2FC=log2FoldChange) %>% right_join(tst) %>% na.omit()



tryme %>% ggplot(aes(x=OTU, y=weight_L2FC, color=Treatment)) + geom_point(aes(size=abs(log2FoldChange))) + coord_flip()

sum(tmpres$OTU$OTU %in% glob_weight$OTU)

tmpres$OTU %>% filter(OTU %in% glob_weight$OTU)

# glob_weight_assoc_0 <- blarg(phyloseq_obj = FS12a_prune, day = 'D0', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
# 
# unique(meta$treatment)

  
###################
FS12a_prune@sam_data$treatment

Cont_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'NC'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
RPS_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'RSpo'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
RCS_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'RScn'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
Acid_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'FAM'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
Bglu_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'BG'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
ZnCu_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'CuZn'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
Phyto_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'PHY'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
Pos_control_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'CT'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
SBP_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'SBP'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')
Malto_weight_23 <- blarg(phyloseq_obj = subset_samples(FS12a_prune, treatment == 'SCF'), day = 'D23', tissue = 'F', covariate = 'nurse_gain', shrink_type = 'normal')


Cont_weight_23[[2]]$treat <- 'NC'
RPS_weight_23[[2]]$treat <- 'RSpo'
RCS_weight_23[[2]]$treat <- 'RScn'
Acid_weight_23[[2]]$treat <- 'FAM'
Bglu_weight_23[[2]]$treat <- 'BG'
ZnCu_weight_23[[2]]$treat <- 'CuZn'
Phyto_weight_23[[2]]$treat <- 'PHY'
Pos_control_weight_23[[2]]$treat <- 'CT'
SBP_weight_23[[2]]$treat <- 'SBP'
Malto_weight_23[[2]]$treat <- 'SCF'


all_treat_weight <- bind_rows(list(Cont_weight_23[[2]],
                                   RPS_weight_23[[2]],
                                   RCS_weight_23[[2]],
                                   Acid_weight_23[[2]],
                                   ZnCu_weight_23[[2]],
                                   Phyto_weight_23[[2]],
                                   Pos_control_weight_23[[2]],
                                   SBP_weight_23[[2]], 
                                   Malto_weight_23[[2]]))


#####

### I THINK THIS IS WHERE THE TREATMENTS NEED TO BE FIXED #
sig_OTUs <- tmpres$OTU

sig_OTUs$treat2 <- sub('down_treatment_','',sig_OTUs$Treatment)
sig_OTUs$treat2 <- sub('treatment_','',sig_OTUs$treat2)
sig_OTUs$treat2 <- sub('_vs_NC','',sig_OTUs$treat2)


nested_sigs <- sig_OTUs %>% group_by(treat2) %>% nest()

sigs_n_weight <- all_treat_weight %>% mutate(treat2=treat) %>%
  group_by(treat2) %>% nest(.key = 'weight_data') %>% right_join(nested_sigs)


sigs_n_weight[3,]$weight_data
sigs_n_weight[3,]$data


# sigs_n_weight[3,]$weight_data[[1]]$OTU %in% sigs_n_weight[3,]$data[[1]]$OTU


testtt <- function(cont_cov_dat, diffabund_dat){
  these <- cont_cov_dat$OTU %in% diffabund_dat$OTU
  tmp1 <- cont_cov_dat %>% filter(these) %>%
    transmute(OTU=OTU, 
              cont_cov_L2FC = log2FoldChange) %>% 
    right_join(diffabund_dat) %>% na.omit()
  
  
}



sigs_n_weight <- sigs_n_weight %>% filter(!map_lgl(.x = weight_data, .f = is.null))

sigs_n_weight %>%
  mutate(good_OTUs=map2(.x = weight_data, .y = data, .f=testtt)) %>% 
  select(treat2, good_OTUs) %>% unnest() %>%
  mutate(treatment=factor(treat2, levels = names(new_names)[-1])) %>% 
  ggplot(aes(x=Genus, y=cont_cov_L2FC, fill=treatment))+ 
  geom_hline(yintercept = 0) +
  geom_point(shape=21, size=4) + scale_fill_brewer(palette = 'Set1', drop=FALSE, name='Treatment') +
  coord_flip() + ylim(-1.75, 1.75) + 
  geom_text(aes(y=1.25, label=OTU))+
  theme(panel.grid.major.x = element_line(color='grey'))+
  ylab('Strength of relationship with nursery preformance') #+
  ggtitle('OTUs significantly associated with a treatment\n and nursery preformance in that same treatment')



# 
# 
# all_treat_weight <- bind_rows(list(Cont_weight_23[[2]],
#           RPS_weight_23[[2]],
#           RCS_weight_23[[2]],
#           Acid_weight_23[[2]],
#           ZnCu_weight_23[[2]],
#           Phyto_weight_23[[2]],
#           Pos_control_weight_23[[2]],
#           SBP_weight_23[[2]], 
#           Malto_weight_23[[2]]))



### using these two dataframes, I want to identify otus that are both enriched in a diet relative to the negative controls
### AND have a linear relationship with nursery preformance

sig_otus_by_treat<- tmpres$OTU
sig_otus_by_treat <- sig_otus_by_treat %>% filter(day != 'D0')


ints <- all_treat_weight %>% group_by(OTU) %>% tally() %>% filter(n>1)

pos_nurse_gain <- all_treat_weight %>% filter(log2FoldChange > 0)


### One option...

all_treat_weight %>% filter(log2FoldChange > 0) %>% 
  ggplot(aes(x=OTU, y=log2FoldChange, fill=treat)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label=Genus, y=0), hjust = 'left')


###

all_treat_weight %>% filter(OTU %in%  sig_otus_by_treat$OTU) %>% filter(log2FoldChange > 0) %>% 
  ggplot(aes(x=OTU, y=log2FoldChange, fill=treat)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label=Genus, y=0), hjust = 'left')


all_treat_weight %>% filter(OTU %in%  sig_otus_by_treat$OTU) %>% filter(log2FoldChange < 0) %>% 
  ggplot(aes(x=OTU, y=log2FoldChange, fill=treat)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label=Genus, y=0), hjust = 'right')



all_treat_weight %>% group_by(OTU, nurse_gain) %>% tally() %>% filter(n>1)

all_treat_weight %>% filter(log2FoldChange >0.25) %>% filter(treat != 'Control') %>% 
  filter(OTU %in% sig_otus_by_treat$OTU) %>% 
  ggplot(aes(x=OTU, y=log2FoldChange, fill=treat)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label=Genus, y=0), hjust = 'left')



all_treat_weight %>% filter(log2FoldChange < -0.25) %>% filter(treat != 'Control') %>% 
  filter(OTU %in% sig_otus_by_treat$OTU) %>% 
  ggplot(aes(x=OTU, y=log2FoldChange, fill=treat)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label=Genus, y=0), hjust = 'right')







#### new zone #### merge?


tests <- all_treat_weight %>% transmute(OTU = OTU,
                               weight_L2FC=log2FoldChange,
                               weight_pval=padj, 
                               weight_treat=treat) %>% 
  right_join(sig_otus_by_treat) %>% na.omit()

###

all_treat_weight[all_treat_weight$OTU %in% sig_otus_by_treat$OTU,] %>% 
  ggplot(aes(x=OTU, y=log2FoldChange, fill=treat)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label=Genus, y=0)) + facet_wrap(~treat, scales = 'free')


glob_weight$OTU %in% all_treat_weight$OTU

# does it make sense to look at linear relationships between OTU abundance prior to weaning and nursery preformance?
 

















##############
# Control Pen effect stuff


# NO DONT DO IT





Pentrols <- prune_samples(FS12@sam_data$day == 'D23' & FS12@sam_data$treatment == 'Control' & FS12@sam_data$pen != 25, x = FS12)


Pen_rare <- rarefy_even_depth(Pentrols)
unique(Pen_rare@sam_data$pen)

Pen_rare@sam_data %>% group_by(pen) %>% tally()


penmds <- NMDS_ellipse(Pen_rare@otu_table, metadata = Pen_rare@sam_data, grouping_set = 'pen')



penmds[[1]] %>% ggplot(aes(x=MDS1, y=MDS2, xend=centroidX, yend=centroidY, color=pen)) + geom_point(alpha=.3) + geom_segment(alpha=.3) + 
  geom_path(data = penmds[[2]], aes(x=NMDS1, y=NMDS2, color=group), inherit.aes = FALSE)


# global_adon_D23 <- adonis(FS12a_D23@otu_table~treatment,
#                           data = data.frame(FS12a_D23@sam_data), method = 'bray')
# 

adonis(Pen_rare@otu_table~pen, data = data.frame(Pen_rare@sam_data))

pairwise.adonis(x = Pen_rare@otu_table, Pen_rare@sam_data$pen, p.adjust.m = 'fdr')

# no strong detectable pen effect in terms of total community structure
# just because penmates, doesn't suggest they will have similar fecal bacterial community structures

