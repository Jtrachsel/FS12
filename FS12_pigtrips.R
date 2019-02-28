# library(dada2)
# 
# 
# list.files('../FS12/16S_raw/')
# path <- '../FS12/16S_raw/'
# list.files()
# fnFs <- sort(list.files(pattern="_R1_001.fastq", full.names = TRUE))
# fnRs <- sort(list.files(pattern="_R2_001.fastq", full.names = TRUE))
# # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# 
# 
# plotQualityProfile(fnFs[21:22])
# plotQualityProfile(fnRs[21:22])
# 
# 
# filt_path <- file.path("filtered") # Place filtered files in filtered/ subdirectory
# filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
# 
# 
# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
#                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                      compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
# head(out)
# 
# 
# errF <- learnErrors(filtFs, multithread=TRUE)
# errR <- learnErrors(filtRs, multithread=TRUE)
# 
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ = TRUE)
# 
# derepFs <- derepFastq(filtFs, verbose=TRUE)
# derepRs <- derepFastq(filtRs, verbose=TRUE)
# # Name the derep-class objects by the sample names
# names(derepFs) <- sample.names
# names(derepRs) <- sample.names
# 
# dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
# dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# 
# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# seqtab <- makeSequenceTable(mergers)
# table(nchar(getSequences(seqtab)))
# seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# sum(seqtab.nochim)/sum(seqtab)
# 
# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
# colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
# rownames(track) <- sample.names
# head(track)
# 
# write.csv(track, file = 'dada_track.csv')
# library(dada2)
# taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# taxas <- assignSpecies(seqtab.nochim, "silva_species_assignment_v132.fa.gz")
# 
# rm(derepFs)
# rm(derepRs)
# rm(dadaFs)
# rm(dadaRs)
# 
# write.csv(seqtab.nochim, file='FS12_dadaseqtab.csv')
# write.csv(taxa, file = 'FS12_dada_taxa.csv')
# write.csv(taxas, file = 'FS12_dada_species.csv')

######## here from load workspace #########

head(seqtab.nochim)


colnames(seqtab.nochim)

ASVs <- paste('ASV', seq(length(colnames(seqtab.nochim))), sep = '')
ASV_seqs <- data.frame('ASV' = ASVs, 'seq' = colnames(seqtab.nochim))


taxa_df <- as.data.frame(taxa)

taxa_df$seq <- rownames(taxa)

taxa_df <- merge(ASV_seqs, taxa_df, by = 'seq' )



ASV_table <- seqtab.nochim
colnames(ASV_table) <- ASVs

head(ASV_table)
sum(rowSums(ASV_table))
ASV_table_trim <- ASV_table[,colSums(ASV_table) > 5] # at least 5 observations globally
ASV_table_trim2 <- ASV_table_trim[rowSums(ASV_table_trim) > 1000,] # remove samples with less than 1000 reads
ASV_table_trim3 <- ASV_table_trim2[,colSums(ASV_table_trim2 > 0) > 3] # removes otus that are detected in fewer than 4 samples globally (including timecourse data)

length(colnames(ASV_table_trim3)) # 2159 total ASVs detected

gsub('([X12ab]+)([NP]+)([0-9]+)([dWDi]+)([0-9]+)([A-Z]?)', '\\1', rownames(ASV_table_trim3))

meta <- data.frame(sample_ID = rownames(ASV_table), 
           experiment = gsub('([X12ab]+)([NP]+)([0-9]+)([dWDi]+)([0-9]+)([A-Z]?)', '\\1', rownames(ASV_table)), 
           pig_pen = gsub('([X12ab]+)([NP]+[0-9]+)([dWDi]+)([0-9]+)([A-Z]?)', '\\2', rownames(ASV_table)), 
           day = gsub('([X12ab]+)([NP]+[0-9]+)([dWDi]+[0-9]+)([A-Z]?)', '\\3', rownames(ASV_table)), 
           tissue = gsub('([X12ab]+)([NP]+[0-9]+)([dWDi]+[0-9]+)([A-Z]?)', '\\4', rownames(ASV_table)))


#write.csv(meta, 'FS12meta.csv')
#write.csv(ASV_table, 'FS12_ASV.csv')



mocks <- ASV_table[grep('[Mm]ock', rownames(ASV_table)),]
NTCs <-  ASV_table[grepl('NTC', rownames(ASV_table)) | grepl('empty', rownames(ASV_table)),]

#write.csv(ASV_table, 'ASV_table.csv')
#write.csv(taxa_df, 'dada_taxa.csv')

meta %>% filter(experiment =='X12b') %>% ggplot(aes(x=num_seqs)) + geom_histogram(color='black', fill='green4', bins = 50) + ggtitle('Frequency of read depth', subtitle = 'MiSeq V2 kit, 30 cycles, 4 plates, 378 samples, 10 Million reads passing all QC, median depth = 27617.5 reads')
hist(meta[meta$experiment =='X12b',]$num_seqs, breaks = 50)
length(meta[meta$experiment =='X12b',]$sample_ID)

meta$num_seqs <- rowSums(ASV_table)
jenn <- meta[grep('iW', meta$sample_ID),]
library(tidyverse)

jenn <- jenn %>% select(sample_ID, num_seqs)

jenn_counts <- ASV_table[rownames(ASV_table) %in% jenn$sample_ID,]

jenn_counts <- jenn_counts[,colSums(jenn_counts) > 1] # at least 10 observations globally
jenn$num_seqs <- rowSums(jenn_counts)

jenn_taxa <- taxa_df[taxa_df$ASV %in% colnames(jenn_counts),]

rownames(jenn_counts)

#write.csv(jenn, 'jenn.csv')
#write.csv(jenn_counts, 'jenn_counts.csv')
#write.csv(jenn_taxa, 'jenn_taxa.csv')


jenn_counts[,grep('ASV246', colnames(jenn_counts))]
jenn_counts[,grep('ASV4173', colnames(jenn_counts))]
jenn_counts[,grep('ASV3337', colnames(jenn_counts))]
jenn_counts[,grep('ASV2455', colnames(jenn_counts))]
jenn_counts[,grep('ASV8432', colnames(jenn_counts))]

print(c(jenn$num_seqs, jenn$sample_ID))
#ASV_table_trim2 <- ASV_table_trim[rowSums(ASV_table_trim) > 1000,] # remove samples with less than 1000 reads
#ASV_table_trim3 <- ASV_table_trim2[,colSums(ASV_table_trim2 > 0) > 4] # removes otus that are detected in fewer than 5 samples globally (including timecourse data)





rowSums(jenn_counts)





mocks_trim <- mocks[,colSums(mocks) > 0]
rowSums(mocks_trim)
library(vegan)
mocks_trim <- mocks_trim[,colSums(mocks_trim > 0) > 2]
mocks_rare <- rrarefy(mocks_trim, min(rowSums(mocks_trim)))

mocks_rare <- mocks_rare[,colSums(mocks_rare) > 0]
diversity(mocks_rare)
mock_taxa <- taxa_df[taxa_df$ASV %in% colnames(mocks_rare),]

unique(mock_taxa$Genus)

ASV_table_trim2 <- ASV_table_trim[,colSums(ASV_table_trim > 0) > 4]
length(colnames(ASV_table_trim2))

ASV_table_trim3 <- ASV_table_trim2[rowSums(ASV_table_trim2) > 1000,]


# foolin around with total disruption over experiment



rownames(meta) <- meta$sample_ID
meta_trim <- meta[meta$sample_ID %in% rownames(ASV_table_trim3),]

meta_test <- meta_trim[meta_trim$experiment == 'X12b' & meta_trim$tissue =='F',]

ASVs_test <- ASV_table_trim3[rownames(ASV_table_trim3) %in% meta_test$sample_ID,]
rownames(ASVs_test) == meta_test$sample_ID

ASVs_test <- rrarefy(ASVs_test, sample = min(rowSums(ASVs_test)))
ASVs_test <- ASVs_test[,colSums(ASVs_test) > 0]
## Gower works well.... unrarrefied
## jacard binary FALSE works well...



bray.dist <- vegdist(ASVs_test,method = 'bray', binary = FALSE)

meta_test$sample_ID == rownames(ASVs_test)

meta_test$shan2 <- diversity(ASVs_test, base = 2)
meta_test$shan <- diversity(ASVs_test)

#ggplot(meta_test, aes(x=treatment, y=shan, group=set)) + geom_boxplot()

min(rowSums(ASVs_test))
hist(rowSums(ASVs_test))
sort(rowSums(ASVs_test))

dist.data <- as.data.frame(as.matrix(bray.dist))
dist.data$from <- rownames(dist.data)

dist.gather <- gather(data = dist.data, key = 'to', value = 'distance', -from)

#


dist.gather$fromPig <- gsub('([X12ab]+)([NP]+)([0-9]+)([dWDi]+)([0-9]+)([A-Z]?)', '\\3', dist.gather$from)

#

dist.gather$FT <- paste(dist.gather$from, dist.gather$to, sep = ' ')



#dist.gather$TF <- paste(dist.gather$to, dist.gather$from, sep = ' ')

######
# all pig pairwise #

total_ground_covered <- dist.gather[grep('X12bP([0-9]+)D[0-9]+F X12bP\\1D[0-9]+F', dist.gather$FT),] %>% group_by(fromPig) %>% summarise(allpw=sum(distance),
                                                                                                                                         num=n())


total_ground_covered$treatment <- ifelse(total_ground_covered$fromPig %in% rooms$X6, 'control',
                              ifelse(total_ground_covered$fromPig %in% rooms$X7, 'RPS', 
                                     ifelse(total_ground_covered$fromPig %in% rooms$X8, 'Acid', 
                                            ifelse(total_ground_covered$fromPig %in% rooms$X9, 'Zn+Cu',
                                                   ifelse(total_ground_covered$fromPig %in% rooms$X10, 'RCS',
                                                          ifelse(total_ground_covered$fromPig %in% rooms$X11, 'Bglu', 'asdfsa'))))))



total_ground_covered <- total_ground_covered %>% filter(num == 25)

boxplot(total_ground_covered$allpw~total_ground_covered$treatment)
total_ground_covered <- merge(total_ground_covered, sum_sal, by = 'fromPig')

sum_sal
#cor.test(total_ground_covered$allpw, total_ground_covered$AULC)


total_ground_covered %>% group_by(treatment) %>% summarise(AULCvTRIP_P=cor.test(AULC, allpw)$p.value,
                                                AULCvTRIP_T=cor.test(AULC, allpw)$statistic)



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

#rowSums(pig_trips)
#pig_trips <- na.omit(pig_trips)
pig_trips$trip <- rowSums(pig_trips[,c(2:5)])
hist(pig_trips$trip, breaks = 10)



rooms <- read.csv('../FS12/Rooms.csv')

# add treatment data.  This probably isn't the best way to do this...

colnames(sum_sal)[1] <- 'fromPig'

library(funfuns)

#NMDS_ellipse(metadata = meta_test, OTU_table = ASVs_test, grouping_set = 'pig_pen')
###############################

pig_trips$treatment <- ifelse(pig_trips$fromPig %in% rooms$X6, 'control',
                             ifelse(pig_trips$fromPig %in% rooms$X7, 'RPS', 
                                    ifelse(pig_trips$fromPig %in% rooms$X8, 'Acid', 
                                           ifelse(pig_trips$fromPig %in% rooms$X9, 'Zn+Cu',
                                                  ifelse(pig_trips$fromPig %in% rooms$X10, 'RCS',
                                                         ifelse(pig_trips$fromPig %in% rooms$X11, 'Bglu', 'asdfsa'))))))




boxplot(pig_trips$trip~pig_trips$treatment)
pairwise.wilcox.test(x=pig_trips$trip, g=pig_trips$treatment, p.adjust.method = 'none')

colnames(sum_sal)[1] <- 'fromPig'
pig_trips <- merge(pig_trips, sum_sal, by = 'fromPig')
pig_trips %>% filter(treatment == "RPS") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

pig_trips %>% filter(treatment == "control") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')



pig_trips %>% group_by(treatment) %>% summarise(AULCvTRIP_P=cor.test(AULC, trip)$p.value,
                                                AULCvTRIP_T=cor.test(AULC, trip)$statistic,
                                                AULCv02_P=cor.test(AULC, D0_2)$p.value,
                                                AULCv02_T=cor.test(AULC, D0_2)$statistic,
                                                AULCv27_P=cor.test(AULC, D2_7)$p.value,
                                                AULCv27_T=cor.test(AULC, D2_7)$statistic,
                                                AULCv714_P=cor.test(AULC, D7_14)$p.value,
                                                AULCv714_T=cor.test(AULC, D7_14)$statistic,
                                                AULCv1421_P=cor.test(AULC, D14_21)$p.value,
                                                AULCv1421_T=cor.test(AULC, D14_21)$statistic)



pig_trips %>% filter(treatment == "control") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
pig_trips %>% filter(treatment == "RPS") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
pig_trips %>% filter(treatment == "Bglu") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
pig_trips %>% filter(treatment == "Zn+Cu") %>% ggplot(aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

pig_trips %>% 
  ggplot(aes(x=trip, y=AULC, fill=treatment, color=treatment)) +
  geom_point(size=2, shape=21) + geom_smooth(method = 'lm', se=FALSE) +
  ggtitle('Correlation between cumulative community membership change and cumulative shedding', 
          subtitle = 'correlation stats: RPS pval = 0.02, control pval = 0.31, Bglu pval = 0.42') +
  xlab('Cumulative Bray-Curtis distance (presence/abscence)')



#testse <- cor.test(pig_trips$trip, pig_trips$AULC)
#testse$p.value
#testse$statistic


ggplot(pig_trips, aes(x=trip, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
ggplot(pig_trips, aes(x=D0_2, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
ggplot(pig_trips, aes(x=D2_7, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
ggplot(pig_trips, aes(x=D7_14, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')
ggplot(pig_trips, aes(x=D14_21, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

#ggplot(pig_trips, aes(x=sum, y=AULC, color=treatment)) + geom_point() + geom_smooth(method = 'lm')

pig_trips %>% group_by(treatment) %>% summarise(num=n())



# looking for missing samples


sum(ASV_table[grep('P50D0', rownames(ASV_table)),])
sum(ASV_table[grep('P181D7F', rownames(ASV_table)),])



ggplot(pig_trips, aes(x=treatment, y=trip, fill=treatment)) +
  geom_boxplot() + ylab('Cumulative bray-curtis dissimilarity (each pig)') + geom_jitter(size=2.5,width = 0.2, shape=21)+
  ggtitle('Cumulative change in community structure through Salmonella infection')

