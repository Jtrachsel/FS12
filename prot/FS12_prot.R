setwd('~/FS12/prot/')



library(tidyverse)




p1 <- read_csv('NN_SalN_AD_GLU.csv')
p11 <- read_csv('Plate1_adglu_bdmalt.csv')
p11_met <- read_csv('p1_adglu_bdmalt_met.csv')


p11_met <- merge(p11_met, treats, by = 'pig') %>% mutate(pen=pen.y) %>% select(-pen.x, -pen.y)








treats <- read_csv('FS12_treatments.csv')
d2t <- read_csv('FS12_diet_treat_map.csv')

FS12a_meta <- read_csv('FS12a_meta.csv')
FS12a_meta <- merge(FS12a_meta, d2t, by = 'diet')

treats <- merge(treats, d2t, by = 'diet')


p1_met <- read_csv('NN_SalN_96to100map.csv')

p1_tmp <- merge(p1_met, FS12a_meta, by = 'pen')
NN_meta <- p1_tmp %>% mutate(pig=pig.y) %>% select(-pig.x, -pig.y)
SN_meta <- merge(p1_met, treats, by='pig') %>% mutate(pen=pen.y) %>% select(-pen.x, -pen.y)











p1.g <- p1 %>% gather(key=well, value=abs405, -Time) %>% filter(well !='Blank')
p1.g$substrate <- 'AD_glu'

NN_meta <- NN_meta[,match(colnames(SN_meta) ,colnames(NN_meta))]
NN_meta$exp <- 'FS12a'
SN_meta$exp <- 'FS12b'


N_meta <- rbind(NN_meta, SN_meta)


p1_all <- merge(p1.g, N_meta, by='well')

p1_all$Time
lubridate::seconds(p1_all$Time)


p1_all %>% filter(exp == 'FS12b', treatment != 'Malto') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))
p1_all %>% filter(exp == 'FS12a') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))




p11.g <- p11 %>% gather(key=well, value=abs405, -Time) %>% filter(well !='Blank')
p11_all <- merge(p11.g, p11_met, by = 'well')


p11_all %>% filter(day == 0, substrate =='adglu') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))
p11_all %>% filter(day == 2, substrate =='adglu') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))


p11_all %>% filter(day == 0, substrate =='bdmalt') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))
p11_all %>% filter(day == 2, substrate =='bdmalt') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))

