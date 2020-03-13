setwd('~/FS12/prot/')



library(tidyverse)



pull_rsq <- function(model){
  summary(model)$r.squared
}

mod_fun <- function(df){
  lm(data=df, formula = abs405~secs)
}

extract_secs <- function(model){
  
  model %>%
    broom::tidy() %>% 
    filter(term == 'secs') %>% 
    pull(estimate)
}

treats <- read_csv('FS12_treatments.csv')

p1 <- read_csv('NN_SalN_AD_GLU.csv')
p11 <- read_csv('Plate1_adglu_bdmalt.csv')
p11_met <- read_csv('p1_adglu_bdmalt_met.csv')


p11_met <- merge(p11_met, treats, by = 'pig') %>% mutate(pen=pen.y) %>% select(-pen.x, -pen.y)








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

as.numeric(p1_all$Time)
lubridate::seconds(p1_all$Time)

# p1 is necropsy data, prob cec_cont from both the nursury necropsy as well as the salmonella necropsy

p1_all %>% filter(exp == 'FS12b', treatment != 'Malto') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))
p1_all %>% filter(exp == 'FS12a') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))


unique(p1_all$substrate)
unique(p1_all$treatment)
unique(p1_all$exp)
unique(p1_all$substrate)

# P11 is fecal data from D0 and D2 from FS12b only


p11.g <- p11 %>% gather(key=well, value=abs405, -Time) %>% filter(well !='Blank') 
p11_all <- merge(p11.g, p11_met, by = 'well') %>% filter(pig !=101) %>%
  select(-diet, -birth_weight) %>% 
  left_join(treats) %>% 
  mutate(secs=as.numeric(Time))






TEST <- p11_all %>% group_by(pig, day, treatment, substrate) %>% nest() %>% 
  mutate(model=map(data, mod_fun), 
         rate= map_dbl(model, extract_secs),
         rsq=map_dbl(model, pull_rsq)) %>% 
  select(-data, -model)


TEST %>%filter(day==0) %>%  ggplot(aes(x=treatment, y=rate)) +
  geom_boxplot() + geom_text(aes(label=pig))+
  facet_wrap(~substrate, scales = 'free')
  

TEST %>%filter(day==2) %>%  ggplot(aes(x=treatment, y=rate)) +
  geom_boxplot() + geom_text(aes(label=pig))+
  facet_wrap(~substrate, scales = 'free')


TEST$SHED <- ifelse(TEST$pig %in% c(373,321,181,392,97), 'low', 'high')



TEST %>%filter(day==0 & treatment == 'RPS') %>%  ggplot(aes(x=SHED, y=rate)) +
  geom_boxplot() + geom_text(aes(label=pig))+
  facet_wrap(~substrate, scales = 'free')



TEST %>%filter(day==2 & treatment == 'RPS') %>%  ggplot(aes(x=SHED, y=rate)) +
  geom_boxplot() + geom_text(aes(label=pig))+
  facet_wrap(~substrate, scales = 'free')






TEST %>% ungroup() %>% filter(substrate=='adglu'& treatment =='RPS' & day == 0) %>% 
  mutate(day_shed= paste(day, SHED, sep = '_')) %>% 
  ggplot(aes(x=day, y=rate, fill=SHED, group=day_shed)) +
  geom_boxplot() + geom_jitter(width = .2)+
  facet_wrap(~substrate, scales = 'free')



TEST %>% 
  ungroup() %>%
  filter(substrate=='adglu'& treatment =='RPS' & day == 2) %>% 
  ggplot(aes(x=SHED, y=rate, fill=SHED)) +
  geom_boxplot(outlier.color=NA) + geom_jitter(width = .2)






TEST %>% ungroup() %>% filter(substrate=='adglu') %>% 
  group_by(day) %>% nest()
# 
# TEST$data[[190]] %>% ggplot(aes(x=secs, y=abs405)) + geom_line()

# 
# broom::tidy(lm(data = TEST$data[[190]], formula = abs405~secs) )
# 
# 
# check <- summary(lm(data = TEST$data[[190]], formula = abs405~secs))
# 
# check$coefficients
# 




#
  
  
  
  


colnames(p11_all)



unique(p11_all$pig)
p11_all$diet


p11_all %>% filter(day == 0, substrate =='adglu') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))
p11_all %>% filter(day == 2, substrate =='adglu') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))


p11_all %>% filter(day == 0, substrate =='bdmalt') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))
p11_all %>% filter(day == 2, substrate =='bdmalt') %>% ggplot(aes(x=Time, y=abs405, group=well)) + geom_line(aes(color=treatment))

