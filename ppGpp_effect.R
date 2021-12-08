library(RcppCNPy)
library(tidyverse)
library(dygraphs)
library(genbankr)
library(gridExtra)
library(DescTools)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicFiles)
library(edgeR)
library(limma)
library(EBSeq)
library(ggnewscale)

#Loading npy arrays of effect of ppGpp----
mean_ppGpp_effect = npyLoad('mean_ppGpp_effect.npy')
mean_ppGpp_effect[is.na(mean_ppGpp_effect)]=0

sd_ppGpp_effect = npyLoad('sd_ppGpp_effect.npy')
sd_ppGpp_effect[is.na(sd_ppGpp_effect)]=0

ppGpp_effect_df = tibble(position=1:length(mean_ppGpp_effect)*10,
                         mean = mean_ppGpp_effect,
                         sd = sd_ppGpp_effect)

save(ppGpp_effect_df, file='ppGpp_effect.Rdata')

mean_ppGpp_effect_linear = npyLoad('mean_ppGpp_effect_linear.npy')
mean_ppGpp_effect_linear[is.na(mean_ppGpp_effect_linear)]=0

sd_ppGpp_effect_linear = npyLoad('sd_ppGpp_effect_linear.npy')
sd_ppGpp_effect_linear[is.na(sd_ppGpp_effect_linear)]=0

ppGpp_effect_linear_df = tibble(position=1:length(mean_ppGpp_effect_linear)*10,
                         mean = mean_ppGpp_effect_linear,
                         sd = sd_ppGpp_effect_linear)

save(ppGpp_effect_linear_df, file='ppGpp_effect_linear.Rdata')

#Plotting effect of ppGpp ----

load(file='ppGpp_effect.Rdata')
load("annotated_CDSs.RData")
load("PurBox_locations.RData")
source('helperFunctions.R')

gene_start = CDSs %>% filter(gene=='purE') %>% .$start
gene_end = CDSs %>% filter(gene=='purE') %>% .$end

line_color <- c('mean'='firebrick4')
strand_colors <- c("+"="gray80","-"="gray80")

genome_plot = plotEffect(CDSs, 
                         ppGpp_effect_df,
                         1, 4.215e6,lineAlpha=1, 
                         plotFeatures=F, plotPurBox=T, plotError=F)

pur_plot = plotEffect(CDSs, 
                      ppGpp_effect_df,
                      gene_start-1e3, gene_start+1e3,lineAlpha=1, 
                      plotFeatures=T, plotPurBox=T, plotError=T)

plotEffect(CDSs, 
           ppGpp_effect_df,
           1.995e6, 1.9975e6,lineAlpha=1, 
           plotFeatures=T, plotPurBox=T, plotError=T,
           useLog=T)

#Plotting ppGpp effect at each PurR binding site ----

load(file='ppGpp_effect.Rdata')
load("annotated_CDSs.RData")
load("PurBox_locations.RData")
source('helperFunctions.R')
strand_colors <- c("+"="gray80","-"="gray80")
peaks_pm_750bp = read.csv(file='Peaks_pm_750bp.csv')
peaks_pm_500bp = read.csv(file='Peaks_pm_500bp.csv')

for(i in c(seq_len(nrow(peaks_pm_750bp)))){
  p = plotEffectFigure(CDSs, ppGpp_effect_df, 
                       peaks_pm_750bp$plotStart[i], 
                       peaks_pm_750bp$plotEnd[i], 
                       plotFeatures=T, 
                       plotPurBox=T, 
                       plotError=T)
  pdfName = paste0(peaks_pm_750bp$peak_name[i],'_1500bp_window_effect',".pdf")
  pdf(pdfName, height=1.6, width=1.6)
  print(p)
  dev.off()
}

for(i in c(seq_len(nrow(peaks_pm_750bp)))){
  p = plotEffectFigure(CDSs, ppGpp_effect_df, 
                       peaks_pm_750bp$plotStart[i], 
                       peaks_pm_750bp$plotEnd[i], 
                       plotFeatures=T, 
                       plotPurBox=T, 
                       plotError=T,
                       useLog=F)
  pdfName = paste0(peaks_pm_750bp$peak_name[i],'_1500bp_window_effect_nolog',".pdf")
  pdf(pdfName, height=1.6, width=1.6)
  print(p)
  dev.off()
}

dev.off()


plot = 
  ppGpp_effect_df %>% ggplot() +
  geom_line(aes(x=position, y=mean), size=1) +
  theme_classic()


#Plotting ppGpp effect in linear ----
load(file='ppGpp_effect_linear.Rdata')
load("annotated_CDSs.RData")
load("PurBox_locations.RData")
source('helperFunctions.R')
strand_colors <- c("+"="gray80","-"="gray80")
line_color = 'firebrick'
peaks_pm_750bp = read.csv(file='Peaks_pm_750bp.csv')
peaks_pm_500bp = read.csv(file='Peaks_pm_500bp.csv')

for(i in c(seq_len(nrow(peaks_pm_750bp)))){
  p = plotEffectFigure(CDSs, ppGpp_effect_linear_df, 
                       peaks_pm_750bp$plotStart[i], 
                       peaks_pm_750bp$plotEnd[i], 
                       plotFeatures=T, 
                       plotPurBox=T, 
                       plotError=T)
  pdfName = paste0(peaks_pm_750bp$peak_name[i],'_1500bp_window_effect_nolog',".pdf")
  pdf(pdfName, height=1.6, width=1.6)
  print(p)
  dev.off()
}

for(i in c(seq_len(nrow(peaks_pm_750bp)))){
  p = plotEffectFigure(CDSs, ppGpp_effect_df, 
                       peaks_pm_750bp$plotStart[i], 
                       peaks_pm_750bp$plotEnd[i], 
                       plotFeatures=T, 
                       plotPurBox=T, 
                       plotError=T,
                       useLog=F)
  pdfName = paste0(peaks_pm_750bp$peak_name[i],'_1500bp_window_effect_nolog',".pdf")
  pdf(pdfName, height=1.6, width=1.6)
  print(p)
  dev.off()
}

dev.off()


plot = 
  ppGpp_effect_df %>% ggplot() +
  geom_line(aes(x=position, y=mean), size=1) +
  theme_classic()

#Determining effect of ppGpp from mean enrichment datasets----

load(file='p0_rhx_00_IDR.Rdata')
load(file='p0_rhx_10_IDR.Rdata')
load(file='wt_rhx_00_IDR.Rdata')
load(file='wt_rhx_10_IDR.Rdata')


ppGpp_effect = (wt_rhx_10 - wt_rhx_00) - (p0_rhx_10 - p0_rhx_00)

ppGpp_effect_stuff = ppGpp_effect %>% mutate('position' = (1:nrow(ppGpp_effect))*10) %>%
                                      group_by(position) %>%
                                      mutate('mean' = mean())
                                            
      
plot1 = 
  ppGpp_effect_stuff %>% ggplot() +
                         geom_line(aes(x=position, y=wt_rep1_rhx_10), size=1) +
                         theme_classic()

plot2 = 
  ppGpp_effect_stuff %>% ggplot() +
  geom_line(aes(x=position, y=wt_rep2_rhx_10), size=1) +
  theme_classic()

plot3 = 
  ppGpp_effect_stuff %>% ggplot() +
  geom_line(aes(x=position, y=wt_rep3_rhx_10), size=1) +
  theme_classic()

save(p0_rhx_00, file='p0_rhx_00_IDR.Rdata')
save(p0_rhx_10, file='p0_rhx_10_IDR.Rdata')
save(wt_rhx_00, file='wt_rhx_00_IDR.Rdata')
save(wt_rhx_10, file='wt_rhx_10_IDR.Rdata')

source('helperFunctions.R')

#Pulling mean and sd of log2 enrichment for PurR ChIP peaks ----
peak_locations = read.csv(file="PurR_chip_peak_locations.csv")

peak_locations = peak_locations %>% mutate("peak_number"=c(1:15))
save(peak_locations, file="peak_locations.RData")
load("peak_locations.RData")

peak_data = ppGpp_effect_df %>%
  dplyr::filter(position > peak_locations$start, position < peak_locations$end) %>%
  mutate(peak_number = case_when(position > peak_locations[1,1] & position < peak_locations[1,2] ~ peak_locations[1,3],
                                 position > peak_locations[2,1] & position < peak_locations[2,2] ~ peak_locations[2,3],
                                 position > peak_locations[3,1] & position < peak_locations[3,2] ~ peak_locations[3,3],
                                 position > peak_locations[4,1] & position < peak_locations[4,2] ~ peak_locations[4,3],
                                 position > peak_locations[5,1] & position < peak_locations[5,2] ~ peak_locations[5,3],
                                 position > peak_locations[6,1] & position < peak_locations[6,2] ~ peak_locations[6,3],
                                 position > peak_locations[7,1] & position < peak_locations[7,2] ~ peak_locations[7,3],
                                 position > peak_locations[8,1] & position < peak_locations[8,2] ~ peak_locations[8,3],
                                 position > peak_locations[9,1] & position < peak_locations[9,2] ~ peak_locations[9,3],
                                 position > peak_locations[10,1] & position < peak_locations[10,2] ~ peak_locations[10,3],
                                 position > peak_locations[11,1] & position < peak_locations[11,2] ~ peak_locations[11,3],
                                 position > peak_locations[12,1] & position < peak_locations[12,2] ~ peak_locations[12,3],
                                 position > peak_locations[13,1] & position < peak_locations[13,2] ~ peak_locations[13,3],
                                 position > peak_locations[14,1] & position < peak_locations[14,2] ~ peak_locations[14,3],
                                 position > peak_locations[15,1] & position < peak_locations[15,2] ~ peak_locations[15,3]))

peak_data_max = peak_data %>%
  dplyr::group_by(peak_number) %>%
  dplyr::filter(mean==max(mean))

save(peak_data_max, file="final_peak_locations.RData")
write.csv(peak_data_max, file='ppGpp_effect_at_peak_locations.csv')

peak_names = read.csv(file="peak_names.csv")

peak_data_max = peak_data_max %>%
  left_join(peak_names, by="peak_number")

peak_data_max$peak_name <- factor(peak_data_max$peak_name, levels = c('yaaD-yaaE',
                                                                      'purR-yabJ',
                                                                      'pbuG',
                                                                      'purE',
                                                                      'ykbA/mhqA',
                                                                      'pycA',
                                                                      'yolJ',
                                                                      'xpt-pbuX',
                                                                      'nusB-folD',
                                                                      'pbuO',
                                                                      'guaC',
                                                                      'oppBA-yvaV',
                                                                      'opuCA-yvbF',
                                                                      'glyA',
                                                                      'purA'))

plot = peak_data_max %>%
  ggplot(aes(x=peak_name, y=mean)) +
  geom_bar(stat="identity", position = 'dodge') +
  geom_errorbar(aes(x=peak_name, 
                    ymin=mean-sd, 
                    ymax=mean+sd),
                position = position_dodge(0.9),
                width=0.3) +
  xlab('Gene/operon') +
  ylab('(p)ppGpp-dependent PurR enrichment (log2)') +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, face='bold'),
        axis.line = element_line(size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


ggsave(plot=plot, filename='ppGpp_effect_at_peaks_bargraph.pdf', width=4, heigh=3)

#Pulling mean and sd of linear enrichment for PurR ChIP peaks----
peak_locations = read.csv(file="PurR_chip_peak_locations.csv")

peak_locations = peak_locations %>% mutate("peak_number"=c(1:15))
save(peak_locations, file="peak_locations.RData")
load("peak_locations.RData")

peak_data = ppGpp_effect_linear_df %>%
  dplyr::filter(position > peak_locations$start, position < peak_locations$end) %>%
  mutate(peak_number = case_when(position > peak_locations[1,1] & position < peak_locations[1,2] ~ peak_locations[1,3],
                                 position > peak_locations[2,1] & position < peak_locations[2,2] ~ peak_locations[2,3],
                                 position > peak_locations[3,1] & position < peak_locations[3,2] ~ peak_locations[3,3],
                                 position > peak_locations[4,1] & position < peak_locations[4,2] ~ peak_locations[4,3],
                                 position > peak_locations[5,1] & position < peak_locations[5,2] ~ peak_locations[5,3],
                                 position > peak_locations[6,1] & position < peak_locations[6,2] ~ peak_locations[6,3],
                                 position > peak_locations[7,1] & position < peak_locations[7,2] ~ peak_locations[7,3],
                                 position > peak_locations[8,1] & position < peak_locations[8,2] ~ peak_locations[8,3],
                                 position > peak_locations[9,1] & position < peak_locations[9,2] ~ peak_locations[9,3],
                                 position > peak_locations[10,1] & position < peak_locations[10,2] ~ peak_locations[10,3],
                                 position > peak_locations[11,1] & position < peak_locations[11,2] ~ peak_locations[11,3],
                                 position > peak_locations[12,1] & position < peak_locations[12,2] ~ peak_locations[12,3],
                                 position > peak_locations[13,1] & position < peak_locations[13,2] ~ peak_locations[13,3],
                                 position > peak_locations[14,1] & position < peak_locations[14,2] ~ peak_locations[14,3],
                                 position > peak_locations[15,1] & position < peak_locations[15,2] ~ peak_locations[15,3]))

peak_data_max = peak_data %>%
  dplyr::group_by(peak_number) %>%
  dplyr::filter(mean==max(mean))

save(peak_data_max, file="final_peak_locations.RData")
write.csv(peak_data_max, file='ppGpp_effect_linear_at_peak_locations.csv')

peak_names = read.csv(file="peak_names.csv")

peak_data_max = peak_data_max %>%
  left_join(peak_names, by="peak_number")

peak_data_max$peak_name <- factor(peak_data_max$peak_name, levels = c('yaaD-yaaE',
                                                                      'purR-yabJ',
                                                                      'pbuG',
                                                                      'purE',
                                                                      'ykbA/mhqA',
                                                                      'pycA',
                                                                      'yolJ',
                                                                      'xpt-pbuX',
                                                                      'nusB-folD',
                                                                      'pbuO',
                                                                      'guaC',
                                                                      'oppBA-yvaV',
                                                                      'opuCA-yvbF',
                                                                      'glyA',
                                                                      'purA'))

plot = peak_data_max %>%
  ggplot(aes(x=peak_name, y=mean)) +
  geom_bar(stat="identity", position = 'dodge') +
  geom_errorbar(aes(x=peak_name, 
                    ymin=mean-sd, 
                    ymax=mean+sd),
                position = position_dodge(0.9),
                width=0.3) +
  xlab('Gene/operon') +
  ylab('(p)ppGpp-dependent PurR enrichment (log2)') +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, face='bold'),
        axis.line = element_line(size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


ggsave(plot=plot, filename='ppGpp_effect_linear_at_peaks_bargraph.pdf', width=4, heigh=3)
