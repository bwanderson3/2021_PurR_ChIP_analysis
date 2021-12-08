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

options(stringsAsFactors = F)
source('helperFunctions.R')

set.seed(1234)

# working with bootstrapped count data, reading in the medians and upper/lower confidence limits ----
# we already summarized the 100 count bootstraps with 'summarize_all_sample_counts.py'.

data_direc = "Y:/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data"
data_direc = "/mnt/jadelab/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data"
info_file = file.path(data_direc,"PurR_chip_and_input_pairs.csv")

info_df = read_csv(info_file)
pair_info = info_df %>% 
  filter(!is.na(pair_description)) %>%
  dplyr::select(pair, pair_description)
info_df %>% head
info_df = info_df %>%
  mutate(median_file_path = file.path(data_direc, 
                                      Dir,
                                      'jawang',
                                      paste0("Sample_", Sample),
                                      paste0("Sample_", Sample, "_summary_median.npy")),
         upper_cl_file_path = file.path(data_direc, 
                                      Dir,
                                      'jawang',
                                      paste0("Sample_", Sample),
                                      paste0("Sample_", Sample, "_summary_maxci.npy")),
         median = purrr::map(median_file_path, npyLoad))

# Transforming bootstrapped cpm to log2 fold enrichment----

info_df_nest = info_df %>%
  nest(-pair) %>%
  left_join(pair_info) %>%
  mutate(enrichment = purrr::map(data, get_ChIP_enrichment),
         position = purrr::map(enrichment, function(x){1:length(x)*10}))

enrichment_df = info_df_nest %>%
  mutate(j = str_split(pair_description, "_rep\\d_"), 
         group = purrr::map_chr(j, function(x){paste(x[1],x[2],sep="_")})) %>%
  unnest(enrichment, position) %>%
  mutate(log2_enrich = log2(enrichment),
         log2_enrich_finite = ifelse(is.infinite(log2_enrich), 0, log2_enrich))

#Processing enrichment_df for use with calculate_idr.R

enrichment_df_wide = enrichment_df %>% select(pair_description, log2_enrich_finite, position)

enrichment_df_wide = enrichment_df_wide %>% spread(pair_description, log2_enrich_finite)

enrichment_df_wide = enrichment_df_wide %>% select(-'position')

p0_rhx_00 = enrichment_df_wide %>% select(p0_rep1_rhx_00, p0_rep2_rhx_00, p0_rep3_rhx_00)
p0_rhx_10 = enrichment_df_wide %>% select(p0_rep1_rhx_10, p0_rep2_rhx_10, p0_rep3_rhx_10)
wt_rhx_00 = enrichment_df_wide %>% select(wt_rep1_rhx_00, wt_rep2_rhx_00, wt_rep3_rhx_00)
wt_rhx_10 = enrichment_df_wide %>% select(wt_rep1_rhx_10, wt_rep2_rhx_10, wt_rep3_rhx_10)
wt_acgu = enrichment_df_wide %>% select(wt_rep1_acgu, wt_rep2_acgu, wt_rep3_acgu)

save(p0_rhx_00, file='p0_rhx_00_IDR.Rdata')
save(p0_rhx_10, file='p0_rhx_10_IDR.Rdata')
save(wt_rhx_00, file='wt_rhx_00_IDR.Rdata')
save(wt_rhx_10, file='wt_rhx_10_IDR.Rdata')
save(wt_acgu, file='wt_acgu_IDR.Rdata')

#Note that infinite values have been excluded!

summary_df = enrichment_df %>%
  group_by(group,position) %>%
  summarize(mean_log2_enrich = mean(log2_enrich_finite),
            sd_log2_enrich = sd(log2_enrich))

save(summary_df, file="summary_df.RData")
load("summary_df.RData")

#Summarizing the bootstrapped data (instead of transforming to log2 fold enrichment) ----

#This is the same manipulation that we did in the first section (to create info_df)
info_df = read_csv(info_file)
pair_info = info_df %>% 
  filter(!is.na(pair_description)) %>%
  dplyr::select(pair, pair_description)
info_df %>% head
info_df = info_df %>%
  mutate(median_file_path = file.path(data_direc, 
                                      Dir,
                                      'jawang',
                                      paste0("Sample_", Sample),
                                      paste0("Sample_", Sample, "_summary_median.npy")),
         upper_cl_file_path = file.path(data_direc, 
                                        Dir,
                                        'jawang',
                                        paste0("Sample_", Sample),
                                        paste0("Sample_", Sample, "_summary_maxci.npy")),
         median = purrr::map(median_file_path, npyLoad))

#Converting the info_df to working with cpm
info_df_cpm = info_df %>%
              mutate(position = purrr::map(median, function(x){1:length(x)*10}))

summary_cpm = info_df_cpm %>%
              mutate(replicate = str_extract(info_df$Description, "rep."))

summary_cpm = summary_cpm %>% subset(select = c("Description", "median", "replicate", "position"))

summary_cpm = summary_cpm %>% mutate(Description = str_remove(Description, "_rep."))

summary_cpm_analyzed = summary_cpm %>%
                       unnest() %>%
                       group_by(Description, position) %>%
                       summarize(mean_cpm = mean(median),
                                 sd_cpm = sd(median))

#Plotting cpm from triplicate (no smoothed data); these data were used to make genome-wide figures----
save(summary_cpm_analyzed, file="summary_cpm_analyzed.RData")
load("summary_cpm_analyzed.RData")

#Load CDS file (with gene names, etc)
save(CDSs, file="annotated_CDSs.RData")
load("annotated_CDSs.RData")

#Load dataframe with PurBox locations - so they can be mapped if desired
purbox_df = read.csv(file="purbox_locations.csv", colClasses=c("NULL",NA,NA))
save(purbox_df, file="PurBox_locations.RData")
load("PurBox_locations.RData")

#To plot specific samples, change the matrix below with the appropriate sample names
plot_samples = c("wt_acgu_chip")
plot_samples = c("wt_rhx_00_chip","wt_rhx_10_chip")
plot_samples = c('wt_rhx_00_chip')
plot_samples = c('p0_rhx_00_chip','p0_rhx_10_chip')
plot_samples = c("wt_acgu_chip", "wt_rhx_00_chip")

#To find a given gene, replace gene below and execute these functions
gene_start = CDSs %>% filter(gene=='purA') %>% .$start
gene_end = CDSs %>% filter(gene=='purA') %>% .$end

#In helperFunctions.R, I added lines to manually specify colors (fill the colors below).
#The palettes for strands and samples are defined here.
strand_colors <- c("+"="gray80","-"="gray80")

#manually define sample colors
sample_colors <- c(wt_rhx_00_chip = "gray60",
                   wt_rhx_10_chip = "navy",
                   p0_rhx_00_chip = "plum3",
                   p0_rhx_10_chip = "darkorchid3",
                   wt_acgu_chip = "darkorchid4")

#Plot every 3rd data point to make smaller figures
plot_p0_rhx_3rd = plotLocusCPM_avg(CDSs, 
                            summary_cpm_analyzed %>% dplyr::filter(Description %in% plot_samples,
                                                                   position %% 30 == 0),
                            1, 4.215e6,lineAlpha=1,lineSize=1,plotFeatures=F, plotPurBox=F)

plot_wt_rhx_3rd = plotLocusCPM_avg(CDSs, 
                                   summary_cpm_analyzed %>% dplyr::filter(Description %in% plot_samples,
                                                                          position %% 30 == 0),
                                   1, 4.215e6,lineAlpha=1,lineSize=1,plotFeatures=F, plotPurBox=F)

plot_wt_acgu_3rd = plotLocusCPM_avg(CDSs, 
                                   summary_cpm_analyzed %>% dplyr::filter(Description %in% plot_samples,
                                                                          position %% 30 == 0),
                                   1, 4.215e6,lineAlpha=1,lineSize=1,plotFeatures=F, plotPurBox=F)

#Plot every other data point to save some computing power.
plot_2nd = plotLocusCPM_avg(CDSs, 
                             summary_cpm_analyzed %>% dplyr::filter(Description %in% plot_samples,
                                                                    position %% 20 == 0),
                             1, 4.215e6,lineAlpha=1,lineSize=1,plotFeatures=F, plotPurBox=F)

ggsave(plot=plot_2nd, filename='genome_wide_wt_preRHX.eps', width=3.5, height=1.5)
ggsave(plot=plot_2nd, filename='genome_wide_wt_acgu.eps', width=3.5, height=1.5)
ggsave(plot=plot_3rd, filename='genome_wide_wt_acgu_3rd.eps', width=3.5, height=1.5)
ggsave(plot=plot_2nd, filename='genome_wide_wt_+-acgu.eps', width=3.5, height=1.5)
ggsave(plot=plot_wt_rhx_3rd, filename='genome_wide_wt_preRHX.eps', width=3.5, height=1.5)

#Plot all data
plot_both = plotLocusCPM_avg(CDSs, 
                            summary_cpm_analyzed %>% dplyr::filter(Description %in% plot_samples),
                            1, 4.215e6,lineAlpha=1, plotFeatures=F, plotPurBox=F)

plot_acgu = plotLocusCPM_avg(CDSs, 
                             summary_cpm_analyzed %>% dplyr::filter(Description %in% plot_samples),
                             1, 4.215e6, ylims=c(0:200), lineAlpha=1, plotFeatures=F, plotPurBox=F)

plot_p0_rhx = plotLocusCPM_avg(CDSs, 
                               summary_cpm_analyzed %>% dplyr::filter(Description %in% plot_samples),
                               1, 4.215e6,lineAlpha=1, plotFeatures=F, plotPurBox=F)

plot_beforeRHX = plotLocusCPM_avg(CDSs, 
                            summary_cpm_analyzed %>% dplyr::filter(Description=='wt_rhx_00_chip'),
                            1, 4.215e6,lineAlpha=1,plotFeatures=F, plotPurBox=F)

plot_afterRHX = plotLocusCPM_avg(CDSs, 
                            summary_cpm_analyzed %>% dplyr::filter(Description=='wt_rhx_10_chip'),
                            1, 4.215e6,lineAlpha=1,plotFeatures=F, plotPurBox=F)

plotLocusCPM_avg(CDSs, 
                 summary_cpm_analyzed %>% 
                   dplyr::filter(Description %in% plot_samples),
                 4.1e6, 4.2e6,lineAlpha=1,plotFeatures=T, plotPurBox=T, plotError=T)

source('helperFunctions.R')

ggsave(plot=plot_both, filename='genome_view_both.png', width=5, height=1.5, dpi=300)
ggsave(plot=plot_beforeRHX, filename='genome_view_beforerhx.png', width=5, height=1.5, dpi=300)
ggsave(plot=plot_afterRHX, filename='genome_view_afterrhx.png', width=5, height=1.5, dpi=300)
ggsave(plot=plot_acgu, filename='acgu.eps', width=3.5, heigh=1.5)
ggsave(plot=plot_p0_rhx, filename='p0_rhx.eps', width=3.5, heigh=1.5)
ggsave(plot=plot_p0_rhx_3rd, filename='p0_rhx_3rd.eps', width=3.5, heigh=1.5)
ggsave(plot=plot_wt_rhx_3rd, filename='wt_rhx_3rd.eps', width=3.5, heigh=1.5)
ggsave(plot=plot_wt_acgu_3rd, filename='wt_acgu_3rd.eps', width=3.5, heigh=1.5)

#Trying to plot circular chromosome (unsuccessful)
summary_cpm_analyzed %>% filter(Description==c('wt_rhx_00_chip','wt_rhx_10_chip')) %>%
  ggplot(aes(x=position, y=mean_cpm, fill=Description)) + 
  geom_col(width=0.1) +
  scale_color_manual(values=sample_colors) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  coord_polar(theta='x', start=0, direction = 1, clip = 'on') +
  aes(group=rev(Description))


# Analyzing dataset smoothed with NGS_smoothing.py (smoothed with 50bp windows) ----
#Data smoothed with 50bp windows were used to make enrichment figures
data_direc = "Y:/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data"
data_direc = "/mnt/jadelab/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data"

info_file = file.path(data_direc,"PurR_chip_and_input_pairs.csv")
info_df = read_csv(info_file)
pair_info = info_df %>% 
  filter(!is.na(pair_description)) %>%
  dplyr::select(pair, pair_description)

save(info_df, file="info_df.RData")
load("info_df.RData")

save(pair_info, file="pair_info.RData")
load("pair_info.RData")

smoothed_file = file.path(data_direc, "smoothed_data.csv")
smoothed_df = read_csv(smoothed_file) [, 2:26]

#Contains smoothed bootstrapped cpm for each sample
save(smoothed_df, file="smoothed_df.RData")
load("smoothed_df.RData")

#Make it a long tibble; remove "Sample_"
smoothed_long = smoothed_df %>% tidyr::gather(Sample, median)
smoothed_long = smoothed_long %>% mutate(Sample = stringr::str_remove(smoothed_long$Sample, "Sample_"))
smoothed_long$Sample <- as.numeric(as.character(smoothed_long$Sample))

save(smoothed_long, file="smoothed_long.RData")
load("smoothed_long.RData")

smoothed_long = smoothed_long %>%
  group_by(Sample) %>%
  summarize(median = list(median))

smoothed_info = smoothed_long %>%
  left_join(info_df)

save(smoothed_info, file="smoothed_info.RData")
load("smoothed_info.RData")

smoothed_df_nest = smoothed_info %>%
  nest(-pair) %>%
  left_join(pair_info) %>%
  mutate(enrichment = purrr::map(data, get_ChIP_enrichment),
         position = purrr::map(enrichment, function(x){1:length(x)*10}))

enrichment_smoothed_df = smoothed_df_nest %>%
  mutate(j = str_split(pair_description, "_rep\\d_"), 
         group = purrr::map_chr(j, function(x){paste(x[1],x[2],sep="_")})) %>%
  unnest(enrichment, position) %>%
  #filter(!(enrichment=="NaN")) %>%
  mutate(enrich = log2(enrichment),
         log2_enrich_finite = ifelse(is.infinite(log2_enrich), 0, log2_enrich))

summary_smoothed = enrichment_smoothed_df %>%
  group_by(group,position) %>%
  summarize(mean_log2_enrich = mean(log2_enrich_finite),
            sd_log2_enrich = sd(log2_enrich))

summary_smoothed = summary_smoothed %>% filter(!(mean_log2_enrich == "NaN"))

save(summary_smoothed, file="summary_smoothed.RData")
load("summary_smoothed.RData")

# Analyzing dataset smoothed with NGS_smoothing.py (smoothed with 100bp windows) ----
#These data were not used to make figures - all figures were made with 50bp smoothed data
data_direc = "Y:/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data"
data_direc = "/mnt/jadelab/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data"

info_file = file.path(data_direc,"PurR_chip_and_input_pairs.csv")
info_df = read_csv(info_file)
pair_info = info_df %>% 
  filter(!is.na(pair_description)) %>%
  dplyr::select(pair, pair_description)

save(info_df, file="info_df.RData")
load("info_df.RData")

save(pair_info, file="pair_info.RData")
load("pair_info.RData")

smoothed_file_100bp = file.path(data_direc, "smoothed_data_100bp.csv")
smoothed_df_100bp = read_csv(smoothed_file_100bp) [, 2:26]

#Contains smoothed bootstrapped cpm for each sample
save(smoothed_df, file="smoothed_df.RData")
load("smoothed_df.RData")

#Make it a long tibble; remove "Sample_"
smoothed_long_100bp = smoothed_df_100bp %>% tidyr::gather(Sample, median)
smoothed_long_100bp = smoothed_long_100bp %>% mutate(Sample = stringr::str_remove(smoothed_long_100bp$Sample, "Sample_"))
smoothed_long_100bp$Sample <- as.numeric(as.character(smoothed_long_100bp$Sample))

save(smoothed_long, file="smoothed_long.RData")
load("smoothed_long.RData")

smoothed_long_100bp = smoothed_long_100bp %>%
  group_by(Sample) %>%
  summarize(median = list(median))

smoothed_info_100bp = smoothed_long_100bp %>%
  left_join(info_df)

save(smoothed_info, file="smoothed_info.RData")
load("smoothed_info.RData")

smoothed_df_nest_100bp = smoothed_info_100bp %>%
  nest(-pair) %>%
  left_join(pair_info) %>%
  mutate(enrichment = purrr::map(data, get_ChIP_enrichment),
         position = purrr::map(enrichment, function(x){1:length(x)*10}))

enrichment_smoothed_df_100bp = smoothed_df_nest_100bp %>%
  mutate(j = str_split(pair_description, "_rep\\d_"), 
         group = purrr::map_chr(j, function(x){paste(x[1],x[2],sep="_")})) %>%
  unnest(enrichment, position) %>%
  #filter(!(enrichment=="NaN")) %>%
  mutate(log2_enrich = log2(enrichment),
         log2_enrich_finite = ifelse(is.infinite(log2_enrich), 0, log2_enrich))

summary_smoothed_100bp = enrichment_smoothed_df_100bp %>%
  group_by(group,position) %>%
  summarize(mean_log2_enrich = mean(log2_enrich_finite),
            sd_log2_enrich = sd(log2_enrich))

summary_smoothed_100bp = summary_smoothed_100bp %>% filter(!(mean_log2_enrich == "NaN"))

save(summary_smoothed_100bp, file="summary_smoothed_100bp.RData")
load("summary_smoothed_100bp.RData")

#Analyzing 50bp smoothed windows of enrichment without log2 transform ----

load("pair_info.RData")
load("smoothed_info.RData")

smoothed_df_nest = smoothed_info %>%
  nest(-pair) %>%
  left_join(pair_info) %>%
  mutate(enrichment = purrr::map(data, get_ChIP_enrichment),
         position = purrr::map(enrichment, function(x){1:length(x)*10}))

enrichment_smoothed_50bp_nolog2 = smoothed_df_nest %>%
  mutate(j = str_split(pair_description, "_rep\\d_"), 
         group = purrr::map_chr(j, function(x){paste(x[1],x[2],sep="_")})) %>%
  unnest(enrichment, position) %>%
  #filter(!(enrichment=="NaN")) %>%
  mutate(enrich = enrichment)

summary_smoothed_nolog2 = enrichment_smoothed_50bp_nolog2 %>%
  group_by(group,position) %>%
  summarize(mean_enrich = mean(enrich),
            sd_enrich = sd(enrich))

summary_smoothed_nolog2 = summary_smoothed_nolog2 %>% filter(!(mean_enrich == "NaN"))

save(summary_smoothed_nolog2, file="summary_smoothed_50bp_linear.RData")
load("summary_smoothed_50bp_linear.RData")

#Plotting 50bp smoothed linear enrichment ----

load("summary_smoothed_50bp_linear.RData")
load("annotated_CDSs.RData")
load("PurBox_locations.RData")

#To plot specific samples, change the matrix below with the appropriate sample names
plot_samples = c("wt_acgu")
plot_samples = c("wt_rhx_00",
                 "wt_rhx_10", 
                 "p0_rhx_00", 
                 "p0_rhx_10")
plot_samples = c("wt_rhx_00","wt_rhx_10")
plot_samples = c("p0_rhx_00", "p0_rhx_10")
plot_samples = c("wt_rhx_10", "p0_rhx_10")
plot_samples = c("wt_rhx_00", "p0_rhx_00")

#To find a given gene, replace gene below and execute these functions
gene_start = CDSs %>% filter(gene=='purE') %>% .$start
gene_end = CDSs %>% filter(gene=='purE') %>% .$end

#In helperFunctions.R, I added lines to manually specify colors.
#The palettes for strands and samples are defined here.
strand_colors <- c("+"="gray80","-"="gray80")

sample_colors <- c(wt_rhx_00 = "gray50",
                   wt_rhx_10 = "navy",
                   p0_rhx_00 = "plum3",
                   p0_rhx_10 = "darkorchid3",
                   wt_acgu = "darkorchid4")

plotEnrichLinear(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
           1, 4.2e6, plotFeatures=F, plotPurBox=F, plotError=F)

plotEnrichLinear(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
           gene_end-3e3, gene_end+1.4e4, plotFeatures=TRUE, plotPurBox=F, plotError=T)

plotEnrichLinear(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
           gene_start-0.75e3, gene_start+0.1e3, plotFeatures=T, plotPurBox=T, plotError=T)

#PlotEnrichFigure changes size parameters to create a figure that shows up better in PDF/Illustrator

plotEnrichFigure(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
                 gene_end-1.5e3, gene_end+2e3, useLog=F, plotFeatures=TRUE, plotPurBox=T, plotError=T)

plotEnrichFigure(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
                 1, 4.6e6, useLog=F, plotFeatures=F, plotPurBox=F, plotError=F)

pur_operon_wide = plotEnrichFigure(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
                                   gene_end-2e3, gene_end+1.25e4, plotFeatures=TRUE, plotPurBox=F, plotError=T)

pur_operon_zoom = plotEnrichFigure(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
                                   gene_start-1.1e3, gene_start+1.25e3, plotFeatures=TRUE, plotPurBox=T, plotError=T)

pur_operon_zoom_more = plotEnrichFigure(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
                                        gene_start-0.6e3, gene_start+0.1e3, plotFeatures=TRUE, plotPurBox=T, plotError=T)

WT_ACGU_genome = plotEnrichFigure(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
                                  1, 4.2e6, ylims=c(0,200), plotFeatures=F, plotPurBox=F, plotError=F)

source('helperFunctions.R')

ggsave(plot=pur_operon_wide, filename='pur_operon_wide.pdf', width=3.5, height=1.5)

ggsave(plot=pur_operon_zoom_more, filename='pur_operon_zoom_just_peak.pdf', width=2, height=1.25)

ggsave(plot=WT_ACGU_genome, filename='WT_ACGU_genome.pdf', width=3, height=1.5)

tiff("purbox_offset.tif", width=4, height=4, units='in', res=600)
dev.off()


#For plotting multiple plots with for loop

peaks_pm_750bp = read.csv(file='Peaks_pm_750bp.csv')

#Plot 1.5kb windows with acgu data
for(i in c(seq_len(nrow(peaks_pm_750bp)))){
  p = plotEnrichFigure(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
                       peaks_pm_750bp$plotStart[i], peaks_pm_750bp$plotEnd[i], plotFeatures=T, plotPurBox=T, plotError=T)
  pdfName = paste0(peaks_pm_750bp$peak_name[i],'_1500bp_window_acgu',".pdf")
  pdf(pdfName, height=1.6, width=1.6)
  print(p)
  dev.off()
}

#Plot 1.5kb windows with wt and p0 data
for(i in c(seq_len(nrow(peaks_pm_750bp)))){
  p = plotEnrichFigure(CDSs, summary_smoothed_nolog2 %>% dplyr::filter(group %in% plot_samples), 
                       peaks_pm_750bp$plotStart[i], peaks_pm_750bp$plotEnd[i], plotFeatures=T, plotPurBox=F, plotError=T)
  pdfName = paste0(peaks_pm_750bp$peak_name[i],'_1500bp_window_rhx',".pdf")
  pdf(pdfName, height=1.6, width=1.6)
  print(p)
  dev.off()
}


for(i in c(seq_len(nrow(peaks_pm_500bp)))){
  print(peaks_pm_500bp$plotStart[i])
}


#Plotting mean log2 enrichment ----
#No smoothing of the data - these data were not used to 
#You need: data file, CDSs, PurBox file:

#Data
save(summary_df, file="summary_df.RData")
load("summary_df.RData")

#Enrichment smoothed by 50bp
load("summary_smoothed.RData")

#Enrichment smoothed by 100bp
load("summary_smoothed_100bp.RData")

#Load CDS file (with gene names, etc)
save(CDSs, file="annotated_CDSs.RData")
load("annotated_CDSs.RData")

#Load dataframe with PurBox locations - so they can be mapped if desired
purbox_df = read.csv(file="purbox_locations.csv", colClasses=c("NULL",NA,NA))
save(purbox_df, file="PurBox_locations.RData")
load("PurBox_locations.RData")

#To plot specific samples, change the matrix below with the appropriate sample names
plot_samples = c("wt_acgu", "wt_rhx_00","wt_rhx_10")
plot_samples = c("wt_rhx_00","wt_rhx_10")
plot_samples = c("p0_rhx_00", "p0_rhx_10")
plot_samples = c("wt_rhx_10", "p0_rhx_10")
plot_samples = c("wt_rhx_00", "p0_rhx_00")

#To find a given gene, replace gene below and execute these functions
gene_start = CDSs %>% filter(gene=='yaaD') %>% .$start
gene_end = CDSs %>% filter(gene=='yaaD') %>% .$end

#In helperFunctions.R, I added lines to manually specify colors.
#The palettes for strands and samples are defined here.
strand_colors <- c("+"="gray80","-"="gray80")

sample_colors <- c(wt_rhx_00 = "black",
                   wt_rhx_10 = "firebrick1",
                   p0_rhx_00 = "darkseagreen4",
                   p0_rhx_10 = "darkslateblue",
                   wt_acgu = "darkorchid4")

plotEnrich(CDSs, summary_smoothed_100bp %>% dplyr::filter(group %in% plot_samples), 
           15000, 20000, plotFeatures=T, plotPurBox=F, plotError=F)

plotEnrich(CDSs, summary_smoothed_100bp %>% dplyr::filter(group %in% plot_samples), 
           gene_end-1e3, gene_end+0.75e3, plotFeatures=TRUE, plotPurBox=F, plotError=T)

plotEnrich(CDSs, summary_smoothed_100bp %>% dplyr::filter(group %in% plot_samples), 
           gene_start-1e3, gene_start+2e3, plotFeatures=T, plotPurBox=F, plotError=T)

source('helperFunctions.R')

#Pulling mean and sd of log2 enrichment for PurR ChIP peaks ----
peak_locations = read.csv(file="PurR_chip_peak_locations.csv")

peak_locations = peak_locations %>% mutate("peak_number"=c(1:15))
save(peak_locations, file="peak_locations.RData")
load("peak_locations.RData")

peak_data = summary_df %>% dplyr::group_by(group) %>%
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
  dplyr::group_by(group, peak_number) %>%
  dplyr::filter(mean_log2_enrich==max(mean_log2_enrich))

peak_data_max_linear = peak_data %>%
  dplyr::group_by(group, peak_number) %>%
  dplyr::filter(mean_enrich==max(mean_enrich))

peak_data_max_linear = peak_data_max %>%
  left_join(peak_names, by="peak_number")

save(peak_data_max_linear, file="final_peak_locations.RData")
write.csv(peak_data_max_linear, file='final_peak_locations.csv') 

peak_names = read.csv(file="peak_names.csv")

peak_data_max = peak_data_max %>%
  left_join(peak_names, by="peak_number")

final_peak_data = peak_data_max %>%
  separate(group, c('strain', 'condition', 'time'), sep="_")

save(final_peak_data, file="final_peak_data.RData")
load("final_peak_data.RData")



save(peak_data_max, file="peak_data_max.RData")
load("peak_data_max.RData")

write.csv(peak_data_max, file="log2_enrich_at_peaks.csv")

#Pulling mean and sd of linear enrichment (50bp smoothed) for PurR ChIP peaks ----
load("summary_smoothed_50bp_linear.RData")
peak_locations = read.csv(file="PurR_chip_peak_locations.csv")

peak_locations = peak_locations %>% mutate("peak_number"=c(1:15))
save(peak_locations, file="peak_locations.RData")
load("peak_locations.RData")

peak_data_linear = summary_smoothed_nolog2 %>% dplyr::group_by(group) %>%
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

peak_data_max_linear = peak_data_linear %>%
  dplyr::group_by(group, peak_number) %>%
  dplyr::filter(mean_enrich==max(mean_enrich))

peak_data_max_linear = peak_data_max %>%
  left_join(peak_names, by="peak_number")

save(peak_data_max_linear, file="final_peak_locations.RData")
write.csv(peak_data_max_linear, file='final_peak_locations_linear_enrichment.csv') 



#This section can be ignored - original approach to analyzing initial dataset ----
#Plotting log2 fold change
peak_difference %>% filter(strain=='p0') %>%
  ggplot(aes(x=reorder(peak_name,log2_change), y=log2_change)) +
  geom_bar(stat="identity", position = 'dodge', aes(fill=strain)) +
  scale_fill_manual(values=c('darkslategray','burlywood4')) +
  xlab(NULL) +
  ylab('Log2 fold change after RHX') +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_text(hjust=1, 
                                   angle=90, 
                                   vjust=0.3, 
                                   face='italic',
                                   size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14, face='bold'),
        axis.line.y = element_line(size=1),
        axis.line.x = element_,
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) #+
  geom_hline(yintercept=0, linetype='dashed', size=1) +
  facet_wrap(~strain)

#Processing of the first round of PurR ChIP (Run_2691). This is moot after final processing
# read in the wt 00 min enrichment from the actual data, not from bootstrapping
profile1 = as_tibble(npyLoad('wt_10min_vs_untreated_wt_rep1_rhx_00_chip_actual.npy')) %>%
  dplyr::rename(wt_rhx_00 = V1) %>%
  dplyr::mutate(robust_z = convert_to_robustz(wt_rhx_00))
profile1 %>% head
profile1 %>% nrow

profile1 %>%
  ggplot(aes(x=1:nrow(profile1), y=robust_z)) +
  geom_line(alpha=0.5)

# read in the wt 10 min enrichment from the actual data, not from bootstrapping
profile2 = as_tibble(npyLoad('wt_10min_vs_untreated_wt_rep1_rhx_10_chip_actual.npy')) %>%
  dplyr::rename(wt_rhx_10 = V1) %>%
  dplyr::mutate(robust_z = convert_to_robustz(wt_rhx_10))
profile2 %>% head
profile2 %>% nrow

profile2 %>%
  ggplot(aes(x=1:nrow(profile2), y=robust_z)) +
  geom_line(alpha=0.5)

# read in the wt 10 min minus the wt 00 min enrichment from actual data, not from bootstrapping
changeOverZero_wt = as_tibble(npyLoad('wt_10min_vs_untreated_wt_rep1_rhx_10_chip_Sub_wt_rep1_rhx_00_chip_actual.npy')) %>%
  dplyr::rename(ten_vs_zero = V1) %>%
  dplyr::mutate(robust_z = convert_to_robustz(ten_vs_zero),
                strain = "wt",
                position = 1:nrow(.)*10)

changeOverZero_wt %>%
  ggplot(aes(x=position, y=robust_z)) +
  geom_line(alpha=0.5)

# read in the p0 10 min minus the p0 00 min enrichment from actual data, not from bootstrapping
changeOverZero_p0 = as_tibble(npyLoad('p0_10min_vs_untreated_p0_rep1_rhx_10_chip_Sub_p0_rep1_rhx_00_chip_actual.npy')) %>%
  dplyr::rename(ten_vs_zero = V1) %>%
  dplyr::mutate(robust_z = convert_to_robustz(ten_vs_zero),
                strain = "p0",
                position = 1:nrow(.)*10)

changeOverZero_p0 %>%
  ggplot(aes(x=position, y=robust_z)) +
  geom_line(alpha=0.5)

changeOverZero = rbind(changeOverZero_p0, changeOverZero_wt)

changeOverZero %>%
  ggplot(aes(x=position,y=robust_z,color=strain)) +
  geom_line(alpha=0.3)

changeDueToppGpp = tibble(position = changeOverZero_p0$position, ppGppEffect = changeOverZero_wt$robust_z - changeOverZero_p0$robust_z)
changeDueToppGpp %>%
  ggplot(aes(x=position, y=ppGppEffect)) +
  geom_line(alpha=0.5)

# get bootstrapped estimates
medianSignal = npyLoad('wt_rhx_10_min_bootstrapped_wt_rep1_rhx_00_chip_sort_median.npy')
medianSignalDF = tibble(position=1:length(medianSignal)*10, median=medianSignal)
dygraph(medianSignalDF)

medianDiff = npyLoad('wt_rhx_10_min_bootstrapped_wt_rep1_rhx_00_chip_sort_Sub_wt_rep1_rhx_00_chip_sort_median.npy')
upperCIDiff = npyLoad('wt_rhx_10_min_bootstrapped_wt_rep1_rhx_00_chip_sort_Sub_wt_rep1_rhx_00_chip_sort_maxci.npy')
lowerCIDiff = npyLoad('wt_rhx_10_min_bootstrapped_wt_rep1_rhx_00_chip_sort_Sub_wt_rep1_rhx_00_chip_sort_minci.npy')
medianSignalDF = medianSignalDF %>%
  mutate(diff = medianDiff,
         upperCI = upperCIDiff,
         lowerCI = lowerCIDiff)

dygraph(medianSignalDF %>% dplyr::select(position,upperCI,lowerCI,diff) %>% filter(1:nrow(medianSignalDF) %% 5 == 0), main="log2(PurR enrichment 10 minutes after RHX/untreated)") %>%
  dySeries(c("lowerCI","diff","upperCI"), label="enrichment") %>%
  dyOptions(colors=RColorBrewer::brewer.pal(3,"Set1")) %>%
  dyRangeSelector()

medianGR = GRanges(seqnames='CP020102', 
                   ranges=IRanges(start=medianSignalDF$position-medianSignalDF$position[1]+1, width=diff(medianSignalDF$position)[1]))
medianGR$score = ifelse(is.na(medianDiff), 0, medianDiff)
rtracklayer::export(medianGR, "wt_10_minus_wt_00_PurR.wig")

medianSignalDF %>% dplyr::select(position,diff,upperCI,lowerCI) %>% mutate(Chr="CP020102",position=as.integer(position)) %>%
  write_csv("wt_10_minus_wt_00_PurR.csv")

# get bootstrapped estimates for p0
p0medianDiff = npyLoad('p0_bootstrapped_wt_rep1_rhx_00_chip_sort_Sub_wt_rep1_rhx_00_chip_sort_median.npy')
p0upperCIDiff = npyLoad('p0_bootstrapped_wt_rep1_rhx_00_chip_sort_Sub_wt_rep1_rhx_00_chip_sort_maxci.npy')
p0lowerCIDiff = npyLoad('p0_bootstrapped_wt_rep1_rhx_00_chip_sort_Sub_wt_rep1_rhx_00_chip_sort_minci.npy')
p0medianSignalDF = medianSignalDF %>%
  mutate(diff = p0medianDiff,
         upperCI = p0upperCIDiff,
         lowerCI = p0lowerCIDiff)

dygraph(p0medianSignalDF %>% dplyr::select(position,upperCI,lowerCI,diff) %>% filter(1:nrow(medianSignalDF) %% 5 == 0), main="log2(p0 PurR enrichment 10 minutes after RHX/untreated)") %>%
  dySeries(c("lowerCI","diff","upperCI"), label="enrichment") %>%
  dyOptions(colors=RColorBrewer::brewer.pal(3,"Set1")) %>%
  dyRangeSelector()

p0medianGR = GRanges(seqnames='CP020102', 
                   ranges=IRanges(start=p0medianSignalDF$position-p0medianSignalDF$position[1]+1, width=diff(p0medianSignalDF$position)[1]))
p0medianGR$score = ifelse(is.na(p0medianDiff), 0, p0medianDiff)
rtracklayer::export(p0medianGR, "p0_10_minus_p0_00_PurR.wig")

p0medianSignalDF %>% dplyr::select(position,diff,upperCI,lowerCI) %>% mutate(Chr="CP020102",position=as.integer(position)) %>% 
  write_csv("p0_10_minus_wt_00_PurR.csv")


# calculate cpm
data_path = 'Y:/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data'
data_path = '/mnt/jadelab/lab/current/Users/Brent/ChIP-seq/PurR_ChIP_data'
run_name = 'Run_2691'
run_path = paste(data_path, run_name, sep='/')
run_info_base_name = paste0(run_name, '_jawang.csv')
run_info_path = paste(run_path, run_info_base_name, sep='/')

runInfo = read.csv(run_info_path, header=T, skip=18)

# Here we're just naming each sample's bam file name in a column called "bamFilePath".
#   This is useful later for reading in the actual data in one line of code, rather than
#   in a complicated, sloppy-looking for loop.
run2691Info = runInfo %>%
  filter(!(duplicated(Sample_ID))) %>%
  mutate(bamFilePath = paste0(run_path, '/', 'jawang', '/', Sample_ID, '/', Sample_ID, '_sort_filtered.bam'),
         bamIndexPath = paste0(bamFilePath, '.bai'),
         # bamFileObjects = purrr::map2(bamFilePath, bamIndexPath, BamFile),
         replicate = '1')

# add more run info as replicates come in
sampleInfo = run2691Info #%>%
  # bind_rows(run2519Info)
sampleInfo

# parse CDS information
anno_file = 'Y:/lab/current/NGS/Annotations/Bacillus_subtilis_NCIB_3610_CP020102.1.gbk'
anno_file = '/mnt/jadelab/lab/current/NGS/Annotations/Bacillus_subtilis_NCIB_3610_CP020102.1.gbk'

gb = readGenBank(anno_file)
genesDF = genes(gb)
genesDF = as_tibble(genesDF) %>% dplyr::mutate(midpoint = (start+end)/2,
                                        seqnames = "CP020102")

genes_gr = as(genesDF, 'GRanges')

CDSs = as_tibble(transcripts(gb)) %>% 
  mutate(midpoint=(start+end)/2,
         seqnames = "CP020102")

rtpEnd = CDSs[CDSs$product=="replication termination protein",]$end

CDSs = CDSs %>% mutate(       
  headon = ifelse(start > rtpEnd, 
                  ifelse(strand == '-', 0, 1),
                  ifelse(strand == '+', 0, 1)
  )
)

CDS_gr = as(CDSs, 'GRanges')

# count reads aligning to 10-bp windows
starts = seq(1,4215607,by=10)
ends = seq(10,4215607,by=10)

starts = starts[1:length(ends)]

windows_gr = GRanges(seqnames = "CP020102",
                     ranges = IRanges(start=starts, end=ends),
                     strand = "*")

bams = BamFileList(sampleInfo$bamFilePath, index=sampleInfo$bamIndexPath, yieldSize=25000)
se = summarizeOverlaps(windows_gr, bams, mode="IntersectionNotEmpty", ignore.strand=T, inter.feature=FALSE)

all.counts = assays(se)$counts
# toss out exponential phase data
head(all.counts)
tail(all.counts)
# all.counts = all.counts[,1:28] # exclude the 'exponential phase data'

save(all.counts, CDS_gr, windows_gr, file="all_counts.RData")
load('all_counts.RData')

#all.counts.CPM = applyGetCPM(all.counts)

save(all.counts.CPM, CDS_gr, windows_gr, file="all_counts_CPM.RData")
load('all_counts_CPM.RData')

CPMdf = as.data.frame(all.counts.CPM)
col_names = names(CPMdf)
new_names = gsub('_sort_filtered.bam','',col_names)
names(CPMdf) = new_names
names(CPMdf)

CPMdf$position = starts

CPM_df_long = CPMdf %>% gather(key=Sample_ID, value=cpm, -position)
CPM_df_long %>% head
CPM_df_long = CPM_df_long %>%
  left_join(sampleInfo %>% dplyr::select(Sample_ID, Description))
CPM_df_long %>% head

#What we need to plot CPM: CDS df (CDSs), CPM df (CPM_df_long), and positions

#Load CDS file (with gene names, etc)
save(CDSs, file="annotated_CDSs.RData")
load("annotated_CDSs.RData")

#Load data file with reads/million in 10bp windows
save(CPM_df_long, file="CPM_with_descriptions.RData")
load("CPM_with_descriptions.RData")

#Load dataframe with PurBox locations - so they can be mapped if desired
purbox_df = read.csv(file="purbox_locations.csv", colClasses=c("NULL",NA,NA))

save(purbox_df, file="PurBox_locations.RData")
load("PurBox_locations.RData")

#To plot specific samples, change the matrix below with the appropriate sample names
plot_samples = c("wt_rep1_rhx_00_chip","wt_rep1_rhx_10_chip", "p0_rep1_rhx_00_chip", "p0_rep1_rhx_10_chip")
plot_samples = c("acgu_rep1_chip")
plot_samples = c("wt_rep1_rhx_00_chip","wt_rep1_rhx_10_chip")
plot_samples = c("wt_rep1_rhx_00_chip","p0_rep1_rhx_00_chip")
plot_samples = c("acgu_rep1_chip", "wt_rep1_rhx_10_chip")
plot_samples = c("p0_rep1_rhx_00_chip","p0_rep1_rhx_10_chip", "wt_rep1_rhx_00_chip", "wt_rep1_rhx_10_chip","acgu_rep1_chip")
plot_samples = c("p0_rep1_rhx_00_chip","p0_rep1_rhx_10_chip")

#To find a given gene, replace gene below and execute these functions
gene_start = CDSs %>% filter(gene=='purA') %>% .$start
gene_end = CDSs %>% filter(gene=='purA') %>% .$end

#In helperFunctions.R, I added lines to manually specify colors.
#The palettes for strands and samples are defined here.
strand_colors <- c("+"="gray90","-"="gray90")

#for Jade's Wash U presentation
sample_colors <- c(wt_rep1_rhx_00_chip = "black",
                   wt_rep1_rhx_10_chip = "firebrick1",
                   p0_rep1_rhx_00_chip = "darkseagreen2",
                   p0_rep1_rhx_10_chip = "darkgreen",
                   acgu_rep1_chip = "darkorchid4")

sample_colors <- c(wt_rep1_rhx_00_chip = "orange",
                   wt_rep1_rhx_10_chip = "red",
                   p0_rep1_rhx_00_chip = "darkseagreen2",
                   p0_rep1_rhx_10_chip = "darkgreen",
                   acgu_rep1_chip = "darkorchid4")

#To plot with genomic locations you provide
plotLocusCPM(CDSs, 
             CPM_df_long %>% dplyr::filter(Description %in% plot_samples),
             0.69826e6, 0.698325e6,lineAlpha=1,plotFeatures=TRUE, plotPurBox=FALSE)

#To plot positive strand genes with a given gene name
plotLocusCPM(CDSs, 
             CPM_df_long %>% dplyr::filter(Description %in% plot_samples), 
             gene_start-1e3, gene_start+1.5e3, lineAlpha=1, plotPurBox=TRUE)

#To plot negative strand genes with a given gene name
plotLocusCPM(CDSs, 
             CPM_df_long %>% dplyr::filter(Description %in% plot_samples), 
             gene_end-1.5e3, gene_end+1e3, lineAlpha=1, plotPurBox=TRUE)


plotLocus(CDSs, 
          profile1, 
          1,5e5, lineAlpha=0.8)

plotLocus(CDSs, 
          profile1, 
          gene_start-1e3, gene_start+1e3, lineAlpha=0.8)



#Preparing files for GEO submission ----

#Creating CPM datasets for each type of sample
summary_cpm_p0_rhx_00_chip = summary_cpm_analyzed %>% filter(Description == 'p0_rhx_00_chip')
summary_cpm_p0_rhx_00_input = summary_cpm_analyzed %>% filter(Description == 'p0_rhx_00_input')
summary_cpm_p0_rhx_10_chip = summary_cpm_analyzed %>% filter(Description == 'p0_rhx_10_chip')
summary_cpm_p0_rhx_10_input = summary_cpm_analyzed %>% filter(Description == 'p0_rhx_10_input')
summary_cpm_wt_acgu_chip = summary_cpm_analyzed %>% filter(Description == 'wt_acgu_chip')
summary_cpm_wt_acgu_input = summary_cpm_analyzed %>% filter(Description == 'wt_acgu_input')
summary_cpm_wt_rhx_00_chip = summary_cpm_analyzed %>% filter(Description == 'wt_rhx_00_chip')
summary_cpm_wt_rhx_00_input = summary_cpm_analyzed %>% filter(Description == 'wt_rhx_00_input')
summary_cpm_wt_rhx_10_chip = summary_cpm_analyzed %>% filter(Description == 'wt_rhx_10_chip')
summary_cpm_wt_rhx_10_input = summary_cpm_analyzed %>% filter(Description == 'wt_rhx_10_input')

write.table(summary_cpm_analyzed, 'cpm_all_samples.txt')
write.table(summary_cpm_p0_rhx_00_chip, 'cpm_p0_rhx_00_chip.txt')
write.table(summary_cpm_p0_rhx_00_input, 'cpm_p0_rhx_00_input.txt')
write.table(summary_cpm_p0_rhx_10_chip, 'cpm_p0_rhx_10_chip.txt')
write.table(summary_cpm_p0_rhx_10_input, 'cpm_p0_rhx_10_input.txt')
write.table(summary_cpm_wt_acgu_chip, 'cpm_wt_acgu_chip.txt')
write.table(summary_cpm_wt_acgu_input, 'cpm_wt_acgu_input.txt')
write.table(summary_cpm_wt_rhx_00_chip, 'cpm_wt_rhx_00_chip.txt')
write.table(summary_cpm_wt_rhx_00_input, 'cpm_wt_rhx_00_input.txt')
write.table(summary_cpm_wt_rhx_10_chip, 'cpm_wt_rhx_10_chip.txt')
write.table(summary_cpm_wt_rhx_10_input, 'cpm_wt_rhx_10_input.txt')

#Converting enrichment datasets to .txt
summary_log2_enrich_p0_rhx_00 = summary_df %>% filter(group == 'p0_rhx_00')
summary_log2_enrich_p0_rhx_10 = summary_df %>% filter(group == 'p0_rhx_10')
summary_log2_enrich_wt_acgu = summary_df %>% filter(group == 'wt_acgu')
summary_log2_enrich_wt_rhx_00 = summary_df %>% filter(group == 'wt_rhx_00')
summary_log2_enrich_wt_rhx_10 = summary_df %>% filter(group == 'wt_rhx_10')

write.table(summary_df, 'log2_enrichment_all_samples.txt')
write.table(summary_log2_enrich_p0_rhx_00, 'log2_enrichment_p0_rhx_00.txt')
write.table(summary_log2_enrich_p0_rhx_10, 'log2_enrichment_p0_rhx_10.txt')
write.table(summary_log2_enrich_wt_acgu, 'log2_enrichment_wt_acgu.txt')
write.table(summary_log2_enrich_wt_rhx_00, 'log2_enrichment_wt_rhx_00.txt')
write.table(summary_log2_enrich_wt_rhx_10, 'log2_enrichment_wt_rhx_10.txt')

#Converting 50bp smoothed datasets to .txt
summary_log2_enrich_smoothed_p0_rhx_00 = summary_smoothed %>% filter(group == 'p0_rhx_00')
summary_log2_enrich_smoothed_p0_rhx_10 = summary_smoothed %>% filter(group == 'p0_rhx_10')
summary_log2_enrich_smoothed_wt_acgu = summary_smoothed %>% filter(group == 'wt_acgu')
summary_log2_enrich_smoothed_wt_rhx_00 = summary_smoothed %>% filter(group == 'wt_rhx_00')
summary_log2_enrich_smoothed_wt_rhx_10 = summary_smoothed %>% filter(group == 'wt_rhx_10')

write.table(summary_smoothed, 'log2_enrichment_smoothed_all_samples.txt')
write.table(summary_log2_enrich_smoothed_p0_rhx_00, 'log2_enrichment_smoothed_p0_rhx_00.txt')
write.table(summary_log2_enrich_smoothed_p0_rhx_10, 'log2_enrichment_smoothed_p0_rhx_10.txt')
write.table(summary_log2_enrich_smoothed_wt_acgu, 'log2_enrichment_smoothed_wt_acgu.txt')
write.table(summary_log2_enrich_smoothed_wt_rhx_00, 'log2_enrichment_smoothed_wt_rhx_00.txt')
write.table(summary_log2_enrich_smoothed_wt_rhx_10, 'log2_enrichment_smoothed_wt_rhx_10.txt')

#Converting 50bp smoothed linear datasets to .txt
summary_linear_enrich_smoothed_p0_rhx_00 = summary_smoothed_nolog2 %>% filter(group == 'p0_rhx_00')
summary_linear_enrich_smoothed_p0_rhx_10 = summary_smoothed_nolog2 %>% filter(group == 'p0_rhx_10')
summary_linear_enrich_smoothed_wt_acgu = summary_smoothed_nolog2 %>% filter(group == 'wt_acgu')
summary_linear_enrich_smoothed_wt_rhx_00 = summary_smoothed_nolog2 %>% filter(group == 'wt_rhx_00')
summary_linear_enrich_smoothed_wt_rhx_10 = summary_smoothed_nolog2 %>% filter(group == 'wt_rhx_10')

write.table(summary_smoothed_nolog2, 'linear_enrichment_smoothed_all_samples.txt')
write.table(summary_linear_enrich_smoothed_p0_rhx_00, 'linear_enrichment_smoothed_p0_rhx_00.txt')
write.table(summary_linear_enrich_smoothed_p0_rhx_10, 'linear_enrichment_smoothed_p0_rhx_10.txt')
write.table(summary_linear_enrich_smoothed_wt_acgu, 'linear_enrichment_smoothed_wt_acgu.txt')
write.table(summary_linear_enrich_smoothed_wt_rhx_00, 'linear_enrichment_smoothed_wt_rhx_00.txt')
write.table(summary_linear_enrich_smoothed_wt_rhx_10, 'linear_enrichment_smoothed_wt_rhx_10.txt')

#Renaming .RData files
save(summary_cpm_analyzed, file='summary_cpm_all_samples.RData')
save(summary_df, file='summary_log2_enrich_all_samples.RData')
save(summary_smoothed, file='summary_log2_enrichment_smoothed_all_samples.RData')
save(summary_smoothed_nolog2, file='summary_linear_enrichment_all_samples.RData')

save(summary_cpm_p0_rhx_00_chip, file='summary_cpm_p0_rhx_00_chip.RData')
save(summary_cpm_p0_rhx_00_input, file='summary_cpm_p0_rhx_00_input.RData')
save(summary_cpm_p0_rhx_10_chip, file='summary_cpm_p0_rhx_10_chip.RData')
save(summary_cpm_p0_rhx_10_input, file='summary_cpm_p0_rhx_10_input.RData')
save(summary_cpm_wt_acgu_chip, file='summary_cpm_wt_acgu_chip.RData')
save(summary_cpm_wt_acgu_input, file='summary_cpm_wt_acgu_input.RData')
save(summary_cpm_wt_rhx_00_chip, file='summary_cpm_wt_rhx_00_chip.RData')
save(summary_cpm_wt_rhx_00_input, file='summary_cpm_wt_rhx_00_input.RData')
save(summary_cpm_wt_rhx_10_chip, file='summary_cpm_wt_rhx_10_chip.RData')
save(summary_cpm_wt_rhx_10_input, file='summary_cpm_wt_rhx_10_input.RData')

save(summary_log2_enrich_p0_rhx_00, file='summary_log2_enrich_p0_rhx_00.RData')
save(summary_log2_enrich_p0_rhx_10, file='summary_log2_enrich_p0_rhx_10.RData')
save(summary_log2_enrich_wt_acgu, file='summary_log2_enrich_wt_acgu.RData')
save(summary_log2_enrich_wt_rhx_00, file='summary_log2_enrich_wt_rhx_00.RData')
save(summary_log2_enrich_wt_rhx_10, file='summary_log2_enrich_wt_rhx_10.RData')

save(summary_log2_enrich_smoothed_p0_rhx_00, file='summary_log2_enrich_smoothed_p0_rhx_00.RData')
save(summary_log2_enrich_smoothed_p0_rhx_10, file='summary_log2_enrich_smoothed_p0_rhx_10.RData')
save(summary_log2_enrich_smoothed_wt_acgu, file='summary_log2_enrich_smoothed_wt_acgu.RData')
save(summary_log2_enrich_smoothed_wt_rhx_00, file='summary_log2_enrich_smoothed_wt_rhx_00.RData')
save(summary_log2_enrich_smoothed_wt_rhx_10, file='summary_log2_enrich_smoothed_wt_rhx_10.RData')

save(summary_linear_enrich_smoothed_p0_rhx_00, file='summary_linear_enrich_smoothed_p0_rhx_00.RData')
save(summary_linear_enrich_smoothed_p0_rhx_10, file='summary_linear_enrich_smoothed_p0_rhx_10.RData')
save(summary_linear_enrich_smoothed_wt_acgu, file='summary_linear_enrich_smoothed_wt_acgu.RData')
save(summary_linear_enrich_smoothed_wt_rhx_00, file='summary_linear_enrich_smoothed_wt_rhx_00.RData')
save(summary_linear_enrich_smoothed_wt_rhx_10, file='summary_linear_enrich_smoothed_wt_rhx_10.RData')
