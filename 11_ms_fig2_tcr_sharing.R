## visualization -----------------------------------------------------------
# data resources
# The notebook is 11_ms_figs.ipynb Section 1.4.1
# The data is from metadata of the anndata below
# adata_t.obs.to_csv(wd + 't_cells_metadata.csv')
# adata_gdt.obs.to_csv(wd + 'gdt_cells_metadata.csv')
# adata_b.obs.to_csv(wd + 'b_cells_metadata.csv')

# import packges

library(dplyr)
library(purrr)
#library(DESeq2)
library(Azimuth)

library(cowplot) # ggsave2
library(ggrastr) # geom_point_rast
library(ggplot2)
library(ggalluvial)
library(ggrepel)
library(patchwork)
library(Matrix)

#install.packages("viridis")
library(viridis)

library(stringr) # string manipulation
library(magrittr) # %>% 

setwd("/path/to/your/wd/")
wd_path <- getwd()


donor_colors = c("#990000", "#333366", "#336699", "#e296ad", "#009999", 
                 "#DAA520", "#95a674", "#a98a7b", "#6495ED", "#FFA07A")

t_colors = c('#240771', '#9920BE', '#63557C', '#7c6b9c', '#a89dbd', '#d3cede', # cd4
             '#70441C', '#A3764F', '#f1c27d',  # cd8 Tnv/cm # cd8 Tmem KLRC2 # cd8 Trm
             '#ffdbac', '#c54349', '#e2a1a4', '#f1d0d2',# cd8 Tmem IFN # Tmem GZMK # Tmem_LAYN # Tmem  MKI67
             '#36802d', '#77ab59', '#c9df8a', # cd8 Tmem TIGIT # Temra # NKT_like
             '#007bdb', #gdT
             '#66afe9' #MAIT
)

gdt_colors = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61')

### abTCR -------------------------------------------------------------------
abtcr = read.csv('t_cells_metadata.csv', row.names = 1)

abtcr_clean = abtcr[!abtcr$barcodes %>% is_in(gdtcr$barcodes), ]

alluvium_theme <- theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
                        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                        #panel.grid.major = element_blank(), #remove major gridlines
                        panel.grid.minor = element_blank(),
                        axis.line.y.left = element_line(),
                        axis.line.x.bottom = element_line(),
                        axis.title = element_text(size = 8, hjust = 0.5),
                        axis.text = element_text(size = 8))

#TCR alpha sharing group by phenotypes
# TCR sharing strata is cdr3_anno 
TCRA_anno <- abtcr_clean %>% 
  filter(!str_detect(cdr3_aa_TRA, pattern = '\\?')) %>% 
  filter(func_TRA == "functional")%>% 
  filter(!is.na(v_TRA) & cdr3_freq_TRA > 1) %>%
  group_by(cdr3_aa_TRA, anno_t_8, timepoint) %>%  
  summarise(cdr3_freq_TRA = n(), .groups = "keep") %>% 
  arrange(cdr3_freq_TRA, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRA = factor(cdr3_aa_TRA, unique(cdr3_aa_TRA)))

TCRA_anno

TCRA_anno$anno_strata <- paste0(TCRA_anno$anno_t_8, "_", TCRA_anno$cdr3_aa_TRA)

#is_lodes_form(TCRA, key = timepoint, value = cdr3_aa_TRA, id = donor_strata) # TRUE

TCRA_anno$anno_t_8 <- factor(TCRA_anno$anno_t_8, levels = c(    
  'CD4 Tnv/cm',
  'Th1',
  'Th2',
  'Th17',
  'Tfh',
  'Treg',
  'CD8 Tnv/cm',
  #'CD8 Tmem KLRC2',
  'CD8 Trm',
  'CD8 Tmem IFN',
  'CD8 Tmem GZMK',
  'CD8 Tmem LAYN',
  'CD8 Tmem MKI67',
  'CD8 Tmem TIGIT',
  'CD8 Temra',
  'CD8 NKT-like',
  #'gdT Vd2',
  'MAIT'))

# Create a unique list of anno for custom coloring
unique_anno <- unique(TCRA_anno$anno_t_8)

# Define custom colors using Viridis palette
graft_colors <- c('#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b','#e377c2', '#b5bd61', '#17becf', '#aec7e8')

d30_colors <- c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b','#e377c2', '#b5bd61', '#17becf', '#aec7e8')

d100_colors <- c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#17becf')

custom_colors = c(rev(graft_colors), rev(d30_colors), rev(d100_colors))

(ggplot(TCRA,
        aes(x = timepoint, stratum = donor, alluvium = donor_strata,
            fill = cdr3_aa_TRA, label = donor)) +
    geom_alluvium(show.legend = F) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "grey", show.legend = F) +
    geom_stratum(aes(fill = donor), show.legend = F, 
                 alpha = 0, fill = custom_colors)+
    geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    NoLegend()+ xlab("") + ylab("frequence")+ alluvium_theme )%T>% 
  figsave("/tra_sharing_group_donor.pdf", w = 100, h = 100)

# TCRB
TCRB_donor <- abtcr_clean %>% 
  filter(!str_detect(cdr3_aa_TRB, pattern = '\\?')) %>% 
  filter(func_TRB == "functional")%>% 
  filter(!is.na(v_TRB) & cdr3_freq_TRB > 1) %>%
  group_by(cdr3_aa_TRB, donor, timepoint) %>%  
  summarise(cdr3_freq_TRB = n(), .groups = "keep") %>% 
  arrange(cdr3_freq_TRB, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRB = factor(cdr3_aa_TRB, unique(cdr3_aa_TRB)))

TCRB_donor

TCRB_donor$donor_strata <- paste0(TCRB_donor$donor, "_", TCRB_donor$cdr3_aa_TRB)

TCRB_donor$timepoint = factor(TCRB_donor$timepoint, levels = c('Graft', 'D30', 'D100'))

trb_donor_graft = c(
  #"#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  "#a98a7b", #MRD012
  "#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

trb_donor_d30 = c(
  "#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  "#a98a7b", #MRD012
  "#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

trb_donor_d100 = c(
  "#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  #"#DAA520", #MRD009
  #"#95a674", #MRD010
  #"#a98a7b", #MRD012
  "#6495ED"#, #MRD013
  #"#FFA07A"  #MRD016 
)

(ggplot(data = TCRB_donor, aes(x = timepoint, y = cdr3_freq_TRB, 
                               alluvium = donor_strata, 
                               stratum = donor, 
                               fill = donor_strata,
                               label = donor)) +
    alluvium_theme + 
    #geom_alluvium(show.legend = TRUE) + 
    geom_flow(stat = "alluvium", color = "black", #show.legend = TRUE, 
              alpha = 1, linewidth = 0.1, fill = 'grey') +
    geom_stratum(aes(stratum = donor), width = 1/3, show.legend = TRUE, 
                 fill = c(rev(trb_donor_graft), rev(trb_donor_d30), rev(trb_donor_d100)),
                 #color = donor_colors, 
                 size = 0.1, 
                 alpha = 0.6) + 
    #scale_fill_viridis_d(direction = -1, option = 'D', begin = 0.1, end = 0.9) +
    #NoLegend()+ 
    #theme(legend.position = "right")+
    ggtitle("TCR Beta sharing between timepoints") +
    theme(plot.title = element_text(size = 8, hjust = 0.5)) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "right") +
    #geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    xlab("") + ylab("Frequence")) %T>%
  figsave("/trb_sharing_timepoint_group_donor.pdf", w = 60, h = 50)


# freq > 0
# TCRB
TCRB_tp <- abtcr_clean %>% 
  filter(!str_detect(cdr3_aa_TRB, pattern = '\\?')) %>% 
  filter(func_TRB == "functional")%>% 
  filter(!is.na(v_TRB) & cdr3_freq_TRB > 0) %>%
  group_by(cdr3_aa_TRB, donor, timepoint) %>%  
  summarise(cdr3_freq_TRB = n(), .groups = "keep") %>% 
  arrange(cdr3_freq_TRB, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRB = factor(cdr3_aa_TRB, unique(cdr3_aa_TRB)))

TCRB_tp

TCRB_tp$donor_strata <- paste0(TCRB_tp$donor, "_", TCRB_tp$cdr3_aa_TRB)

TCRB_tp$timepoint = factor(TCRB_tp$timepoint, levels = c('Graft', 'D30', 'D100'))

trb_donor_graft = c(
  #"#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  "#a98a7b", #MRD012
  "#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

trb_donor_d30 = c(
  "#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  "#a98a7b", #MRD012
  "#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

trb_donor_d100 = c(
  "#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  #"#DAA520", #MRD009
  #"#95a674", #MRD010
  #"#a98a7b", #MRD012
  "#6495ED"#, #MRD013
  #"#FFA07A"  #MRD016 
)

(ggplot(data = TCRB_tp, aes(x = timepoint, y = cdr3_freq_TRB, 
                            alluvium = donor_strata, 
                            stratum = donor, 
                            fill = donor_strata,
                            label = donor)) +
    alluvium_theme + 
    #geom_alluvium(show.legend = TRUE) + 
    geom_flow(stat = "alluvium", color = "black", #show.legend = TRUE, 
              alpha = 1, linewidth = 0.1, fill = 'grey') +
    geom_stratum(aes(stratum = donor), width = 1/3, show.legend = TRUE, 
                 fill = c(rev(trb_donor_graft), rev(trb_donor_d30), rev(trb_donor_d100)),
                 #color = donor_colors, 
                 size = 0.1, 
                 alpha = 0.6) + 
    #scale_fill_viridis_d(direction = -1, option = 'D', begin = 0.1, end = 0.9) +
    #NoLegend()+ 
    #theme(legend.position = "right")+
    ggtitle("TCR Beta sharing between timepoints") +
    theme(plot.title = element_text(size = 8, hjust = 0.5)) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "right") +
    #geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    xlab("Timepoints") + ylab("Counts")) %T>%
  figsave("/trb_sharing_timepoint_group_donor_freq_over_0.pdf", w = 50, h = 50)

# Extract persiting clones
result <- TCRB_tp %>%
  group_by(cdr3_aa_TRB) %>%
  summarize(unique_timepoints = n_distinct(timepoint)) %>%
  filter(unique_timepoints >= 2) %>%
  select(cdr3_aa_TRB)

filtered_data <- TCRB_tp %>%
  inner_join(result, by = "cdr3_aa_TRB")

filtered_data %>% arrange(-cdr3_freq_TRB)

write.csv(result, file = '/persisting_tcrb_cdr3s.csv')


### gdTCR -------------------------------------------------------------------
gdtcr = read.csv('gdt_cells_metadata.csv', row.names = 1)

#TCR gamma sharing group by donor
# TCR sharing strata is cdr3_anno 
TCRD_tp <- gdtcr %>% 
  #filter(anno_gdt_1 %in% c('Vd1 Tnv', 'Vd1 T GZMK', 'Vd1 CTL')) %>% 
  filter(!str_detect(cdr3_aa_TRD, pattern = '\\?')) %>% 
  #filter(func_TRD == "functional")%>% 
  filter(!is.na(v_TRD) & cdr3_freq_TRD > 1) %>%
  group_by(cdr3_aa_TRD, donor, timepoint) %>%  
  summarise(cdr3_freq_TRD = n(), .groups = "keep") %>% 
  arrange(cdr3_freq_TRD, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRD = factor(cdr3_aa_TRD, unique(cdr3_aa_TRD)))

TCRD_tp$timepoint = factor(TCRD_tp$timepoint, levels = c('Graft', 'D30', 'D100'))

TCRD_tp$to_strata <- paste0(TCRD_tp$donor, "_", TCRD_tp$cdr3_aa_TRD)

trd_donor_graft = c(
  #"#990000", #MRD002 
  #"#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  #"#a98a7b", #MRD012
  #"#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

trd_donor_d30 = c(
  #"#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  #"#95a674", #MRD010
  #"#a98a7b", #MRD012
  #"#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

trd_donor_d100 = c(
  #"#990000", #MRD002 
  #"#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  #"#DAA520", #MRD009
  #"#95a674", #MRD010
  #"#a98a7b", #MRD012
  "#6495ED"#, #MRD013
  #"#FFA07A"  #MRD016 
)


(ggplot(data = TCRD_tp, aes(x = timepoint, y = cdr3_freq_TRD, 
                            alluvium = to_strata, 
                            stratum = donor, 
                            fill = to_strata,
                            label = donor)) +
    alluvium_theme + 
    #geom_alluvium(show.legend = TRUE) + 
    geom_flow(stat = "alluvium", color = "black", #show.legend = TRUE, 
              alpha = 1, linewidth = 0.1, fill = 'grey') +
    geom_stratum(aes(stratum = donor), width = 1/3, show.legend = TRUE, 
                 fill = c(rev(trd_donor_graft), rev(trd_donor_d30), rev(trd_donor_d100)),
                 #color = donor_colors, 
                 size = 0.1, 
                 alpha = 0.6) + 
    #scale_fill_viridis_d(direction = -1, option = 'D', begin = 0.1, end = 0.9) +
    #NoLegend()+ 
    #theme(legend.position = "right")+
    ggtitle("TCR Delta sharing between timepoints") +
    theme(plot.title = element_text(size = 8, hjust = 0.5)) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "right") +
    #geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    xlab("") + ylab("Frequence")) %T>%
  figsave("abtcr_ms/trd_sharing_timepoint_group_donor.pdf", w = 60, h = 50)

# freq > 0
TCRD_tp <- gdtcr %>% 
  #filter(anno_gdt_1 %in% c('Vd1 Tnv', 'Vd1 T GZMK', 'Vd1 CTL')) %>% 
  filter(!str_detect(cdr3_aa_TRD, pattern = '\\?')) %>% 
  #filter(func_TRD == "functional")%>% 
  filter(!is.na(v_TRD) & cdr3_freq_TRD > 0) %>%
  group_by(cdr3_aa_TRD, donor, timepoint) %>%  
  summarise(cdr3_freq_TRD = n(), .groups = "keep") %>% 
  arrange(cdr3_freq_TRD, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRD = factor(cdr3_aa_TRD, unique(cdr3_aa_TRD)))

TCRD_tp$timepoint = factor(TCRD_tp$timepoint, levels = c('Graft', 'D30', 'D100'))

TCRD_tp$to_strata <- paste0(TCRD_tp$donor, "_", TCRD_tp$cdr3_aa_TRD)

trd_donor_graft = c(
  #"#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  #"#a98a7b", #MRD012
  "#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

trd_donor_d30 = c(
  #"#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  "#a98a7b", #MRD012
  #"#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

trd_donor_d100 = c(
  "#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  #"#DAA520", #MRD009
  #"#95a674", #MRD010
  #"#a98a7b", #MRD012
  "#6495ED"#, #MRD013
  #"#FFA07A"  #MRD016 
)


(ggplot(data = TCRD_tp, aes(x = timepoint, y = cdr3_freq_TRD, 
                            alluvium = to_strata, 
                            stratum = donor, 
                            fill = to_strata,
                            label = donor)) +
    alluvium_theme + 
    #geom_alluvium(show.legend = TRUE) + 
    geom_flow(stat = "alluvium", color = "black", #show.legend = TRUE, 
              alpha = 1, linewidth = 0.1, fill = 'grey') +
    geom_stratum(aes(stratum = donor), width = 1/3, show.legend = TRUE, 
                 fill = c(rev(trd_donor_graft), rev(trd_donor_d30), rev(trd_donor_d100)),
                 #color = donor_colors, 
                 size = 0.1, 
                 alpha = 0.6) + 
    #scale_fill_viridis_d(direction = -1, option = 'D', begin = 0.1, end = 0.9) +
    #NoLegend()+ 
    #theme(legend.position = "right")+
    ggtitle("TCR Delta sharing between timepoints") +
    theme(plot.title = element_text(size = 8, hjust = 0.5)) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "right") +
    #geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    xlab("Timepoints") + ylab("Counts")) %T>%
  figsave("abtcr_ms/trd_sharing_timepoint_group_donor_freq_over_0.pdf", w = 50, h = 50)


### BCR ---------------------------------------------------------------------
bcr_final = read.csv('b_cells_metadata.csv', row.names = 1)

# TCR sharing strata is cdr3_anno 
BCR_tp <- bcr_final %>% 
  #filter(anno_gdt_1 %in% c('Vd1 Tnv', 'Vd1 T GZMK', 'Vd1 CTL')) %>% 
  filter(!str_detect(cdr3_aa_IGH, pattern = '\\?')) %>% 
  #filter(func_TRD == "functional")%>% 
  filter(!is.na(v_IGH) & cdr3_freq_IGH > 0) %>%
  group_by(cdr3_aa_IGH, donor, timepoint) %>%  
  summarise(cdr3_freq_IGH = n(), .groups = "keep") %>% 
  arrange(cdr3_freq_IGH, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_IGH = factor(cdr3_aa_IGH, unique(cdr3_aa_IGH)))

BCR_tp$timepoint = factor(BCR_tp$timepoint, levels = c('Graft', 'D30', 'D100'))

BCR_tp$to_strata <- paste0(BCR_tp$donor, "_", BCR_tp$cdr3_aa_IGH)

bcr_donor_graft = c(
  #"#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  "#a98a7b", #MRD012
  "#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

bcr_donor_d30 = c(
  #"#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  "#DAA520", #MRD009
  "#95a674", #MRD010
  "#a98a7b", #MRD012
  "#6495ED", #MRD013
  "#FFA07A"  #MRD016 
)

bcr_donor_d100 = c(
  "#990000", #MRD002 
  "#333366", #MRD003
  "#336699", #MRD004
  "#e296ad", #MRD007
  "#009999", #MRD008
  #"#DAA520", #MRD009
  #"#95a674", #MRD010
  #"#a98a7b", #MRD012
  "#6495ED"#, #MRD013
  #"#FFA07A"  #MRD016 
)


(ggplot(data = BCR_tp, aes(x = timepoint, y = cdr3_freq_IGH, 
                           alluvium = to_strata, 
                           stratum = donor, 
                           fill = to_strata,
                           label = donor)) +
    alluvium_theme + 
    #geom_alluvium(show.legend = TRUE) + 
    geom_flow(stat = "alluvium", color = "black", #show.legend = TRUE, 
              alpha = 1, linewidth = 0.1, fill = 'grey') +
    geom_stratum(aes(stratum = donor), width = 1/3, show.legend = TRUE, 
                 fill = c(rev(bcr_donor_graft), rev(bcr_donor_d30), rev(bcr_donor_d100)),
                 #color = donor_colors, 
                 size = 0.1, 
                 alpha = 0.6) + 
    #scale_fill_viridis_d(direction = -1, option = 'D', begin = 0.1, end = 0.9) +
    #NoLegend()+ 
    #theme(legend.position = "right")+
    ggtitle("IGH sharing between timepoints") +
    theme(plot.title = element_text(size = 8, hjust = 0.5)) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "right") +
    #geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    xlab("Timepoints") + ylab("Counts")) %T>%
  figsave("/igh_sharing_timepoint_group_donor_over_0.pdf", w = 50, h = 50)


# TCR sharing strata is cdr3_anno 
BCR_tp <- bcr_final %>% 
  #filter(anno_gdt_1 %in% c('Vd1 Tnv', 'Vd1 T GZMK', 'Vd1 CTL')) %>% 
  filter(!str_detect(cdr3_aa_IGH, pattern = '\\?')) %>% 
  #filter(func_TRD == "functional")%>% 
  filter(!is.na(v_IGH) & cdr3_freq_IGH > 1) %>%
  group_by(cdr3_aa_IGH, donor, timepoint) %>%  
  summarise(cdr3_freq_IGH = n(), .groups = "keep") %>% 
  arrange(cdr3_freq_IGH, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_IGH = factor(cdr3_aa_IGH, unique(cdr3_aa_IGH)))

BCR_tp$timepoint = factor(BCR_tp$timepoint, levels = c('Graft', 'D30', 'D100'))

BCR_tp$to_strata <- paste0(BCR_tp$donor, "_", BCR_tp$cdr3_aa_IGH)

bcr_donor_graft = c(#'#1f77b4', #MRD002
  '#ff7f0e', #MRD003
  '#279e68', #MRD004
  '#d62728', #MRD007
  '#aa40fc', #MRD008
  '#8c564b', #MRD009
  '#e377c2', #MRD010
  '#b5bd61', #MRD012
  '#17becf', #MRD013
  '#aec7e8'  #MRD016
)

bcr_donor_d30 = c(#'#1f77b4', #MRD002
  '#ff7f0e', #MRD003
  '#279e68', #MRD004
  '#d62728', #MRD007
  '#aa40fc', #MRD008
  '#8c564b', #MRD009
  '#e377c2', #MRD010
  '#b5bd61', #MRD012
  '#17becf', #MRD013
  '#aec7e8'  #MRD016
)

bcr_donor_d100 = c('#1f77b4', #MRD002
                   '#ff7f0e', #MRD003
                   '#279e68', #MRD004
                   '#d62728', #MRD007
                   '#aa40fc', #MRD008
                   #'#8c564b', #MRD009
                   #'#e377c2', #MRD010
                   #'#b5bd61', #MRD012
                   '#17becf' #MRD013
                   #'#aec7e8'  #MRD016
)


(ggplot(data = BCR_tp, aes(x = timepoint, y = cdr3_freq_IGH, 
                           alluvium = to_strata, 
                           stratum = donor, 
                           fill = to_strata,
                           label = donor)) +
    alluvium_theme + 
    #geom_alluvium(show.legend = TRUE) + 
    geom_flow(stat = "alluvium", color = "black", #show.legend = TRUE, 
              alpha = 1, linewidth = 0.1, fill = 'grey') +
    geom_stratum(aes(stratum = donor), width = 1/3, show.legend = TRUE, 
                 fill = c(rev(bcr_donor_graft), rev(bcr_donor_d30), rev(bcr_donor_d100)),
                 #color = donor_colors, 
                 size = 0.1, 
                 alpha = 0.6) + 
    #scale_fill_viridis_d(direction = -1, option = 'D', begin = 0.1, end = 0.9) +
    #NoLegend()+ 
    #theme(legend.position = "right")+
    ggtitle("IGH sharing between timepoints") +
    theme(plot.title = element_text(size = 8, hjust = 0.5)) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), 
          legend.position = "right") +
    #geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    xlab("") + ylab("Frequence")) %T>%
  figsave("/igh_sharing_timepoint_group_donor_over_0.pdf", w = 60, h = 50)
