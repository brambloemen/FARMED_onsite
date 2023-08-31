# libraries --------------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales)
library(RColorBrewer)
# functions --------------------------------------------------------------------

cat_tsv <- function(dir, pattern){
    files <- list.files(dir, pattern)
    df <- data.frame()
    for (file in files){
        fp_complete <- paste0(dir, file)
        newfile <- fread(fp_complete, sep="\t", integer64 = "numeric")
        newfile <- newfile %>%
            mutate(Experiment = str_extract(file, "[A-Za-z_].+$"),
                   Experiment = str_remove(Experiment, pattern))
        df <- rbind(df, newfile)
    }
    return(df)
}

rarefaction_byExpTime <- function(df, ...){
  df <- df %>%
    arrange(Experiment, Time) %>%
    # take only the first occurence of a species per experiment
    group_by(Experiment, ...) %>% slice_head(n = 1) %>% 
    group_by(Experiment, Time) %>%
    reframe(n_unique = n(),
              .groups = "drop") %>%
    group_by(Experiment) %>%
    reframe(Time = Time,
              n_unique = cumsum(n_unique),
              .groups = "drop") %>%
    complete(Experiment, Time) %>%
    fill(n_unique, .direction = "down")
  return(df)
}


lineplot <- function(df, yvar){
  plot <- df %>% ggplot(aes(x = Time, y = {{yvar}})) +
                  geom_line(aes(color=Experiment)) + ylim(0, NA)
  return(plot)
}

# Import data and generate plots --------------------------------------------------------------------

abfhpv <- cat_tsv("../Plots/Summarized_data/CovByTime/", "_abfhpvByTime.tsv") %>%
  mutate(Organism_QID = n_match_bases/Total_readlength) %>%
  filter(!str_detect(Species, "([P|p]lasmid)|([V|v]irus)|([P|p]hage)")) %>%
  mutate(Experiment = case_when(
           Experiment == "Fecal_background" ~ "UEQL",
           Experiment == "Fecal_GMSspikeI_CB_beads_bc2" ~ "BQR",
           Experiment == "Fecal_GMSspikeI_CB_beads_VSK" ~ "BQV",
           Experiment == "Fecal_GMSspikeI_CB_RAD" ~ "BDR",
           Experiment == "Fecal_GMSspikeI_QuickDNA_LSK" ~ "EQL",
           Experiment == "Fecal_GMSspikeI_QuickDNA_RAD" ~ "EQR",
           .default = "DMC"
         ))

abfhpv_SpeciesByTime <- abfhpv %>% filter(Organism_QID >= 0.80) %>% rarefaction_byExpTime(Species)
abfhpv_SpeciesByTime <- abfhpv %>% filter(Organism_QID >= 0.80, Coverage >= 5) %>% rarefaction_byExpTime(Species)

# reference community data
reference_data <- read.csv("/scratch/brbloemen/FARMED/WP3_T3_onsite_test/FARMED_onsite/GMSspikeI_AMRprofile.csv", sep = ";") %>%
    mutate(Reference_AMR_link=TRUE,
           Species = Organism) %>%
    select(-Organism)
reference_ARGs <- data.frame(reference_data$ARG) %>% mutate(Reference_AMR=TRUE)

reference_comm <- read.csv("/scratch/brbloemen/FARMED/WP3_T3_onsite_test/FARMED_onsite/zymo_GMS.csv", sep = ";", dec=",")
# resistance classes
resistance_classes <- read.csv("/scratch/brbloemen/FARMED/WP3_T3_onsite_test/FARMED_onsite/resistance_classes.tsv", sep="\t") %>%
  mutate(Resistance = str_remove(Resistance, "\\sresistance")) %>% select(Gene, Resistance)


# AMR linking summaries
GMSdb <- cat_tsv("../Plots/Summarized_data/AMRlinking/", "_GMSvResF.tsv") %>%
  mutate(Time = Time_hours) %>%
  mutate(Experiment = case_when(
           Experiment == "Fecal_background" ~ "UEQL",
           Experiment == "Fecal_GMSspikeI_CB_beads_bc2" ~ "BQR",
           Experiment == "Fecal_GMSspikeI_CB_beads_VSK" ~ "BQV",
           Experiment == "Fecal_GMSspikeI_CB_RAD" ~ "BDR",
           Experiment == "Fecal_GMSspikeI_QuickDNA_LSK" ~ "EQL",
           Experiment == "Fecal_GMSspikeI_QuickDNA_RAD" ~ "EQR",
           .default = "DMC"
         ))
GMSdb$Experiment <- ordered(GMSdb$Experiment, levels=c("UEQL", "DMC", "EQL", "EQR", "BDR", "BQR", "BQV"))

# Unique ARGs by time per experiment
GMSdb_unique_ARGs <- GMSdb %>% rarefaction_byExpTime(ARG)

# absolute ARG counts per experiment over time
ARG_byExpTime <- GMSdb %>%
  arrange(Experiment, Time) %>%
    group_by(Experiment) %>%
    reframe(n_reads = cumsum(n_reads), Time = Time) %>%
    group_by(Experiment, Time) %>%
    slice_tail(n=1) %>% ungroup() %>%
    complete(Experiment, Time) %>%
    fill(n_reads, .direction = "down")

combined <- merge(GMSdb_unique_ARGs, abfhpv_SpeciesByTime, by=c("Time", "Experiment"), all=TRUE) %>%
  rename(n_uniq_ARGS = n_unique.x, n_uniq_Spec = n_unique.y)
combined <- merge(combined, ARG_byExpTime, by=c("Time", "Experiment"), all=TRUE)
combined <- combined %>%
    filter(Experiment != "DMC") %>%
    droplevels("DMC") %>%
    complete(Experiment, Time) %>%
    fill(n_uniq_ARGS, n_uniq_Spec, n_reads, .direction = "down") %>%
    pivot_longer(names_to = "Variable", values_to = "Value", cols = c("n_uniq_ARGS", "n_uniq_Spec", "n_reads"))

hex <- hue_pal()(6)

var.labs <- c(n_reads = "ARG counts", n_uniq_ARGS = "Unique ARGs detected", n_uniq_Spec = "Unique species detected")
p_time <- combined %>% ggplot(aes(x = Time, y = Value, col = Experiment)) +
  facet_wrap(facets = ~ Variable, ncol = 1, scales="free", labeller = labeller(Variable = var.labs)) +
  scale_color_manual(values = hex) +
  geom_line(linewidth=1) +
  ylab(NULL) +
  xlim(0, 72) +
  ylim(0, NA) +
  theme_classic() +
  theme(legend.position = "bottom",
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size=14),
        axis.title.x = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size=14))
p_time

# Split the facets into separate ggplot objects
p_list <- lapply(unique(combined$Variable), function(var) {
  p_subset <- p_time %+% subset(combined, Variable == var)
  p_subset
})
# Access individual subplots as separate ggplot objects
p_unARG <- p_list[[1]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 16)) + ylab("Unique species")
p_unSpec <- p_list[[2]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 16)) + ylab('Unique ARGs')
p_ARGcount <- p_list[[3]] + theme(axis.title.y = element_text(size = 16)) + ylab("ARG read count")

p_time2 <- ggarrange(p_unSpec,
                 p_unARG,
                 p_ARGcount,
                 labels = c("a", "b", "c"),
                 nrow = 3, align = "v", common.legend = TRUE, legend = "bottom")
p_time2

# Heatmap with number of resistance genes per experiment
# ARG counts by resistance class and experiment
DMC <- merge(reference_data, reference_comm, by.x="Species", by.y="Organism", all.x=TRUE) %>%
  mutate(Experiment = "DMC", n_reads = Perc_gDNA)
ARG_exp <- rbindlist(list(GMSdb, DMC), fill=TRUE)  
ARG_exp <- merge(ARG_exp, resistance_classes, by.x="ARG", by.y="Gene", all.x=TRUE) 

ARG_exp_byARG <- ARG_exp %>%
  group_by(Experiment, ARG, Resistance) %>%
  summarize(n_reads = sum(n_reads))

ARG_exp_byARG <- ARG_exp_byARG %>% filter(Experiment != "DMC")
midpoint <-  log10(median(ARG_exp_byARG$n_reads))
ARG_exp_byARG$Resistance <- factor(ARG_exp_byARG$Resistance)
order_ARGs <- ARG_exp_byARG %>% filter(Experiment != "DMC") %>% group_by(ARG) %>% 
    mutate(med_reads = median(n_reads)) %>% 
    arrange(desc(Resistance), med_reads) %>% select(ARG) %>% unique()
order_ARGs <- order_ARGs$ARG


ARG_exp_byARG$ARG <- ordered(ARG_exp_byARG$ARG, levels = order_ARGs)
heatmap <- ARG_exp_byARG %>%
 ggplot(aes(x = Experiment, y = ARG, fill= log10(n_reads))) +
  geom_tile() +
  geom_text(aes(label = n_reads), size = 3) +
   scale_fill_gradient2(name="Log10 of ARG counts",
                       low="blue", mid="#dbdbdb", high="red", midpoint = midpoint) +
   scale_y_discrete(guide = guide_axis(position = "right")) +
  theme_classic() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)
heatmap

ARGclasses <- ggplot(ARG_exp_byARG, aes(x=1, y=ARG, fill=Resistance)) +
geom_tile() +
theme_void() +
guides(fill = guide_legend(position = "bottom", ncol=3, nrow=3)) +
theme(legend.title = element_blank(),
        legend.text = element_text(size=14)) +
scale_fill_brewer(palette = "Set3")

legend_heatmap <- get_legend(heatmap)
legend_ARGclass <- get_legend(ARGclasses)
arrange_legends <- cowplot::plot_grid(legend_heatmap, legend_ARGclass, ncol=2, rel_widths = c(0.5,1))

heatmap_ARG <- ggarrange(heatmap + theme(legend.position = "none"), ARGclasses + theme(legend.position = "none"),
                         ncol=2, align = "h", widths = c(20, 1))
heatmap_ARG <- ggarrange(heatmap_ARG, arrange_legends, nrow=2, align = "v", heights = c(10, 1))
heatmap_ARG

fig <- ggarrange(p_time2, heatmap_ARG, ncol = 2, align = "v", widths = c(0.7, 1), labels = c("", "d"))
fig

# absolute numbers of ARGs detected throughout experiments

# barchart with ARGs detected, separated by exp and colored by resistance class
# absolute numbers
ARG_exp_byclass <- ARG_exp %>%
  group_by(Experiment, Resistance) %>%
  summarize(n_reads = sum(n_reads))

ARG_exp_byclass %>% filter(Experiment != "DMC") %>% 
  ggplot(aes(x = Experiment, y = n_reads, fill= Resistance)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90))

# relative numbers
ARG_exp_byclass %>% ggplot(aes(x = Experiment, y = n_reads, fill= Resistance)) +
  geom_col(position = "fill") +
  theme(axis.text.x = element_text(angle = 90))