# libraries --------------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales)

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
  filter(!str_detect(Species, "([P|p]lasmid)|([V|v]irus)|([P|p]hage)"))

abfhpv_SpeciesByTime <- abfhpv %>% filter(Organism_QID >= 0.80) %>% rarefaction_byExpTime(Species)
abfhpv_SpeciesByTime <- abfhpv %>% filter(Organism_QID >= 0.80, Coverage >= 5) %>% rarefaction_byExpTime(Species)

p_abfhpv_uniqueSpecies <- ggplot(abfhpv_SpeciesByTime, aes(x=Time, y=n_unique)) +
    geom_line(aes(color = Experiment)) + ylim(0,NA)

p_abfhpv_uniqueSpecies

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
  mutate(Time = Time_hours)

# Unique ARGs by time per experiment
GMSdb_unique_ARGs <- GMSdb %>% rarefaction_byExpTime(ARG)
p_GMSdb_uniqueARGs <- lineplot(GMSdb_unique_ARGs, n_unique)
p_GMSdb_uniqueARGs

# ARG counts by resistance class and experiment
RC <- merge(reference_data, reference_comm, by.x="Species", by.y="Organism", all.x=TRUE) %>%
  mutate(Experiment = "RC", n_reads = Perc_gDNA)
ARG_exp <- rbindlist(list(GMSdb, RC), fill=TRUE)  %>%
  mutate(Experiment = case_when(
           Experiment == "Fecal_background" ~ "UEQL",
           Experiment == "Fecal_GMSspikeI_CB_beads_bc2" ~ "BQR",
           Experiment == "Fecal_GMSspikeI_CB_beads_VSK" ~ "BQV",
           Experiment == "Fecal_GMSspikeI_CB_RAD" ~ "BDR",
           Experiment == "Fecal_GMSspikeI_QuickDNA_LSK" ~ "EQL",
           Experiment == "Fecal_GMSspikeI_QuickDNA_RAD" ~ "EQR",
           .default = "RC"
         ))
ARG_exp$Experiment <- ordered(ARG_exp$Experiment, levels=c("UEQL", "RC", "EQL", "EQR", "BDR", "BQR", "BQV"))

ARG_exp <- merge(ARG_exp, resistance_classes, by.x="ARG", by.y="Gene", all.x=TRUE) 
ARG_exp_byclass <- ARG_exp %>%
  group_by(Experiment, Resistance) %>%
  summarize(n_reads = sum(n_reads))
ARG_exp_byARG <- ARG_exp %>%
  group_by(Experiment, ARG, Resistance) %>%
  summarize(n_reads = sum(n_reads))

# barchart with ARGs detected, separated by exp and colored by resistance class
# absolute numbers
ARG_exp_byclass %>% filter(Experiment != "RC") %>% 
  ggplot(aes(x = Experiment, y = n_reads, fill= Resistance)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90))

# relative numbers
ARG_exp_byclass %>% ggplot(aes(x = Experiment, y = n_reads, fill= Resistance)) +
  geom_col(position = "fill") +
  theme(axis.text.x = element_text(angle = 90))

  # Heatmap with number of resistance genes per experiment
midpoint <-  log10(median(ARG_exp_byARG$n_reads))
ARG_exp_byARG %>% filter(Experiment != "RC") %>% 
  ggplot(aes(x = Experiment, y = ARG, fill= log10(n_reads))) +
  geom_tile() +
  geom_text(aes(label = n_reads)) +
   scale_fill_gradient2(name="Log10 of reads",
                       low="blue", mid="#e6e6e6", high="red", midpoint = midpoint) +
  theme(axis.text.x = element_text(angle = 90))


# absolute numbers of ARGs detected throughout experiments
ARG_byExpTime <- GMSdb %>%
  arrange(Experiment, Time) %>%
    group_by(Experiment) %>%
    reframe(n_reads = cumsum(n_reads), Time = Time) %>%
    group_by(Experiment, Time) %>%
    slice_tail(n=1) %>% ungroup() %>%
    complete(Experiment, Time) %>%
    fill(n_reads, .direction = "down")

ARG_byExpTime %>%
    ggplot(aes(x=Time, y=n_reads, col=Experiment)) +
    geom_line()