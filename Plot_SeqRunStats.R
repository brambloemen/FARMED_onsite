# libraries ---------------------------------------------------------------
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

# functions ---------------------------------------------------------------

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

calculate_N50 <- function(lengths) {
  # Sort the read lengths in descending order
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  
  # Calculate the cumulative sum of the read lengths
  cumulative_sum <- cumsum(as.numeric(sorted_lengths))
  
  # Find the index of the read length that is greater than or equal to half the total length
  n50_index <- which(cumulative_sum >= (sum(lengths) / 2))[1]
  
  # Return the N50 value
  return(sorted_lengths[n50_index])
}

# data processing ---------------------------------------------------------

readstats <- cat_tsv("results/readstats/", ".tsv")

readstats$Experiment <- ordered(readstats$Experiment, levels=c("UEQL", "EQL", "EQR", "BDR", "BQR", "BQV"))

# reference community data
ref_community <- read.csv("./zymo_GMS.csv", sep = ";", dec = ",")
ref_community <- ref_community %>%
  arrange(desc(Perc_gDNA), Organism)

# sort by descending zymo GMS gDNA% -> make organisms into ordered factors
ref_community$Organism <- factor(ref_community$Organism, levels = ref_community$Organism,
                                 ordered = TRUE)

# Process combined KMA results summaries
KMAdata <- fread("results/mapping/kma_GMS_summary.tsv", sep = "\t", integer64 = "numeric")
KMAdata <- KMAdata %>%
  mutate(Experiment = str_remove(Experiment, "_GMSdb"))

KMAdata <- merge(KMAdata, ref_community, by.x="Species", by.y = "Organism", all = TRUE)
KMAdata$Experiment <- ordered(KMAdata$Experiment, levels=c("UEQL", "EQL", "EQR", "BDR", "BQR", "BQV"))

# filtering and sorting to plot KMA classification results on heatmap
KMAdata <- KMAdata %>%
  filter(p_bpTotal > 0, !is.na(Perc_gDNA)) %>%
  arrange(desc(Perc_gDNA)) %>%
  group_by(Experiment) %>%
  mutate(p_bpTotal = p_bpTotal/sum(p_bpTotal),
            Perc_gDNA = Perc_gDNA) %>%
  mutate(norm_abundance=log(100*p_bpTotal/Perc_gDNA))
KMAdata$perc <- signif(KMAdata$Perc_gDNA, digits=2)
KMAdata$format_perc <- ifelse(KMAdata$perc <0.05, sprintf("%.2e",KMAdata$perc), KMAdata$perc)
KMAdata$Species_ab <- str_replace(KMAdata$Species, "(\\w)\\w*\\s(\\w)\\w*", "\\1.\\2.")
KMAdata$Species_ab <- paste0(KMAdata$Species_ab, ": ", KMAdata$format_perc, "%")
KMAdata$Species_ab <- ordered(KMAdata$Species_ab, levels=unique(KMAdata$Species_ab))


# Plots -------------------------------------------------------------------
N50s <- readstats %>% group_by(Experiment) %>% reframe(N50 = calculate_N50(lengths))

readlengths <- readstats %>% ggplot(aes(x = Experiment, y = lengths)) +
  geom_violin(position = position_dodge(0.7), alpha=0.75, fill = "#00BFC4") +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.7), width=0.2, alpha=0.75, color = "black") +
  geom_text(data = N50s, aes(x = Experiment, y = 10^5, label = N50)) +
  ylab("Read length (bp)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=14), axis.title.y = element_text(size=14),
      plot.margin = margin(b=0.5, unit = 'cm')) +
  scale_y_log10()

seq_yield <- readstats %>% select(-quals) %>% group_by(Experiment) %>% reframe(Totalbp = sum(lengths))
seq_yield_plot <- seq_yield %>% ggplot(aes(x = Experiment, y = Totalbp)) +
  geom_col(position = position_dodge(0.7), width = 0.5, fill = "#00BFC4") +
  ylab("Total sequenced basepairs") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size = 14),
    plot.margin = margin(b=0.5, unit = 'cm'))


fig1 <- ggarrange(readlengths,
                 seq_yield_plot,
                 labels = c("a", "b"),
                 nrow = 2, align = "v", common.legend = TRUE)
fig1


# heatmap
heatmap <- KMAdata %>%
  ggplot(aes(x=Experiment, y=Species_ab)) +
  geom_tile(aes(fill=norm_abundance)) +
  scale_fill_gradient2(name="Log ratio \nobserved/expected",
                       low="blue", mid="#dbdbdb", high="red", midpoint = 0) +
  # reverse organism order on heatmap y axis (i.e. top=high abundance, bottom=low)
  scale_y_discrete(limits=rev, name=NULL, position = "right") + 
  scale_x_discrete(position = "bottom", name=NULL) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust=1, size=14),
        axis.text.x = element_text(size=14),
        legend.position="top", 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.spacing.x = unit(1, "cm"),
        plot.margin = margin(b=0.5, r=1, l=0.5, unit = 'cm'))

fig2 <- ggarrange(fig1, heatmap, ncol = 2, labels = c("", "c"))
fig2


pdf(NULL)
exportpath <- paste0(getwd(), "/results/plots/SeqRunStats.png")  # !!! always same filename; will overwrite previous ones
print(fig2)
ggsave(exportpath, width = 16, height = 9, dpi = 300)



# KMA against abfhpv database
# Process combined KMA results summaries
KMAdata2 <- fread("results/mapping/kma_abfhpv_summary.tsv", sep = "\t", integer64 = "numeric")
KMAdata2 <- KMAdata2 %>%
  mutate(Experiment = str_remove(Experiment, "_abfhpv"))

KMAdata2 <- merge(KMAdata2, ref_community, by.x="Species", by.y = "Organism", all = TRUE)
KMAdata2$Experiment <- ordered(KMAdata2$Experiment, levels=c("UEQL", "EQL", "EQR", "BDR", "BQR", "BQV"))

# filtering and sorting to plot KMA classification results on heatmap
KMAdata2_ref <- KMAdata2 %>%
  filter(p_bpTotal > 0, !is.na(Perc_gDNA)) %>%
  arrange(desc(Perc_gDNA)) %>%
  mutate(p_bpTotal = p_bpTotal/sum(p_bpTotal),
            Perc_gDNA = Perc_gDNA) %>%
  mutate(norm_abundance=log(100*p_bpTotal/Perc_gDNA))
KMAdata2_ref$Species_ab <- paste0(KMAdata2_ref$Species, " - ", signif(KMAdata2_ref$Perc_gDNA, digits = 2), "%")
KMAdata2_ref$Species_ab <- ordered(KMAdata2_ref$Species_ab, levels=unique(KMAdata2_ref$Species_ab))


# heatmap
heatmap_abfpv <- KMAdata2_ref %>%
  ggplot(aes(x=Experiment, y=Species_ab)) +
  geom_tile(aes(fill=norm_abundance)) +
  scale_fill_gradient2(name="Log10 observed/expected",
                       low="blue", mid="#dbdbdb", high="red", midpoint = 0) +
  # reverse organism order on heatmap y axis (i.e. top=high abundance, bottom=low)
  scale_y_discrete(limits=rev, name=NULL) + 
  scale_x_discrete(position = "bottom", name=NULL) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position="top")


# heatmap all, no ref community
KMAdata2_all <- KMAdata2 %>%
  filter(!str_detect(Species, "[P|p]lasmid") & !str_detect(Species, "[P|p]hage") & !str_detect(Species, "[V|v]irus")) %>%
  filter(mean_queryID > 80, mean_template_coverage > 5, is.na(Perc_gDNA))

midpoint <- median(log(KMAdata2_all$p_bpTotal), na.rm=TRUE)


KMAdata2_all <- KMAdata2 %>%
  filter(!str_detect(Species, "[P|p]lasmid") & !str_detect(Species, "[P|p]hage") & !str_detect(Species, "[V|v]irus")) %>%
  filter(is.na(Perc_gDNA), mean_queryID > 80, mean_template_coverage>5)

UEQL_species <- KMAdata2_all %>% filter(Experiment=="UEQL", mean_queryID > 80, mean_template_coverage>5) %>% select(Species)


heatmap_abfpv_all <- KMAdata2_all %>%
  filter(Species %in% UEQL_species$Species) %>%
  ggplot(aes(x=Experiment, y=Species)) +
  geom_tile(aes(fill=log(p_bpTotal))) +
  scale_fill_gradient2(name="Log10 of proportion of mapped bases",
                       low="blue", mid="#dbdbdb", high="red", midpoint = midpoint) +
  # reverse organism order on heatmap y axis (i.e. top=high abundance, bottom=low)
  scale_y_discrete(limits=rev, name=NULL) + 
  scale_x_discrete(position = "bottom", name=NULL) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position="top")
heatmap_abfpv_all


pdf(NULL)
exportpath <- paste0(getwd(), "/results/plots/heatmap_background_species.png")  # !!! always same filename; will overwrite previous ones
print(heatmap_abfpv_all)
ggsave(exportpath, width = 16, height = 9, dpi = 300)