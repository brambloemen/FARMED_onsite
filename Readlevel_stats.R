library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

readstats <- list.files("../Plots/Summarized_data")
readstats <- readstats[str_detect(readstats,"\\w+.readstats.class.tsv")]


ref_com <- read.csv("zymo_GMS.csv", sep=";", dec=",")

# Parse list of readstat files, extract experiment name from filename, add to each readstat dataset
exp_names <- lapply(readstats, str_extract, "\\w+.readstats.class.tsv$")
exp_names <- lapply(exp_names, str_remove, ".readstats.class.tsv")
readstats <- lapply(readstats, function(x) paste0("../Plots/Summarized_data/", x))
readstats <- lapply(readstats, read.csv, sep="\t")
for (k in 1:length(readstats)){
  readstats[[k]]$Experiment <- exp_names[[k]]
}
readstats <- rbindlist(readstats)

# Parse reference community data
ref_com$Species <- ref_com$Organism
ref_com <- ref_com %>% select(Species, Perc_gDNA) %>% 
  distinct() %>%
  arrange(desc(Perc_gDNA))
ref_Species <- ordered(ref_com$Species)

# Function to clean organism names (to species level)
# clean_org_name <- function(organism){
#   organism <- str_remove_all(organism, "Synthetic|\\[|\\]|\\scf.")S 
#   organism <- str_replace_all(organism, "_", " ")
#   organism <- str_replace_all(organism, "E.coli.*", "Escherichia coli")
#   organism <- str_replace_all(organism, "veillonella", "Veillonella")
#   organism <- str_replace_all(organism, "\\.", " ")
#   organism <- str_extract(organism, "[:upper:]{1}[:lower:]+\\s[:alpha:]+\\.?")
#   organism <- str_replace(organism, "Clostridium difficil+e", "Clostridioides difficile")
#   return(organism)
# }

# preprocess data
readstats <- readstats %>%
  mutate(Species = Template) %>%
  filter(Species %in% ref_com$Species) %>%
  mutate(Species = ordered(Species, levels=ref_Species)) %>%
  group_by(Species, Experiment)

# barplot mean readlengths
plot <- readstats %>%
  summarize(mean_readlength = mean(Length)) %>%
  ggplot(aes(x=Species, y=mean_readlength, fill=Experiment)) +
  geom_col(position = "dodge", color="black") +
  coord_cartesian(ylim=c(0,15000)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

# pdf(NULL)
# exportpath <- paste0(getwd(), "/Mean_readlengths_by_species.png")  # !!! always same filename; will overwrite previous ones
# print(exportpath)
# print(plot)
# ggsave(exportpath, width = 16, height = 9, dpi = 100)

# boxplot
boxplot <- readstats %>%
  ggplot(aes(x=Species, y=Length, fill=Experiment)) +
  geom_boxplot(position = "dodge", outlier.shape = NA) +
  # coord_cartesian(ylim=c(0,20000)) +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

# pdf(NULL)
# exportpath <- paste0(getwd(), "/Boxplot_readlengths_by_species.png")  # !!! always same filename; will overwrite previous ones
# print(exportpath)
# print(boxplot)
# ggsave(exportpath, width = 16, height = 9, dpi = 100)