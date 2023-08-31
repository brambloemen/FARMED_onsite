# libraries --------------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales)
library(circlize)
# library(networkD3)


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

cumsum_reads <- function(df, ...){
  df <- df %>% 
    arrange(Time_hours, ...) %>%
    group_by(...) %>%
    reframe(Time_hours = Time_hours,
              ... = ...,
              n_reads = cumsum(n_reads),
              n_match_bases1 = cumsum(n_match_bases1))
}

# read data --------------------------------------------------------------------

# reference community data
reference_data <- read.csv("./GMSspikeI_AMRprofile.csv", sep = ";") %>%
    mutate(Reference_AMR_link=TRUE,
           Species = Organism) %>%
    select(-Organism)
reference_ARGs <- data.frame(reference_data$ARG) %>% mutate(Reference_AMR=TRUE)

reference_comm <- read.csv("./zymo_GMS.csv", sep = ";", dec=",")
# resistance classes
resistance_classes <- read.csv("./resistance_classes.tsv", sep="\t") %>%
  mutate(Resistance = str_remove(Resistance, "\\sresistance")) %>% select(Gene, Resistance)


abfhpv <- cat_tsv("../Plots/Summarized_data/AMRlinking/", "_abfhpvvResF.tsv") %>%
  mutate(Experiment = case_when(
           Experiment == "Fecal_background" ~ "FbEQL",
           Experiment == "Fecal_GMSspikeI_CB_beads_bc2" ~ "FBQR",
           Experiment == "Fecal_GMSspikeI_CB_beads_VSK" ~ "FBQV",
           Experiment == "Fecal_GMSspikeI_CB_RAD" ~ "FBDR",
           Experiment == "Fecal_GMSspikeI_QuickDNA_LSK" ~ "FEQL",
           Experiment == "Fecal_GMSspikeI_QuickDNA_RAD" ~ "FEQR"
         )) %>%
  mutate(Species = ifelse(str_detect(Species, "[P|p]lasmid"), "Plasmid", Species))


GMSdb <- cat_tsv("../Plots/Summarized_data/AMRlinking/", "_GMSvResF.tsv") %>%
  mutate(Experiment = case_when(
           Experiment == "Fecal_background" ~ "FbEQL",
           Experiment == "Fecal_GMSspikeI_CB_beads_bc2" ~ "FBQR",
           Experiment == "Fecal_GMSspikeI_CB_beads_VSK" ~ "FBQV",
           Experiment == "Fecal_GMSspikeI_CB_RAD" ~ "FBDR",
           Experiment == "Fecal_GMSspikeI_QuickDNA_LSK" ~ "FEQL",
           Experiment == "Fecal_GMSspikeI_QuickDNA_RAD" ~ "FEQR"
         )) %>%
  mutate(Species = ifelse(str_detect(Species, "[P|p]lasmid"), "Plasmid", Species))


chordplot_AMRlinks <- function(df, Exp, readsThreshold=1, bpThreshold=50000){
 # reformat dataframe: keep only the cumulative counts at final time;
  # filter with threshold on how much of the ARG reads need to be aligned to the taxon template
  # !TODO: fix so that ARG-host under treshold are cateogrized as unclassified hosts
  df_ARGlink <- df %>%
      filter(Experiment == Exp) %>%
      cumsum_reads(Species, ARG) %>%
      mutate(ARGlink = paste(Species, ARG, sep = " - ")) %>%
      arrange(Time_hours) %>% group_by(Species, ARG) %>% slice_tail(n = 1) %>% 
      # filter(Species == "Unmapped" | n_match_bases1>bpThreshold, n_reads > readsThreshold) %>%
      mutate(Species = ifelse(n_match_bases1>bpThreshold & n_reads > readsThreshold, Species, "Unmapped")) %>%
      group_by(Species, ARG) %>% reframe(n_reads = sum(n_reads)) %>%
      mutate(Species = ifelse(Species %in% reference_comm$Organism, paste0(Species, "*"), Species))
  
  # reformat to matrix for chordDiagram
  m_ARGlink <- df_ARGlink %>% select(Species, ARG, n_reads) %>% pivot_wider(names_from="ARG", values_from="n_reads") %>% as.matrix()
  m_ARGlink2 <- m_ARGlink[1:nrow(m_ARGlink),2:ncol(m_ARGlink)]
  # coerce to numeric
  m_ARGlink2 <- matrix(as.numeric(m_ARGlink2), ncol = ncol(m_ARGlink2))
  rownames(m_ARGlink2) <- m_ARGlink[1:nrow(m_ARGlink),1]
  colnames(m_ARGlink2) <- colnames(m_ARGlink)[2:length(colnames(m_ARGlink))]
  # chordDiagram only takes table or df as input, not matrix
  m_ARGlink2 <- as.table(m_ARGlink2)

  RC_spec_star <- paste0(reference_comm$Organism, "*")
  refspecs <- rownames(m_ARGlink2)[rownames(m_ARGlink2) %in% RC_spec_star]
  border_mat <- m_ARGlink2
  border_mat[] <- NA
  border_mat[refspecs,1:ncol(m_ARGlink2)] <- "black"

  # chordDiagram(m_ARGlink2, link.border = border_mat, link.lwd = 2)
  chordDiagram(m_ARGlink2, link.border = border_mat, link.lwd = 2,annotationTrack = "grid", preAllocateTracks = 1)
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.par(canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2))
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  
  
  circos.clear()
}

chordplot_AMRlinks(abfhpv, "FbEQL")
chordplot_AMRlinks(abfhpv, "FEQL")
chordplot_AMRlinks(abfhpv, "FEQR")
chordplot_AMRlinks(abfhpv, "FBDR")
chordplot_AMRlinks(abfhpv, "FBQR")
chordplot_AMRlinks(abfhpv, "FBQV")

chordplot_AMRlinks(GMSdb, "FBQV")



# sankey_AMRlinks <- function(df, Exp, kmares, QID = 90, cov = 20){
#   links <- data.frame(
#     source = df$Species,
#     target = df$ARG,
#     value = df$n_reads,
#     link = df$ARGlink)
#   nodes <- data.frame(
#     name = unique(c(df$Species, df$ARG)))

#   # sankeyNetwork() requires integers (0-indexed) -> match source vs nodes_name to get integer ID
#   links$IDsource <- match(links$source, nodes$name)-1 
#   links$IDtarget <- match(links$target, nodes$name)-1

#   sankey <- sankeyNetwork(
#     Links = links,
#     Nodes = nodes,
#     Source = "IDsource",
#     Target = "IDtarget",
#     Value = "value",
#     NodeID = "name",
#     sinksRight = TRUE,
#     width = 500,
#     height = 1600
#   )
#   saveNetwork(sankey, "test.html", selfcontained = FALSE)
#   webshot::webshot("test.html","test.png", vwidth = 500, vheight = 1600)
# }
# sankey_AMRlinks(abfhpv, "FBQV", KMAdata_abfhpv)


