library("tidyverse")

rm(list = ls())

df <- read.csv("MalaiseSamples_2024.csv", 
               sep=";", header=TRUE, as.is=TRUE)

ggplot(data = df, aes(x = Week, y = Count, by = Taxon)) + 
  geom_point(aes(color = Taxon))

#1. Look at number of active taxa

active_taxa <- df %>%
  filter(Count > 0) %>%
  distinct(Week, Taxon) %>%
  count(Week, name = "active_families")

taxa_active <- ggplot(data = active_taxa, aes(x = Week, y = active_families)) +
  geom_point()+
  geom_smooth(span= 1.2) +
  scale_x_continuous(limits = c(23, 35), breaks = seq(23, 35, by = 2)) +
  theme_minimal() +
  ylab("Number of active groups") +
  xlab("Week") +
  geom_vline(xintercept = 27, linetype = "dashed", color = "black", alpha = 0.6) +
  geom_vline(xintercept = 31, linetype = "dashed", color = "black", alpha = 0.6) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )

taxa_active

ggsave(plot = taxa_active, filename = "Active_taxa.jpg", width = 6.5,
       height = 5.26, dpi = 450)

#2. Look at abundances of most frequent taxa across entire period

weekly_family <- df %>%
  group_by(Week, Taxon) %>%
  summarise(
    total_count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

top5_week <- weekly_family %>%
  group_by(Week) %>%
  slice_max(total_count, n = 5, with_ties = FALSE) %>%
  ungroup()
active_taxa <- df %>%
  filter(Count > 0) %>%
  distinct(Week, Taxon) %>%
  count(Week, name = "active_families")

top5_active <- ggplot(data = top5_week,
                                    aes(x = Week, y = log(total_count), 
                                            by = Taxon)) +
  geom_point(aes(color = Taxon), size = 2) +
  geom_line(aes(color = Taxon), linewidth = 1) +
  scale_x_continuous(limits = c(23, 35), breaks = seq(23, 35, by = 2)) +
  theme_minimal() +
  ylab("ln(Abundance)") +
  xlab("Week") +
  geom_vline(xintercept = 27, linetype = "dashed", color = "black", alpha = 0.6) +
  geom_vline(xintercept = 31, linetype = "dashed", color = "black", alpha = 0.6) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )

top5_active

ggsave(plot = top5_active, filename = "top5_active.jpg", width = 6.5,
       height = 5.26, dpi = 450)

#3. All taxa but with total abundance > 5 to remove rare ones

all5_taxa <- df %>%
  group_by(Week, Taxon) %>%
  summarise(
    total_count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_count >= 5)

taxa_all <- ggplot(data = all5_taxa,
                      aes(x = Week, y = log(total_count), 
                          by = Taxon)) +
  geom_point(aes(color = Taxon), size = 2) +
  geom_line(aes(color = Taxon), linewidth = 1) +
  scale_x_continuous(limits = c(23, 35), breaks = seq(23, 35, by = 2)) +
  theme_minimal() +
  ylab("ln(Abundance)") +
  xlab("Week") +
  geom_vline(xintercept = 27, linetype = "dashed", color = "black", alpha = 0.6) +
  geom_vline(xintercept = 31, linetype = "dashed", color = "black", alpha = 0.6) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )

taxa_all

ggsave(plot = taxa_all, filename = "all5_taxa.jpg", width = 6.5,
       height = 5.26, dpi = 450)

#Total abundance across season

abundance_total <- df %>%
  group_by(Week) %>%
  summarise(
    total_count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

total_abu <- ggplot(data = abundance_total, aes(x = Week, 
                                                y = total_count)) +
  geom_point() +
  geom_smooth(span = 1.25) +
  scale_x_continuous(limits = c(23, 35), breaks = seq(23, 35, by = 2)) +
  theme_minimal() +
  ylab("Total abundance") +
  xlab("Week") +
  geom_vline(xintercept = 27, linetype = "dashed", color = "black", alpha = 0.6) +
  geom_vline(xintercept = 31, linetype = "dashed", color = "black", alpha = 0.6) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )

total_abu

ggsave(plot = total_abu, filename = "total_abundance.jpg", width = 6.5,
       height = 5.26, dpi = 450)  

###############
#NETWORK
##############
#PossiblePackages
#library(cheddar) Load from Git
#library(enaR) Load from Git
library(bipartite)
library(igraph)

matrix <- read.csv("mergedData_wide_PresenceAbsence.csv", sep = ";")
#ISSUE WITH ZOTU LISTS NEEDS TO BE ADRESSED BEFORE WE CONTINUE

diet_matrix <- matrix[, 10:487]

Insect_ID <- read.csv("mergedData_long_format.csv", sep = ";")

Base_ID <- read.csv("Foodweb_samples_2024.csv", sep = ";")

Bird_ID <- Base_ID %>%
  filter(Group == "Bird") %>%
  select(ID, Species) %>%
  rename(
    Sample = ID,
    Family   = Species
  )

#Exchange SampleID for species (predators)

#Remove the bloodmeal samples for now
excluded_primers <- c("P16S", "12SV5")

Insect_Filter <- Insect_ID %>%
  filter(!Primerset %in% excluded_primers)

#Sum reads across remaining primers
summed_reads <- Insect_Filter %>%
  group_by(Sample, Family) %>%
  summarise(total_reads = sum(Count), .groups = "drop")

#Identify families that cannot be predators
non_predator_species <- c(
  "Chironomidae",
  "Geometridae",
  "Carabidae",
  "Unknown",
  "Noctuidae",
  "Simuliidae"
)

dominant_species <- summed_reads %>%
  filter(!Family %in% non_predator_species) %>% 
  group_by(Sample) %>%
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Sample, Family)

#Fill in the bird samples

dominant_no_birds <- dominant_species %>%
  filter(!Sample %in% Bird_ID$Sample)

predator_lookup <- bind_rows(
  dominant_no_birds,
  Bird_ID %>% select(Sample, Family)
)

#Name lookup vector
id_to_species <- setNames(dominant_species$Family,
                          dominant_species$Sample)

#Exchange names
colnames(diet_matrix) <- id_to_species[colnames(diet_matrix)]

#Check missing IDs
setdiff(colnames(diet_mat), names(id_to_species))

#Check duplicates
any(duplicated(colnames(diet_mat)))

#Exclude mosquitoes for bloodmeal analysis

