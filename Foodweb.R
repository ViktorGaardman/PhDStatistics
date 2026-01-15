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

library(readxl)
library(readr)

matrix <- read_xlsx("mergedData_wide_PresenceAbsence_no_zotu.xlsx")

diet_matrix <- matrix[, 9:489]

prey_data <- matrix[, 1:8]

Insect_ID <- read_xlsx("mergedData_long_format_no_zotu.xlsx")

Base_ID <- read.csv("Foodweb_samples_2024.csv", sep = ";")

Base_ID <- Base_ID %>%
  mutate(
    Date = str_pad(Date, width = 6, side = "left", pad = "0")  # make all dates 6 digits
  )

# Create month column
Base_ID <- Base_ID %>%
  mutate(
    Month = str_sub(Date, 3, 4),  # extract MM from DDMMYY
    Month = case_when(
      Month == "06" ~ "June",
      Month == "07" ~ "July",
      Month %in% c("08", "09") ~ "August",
      TRUE ~ NA_character_  # catch other months if present
    )
  )

Bird_ID <- Base_ID %>%
  filter(Group == "Bird") %>%
  select(ID, Species) %>%
  rename(
    Sample = ID,
    Genus   = Species
  )

#Remove the bloodmeal samples for now
excluded_primers <- c("P16S", "12S-V5")

Insect_Filter <- Insect_ID %>%
  filter(!PrimerPair %in% excluded_primers)

#Sum reads across remaining primers
summed_reads <- Insect_Filter %>%
  group_by(Sample, Family, Genus) %>%
  summarise(total_reads = sum(Count), .groups = "drop")

#Identify families that cannot be predators
non_predator_species <- c(
  "Chironomidae",
  "Geometridae",
  "Carabidae",
  "Unknown",
  "Noctuidae",
  "Simuliidae",
  "Byrrhidae",
  "Culicidae",
  "Curculionidae",
  "Ichneumonidae",
  "Empididae",
  "Tipulidae",
  "Syrphidae",
  "Oribatulidae",
  "Tortricidae"
)

dominant_species <- summed_reads %>%
  filter(!Family %in% non_predator_species) %>% 
  group_by(Sample) %>%
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Sample, Genus)

#Fill in the bird samples

dominant_no_birds <- dominant_species %>%
  filter(!Sample %in% Bird_ID$Sample)

predator_lookup <- bind_rows(
  dominant_no_birds,
  Bird_ID %>% select(Sample, Genus)
)

#Remove the bloodmeal samples from diel matrix
diet_matrix_clean <- diet_matrix[, colnames(diet_matrix) %in% predator_lookup$Sample]

# 1. Store original column names
original_colnames <- colnames(diet_matrix_clean)

#Replace names using the lookup vector
lookup_vector <- setNames(predator_lookup$Genus, predator_lookup$Sample)
new_colnames <- lookup_vector[original_colnames]

#Assign new column names (some may be NA)
colnames(diet_matrix_clean) <- new_colnames

#Identify which columns became NA and backtrack their original Sample IDs
na_columns <- which(is.na(new_colnames))
if(length(na_columns) > 0){
  missing_samples <- original_colnames[na_columns]
  warning("These columns did not match any Sample ID in predator_lookup: ", 
          paste(missing_samples, collapse = ", "))
}



#No NAs - we are good to go! :D

write_csv(diet_matrix_clean, "diet_matrix_renamed.csv")

#Add prey names to the rows
prey_names <- prey_data$Assigned

#Store original column names
original_colnames <- colnames(diet_matrix_clean)

# Make them unique internally
colnames(diet_matrix_clean) <- make.unique(original_colnames, sep = "_sample_")

#Add prey names as a column
diet_long <- diet_matrix_clean %>%
  as.data.frame() %>%         # ensure it's a dataframe
  mutate(Prey = prey_names) %>%
  pivot_longer(
    cols = -Prey,             # keep Prey column
    names_to = "PredatorSample",
    values_to = "Presence"
  )

#Add original names back
diet_long <- diet_long %>%
  mutate(Predator = sub("_sample_.*", "", PredatorSample))

#Sum presence per Predator-Prey combination
interaction_counts <- diet_long %>%
  group_by(Predator, Prey) %>%
  summarise(TotalPresence = sum(Presence), .groups = "drop")

# 3. Calculate total number of samples per predator
total_samples <- diet_long %>%
  group_by(Predator) %>%
  summarise(TotalSamples = n(), .groups = "drop")

# 4. Join totals and calculate interaction strength
edgelist <- interaction_counts %>%
  left_join(total_samples, by = "Predator") %>%
  mutate(InteractionStrength = TotalPresence / TotalSamples) %>%
  select(Predator, Prey, InteractionStrength)

# 5. Optional: remove zero interactions (if you only want non-zero links)
edgelist <- edgelist %>% filter(InteractionStrength > 0)

#Check list
head(edgelist)

edgelist <- edgelist %>%
  filter(!Predator %in% "Unknown")

# Create a new column with prey genus (everything before the first underscore)
edgelist <- edgelist %>%
  mutate(PreyGenus = sub("_.*", "", Prey))

#Remove predator DNA/cannibalism
edgelist_clean <- edgelist %>%
  filter(Predator != PreyGenus) %>%  
  select(-PreyGenus)                

write.csv(edgelist_clean, "edgelist_2024.csv")

###########
#Actual fooweb analysis!!
#Basics for presentation

library(bipartite)
library(igraph)

#Rename columns for bipartite
edgelist_for_bipartite <- edgelist_clean %>%
  rename(
    higher = Predator,  # predator
    lower  = Prey,      # prey
    freq   = InteractionStrength
  )


weblist <- frame2webs(edgelist_for_bipartite, 
                      varnames = c("higher", "lower", "freq"))

# Convert long edgelist to matrix format
network_matrix <- edgelist_clean %>%
  pivot_wider(
    names_from = Predator,
    values_from = InteractionStrength,
    values_fill = 0
  ) %>%
  column_to_rownames("Prey") %>%  # rows = prey
  as.matrix()

# Check
dim(network_matrix)
head(network_matrix)

basicweb <- plotweb(
  network_matrix,
  srt = 45,
  text_size = 0.5,
  sorting = "dec"
)



low.abun = TRUE,       # scale prey bars by total interactions
high.abun = TRUE,      # scale predator bars by total interactions
labsize = 0.7,         # reduce label size
low.col = "orange",    # prey bar color
high.col = "skyblue",  # predator bar color
col.interaction = "gray" # links color



#LEFTOVER CODE
