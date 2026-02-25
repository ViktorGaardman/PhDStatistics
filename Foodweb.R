library(readxl)
library(readr)
library(tidyverse)
library(bipartite)

rm(list = ls())

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
  )%>%
  filter(Comment %in% "")

Bird_ID <- Base_ID %>%
  filter(Group == "Bird") %>%
  select(ID, Species, Month) %>%
  rename(
    Sample = ID,
    Family   = Species,
    Month = Month
  )

NonSpiderPred_ID <- Base_ID %>%
  filter(Group %in% c("Scatophagidae", "Aquatic predator")) %>%
  select(ID, Group, Species, Month) %>%
  rename(
    Sample = ID,
    Group2   = Species,
    Month = Month
  )

NonSpiderPred_ID$Family <- fct_collapse(
  NonSpiderPred_ID$Group2,
  Dytiscidae = c("C. dolabratus larvae", "C. dolabratus", "Colymbetes dolabratus",
                 "C. Dolabratus"),
  Chironomidae = "Chironomid larvae",
  Scathophagidae = "")

Insect_ID <- Insect_ID %>%
  mutate(
    Date = str_pad(Date, width = 6, side = "left", pad = "0")  # make all dates 6 digits
  )

# Create month column
Insect_ID <- Insect_ID %>%
  mutate(
    Month = str_sub(Date, 3, 4),  # extract MM from DDMMYY
    Month = case_when(
      Month == "06" ~ "June",
      Month == "07" ~ "July",
      Month %in% c("08", "09") ~ "August",
      TRUE ~ NA_character_  # catch other months if present
    )
  )

#Remove the bloodmeal samples for now
excluded_primers <- c("P16S", "12S-V5")

Insect_Filter <- Insect_ID %>%
  filter(!PrimerPair %in% excluded_primers)

#Sum reads across remaining primers
summed_reads <- Insect_Filter %>%
  group_by(Sample, Family, Genus, Month) %>%
  summarise(total_reads = sum(Count), .groups = "drop")

#Remove non-spider arthropod predators from field data
dominant_spidersbirds <- summed_reads %>%
  filter(!Sample %in% NonSpiderPred_ID$Sample)

#Remove bird predators, leaving only spiders
dominant_spiders <- dominant_spidersbirds %>%
  filter(!Sample %in% Bird_ID$Sample)

#Identify families that cannot be spider predators (all non-spiders)
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
  "Tortricidae",
  "Dytiscidae",
  "Scathophagidae"
)

dominant_spider_sp <- dominant_spiders %>%
  filter(!Family %in% non_predator_species) %>% 
  group_by(Sample) %>%
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Sample, Family, Month)

birdsspiders <- bind_rows(
  dominant_spider_sp,
  Bird_ID %>% select(Sample, Family, Month)
)

predator_lookup <- bind_rows(
  birdsspiders,
  NonSpiderPred_ID %>% select(Sample, Family, Month)
)

predator_lookup <- predator_lookup %>%
  filter(!Family %in% "Unknown")

#Split matrix by month

make_monthly_diet_matrices <- function(diet_matrix, predator_lookup) {
  
  # Helper function for a single month
  make_one_month <- function(month_name) {
    
    month_df <- predator_lookup %>%
      filter(Month == month_name) %>%
      select(Sample, Family)
    
    # Filter diet matrix
    mat <- diet_matrix[, colnames(diet_matrix) %in% month_df$Sample, drop = FALSE]
    
    # Rename columns Sample -> Genus
    lookup <- setNames(month_df$Family, month_df$Sample)
    colnames(mat) <- lookup[colnames(mat)]
    
    return(mat)
  }

    # Create matrices
  list(
    June   = make_one_month("June"),
    July   = make_one_month("July"),
    August = make_one_month("August")
  )
}

monthly_matrices <- make_monthly_diet_matrices(
  diet_matrix = diet_matrix,
  predator_lookup = predator_lookup
)

diet_june_Fam   <- monthly_matrices$June
diet_july_Fam   <- monthly_matrices$July
diet_august_Fam <- monthly_matrices$August



#Remove the bloodmeal samples from diel matrix
diet_matrix_clean <- diet_matrix[, colnames(diet_matrix) %in% predator_lookup$Sample]

# 1. Store original column names
original_colnames <- colnames(diet_matrix_clean)

#Replace names using the lookup vector
lookup_vector <- setNames(predator_lookup$Family, predator_lookup$Sample)
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

write_csv(diet_matrix_clean, "diet_matrix_Fam.csv")
write_csv(diet_june_Fam, "diet_matrix_june_Fam.csv")
write_csv(diet_july_Fam, "diet_matrix_july_Fam.csv")
write_csv(diet_august_Fam, "diet_matrix_august_Fam.csv")

diet_matrix_clean <- read.csv("diet_Matrix_Fam.csv")

#Add prey names to the rows (here genera)
prey_names <- prey_data$Family

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
  filter(!Prey %in% "NA")

edgelist <- edgelist %>%
  filter(!Prey %in% "Fringillidae")

# Create a new column with prey genus (everything before the first underscore)
#edgelist <- edgelist %>%
#  mutate(PreyGenus = sub("_.*", "", Prey))

#Remove predator DNA/cannibalism
#edgelist_clean <- edgelist %>%
#  filter(Predator != Prey) %>%  
#  select(-Prey)                

#Remove predator DNA/cannibalism (Family level)
edgelist_clean <- edgelist %>%
  filter(Predator != Prey)

write.csv(edgelist_clean, "edgelist_2024_Fam.csv")

#Monthly edgeslists
diet_matrix_to_edgelist <- function(diet_matrix, prey_names) {
  
  # 1. Store original column names
  original_colnames <- colnames(diet_matrix)
  
  # 2. Make column names unique internally
  colnames(diet_matrix) <- make.unique(original_colnames, sep = "_sample_")
  
  # 3. Convert to long format
  diet_long <- diet_matrix %>%
    as.data.frame() %>%
    mutate(Prey = prey_names) %>%
    pivot_longer(
      cols = -Prey,
      names_to = "PredatorSample",
      values_to = "Presence"
    )
  
  # 4. Recover predator genus from column names
  diet_long <- diet_long %>%
    mutate(Predator = sub("_sample_.*", "", PredatorSample))
  
  # 5. Count prey occurrences per predator
  interaction_counts <- diet_long %>%
    group_by(Predator, Prey) %>%
    summarise(TotalPresence = sum(Presence), .groups = "drop")
  
  # 6. Count number of samples per predator
  total_samples <- diet_long %>%
    group_by(Predator) %>%
    summarise(TotalSamples = n_distinct(PredatorSample), .groups = "drop")
  
  # 7. Calculate interaction strength
  edgelist <- interaction_counts %>%
    left_join(total_samples, by = "Predator") %>%
    mutate(InteractionStrength = TotalPresence / TotalSamples) %>%
    select(Predator, Prey, InteractionStrength) %>%
    filter(InteractionStrength > 0)
  
  # 8. Remove predator DNA / cannibalism (family-level match) and NA
  edgelist_clean <- edgelist %>%
    filter(Predator != Prey) %>%
    filter(!Prey %in% "NA")
  
  return(edgelist_clean)
}

edgelist_june_Fam   <- diet_matrix_to_edgelist(diet_june_Fam, prey_data$Family)
edgelist_july_Fam   <- diet_matrix_to_edgelist(diet_july_Fam, prey_data$Family)
edgelist_august_Fam <- diet_matrix_to_edgelist(diet_august_Fam, prey_data$Family)


write.csv(edgelist_june_Fam, "edgelist_june_2024_Fam.csv")
write.csv(edgelist_july_Fam, "edgelist_july_2024_Fam.csv")
write.csv(edgelist_august_Fam, "edgelist_august_2024_Fam.csv")




###########
#Actual foodweb analysis!!

edgelist_clean <- read.csv("edgelist_2024_Fam.csv")
edgelist_june_Fam <- read.csv("edgelist_june_2024_Fam.csv")
edgelist_july_Fam <- read.csv("edgelist_July_2024_Fam.csv")
edgelist_august_Fam <- read.csv("edgelist_august_2024_Fam.csv")

edgelist_clean <- edgelist_clean[,2:4]
edgelist_june_Fam <- edgelist_june_Fam[,2:4]
edgelist_july_Fam <- edgelist_july_Fam[,2:4]
edgelist_august_Fam <- edgelist_august_Fam[,2:4]

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


png("foodweb_basic.png", width = 3000, height = 1500, res = 300)

plotweb(
  network_matrix,
  srt = 45,
  text_size = 0.5,
  sorting = "dec"
)

dev.off()

#Monthly

edgelist_to_matrix <- function(edgelist) {
  
  network_matrix <- edgelist %>%
    pivot_wider(
      names_from  = Predator,
      values_from = InteractionStrength,
      values_fill = 0
    ) %>%
    column_to_rownames("Prey") %>%
    as.matrix()
  
  return(network_matrix)
}

edgelist_june_Fam <- edgelist_june_Fam %>%
  filter(!Prey %in% "NA")
edgelist_july_Fam <- edgelist_july_Fam %>%
  filter(!Prey %in% "NA")
edgelist_august_Fam <- edgelist_august_Fam %>%
  filter(!Prey %in% "NA")
edgelist_august_Fam <- edgelist_august_Fam %>%
  filter(!Prey %in% "Fringillidae")

matrix_june_Fam   <- edgelist_to_matrix(edgelist_june_Fam)
matrix_july_Fam   <- edgelist_to_matrix(edgelist_july_Fam)
matrix_august_Fam <- edgelist_to_matrix(edgelist_august_Fam)

plotweb(
  matrix_june_Fam,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE
)

visweb(
  matrix_june_Fam
)


foodwebs <- list(
  June   = matrix_june_Fam,
  July   = matrix_july_Fam,
  August = matrix_august_Fam
)

for (month in names(foodwebs)) {
  
  png(
    filename = paste0("Fam_foodweb_", month, ".png"),
    width = 3000,
    height = 1500,
    res = 300
  )
  
  plotweb(
    foodwebs[[month]],
    srt = 45,
    text_size = 0.5,
    sorting = "dec"
  )
  
  title(paste("Fam food web –", month))
  
  dev.off()
}


#Weighted analyses
#1. Load malaise and pitfall data
flyingprey <- read.csv("FlyingAbundances2024_edit.csv", sep = ";")
groundprey <- read.csv("Pitfall2024_edited.csv", sep = ";")
groundprey <- groundprey[,1:5]

flying_prey_june <- flyingprey %>%
  filter(Month == "June")

ground_prey_june <- groundprey %>%
  filter(Month == "6")

flying_prey_june_matrix <- flying_prey_june %>%
  select(Taxon, Total) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Total,
    values_fill = 0
  ) %>%
  as.matrix()

#Drop flying arthropods from pitfall traps
ground_prey_june <- ground_prey_june %>%
  filter(!Family %in% c("Scelionidae", "Ichneumonidae",
         "Syrphidae", "Sciaridae", "Phoridae",
         "Muscidae", "Chironomidae", "Anthomyiidae",
         "Psylloidea", "Coccoidea", "Aphidoidea"))



ground_prey_june_matrix <- ground_prey_june %>%
  select(Family, Count) %>%
  pivot_wider(
    names_from = Family,
    values_from = Count,
    values_fill = 0
  ) %>%
  as.matrix()

#2. Split diet into ground dwelling prey and flying prey
edgelist_june_flying <- edgelist_june_Fam %>%
  filter(Prey %in% c("Braconidae", "Chironomidae", "Ichneumonidae",
         "Noctuidae", "Vespidae", "Culicidae"))

edgelist_june_ground <- edgelist_june_Fam %>%
  filter(!Prey %in% c("Braconidae", "Chironomidae", "Ichneumonidae",
                     "Noctuidae", "Vespidae", "Culicidae"))

matrix_june_flying <- edgelist_to_matrix(edgelist_june_flying)
matrix_june_ground <- edgelist_to_matrix(edgelist_june_ground)

# Transpose your matrix so predators are rows and prey are columns
matrix_june_flying_transposed <- t(matrix_june_flying)
matrix_june_ground_transposed <- t(matrix_june_ground)

# Subset flying_prey_june_matrix to only include prey in matrix_june_flying
flying_prey_june_matrix2 <- flying_prey_june_matrix[, colnames(matrix_june_flying_transposed)]
ground_prey_june_matrix2 <- ground_prey_june_matrix[, colnames(matrix_june_ground_transposed)]

flying_prey_june_matrix2 <- matrix(flying_prey_june_matrix2, nrow = 1)
colnames(flying_prey_june_matrix2) <- colnames(matrix_june_flying_transposed)

ground_prey_june_matrix2 <- matrix(ground_prey_june_matrix2, nrow = 1)
colnames(ground_prey_june_matrix2) <- colnames(matrix_june_ground_transposed)

# Get prey names in each matrix
prey_in_matrix <- colnames(flying_prey_june_matrix2)
prey_in_abundance <- colnames(flying_prey_june_matrix)

prey_in_matrix_g <- colnames(ground_prey_june_matrix2)
prey_in_abundance_g <- colnames(ground_prey_june_matrix)


# Find prey in abundance but not in matrix
missing_prey <- setdiff(prey_in_abundance, prey_in_matrix)
missing_prey_g <- setdiff(prey_in_abundance_g, prey_in_matrix_g)

# Print the result
print(missing_prey)
print(missing_prey_g)

# Add new prey columns with all zeros for prey in Malaise or pitfall but not diet
new_prey <- c("Anthomyiidae", "Calliophoridae", "Mycetophilidae",
              "Phoridae", "Scathophagidae", "Syrphidae", "Encyrtidae",
              "Geometridae")  # names of prey not in diet
matrix_june_flying_expanded <- cbind(matrix_june_flying_transposed,
                                     matrix(0, nrow = nrow(matrix_june_flying_transposed),
                                            ncol = length(new_prey),
                                            dimnames = list(rownames(matrix_june_flying_transposed), new_prey)))

matrix_june_flying_expanded_t <- t(matrix_june_flying_expanded)

# Get the column names of the expanded matrix (prey)
prey_order <- rownames(matrix_june_flying_expanded_t)

# Reorder flying_prey_june_matrix to match
flying_prey_june_ordered <- flying_prey_june_matrix[, prey_order]

# Add new prey columns with all zeros for prey in Malaise or pitfall but not diet
new_prey_g <- c("Linyphiidae", "Coccinellidae", "Thysanoptera")  # names of prey not in diet
matrix_june_ground_expanded <- cbind(matrix_june_ground_transposed,
                                     matrix(0, nrow = nrow(matrix_june_ground_transposed),
                                            ncol = length(new_prey_g),
                                            dimnames = list(rownames(matrix_june_ground_transposed), new_prey_g)))

matrix_june_ground_expanded_t <- t(matrix_june_ground_expanded)

# Get the column names of the expanded matrix (prey)
prey_order_g <- rownames(matrix_june_ground_expanded_t)

# Reorder flying_prey_june_matrix to match
ground_prey_june_ordered <- ground_prey_june_matrix[, prey_order_g]


plotweb(
  matrix_june_flying_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = flying_prey_june_ordered
)

plotweb(
  matrix_june_ground_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = ground_prey_june_ordered
)

###JULY
flying_prey_july <- flyingprey %>%
  filter(Month == "July")

ground_prey_july <- groundprey %>%
  filter(Month == "7")

#Drop spiders from malaise samples
flying_prey_july <- flying_prey_july %>%
  filter(!Taxon %in% c("Araneidae", "Phalangiidae"))

flying_prey_july_matrix <- flying_prey_july %>%
  select(Taxon, Total) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Total,
    values_fill = 0
  ) %>%
  as.matrix()

#Drop flying arthropods from pitfall traps
ground_prey_july <- ground_prey_july %>%
  filter(!Family %in% c("Scelionidae", "Ichneumonidae",
                        "Syrphidae", "Sciaridae", "Phoridae",
                        "Muscidae", "Chironomidae", "Anthomyiidae",
                        "Psylloidea", "Coccoidea", "Aphidoidea",
                        "Calliphoridae", "Cecidomyiidae", "Chalcidoidea",
                        "Dolichopodidae", "Mycetophilidae",
                        "Syrphidae", "Simuliidae", "Tachinidae"))


ground_prey_july_matrix <- ground_prey_july %>%
  select(Family, Count) %>%
  pivot_wider(
    names_from = Family,
    values_from = Count,
    values_fill = 0
  ) %>%
  as.matrix()

#2. Split diet into ground dwelling prey and flying prey
edgelist_july_flying <- edgelist_july_Fam %>%
  filter(Prey %in% c("Anthomyiidae", "Braconidae", "Chironomidae", 
                     "Ichneumonidae", "Tipulidae", "Tortricidae",
                     "Noctuidae", "Vespidae", "Culicidae",
                     "Empididae", "Geometridae", 
                     "Scathophagidae", "Sciaridae", "Tachinidae",
                     "Syrphidae", "Simuliidae"))

edgelist_july_ground <- edgelist_july_Fam %>%
  filter(!Prey %in% c("Anthomyiidae", "Braconidae", "Chironomidae", 
                      "Ichneumonidae", "Tipulidae", "Tortricidae",
                      "Noctuidae", "Vespidae", "Culicidae",
                      "Empididae", "Geometridae", 
                      "Scathophagidae", "Sciaridae", "Tachinidae",
                      "Syrphidae", "Simuliidae"))

matrix_july_flying <- edgelist_to_matrix(edgelist_july_flying)
matrix_july_ground <- edgelist_to_matrix(edgelist_july_ground)

# Transpose your matrix so predators are rows and prey are columns
matrix_july_flying_transposed <- t(matrix_july_flying)
matrix_july_ground_transposed <- t(matrix_july_ground)

# Subset flying_prey_june_matrix to only include prey in matrix_june_flying
flying_prey_july_matrix2 <- flying_prey_july_matrix[, colnames(matrix_july_flying_transposed)]
ground_prey_july_matrix2 <- ground_prey_july_matrix[, colnames(matrix_july_ground_transposed)]

flying_prey_july_matrix2 <- matrix(flying_prey_july_matrix2, nrow = 1)
colnames(flying_prey_july_matrix2) <- colnames(matrix_july_flying_transposed)

ground_prey_july_matrix2 <- matrix(ground_prey_july_matrix2, nrow = 1)
colnames(ground_prey_july_matrix2) <- colnames(matrix_july_ground_transposed)

# Get prey names in each matrix
prey_in_matrix_july <- colnames(flying_prey_july_matrix2)
prey_in_abundance_july <- colnames(flying_prey_july_matrix)

prey_in_matrix_july_g <- colnames(ground_prey_july_matrix2)
prey_in_abundance_july_g <- colnames(ground_prey_july_matrix)


# Find prey in abundance but not in matrix
missing_prey_july <- setdiff(prey_in_abundance_july, prey_in_matrix_july)
missing_prey_july_g <- setdiff(prey_in_abundance_july_g, prey_in_matrix_july_g)

# Print the result
print(missing_prey_july)
print(missing_prey_july_g)

# Add new prey columns with all zeros for prey in Malaise or pitfall but not diet
new_prey_july <- c("Agromyzidae", "Ceratopogonidae", "Dolichopodidae",
              "Hypodermatidae", "Muscidae", "Mycetophilidae", "Encyrtidae",
              "Psyllidae", "Braconidae", "Diapriidae", "Figitidae",
              "Hemerobiidae", "Limnephilidae")  # names of prey not in diet
matrix_july_flying_expanded <- cbind(matrix_july_flying_transposed,
                                     matrix(0, nrow = nrow(matrix_july_flying_transposed),
                                            ncol = length(new_prey_july),
                                            dimnames = list(rownames(matrix_july_flying_transposed), new_prey_july)))

matrix_july_flying_expanded_t <- t(matrix_july_flying_expanded)

# Get the column names of the expanded matrix (prey)
prey_order_july <- rownames(matrix_july_flying_expanded_t)

# Reorder flying_prey_june_matrix to match
flying_prey_july_ordered <- flying_prey_july_matrix[, prey_order_july]

# Add new prey columns with all zeros for prey in Malaise or pitfall but not diet
new_prey_july_g <- c("Lygaeidae", "Coccinellidae",
                     "Thysanoptera", "Cicadellidae")  # names of prey not in diet
matrix_july_ground_expanded <- cbind(matrix_july_ground_transposed,
                                     matrix(0, nrow = nrow(matrix_july_ground_transposed),
                                            ncol = length(new_prey_july_g),
                                            dimnames = list(rownames(matrix_july_ground_transposed), new_prey_july_g)))

matrix_july_ground_expanded_t <- t(matrix_july_ground_expanded)

# Get the column names of the expanded matrix (prey)
prey_order_july_g <- rownames(matrix_july_ground_expanded_t)

# Reorder flying_prey_june_matrix to match
ground_prey_july_ordered <- ground_prey_july_matrix[, prey_order_july_g]


plotweb(
  matrix_july_flying_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = flying_prey_july_ordered
)

plotweb(
  matrix_july_ground_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = ground_prey_july_ordered
)

###JULY
flying_prey_aug <- flyingprey %>%
  filter(Month == "August")

ground_prey_aug <- groundprey %>%
  filter(Month == "8")

#Drop spiders from malaise samples
flying_prey_aug <- flying_prey_aug %>%
  filter(!Taxon %in% c("Araneidae", "Phalangiidae", "Curculionidae",
                       "Lygaeidae"))

flying_prey_aug_matrix <- flying_prey_aug %>%
  select(Taxon, Total) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = Total,
    values_fill = 0
  ) %>%
  as.matrix()

#Drop flying arthropods from pitfall traps
ground_prey_aug <- ground_prey_aug %>%
  filter(!Family %in% c("Scelionidae", "Ichneumonidae",
                        "Syrphidae", "Sciaridae", "Phoridae",
                        "Muscidae", "Chironomidae", "Anthomyiidae",
                        "Psylloidea", "Coccoidea", "Aphidoidea",
                        "Calliphoridae", "Cecidomyiidae", "Chalcidoidea",
                        "Dolichopodidae", "Mycetophilidae",
                        "Syrphidae", "Simuliidae", "Tachinidae",
                        "Cynipoidea", "Scatophagidae", "Vespidae"))


ground_prey_aug_matrix <- ground_prey_aug %>%
  select(Family, Count) %>%
  pivot_wider(
    names_from = Family,
    values_from = Count,
    values_fill = 0
  ) %>%
  as.matrix()

#2. Split diet into ground dwelling prey and flying prey
edgelist_august_Fam <- edgelist_august_Fam %>%
  filter(!Prey %in% "Corvidae")

edgelist_aug_flying <- edgelist_august_Fam %>%
  filter(Prey %in% c("Chironomidae", 
                     "Ichneumonidae", "Empididae", 
                     "Scathophagidae","Simuliidae"))

edgelist_aug_ground <- edgelist_august_Fam %>%
  filter(!Prey %in% c("Chironomidae", 
                      "Ichneumonidae", "Empididae", 
                      "Scathophagidae","Simuliidae"))

matrix_aug_flying <- edgelist_to_matrix(edgelist_aug_flying)
matrix_aug_ground <- edgelist_to_matrix(edgelist_aug_ground)

# Transpose your matrix so predators are rows and prey are columns
matrix_aug_flying_transposed <- t(matrix_aug_flying)
matrix_aug_ground_transposed <- t(matrix_aug_ground)

# Subset flying_prey_june_matrix to only include prey in matrix_june_flying
flying_prey_aug_matrix2 <- flying_prey_aug_matrix[, colnames(matrix_aug_flying_transposed)]
ground_prey_aug_matrix2 <- ground_prey_aug_matrix[, colnames(matrix_aug_ground_transposed)]

flying_prey_aug_matrix2 <- matrix(flying_prey_aug_matrix2, nrow = 1)
colnames(flying_prey_aug_matrix2) <- colnames(matrix_aug_flying_transposed)

ground_prey_aug_matrix2 <- matrix(ground_prey_aug_matrix2, nrow = 1)
colnames(ground_prey_aug_matrix2) <- colnames(matrix_aug_ground_transposed)

# Get prey names in each matrix
prey_in_matrix_aug <- colnames(flying_prey_aug_matrix2)
prey_in_abundance_aug <- colnames(flying_prey_aug_matrix)

prey_in_matrix_aug_g <- colnames(ground_prey_aug_matrix2)
prey_in_abundance_aug_g <- colnames(ground_prey_aug_matrix)


# Find prey in abundance but not in matrix
missing_prey_aug <- setdiff(prey_in_abundance_aug, prey_in_matrix_aug)
missing_prey_aug_g <- setdiff(prey_in_abundance_aug_g, prey_in_matrix_aug_g)

# Print the result
print(missing_prey_aug)
print(missing_prey_aug_g)

# Add new prey columns with all zeros for prey in Malaise or pitfall but not diet
new_prey_aug <- c("Agromyzidae", "Anthomyiidae",
                  "Calliphoridae", "Ceratopogonidae", 
                  "Culicidae", "Dolichopodidae",
                   "Hypodermatidae", "Muscidae", "Mycetophilidae", "Phoridae",
                  "Sciaridae", "Syrphidae", "Tachinidae",
                  "Aphididae", "Psyllidae", "Braconidae",
                  "Diapriidae", "Encyrtidae", "Eulophidae",
                  "Pteromalidae", "Vespidae", "Geometridae",
                  "Microlepidoptera", "Hemerobiidae", "Limnephilidae")  # names of prey not in diet
matrix_aug_flying_expanded <- cbind(matrix_aug_flying_transposed,
                                     matrix(0, nrow = nrow(matrix_aug_flying_transposed),
                                            ncol = length(new_prey_aug),
                                            dimnames = list(rownames(matrix_aug_flying_transposed),
                                                            new_prey_aug)))

matrix_aug_flying_expanded_t <- t(matrix_aug_flying_expanded)

# Get the column names of the expanded matrix (prey)
prey_order_aug <- rownames(matrix_aug_flying_expanded_t)

# Reorder flying_prey_june_matrix to match
flying_prey_aug_ordered <- flying_prey_aug_matrix[, prey_order_aug]

# Add new prey columns with all zeros for prey in Malaise or pitfall but not diet
new_prey_aug_g <- c("Gnaphosidae", "Hahniidae",
                     "Linyphiidae", "Thomisidae",
                    "Curculionidae", "Phalangiidae",
                    "Thysanoptera")  # names of prey not in diet
matrix_aug_ground_expanded <- cbind(matrix_aug_ground_transposed,
                                     matrix(0, nrow = nrow(matrix_aug_ground_transposed),
                                            ncol = length(new_prey_aug_g),
                                            dimnames = list(rownames(matrix_aug_ground_transposed), 
                                                            new_prey_aug_g)))

matrix_aug_ground_expanded_t <- t(matrix_aug_ground_expanded)

# Get the column names of the expanded matrix (prey)
prey_order_aug_g <- rownames(matrix_aug_ground_expanded_t)

# Reorder flying_prey_june_matrix to match
ground_prey_aug_ordered <- ground_prey_aug_matrix[, prey_order_aug_g]


plotweb(
  matrix_aug_flying_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = flying_prey_aug_ordered
)

plotweb(
  matrix_aug_ground_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = ground_prey_aug_ordered
)

#Divide into birds, dungflies, and spiders.



#We should take prey eating prey into account, but 
#not enough data for that yet I would say
#Ignore it for now and calculate some basic metrics

#Degree

#Dependence

#Species strength
strength(network_matrix, type = "Bascompte")

#Modularity



#####################
#Piechart data

june_df <- read.csv("edgelist_june_2024_Fam.csv")

june_df$Pred_Group <- fct_collapse(
  june_df$Predator,
  Birds = c("Lapland bunting", "Snow bunting", "Wheatear"),
  GroundSpiders = c("Phalangiidae", "Lycosidae"),
  Orbweavers = c("Araneidae", "Thomisidae"),
  Divingbeetle = "Dytiscidae",
  Midges = "Chironomidae"
)

june_df <- june_df %>%
  group_by(Pred_Group, Prey) %>%
  summarise(
    TotalInteraction = sum(InteractionStrength),
    .groups = "drop"
  )

june_df$Order <- fct_collapse(
  june_df$Prey,
  Aranea = c("Araneidae", "Lycosidae", "Philodromidae", "Thomisidae", "Gnaphosidae"),
  Opiliones = "Phalangiidae",
  Coleoptera = c("Dytiscidae", "Carabidae", "Curculionidae", "Byrrhidae"),
  Diptera = c("Chironomidae", "Culicidae"),
  Hymenoptera = c("Ichneumonidae", "Braconidae", "Vespidae"),
  Lepidoptera = "Noctuidae"
)

july_df <- read.csv("edgelist_july_2024_Fam.csv")

july_df <- july_df %>%
  filter(!is.na(Prey))

july_df$Pred_Group <- fct_collapse(
  july_df$Predator,
  Birds = c("Lapland bunting", "Snow bunting", "Wheatear"),
  Groundspiders = c("Phalangiidae", "Lycosidae", "Philodromidae"),
  Orbweavers = c("Araneidae", "Thomisidae"),
  Divingbeetle = "Dytiscidae",
  Midges = "Chironomidae"
)

july_df <- july_df %>%
  group_by(Pred_Group, Prey) %>%
  summarise(
    TotalInteraction = sum(InteractionStrength),
    .groups = "drop"
  )

july_df$Order <- fct_collapse(
  july_df$Prey,
  Aranea = c("Araneidae", "Lycosidae", "Linyphiidae",
             "Philodromidae", "Thomisidae", "Gnaphosidae"),
  Opiliones = "Phalangiidae",
  Coleoptera = c("Dytiscidae", "Carabidae", "Curculionidae", "Byrrhidae"),
  Diptera = c("Scathophagidae", "Empididae", "Simuliidae",
              "Anthomyiidae", "Chironomidae", "Culicidae",
              "Sciaridae", "Syrphidae", "Tachinidae", "Tipulidae"),
  Hymenoptera = "Ichneumonidae",
  Lepidoptera = c("Geometridae", "Noctuidae", "Tortricidae"),
  Sarcoptiformes = "Oribatulidae"
)

aug_df <- read.csv("edgelist_august_2024_Fam.csv")

aug_df <- aug_df %>%
  filter(!Prey %in% "Fringillidae") %>%
  filter(!Prey %in% "Corvidae")

aug_df$Pred_Group <- fct_collapse(
  aug_df$Predator,
  Birds = "Wheatear",
  Groundspiders = c("Lycosidae", "Philodromidae"),
  Orbweavers = "Araneidae",
  Divingbeetle = "Dytiscidae",
  Predatoryflies = "Scathophagidae"
)

aug_df <- aug_df %>%
  group_by(Pred_Group, Prey) %>%
  summarise(
    TotalInteraction = sum(InteractionStrength),
    .groups = "drop"
  )

aug_df$Order <- fct_collapse(
  aug_df$Prey,
 Aranea = c("Araneidae", "Lycosidae"),
 Coleoptera = "Dytiscidae",
 Diptera = c("Scathophagidae", "Empididae", "Simuliidae",
             "Chironomidae"),
 Hymenoptera = "Ichneumonidae"
)

writexl::write_xlsx(june_df, "Piechart_june24.xlsx")
writexl::write_xlsx(july_df, "Piechart_july24.xlsx")
writexl::write_xlsx(aug_df, "Piechart_aug24.xlsx")


######
#BLOODMEAL SAMPLES
#####
#Include only bloodmeal primers
included_primers <- c("P16S", "12S-V5")

Bloodmeal_filter <- Insect_ID %>%
  filter(PrimerPair == included_primers)

Bloodmeal_filter <- Bloodmeal_filter %>%
  filter(!Order %in% "NA")

#Only include bloodmeal samples from diel matrix
diet_matrix_blood <- diet_matrix[, colnames(diet_matrix) %in% Bloodmeal_filter$Sample]

##Add prey names to the rows
prey_names_blood <- prey_data %>%
  pull(Order)

#Add prey names as a column
diet_long_blood <- diet_matrix_blood %>%
  as.data.frame() %>%         # ensure it's a dataframe
  mutate(Prey = prey_names_blood) %>%
  pivot_longer(
    cols = -Prey,             # keep Prey column
    names_to = "PredatorSample",
    values_to = "Presence"
  )

#Keep only bloodmeal prey
diet_long_blood <- diet_long_blood %>%
  filter(Prey %in% c("Passeriformes", "Artiodactyla",
                   "Carnivora", "Primates"))

#Make all samples Aedes
diet_long_blood <- diet_long_blood %>%
  mutate(PredatorSample = "Aedes")

#Sum presence per Predator-Prey combination
interaction_counts <- diet_long_blood %>%
  group_by(PredatorSample, Prey) %>%
  summarise(TotalPresence = sum(Presence), .groups = "drop")

# 3. Calculate total number of samples per predator
total_samples <- diet_long_blood %>%
  group_by(PredatorSample) %>%
  summarise(TotalSamples = n(), .groups = "drop")

# 4. Join totals and calculate interaction strength
edgelist_blood <- interaction_counts %>%
  left_join(total_samples, by = "PredatorSample") %>%
  mutate(InteractionStrength = TotalPresence / TotalSamples) %>%
  select(PredatorSample, Prey, InteractionStrength)

# Convert to matrix format
network_blood <- edgelist_blood %>%
  pivot_wider(
    names_from = PredatorSample,
    values_from = InteractionStrength,
    values_fill = 0
  ) %>%
  column_to_rownames("Prey") %>%  # rows = prey
  as.matrix()

rownames(network_blood)
colnames(network_blood)

plotweb(network_blood,
        text_size = 0.5)

png("bloodmeal_web.png", width = 3000, height = 1500, res = 300)

plotweb(
  network_blood,
  text_size = 0.5
)

dev.off()

###################

#Tri-trophic network visualization
library(igraph)
library(tidygraph)
library(ggraph)

edgelist_clean <- read.csv("edgelist_2024_Fam.csv")

edgelist_clean <- edgelist_clean[,2:4]


#Create basic graph
#Add the interactions
g <- graph_from_data_frame(
  edgelist_clean,
  directed = TRUE,
  vertices = unique(c(edgelist_clean$Predator, edgelist_clean$Prey))
)

#Add interaction strength as weights
E(g)$weight <- edgelist_clean$InteractionStrength

#Check if it worked
g
is_weighted(g)

#Weighted degree metrics
strength_in  <- strength(g, mode = "in",  weights = E(g)$weight)
strength_out <- strength(g, mode = "out", weights = E(g)$weight)

#Add trophic level to each predator and prey
A <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)

# Normalize rows to diet proportions
P <- A / rowSums(A)
P[is.na(P)] <- 0

n <- nrow(P)
I <- diag(n)

TL <- solve(I - P, rep(1, n))
names(TL) <- V(g)$name

V(g)$trophic_level <- TL

summary(V(g)$trophic_level)

#Define three trophic levels
V(g)$trophic_class <- cut(
  V(g)$trophic_level,
  breaks = c(-Inf, 1.5, 3.5, Inf),
  labels = c("Basal", "Intermediate", "Top")
)



table(V(g)$trophic_class)
tapply(V(g)$trophic_level, V(g)$trophic_class, range)

V(g)$rank <- as.numeric(V(g)$trophic_class)
V(g)$rank

#Plot using tidygraph and ggraph
lay <- layout_with_sugiyama(
  g,
  layers = V(g)$rank,
  reorder = TRUE,  # reorders nodes horizontally for fewer crossings
  hgap = 2,   # horizontal gap between nodes
  vgap = 3    # vertical gap between layers
)

V(g)$trophic_level[V(g)$name %in% c("Snow bunting", "Lapland bunting",
                                    "Wheatear")] <- max(V(g)$trophic_level) + 0.1

tg <- as_tbl_graph(g)

ggraph(tg, layout = "manual",
       x = lay$layout[,1],
       y = -lay$layout[,2]) +   # top predators are at the top
  geom_edge_link(aes(width = weight), alpha = 0.6) +
  geom_node_text(aes(label = name, color = trophic_class), size = 3, repel = TRUE) +
  scale_edge_width(range = c(0.3, 2)) +
  theme_graph(base_family = "Arial")

species_trophic <- data.frame(
  Species = V(g)$name,
  TrophicLevel = V(g)$trophic_level
)

species_trophic


# Show TL and prey of a species
data.frame(
  Species = V(g)$name,
  TrophicLevel = V(g)$trophic_level,
  PreyCount = degree(g, mode = "out"),
  PredatorCount = degree(g, mode = "in")
)

