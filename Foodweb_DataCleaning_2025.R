library(readxl)
library(readr)
library(tidyverse)

rm(list = ls())

set.seed(123)

matrix <- read_xlsx("mergedData_wide_ReadCounts_2025.xlsx")

levels(as.factor(matrix$Phylum))

#Remove the plants and worms for now
matrix <- matrix |> 
  filter(!Phylum %in% c("Streptophyta", "Annelida"))

diet_matrix <- matrix[, 8:1018]

prey_data <- matrix[, 1:7]

diet_matrix <- ifelse(diet_matrix > 0, 1, 0)

Insect_ID <- read.csv("mergedData_long_format_2025.csv", sep = ";")

Base_ID <- read.csv("FoodwebSamples_2025.csv", sep = ";")

Base_ID <- Base_ID %>%
  mutate(
    Date = str_pad(Date, width = 6, side = "left", pad = "0")  # make all dates 6 digits
  ) |> 
  filter(!Group %in% "")

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
  filter(Group == "Birds") %>%
  select(ID, Species, Month) %>%
  rename(
    Sample = ID,
    Family   = Species,
    Month = Month
  )

levels(as.factor(Base_ID$Group))

NonSpiderPred_ID <- Base_ID %>%
  filter(Group %in% c("Predatory_flies", "Aquatic_predators")) %>%
  select(ID, Group, Species, Month) %>%
  rename(
    Sample = ID,
    Group2   = Species,
    Month = Month
  )

levels(as.factor(NonSpiderPred_ID$Group2))

NonSpiderPred_ID$Family <- fct_collapse(
  NonSpiderPred_ID$Group2,
  Dytiscidae = "Colymbetes_dolabratus",
  Gyrinidae = "Gyrinus_opacus",
  Scathophagidae = c("Scatophagidae", "Unknown", ""))

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

levels(as.factor(Insect_ID$PrimerPair))

#Remove the bloodmeal and plant samples for now
excluded_primers <- c("ChoP16S", "12S-V5", "rbcLa",
                      "trnLch")

Insect_Filter <- Insect_ID %>%
  filter(!PrimerPair %in% excluded_primers)

levels(as.factor(Insect_Filter$PrimerPair))

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

levels(as.factor(dominant_spiders$Family))

#Identify families that cannot be spider predators (all non-spiders)
non_predator_species <- c(
  "Braconidae",
  "Carabidae",
  "Chironomidae",
  "Coccinellidae",
  "Culicidae",
  "Dytiscidae",
  "Geometridae",
  "Gyrinidae",
  "Ichneumonidae",
  "Muscidae",
  "Noctuidae",
  "Scathophagidae",
  "Sciaridae",
  "Syrphidae",
  "Tortricidae"
)

dominant_spider_sp <- dominant_spiders %>%
  filter(!Family %in% non_predator_species) %>% 
  group_by(Sample) %>%
  slice_max(order_by = total_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Sample, Family, Month)

#Sample 285 is lacking a family - should be Lycosidae
dominant_spider_sp[dominant_spider_sp$Sample == "KF25-285", "Family"] <- "Lycosidae"

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

write_csv(as.data.frame(diet_matrix_clean), "diet_matrix_Fam25.csv")
write_csv(as.data.frame(diet_june_Fam), "diet_matrix_june_Fam25.csv")
write_csv(as.data.frame(diet_july_Fam), "diet_matrix_july_Fam25.csv")
write_csv(as.data.frame(diet_august_Fam), "diet_matrix_august_Fam25.csv")

diet_matrix_clean <- read.csv("diet_Matrix_Fam25.csv")

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

levels(as.factor(edgelist$Prey))

edgelist <- edgelist %>%
  filter(!Prey %in% "NA")

#Remove unlikely families
edgelist <- edgelist %>%
  filter(!Prey %in% c("Blattidae", "Formicidae", "Elateridae", "Cantharidae",
                      "Lasiocampidae"))

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

write.csv(edgelist_clean, "edgelist_2025_Fam.csv")

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


write.csv(edgelist_june_Fam, "edgelist_june_2025_Fam.csv")
write.csv(edgelist_july_Fam, "edgelist_july_2025_Fam.csv")
write.csv(edgelist_august_Fam, "edgelist_august_2025_Fam.csv")


