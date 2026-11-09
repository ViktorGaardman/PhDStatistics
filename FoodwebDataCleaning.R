library(readxl)
library(readr)
library(tidyverse)

rm(list = ls())

set.seed(123)

matrix <- read_xlsx("mergedData_wide_ReadCounts_2024.xlsx")

#Keep only arthropods
levels(as.factor(matrix$Phylum))

matrix <- matrix |> 
  filter(!Phylum == "Chordata")

diet_matrix <- matrix[, 9:489]

prey_data <- matrix[, 1:8]

diet_matrix <- diet_matrix %>%
  mutate(Family = prey_data$Family) %>%
  filter(!is.na(Family)) %>%
  group_by(Family) %>%
  summarise(across(everything(), max),
            .groups = "drop")

Insect_ID <- read.csv("mergedData_long_format_no_zotu_2024.csv", sep = ";")

Base_ID <- read.csv("Foodweb_samples_2024.csv", sep = ";")

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
  filter(Group == "Bird") %>%
  select(ID, Species, Month) %>%
  rename(
    Sample = ID,
    Family   = Species,
    Month = Month
  )

levels(as.factor(Base_ID$Group))

NonSpiderPred_ID <- Base_ID %>%
  filter(Group %in% c("Scatophagidae", "Aquatic_predators")) %>%
  select(ID, Group, Species, Month) %>%
  rename(
    Sample = ID,
    Group2   = Species,
    Month = Month
  )

levels(as.factor(NonSpiderPred_ID$Group2))

NonSpiderPred_ID$Family <- fct_collapse(
  NonSpiderPred_ID$Group2,
  Dytiscidae = c("Colymbetes_dolabratus", "C._dolabratus",
  "C._Dolabratus", "C._dolabratus_larvae"),
  Chironomidae = "Chironomid_larvae",
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

levels(as.factor(Insect_ID$PrimerPair))

#Remove the bloodmeal and plant samples for now
excluded_primers <- c("P16S", "12S-V5")

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
  "Chironomidae",
  "Culicidae",
  "Dytiscidae",
  "Geometridae",
  "Scathophagidae",
  "Simuliidae"
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

# Store prey family names separately
prey_names <- diet_matrix$Family

# Remove family column, leaving only sample columns
diet_matrix_clean <- diet_matrix %>%
  select(-Family)

# Remove samples not identified as predators
diet_matrix_clean <- diet_matrix_clean[,
                                       colnames(diet_matrix_clean) %in% predator_lookup$Sample
]

# Store original column names
original_colnames <- colnames(diet_matrix_clean)

# Replace sample IDs with predator families
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

#Store original column names
original_colnames <- colnames(diet_matrix_clean)

# Make them unique internally
colnames(diet_matrix_clean) <- make.unique(original_colnames, sep = "_sample_")

#Turn read counts into pres abs
diet_matrix_clean <- ifelse(diet_matrix_clean > 0, 1, 0) 

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

#remove predator DNA
diet_long_clean <- diet_long %>%
  filter(Prey != Predator) %>%
  filter(!Prey %in% "NA") |> 
  filter(!Prey %in% "Fringillidae")

edgelist <- diet_long_clean %>%
  group_by(Predator, Prey) %>%
  summarise(
    TotalPresence = sum(Presence),
    .groups = "drop"
  ) %>%
  left_join(
    diet_long %>%
      distinct(Predator, PredatorSample) %>%
      count(Predator, name = "N"),
    by = "Predator"
  ) %>%
  mutate(
    FOO = TotalPresence / N
  )

write.csv(edgelist, "edgelist_2024_Fam.csv")

#Monthly edgeslists
diet_matrix_to_edgelist <- function(diet_matrix, prey_names) {
  #Store original column names
  original_colnames <- colnames(diet_matrix)
  
  # Make them unique internally
  colnames(diet_matrix) <- make.unique(original_colnames, sep = "_sample_")
  
  #Turn read counts into pres abs
  diet_matrix <- ifelse(diet_matrix > 0, 1, 0) 
  
  #Add prey names as a column
  diet_long <- diet_matrix %>%
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
  
  #remove predator DNA
  diet_long_clean <- diet_long %>%
    filter(Prey != Predator) %>%
    filter(!Prey %in% "NA")
  
  edgelist <- diet_long_clean %>%
    group_by(Predator, Prey) %>%
    summarise(
      TotalPresence = sum(Presence),
      .groups = "drop"
    ) %>%
    left_join(
      diet_long %>%
        distinct(Predator, PredatorSample) %>%
        count(Predator, name = "N"),
      by = "Predator"
    ) %>%
    mutate(
      FOO = TotalPresence / N
    )
  return(edgelist)
}

prey_names <- diet_matrix$Family

edgelist_june_Fam   <- diet_matrix_to_edgelist(diet_june_Fam, prey_names)
edgelist_july_Fam   <- diet_matrix_to_edgelist(diet_july_Fam, prey_names)
edgelist_august_Fam <- diet_matrix_to_edgelist(diet_august_Fam, prey_names)


write.csv(edgelist_june_Fam, "edgelist_june_2024_Fam.csv")
write.csv(edgelist_july_Fam, "edgelist_july_2024_Fam.csv")
write.csv(edgelist_august_Fam, "edgelist_august_2024_Fam.csv")



#RRA VERSION

rm(list = ls())

set.seed(123)

matrix <- read_xlsx("mergedData_wide_ReadCounts_2024.xlsx")

#Keep only arthropods
levels(as.factor(matrix$Phylum))

matrix <- matrix |> 
  filter(!Phylum == "Chordata")

diet_matrix_rel <- matrix[, 9:489]

prey_data <- matrix[, 1:8]

#Aggregate dataset to family level
diet_matrix_rel <- diet_matrix_rel %>%
  mutate(Family = prey_data$Family) %>%
  filter(!is.na(Family)) %>%
  group_by(Family) %>%
  summarise(across(everything(), sum),
            .groups = "drop")

Insect_ID <- read.csv("mergedData_long_format_no_zotu_2024.csv", sep = ";")

Base_ID <- read.csv("Foodweb_samples_2024.csv", sep = ";")

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
  filter(Group == "Bird") %>%
  select(ID, Species, Month) %>%
  rename(
    Sample = ID,
    Family   = Species,
    Month = Month
  )

levels(as.factor(Base_ID$Group))

NonSpiderPred_ID <- Base_ID %>%
  filter(Group %in% c("Scatophagidae", "Aquatic_predators")) %>%
  select(ID, Group, Species, Month) %>%
  rename(
    Sample = ID,
    Group2   = Species,
    Month = Month
  )

levels(as.factor(NonSpiderPred_ID$Group2))

NonSpiderPred_ID$Family <- fct_collapse(
  NonSpiderPred_ID$Group2,
  Dytiscidae = c("Colymbetes_dolabratus", "C._dolabratus",
                 "C._Dolabratus", "C._dolabratus_larvae"),
  Chironomidae = "Chironomid_larvae",
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

levels(as.factor(Insect_ID$PrimerPair))

#Remove the bloodmeal and plant samples for now
excluded_primers <- c("P16S", "12S-V5")

Insect_Filter <- Insect_ID %>%
  filter(!PrimerPair %in% excluded_primers)

levels(as.factor(Insect_Filter$PrimerPair))

#Sum reads across remaining primers
summed_reads <- Insect_Filter %>%
  group_by(Sample, Family, Genus, Month) %>%
  summarise(total_reads = sum(Count),
            .groups = "drop")


#Remove non-spider arthropod predators from field data
dominant_spidersbirds <- summed_reads %>%
  filter(!Sample %in% NonSpiderPred_ID$Sample)

#Remove bird predators, leaving only spiders
dominant_spiders <- dominant_spidersbirds %>%
  filter(!Sample %in% Bird_ID$Sample)

levels(as.factor(dominant_spiders$Family))

#Identify families that cannot be spider predators (all non-spiders)
non_predator_species <- c(
  "Chironomidae",
  "Culicidae",
  "Dytiscidae",
  "Geometridae",
  "Scathophagidae",
  "Simuliidae"
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

make_monthly_diet_matrices <- function(diet_matrix_rel, predator_lookup) {
  
  # Helper function for a single month
  make_one_month <- function(month_name) {
    
    month_df <- predator_lookup %>%
      filter(Month == month_name) %>%
      select(Sample, Family)
    
    # Filter diet matrix
    mat <- diet_matrix_rel[, colnames(diet_matrix_rel) %in% month_df$Sample, drop = FALSE]
    
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
  diet_matrix_rel = diet_matrix_rel,
  predator_lookup = predator_lookup
)

diet_june_Fam   <- monthly_matrices$June
diet_july_Fam   <- monthly_matrices$July
diet_august_Fam <- monthly_matrices$August

# Store prey family names separately
prey_names <- diet_matrix_rel$Family

# Remove family column, leaving only sample columns
diet_matrix_clean <- diet_matrix_rel %>%
  select(-Family)

# Remove samples not identified as predators
diet_matrix_clean <- diet_matrix_clean[,
 colnames(diet_matrix_clean) %in% predator_lookup$Sample
]

# Store original column names
original_colnames <- colnames(diet_matrix_clean)

# Replace sample IDs with predator families
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

write_csv(diet_matrix_clean, "diet_matrix_Fam_RRA.csv")
write_csv(diet_june_Fam, "diet_matrix_june_Fam_RRA.csv")
write_csv(diet_july_Fam, "diet_matrix_july_Fam_RRA.csv")
write_csv(diet_august_Fam, "diet_matrix_august_Fam_RRA.csv")

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
    values_to = "TotalReads"
  )

#Add original names back
diet_long <- diet_long %>%
  mutate(Predator = sub("_sample_.*", "", PredatorSample))

#remove predator DNA
diet_long_clean <- diet_long %>%
  filter(Prey != Predator) %>%
  filter(!Prey %in% "NA") %>%
  filter(TotalReads > 0)

#Calculate RRA within each sample
diet_long_clean <- diet_long_clean %>%
  group_by(PredatorSample) %>%
  mutate(
    SampleRRA = TotalReads / sum(TotalReads, na.rm = TRUE)
  ) %>%
  ungroup()

n_pred_samples <- diet_long %>%
  distinct(Predator, PredatorSample) %>%
  count(Predator, name = "N")

edgelist <- diet_long_clean %>%
  group_by(Predator, Prey) %>%
  summarise(
    SumRRA = sum(SampleRRA),
    .groups = "drop"
  ) %>%
  left_join(n_pred_samples, by = "Predator") %>%
  mutate(
    MeanRRA = SumRRA / N
  )

write.csv(edgelist, "edgelist_2024_Fam_RRA.csv")

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
      values_to = "TotalReads"
    )
  
  # 4. Recover predator genus from column names
  diet_long <- diet_long %>%
    mutate(Predator = sub("_sample_.*", "", PredatorSample))
  
  #remove predator DNA
  diet_long_clean <- diet_long %>%
    filter(Prey != Predator) %>%
    filter(!Prey %in% "NA") %>%
    filter(TotalReads > 0)
  
  #Calculate RRA within each sample
  diet_long_clean <- diet_long_clean %>%
    group_by(PredatorSample) %>%
    mutate(
      SampleRRA = TotalReads / sum(TotalReads, na.rm = TRUE)
    ) %>%
    ungroup()
  
  n_pred_samples <- diet_long %>%
    distinct(Predator, PredatorSample) %>%
    count(Predator, name = "N")
  
  edgelist <- diet_long_clean %>%
    group_by(Predator, Prey) %>%
    summarise(
      SumRRA = sum(SampleRRA),
      .groups = "drop"
    ) %>%
    left_join(n_pred_samples, by = "Predator") %>%
    mutate(
      MeanRRA = SumRRA / N
    )
  
  return(edgelist)
}

prey_names <- diet_matrix_rel$Family

edgelist_june_Fam   <- diet_matrix_to_edgelist(diet_june_Fam, prey_names)
edgelist_july_Fam   <- diet_matrix_to_edgelist(diet_july_Fam, prey_names)
edgelist_august_Fam <- diet_matrix_to_edgelist(diet_august_Fam, prey_names)


write.csv(edgelist_june_Fam, "edgelist_june_2024_Fam_RRA.csv")
write.csv(edgelist_july_Fam, "edgelist_july_2024_Fam_RRA.csv")
write.csv(edgelist_august_Fam, "edgelist_august_2024_Fam_RRA.csv")
