





flying_prey_june_matrix <- flying_prey_june %>%
  select(Family, Count) %>%
  pivot_wider(
    names_from = Family,
    values_from = Count,
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

matrix_june <- edgelist_to_matrix(edgelist_june_Fam)
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

#Add color to Culicidae
lower_color <- rep("grey20", nrow(matrix_june_flying_expanded_t))
names(lower_color) <- rownames(matrix_june_flying_expanded_t)
lower_color["Culicidae"] <- "goldenrod"

#Add sample size of predators
sample_size <- read.csv("SampleSizePred2024.csv", sep = ";")
colnames(sample_size) <- gsub("\\.", " ", colnames(sample_size))

sample_size_june <- sample_size %>%
  filter(Month == "June")
sample_size_june <- sample_size_june[,2:13]

sample_size_july <- sample_size %>%
  filter(Month == "July")
sample_size_july <- sample_size_july[,2:13]

sample_size_aug <- sample_size %>%
  filter(Month == "August")
sample_size_aug <- sample_size_aug[,2:13]

predator_order_june <- colnames(matrix_june_flying_expanded_t)
sample_size_june_ord <- sample_size_june[, predator_order_june, drop = FALSE]

#Extract row values as a vector
higher_abundances <- as.numeric(sample_size_june_ord[1, ])
names(higher_abundances) <- colnames(sample_size_june_ord)

#Add color to Culicidae
lower_color <- rep("grey20", nrow(matrix_june_flying_expanded_t))
names(lower_color) <- rownames(matrix_june_flying_expanded_t)
lower_color["Culicidae"] <- "goldenrod"

png("foodweb_june_24.png", width = 3000, height = 1500, res = 300)

plotweb(
  matrix_june,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  #  curved_links = TRUE,
  #  lower_color = lower_color,
  #  link_color = "lower"
)
dev.off()


png("foodweb_june_flying_2024.png", width = 3000, height = 1000, res = 300)

plotweb(
  matrix_june_flying_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = flying_prey_june_ordered,
  higher_abundances = higher_abundances,
  lower_color = lower_color,
  link_color = "lower"
)

dev.off()

predator_order_june_g <- colnames(matrix_june_ground_expanded_t)
sample_size_june_ord_g <- sample_size_june[, predator_order_june_g, drop = FALSE]


higher_abundances_june_g <- as.numeric(sample_size_june_ord_g[1, ])
names(higher_abundances_june_g) <- colnames(sample_size_june_ord_g)

png("foodweb_june_ground_2024.png", width = 3000, height = 1250, res = 300)

plotweb(
  matrix_june_ground_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  higher_abundances = higher_abundances_june_g,
  lower_abundances = ground_prey_june_ordered
)

dev.off()



#Drop spiders from malaise samples
flying_prey_july <- flying_prey_july %>%
  filter(!Family %in% c("Araneidae", "Phalangiidae"))

flying_prey_july_matrix <- flying_prey_july %>%
  select(Family, Count) %>%
  pivot_wider(
    names_from = Family,
    values_from = Count,
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

matrix_july <- edgelist_to_matrix(edgelist_july_Fam)
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

predator_order_july <- colnames(matrix_july_flying_expanded_t)
sample_size_july_ord <- sample_size_july[, predator_order_july, drop = FALSE]

#Extract row values as a vector
higher_abundances <- as.numeric(sample_size_july_ord[1, ])
names(higher_abundances) <- colnames(sample_size_july_ord)


png("foodweb_july_flying_2024.png", width = 3000, height = 1250, res = 300)

plotweb(
  matrix_july_flying_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = flying_prey_july_ordered,
  higher_abundances = higher_abundances,
  lower_color = lower_color_jul,
  link_color = "lower"
)

dev.off()


predator_order_july_g <- colnames(matrix_july_ground_expanded_t)
sample_size_july_ord_g <- sample_size_july[, predator_order_july_g, drop = FALSE]


higher_abundances_july_g <- as.numeric(sample_size_july_ord_g[1, ])
names(higher_abundances_july_g) <- colnames(sample_size_july_ord_g)


png("foodweb_july_ground_2024.png", width = 3000, height = 1250, res = 300)

plotweb(
  matrix_july_ground_expanded_t,
  srt = 45,
  text_size = 0.5,
  sorting = "dec",
  curved_links = TRUE,
  lower_abundances = ground_prey_july_ordered,
  higher_abundances = higher_abundances_july_g
)

dev.off()
