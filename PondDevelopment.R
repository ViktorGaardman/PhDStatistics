library(tidyverse)
library(ggcorrplot)
library(car)
library(patchwork)
library(lme4)
library(performance)
library(DHARMa)
library(ggeffects)
library(nlme)

set.seed(123)

rm(list = ls())

setwd("C:/Users/vikto/OneDrive/Skrivbord/PhD/Fieldwork/Pond_Data")

#Load pond data

ponds2024 <- read.csv("Ponds2024_Clean.csv", 
                      sep=";", header=TRUE, as.is=TRUE)

ponds2025 <- read.csv("Ponds2025_Clean.csv", 
                      sep=";", header=TRUE, as.is = TRUE)

#Combine pond datasets

ponds2024$Temp_PE <- as.numeric(ponds2024$Temp_PE)

pond_data <- bind_rows(ponds2024, ponds2025)


#Standardize data
pond_data <- pond_data %>% 
  mutate(across(where(is.character) & c("Average_Volume",
                      "Conductivity", "PH", "Average_Temp", "Temp_PE",
                      "Arthropod_Diversity", "Plant_Diversity",
                      "Veg_Cov", "Perimeter"
  ), ~ as.numeric(gsub(",", ".", .))))
                                        

pond_data <- pond_data %>% 
  mutate(across(c("Mosses_PresAbs", "Eriophorum_PresAbs", "Pred_Den",
                  "Bottom_substrate", "Habitat", "Year"), as.factor))

str(pond_data)

#Look at multicollinearity of environmental variables
# Select numerical columns and compute correlation matrix
#Already removed Perimeter, Arhtropod_Diversity and 
#Plant_Div due to high correlations
pond_corr <- pond_data %>% 
  select(Average_Volume, Conductivity, PH, 
         Perimeter, Arthropod_Diversity, Plant_Diversity, 
         Veg_Cov, Average_Temp) %>%
  cor(use = "pairwise.complete.obs")  # Avoids NA issues

# Plot the correlation matrix
ggcorrplot(pond_corr, lab = TRUE, type = "lower", hc.order = TRUE)

#Distribution of factor levels
table(pond_data$Habitat, pond_data$Bottom_substrate)
#Remove habitat, use only bottom substrate.

#Standardize numeric variables(?)
pond_data <- pond_data %>%
  mutate(across(where(is.numeric) & !all_of(c("PondID", "Mean_Larvae1",
                                              "Mean_Larvae2", "Larvae_Sampled1",
                                              "Larvae_Sampled2")), ~ as.numeric(scale(.))))

#Interestig correlations: 

ggplot(data = pond_data, aes (x = Pred_Den, y =Arthropod_Diversity)) + 
  geom_point() + geom_smooth(method = 'lm')

ggplot(data = pond_data, aes (x = Average_Temp, y =Arthropod_Diversity)) + 
  geom_point() + geom_smooth(method = 'lm')

ggplot(data = pond_data, aes (x = Average_Temp, y =Average_Volume)) + 
  geom_point() + geom_smooth(method = 'lm')

#Load emergence trap data

emg_2024 <- read.csv("Emergence2024.csv", 
                     sep=";", header=TRUE, as.is=TRUE)

emg_2025 <- read.csv("Emergence2025.csv", 
                    sep=";", header=TRUE, as.is=TRUE)

long_2025 <- read.csv("Longevity2025_Clean.csv", 
                      sep=";", header=TRUE, as.is=TRUE)

long_2025 <- long_2025[,2:8]

emg_dat <- bind_rows(emg_2024, emg_2025)

#Drop comment column
emg_data <- emg_data[,1:9]

emg_data <- bind_rows(emg_dat, long_2025)

#Prep data
emg_data <- emg_data %>% 
  mutate(across(where(is.character) & c("Weight", "Wing_Length",
                                        "Body_Length",
                                        "NumericDate", "Longevity"),
                ~ as.numeric(gsub(",", ".", .))))

emg_data <- emg_data %>%
  mutate(across(c("Species", "Sex", "PondID"), as.factor))

str(emg_data)

# Compute first, last, and peak emergence date per pond and species
emergence_summary <- emg_data %>%
  group_by(PondID, Species) %>%
  summarise(
    first_emergence = min(NumericDate, na.rm = TRUE),   # First recorded emergence
    last_emergence = max(NumericDate, na.rm = TRUE),    # Last recorded emergence
    peak_emergence = NumericDate[which.max(table(NumericDate))], # Date with most observations
    .groups = "drop"
  )

#Remove NA
emergence_summary <- emergence_summary %>% 
  filter(!is.na(Species))

#Plot the relationships between wing lenght, body length and weight

p1 <- ggplot(data= emg_data, aes(x = Weight, y = Wing_Length)) + geom_point() +
  geom_smooth(method = 'lm'
  )

p2 <- ggplot(data= emg_data, aes(x = Weight, y = Body_Length)) + geom_point() +
  geom_smooth(method = 'lm'
  )

p3 <- ggplot(data= emg_data, aes(x = Body_Length, y = Wing_Length)) + geom_point() +
  geom_smooth(method = 'lm'
  )

(p1+p2+p3)

#Use wing length and body length as responses

# Step 1: Merge "emg_data" with "emergence_summary" using "PondID"
merged_data <- emg_data %>%
  left_join(emergence_summary, by = c("PondID", "Species"))

temporal_data <- emergence_summary %>%
  left_join(pond_data, by ="PondID")

# Step 2: Merge the result with "pond_data" to include environmental variables
full_dataset <- merged_data %>%
  left_join(pond_data, by = "PondID")

#NOTE: THIS DATASET HAS ONLY DATA FOR THE PONDS WITH EMERGENCE DATA


# Save to CSV file
write_csv2(full_dataset, "full_combined_data24-25.csv")


#Divide into one dataframe for nigripes and one for impiger
full_dataset_N <- full_dataset %>% filter(Species == "N")  # Subset for species N
full_dataset_I <- full_dataset %>% filter(Species == "I")  # Subset for species I

full_dataset_N$PondID <- as.factor(full_dataset_N$PondID)

#Model the data
#Model1:: What factors influence the mean number of larvae in the pond?
#Since we don't know exactly when in their development we reach the larvae,
#this is a very shaky test in my opinion.

larv_mod1 <- lmer(data = pond_data, 
                Mean_Larvae1 ~ 
                  Average_Volume +
                  Habitat +
                  Pred_Den +
                  Conductivity + 
                  PH + 
                  Temp_PE + 
                  Veg_Cov + 
                  Bottom_substrate +
                  (1|Year))



summary(larv_mod1)
Anova(larv_mod1, type = c('2'))

check_normality(larv_mod1)  # Tests if residuals are normally distributed
check_heteroscedasticity(larv_mod1)  # Checks if residual variance is consistent
check_collinearity(larv_mod1)
check_outliers(larv_mod1)
vif(larv_mod1)


#Model2:: What factors influence the weight of emerging mosquitoes?
##Values of EriophorumPresAbs and Veg_Cov were identical, so Eriophorum removed
#Habitat removed because it fucked the model completely. Still unsure why...
#Plant_Diversity removed  because it was moderately correlated with Veg_Cov
#Arthropod diversity removed due to correlation with temperature, volume and predator density
#Likely to be more predators where there is higher diversity!

#Body length
Emergence_Body_N <- lmer(data = full_dataset_N, 
                           Body_Length ~ 
                           Pred_Den +
                           Average_Volume +
                           Habitat +
                           Bottom_substrate +
                           Conductivity +
                           PH + 
                          Temp_PE + 
                             Veg_Cov + 
                           Sex +
                           NumericDate +
                           (1|PondID) +
                             (1|Year))


summary(Emergence_Body_N)
Anova(Emergence_Body_N, type = c('2')) #PH, habitat and Veg_Cov significant
#nope not when year is included xD
AIC(Emergence_Body_N)

check_normality(Emergence_Body_N)  # Tests if residuals are normally distributed
check_heteroscedasticity(Emergence_Body_N)  # Checks if residual variance is consistent
check_collinearity(Emergence_Body_N)
sim_res <- simulateResiduals(Emergence_Body_N)
plot(sim_res)
check_outliers(Emergence_Body_N)
plot(Emergence_Body_N)      # residuals vs fitted
qqnorm(resid(Emergence_Body_N))
qqline(resid(Emergence_Body_N))
vif(Emergence_Body_N)

#Wing length. 
#AKA nothing influences wing length. But body length does, so let's focus on that
Emergence_Wing_N <- lmer(data = full_dataset_N, 
                         Wing_Length ~ 
                           Average_Volume + 
                           Pred_Den +
                           Habitat +
                           Conductivity + 
                           PH + 
                           Temp_PE + 
                           Veg_Cov + 
                           Bottom_substrate +
                           NumericDate +
                           Sex +
                           (1|PondID) +
                           (1|Year))


plot(Emergence_Wing_N)      # residuals vs fitted
qqnorm(resid(Emergence_Wing_N))
qqline(resid(Emergence_Wing_N))

summary(Emergence_Wing_N)
Anova(Emergence_Wing_N, type = c('2')) #PH, habitat and Veg_Cov significant
AIC(Emergence_Wing_N)

check_normality(Emergence_Wing_N)  # Tests if residuals are normally distributed
check_heteroscedasticity(Emergence_Wing_N)  # Checks if residual variance is consistent
check_collinearity(Emergence_Wing_N)
check_outliers(Emergence_Wing_N)
vif(Emergence_Wing_N)

#What factors influence longevity

long_mod <- lmer(data = full_dataset_N,
                 Longevity ~ 
                   Average_Volume + 
                   Pred_Den +
                   Habitat +
                   Conductivity + 
                   PH + 
                   Temp_PE + 
                   Veg_Cov + 
                   Bottom_substrate +
                   Sex +
                   (1|PondID))

Anova(long_mod, type = c('2')) #pH, conductivity and temp significant! 
#But really though or just phacking artifact??

plot(long_mod)      # residuals vs fitted
qqnorm(resid(long_mod))
qqline(resid(long_mod))

ggplot(data = full_dataset_N, aes(x = Longevity, fill = Sex)) + 
  geom_bar()


#Plotting longevity against pH, conductivity, temp

plot_habitat <- ggplot(full_dataset_N, aes(x = Habitat, y = Body_Length)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE, aes(color = Habitat)) +
  geom_point(aes(color= Habitat)) +
  theme_minimal() +
  ylab("Body length") +
  xlab("Habitat") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )

plot_habitat


#Predict for continuous variables

conductivity_pred<- ggpredict(Emergence_Body_N, terms = "Veg_Cov")

plotpred_veg <- ggplot(vegcov_pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 0.8, color = "firebrick") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "red", alpha = 0.3) +
  geom_smooth() +
  geom_point() +
  theme_minimal() +
  ylab("Predicted Body length") +
  xlab("Vegetation cover") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )

plotpred_veg

plot_veg <- ggplot(full_dataset_N, aes(x = Veg_Cov, y = Body_Length)) +
  stat_smooth(method= 'lm', linewidth = 0.8, color = "black") +
  geom_point() +
  theme_minimal() +
  ylab("Body length") +
  xlab("Vegetation cover") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black"))

plot_veg

ph_pred<- ggpredict(Emergence_Body_N, terms = "PH")

plotpred_ph <- ggplot(ph_pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 0.8, color = "firebrick") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "red", alpha = 0.3) +
  geom_smooth() +
  geom_point() +
  theme_minimal() +
  ylab("Predicted Body length") +
  xlab("PH") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )

plotpred_ph

plot_ph <- ggplot(full_dataset_N, aes(x = PH, y = Body_Length)) +
  stat_smooth(method= 'lm', linewidth = 0.8, color = "black") +
  geom_point() +
  theme_minimal() +
  ylab("Body length") +
  xlab("PH") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black"))

plot_ph

observed <- plot_ph + plot_veg + plot_habitat

ggsave(plot = observed, filename = "Observed_PondEffects_Body.png",
       height = 5.26, width = 13, dpi = 300)

#Check so rocky habitat (low body wieght) is not correlated with
#Veg cover or pH

hab_veg_plot <- ggplot(full_dataset_N, aes(x = Habitat, y = Veg_Cov))+
  geom_boxplot()+
  geom_point()

hab_pH_plot <- ggplot(full_dataset_N, aes(x = Habitat, y = PH))+
  geom_boxplot()+
  geom_point()

pH_veg_plot <- ggplot(full_dataset_N, aes(x = PH, y = Veg_Cov))+
  geom_point()

correlations <- pH_veg_plot + hab_pH_plot + hab_veg_plot

ggsave(plot = correlations, filename = "pred_corrs.png",
       height = 5.26, width = 13, dpi = 300)

testplot <- ggplot(full_dataset_N, aes(x = Veg_Cov, y = Body_Length))+
  geom_point(aes(color = Habitat))+
  stat_smooth(aes(color = Habitat), method = "lm")

ggsave(plot = testplot, filename = "corr_VegBodHab.png",
       height = 5.26, width = 6.5, dpi = 300)

#We need more data on Impiger to make a good model.
Emergence_Wing_I <- lmer(data = full_dataset_I, 
                           Body_Length ~ 
                             Habitat +
                             Conductivity + 
                             PH + 
                           Veg_Cov +
                             (1|PondID) +
                             (1|NumericDate))

summary(Emergence_Wing_I)
Anova(Emergence_Wing_I, type = c('3'))
AIC(Emergence_Wing_I)

check_model(Emergence_Weight_I)
check_normality(Emergence_Weight_N)  # Tests if residuals are normally distributed
check_heteroscedasticity(Emergence_Weight_N)  # Checks if residual variance is consistent
check_collinearity(Emergence_Weight_N)
check_outliers(Emergence_Weight_N)
vif(Emergence_Weight_N)

#Model3&4:: What factors influence peak emergence date and end emergence date?
#Using only the variables from both years

temporal_N <- temporal_data %>% filter(Species == "N")  # Subset for species N
temporal_I <- temporal_data %>% filter(Species == "I")  # Subset for species I

EmergencePeak_N <- lmer(data = temporal_N, 
                       peak_emergence ~ 
                         Average_Volume + 
                         Habitat +
                         Bottom_substrate +
                         Conductivity + 
                         PH + 
                         Temp_PE + 
                         Pred_Den +
                         Veg_Cov +
                         PondID +
                        (1|Year))


summary(EmergencePeak_N)
Anova(EmergencePeak_N, type = c('2'))

check_normality(EmergencePeak_N)  # Tests if residuals are normally distributed
check_heteroscedasticity(EmergencePeak_N)  # Checks if residual variance is consistent
sim_res <- simulateResiduals(EmergencePeak_N)
plot(sim_res)
plot(EmergencePeak_N)      # residuals vs fitted
qqnorm(resid(EmergencePeak_N))
qqline(resid(EmergencePeak_N))
check_collinearity(EmergencePeak_N)
check_outliers(EmergencePeak_N)

#Habitat and PH important!

#Impiger has too little data to make any good models
#We can't really trust end emergence or start emergence so those are skipped

ggplot() +
  geom_boxplot(data = full_dataset_N, aes(x = as.factor(PondID), y = Wing_Length))

ggplot()+
  geom_boxplot(data = full_dataset_N, aes(x = as.factor(PondID), y = Body_Length))
