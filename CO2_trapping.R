library(tidyverse)
library(readr)
library(lme4)
library(car)
library(emmeans)
library(performance)
library(lubridate)
library(DHARMa)
library(patchwork)


rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

set.seed(123)

df <- read.csv("TrappingCombinations.csv", sep =";", h = T)

#Sampled individuals per hour
df <- df %>%
  filter(!is.na(Family)) %>%
  group_by(ID, Family, Species, TrapType, Date_Num) %>%
  summarize(
    CountPerHour = Count / Hours_sampling
  )


df$TrapType <- as_factor(df$TrapType)

df_fam <- df %>% 
  group_by(Family, TrapType, Date_Num) %>% 
  summarize(
    SumCountPerHour = sum(CountPerHour)
  )

df <- df %>% 
  mutate(
    Species = fct_recode(
      Species,
      "Oc. impiger" = "Oc_impiger",
      "Oc. nigripes" = "Oc_nigripes",
      "S. vittatum" = "S_vittatum",
      "S. rostratum" = "S_rostratum"
    )) %>% 
  mutate(Species = fct_relevel(
    Species,
    "Oc. impiger",
    "S. rostratum",
    "Oc. nigripes",
    "S. vittatum"))

species_boxplot <- ggplot(df, aes(x = TrapType, y = CountPerHour, color=
                                         Family)) +
  facet_wrap(~ Species, scales = "free_y") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_classic()+
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  ylab("Individuals per hour")+
#  scale_y_continuous(limits = c(0,20)) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )

species_boxplot


family_boxplot <- ggplot(df, aes(x = TrapType, y = CountPerHour, color=
                                    Family)) +
  facet_wrap(~ Family, scales = "free_y") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_classic()+
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  ylab("Individuals per hour")+
#  scale_y_continuous(limits = c(0,20)) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )

family_boxplot

family_boxplot_2 <- family_boxplot / plot_spacer()



df$Family <- as.factor(df$Family)
df$TrapType <- as.factor(df$TrapType)
df$CountPerHour <- as.numeric(df$CountPerHour)


compmod_Fam <- glmmTMB::glmmTMB(
  CountPerHour ~ Family * TrapType + (1|Date_Num),
  family = Gamma(link = "log"),
  data = df
)

simulationOutput <- simulateResiduals(fittedModel = compmod_Fam)


png("ModDiag_Bait_Fam.png", width = 3000, height = 1500, res = 300)



plot(simulationOutput)

dev.off()

testOutliers(simulationOutput)
testDispersion(simulationOutput)
testUniformity(simulationOutput)

summary(compmod_Fam)

Anova(compmod_Fam, type = "3")

pair <- emmeans(compmod_Fam, ~ TrapType | Family, 
                type = "response", adjust = "none")

pairs(pair)
print(pair)

comp_pred <- emmeans(
  compmod_Fam,
  ~ Family | TrapType,
  type = "response"
)

plot_pred_comp <- as.data.frame(comp_pred)

predicted_comppred <- ggplot() +
  geom_jitter(data = df, aes(x = TrapType, y = CountPerHour, color = Family), alpha = 0.5,
              width = 0.1) +
  geom_errorbar(
    data = plot_pred_comp,
    aes(
      x = TrapType,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = Family
    ),
    width = 0.2,
    linewidth = 1
  ) +
  geom_point(data = plot_pred_comp, aes(x = TrapType, y = response, color = Family), size = 4) +
  #  scale_y_continuous(limits = c(-2, 10), breaks = scales::pretty_breaks(n = 11)) +
  theme_classic() +
  ylab("Individuals per hour") +
  facet_wrap(~ Family, scales = "free_y") +
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14)
  )

predicted_comppred


boxplots <- (predicted_comppred / species_boxplot) + plot_annotation(tag_levels = "A")

boxplots_h <- boxplots + plot_layout(heights = c(1, 2))

ggsave(boxplots_h, filename = "Fig2_CO2.png", dpi = 450,
       height = 10.56, width = 10)



type3 <- list(TrapType = contr.sum, Family = contr.sum)

mod <- lmer(log(SumCountHour+0.1) ~
            TrapType * 
              Family +
            (1|Date_Num),
            contrasts = type3,
          data = df_fam)

plot(mod)
Anova(mod, type = 'III')

emm <- emmeans(mod, specs = pairwise ~ TrapType)
emm$contrasts

##C-flux data

files <- list.files(
  "Viktor_c_flux",
  pattern = "\\.txt$",
  full.names = TRUE,
  recursive = TRUE
)

#Find header line
lines <- readLines(files[1], warn = FALSE)
header_line <- grep("DATE/TIME, DATA FORMAT", lines)
header_line #3


read_logger_file <- function(f) {
  lines <- readLines(f, warn = FALSE)
  
  # Find the real header
  header_line <- grep("DATE/TIME, DATA FORMAT", lines)
  
  # Keep everything from the header line onwards
  lines <- lines[header_line:length(lines)]
  
  # Remove wrapping quotes if any
  lines <- ifelse(
    str_starts(lines, '"') & str_ends(lines, '"'),
    str_sub(lines, 2, -2),
    lines
  )
  
  data <- read_delim(
    I(lines),
    delim = ",",
    trim_ws = TRUE,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  # Add file name column
  data <- data %>% mutate(SourceFile = basename(f))
  
  return(data)
}


CO2_data <- map_dfr(files, read_logger_file) 

CO2_data$CO2 <- as.numeric(CO2_data$CO2)

CO2_data <- CO2_data %>%
  filter(!is.na(CO2)) 

CO2_data <- CO2_data %>%
  mutate(TAIR = as.numeric(TAIR)) %>%
  mutate(AIRPRES = as.numeric(`AIR PRESSURE`)) %>%
  mutate(
    EventSec = period_to_seconds(hms(`EVENT TIME`))
  )

CO2_data <- CO2_data %>%
  mutate(
    Hour = str_extract(SourceFile, "(?<=_)[0-9]+(?=\\.txt)")
  ) %>%
  mutate(Hour = as.numeric(Hour))  # optional, makes it numeric

CO2_data$Hour <- as.factor(CO2_data$Hour)
CO2_data$SourceFile <- as.factor(CO2_data$SourceFile)


CO2_data$`EVENT DATE` <- as.factor(CO2_data$`EVENT DATE`)
levels(CO2_data$`EVENT DATE`)
#%>%
#  filter(CO2 >= 600, CO2 <= 2010)

write.csv(CO2_data, "CO2_data_cleaned.csv")

CO2_data <- read_csv("CO2_data_cleaned.csv")

ggplot(subset(CO2_data, SourceFile == "bag1_10.txt"),
       aes(x = EventSec, y = CO2)) +
  geom_point()+
  geom_smooth()+
facet_wrap(~ `EVENT DATE`, scales = "free_x")

ggplot(subset(CO2_filtered, `EVENT DATE` == "10/06/25"),
       aes(x = EventSec, y = CO2)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~ SourceFile, scales = "free_x")

#Remove wonky measurements:
CO2_cleaned <- CO2_data %>%
  filter(!(`EVENT DATE` == "13/05/25" & SourceFile %in% c("bag1_9.txt",
                                                        "bag2_9.txt",
                                                        "bag3_13.txt",
                                                        "bag2_15.txt",
                                                        "bag1_16.txt",
                                                        "bag2_19.txt")))%>%
  filter(!(`EVENT DATE` == "25/05/25" & SourceFile %in% c("bag2_10.txt",
                                                          "bag3_11.txt",
                                                        "bag1_12.txt",
                                                        "bag2_12.txt",
                                                        "bag3_12.txt",
                                                        "bag2_13.txt",
                                                        "bag3_13.txt",
                                                        "bag3_14.txt",
                                                        "bag1_15.txt",
                                                        "bag2_15.txt",
                                                        "bag3_15.txt",
                                                        "bag1_16.txt",
                                                        "bag3_16.txt",
                                                        "bag2_17.txt",
                                                        "bag3_17.txt",
                                                        "bag1_18.txt",
                                                        "bag2_18.txt",
                                                        "bag3_18.txt",
                                                        "bag1_20.txt",
                                                        "bag1_22.txt",
                                                        "bag2_22.txt",
                                                        "bag3_22.txt")))%>%
  filter(!(`EVENT DATE` == "10/06/25" & SourceFile %in% c("bag2_10.txt",
                                                        "bag1_14.txt",
                                                        "bag2_15.txt",
                                                        "bag1_16.txt",
                                                        "bag3_18.txt")))



ggplot(subset(CO2_cleaned, SourceFile == "bag3_16.txt"),
       aes(x = EventSec, y = CO2)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~ `EVENT DATE`, scales = "free_x")

#Trim some measurements to start at the right time:
CO2_cleaned <- CO2_cleaned %>%
  filter(
    !(`EVENT DATE` == "10/06/25" & SourceFile == "bag1_12.txt") |
      EventSec >= 36690
  )

CO2_cleaned <- CO2_cleaned %>%
  filter(
    !(`EVENT DATE` == "13/05/25" & SourceFile == "bag2_16.txt") |
      EventSec >= 50306
  )

CO2_cleaned <- CO2_cleaned %>%
  filter(
    !(`EVENT DATE` == "25/05/25" & SourceFile == "bag3_20.txt") |
      EventSec >= 65275
  )

CO2_cleaned <- CO2_cleaned %>%
  filter(
    !(`EVENT DATE` == "25/05/25" & SourceFile == "bag3_11.txt") |
      EventSec >= 32030
  )

write.csv(CO2_cleaned, "CO2_data_trimmed.csv")

CO2_cleaned <- read_csv("CO2_data_trimmed.csv")
#Calculate ml CO2 per minute
#First per second, then times 60
#USe only values between 1000-2000ppm

#Get time in seconds
CO2_filtered <- CO2_cleaned %>%
  filter(CO2 >= 1000, CO2 <= 2000) %>%
  mutate(TAIR = as.numeric(TAIR)) %>%
  mutate(AIRPRES = as.numeric(`AIR PRESSURE`)) %>%
  mutate(
    EventSec = period_to_seconds(hms(`EVENT TIME`))
  )

ggplot(subset(CO2_filtered, `EVENT DATE` == "25/05/25"),
       aes(x = EventSec, y = CO2)) +
  geom_point()+
  geom_smooth()+
  facet_wrap(~ SourceFile, scales = "free_x")

CO2_min <- CO2_filtered %>%
  group_by(`EVENT DATE`, SourceFile) %>%
  summarize(
    TotalRelease = max(CO2) - min(CO2), #How much was released
    TimeSec = max(EventSec) - min(EventSec), #During how many seconds
    ReleasePerMin = TotalRelease/TimeSec * 60,
    Temp = mean(TAIR),
    Pressure = mean(AIRPRES)
  )

CO2_min <- CO2_min %>%
  mutate(
    Hour = str_extract(SourceFile, "(?<=_)[0-9]+(?=\\.txt)")
  ) %>%
  mutate(Hour = as.numeric(Hour))  # optional, makes it numeric

#Calculate ppm to ml
#ml = Volume of chamber in mm * ppm/10^6 * Pressure * 273.15/Temp in Kelvin

CO2_min <- CO2_min %>%
  mutate(
    Kelvin = Temp + 273.15
  ) %>%
  mutate(Kelvin = as.numeric(Kelvin))

#Calculate frame size
#33x33x34cm = 0.33*0.33*0.34m 
0.33* 0.33* 0.34 #0.037026 m^3 = 37 litres = 37 000 ml
#Standard pressure 1013.25 mbar
#OBS: Pressure and temp already accounted for in CO2 measurements!

CO2_min <- CO2_min %>%
  mutate(
    mlpermin = 37000 * ReleasePerMin/10^6
  ) 

filtered <- CO2_min %>%
  filter(Hour >= 13, Hour <= 16)

max(CO2_min$mlpermin)
min(CO2_min$mlpermin)
mean(CO2_min$mlpermin)

write.csv(CO2_min, "CO2_per_minute.csv")

CO2_min <- read_csv("CO2_per_minute.csv")

filter <- CO2_min %>%
  filter(Hour == 22)
min(filter$mlpermin)
max(filter$mlpermin)
mean(filter$mlpermin)



#Add dummy 0s at 8:00
CO2_df <- CO2_min |> 
  mutate(
    Year = str_sub(`EVENT DATE`, 7,8),
    Month = str_sub(`EVENT DATE`, 4,5),
    Day = str_sub(`EVENT DATE`, 1, 2)# extract MM from DDMMYY
  )

CO2_df <- CO2_df |> 
  mutate(
    Day = fct_recode(
      Day,
      "03" = "13",
      "15" = "25",
      "31" = "10",
      
    ))

CO2_filt <- CO2_df |> 
  select(c(Hour, mlpermin, Day))

CO2_min0 <- rbind(CO2_filt, data.frame(Hour = 8, mlpermin = 0, Day = c("03", "15", "31")))


CO2_plot <- CO2_min0 %>%
  mutate(
    HourCult = Hour - 8
  )

#Create a linear graph of release rates
CO2_mean <- CO2_plot |>
  group_by(HourCult) |> 
  summarize(
    mean_rate = mean(mlpermin),
    sd_rate = sd(mlpermin)
  )

CO2_plot <- CO2_plot |> 
  mutate(
    Day = fct_relevel(Day, "03", "15", "31")
  )

min(CO2_plot$mlpermin)

#Look at data
mlpermin_plotraw <- ggplot(CO2_plot, aes(x = HourCult, y = mlpermin)) + 
  theme_bw()+
  xlab("Hours of fermentation")+
  ylab("ml CO2 / min")+
  geom_point(color = "grey20", size = 2, alpha = 0.7)+
#  scale_color_manual(values = c("#D55E00", "#56B4e9", "#CC79A7")) +
  geom_hline(yintercept = 185, col = "black", linetype = "dashed", lwd = 0.5) +
  geom_line(
    data = CO2_mean,
    aes(y = mean_rate,
    x = HourCult),
    linewidth = 0.8
  ) +
  geom_errorbar(
    data = CO2_mean,
    aes(
    x = HourCult,
    ymax = mean_rate + sd_rate,
    ymin = mean_rate - sd_rate),
    width = 0.2,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +

#  geom_smooth(
#    method = "gam",
#    formula = y ~ poly(x, 4),
#    se = FALSE,
#    color = "black"
#  ) +
  scale_x_continuous(limits = c(0,14), n.breaks = 8)+
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
        )

mlpermin_plotraw

mlpermin_plots <- (mlpermin_plotraw / temp_timeplot) + plot_annotation(tag_levels = 'A')

ggsave(plot = mlpermin_plotraw, filename = "Fig3_Co2.png",
       dpi = 450, height= 4, width = 6.5)

#model
#First calculate mean temp per day
CO2_temp <- CO2_min %>%
  group_by(`EVENT DATE`) %>%
  summarise(
    meantemp = mean(Temp),
    mltotal = sum(mlpermin)
  )


mlmod <- lmer(log(mlpermin) ~
              Temp +
                (1|`EVENT DATE`),
            data = CO2_min)

plot(mlmod)
summary(mlmod)



tempml_plot <- ggplot(CO2_min, aes(x = Temp, y = mlpermin)) + 
  theme_bw()+
  stat_smooth(method = 'gam', 
              formula = y ~ poly(x, 3),
              se = FALSE,
              color = "cornflowerblue"
              )+
  geom_point()+
  xlab("Temperature (°C)")+
  ylab("ml CO2 / min")+
  scale_x_continuous(limits = c(12,32), n.breaks = 12)+
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())

tempml_plot

ggsave(plot = tempml_plot, filename = "Temp_CO2.TIFF",
       dpi = 300, height = 3.5, width = 6.5)

temptime_plot <- ggplot(CO2_min, aes(x = Hour, y = Temp)) + 
  theme_bw()+
  stat_smooth(method = 'gam', 
              formula = y ~ poly(x, 3),
              se = FALSE,
              color = "cornflowerblue"
  )+
  geom_point()+
  xlab("Hour of day")+
  ylab("Temperature (°C)")+
  scale_x_continuous(limits = c(9,22), n.breaks = 8)+
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())

temptime_plot

ggsave(plot = temptime_plot, filename = "Temp_Time.TIFF",
       dpi = 300, height = 3.5, width = 6.5)

head(CO2_min)


####Cannister comparison

comp_df <- read.csv("YeastVsCannister.csv", sep =";", h = T)

comp_df <- comp_df %>% 
  mutate_at(c("Site", "Type",
            "Date", "Species", "Group"), as.factor)

summ <- comp_df |> 
  group_by(Group, Type) |> 
  summarize(
    sum = sum(Count)
  )

#Make a gamma model instead of lmer(?)
comp_mod <- glmmTMB::glmmTMB(
  CountPerHour ~ Type * Species + (1|Date) + (1|Site),
    family = Gamma(link = "log"),
    data = comp_df
)


simulationOutput <- simulateResiduals(fittedModel = comp_mod)


png("ModDiag_CO2_Species.png", width = 3000, height = 1500, res = 300)



plot(simulationOutput)

dev.off()

testOutliers(simulationOutput)
testDispersion(simulationOutput)
testUniformity(simulationOutput)

summary(comp_mod)

Anova(comp_mod, type = "3")

pair <- emmeans(comp_mod, ~ Type | Species, 
                type = "response", adjust = "none")

pairs(pair)
print(pair)

emm <- emmeans(
  comp_mod,
  ~ Type | Species,
  type = "response"
)

plot_pred <- as.data.frame(emm)

plot_pred <- plot_pred %>% 
  mutate(
    Species = fct_recode(
    Species,
    "Oc. impiger" = "Oc_impiger",
    "Oc. nigripes" = "Oc_nigripes",
    "S. vittatum" = "S_vittatum",
    "S. rostratum" = "S_rostratum"
  )) %>% 
  mutate(Species = fct_relevel(
    Species,
    "Oc. impiger",
    "S. rostratum",
    "Oc. nigripes",
    "S. vittatum"))

comp_df <- comp_df %>% 
  mutate(
    Species = fct_recode(
      Species,
      "Oc. impiger" = "Oc_impiger",
      "Oc. nigripes" = "Oc_nigripes",
      "S. vittatum" = "S_vittatum",
      "S. rostratum" = "S_rostratum"
    )) %>%  mutate(Species = fct_relevel(
         Species,
         "Oc. impiger",
         "S. rostratum",
         "Oc. nigripes",
         "S. vittatum"))



predicted_pred <- ggplot() +
  geom_jitter(data = comp_df, aes(x = Type, y = CountPerHour, color = Type), alpha = 0.5,
              width = 0.1) +
  geom_errorbar(
    data = plot_pred,
    aes(
      x = Type,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = Type
    ),
    width = 0.2,
    linewidth = 1
  ) +
  geom_point(data = plot_pred, aes(x = Type, y = response, color = Type), size = 4) +
  #  scale_y_continuous(limits = c(-2, 10), breaks = scales::pretty_breaks(n = 11)) +
  theme_classic() +
  ylab("Individuals per hour") +
  facet_wrap(~ Species, scales = "free_y") +
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.ticks.x = element_blank()
  )

predicted_pred


#Gamma model groups

group_df <- comp_df %>% 
  group_by(Group, Date, Type, Site) %>% 
  summarize(
    Sum = sum(CountPerHour)
  )

group_mod <- glmmTMB::glmmTMB(
  Sum ~ Type * Group + (1|Date) + (1|Site),
  family = Gamma(link = "log"),
  data = group_df
)

simulationOutput <- simulateResiduals(fittedModel = group_mod)


png("ModDiag_CO2_Fam.png", width = 3000, height = 1500, res = 300)



plot(simulationOutput)

dev.off()

testOutliers(simulationOutput)
testDispersion(simulationOutput)
testUniformity(simulationOutput)

summary(group_mod)

Anova(group_mod, type = "3")

pair <- emmeans(group_mod, ~ Type | Group, 
                type = "response", adjust = "none")

pairs(pair)
print(pair)

emm_group <- emmeans(
  group_mod,
  ~ Group | Type,
  type = "response"
)

group_pred <- as.data.frame(emm_group)

predicted_group <- ggplot() +
  geom_jitter(data = group_df, aes(x = Type, y = Sum, color = Type), alpha = 0.5,
              width = 0.1) +
  geom_errorbar(
    data = group_pred,
    aes(
      x = Type,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = Type
    ),
    width = 0.2,
    linewidth = 1
  ) +
  geom_point(data = group_pred, aes(x = Type, y = response, color = Type), size = 4) +
  #  scale_y_continuous(limits = c(-2, 10), breaks = scales::pretty_breaks(n = 11)) +
  theme_classic() +
  ylab("Individuals per hour") +
  facet_wrap(~ Group, scales = "free_y") +
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 14),
    legend.position = c(0.95, 0.8)
  )

predicted_group

predicted_catch <- predicted_group / predicted_pred

pred_catch <- predicted_catch + plot_annotation(tag_levels = "A") + plot_layout(heights = c(1,2))

ggsave(pred_catch, filename = "Fig4_CO2.png", dpi = 450,
       height = 10.56, width = 10)


##Calculate mean temp

Tempdf_raw <- read.delim("TempData_Kobb.txt")

Tempdf <- Tempdf_raw |> 
  mutate(
    Year = str_sub(Date, 3,4),
    Month = str_sub(Date, 6,7),
    Day = str_sub(Date, 9, 10)# extract MM from DDMMYY
    )

Tempdf <- Tempdf |> 
  filter(Year == "25",
         Month == "07", 
         Day %in% c("03", "15" ,"31")) 

Tempdf <- Tempdf |>
  filter(
    Time >= "08:00:00",
    Time <= "22:00:00"
  )

Tempdf <- Tempdf |> 
  mutate(time_sec = period_to_seconds(hms(Time)) / 3600)

Tempdf <- Tempdf |> 
  mutate(
  Day = fct_relevel(Day, "03", "15", "31")
  )

temp_timeplot <- ggplot(Tempdf, aes(x = (time_sec - 8), y = air_temperature_2m,
                   color = Day)) +
  theme_bw()+
#  stat_smooth(method = 'gam', 
#              formula = y ~ poly(x, 3),
#              se = FALSE,
#              color = "black"
#  )+
  scale_color_manual(values = c("#D55E00", "#56B4e9", "#CC79A7")) +
  geom_point(size = 2)+
  xlab("Hours of fermentation")+
  ylab("Temperature (°C)")+
  scale_x_continuous(limits = c(0,14), n.breaks = 8)+
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())

temp_timeplot

#Mean temp during CO2 release rate
CO2_df <- CO2_min |> 
  mutate(
    Year = str_sub(`EVENT DATE`, 7,8),
    Month = str_sub(`EVENT DATE`, 4,5),
    Day = str_sub(`EVENT DATE`, 1, 2)# extract MM from DDMMYY
  )

CO2_df <- CO2_df |> 
  mutate(
  Day = fct_recode(
    Day,
    "03" = "13",
              "15" = "25",
             "31" = "10",
              
  ))

CO2_filt <- CO2_df |> 
  select(c(Hour, mlpermin, Day))

Tempdf$Hour <- Tempdf$time_sec

Tempdf_filt <- Tempdf |> 
  select(c(air_temperature_2m, Hour, Day))

Release_df <- CO2_filt |> 
  left_join(Tempdf_filt, by = c("Hour", "Day"))

release_mod <- lm(mlpermin ~ air_temperature_2m,
                  data = Release_df)

summary(release_mod)
car::Anova(release_mod, type = '2')

ml_timeplot <- ggplot(Release_df, aes(x = air_temperature_2m, y = mlpermin, color = Day)) +
  theme_bw()+
  scale_color_manual(values = c("#D55E00", "#56B4e9", "#CC79A7")) +
#  stat_smooth(method = 'lm', color = "black") +
  geom_point()+
  ylab("ml CO2 / min")+
  xlab("Temperature (°C)")+
  scale_x_continuous(
    n.breaks = 6)+
  facet_wrap(~Day, scales = "free_x") +
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())
ml_timeplot
