library(tidyverse)
library(readr)
library(lme4)
library(car)
library(emmeans)
library(performance)
library(lubridate)

df <- read.csv("TrappingCombinations.csv", sep =";", h = T)

#Sampled individuals per hour
hours <- df %>%
  group_by(ID) %>%
  summarize(
    CountPerHour = Count / Hours_sampling
  )

df <- df %>%
  filter(!is.na(Family)) %>%
  left_join(hours, by = "ID")

df_fam <- df %>%
  group_by(TrapType, Family,
           Date_Num) %>%
  summarize(
    SumCountHour = sum(CountPerHour)
  )

df_fam <- df_fam %>%
  filter(!Date_Num %in% "")

raw_boxplot <- ggplot(df_fam, aes(x = TrapType, y = SumCountHour, color=
                                         Family)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_bw()+
  ylab("Capture per hour")+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
    axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 


raw_boxplot

ggsave(plot= raw_boxplot, filename = "Attractant_boxplot.TIFF",
       dpi = 300, width = 6.5, height = 3.5)

df_tot <- df_fam %>%
  group_by(Family,
           Date_Num) %>%
  summarize(
    TotalCount = sum(SumCountHour)
  )

time_plot <- ggplot(df_tot, aes(x = Date_Num, y = TotalCount,
                                       color = Family)) +
  stat_smooth(method = 'lm', 
              se = FALSE, 
              aes(color = Family)) +
  geom_point() +
  theme_bw()+
  ylab("Capture per hour")+
  scale_y_continuous(limits = c(0,70))+
  theme(legend.position="right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        legend.direction='vertical',axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5)) 

time_plot

ggsave(plot = time_plot, filename = "Capture_time.TIFF",
       dpi = 300, height = 3.5, width = 6.5)

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

CO2_filtered <- CO2_data %>%
  filter(!is.na(CO2)) %>%
  filter(CO2 >= 600, CO2 <= 2010)

write.csv(CO2_filtered, "CO2_data_cleaned.csv")

#Calculate ml CO2 per minute
#First per second, then times 60

#Get time in seconds
CO2_filtered <- CO2_filtered %>%
  mutate(TAIR = as.numeric(TAIR)) %>%
  mutate(AIRPRES = as.numeric(`AIR PRESSURE`)) %>%
  mutate(
    EventSec = period_to_seconds(hms(`EVENT TIME`))
  )

CO2_min <- CO2_filtered %>%
  group_by(`EVENT DATE`, SourceFile) %>%
  summarize(
    TotalRelease = max(CO2) - min(CO2), #How much was released
    TimeSec = max(EventSec) - min(EventSec), #During how many seconds
    ReleasePerMin = (TotalRelease/TimeSec) * 60,
    Temp = mean(TAIR),
    Pressure = mean(AIRPRES)
  )

CO2_min <- CO2_min %>%
  mutate(
    Hour = str_extract(SourceFile, "(?<=_)[0-9]+(?=\\.txt)")
  ) %>%
  mutate(Hour = as.numeric(Hour))  # optional, makes it numeric

CO2_min <- CO2_min %>%
  filter(ReleasePerMin <= 10000)

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

write.csv(CO2_min, "CO2_per_minute.csv")

#Look at data
mlpermin_plotraw <- ggplot(CO2_min, aes(x = Hour, y = mlpermin)) + 
  geom_point()+
  theme_bw()+
  xlab("Hour of day")+
  ylab("ml CO2 / min")+
  scale_x_continuous(limits = c(9,22), n.breaks = 8)+
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())

mlpermin_plotraw

ggsave(plot = mlpermin_plotraw, filename = "ml_per_min.TIFF",
       dpi = 300, height= 3.5, width = 6.5)

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
  stat_smooth(method = 'lm', color = "cornflowerblue",
                                 fill = "cornflowerblue")+
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
