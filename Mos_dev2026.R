library(tidyverse)
library(readr)
library(patchwork)


rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

df <- read.csv("Mos_dev2026.csv", sep =";", h = T)

head(df)

df <- df %>% 
  mutate(across(c(Temp, ID, State), as.factor)) %>% 
  filter(!State %in% "")


ggplot(data = df, aes(x = Exp_day, fill = State)) +
  geom_bar() +
  facet_wrap(~Temp, scale = "free") +
  xlab("Experimental day") +
  ylab("Count") +
  theme_classic()
         