library(tidyverse)
library(readr)
library(lubridate)

rm(list = ls())

df <- read_delim("Biobasis_Nuuk_Arthropod_Raw.txt")


df_2024 <- df %>%
  filter(year(Date) == 2024)

df_2024 <- df_2024 %>%
  mutate(
    across(Species_A:Species_D, ~ na_if(., -9999))
  )

df_2024 <- df_2024 %>%
  filter(!is.na(Species_A))

df2024_sum <- df_2024 %>%
  group_by(Date, Plot, Order, Family) %>%
    mutate(
      Count = rowSums(across(Species_A:Species_D), na.rm = TRUE)
    )

sum_2024 <- df2024_sum %>%
  select(Date, Plot, Order, Family, Count) %>%
  filter(!Count  == 0)

df2024_final <- sum_2024 %>%
  group_by(Date, Order, Family) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

df2024_final_month <- df2024_final %>%
  mutate(
    Month = month(Date)
  ) %>%
  filter(!Month %in% c(9, 10)) %>%
  group_by(year(Date), Month, Order, Family) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

writexl::write_xlsx(df2024_final_month, "Pitfall2024_clean.xlsx")
