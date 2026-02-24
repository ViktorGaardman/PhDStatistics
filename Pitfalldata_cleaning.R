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


#Malaise summary
df <- read_delim("MalaiseSamples_2024.csv")

#Exchange week to date
df <- df %>%
  mutate(
    Month = cut(
      Week,
      breaks = c(-Inf, 26, 30, Inf),
      labels = c("June", "July", "August"),
      right = TRUE
    )
  )

#Sum across plots
df2 <- df %>%
  select(TrapNumber, Month, Order, Taxon, Count)

df_sum <- df2 %>%
  group_by(Month, Order, Taxon) %>%
  summarise(
    Total = sum(Count, na.rm = TRUE),
    .groups = "drop"
  )

writexl::write_xlsx(df_sum, "Malaise_pie.xlsx")


###

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

