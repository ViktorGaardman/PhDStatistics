library("tidyverse")
library('forcats')
library('glmmTMB')
library('car')
library('ggeffects')
library('DHARMa')
library('emmeans')

rm(list = ls())

setwd("C:/Users/vikto/OneDrive/Skrivbord/PhD/R/PhDStatistics")


options(contrasts = c("contr.sum", "contr.poly"))

set.seed(123)

df <- read.csv("DirectPredation_CombinedClean.csv", 
               sep=";", header=TRUE, as.is=TRUE)

head(df)

df <- df %>%
  mutate_at(c("ID", "Sex",
              "Beetles"), as.factor) %>% 
  mutate_at(c("Temp", "Added", "Dead",
              "Alive", "Total", "Exp_Day"), as.numeric)

df <- df %>%
  filter(!Exp_Day %in% 0)

df <- df %>%
  filter(!Beetles %in% 2)

df <- df %>%
  filter(!Sex %in% "Control")

#Remove 9 degrees, it was made so late in season that data is not comparable
df <- df %>% 
  filter(!Temp %in% 9)

#Add obs level ID
df$obs <- 1:nrow(df)

#Basic plot
pred_plot <- ggplot(data = df, aes(Temp, Dead, color = Sex))+
  geom_jitter(alpha = 0.6,
              width = 0.1) +
  geom_smooth(linewidth = 0.8) +
  theme_minimal() +
  ylab("Daily consumption") +
  xlab("Temperature") +
  scale_x_continuous(limits = c(2, 17), breaks = scales::pretty_breaks(n = 8)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black"))

pred_plot

ggsave("Observed_Direct.png", plot = pred_plot, 
       width = 6.5, height = 5.26, dpi = 450)

#fit model
predation_mod <- glmmTMB(cbind(Dead, Alive) 
                     ~ Temp + Sex + (1|ID) + (1|Exp_Day) + (1|obs),
                     data=df, family = binomial(link = "logit"))

car::Anova(predation_mod, type = 'III')

#Check assumptions
testResiduals(predation_mod, plot = TRUE)

plot(residuals(predation_mod) ~
       predict(predation_mod,type="link"),xlab=expression(hat(eta)),
     ylab="Deviance residuals",pch=20,col="blue")

performance::check_overdispersion(predation_mod) 

#Predict feeding rates per sex

emm <- emmeans(
  predation_mod,
  ~ Temp | Sex,
  at = list(Temp = seq(min(df$Temp), max(df$Temp), length.out = 100)),
  type = "response"
)

plot_pred <- as.data.frame(emm)

predicted_pred <- ggplot(plot_pred, aes(x = Temp, y = prob, color = Sex,
                                        fill = Sex)) +
  geom_ribbon(data = plot_pred,
              aes(x = Temp, ymin = asymp.LCL, ymax = asymp.UCL),
              alpha = 0.1) +
  geom_jitter(data = df, aes(x = Temp, y = (Dead/20), color = Sex), alpha = 0.6,
             width = 0.3) +
  geom_line(data = plot_pred, aes(x = Temp, y = prob, color = Sex), linewidth = 1) +

    scale_color_manual(values = c("#009E73", "#E69F00")) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  scale_x_continuous(limits = c(2, 17), breaks = scales::pretty_breaks(n = 8)) +
  theme_classic() +
  ylab("Predicted proportion consumed") +
  xlab("Temperature (°C)") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

predicted_pred

ggsave("Predicted_Direct.png", plot = predicted_pred, 
       width = 6.5, height = 5.26, dpi = 450)
