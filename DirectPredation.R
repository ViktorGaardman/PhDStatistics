library("tidyverse")
library('forcats')
library('glmmTMB')
library('car')
library('ggeffects')
library('DHARMa')

rm(list = ls())

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

#fit model
predation_mod <- glmmTMB(cbind(Dead, Alive) 
                     ~ Temp + (1|ID) + (1|Exp_Day),
                     data=df, family = binomial(link = "logit"))

car::Anova(predation_mod, type = 'III')

#Check assumptions
testResiduals(predation_mod, plot = TRUE)

plot(residuals(predation_mod) ~
       predict(predation_mod,type="link"),xlab=expression(hat(eta)),
     ylab="Deviance residuals",pch=20,col="blue")

performance::check_overdispersion(predation_mod) 

#Predict feeding rates and plot
newdat <- data.frame(
  Temp = seq(min(df$Temp),
                    max(df$Temp),
             length.out = 110),
             ID = df$ID[9],
             Exp_Day = df$Exp_Day[1]
)

pred <- predict(predation_mod,
                newdata = newdat,
                type = "response",
                se.fit = TRUE,
                re.form = NA)

newdat$fit <- pred$fit
newdat$lwr <- plogis(qlogis(pred$fit) - 1.96 * pred$se.fit)
newdat$upr <- plogis(qlogis(pred$fit) + 1.96 * pred$se.fit)

predicted_pred <- ggplot() +
  geom_jitter(data = df, aes(Temp, Dead / 20), alpha = 0.6, color = "black",
             width = 0.2) +
  geom_line(data = newdat, aes(x = Temp, y = fit), linewidth = 1,
            color = "cornflowerblue") +
  geom_ribbon(data = newdat,
              aes(x = Temp, ymin = lwr, ymax = upr),
              alpha = 0.2, fill = "cornflowerblue") +
  scale_x_continuous(limits = c(3, 15), breaks = seq(3, 15, by = 3)) +
  theme_minimal() +
  ylab("Predicted predation rate
       (Proportion eaten)") +
  xlab("Temperature") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black")
  )


ggsave("Predicted_Direct.png", plot = predicted_pred, 
       width = 6.5, height = 5.26, dpi = 450)
