rm(list = ls())

library(ggplot2)
library(car)
library(readxl)
library(GGally)
library(lme4)
library(bbmle)
library(dplyr)
library(tidyr)
library(nlme)
library(MuMIn)
library(stringr)
library(effects)
library(heplots)


################################################################################################################################
#
#
#     Input data and filtering
#
#
#################################################################################################################################

raw.ejac <- read_xlsx(path = "~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/02_Raw_Data/ejaculate_size_data.xlsx")
raw.ejac$SL.mm <- as.numeric(raw.ejac$SL.mm)
raw.ejac$AreaEjaculate.pixels.ROUNDDOWN <- as.numeric(raw.ejac$AreaEjaculate.pixels.ROUNDDOWN)
raw.ejac$pop <- as.factor(raw.ejac$pop)
raw.ejac$food <- as.factor(raw.ejac$food)
raw.ejac$tank <- as.factor(raw.ejac$tank)
raw.ejac$DateFoodStart <- as.Date(raw.ejac$DateFoodStart, format = "%m/%d/%y")
raw.ejac$TrialDate <- as.Date(raw.ejac$TrialDate, format = "%m/%d/%y")
raw.ejac$TimeFoodTreat <- as.numeric(raw.ejac$TrialDate - raw.ejac$DateFoodStart)
#ejac.complete.cases <- raw.ejac[1:32,] # this will need to be changed as more data is added

# ggplot(data = raw.ejac, aes(x = interaction(pop,food), y = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
ggplot(data = raw.ejac, aes(x = pop, y = AreaEjaculate.pixels.ROUNDDOWN, fill = food)) +
  geom_boxplot() +
  theme_classic()
# 
# ggplot(data = raw.ejac, aes(x = pop, y = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
# ggplot(data = raw.ejac, aes(x = food, y = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
# ggplot(data = raw.ejac, aes(x = SL.mm, y = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_point() + 
#   theme_classic()
# 
ejac.size.plot <- ggplot() +
  geom_boxplot(data = raw.ejac, aes(x = pop, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10), color = food), outlier.alpha = 0) +
  geom_point(data = raw.ejac, aes(x = pop, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10), color = food), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("Ejaculate Size") + 
  theme_classic()
ejac.size.plot

pdf("~/Downloads/FigureS1_sperm_comp.pdf", width = 4, height = 4)
ejac.size.plot
dev.off()



# Calculate sample sizes:
sum(interaction(raw.ejac$pop, raw.ejac$food) == "bon.foodH")
sum(interaction(raw.ejac$pop, raw.ejac$food) == "pso.foodH")
sum(interaction(raw.ejac$pop, raw.ejac$food) == "bon.foodL")
sum(interaction(raw.ejac$pop, raw.ejac$food) == "pso.foodL")

# ggplot(data = raw.ejac, aes(x = pop, y = AreaEjaculate.pixels.ROUNDDOWN / SL.mm, fill = food)) + 
#   geom_boxplot() + 
#   theme_classic()


#lmer(lmer(AreaEjaculate.pixels.ROUNDDOWN ~ pop*food + SL.mm +asdf , na.action = "na.fail" , data = data))
# ejac.mod <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop*food, data = raw.ejac)
# Anova(ejac.mod, type = "III")

# ejac.mod1 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop, data = raw.ejac)
# ejac.mod2 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ food, data = raw.ejac)
# ejac.mod3 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop + food, data = raw.ejac)
# ejac.mod4 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ SL.mm, data = raw.ejac)
# ejac.mod5 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop + SL.mm, data = raw.ejac)
# ejac.mod6 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ food + SL.mm, data = raw.ejac)
# ejac.mod7 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop + food + SL.mm, data = raw.ejac)
# ejac.mod8 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop*food, data = raw.ejac)
# ejac.mod9 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop*food + SL.mm, data = raw.ejac)
sum(interaction(raw.ejac$pop, raw.ejac$food)=="pso.foodH")
sum(interaction(raw.ejac$pop, raw.ejac$food)=="pso.foodL")
sum(interaction(raw.ejac$pop, raw.ejac$food)=="bon.foodH")
sum(interaction(raw.ejac$pop, raw.ejac$food)=="bon.foodL")


glob.mod <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop*food + SL.mm + TimeFoodTreat, data = raw.ejac, na.action = "na.fail")
glob.mod.dredge <- dredge(glob.mod, rank = "AICc", m.lim = c(0, 4))
print(glob.mod.dredge, abbrev.names = TRUE)
ejac.avg <- model.avg(glob.mod.dredge, subset = cumsum(weight) <= .95)
summary(ejac.avg)
confint(ejac.avg)


raw.ejac$std.ejac <- raw.ejac$AreaEjaculate.pixels.ROUNDDOWN / log(raw.ejac$SL.mm,10)
ggplot(data = raw.ejac, aes(x = std.ejac)) + 
  geom_histogram()
ggplot(data = raw.ejac, aes(x = log(SL.mm, 10), y = std.ejac / 500000)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_classic()

ggplot(data = raw.ejac, aes(x = pop, y = std.ejac)) + 
  geom_boxplot() + 
  theme_classic()
ggplot(data = raw.ejac, aes(x = food, y = std.ejac)) + 
  geom_boxplot() + 
  theme_classic()

raw.ejac$z.ejac <- (raw.ejac$std.ejac - mean(raw.ejac$std.ejac)) / sd(raw.ejac$std.ejac)
ggplot(data = raw.ejac, aes(x = log(SL.mm, 10), y = z.ejac)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  ylab("Size-Corrected Ejaculate Size (zscore)") + 
  theme_classic()
ggplot(data = raw.ejac, aes(x = pop, y = std.ejac)) + 
  geom_boxplot()
ggplot(data = raw.ejac, aes(x = food, y = std.ejac)) + 
  geom_boxplot()

glob.std.mod <- lm(std.ejac ~ pop*food + log(SL.mm,10) + TimeFoodTreat, data = raw.ejac, na.action = "na.fail")
glob.std.mod.dredge <- dredge(glob.std.mod, rank = "AICc", m.lim = c(0, 4))
print(glob.std.mod.dredge, abbrev.names = TRUE)
ejac.std.avg <- model.avg(glob.std.mod.dredge, subset = cumsum(weight) <= .95)
summary(ejac.std.avg)
confint(ejac.std.avg)

# ggplot(data = raw.ejac, aes(x = pop, y = AreaEjaculate.pixels.ROUNDDOWN, fill = food)) + 
#   geom_boxplot()
# ggplot(data = raw.ejac, aes(x = pop, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10), fill = food)) + 
#   geom_boxplot()
# ggplot(data = raw.ejac, aes(x = AreaEjaculate.pixels.ROUNDDOWN, color = interaction(pop, food))) + 
#   geom_histogram()
# ggplot(data = raw.ejac, aes(x = log(AreaEjaculate.pixels.ROUNDDOWN, 10), color = interaction(pop,food))) + 
#   geom_histogram()


#write.csv(print(dredge.out, abbrev.names = TRUE), "~/Downloads/explore.mod.select.csv")
# m.avg.out <- model.avg(dredge.out)
# m.avg.out
# summary(m.avg.out)

#interaction is non-significant, so using additive model instead
best.ejac.mod <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop + SL.mm, data = raw.ejac, na.action = "na.fail")
etasq(best.ejac.mod, type = "III", partial = T, anova = T)
#best.explore.mod <- lmer(PC1 ~ pop * food.treat + (1 | tank), data = new.explore, na.action = "na.fail")
Anova(best.ejac.mod, type = "III")
confint(best.ejac.mod)

# ejac.mod1 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop, data = raw.ejac)
# ejac.mod2 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ food, data = raw.ejac)
# ejac.mod3 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop + food, data = raw.ejac)
# ejac.mod4 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ SL.mm, data = raw.ejac)
# ejac.mod5 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop + SL.mm, data = raw.ejac)
# ejac.mod6 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ food + SL.mm, data = raw.ejac)
# ejac.mod7 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop + food + SL.mm, data = raw.ejac)
# ejac.mod8 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop*food, data = raw.ejac)
# ejac.mod9 <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop*food + SL.mm, data = raw.ejac)

model.sel(ejac.mod1, ejac.mod2, ejac.mod3, ejac.mod4, ejac.mod5, ejac.mod6, ejac.mod7, ejac.mod8, ejac.mod9)

plot(Effect("SL.mm", ejac.mod5))
plot(Effect("pop", ejac.mod5))

sl.effect <- plot(Effect("SL.mm", ejac.mod5))
pop.effect <- plot(Effect("pop", ejac.mod5))

sl.effect.data <- data.frame(sl.x = sl.effect[["panel.args"]][[1]]$x, sl.y = sl.effect[["panel.args"]][[1]]$y)
pop.effect.data <- data.frame(pop.x = pop.effect[["panel.args"]][[1]]$x, pop.y = pop.effect[["panel.args"]][[1]]$y)
Anova(ejac.mod5, type = "III")
# Anova(ejac.mod4, type = "III")
# Anova(ejac.mod6, type = "III")
# Anova(ejac.mod7, type = "III")
# plot(Effect(c("pop", "SL.mm"), ejac.mod5))
# ejac.global.lm <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop*food + SL.mm, data = raw.ejac[which(raw.ejac$SL.mm != "NA"),], na.action = "na.fail")
# ejac.global.lm <- lm(PercentTotalAreaSpermBundles ~ pop*food + SL.mm, data = raw.ejac[which(raw.ejac$SL.mm != "NA"),], na.action = "na.fail")
# Anova(ejac.global.lm)
# ejac.dredge <- dredge(global.model = ejac.global.lm, rank = "AICc", m.max = 4)
# print(mot.dredge, abbrev.names = TRUE)
# ggplot(data = raw.ejac, aes(x = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_density()
# ggplot(data = raw.ejac, aes(x = log(AreaEjaculate.pixels.ROUNDDOWN, 10))) + 
#   geom_density()

# ggplot(data = raw.ejac, aes(x = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_histogram(bins = )
# ggplot(data = raw.ejac, aes(x = log(AreaEjaculate.pixels.ROUNDDOWN, 10))) + 
#   geom_histogram()

ggplot(data = raw.ejac, aes(x = SL.mm, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10))) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_classic()

ggplot(data = raw.ejac, aes(x = pop, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10), fill = food)) + 
  geom_boxplot() + 
  theme_classic()

ggplot(data = raw.ejac, aes(x = food, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10), fill = food)) + 
  geom_boxplot() + 
  theme_classic()

ggplot(data = raw.ejac, aes(x = pop, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10), fill = food)) + 
  geom_boxplot() + 
  theme_classic()

# ggplot(data = raw.ejac, aes(x = SL.mm, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10), color = pop)) + 
#   geom_point() + 
#   scale_color_manual(values = c("steelblue", "gold")) + 
#   geom_segment(aes(x = sl.effect.data$sl.x[1], xend = sl.effect.data$sl.x[5], y = sl.effect.data$sl.y[1], yend = sl.effect.data$sl.y[5])) + 
#   theme_classic()
ggplot(data = raw.ejac, aes(x = SL.mm, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10))) + 
  geom_point(aes(fill = pop), size = 7, pch = 21) + 
  scale_fill_manual(values = c("steelblue", "gold")) + 
  geom_segment(aes(x = sl.effect.data$sl.x[1], xend = sl.effect.data$sl.x[5], y = sl.effect.data$sl.y[1], yend = sl.effect.data$sl.y[5]), size = 2) + 
  theme_classic()

# ggplot(data = raw.ejac, aes(x = TimeFoodTreat, y = log(AreaEjaculate.pixels.ROUNDDOWN, 10), color = pop)) + 
#   geom_point() + 
#   theme_classic()


# ggplot(data = ejac.complete.cases, aes(x = interaction(pop,food), y = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
# ggplot(data = ejac.complete.cases, aes(x = pop, y = AreaEjaculate.pixels.ROUNDDOWN, fill = food)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
# ggplot(data = ejac.complete.cases, aes(x = pop, y = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
# ggplot(data = ejac.complete.cases, aes(x = food, y = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
# ggplot(data = ejac.complete.cases, aes(x = SL.mm, y = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_point() + 
#   theme_classic()
# 
# # Calculate sample sizes:
# sum(interaction(ejac.complete.cases$pop, ejac.complete.cases$food) == "bon.foodH")
# sum(interaction(ejac.complete.cases$pop, ejac.complete.cases$food) == "pso.foodH")
# sum(interaction(ejac.complete.cases$pop, ejac.complete.cases$food) == "bon.foodL")
# sum(interaction(ejac.complete.cases$pop, ejac.complete.cases$food) == "pso.foodL")
# 
# ggplot(data = ejac.complete.cases, aes(x = pop, y = AreaEjaculate.pixels.ROUNDDOWN / SL.mm, fill = food)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
# 
# #lmer(lmer(AreaEjaculate.pixels.ROUNDDOWN ~ pop*food + SL.mm +asdf , na.action = "na.fail" , data = data))
# ejac.mod <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop*food, data = ejac.complete.cases)
# Anova(ejac.mod, type = "III")
# 
# ejac.mod1 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop, data = ejac.complete.cases)
# ejac.mod2 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ food, data = ejac.complete.cases)
# ejac.mod3 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop + food, data = ejac.complete.cases)
# ejac.mod4 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ SL.mm, data = ejac.complete.cases)
# ejac.mod5 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop + SL.mm, data = ejac.complete.cases)
# ejac.mod6 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ food + SL.mm, data = ejac.complete.cases)
# ejac.mod7 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop + food + SL.mm, data = ejac.complete.cases)
# ejac.mod8 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop*food, data = ejac.complete.cases)
# ejac.mod9 <- lm(log(AreaEjaculate.pixels.ROUNDDOWN,10) ~ pop*food + SL.mm, data = ejac.complete.cases)
# 
# model.sel(ejac.mod1, ejac.mod2, ejac.mod3, ejac.mod4, ejac.mod5, ejac.mod6, ejac.mod7, ejac.mod8, ejac.mod9)
# 
# plot(Effect("SL.mm", ejac.mod5))
# plot(Effect("pop", ejac.mod5))
# Anova(ejac.mod5, type = "III")
# # ejac.global.lm <- lm(AreaEjaculate.pixels.ROUNDDOWN ~ pop*food + SL.mm, data = ejac.complete.cases[which(ejac.complete.cases$SL.mm != "NA"),], na.action = "na.fail")
# # ejac.global.lm <- lm(PercentTotalAreaSpermBundles ~ pop*food + SL.mm, data = ejac.complete.cases[which(ejac.complete.cases$SL.mm != "NA"),], na.action = "na.fail")
# # Anova(ejac.global.lm)
# # ejac.dredge <- dredge(global.model = ejac.global.lm, rank = "AICc", m.max = 4)
# # print(mot.dredge, abbrev.names = TRUE)
# ggplot(data = ejac.complete.cases, aes(x = AreaEjaculate.pixels.ROUNDDOWN)) + 
#   geom_density()











