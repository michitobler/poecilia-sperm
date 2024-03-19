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
library(effectsize)
library(plotly)


################################################################################################################################
#
#
#     Input data and filtering
#
#
#################################################################################################################################

raw.sperm.comp <- read_excel(path = "~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/02_Raw_Data/sperm_competition_thresholds_and_CASA_results_t0.5_t10.5.xlsx")
#discard.ids <- raw.sperm.comp$VideoID[which(is.na(raw.sperm.comp$AvgVCL))]
discard.ids <- raw.sperm.comp$VideoID[which(raw.sperm.comp$AvgVCL == "NaN")]
sperm.competition <- raw.sperm.comp[which(!(raw.sperm.comp$VideoID %in% discard.ids) & raw.sperm.comp$Sulfide != 250),] %>% dplyr::select(VideoID, Tank, Pop, Food, Treatment, Sulfide, MaleSL, Timepoint, PercentMotility, AvgVCL, AvgVAP, AvgVSL, AvgLIN, NumSpermTracked)

numerics <- c("PercentMotility", "MaleSL", "AvgVCL", "AvgVAP", "AvgVSL", "AvgLIN", "NumSpermTracked")
sperm.competition[,which(colnames(sperm.competition) %in% numerics)] <- lapply(sperm.competition[,which(colnames(sperm.competition) %in% numerics)],as.numeric)
factors <- c("VideoID", "Tank", "Pop", "Food", "Treatment", "Sulfide", "Timepoint")
sperm.competition[,which(colnames(sperm.competition) %in% factors)] <- lapply(sperm.competition[,which(colnames(sperm.competition) %in% factors)],as.factor)
#ggpairs(sperm.competition[,c("PercentMotility", "MaleSL", "AvgVAP", "AvgLIN", "NumSpermTracked")])

################################################################################################################################
#
#
#     Exploratory plotting with raw data
#
#
#################################################################################################################################

### VAP
# at 0 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 0),], aes(x = Timepoint, y = AvgVAP, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

# at 50 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 50),], aes(x = Timepoint, y = AvgVAP, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

### at 100 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 100),], aes(x = Timepoint, y = AvgVAP, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

# Across timepoints and sulfide concentrations
ggplot(data = sperm.competition, aes(x = Pop, y = AvgVAP, fill = Treatment)) + 
  theme_classic() + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + 
  ylim(c(25,125)) + 
  scale_fill_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"))

ggplot(data = sperm.competition, aes(x = Sulfide, y = AvgVAP, fill = Treatment)) + 
  theme_classic() + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + 
  ylim(c(25,125)) + 
  scale_fill_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"))


### Motility
# at 0 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 0),], aes(x = Timepoint, y = PercentMotility, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

# at 50 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 50),], aes(x = Timepoint, y = PercentMotility, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

### at 100 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 100),], aes(x = Timepoint, y = PercentMotility, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

ggplot(data = sperm.competition, aes(x = Sulfide, y = PercentMotility, fill = Treatment)) + 
  theme_classic() + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + 
  #ylim(c(0.8,0.95)) + 
  scale_fill_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"))


### Linearity
# at 0 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 0),], aes(x = Timepoint, y = AvgLIN, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

# at 50 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 50),], aes(x = Timepoint, y = AvgLIN, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

### at 100 sulfide
ggplot(data = sperm.competition[which(sperm.competition$Sulfide == 100),], aes(x = Timepoint, y = AvgLIN, color = Treatment)) + 
  theme_classic() + 
  stat_summary(fun = mean) + 
  stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
  stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
  scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")

ggplot(data = sperm.competition, aes(x = Pop, y = AvgLIN, fill = Treatment)) + 
  theme_classic() + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + 
  ylim(c(0.8,0.95)) + 
  scale_fill_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"))

ggplot(data = sperm.competition, aes(x = Sulfide, y = AvgLIN, fill = Treatment)) + 
  theme_classic() + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + 
  ylim(c(0.8,0.95)) + 
  scale_fill_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"))


################################################################################################################################
#
#
#     Modelling with raw data
#
#
#################################################################################################################################

### Raw data
sperm.competition.noNA <- sperm.competition[complete.cases(sperm.competition),]
#numSamp <- nrow(sperm.competition.noNA)

# ggplot(data = sperm.competition.noNA[which(sperm.competition.noNA$NumSpermTracked * sperm.competition.noNA$PercentMotility > 10),], aes(x = NumSpermTracked * PercentMotility, y = AvgVAP, color = Food)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# tdf <- sperm.competition.noNA
# tdf$NumMoving <- tdf$NumSpermTracked * tdf$PercentMotility
# summary(lm(AvgVAP ~ NumMoving, data = tdf))
# ggplot(data = tdf, aes(x = NumSpermTracked, y = NumMoving)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")
# ggplot(data = tdf, aes(x = asin(sqrt(PercentMotility)))) + 
#   geom_histogram()
# ggplot(data = tdf, aes(x = asin(sqrt(PercentMotility)))) + 
#   geom_density()
# ggplot(data = tdf, aes(x = PercentMotility)) + 
#   geom_density()
ggplot(data = sperm.competition.noNA[which(sperm.competition.noNA$Timepoint == 0.5),], aes(x = asin(sqrt(PercentMotility)), y = AvgVAP, color = interaction(Pop, Food))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition.noNA[which(sperm.competition.noNA$Timepoint == 10.5),], aes(x = asin(sqrt(PercentMotility)), y = AvgVAP, color = interaction(Pop, Food))) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(data = sperm.competition.noNA[which(sperm.competition.noNA$Timepoint == 0.5),], aes(x = asin(sqrt(PercentMotility)), y = asin(sqrt(AvgLIN)), color = interaction(Pop, Food))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition.noNA[which(sperm.competition.noNA$Timepoint == 10.5),], aes(x = asin(sqrt(PercentMotility)), y = asin(sqrt(AvgLIN)), color = interaction(Pop, Food))) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(data = sperm.competition.noNA, aes(x = Timepoint, y = asin(sqrt(PercentMotility)))) + 
  geom_boxplot()
ggplot(data = sperm.competition.noNA, aes(x = Timepoint, y = asin(sqrt(AvgLIN)))) + 
  geom_boxplot()
ggplot(data = sperm.competition.noNA, aes(x = Timepoint, y = AvgVAP)) + 
  geom_boxplot()

ggplot(data = sperm.competition.noNA, aes(x = AvgVAP, color = Pop)) + 
  geom_histogram()
ggplot(data = sperm.competition.noNA, aes(x = AvgVCL, color = Pop)) + 
  geom_histogram()
ggplot(data = sperm.competition.noNA, aes(x = AvgVSL, color = Pop)) + 
  geom_histogram()
ggplot(data = sperm.competition.noNA, aes(x = AvgLIN, color = Pop)) + 
  geom_histogram()
ggplot(data = sperm.competition.noNA, aes(x = PercentMotility, color = Pop)) + 
  geom_histogram()

plot_ly(sperm.competition.noNA, x = ~asin(sqrt(PercentMotility)), y = ~asin(sqrt(AvgLIN)), z = ~AvgVAP, color = ~interaction(Pop, Food))#, colors = c("steelblue", "#CB6E37"))

dv <- data.frame(trans.mot = asin(sqrt(sperm.competition.noNA$PercentMotility)), trans.lin = asin(sqrt(sperm.competition.noNA$AvgLIN)), vap = sperm.competition.noNA$AvgVAP)

ggpairs(dv)

sc.pca <- prcomp(dv, center = T, scale. = T)
sc.pca
summary(sc.pca)
test <- sperm.competition.noNA
test$PC1 <- sc.pca$x[,1]
test$PC2 <- sc.pca$x[,2]
test$PC3 <- sc.pca$x[,3]

ggplot(data = test[which(test$Timepoint == 0.5),], aes(x = PC1, y = PC2, color = interaction(Pop, Food))) + 
  geom_point()
ggplot(data = test[which(test$Timepoint == 10.5),], aes(x = PC1, y = PC2, color = interaction(Pop, Food))) + 
  geom_point()
ggplot(data = test, aes(x = Timepoint, y = PC1)) + 
  geom_boxplot()
ggplot(data = test, aes(x = Timepoint, y = PC1, fill = interaction(Pop, Food))) + 
  geom_boxplot()
ggplot(data = test, aes(x = interaction(Pop, Food), y = PC1, fill = Timepoint)) + 
  geom_boxplot()
ggplot(data = test, aes(x = interaction(Pop, Food), y = PC2, fill = Timepoint)) + 
  geom_boxplot()
ggplot(data = test, aes(x = interaction(Pop, Food), y = PC3, fill = Timepoint)) + 
  geom_boxplot()
ggplot(data = test, aes(x = PC1, y = PC2, color = Food)) + 
  geom_point()
ggplot(data = test, aes(x = PC2, y = PC3, color = interaction(Pop, Food))) + 
  geom_point()
ggplot(data = test, aes(x = PC1, y = PC3, color = interaction(Pop, Food))) + 
  geom_point()
ggplot(data = test[which(test$Timepoint == 0.5),], aes(x = MaleSL, y = PC1, color = interaction(Pop, Food))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = test[which(test$Timepoint == 10.5),], aes(x = MaleSL, y = PC1, color = interaction(Pop, Food))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = test, aes(x = MaleSL, y = PC1, color = Timepoint)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = test, aes(x = MaleSL, y = PC1)) + 
  geom_point() + 
  geom_smooth(method = "lm")

glob.pc1 <- lmer(PC1 ~ Pop * Food * Sulfide + Timepoint + MaleSL + (1 | Tank), na.action = "na.fail", data = test)
dredge.pc1 <- dredge(glob.pc1, rank = "AICc", m.lim = c(1,4))
print(dredge.pc1)
best.pc1.mod <- lmer(PC1 ~ Food + (1 | Tank), data = test)
Anova(best.pc1.mod)
plot(Effect("Food", best.pc1.mod))
eta_squared(best.pc1.mod, ci = .95)


glob.pc2 <- lmer(PC2 ~ Pop * Food * Sulfide + Timepoint + MaleSL + (1 | Tank), na.action = "na.fail", data = test)
dredge.pc2 <- dredge(glob.pc2, rank = "AICc", m.lim = c(1,4))
print(dredge.pc2)
best.pc2.mod <- lmer(PC2 ~ Pop + (1 | Tank), data = test)
Anova(best.pc2.mod, type = "III")
plot(Effect("Pop", best.pc2.mod))
eta_squared(best.pc2.mod, ci = .95)


glob.pc3 <- lmer(PC3 ~ Pop * Food * Sulfide + Timepoint + MaleSL + (1 | Tank), na.action = "na.fail", data = test)
dredge.pc3 <- dredge(glob.pc3, rank = "AICc", m.lim = c(1,4))
print(dredge.pc3)
best.pc3.mod <- lmer(PC3 ~ Food + (1 | Tank), data = test)
Anova(best.pc3.mod, type = "III")
plot(Effect("Food", best.pc3.mod))
eta_squared(best.pc3.mod, ci = .95)

# Motility
mot.global.lmer <- lmer(asin(sqrt(PercentMotility)) ~ Pop * Food * Sulfide * Timepoint + MaleSL + (1 | Tank), na.action = "na.fail", data = sperm.competition.noNA)
#mot.global.glmer <- glmer(PercentMotility ~ (Pop + Food + Sulfide + Timepoint)^2 + MaleSL + (1 | Tank) + (1 | VideoID), family = "binomial", na.action = "na.fail" , data = sperm.competition.noNA)
mot.dredge <- dredge(global.model = mot.global.lmer, rank = "AICc", m.lim = c(1,4))
print(mot.dredge, abbrev.names = TRUE)
#mot.dredge <- dredge(global.model = mot.global.lmer, rank = "AICc", m.lim = c(1,4))
#print(mot.dredge, abbrev.names = TRUE)

#mot.avg.out <- model.avg(get.models(mot.dredge, subset = delta <= 2))
#summary(mot.avg.out)

# ggplot(sperm.competition.noNA, aes(MaleSL, PercentMotility, color=Sulfide))+
#   geom_point()+
#   geom_smooth(method="lm")
# 
# ggplot(sperm.competition.noNA, aes(Sulfide, PercentMotility, color=Pop))+
#   geom_boxplot()

#mot.glmer <- glmer((PercentMotility * NumSpermTracked) / NumSpermTracked~ Food + Sulfide + Timepoint + MaleSL + (1 | Tank) + (1 | VideoID), family = "binomial", na.action = "na.fail", weights = NumSpermTracked, data = sperm.competition.noNA)
# mot.glmer <- glmer(PercentMotility ~ Sulfide + MaleSL + (1 | Tank) + (1 | VideoID), family = "binomial", na.action = "na.fail", data = sperm.competition.noNA)
# mot.glmer <- glmer(PercentMotility ~ Food + Sulfide + Timepoint + MaleSL + (1 | Tank) + (1 | VideoID), family = "binomial", na.action = "na.fail", data = sperm.competition.noNA)
best.mot.lmer <- lmer(asin(sqrt(PercentMotility)) ~ Timepoint + (1 | Tank), na.action = "na.fail", data = sperm.competition.noNA)
#mot.lmer <- lmer(asin(sqrt(PercentMotility)) ~ Pop * Food * Sulfide + MaleSL + (1 | Tank) + (1 | VideoID), na.action = "na.fail", data = sperm.competition.noNA)
#mot.lmer <- lmer(asin(sqrt(PercentMotility)) ~ Pop * Food * Sulfide + Timepoint +  MaleSL + (1 | Tank), na.action = "na.fail", data = sperm.competition.noNA)
#mot.lmer <- lmer(PercentMotility ~ Pop * Food * Sulfide + MaleSL + (1 | Tank) + (1 | VideoID), na.action = "na.fail", data = sperm.competition.noNA)
#mot.glmer <- glmer(PercentMotility ~ Sulfide + MaleSL + (1 | Tank) + (1 | VideoID), family = "binomial", na.action = "na.fail" , data = sperm.competition.noNA)
summary(best.mot.lmer)
Anova(best.mot.lmer, type = "III")
eta_squared(best.mot.lmer, ci = .95)

ggplot(data = sperm.competition.noNA, aes(x = Timepoint, y = asin(sqrt(PercentMotility)))) + 
  geom_boxplot()

plot(Effect("Timepoint", best.mot.lmer))



# mot.lmer.select <- lmer(asin(sqrt(PercentMotility)) ~ Pop + Sulfide + Timepoint + Pop*Sulfide + (1 | Tank) + (1 | VideoID), na.action = "na.fail", data = sperm.competition.noNA)
# mot.lmer.select.fit <- Anova(mot.lmer.select, type = "III")
# mot.lmer.select.sum <- summary(mot.lmer.select.fit)
# dredge(mot.lmer.select)
# eta_squared(mot.lmer.select, ci = .95)

plot(Effect("Timepoint", best.mot.lmer))

ggplot(data = sperm.competition.noNA, aes(x = MaleSL, y = PercentMotility, color = Sulfide)) + 
  geom_point() + 
  scale_color_manual(values = c("#FFFDA6", "#FFD700", "#FFB90F")) + 
  geom_smooth(aes(group = Sulfide), method = "lm", se = F)
ggplot(data = sperm.competition.noNA[which(sperm.competition.noNA$Timepoint == 10.5),], aes(x = MaleSL, y = PercentMotility, color = interaction(Pop, Food))) + 
  geom_point() + 
  geom_smooth(aes(group = interaction(Pop, Food)), method = "lm", se = F)
ggplot(data = sperm.competition.noNA, aes(x = Pop, y = PercentMotility, fill = Pop)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = c("steelblue", "gold"))
ggplot(data = sperm.competition.noNA, aes(x = Sulfide, y = PercentMotility, fill = Sulfide)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = c("#FFFDA6", "#FFD700", "#FFB90F"))
ggplot(data = sperm.competition.noNA, aes(x = Timepoint, y = PercentMotility, fill = Timepoint)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = c("lightgreen", "red"))
ggplot(data = sperm.competition.noNA, aes(x = Food, y = PercentMotility, fill = Food)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = c("tan", "grey"))
# ggplot(data = sperm.competition.noNA, aes(x = Pop, y = PercentMotility, fill = Food)) + 
#   geom_boxplot() + 
#   theme_classic() + 
#   scale_fill_manual(values = c("tan", "grey"))
ggplot(data = sperm.competition.noNA, aes(x = interaction(Food, Pop), y = PercentMotility, fill = interaction(Food, Pop))) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_fill_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"))

ggplot(data = sperm.competition.noNA, aes(x = Sulfide, y = PercentMotility, fill = Sulfide)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#FFFDA6", "#FFD700", "#FFB90F"))

ggplot(data = sperm.competition, aes(x = Sulfide, y = PercentMotility, fill = Sulfide)) + 
  theme_classic() + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + 
  #ylim(c(0.8,0.95)) + 
  scale_fill_manual(values = c("#FFFDA6", "#FFD700", "#FFB90F"))

sl.mot.effect <- plot(Effect("MaleSL",mot.glmer))
sulf.mot.effect <- plot(Effect("Sulfide",mot.glmer))
sl.mot.effect.data <- data.frame(sl.x = sl.mot.effect[["panel.args"]][[1]]$x, sl.y = sl.mot.effect[["panel.args"]][[1]]$y)
sulf.mot.effect.data <- data.frame(sulf.x = sulf.mot.effect[["panel.args"]][[1]]$x, sulf.y = sulf.mot.effect[["panel.args"]][[1]]$y)
plot(Effect(c("MaleSL", "Sulfide"),mot.glmer))


# VAP
vap.global.lmer <- lmer(AvgVAP ~ Pop * Food * Sulfide + Timepoint + MaleSL + (1 | Tank), na.action = "na.fail" , data = sperm.competition.noNA)
vap.dredge <- dredge(global.model = vap.global.lmer, rank = "AICc", m.lim = c(1,4))
print(vap.dredge, abbrev.names = TRUE)
vap.lmer <- lmer(AvgVAP ~ Food + Pop + Sulfide + Food*Sulfide + (1 | Tank), na.action = "na.fail" , data = sperm.competition.noNA)
#vap.lme <- lme(AvgVAP ~ Pop*Sulfide + Food, random = list(~ 1 | Tank, ~1 | VideoID), data = sperm.competition.noNA)
#vap.lmer <- lmer(AvgVAP ~ Food*Sulfide + Pop + (1 | Tank), na.action = "na.fail" , data = sperm.competition.noNA, REML=TRUE)
summary(vap.lmer)
Anova(vap.lmer, type = "III")
confint(vap.lmer)
eta_squared(vap.lmer, ci = .95)

plot(Effect("Pop", vap.lmer))
plot(Effect("Food", vap.lmer))
plot(Effect("Sulfide", vap.lmer))
plot(Effect(c("Food", "Sulfide"), vap.lmer))


# plot(Effect(c("Pop", "Food", "Sulfide"), vap.lmer))
# plot(Effect(c("Pop", "Food", "Sulfide"), vap.lme))
# plot(Effect(c("Food", "Sulfide"), vap.lmer))
# plot(Effect(c("Food", "Sulfide"), vap.lme))
# plot(Effect(c("Pop", "Sulfide"), vap.lmer))
# plot(Effect(c("Pop", "Sulfide"), vap.lme))
# plot(Effect(c("Food", "Pop"), vap.lmer))
# plot(Effect(c("Food", "Pop"), vap.lme))




# Linearity
lin.global.lmer <- lmer(asin(sqrt(AvgLIN)) ~ Pop * Food * Sulfide + Timepoint + MaleSL + (1 | Tank), na.action = "na.fail" , data = sperm.competition.noNA)
lin.dredge <- dredge(global.model = lin.global.lmer, rank = "AICc", m.lim = c(1,4))
print(lin.dredge, abbrev.names = TRUE)

best.lin.mod <- lmer(asin(sqrt(AvgLIN)) ~ Food + (1 | Tank), na.action = "na.fail" , data = sperm.competition.noNA)
Anova(best.lin.mod, type = "III")
confint(best.lin.mod)
eta_squared(best.lin.mod)
plot(Effect("Food", best.lin.mod))

ggplot(data = sperm.competition.noNA, aes(x = Food, y = asin(sqrt(AvgLIN)))) + 
  geom_boxplot()
# lin.glmer.sl<- glmer(AvgLIN ~ MaleSL + (1 | Tank) + (1 | VideoID), na.action = "na.fail" , family = "binomial", data = sperm.competition.noNA)
# plot(Effect("MaleSL", lin.glmer.sl))
# lin.glmer.pop <- glmer(AvgLIN ~ Pop + (1 | Tank) + (1 | VideoID), na.action = "na.fail" , family = "binomial", data = sperm.competition.noNA)
# plot(Effect("Pop", lin.glmer.pop))
# lin.glmer.time<- glmer(AvgLIN ~ Timepoint + (1 | Tank) + (1 | VideoID), na.action = "na.fail" , family = "binomial", data = sperm.competition.noNA)
# plot(Effect("Timepoint", lin.glmer.time))


# lin.avg.out <- model.avg(lin.dredge)
# summary(lin.avg.out)
# summary(lin.glmer.time)

ggplot(data = sperm.competition.noNA, aes(x = interaction(Pop, Food), y = AvgLIN, fill = interaction(Pop, Food))) + 
  geom_boxplot() + 
  theme_classic()




################################################################################################################################
#
#
#     Modelling with data normalized to t=0 and zero sulfide within each male.
#
#
#################################################################################################################################

data <- sperm.competition.noNA
ids <- as.character(data$VideoID)
id.list <- str_split(ids, "_")
id.df <- data.frame(matrix(unlist(id.list), nrow = length(id.list), byrow=TRUE), stringsAsFactors=FALSE)
data$fishID <- id.df$X3 
data[-c(1:8, 10, 12, 14, 15)] <- as.data.frame(sapply(data[-c(1:8, 10, 12, 14, 15)], function(col) ave(col, data$fishID, FUN=function(x) x / x[1]))) # need to do some hard-coding here unfortunately......

# ### VAP
# # at 0 sulfide
# ggplot(data = data[which(data$Sulfide == 0),], aes(x = Timepoint, y = AvgVAP, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")
# 
# # at 50 sulfide
# ggplot(data = data[which(data$Sulfide == 50),], aes(x = Timepoint, y = AvgVAP, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")
# 
# # at 100 sulfide
# ggplot(data = data[which(data$Sulfide == 100),], aes(x = Timepoint, y = AvgVAP, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")
# 
# ### Motility
# # at 0 sulfide
# ggplot(data = data[which(data$Sulfide == 0),], aes(x = Timepoint, y = PercentMotility, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")
# 
# # at 50 sulfide
# ggplot(data = data[which(data$Sulfide == 50),], aes(x = Timepoint, y = PercentMotility, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")
# 
# # at 100 sulfide
# ggplot(data = data[which(data$Sulfide == 100),], aes(x = Timepoint, y = PercentMotility, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")
# 
# ### Linearity
# # at 0 sulfide
# ggplot(data = data[which(data$Sulfide == 0),], aes(x = Timepoint, y = AvgLIN, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")
# 
# # at 50 sulfide
# ggplot(data = data[which(data$Sulfide == 50),], aes(x = Timepoint, y = AvgLIN, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")
# 
# # at 100 sulfide
# ggplot(data = data[which(data$Sulfide == 100),], aes(x = Timepoint, y = AvgLIN, color = Treatment)) + 
#   theme_classic() + 
#   stat_summary(fun = mean) + 
#   stat_summary(aes(color = Treatment), fun.data = mean_cl_boot, geom = "errorbar", width = 0.1) +
#   stat_summary(aes(group = Treatment), geom = "line", fun = mean) + 
#   scale_color_manual(values = c("steelblue", "deepskyblue2", "darkgoldenrod1", "gold"), name = "Treatment")




mot.norm.global.lmer <- lmer(PercentMotility ~ Pop * Food * Sulfide + Timepoint + MaleSL + (1 | Tank), na.action = "na.fail" , data = data)
mot.norm.dredge <- dredge(global.model = mot.norm.global.lmer, rank = "AICc", m.lim = c(1,4))
print(mot.norm.dredge, abbrev.names = TRUE)

mot.norm.lmer <- lmer(PercentMotility ~ Timepoint + (1 | Tank), na.action = "na.fail", data = data)
#mot.norm.lme <- lme(PercentMotility ~ Food + Pop + Timepoint, random = list(~1 | Tank, ~1 | VideoID), data = data)
#summary(mot.norm.lme)
Anova(mot.norm.lmer, type = "III")
plot(Effect("Timepoint", mot.norm.lmer))

ggplot(data = data, aes(x = Timepoint, y = PercentMotility)) + 
  geom

# VAP
vap.norm.global.lmer <- lmer(AvgVAP ~ Pop * Food * Sulfide + Timepoint + MaleSL + (1 | Tank), data = data, na.action = "na.fail")
#vap.norm.global.lme <- lme(AvgVAP ~ (Pop + Food + Sulfide + Timepoint)^2 + MaleSL, random = list(~1 | Tank, ~1 | VideoID), data = data)
vap.norm.dredge <- dredge(global.model = vap.norm.global.lmer, rank = "AICc", m.lim = c(1,4))
print(vap.norm.dredge, abbrev.names = TRUE)
best.norm.vap.mod <- lmer(AvgVAP ~ Pop + (1 | Tank), data = data, na.action = "na.fail")
Anova(best.norm.vap.mod, type = "III")
plot(Effect("Pop", best.norm.vap.mod))
confint(best.norm.vap.mod)
eta_squared(best.norm.vap.mod)
#vap.norm.lme <- lme(AvgVAP ~ Pop, random = list(~1 | Tank, ~1 | VideoID), data = data)
#vap.lme <- lme(AvgVAP ~ Food*Sulfide + Pop, random = list(~ 1 | Tank, ~1 | VideoID), data = sperm.competition.noNA)
#vap.lmer <- lmer(AvgVAP ~ Food*Sulfide + Pop + (1 | Tank), na.action = "na.fail" , data = sperm.competition.noNA, REML=TRUE)
summary(vap.norm.lme)

plot(Effect("Pop", vap.norm.lme))





# Linearity
lin.norm.global.lmer <- lmer(AvgLIN ~ Pop * Food * Sulfide + Timepoint + MaleSL + (1 | Tank), na.action = "na.fail", data = data)
#lin.norm.global.lme <- lme(AvgLIN ~ (Pop + Food + Sulfide + Timepoint)^2 + MaleSL, random = list(~1 | Tank, ~1 | VideoID), data = data)
lin.norm.dredge <- dredge(global.model = lin.norm.global.lmer, rank = "AICc", m.lim = c(1,4))
print(lin.norm.dredge, abbrev.names = TRUE)

best.lin.norm.mod <- lmer(AvgLIN ~ Pop + (1 | Tank), na.action = "na.fail", data = data)
Anova(best.lin.norm.mod, type = "III")
confint(best.lin.norm.mod)
plot(Effect("Pop", best.lin.norm.mod))
eta_squared(best.lin.norm.mod, ci = .95)




