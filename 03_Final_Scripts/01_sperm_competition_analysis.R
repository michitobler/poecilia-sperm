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
library(EnvStats)


################################################################################################################################
#
#
#     Input data and filtering
#
#
#################################################################################################################################

raw.sperm.comp <- read_excel(path = "~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/02_Raw_Data/sperm_competition_taco_puya_sulfide_gradient.xlsx")
#discard.ids <- raw.sperm.comp$VideoID[which(is.na(raw.sperm.comp$AvgVCL))]
discard.ids <- raw.sperm.comp$VideoID[which(raw.sperm.comp$AvgVCL == "NaN")]
sperm.competition <- raw.sperm.comp[which(!(raw.sperm.comp$VideoID %in% discard.ids) & raw.sperm.comp$Sulfide != 250),] %>% dplyr::select(VideoID, Timepoint, Tank, IndivID, DateFoodStart, TrialDateTime, Drainage, Sulfidic, Pop, Food, Treatment, Sulfide, MaleSL, PercentMotility, AvgVCL, AvgVAP, AvgVSL, AvgLIN, NumSpermTracked, NumSpermMotile)
raw.cfc <- read_excel(path = "~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/02_Raw_Data/cryptic_female_choice_metadata.xlsx")[,1:26]
colnames(raw.cfc)[17:26] <- c("PercentMotility", "AvgVCL", "AvgVAP", "AvgVSL", "AvgLIN", "AvgWOB", "AvgPROG", "AvgBCF", "NumSpermTracked", "NumSpermMotile")
numerics <- c("PercentMotility", "MaleSL", "AvgVCL", "AvgVAP", "AvgVSL", "AvgLIN", "NumSpermTracked", "NumSpermMotile")
sperm.competition[,which(colnames(sperm.competition) %in% numerics)] <- lapply(sperm.competition[,which(colnames(sperm.competition) %in% numerics)],as.numeric)
raw.cfc[,which(colnames(raw.cfc) %in% numerics)] <- lapply(raw.cfc[,which(colnames(raw.cfc) %in% numerics)],as.numeric)
factors <- c("VideoID", "Tank", "Pop", "Food", "Treatment", "Sulfide", "Timepoint", "Sulfidic", "Drainage")
sperm.competition[,which(colnames(sperm.competition) %in% factors)] <- lapply(sperm.competition[,which(colnames(sperm.competition) %in% factors)],as.factor)
raw.cfc[,which(colnames(raw.cfc) %in% factors)] <- lapply(raw.cfc[,which(colnames(raw.cfc) %in% factors)],as.factor)
sperm.competition$d.t.food.start <- strptime(sperm.competition$DateFoodStart, format = "%m/%d/%y %H:%M")
sperm.competition$d.t.trial <- strptime(sperm.competition$TrialDateTime, format = "%m/%d/%y %H:%M")
sperm.competition$treat.length <- sperm.competition$d.t.trial - sperm.competition$d.t.food.start
sperm.competition$treat.length <- as.numeric(sperm.competition$treat.length)

zero.sulf.sc <- sperm.competition[which(sperm.competition$Sulfide == 0), c(1,2,4,6,7,8,9,13:20)]
colnames(zero.sulf.sc)[1] <- "TrialID"

cfc.sub <- raw.cfc[,c(1,15,2,8,9,11,13,6,17,18,19,20,21,25,26)]
colnames(cfc.sub)[c(3,4,5,6,7)] <- c("IndivID", "TrialDateTime", "Drainage", "Sulfidic", "Pop")
sc.cfc.mix <- as.data.frame(rbind(zero.sulf.sc, cfc.sub))
ggplot(data = sc.cfc.mix, aes(x = Drainage, y = AvgVSL, color = Sulfidic)) + 
  geom_boxplot()

################################################################################################################################
#
#
#     Exploratory Plotting
#
#
#################################################################################################################################

#### Plotting distributions of DVs
# VAP
ggplot(data = sperm.competition, aes(x = AvgVAP)) + 
  geom_histogram(binwidth = 1)
ggplot(data = sperm.competition, aes(x = AvgVAP, color = Drainage)) + 
  geom_density() 
ggplot(data = sperm.competition, aes(x = Drainage, y = AvgVAP, color = Sulfidic)) + 
  geom_boxplot() 
ggplot(data = sperm.competition, aes(x = Drainage, y = AvgVAP, color = Sulfidic, fill = Food)) + 
  geom_boxplot() 

# MOT
ggplot(data = sperm.competition, aes(x = asin(sqrt(PercentMotility)))) + 
  geom_histogram(binwidth = .01)
ggplot(data = sperm.competition, aes(x = asin(sqrt(PercentMotility)), color = Drainage)) + 
  geom_density() 
ggplot(data = sperm.competition, aes(x = Sulfidic, y = asin(sqrt(PercentMotility)), fill = Food)) + 
  geom_boxplot() 
ggplot(data = sperm.competition, aes(x = Drainage, y = asin(sqrt(PercentMotility)))) + 
  geom_boxplot() 
ggplot(data = sperm.competition, aes(x = Sulfidic, y = AvgVSL, fill = Food)) + 
  geom_boxplot() 

# LIN
ggplot(data = sperm.competition, aes(x = asin(sqrt(AvgLIN)))) + 
  geom_histogram(binwidth = .01)
ggplot(data = sperm.competition, aes(x = asin(sqrt(AvgLIN)), color = Drainage)) + 
  geom_density() 
ggplot(data = sperm.competition, aes(x = Drainage, y = asin(sqrt(AvgLIN)), color = Sulfidic)) + 
  geom_boxplot() 

# VAP vs MOT
ggplot(data = sperm.competition, aes(x = AvgVAP, y = asin(sqrt(PercentMotility)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# VAP vs LIN
ggplot(data = sperm.competition, aes(x = AvgVAP, y = asin(sqrt(AvgLIN)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# MOT vs LIN
ggplot(data = sperm.competition, aes(x = asin(sqrt(PercentMotility)), y = asin(sqrt(AvgLIN)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")



#### Now some plotting of covariates
# Effects of treatment length
ggplot(data = sperm.competition, aes(x = treat.length)) + 
  geom_histogram(binwidth = 1)
ggplot(data = sperm.competition, aes(x = treat.length, y = AvgVAP, color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition, aes(x = treat.length, y = asin(sqrt(PercentMotility)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition, aes(x = treat.length, y = asin(sqrt(AvgLIN)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# Effects of SL
ggplot(data = sperm.competition, aes(x = MaleSL)) + 
  geom_histogram(binwidth = 1)
ggplot(data = sperm.competition, aes(x = MaleSL, y = AvgVAP, color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition, aes(x = MaleSL, y = asin(sqrt(PercentMotility)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition, aes(x = MaleSL, y = asin(sqrt(PercentMotility)))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition, aes(x = MaleSL, y = asin(sqrt(AvgLIN)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# Effects of sperm numbers
ggplot(data = sperm.competition, aes(x = NumSpermMotile, color = Drainage)) + 
  geom_histogram(binwidth = 1)
ggplot(data = sperm.competition, aes(x = NumSpermMotile, y = AvgVAP, color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition[which(sperm.competition$NumSpermMotile >= 10),], aes(x = NumSpermMotile, y = AvgVAP, color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition, aes(x = NumSpermMotile, y = asin(sqrt(PercentMotility)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = sperm.competition, aes(x = NumSpermMotile, y = asin(sqrt(AvgLIN)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")


# make bins and calculate sd 
binned.motile <- c()
for (i in 1:nrow(sperm.competition)){
  if (sperm.competition$NumSpermMotile[i] <= 5){
    binned.motile[i] <- 1
  } else if (sperm.competition$NumSpermMotile[i] > 5 & sperm.competition$NumSpermMotile[i] <= 10){
    binned.motile[i] <- 2
  } else if (sperm.competition$NumSpermMotile[i] > 10 & sperm.competition$NumSpermMotile[i] <= 15){
    binned.motile[i] <- 3
  } else if (sperm.competition$NumSpermMotile[i] > 15 & sperm.competition$NumSpermMotile[i] <= 20){
    binned.motile[i] <- 4
  } else if (sperm.competition$NumSpermMotile[i] > 20 & sperm.competition$NumSpermMotile[i] <= 25){
    binned.motile[i] <- 5
  } else if (sperm.competition$NumSpermMotile[i] > 25 & sperm.competition$NumSpermMotile[i] <= 30){
    binned.motile[i] <- 6
  } else if (sperm.competition$NumSpermMotile[i] > 30 & sperm.competition$NumSpermMotile[i] <= 35){
    binned.motile[i] <- 7
  } else if (sperm.competition$NumSpermMotile[i] > 35 & sperm.competition$NumSpermMotile[i] <= 40){
    binned.motile[i] <- 8
  } else if (sperm.competition$NumSpermMotile[i] > 40 & sperm.competition$NumSpermMotile[i] <= 45){
    binned.motile[i] <- 9
  } else if (sperm.competition$NumSpermMotile[i] > 45 & sperm.competition$NumSpermMotile[i] <= 50){
    binned.motile[i] <- 10
  } else if (sperm.competition$NumSpermMotile[i] > 50 & sperm.competition$NumSpermMotile[i] <= 55){
    binned.motile[i] <- 11
  } else if (sperm.competition$NumSpermMotile[i] > 55 & sperm.competition$NumSpermMotile[i] <= 60){
    binned.motile[i] <- 12
  } else if (sperm.competition$NumSpermMotile[i] > 60 & sperm.competition$NumSpermMotile[i] <= 65){
    binned.motile[i] <- 13
  } else if (sperm.competition$NumSpermMotile[i] > 65 & sperm.competition$NumSpermMotile[i] <= 70){
    binned.motile[i] <- 14
  } else if (sperm.competition$NumSpermMotile[i] > 70 & sperm.competition$NumSpermMotile[i] <= 75){
    binned.motile[i] <- 15
  } else if (sperm.competition$NumSpermMotile[i] > 75 & sperm.competition$NumSpermMotile[i] <= 80){
    binned.motile[i] <- 16
  } else if (sperm.competition$NumSpermMotile[i] > 80 & sperm.competition$NumSpermMotile[i] <= 85){
    binned.motile[i] <- 17
  } else if (sperm.competition$NumSpermMotile[i] > 85 & sperm.competition$NumSpermMotile[i] <= 90){
    binned.motile[i] <- 18
  } else if (sperm.competition$NumSpermMotile[i] > 90 & sperm.competition$NumSpermMotile[i] <= 95){
    binned.motile[i] <- 19
  } else if (sperm.competition$NumSpermMotile[i] > 95 & sperm.competition$NumSpermMotile[i] <= 100){
    binned.motile[i] <- 20
  } else if (sperm.competition$NumSpermMotile[i] > 100){
    binned.motile[i] <- 21
  }
}


sperm.competition$binned.motile <- binned.motile
ggplot(data = sperm.competition, aes(x = binned.motile, y = AvgVAP)) + 
  geom_point()

sd.bins <- aggregate(cbind(PercentMotility, AvgVAP, AvgLIN) ~ binned.motile, data = sperm.competition, FUN = sd)
cv.bins <- aggregate(cbind(PercentMotility, AvgVAP, AvgLIN) ~ binned.motile, data = sperm.competition, FUN = EnvStats::cv)
num.bins <- aggregate(cbind(PercentMotility, AvgVAP, AvgLIN) ~ binned.motile, data = sperm.competition, FUN = length)

ggplot(data = sd.bins, aes(x = binned.motile, y = PercentMotility)) + 
  geom_point()
ggplot(data = cv.bins, aes(x = binned.motile, y = PercentMotility)) + 
  geom_point()
ggplot(data = sd.bins, aes(x = binned.motile, y = AvgLIN)) + 
  geom_point()
ggplot(data = cv.bins, aes(x = binned.motile, y = AvgLIN)) + 
  geom_point()
ggplot(data = sd.bins, aes(x = binned.motile, y = AvgVAP)) + 
  geom_point()
ggplot(data = cv.bins, aes(x = binned.motile, y = AvgVAP)) + 
  geom_point()

for.analysis.sc <- sperm.competition[which(sperm.competition$NumSpermMotile >= 5),]

################################################################################################################################
#
#
#     And now some data analysis
#
#
#################################################################################################################################

# MOT
mot.global.lmer <- lmer(asin(sqrt(PercentMotility)) ~ Sulfidic * Food * Sulfide + Timepoint + MaleSL + treat.length + (1 | Drainage) + (1 | IndivID), na.action = "na.fail", data = for.analysis.sc)
mot.dredge <- dredge(global.model = mot.global.lmer, rank = "AICc", m.lim = c(1,4))
print(mot.dredge, abbrev.names = TRUE)
#best.mot.lmer <- lmer(asin(sqrt(PercentMotility)) ~ Timepoint + (1 | IndivID), data = for.analysis.sc)
#summary(best.mot.lmer)
#Anova(best.mot.lmer, type = "III")
#plot(Effect("Timepoint", best.mot.lmer))
#eta_squared(best.mot.lmer)
#mot.avg <- model.avg(mot.dredge)#, subset = delta < 7)
mot.avg <- model.avg(mot.dredge, subset = cumsum(weight) <= .95)
summary(mot.avg)
confint(mot.avg)
coefTable(mot.avg, full = F)
confint(mot.avg)

#Weights(mot.dredge)
# full.mot.mod <- lmer(asin(sqrt(PercentMotility)) ~ Timepoint + Sulfidic + Food + MaleSL + (1 | Drainage) + (1 | IndivID), data = for.analysis.sc)
# summary(full.mot.mod)
# Anova(full.mot.mod)
# plot(Effect("MaleSL", full.mot.mod))
# plot(Effect("Timepoint", full.mot.mod))
# plot(Effect("Food", full.mot.mod))
# plot(Effect("Sulfidic", full.mot.mod))

# VAP
vap.global.lmer <- lmer(AvgVAP ~ Sulfidic * Food * Sulfide + Timepoint + MaleSL + treat.length + (1 | Drainage) + (1 | IndivID), na.action = "na.fail", data = for.analysis.sc)
vap.dredge <- dredge(global.model = vap.global.lmer, rank = "AICc", m.lim = c(1,4))
print(vap.dredge, abbrev.names = TRUE)
vap.avg <- model.avg(vap.dredge, subset = cumsum(weight) <= .95)
summary(vap.avg)
confint(vap.avg)
coefTable(vap.avg, full = F)

ggplot(data = for.analysis.sc, aes(x = Sulfide, y = AvgVAP, fill = Food)) + 
  geom_boxplot()
ggplot(data = for.analysis.sc, aes(x = Food, y = AvgVAP, fill = Sulfide)) + 
  geom_boxplot()
ggplot(data = for.analysis.sc, aes(x = Sulfidic, y = AvgVAP, fill = Sulfide)) + 
  geom_boxplot()
ggplot(data = for.analysis.sc, aes(x = Sulfide, y = AvgVAP, fill = Sulfidic)) + 
  geom_boxplot()



# VSL
vsl.global.lmer <- lmer(AvgVSL ~ Sulfidic * Food * Sulfide + Timepoint + MaleSL + treat.length + (1 | Drainage) + (1 | IndivID), na.action = "na.fail", data = for.analysis.sc)
vsl.dredge <- dredge(global.model = vsl.global.lmer, rank = "AICc", m.lim = c(1,4))
print(vsl.dredge, abbrev.names = TRUE)
vsl.avg <- model.avg(vsl.dredge, subset = cumsum(weight) <= .95)
summary(vsl.avg)
confint(vsl.avg)
coefTable(vsl.avg, full = F)

ggplot(data = for.analysis.sc, aes(x = Sulfidic, y = AvgVSL, fill = Food)) + 
  geom_boxplot()
ggplot(data = for.analysis.sc, aes(x = Sulfide, y = AvgVSL, fill = Food)) + 
  geom_boxplot()
ggplot(data = for.analysis.sc, aes(x = MaleSL, y = AvgVSL)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# best.vap.lmer <- lmer(AvgVAP ~ Food * Sulfide + treat.length + (1 | Drainage) + (1 | IndivID), data = sperm.competition)
# full.vap.lmer <- lmer(AvgVAP ~ Food * Sulfide + Sulfide * Sulfidic + Timepoint + treat.length + (1 | Drainage) + (1 | IndivID), data = sperm.competition)
# summary(full.vap.lmer)
# Anova(full.vap.lmer, type = "III")
# eta_squared(full.vap.lmer)
# plot(Effect(c("Sulfide", "Food"), full.vap.lmer))
# plot(Effect(c("Sulfidic", "Sulfide"), full.vap.lmer))
# plot(Effect("Food", full.vap.lmer))
# plot(Effect("Sulfide", full.vap.lmer))
# plot(Effect("Sulfidic", full.vap.lmer))
# plot(Effect("treat.length", full.vap.lmer))
# plot(Effect("Timepoint", full.vap.lmer))


# LIN
lin.global.lmer <- lmer(asin(sqrt(AvgLIN)) ~ Sulfidic * Food * Sulfide + Timepoint + MaleSL + treat.length + (1 | Drainage) + (1 | IndivID), na.action = "na.fail", data = for.analysis.sc)
lin.dredge <- dredge(global.model = lin.global.lmer, rank = "AICc", m.lim = c(1,4))
print(lin.dredge, abbrev.names = TRUE)
lin.avg <- model.avg(lin.dredge, subset = cumsum(weight) <= .95)
summary(lin.avg)
confint(lin.avg)
coefTable(lin.avg, full = F)

best.lin.lmer <- lmer(asin(sqrt(AvgLIN)) ~ Food + (1 | Drainage) + (1 | IndivID), data = sperm.competition)
summary(best.lin.lmer)
Anova(best.lin.lmer, type = "III")
plot(Effect("Food", best.lin.lmer))
eta_squared(best.lin.lmer)



ggplot(data = for.analysis.sc, aes(x = AvgVAP, y = asin(sqrt(PercentMotility)))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = for.analysis.sc, aes(x = AvgVAP, y = asin(sqrt(AvgLIN)))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = for.analysis.sc, aes(x = asin(sqrt(PercentMotility)), y = asin(sqrt(AvgLIN)))) + 
  geom_point() + 
  geom_smooth(method = "lm")


vsl.plot <- ggplot() +
  geom_boxplot(data = for.analysis.sc, aes(x = Sulfidic, y = AvgVSL, color = Food), outlier.alpha = 0) +
  geom_point(data = for.analysis.sc, aes(x = Sulfidic, y = AvgVSL, color = Food), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("Sperm Velocity (VSL)") + 
  theme_classic() + 
  theme(aspect.ratio = 1) + 
  theme(legend.position = "none")
vsl.plot

mot.plot <- ggplot() +
  geom_boxplot(data = for.analysis.sc, aes(x = Sulfidic, y = asin(sqrt(PercentMotility)), color = Food), outlier.alpha = 0) +
  geom_point(data = for.analysis.sc, aes(x = Sulfidic, y = asin(sqrt(PercentMotility)), color = Food), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("Ejaculate Motility (MOT)") + 
  theme_classic() + 
  theme(aspect.ratio = 1) + 
  theme(legend.position = "none")
mot.plot

pdf("~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/04_Manuscript/chapter3_supplemental_figures_first_draft/FigureS2_sperm_comp.pdf", width = 9, height = 3)
plot_grid(mot.plot, vsl.plot, labels = c("A", "B"), ncol = 2, nrow = 1)
dev.off()


mot.sl.plot <- ggplot(data = for.analysis.sc, aes(x = MaleSL, y = asin(sqrt(PercentMotility)))) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_classic() + 
  ylab("MOT (in Sulfide Gradient)") + 
  theme(aspect.ratio = 1)
mot.sl.plot

vsl.food.sulf.plot <- ggplot(data = for.analysis.sc, aes(x = Food, y = AvgVSL, fill = Sulfide)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab("VSL (in Sulfide Gradient)") + 
  #scale_fill_manual(values = c("#FFF59D", "#FFF176", "#FFE838"))
  scale_fill_manual(values = c("#FFF9C4", "#FFF176", "#FFE838")) + 
  theme(aspect.ratio = 1) + 
  theme(legend.position = "none")
vsl.food.sulf.plot
pdf("~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/04_Manuscript/Figure3_first.draft.pdf", height = 3, width = 9)
plot_grid(mot.sl.plot, vsl.food.sulf.plot, labels=c("A", "B"), ncol = 2, nrow = 1)
dev.off()

################################################################################################################################
#
#
#     Multivariate analysis
#
#
#################################################################################################################################

dv <- data.frame(trans.mot = asin(sqrt(for.analysis.sc$PercentMotility)), trans.lin = asin(sqrt(for.analysis.sc$AvgLIN)), vap = for.analysis.sc$AvgVAP)

ggpairs(dv)

sc.pca <- prcomp(dv, center = T, scale. = T)
sc.pca
summary(sc.pca)


pc.data <- as.data.frame(cbind(for.analysis.sc, sc.pca$x))
evec <- sc.pca$rotation
eval <- sc.pca$sdev^2
loadings <- data.frame(matrix(nrow = nrow(evec), ncol = ncol(evec)))
for (i in 1:length(colnames(evec))){
  loadings[,i] <- evec[,i] * sqrt(eval[i])
}

colnames(loadings) <- colnames(evec)
rownames(loadings) <- rownames(evec)

ggplot(data = pc.data, aes(x = PC1, y = PC2, color = Sulfidic)) + 
  geom_point(aes(shape = Food)) + 
  scale_color_manual(values = c("steelblue", "darkgoldenrod1"), name = NULL) + 
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(1/2, "picas")), 
               color = "black", 
               size = .5) + 
  annotate("text", 
           x = loadings$PC1, 
           y = loadings$PC2, 
           label = rownames(loadings), 
           size = 2) + 
  theme_classic()

ggplot(data = pc.data, aes(x = PC2, y = PC3, color = Sulfidic)) + 
  geom_point(aes(shape = Food)) + 
  scale_color_manual(values = c("steelblue", "darkgoldenrod1"), name = NULL) + 
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC2, yend = PC3), 
               arrow = arrow(length = unit(1/2, "picas")), 
               color = "black", 
               size = .5) + 
  annotate("text", 
           x = loadings$PC2, 
           y = loadings$PC3, 
           label = rownames(loadings), 
           size = 2) + 
  theme_classic()

### Analyze PC1
pc1.glob <- lmer(PC1 ~ Sulfidic * Food * Sulfide + Timepoint + MaleSL + treat.length + (1 | Drainage) + (1 | IndivID), na.action = "na.fail", data = pc.data)
pc1.dredge <- dredge(global.model = pc1.glob, rank = "AICc", m.lim = c(1,4))
print(pc1.dredge, abbrev.names = TRUE)
pc1.avg <- model.avg(pc1.dredge, subset = cumsum(weight) <= .95)
summary(pc1.avg)
confint(pc1.avg)

ggplot(data = pc.data, aes(x = treat.length, y = PC1)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(data = pc.data, aes(x = Timepoint, y = PC1)) + 
  geom_boxplot()



### Analyze PC2
pc2.glob <- lmer(PC2 ~ Sulfidic * Food * Sulfide + Timepoint + MaleSL + treat.length + (1 | Drainage) + (1 | IndivID), na.action = "na.fail", data = pc.data)
pc2.dredge <- dredge(global.model = pc2.glob, rank = "AICc", m.lim = c(1,4))
print(pc2.dredge, abbrev.names = TRUE)
pc2.avg <- model.avg(pc2.dredge, subset = cumsum(weight) <= .95)
summary(pc2.avg)
confint(pc2.avg)

ggplot(data = pc.data, aes(x = MaleSL, y = PC2)) + 
  geom_point() + 
  geom_smooth(method = "lm")



### Analyze PC3
pc3.glob <- lmer(PC3 ~ Sulfidic * Food * Sulfide + Timepoint + MaleSL + treat.length + (1 | Drainage) + (1 | IndivID), na.action = "na.fail", data = pc.data)
pc3.dredge <- dredge(global.model = pc3.glob, rank = "AICc", m.lim = c(1,4))
print(pc3.dredge, abbrev.names = TRUE)
pc3.avg <- model.avg(pc3.dredge, subset = cumsum(weight) <= .95)
summary(pc3.avg)
confint(pc3.avg)

ggplot(data = pc.data, aes(x = Sulfide, y = PC3)) + 
  geom_boxplot() 











