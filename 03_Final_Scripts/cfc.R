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
library(ggpubr)


################################################################################################################################
#
#
#     Input data and filtering
#
#
#################################################################################################################################

raw.cfc <- read_excel(path = "~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/02_Raw_Data/cryptic_female_choice_metadata.xlsx")
#discard.ids <- raw.sperm.comp$VideoID[which(is.na(raw.sperm.comp$AvgVCL))]
discard.ids <- raw.cfc$TrialID[which(raw.cfc$ZeroOF.AvgVCL == "NaN" | raw.cfc$AvgVCL == "NaN")]
zero.of.raw <- raw.cfc[which(!(raw.cfc$TrialID %in% discard.ids)), c(1:26, 38)]
zero.of <- zero.of.raw[which(!(duplicated(zero.of.raw$ZeroOF.PercentMotility))),]
zero.of$of.treatment <- NA
for (i in 1:nrow(zero.of)){
  zero.of$of.treatment[i] <- "A.ZeroOF"
}

of.treatment <- raw.cfc[which(!(raw.cfc$TrialID %in% discard.ids)),c(1:15, 27:38)]
of.treatment$of.treatment <- NA
for (j in 1:nrow(of.treatment)){
  if (interaction(of.treatment$MaleSulfidic[j], of.treatment$FemaleSulfidic[j]) == "NS.S") {
    of.treatment$of.treatment[j] <- "C.NSMale.SFemale.OF"
  } else if (interaction(of.treatment$MaleSulfidic[j], of.treatment$FemaleSulfidic[j]) == "NS.NS") {
    of.treatment$of.treatment[j] <- "B.NSMale.NSFemale.OF"
  } else if (interaction(of.treatment$MaleSulfidic[j], of.treatment$FemaleSulfidic[j]) == "S.NS") {
    of.treatment$of.treatment[j] <- "E.SMale.NSFemale.OF"
  } else if (interaction(of.treatment$MaleSulfidic[j], of.treatment$FemaleSulfidic[j]) == "S.S") {
      of.treatment$of.treatment[j] <- "D.SMale.SFemale.OF"
  }
}

names(zero.of)[16:26] <- names(of.treatment[16:26])

cfc <- as.data.frame(rbind(zero.of, of.treatment))
names(cfc)[16] <- "VideoID"

#cfc <- raw.cfc[which(!(raw.cfc$TrialID %in% discard.ids)),] #%>% dplyr::select(VideoID, Timepoint, Tank, IndivID, DateFoodStart, TrialDateTime, Drainage, Sulfidic, Pop, Food, Treatment, Sulfide, MaleSL, PercentMotility, AvgVCL, AvgVAP, AvgVSL, AvgLIN, NumSpermTracked, NumSpermMotile)

numerics <- c("MaleSL", "FemaleSL", "PercentMotility", "AvgVCL", "AvgVAP", "AvgVSL", "AvgLIN", "NumSpermTracked", "NumSpermMotile")
cfc[,which(colnames(cfc) %in% numerics)] <- lapply(cfc[,which(colnames(cfc) %in% numerics)],as.numeric)
factors <- c("TrialID", "MaleID", "FemaleID", "Block", "MaleDrainage", "FemaleDrainage", "MaleSulfidic", "FemaleSulfidic", "MalePop", "FemalePop")
cfc[,which(colnames(cfc) %in% factors)] <- lapply(cfc[,which(colnames(cfc) %in% factors)],as.factor)
#cfc$d.t.food.start <- strptime(sperm.competition$DateFoodStart, format = "%m/%d/%y %H:%M")
cfc$d.t.trial <- strptime(cfc$TrialDate.Time, format = "%m/%d/%y %H:%M")
#sperm.competition$treat.length <- sperm.competition$d.t.trial - sperm.competition$d.t.food.start
#sperm.competition$treat.length <- as.numeric(sperm.competition$treat.length)

cfc.for.analysis <- cfc[which(cfc$NumSpermMotile >= 5),]

ggplot(data = cfc.for.analysis, aes(x = Block, y = AvgVSL)) + 
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
ggplot(data = cfc.for.analysis, aes(x = AvgVAP)) + 
  geom_histogram(binwidth = 3)
ggplot(data = cfc.for.analysis, aes(x = AvgVAP, color = FemaleSulfidic)) + 
  geom_histogram() 
ggplot(data = cfc.for.analysis, aes(x = interaction(MaleSulfidic, FemaleSulfidic), y = AvgVAP)) + 
  geom_boxplot() 
ggplot(data = cfc.for.analysis, aes(x = NumSpermMotile, y = AvgVAP)) + 
  geom_point() + 
  geom_smooth(method = "lm")


# VSL
ggplot(data = cfc.for.analysis, aes(x = AvgVSL)) + 
  geom_histogram(binwidth = 3)
ggplot(data = cfc.for.analysis, aes(x = AvgVSL, color = FemaleSulfidic)) + 
  geom_histogram() 
ggplot(data = cfc.for.analysis, aes(x = interaction(MaleSulfidic, FemaleSulfidic), y = AvgVSL)) + 
  geom_boxplot()
ggplot(data = cfc.for.analysis, aes(x = MaleSulfidic, y = AvgVSL, fill = FemaleSulfidic)) + 
  geom_boxplot()
ggplot(data = cfc.for.analysis, aes(x = NumSpermMotile, y = AvgVSL)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc.for.analysis[which(cfc.for.analysis$of.treatment == "A.ZeroOF"),], aes(x = MaleSulfidic, y = AvgVSL)) + 
  geom_boxplot() 
# VSL in zero OF
ggplot(data = cfc.for.analysis[which(cfc.for.analysis$of.treatment == "A.ZeroOF"),], aes(x = MaleSulfidic, y = AvgVSL, color = Timepoint)) + 
  geom_boxplot() 
ggplot(data = cfc.for.analysis[which(cfc.for.analysis$Timepoint == "0.5-2"),], aes(x = of.treatment, y = AvgVSL)) + 
  geom_boxplot() 
ggplot(data = cfc.for.analysis[which(cfc.for.analysis$Timepoint == "10.5-12"),], aes(x = of.treatment, y = AvgVSL)) + 
  geom_boxplot() 
ggplot(data = cfc.for.analysis, aes(x = of.treatment, y = AvgVSL)) + 
  geom_boxplot() 

ggplot(data = cfc.for.analysis, aes(x = AvgVSL, y = AvgVAP)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  stat_regline_equation(label.x=30, label.y=110) + 
  stat_cor(aes(label=..rr.label..), label.x=30, label.y=100)
ggplot(data = cfc.for.analysis, aes(x = AvgVCL, y = AvgVAP)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  stat_regline_equation(label.x=60, label.y=110) + 
  stat_cor(aes(label=..rr.label..), label.x=60, label.y=100)
ggplot(data = cfc.for.analysis, aes(x = AvgVCL, y = AvgVSL)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  stat_regline_equation(label.x=60, label.y=110) + 
  stat_cor(aes(label=..rr.label..), label.x=60, label.y=100)


# MOT
ggplot(data = cfc.for.analysis, aes(x = asin(sqrt(PercentMotility)))) + 
  geom_histogram(binwidth = .05)
# ggplot(data = cfc.for.analysis, aes(x = asin(sqrt(PercentMotility)), color = Drainage)) + 
#   geom_density() 
ggplot(data = cfc.for.analysis, aes(x = of.treatment, y = asin(sqrt(PercentMotility)))) + 
  geom_boxplot() 
ggplot(data = cfc.for.analysis[which(cfc.for.analysis$of.treatment == "A.ZeroOF"),], aes(x = MaleSulfidic, y = asin(sqrt(PercentMotility)), color = Timepoint)) + 
  geom_boxplot() 
ggplot(data = cfc.for.analysis[which(cfc.for.analysis$of.treatment != "A.ZeroOF"),], aes(x = of.treatment, y = asin(sqrt(PercentMotility)))) + 
  geom_boxplot() 


# LIN
ggplot(data = cfc.for.analysis, aes(x = asin(sqrt(AvgLIN)))) + 
  geom_histogram(binwidth = .01)

ggplot(data = cfc.for.analysis, aes(x = of.treatment, y = asin(sqrt(AvgLIN)))) + 
  geom_boxplot() 

# VSL vs MOT
ggplot(data = cfc.for.analysis, aes(x = AvgVSL, y = asin(sqrt(PercentMotility)))) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  stat_regline_equation(label.x=90, label.y=1) + 
  stat_cor(aes(label=..rr.label..), label.x=90, label.y=.9)

# VSL vs LIN
ggplot(data = cfc.for.analysis, aes(x = AvgVSL, y = asin(sqrt(AvgLIN)))) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  stat_regline_equation(label.x=90, label.y=1) + 
  stat_cor(aes(label=..rr.label..), label.x=90, label.y=.9)

# MOT vs LIN
ggplot(data = cfc.for.analysis, aes(x = asin(sqrt(PercentMotility)), y = asin(sqrt(AvgLIN)))) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  stat_regline_equation(label.x=.3, label.y=1) + 
  stat_cor(aes(label=..rr.label..), label.x=.3, label.y=.9)



#### Now some plotting of covariates
# Effects of treatment length
# ggplot(data = cfc, aes(x = treat.length)) + 
#   geom_histogram(binwidth = 1)
# ggplot(data = cfc, aes(x = treat.length, y = AvgVAP, color = Drainage)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")
# ggplot(data = cfc, aes(x = treat.length, y = asin(sqrt(PercentMotility)), color = Drainage)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")
# ggplot(data = cfc, aes(x = treat.length, y = asin(sqrt(AvgLIN)), color = Drainage)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")

# Effects of SL
ggplot(data = cfc, aes(x = MaleSL)) + 
  geom_histogram(binwidth = 1)
ggplot(data = cfc, aes(x = FemaleSL)) + 
  geom_histogram(binwidth = 1)
ggplot(data = cfc, aes(x = MaleSL, y = AvgVSL)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc, aes(x = FemaleSL, y = AvgVSL)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc, aes(x = MaleSL, y = asin(sqrt(PercentMotility)))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc, aes(x = FemaleSL, y = asin(sqrt(PercentMotility)))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc, aes(x = MaleSL, y = asin(sqrt(AvgLIN)))) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc, aes(x = FemaleSL, y = asin(sqrt(AvgLIN)))) + 
  geom_point() + 
  geom_smooth(method = "lm")

# Effects of sperm numbers
ggplot(data = cfc, aes(x = NumSpermMotile, color = Drainage)) + 
  geom_histogram(binwidth = 1)
ggplot(data = cfc, aes(x = NumSpermMotile, y = AvgVSL)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc[which(cfc$NumSpermMotile >= 10),], aes(x = NumSpermMotile, y = AvgVAP, color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc, aes(x = NumSpermMotile, y = asin(sqrt(PercentMotility)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = cfc, aes(x = NumSpermMotile, y = asin(sqrt(AvgLIN)), color = Drainage)) + 
  geom_point() + 
  geom_smooth(method = "lm")

final.cfc <- cfc.for.analysis[which(cfc.for.analysis$of.treatment != "A.ZeroOF"),]
final.cfc$foreign.native <- NA
for (k in 1:nrow(final.cfc)){
  if (final.cfc$of.treatment[k] == "B.NSMale.NSFemale.OF") {
    final.cfc$foreign.native[k] <- "native"
  } else if (final.cfc$of.treatment[k] == "C.NSMale.SFemale.OF") {
    final.cfc$foreign.native[k] <- "foreign"
  } else if (final.cfc$of.treatment[k] == "D.SMale.SFemale.OF") {
    final.cfc$foreign.native[k] <- "native"
  } else if (final.cfc$of.treatment[k] == "E.SMale.NSFemale.OF") {
    final.cfc$foreign.native[k] <- "foreign"
  }
}

ggplot(data = final.cfc[which(final.cfc$Timepoint == "0.5-2"),], aes(x = foreign.native, y = AvgVSL)) + 
  geom_boxplot()
ggplot(data = final.cfc[which(final.cfc$Timepoint == "10.5-12"),], aes(x = foreign.native, y = AvgVSL)) + 
  geom_boxplot()
ggplot(data = final.cfc, aes(x = foreign.native, y = AvgVSL, color = MaleSulfidic)) + 
  geom_boxplot()
ggplot(data = final.cfc, aes(x = foreign.native, y = AvgVSL, color = Timepoint)) + 
  geom_boxplot()
ggplot(data = final.cfc, aes(x = MaleSulfidic, y = AvgVSL, fill = foreign.native)) + 
  geom_boxplot()
ggplot(data = final.cfc, aes(x = MaleSulfidic, y = asin(sqrt(PercentMotility)), fill = foreign.native)) + 
  geom_boxplot()


# pdf("~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/04_Manuscript/Figure4_first_draft.pdf", height = 4, width = 4)
# ggplot(data = final.cfc, aes(x = MaleSL, y = AvgVSL)) + 
#   geom_point() + 
#   theme_classic() + 
#   ylab("VSL in OF") + 
#   geom_smooth(method = "lm")
# dev.off()




#global.cfc.vsl.mod <- lmer(formula = AvgVSL ~ MaleSL + FemaleSL + MaleSulfidic * FemaleSulfidic * Timepoint + (1 | Block) + (1 | TrialID), data = cfc.for.analysis, na.action = "na.fail")
global.cfc.vsl.mod <- lmer(formula = AvgVSL ~ foreign.native * MaleSulfidic * Timepoint + MaleSL + FemaleSL + (1 | MaleID) + (1 | FemaleID), data = final.cfc, na.action = "na.fail")
#global.cfc.vsl.mod <- lmer(formula = AvgVSL ~ of.treatment + Timepoint + (1 | MaleID) + (1 | FemaleID), data = cfc, na.action = "na.fail")
cfc.vsl.dredge <- dredge(global.cfc.vsl.mod, m.lim = c(0,4))
print(cfc.vsl.dredge, abbrev.names = TRUE)
cfc.vsl.avg <- model.avg(cfc.vsl.dredge, subset = cumsum(weight) <= .95)
summary(cfc.vsl.avg)
confint(cfc.vsl.avg)

emmeans::emmeans(global.cfc.vsl.mod, pairwise ~ of.treatment)

# top.cfc.vsl.mod <- lmer(AvgVSL ~ foreign.native * MaleSulfidic + MaleSL + (1 | MaleID) + (1 | FemaleID), data = final.cfc)
# Anova(top.cfc.vsl.mod, type = "III")
global.cfc.vap.mod <- lmer(formula = AvgVAP ~ foreign.native * MaleSulfidic * Timepoint + MaleSL + FemaleSL + (1 | MaleID) + (1 | FemaleID), data = final.cfc, na.action = "na.fail")
cfc.vap.dredge <- dredge(global.cfc.vap.mod, m.lim = c(0,4))
print(cfc.vap.dredge, abbrev.names = TRUE)
cfc.vap.avg <- model.avg(cfc.vap.dredge, subset = cumsum(weight) <= .95)
summary(cfc.vap.avg)
confint(cfc.vap.avg)

global.cfc.mot.mod <- lmer(formula = asin(sqrt(PercentMotility)) ~ foreign.native * MaleSulfidic * Timepoint + MaleSL + FemaleSL + (1 | MaleID) + (1 | FemaleID), data = final.cfc, na.action = "na.fail")
cfc.mot.dredge <- dredge(global.cfc.mot.mod, m.lim = c(0,4))
print(cfc.mot.dredge, abbrev.names = TRUE)
cfc.mot.avg <- model.avg(cfc.mot.dredge, subset = cumsum(weight) <= .95)
summary(cfc.mot.avg)
confint(cfc.mot.avg)


x <- lmer(AvgVSL ~ MaleSL + Timepoint + of.treatment + (1 + of.treatment | Block), data = final.cfc)
x2 <- lmer(AvgVSL ~ MaleSL + Timepoint + of.treatment + (1 | Block / TrialID), data = final.cfc)
x.null <- lmer(AvgVSL ~ MaleSL + Timepoint + (1 + of.treatment | Block), data = final.cfc)
anova(x, x2)
summary(x2)

coef(x2)

y <- lmerTest::lmer(AvgVSL ~ MaleSL + Timepoint + of.treatment + (1 | Block / TrialID), data = final.cfc)
y <- lmerTest::lmer(AvgVSL ~ MaleSL + Timepoint + of.treatment + (1 | MaleID) + (1 | FemaleID) + (1 | TrialID), data = final.cfc)

anova(y, type = "III")
coef(y)
summary(y)
ggplot(data = final.cfc, aes(x = interaction(MaleID, FemaleID), y = AvgVSL, color = Block)) + 
  geom_point()

# diff.df <- data.frame(NA)
# for (x in 1:length(unique(final.cfc$TrialID))){
#   curr.data <- as.data.frame(final.cfc[which(final.cfc$TrialID == as.character(unique(final.cfc$TrialID)[x])),])
# }



# Relative data
cfc.rel <- raw.cfc[which(!(raw.cfc$TrialID %in% discard.ids)), c(1:15, 39:46, 26, 37)]

rel.numerics <- c("MaleSL", "FemaleSL", "relMOT", "relVCL", "relVAP", "relVSL", "relLIN", "relWOB", "relPROG", "relBCF", "ZeroOF.NumSpermMotile", "NumSpermMotile")
cfc.rel[,which(colnames(cfc.rel) %in% rel.numerics)] <- lapply(cfc.rel[,which(colnames(cfc.rel) %in% rel.numerics)],as.numeric)
rel.factors <- c("TrialID", "MaleID", "FemaleID", "Block", "MaleDrainage", "FemaleDrainage", "MaleSulfidic", "FemaleSulfidic", "MalePop", "FemalePop")
cfc.rel[,which(colnames(cfc.rel) %in% rel.factors)] <- lapply(cfc.rel[,which(colnames(cfc.rel) %in% rel.factors)],as.factor)

cfc.rel$foreign.native <- NA
for (y in 1:nrow(cfc.rel)){
  if (interaction(cfc.rel$MalePop[y], cfc.rel$FemalePop[y]) == "BON.BON") {
    cfc.rel$foreign.native[y] <- "native"
  } else if (interaction(cfc.rel$MalePop[y], cfc.rel$FemalePop[y]) == "PSO.PSO") {
    cfc.rel$foreign.native[y] <- "native"
  } else if (interaction(cfc.rel$MalePop[y], cfc.rel$FemalePop[y]) == "PSO.BON") {
    cfc.rel$foreign.native[y] <- "foreign"
  } else if (interaction(cfc.rel$MalePop[y], cfc.rel$FemalePop[y]) == "BON.PSO") {
    cfc.rel$foreign.native[y] <- "foreign"
  }
}


ggplot(data = cfc.rel, aes(x = interaction(MaleSulfidic, FemaleSulfidic), y = relVSL)) + 
  geom_boxplot()

ggplot(data = cfc.rel, aes(x = interaction(MaleSulfidic, FemaleSulfidic), y = relMOT)) + 
  geom_boxplot()

ggplot(data = cfc.rel, aes(x = relVSL)) + 
  geom_histogram()
ggplot(data = cfc.rel, aes(x = sqrt(relMOT))) + 
  geom_histogram()
ggplot(data = cfc.rel, aes(x = MaleSulfidic, y = relVSL, color = foreign.native)) + 
  geom_boxplot()
ggplot(data = cfc.rel, aes(x = MaleSulfidic, y = relMOT, color = foreign.native)) + 
  geom_boxplot()


mot.plot <- ggplot() +
  geom_boxplot(data = final.cfc, aes(x = MaleSulfidic, y = asin(sqrt(PercentMotility)), color = foreign.native), outlier.alpha = 0) +
  geom_point(data = final.cfc, aes(x = MaleSulfidic, y = asin(sqrt(PercentMotility)), color = foreign.native), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("Ejaculate Motility (MOT)") + 
  theme_classic() + 
  theme(aspect.ratio = 1) + 
  theme(legend.position = "none")
mot.plot

vsl.plot <- ggplot() +
  geom_boxplot(data = final.cfc, aes(x = MaleSulfidic, y = AvgVSL, color = foreign.native), outlier.alpha = 0) +
  geom_point(data = final.cfc, aes(x = MaleSulfidic, y = AvgVSL, color = foreign.native), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("VSL") + 
  theme_classic() + 
  theme(aspect.ratio = 1) + 
  theme(legend.position = "none")
vsl.plot

pdf("~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/04_Manuscript/chapter3_supplemental_figures_first_draft/FigureS4_sperm_comp.pdf", width = 9, height = 3)
plot_grid(mot.plot, vsl.plot, labels = c("A", "B"), ncol = 2, nrow = 1)
dev.off()

ggplot() +
  geom_boxplot(data = final.cfc, aes(x = MaleSulfidic, y = AvgVSL, color = foreign.native), outlier.alpha = 0) +
  geom_point(data = final.cfc, aes(x = MaleSulfidic, y = AvgVSL, color = foreign.native), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("VSL") + 
  theme_classic() + 
  theme(aspect.ratio = 1)
ggplot() +
  geom_boxplot(data = final.cfc, aes(x = MaleSulfidic, y = asin(sqrt(PercentMotility)), color = foreign.native), outlier.alpha = 0) +
  geom_point(data = final.cfc, aes(x = MaleSulfidic, y = asin(sqrt(PercentMotility)), color = foreign.native), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("Ejaculate Motility (MOT)") + 
  theme_classic() + 
  theme(aspect.ratio = 1)
