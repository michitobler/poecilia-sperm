rm(list = ls())

library(ggplot2)
library(stringr)
library(plyr)
library(dplyr)
library(readxl)
library(lubridate)
library(data.table)
library(lme4)
library(car)


#raw.metadata <- read_xlsx("~/Downloads/viability_test_2-8-22/sperm_viability_metadata.xlsx")
raw.metadata <- read_xlsx("~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/02_Raw_Data/viability/sperm_viability_metadata.xlsx")

# remove problematic samples that have a "0" in the "keep" column
#metadata <- raw.metadata[is.na(raw.metadata$Notes),]
metadata <- raw.metadata[which(raw.metadata$Keep == 1),]
#metadata <- metadata[1:60,] #remove this once more data is input
# metadata$ImageID <- as.character(metadata$ImageID)
# metadata$OutputResultsC0 <- as.character(metadata$OutputResultsC0)
# metadata$OutputResultsC1 <- as.character(metadata$OutputResultsC1)
metadata$DateFoodStart <- as.Date(metadata$DateFoodStart, format = "%m-%d-%y")
metadata$TrialDate <- as.Date(metadata$TrialDate, format = "%m-%d-%y")
metadata$SL.mm <- as.numeric(metadata$SL.mm)
# metadata$TankID <- as.factor(metadata$TankID)
# metadata$Pop <- as.factor(metadata$Pop)
# metadata$Indiv <- as.factor(metadata$Indiv)
# metadata$TrialID <- as.factor(metadata$TrialID)
# metadata$IndivID <- as.factor(metadata$IndivID)
# metadata$Food <- as.factor(metadata$Food)
# metadata$Timpoint <- as.factor(metadata$Timpoint)
# metadata$Replicate <- as.factor(metadata$Replicate)


metadata$extractT <- hm(metadata$ExtractionTime)
metadata$activT0 <- hm(metadata$t0ActivationTime)
metadata$activT4 <- hm(metadata$t4ActivationTime)
metadata$DaysFoodTreat <- metadata$TrialDate - metadata$DateFoodStart

liveNum <- c()
liveArea <- c()
for (i in 1:nrow(metadata)){
  live <- read.csv(paste("~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/02_Raw_Data/viability/raw_viability_data/", metadata$OutputResultsC0[i], sep = ""))
  liveNum <- append(liveNum, nrow(live))
  liveArea <- append(liveArea, sum(live$Area))
}

deadNum <- c()
deadArea <- c()
for (i in 1:nrow(metadata)){
  dead <- read.csv(paste("~/Desktop/Projects/2020_Sperm_Competition_Pmex_Reproductive_Isolation/02_Raw_Data/viability/raw_viability_data/", metadata$OutputResultsC1[i], sep = ""))
  deadNum <- append(deadNum, nrow(dead))
  deadArea <- append(deadArea, sum(dead$Area))
}

metadata$liveNum <- liveNum
metadata$liveArea <- liveArea
metadata$deadNum <- deadNum
metadata$deadArea <- deadArea

# avgPropLiveNum <- c()
# avgPropLiveArea <- c()
# avg.of.props.num <- c()
# avg.of.props.area <- c()

# for (j in seq(2, nrow(metadata), by = 2)){
#   avgPropLiveNum <- append(avgPropLiveNum, 
#                            rep((metadata$liveNum[j] + metadata$liveNum[j - 1]) / (metadata$liveNum[j] + metadata$liveNum[j - 1] + metadata$deadNum[j] + metadata$deadNum[j - 1]), 
#                                2))
#   avgPropLiveArea <- append(avgPropLiveArea, 
#                             rep((metadata$liveArea[j] + metadata$liveArea[j - 1]) / (metadata$liveArea[j] + metadata$liveArea[j - 1] + metadata$deadArea[j] + metadata$deadArea[j - 1]), 
#                                 2))
#   avg.of.props.num <- append(avg.of.props.num, 
#                              rep(mean(c((metadata$liveNum[j] / (metadata$liveNum[j] + metadata$deadNum[j])), (metadata$liveNum[j-1] / (metadata$liveNum[j-1] + metadata$deadNum[j-1])))), 
#                                  2))
#   avg.of.props.area <- append(avg.of.props.area, 
#                               rep(mean(c((metadata$liveArea[j] / (metadata$liveArea[j] + metadata$deadArea[j])), (metadata$liveArea[j-1] / (metadata$liveArea[j-1] + metadata$deadArea[j-1])))),
#                                   2))
# }
metadata$propNum <- metadata$liveNum / (metadata$liveNum + metadata$deadNum)
avg.of.props.num <- aggregate(propNum ~ TrialID, data = metadata, FUN = mean)
avgLiveNum <- aggregate(liveNum ~ TrialID, data = metadata, FUN = mean)
avgDeadNum <- aggregate(deadNum ~ TrialID, data = metadata, FUN = mean)
avg.num.df <- merge(avgLiveNum, avgDeadNum, by = "TrialID")
avg.num.df$avgPropLiveNum <- avg.num.df$liveNum / (avgLiveNum$liveNum + avgDeadNum$deadNum)

#avgPropLiveNum <- avgLiveNum$liveNum / (avgLiveNum$liveNum + avgDeadNum$deadNum)


metadata$propArea <- metadata$liveArea / (metadata$liveArea + metadata$deadArea)
avg.of.props.area <- aggregate(propArea ~ TrialID, data = metadata, FUN = mean)
avgLiveArea <- aggregate(liveArea ~ TrialID, data = metadata, FUN = mean)
avgDeadArea <- aggregate(deadArea ~ TrialID, data = metadata, FUN = mean)
#avgPropLiveArea <- avgLiveArea$liveArea / (avgLiveArea$liveArea + avgDeadArea$deadArea)
avg.area.df <- merge(avgLiveArea, avgDeadArea, by = "TrialID")
avg.area.df$avgPropLiveArea <- avg.area.df$liveArea / (avgLiveArea$liveArea + avgDeadArea$deadArea)

avg.data <- merge(avg.num.df, avg.area.df, by = "TrialID")
avg.df <- join_all(list(avgLiveNum, avgDeadNum, avgLiveArea, avgDeadArea, avg.of.props.num, avg.of.props.area), by = "TrialID")
#avg <- as.data.frame(matrix(nrow = nrow(avg.data), ncol = ncol(metadata)))
#avg <- as.data.frame(matrix(nrow = nrow(avg.data), ncol = 17))
avg <- NULL

#colnames(avg) <- colnames(metadata)

# for (i in 1:nrow(avg.data)) {
#   ids <- which(as.character(metadata$TrialID) == as.character(avg.data$TrialID[i]))
#   row <- metadata[ids[1], as.character(c("ImageID", "DateFoodStart", "TrialDate", "SL.mm", "ExtractionTime", "t0ActivationTime", 
#                      "t4ActivationTime", "TankID", "Pop", "Indiv", "IndivID", "Food", "Timpoint", "extractT", 
#                      "activT0", "activT4", "DaysFoodTreat"))]
#   #avg[i,] <- row
#   avg <- as.data.frame(rbind(avg, row))
# }
for (i in 1:nrow(avg.data)) {
  ids <- which(as.character(metadata$TrialID) == as.character(avg.data$TrialID[i]))
  row <- as.data.frame(metadata[ids[1], c("ImageID", "DateFoodStart", "TrialDate", "SL.mm", "ExtractionTime", "t0ActivationTime", 
                                         "t4ActivationTime", "TankID", "Pop", "Indiv", "IndivID", "Food", "Timpoint", "extractT", 
                                         "activT0", "activT4", "DaysFoodTreat")])
  avg <- rbind(avg, row)
  #avg[i,] <- as.character(row)
  #avg <- as.data.frame(rbind(avg, row))
}


#colnames(avg) <- c("ImageID", "DateFoodStart", "TrialDate", "SL.mm", "ExtractionTime", "t0ActivationTime", 
#                   "t4ActivationTime", "TankID", "Pop", "Indiv", "IndivID", "Food", "Timpoint", "extractT", 
#                   "activT0", "activT4", "DaysFoodTreat")
#cbind(avg.of.props.num, avg.of.props.area, avgPropLiveNum, avgPropLiveArea)
# avgPropLiveNum <- A
# metadata$avgPropLiveNum <- avgPropLiveNum
# metadata$avgPropLiveArea <- avgPropLiveArea
# metadata$avg.of.props.num <- avg.of.props.num
# metadata$avg.of.props.area <- avg.of.props.area

#avg <- metadata[seq(1, nrow(metadata), by = 2),c("ImageID", "DateFoodStart", "TrialDate", "SL.mm", "ExtractionTime", "t0ActivationTime", 
#                                                 "t4ActivationTime", "TankID", "Pop", "Indiv", "IndivID", "Food", "Timpoint", "extractT", 
#                                                 "activT0", "activT4", "DaysFoodTreat", "avgPropLiveNum", "avgPropLiveArea", "avg.of.props.num", "avg.of.props.area")]

avg$date.time0 <- paste(avg$TrialDate, avg$t0ActivationTime)
avg$date.time4 <- paste(avg$TrialDate, avg$t4ActivationTime)
avg$timeBetweenAssays <- as.numeric(difftime(avg$date.time4, avg$date.time0))
names.trial <- as.data.frame(str_split(avg$ImageID, pattern = "_", simplify = T)[,1:6])
names(names.trial) <- c("date", "tank", "pop", "food", "indiv", "timepoint")
#names.all <- do.call(paste, c(names.trial[col.names.trial], sep = "_"))
avg$TrialID <- paste(names.trial$date, names.trial$tank, names.trial$pop, names.trial$food, names.trial$indiv, names.trial$timepoint, sep = "_")


final.dataset <- merge(avg, avg.data, by = "TrialID")




ggplot(data = final.dataset, aes(x = SL.mm, y = avgPropLiveArea)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = final.dataset, aes(x = TrialDate, y = avgPropLiveArea)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = final.dataset, aes(x = DaysFoodTreat, y = avgPropLiveArea)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = final.dataset, aes(x = timeBetweenAssays, y = avgPropLiveArea)) + 
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(data = final.dataset, aes(x = timeBetweenAssays)) + 
  geom_histogram()
# ggplot(data = final.dataset, aes(x = Pop, y = avgPropLiveArea, fill = Food, color = Timpoint)) + 
#   geom_boxplot() + 
#   geom_point(data = final.dataset, aes(x = Pop, y = avgPropLiveArea, color = Food), position = position_jitterdodge(.15), alpha = 0.4, size = 2.5) + 
#   scale_color_manual(values = c("black", "grey")) + 
#   scale_fill_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL)

ggplot() + 
  geom_boxplot(data = final.dataset, aes(x = factor(interaction(Pop, Food, Timpoint), levels = c("BON.H.0hr", "BON.L.0hr", "PSO.H.0hr", "PSO.L.0hr", "BON.H.4hr", "BON.L.4hr", "PSO.H.4hr", "PSO.L.4hr"), ordered = T), y = avgPropLiveArea, color = factor(interaction(Pop, Food, Timpoint), levels = c("BON.H.0hr", "BON.L.0hr", "PSO.H.0hr", "PSO.L.0hr", "BON.H.4hr", "BON.L.4hr", "PSO.H.4hr", "PSO.L.4hr"), ordered = T))) + 
  geom_point(data = final.dataset, aes(x = factor(interaction(Pop, Food, Timpoint), levels = c("BON.H.0hr", "BON.L.0hr", "PSO.H.0hr", "PSO.L.0hr", "BON.H.4hr", "BON.L.4hr", "PSO.H.4hr", "PSO.L.4hr"), ordered = T), y = avgPropLiveArea, color = factor(interaction(Pop, Food, Timpoint), levels = c("BON.H.0hr", "BON.L.0hr", "PSO.H.0hr", "PSO.L.0hr", "BON.H.4hr", "BON.L.4hr", "PSO.H.4hr", "PSO.L.4hr"), ordered = T)), position = position_jitterdodge(.75), alpha = 0.4, size = 2.5) + 
  scale_color_manual(values = c("steelblue", "steelblue", "darkgoldenrod1", "darkgoldenrod1", "steelblue", "steelblue", "darkgoldenrod1", "darkgoldenrod1"), name = NULL) 

ggplot() + 
  geom_boxplot(data = final.dataset, aes(x = Pop, y = avgPropLiveArea, color = Food)) + 
  geom_point(data = final.dataset, aes(x = Pop, y = avgPropLiveArea, color = Food)) + 
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) 


ggplot(data = final.dataset, aes(x = Pop, y = avgPropLiveArea, fill = Pop)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("steelblue", "darkgoldenrod1"), name = NULL)
ggplot(data = final.dataset, aes(x = Food, y = avgPropLiveArea, fill = Food)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL)

ggplot(data = final.dataset[which(final.dataset$Timpoint == "0hr"),], aes(x = Pop, y = avgPropLiveArea, fill = Food)) + 
  geom_boxplot() + 
  scale_color_manual(values = c("black", "grey")) + 
  ylim(c(0, 1))
ggplot(data = final.dataset[which(final.dataset$Timpoint == "4hr"),], aes(x = Pop, y = avgPropLiveArea, fill = Food)) + 
  geom_boxplot() + 
  scale_color_manual(values = c("black", "grey")) + 
  ylim(c(0, 1))

t0.plot <- ggplot() +
  geom_boxplot(data = final.dataset[which(final.dataset$Timpoint == "0hr"),], aes(x = Pop, y = avgPropLiveArea, color = Food), outlier.alpha = 0) +
  geom_point(data = final.dataset[which(final.dataset$Timpoint == "0hr"),], aes(x = Pop, y = avgPropLiveArea, color = Food), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("Sperm Viability 0 hr After Activation") + 
  theme_classic() + 
  ylim(c(0,1))
t0.plot

t4.plot <- ggplot() +
  geom_boxplot(data = final.dataset[which(final.dataset$Timpoint == "4hr"),], aes(x = Pop, y = avgPropLiveArea, color = Food), outlier.alpha = 0) +
  geom_point(data = final.dataset[which(final.dataset$Timpoint == "4hr"),], aes(x = Pop, y = avgPropLiveArea, color = Food), position = position_jitterdodge(0.15), alpha = 0.4, size = 2.5) +
  scale_color_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL) + 
  ylab("Sperm Viability 4 hr After Activation") + 
  theme_classic() + 
  ylim(c(0,1))
t4.plot

pdf("~/Downloads/FigureS3_sperm_comp.pdf", width = 9, height = 3)
plot_grid(t0.plot, t4.plot, labels = c("A", "B"), ncol = 2, nrow = 1)
dev.off()



df <- final.dataset
doubles.df <- df %>% 
  dplyr::group_by(IndivID) %>%
  dplyr::filter(n_distinct(avgPropLiveArea) == 2)

# diff.df <- data.frame(matrix(nrow = nrow(doubles.df) / 2, ncol = ncol(doubles.df)))
# names(diff.df) <- names(doubles.df)

# for (i in unique(doubles.df$IndivID)){
#   curr.data <- doubles.df[which(doubles.df$IndivID == i),]
#   diff.prop.live <- curr.data$avgPropLiveArea[2] - curr.data$avgPropLiveArea[1]
# }
diff.df <- data.frame()
diff.prop.live <- c()
for (i in 1:length(unique(doubles.df$IndivID))){
  curr.data <- as.data.frame(doubles.df[which(doubles.df$IndivID == unique(doubles.df$IndivID)[i]),])
  single.data <- curr.data[1,]
  single.data$t0.area <- single.data$avgPropLiveArea
  second.data <- curr.data[2,]
  single.data$t4.area <- second.data$avgPropLiveArea
  diff.df <- rbind(diff.df, single.data)
  #diff.df[i,] <- curr.data[1,]  
  
  curr.diff.prop.live <- as.numeric(curr.data$avgPropLiveArea[2]) - as.numeric(curr.data$avgPropLiveArea[1])
  diff.prop.live <- append(diff.prop.live, curr.diff.prop.live)
  #diff.prop.live <- as.numeric(diff.prop.live)
  #diff.df$diff.prop.live[i] <- as.numeric(diff.prop.live)
}
diff.df$diff.prop.live <- diff.prop.live

ggplot(data = diff.df, aes(x = diff.prop.live, color = Pop)) + 
  geom_density()

#global.model <- lmer(asin(sqrt(avgPropLiveArea)) ~ Pop*Food + SL.mm  + DaysFoodTreat + Timpoint + (1 | IndivID), data = final.dataset, na.action = "na.fail")
t0.df <- final.dataset[which(final.dataset$Timpoint == "0hr"),]
global.model <- lmer(asin(sqrt(avgPropLiveArea)) ~ Pop*Food + Timpoint + DaysFoodTreat + SL.mm + (1 | IndivID), data = final.dataset, na.action = "na.fail")

dredge.mod <- dredge(global.model, rank = "AICc", m.lim = c(0,4))
print(dredge.mod, abbrev.names = TRUE)
v.avg <- model.avg(dredge.mod, subset = cumsum(weight) <= .95)
summary(v.avg)

best.v.mod <- lmer(asin(sqrt(avgPropLiveArea)) ~ Timpoint + (1 | IndivID), data = final.dataset)
Anova(best.v.mod)

best.viability.model <- lmerTest::lmer(asin(sqrt(avgPropLiveArea)) ~ Timpoint + (1 | IndivID), data = final.dataset)
anova(best.viability.model)
confint(best.v.mod)
confint(best.viability.model)
eta_squared(best.v.mod)
coefTable(best.v.mod)



global.model.diff <- lm(asin(sqrt(diff.prop.live)) ~ Pop*Food + SL.mm + DaysFoodTreat, data = diff.df, na.action = "na.fail")
dredge.diff.mod <- dredge(global.model.diff, rank = "AICc", m.lim = c(1,4))
print(dredge.diff.mod)
v.diff.avg <- model.avg(dredge.diff.mod, subset = cumsum(weight) <= .95)
summary(v.diff.avg)


ggplot(data = diff.df, aes(x = DaysFoodTreat, y = diff.prop.live)) + 
  geom_point() + 
  geom_smooth(method = "lm")



Anova(lmer(asin(sqrt(avgPropLiveArea)) ~ Timpoint + (1 | IndivID), data = final.dataset, na.action = "na.fail"))#write.csv(print(dredge.out, abbrev.names = TRUE), "~/Downloads/explore.mod.select.csv")
confint(global.model)

#interaction is non-significant, so using additive model instead
best.explore.mod <- lmer(PC1 ~ pop + food.treat + (1 | tank), data = new.explore, na.action = "na.fail")

#best.explore.mod <- lmer(PC1 ~ pop * food.treat + (1 | tank), data = new.explore, na.action = "na.fail")
Anova(best.explore.mod, type = "III")
confint(best.explore.mod)


summary(glmer(avgPropLiveArea ~ Pop*Food + DaysFoodTreat + Timpoint + (1 | IndivID), data = avg, family = "binomial", nAGQ = 0))
summary(glmer(avgPropLiveArea ~ Pop*Food + Timpoint + (1 | IndivID), data = avg, family = "binomial", nAGQ = 0))



# avg$diffActivTime <- avg$activT4 - avg$activT0
# ggplot(data = avg, aes(x = SL.mm, y = avgPropLiveArea)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")
# ggplot(data = avg, aes(x = TrialDate, y = avgPropLiveArea)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")
# ggplot(data = avg, aes(x = DaysFoodTreat, y = avgPropLiveArea)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")
# ggplot(data = avg, aes(x = timeBetweenAssays, y = avgPropLiveArea)) + 
#   geom_point() + 
#   geom_smooth(method = "lm")
# ggplot(data = avg, aes(x = timeBetweenAssays)) + 
#   geom_histogram()
# ggplot(data = avg, aes(x = Pop, y = avgPropLiveArea, fill = Food, color = Timpoint)) + 
#   geom_boxplot() + 
#   scale_color_manual(values = c("black", "grey")) + 
#   scale_fill_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL)
# ggplot(data = avg, aes(x = Pop, y = avgPropLiveArea, fill = Food, color = Timpoint)) + 
#   geom_boxplot() + 
#   scale_color_manual(values = c("black", "grey")) + 
#   scale_fill_manual(values = c("#5DA5DA", "#FAA43A"), name = NULL)
# ggplot(data = avg, aes(x = interaction(Pop, Food), y = avgPropLiveArea, color = Timpoint)) + 
#   geom_boxplot()
# ggplot(data = avg, aes(x = factor(interaction(Pop, Food), levels = c("BON.H", "BON.L", "PSO.H", "PSO.L"), ordered = T), y = avgPropLiveArea, color = Timpoint, fill = Pop)) +
#   theme_classic() + 
#   geom_boxplot() + 
#   scale_color_manual(values = c("grey", "black")) + 
#   scale_fill_manual(values = c("steelblue", "gold"))
#   
# ggplot(data = metadata, aes(x = SL.mm, y = avgPropLiveArea, color = Timpoint)) + 
#   geom_point()
# ggplot(data = metadata, aes(x = DaysFoodTreat, y = avgPropLiveArea, color = Timpoint)) + 
#   geom_point()
# 
# ggplot(data = metadata, aes(x = Pop, y = avgPropLiveArea, fill = Food)) + 
#   geom_boxplot()
# ggplot(data = metadata, aes(x = Pop, y = avgPropLiveNum, fill = Food)) + 
#   geom_boxplot()
# ggplot(data = metadata, aes(x = Pop, y = avg.of.props.area, fill = Food)) + 
#   geom_boxplot()
# ggplot(data = metadata, aes(x = Pop, y = avg.of.props.num, fill = Food)) + 
#   geom_boxplot()

# T0 <- avg[which(avg$Timpoint == "0hr"),]
# T4 <- avg[which(avg$Timpoint == "4hr"),]
# 
# ggplot(data = T0, aes(x = Pop, y = avgPropLiveArea, fill = Food)) + 
#   geom_boxplot()
# ggplot(data = T4, aes(x = Pop, y = avgPropLiveArea, fill = Food)) + 
#   geom_boxplot()

# summary(glmer(avgPropLiveArea ~ Pop*Food + SL.mm  + DaysFoodTreat + Timpoint + (1 | IndivID), data = avg, family = "binomial", nAGQ = 0))
# summary(glmer(avgPropLiveArea ~ Pop*Food + DaysFoodTreat + Timpoint + (1 | IndivID), data = avg, family = "binomial", nAGQ = 0))
# summary(glmer(avgPropLiveArea ~ Pop*Food + Timpoint + (1 | IndivID), data = avg, family = "binomial", nAGQ = 0))

# metadata$propAliveNum <- metadata$liveNum / (metadata$liveNum + metadata$deadNum)
# metadata$propAliveArea <- metadata$liveArea / (metadata$liveArea + metadata$deadArea)
# 
# avgAliveNum <- c()
# avgAliveArea <- c()
# for (j in seq(2, nrow(metadata), by = 2)){
#   numVector <- rep(mean(c(metadata$propAliveNum[j], metadata$propAliveNum[j-1])), 2)
#   areaVector <- rep(mean(c(metadata$propAliveArea[j], metadata$propAliveArea[j-1])), 2)
#   append(avgAliveNum, numVector)
#   append(avgAliveArea, areaVector)
# }

# data <- setDT(metadata)
# avg <- data[,.(avgPropNumAlive = mean(propAliveNum), 
#                avgPropAreaAlive = mean(propAliveArea), 
#                Pop = unique(Pop), 
#                by = TrialID)]


# res.folder <- "/Volumes/NO NAME/2021_SpermViability/02_fiji_results/"
# 
# 
# raw.inputs.live <- list.files(path = raw.folder, pattern = "C0greenChannel")
# raw.inputs.dead <- list.files(path = raw.folder, pattern = "C1redChannel")
# 











