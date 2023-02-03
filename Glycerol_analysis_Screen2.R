##### Libraries #####

library(tidyverse)
library(cytominer)
library(pheatmap)
library(correlation)
library(missMDA)
library(FactoMineR)
library("corrplot")
library(factoextra)
library(cluster)
library(caret)
library(data.table)
library(DataExplorer)
library(tidymodels)
library(AICcmodavg)
library(multcomp)
library(devtools)
library(rstatix)
library(DescTools)
library(ggrepel)

##### Import Data ####

data_file <- file.choose()
data_samples <- read.csv(data_file) 
data_samples$Signal <- as.double(as.character(data_samples$Signal))
head(data_samples)
plot_intro(data_samples)
write.csv(data_samples, "data_samples_screen2.csv")

data_std <- file.choose()
std_all <- read.csv(data_std)
head(std_all)
std_all <- na.omit(std_all)
plot_intro


data_samples_re_file <- file.choose()
data_samples_re <- read.csv(data_samples_re_file)
data_samples_re$Signal <- as.double(as.character(data_samples_re$Signal))
head(data_samples_re)
data_samples_re <- na.omit(data_samples_re)
plot_intro(data_samples_re)
data_samples_re <- data_samples_re %>% rename(Image_Metadata_Well = Destination.Well)
data_samples_re <- data_samples_re[!duplicated(data_samples_re$Image_Metadata_Well), ]
data_samples_re <- data_samples_re %>%
  mutate(Image_Metadata_GeneID = recode(Image_Metadata_GeneID, PLIN1 = "siPLIN1"))


std_all_re_file <- file.choose()
std_all_re <- read.csv(std_all_re_file)
head(std_all_re)
std_all_re <- na.omit(std_all_re)
plot_intro(std_all_re)

ldh_file <- file.choose()
ldh <- read.csv(ldh_file)
ldh$Abs <- as.numeric(as.character(ldh$Abs))


##### Linear regression #####

lr <- std_all %>% group_by(Image_Metadata_Plate, Conc) %>% summarise(Signal = median(Signal))

lm <- lr %>% nest(data = -Image_Metadata_Plate) %>%
  mutate(lm = map(data, ~lm(Conc ~ 0 + Signal, data = .)),
         results = map(lm, glance)) %>% 
  unnest(results)  

r_sq <-  ggplot(lm, aes(x = factor(Image_Metadata_Plate), y = r.squared)) +
  geom_bar(stat = "identity") +
  labs(x = "Plate", y = expression(R^{2})) + 
  labs(title = "Glycerol Standard Curve - 2nd Screen") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) 
r_sq


##Remeasured

lr_re <- std_all_re %>% group_by(Image_Metadata_Plate, Conc) %>% summarise(Signal = median(Signal))

lm_re <- lr_re %>% nest(data = -Image_Metadata_Plate) %>%
  mutate(lm = map(data, ~lm(Conc ~ 0 + Signal, data = .)),
         results = map(lm, glance)) %>% 
  unnest(results)  

r_sq_re <-  ggplot(lm_re, aes(x = factor(Image_Metadata_Plate), y = r.squared)) +
  geom_bar(stat = "identity") +
  labs(x = "Plate", y = expression(R^{2})) + 
  labs(title = "Glycerol Standard Curve - 2nd Screen - Remeasured") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) 
r_sq_re


##### Data fit #####

data_fit <- data_samples %>% nest(data = -Image_Metadata_Plate)
data_fit <- data_fit %>% rename(data_samples = data)
data_fit <- left_join(lm, data_fit)
data_fit <- data_fit[,-c(4:15)]
fit <- data_fit %>%  mutate(Pred = map2_df(lm, data_samples, predict)) 

names_plate <- data_samples %>% group_by(Image_Metadata_Plate) %>% summarise()
Image_Metadata_Plate <- names_plate$Image_Metadata_Plate
names_well <- data_samples %>% group_by(Image_Metadata_Well) %>% summarise()
Image_Metadata_Well <- names_well$Image_Metadata_Well
annotation <- data_samples[,c(1:4)]

pred <- as.data.frame(t(fit[[5]]))
colnames(pred) <- c(Image_Metadata_Plate)
pred <- cbind(names_well, pred)
pred <- pred %>% pivot_longer(cols = 2:11,
                     names_to = "Image_Metadata_Plate",
                     values_to = "Conc")

data_fitted <- left_join(annotation, pred)
data_fitted <- right_join(data_samples, pred)


data_fitted$patient <- "P4"

data_fitted <- data_fitted %>% select(patient, everything())
data_fitted <- na.omit(data_fitted)



## Remeasured

data_fit_re <- data_samples_re %>% nest(data = -Image_Metadata_Plate)
data_fit_re <- data_fit_re %>% rename(data_samples_re = data)
data_fit_re <- left_join(lm_re, data_fit_re)
data_fit_re <- data_fit_re[,-c(4:15)]
fit_re <- data_fit_re %>%  mutate(Pred = map2_df(lm, data_samples_re, predict))
names_plate_re <- data_samples_re %>% group_by(Image_Metadata_Plate) %>% summarise()
Image_Metadata_Plate_re <- names_plate_re$Image_Metadata_Plate
names_well_re <- data_samples_re %>% group_by(Image_Metadata_Well) %>% summarise()
Image_Metadata_Well_re <- names_well_re$Image_Metadata_Well
annotation_re <- data_samples_re[,c(1:4)]
pred_re <- as.data.frame(t(fit_re[[5]]))
colnames(pred_re) <- c(Image_Metadata_Plate_re)
pred_re <- cbind(names_well_re, pred_re)
pred_re <- pred_re %>% pivot_longer(cols = 2,
                              names_to = "Image_Metadata_Plate",
                              values_to = "Conc")

data_fitted_re <- left_join(annotation_re, pred_re)
data_fitted_re <- right_join(data_samples_re, pred_re)
data_fitted_re$patient <- "P4"
data_fitted_re <- data_fitted_re %>% select(patient, everything())
data_fitted_re <- na.omit(data_fitted_re)


##### Compare 2 glycerol measurements #####

## Select second measurement data from both dfs

genes_1 <- data_samples %>% group_by(Image_Metadata_GeneID) %>% summarise
genes_1 <- genes_1$Image_Metadata_GeneID
genes_2 <- data_samples_re %>% group_by(Image_Metadata_GeneID) %>% summarise
genes_2 <- genes_2$Image_Metadata_GeneID

common_1 <- subset(data_fitted[,-c(1,2)], Image_Metadata_GeneID %in% genes_2)
common_1$measurement <- "1"
common_2 <- data_fitted_re[,2:6]
common_2$measurement <- "2"
common <- rbind(common_1, common_2)

## Correlation between 2 measurements

stats_common <- common %>% group_by(measurement, Image_Metadata_GeneID)  %>% 
  summarise(Min = min(Conc,na.rm = TRUE),
            Q1 = quantile(Conc,probs = .25,na.rm = TRUE),
            Median = median(Conc, na.rm = TRUE),
            Q3 = quantile(Conc,probs = .75,na.rm = TRUE),
            Max = max(Conc,na.rm = TRUE),
            Mean = mean(Conc, na.rm = TRUE),
            SD = sd(Conc, na.rm = TRUE),
            CV = abs((sd(Conc, na.rm = TRUE)/mean(Conc, na.rm = TRUE)*100)),
            MAD = mad(Conc, na.rm = TRUE),
            iqr = IQR(Conc, na.rm = TRUE),
            n = n(),
            Missing = sum(is.na(Conc)))

common$Image_Metadata_GeneID <- factor(common$Image_Metadata_GeneID)
common$Image_Metadata_GeneID <- relevel(common$Image_Metadata_GeneID, ref = "Neu. ctrl.")
common_bar <- common %>% group_by(measurement, Image_Metadata_GeneID) %>%
  ggplot(aes(x = Image_Metadata_GeneID, y = Conc, fill=measurement)) + geom_boxplot(position = "dodge") +
  scale_fill_manual(values = c("#7aa6c2", "#de425b")) +
  coord_flip() +
  labs(title = "Glycerol Concentration", subtitle = "Measurement 1 Vs. Measurement 2", fill = "Measurement") +
  xlab("Gene") + ylab("Glycerol Concentration (uM)") +
  theme_light() +
  theme(plot.title = element_text(size = 20),
        plot.subtitle =  element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
common_bar

##### Final df with measurements for 2nd meas #####

genes_3 <- genes_2[!genes_2 %in% c("Neu. ctrl.", "siPLIN1")]
data_fitted_2 <- data_fitted[ ! data_fitted$Image_Metadata_GeneID %in% genes_3, ]
data_fitted_2 <- rbind(data_fitted_2[,-2], data_fitted_re)


##### Descriptive Statistics #####

stats <- data_fitted_2 %>% group_by(Image_Metadata_GeneID)  %>% 
  summarise(Min = min(Conc,na.rm = TRUE),
            Q1 = quantile(Conc,probs = .25,na.rm = TRUE),
            Median = median(Conc, na.rm = TRUE),
            Q3 = quantile(Conc,probs = .75,na.rm = TRUE),
            Max = max(Conc,na.rm = TRUE),
            Mean = mean(Conc, na.rm = TRUE),
            SD = sd(Conc, na.rm = TRUE),
            CV = abs((sd(Conc, na.rm = TRUE)/mean(Conc, na.rm = TRUE)*100)),
            MAD = mad(Conc, na.rm = TRUE),
            iqr = IQR(Conc, na.rm = TRUE),
            n = n(),
            Missing = sum(is.na(Conc)))

ggplot(data_fitted_2, aes(x=Conc)) + geom_histogram(color="black", fill="grey") + 
  geom_vline(aes(xintercept=mean(Conc)), color="blue", linetype="dashed") +
  labs(title = "Raw data") + xlab("Count") + ylab("Concentration") +
  theme_light() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

data_fitted_iqr <- data_fitted_2 %>% left_join(stats) %>%
  group_by(patient, Image_Metadata_GeneID) %>%
  subset(Conc > (Q1 - 1.5*iqr) & Conc < (Q3 + 1.5*iqr))
data_fitted_iqr <- data_fitted_iqr[,-c(5,7:18)]

stats_iqr <- data_fitted_iqr %>% group_by(Image_Metadata_GeneID)  %>% 
  summarise(Min = min(Conc,na.rm = TRUE),
            Q1 = quantile(Conc,probs = .25,na.rm = TRUE),
            Median = median(Conc, na.rm = TRUE),
            Q3 = quantile(Conc,probs = .75,na.rm = TRUE),
            Max = max(Conc,na.rm = TRUE),
            Mean = mean(Conc, na.rm = TRUE),
            SD = sd(Conc, na.rm = TRUE),
            CV = abs((sd(Conc, na.rm = TRUE)/mean(Conc, na.rm = TRUE)*100)),
            MAD = mad(Conc, na.rm = TRUE),
            iqr = IQR(Conc, na.rm = TRUE),
            n = n(),
            Missing = sum(is.na(Conc)))
stats_iqr <- left_join(data_fitted_iqr[,1:5], stats_iqr)

ggplot(data_fitted_iqr, aes(x=Conc)) + geom_histogram(color="black", fill="grey") + 
  geom_vline(aes(xintercept=mean(Conc)), color="blue", linetype="dashed") +
  labs(title = "Raw data | IQR outlier elimination") + xlab("Count") + ylab("Concentration") +
  theme_light() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

high_cv <- data_fitted_iqr %>% left_join(stats_iqr) %>% filter(CV > 20) 
genes_highCV <- unique(high_cv$Image_Metadata_GeneID) 
genes_highCV

genes_highcv_plot <- ggplot(high_cv, aes(x=Image_Metadata_GeneID, y=CV)) + geom_col() +
  labs(title = "Screen 2  | Genes with CV > 20%") +
  ylab("CV%") + xlab("Gene") +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30),
        plot.subtitle =  element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10))
genes_highcv_plot


##### Normalization  #####

data_norm <- data_fitted_2[,-5] 
plot_intro(data_norm)

neg <- data_norm %>% filter(Image_Metadata_GeneID == "Neu. ctrl.")

data_norm <-normalize(population = data_norm ,variables = colnames(data_norm[5]), 
                      strata =  "Image_Metadata_Plate", 
                      sample = neg, operation = "robustize",
                      na.rm = T)

stats_norm <- data_norm %>% group_by(Image_Metadata_GeneID)  %>% 
  summarise(Min = min(Conc,na.rm = TRUE),
            Q1 = quantile(Conc,probs = .25,na.rm = TRUE),
            Median = median(Conc, na.rm = TRUE),
            Q3 = quantile(Conc,probs = .75,na.rm = TRUE),
            Max = max(Conc,na.rm = TRUE),
            Mean = mean(Conc, na.rm = TRUE),
            SD = sd(Conc, na.rm = TRUE),
            CV = abs((sd(Conc, na.rm = TRUE)/mean(Conc, na.rm = TRUE)*100)),
            MAD = mad(Conc, na.rm = TRUE),
            iqr = IQR(Conc, na.rm = TRUE),
            n = n(),
            Missing = sum(is.na(Conc)))

data_norm_iqr <- data_norm %>% left_join(stats_norm) %>%
  group_by(Image_Metadata_GeneID) %>% 
  subset(Conc > (Q1 - 1*iqr) & Conc < (Q3 + 1*iqr))

stats_norm_iqr<- data_norm_iqr %>% group_by(Image_Metadata_GeneID, Image_Metadata_Plate)  %>% 
  summarise(Min = min(Conc,na.rm = TRUE),
            Q1 = quantile(Conc,probs = .25,na.rm = TRUE),
            Median = median(Conc, na.rm = TRUE),
            Q3 = quantile(Conc,probs = .75,na.rm = TRUE),
            Max = max(Conc,na.rm = TRUE),
            Mean = mean(Conc, na.rm = TRUE),
            SD = sd(Conc, na.rm = TRUE),
            CV = abs((sd(Conc, na.rm = TRUE)/mean(Conc, na.rm = TRUE)*100)),
            MAD = mad(Conc, na.rm = TRUE),
            iqr = IQR(Conc, na.rm = TRUE),
            n = n(),
            Missing = sum(is.na(Conc)))
#saveRDS(data_norm_iqr, "data_norm_iqr_screen2.RDS")

high_cv_norm <- stats_norm_iqr %>% filter(CV > 20)
genes_highCV_norm <- unique(high_cv_norm$Image_Metadata_GeneID) 
genes_highCV_norm

genes_highcv_norm_plot <- ggplot(high_cv_norm, aes(x=Image_Metadata_GeneID, y=CV)) + geom_col() +
  labs(title = "Screen 1  | Genes with CV > 20%") +
  ylab("CV%") + xlab("Gene") +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 30),
        plot.subtitle =  element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10))
genes_highcv_norm_plot


data_norm_cv <- data_norm_iqr[ !data_norm_iqr$Image_Metadata_GeneID %in% genes_highCV, ]

##### Plots normalized data #####

## Batch effect

bef_norm <- data_norm  %>% ggplot(aes(x = patient, y = Conc, fill = patient)) + geom_violin() + geom_boxplot(width = 0.1) +
  ylim(-10,40) +
  labs(title = "Batch effect", subtitle = "Normalized data") +
  ylab("z-score") + xlab("Patient") +
  scale_fill_manual(values = c("#004b6e", "#bc5090", "#ffa600")) +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
bef_norm


bef_norm_iqr <- data_norm_iqr %>% ggplot(aes(x = patient, y = Conc, fill = patient)) + geom_violin() + geom_boxplot(width = 0.1) +
  ylim(-10,40) + 
  labs(title = "Batch effect", subtitle = "Normalized data | IQR outlier elimination") +
  ylab("z-score") + xlab("Patient") +
  scale_fill_manual(values = c("#004b6e", "#bc5090", "#ffa600")) +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
bef_norm_iqr


## Controls

c_norm <-  bind_rows(
  data_norm %>% filter(Image_Metadata_GeneID == "Neu. ctrl."),
  data_norm %>% filter(Image_Metadata_GeneID == "siPLIN1"),
  data_norm %>% filter(Image_Metadata_GeneID == "EBF1")
) 
c_norm$Image_Metadata_GeneID <- factor(c_norm$Image_Metadata_GeneID, levels = c("Neu. ctrl.", "siPLIN1", "EBF1"))
c_pl_norm <-  ggplot(c_norm, aes(x = Image_Metadata_GeneID, y = Conc, fill = Image_Metadata_GeneID)) + 
  geom_violin(aes(fill = Image_Metadata_GeneID)) + geom_boxplot(aes(fill = Image_Metadata_GeneID), width = 0.1) +
  ylim(-10,55) +
  labs(title = "Normalized data | Raw", subtitle = "Normalization per plate", fill = "Gene") +
  ylab("z-score") + xlab("Gene") +
  scale_x_discrete(labels = c("Negative control", "PLIN1", "EBF1")) +
  scale_fill_manual(values = c("#ebebeb", "#004b6e", "#5f7d96")) +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
c_pl_norm


c_norm_iqr <-  bind_rows(
  data_norm_iqr %>% filter(Image_Metadata_GeneID == "Neu. ctrl."),
  data_norm_iqr %>% filter(Image_Metadata_GeneID == "siPLIN1"),
  data_norm_iqr %>% filter(Image_Metadata_GeneID == "EBF1")
) 
c_norm_iqr$Image_Metadata_GeneID <- factor(c_norm_iqr$Image_Metadata_GeneID, levels = c("Neu. ctrl.", "siPLIN1", "EBF1"))
c_pl_norm_iqr <-  ggplot(c_norm_iqr, aes(x = Image_Metadata_GeneID, y = Conc)) + 
  geom_violin(aes(fill = Image_Metadata_GeneID)) + geom_boxplot(aes(fill = Image_Metadata_GeneID), width = 0.1) +
  ylim(-10,55) +
  labs(title = "Normalized data | IQR outlier elimination", subtitle = "Normalization per plate", fill = "Gene") +
  ylab("z-score") + xlab("Patient") +
  scale_x_discrete(name = "Gene", labels = c("Negative control", "PLIN1", "Ebf1")) +
  scale_fill_manual(values = c("#ebebeb", "#004b6e", "#5f7d96")) +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
c_pl_norm_iqr


## Histograms

ggplot(data_norm, aes(x=Conc)) + geom_histogram(aes(y=..density..), color="black", fill="grey") + 
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(Conc)), color="blue", linetype="dashed") +
  labs(title = "Normalized data") + xlab("Concentration") + ylab("Count") +
  theme_light() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))


ggplot(data_norm_iqr, aes(x=Conc)) + geom_histogram(aes(y=..density..), color="black", fill="grey") + 
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(Conc)), color="blue", linetype="dashed") +
  labs(title = "Normalized data | IQR outlier elimination") + xlab("Concentration") + ylab("Count") +
  theme_light() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))



##### LDH data ####

neg_ldh <- ldh %>% filter(Image_Metadata_GeneID == "Neu. ctrl.") %>% na.omit()

ldh_norm <-normalize(population = ldh ,variables = colnames(ldh[6]), 
                      strata =  "Image_Metadata_Plate", 
                      sample = neg_ldh, operation = "robustize",
                      na.rm = TRUE)

stats_norm_ldh <- ldh_norm %>% group_by(Image_Metadata_GeneID)  %>% 
  summarise(Min = min(Abs,na.rm = TRUE),
            Q1 = quantile(Abs,probs = .25,na.rm = TRUE),
            Median = median(Abs, na.rm = TRUE),
            Q3 = quantile(Abs,probs = .75,na.rm = TRUE),
            Max = max(Abs,na.rm = TRUE),
            Mean = mean(Abs, na.rm = TRUE),
            SD = sd(Abs, na.rm = TRUE),
            CV = abs((sd(Abs, na.rm = TRUE)/mean(Abs, na.rm = TRUE)*100)),
            MAD = mad(Abs, na.rm = TRUE),
            iqr = IQR(Abs, na.rm = TRUE),
            n = n(),
            Missing = sum(is.na(Abs)))

ldh_norm_iqr<- ldh_norm %>% left_join(stats_norm_ldh) %>%
  group_by(Image_Metadata_GeneID) %>% 
  subset(Abs > (Q1 - 0.5*iqr) & Abs < (Q3 + 0.5*iqr))

stats_norm_ldh_iqr<- ldh_norm_iqr %>% group_by(Image_Metadata_GeneID)  %>% 
  summarise(Min = min(Abs,na.rm = TRUE),
            Q1 = quantile(Abs,probs = .25,na.rm = TRUE),
            Median = median(Abs, na.rm = TRUE),
            Q3 = quantile(Abs,probs = .75,na.rm = TRUE),
            Max = max(Abs,na.rm = TRUE),
            Mean = mean(Abs, na.rm = TRUE),
            SD = sd(Abs, na.rm = TRUE),
            CV = abs((sd(Abs, na.rm = TRUE)/mean(Abs, na.rm = TRUE)*100)),
            MAD = mad(Abs, na.rm = TRUE),
            iqr = IQR(Abs, na.rm = TRUE),
            n = n(),
            Missing = sum(is.na(Abs)))


ggplot(ldh_norm_iqr, aes(x=Abs)) + geom_histogram(color="black", fill="grey") + 
  geom_vline(aes(xintercept=mean(Abs)), color="blue", linetype="dashed") +
  labs(title = "Raw data") + xlab("Count") + ylab("LDH") +
  theme_light() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))


ldh_norm_iqr$Image_Metadata_GeneID <- factor(ldh_norm_iqr$Image_Metadata_GeneID)
ldh_norm_iqr$Image_Metadata_GeneID <- relevel(ldh_norm_iqr$Image_Metadata_GeneID, ref = "Neu. ctrl.")
ldh_norm_col <- ldh_norm_iqr[ldh_norm_iqr$Image_Metadata_GeneID %in% c("Neu. ctrl.", top_genes_norm$Image_Metadata_GeneID), ] %>% 
  ggplot(aes(x = Image_Metadata_GeneID, y = Abs)) + geom_boxplot(fill = "#21918c") +
  xlab("Gene") + ylab("LDH") +
  geom_hline(yintercept=-3, linetype="dashed", color = "black") +
  geom_hline(yintercept=3, linetype="dashed", color = "black") +
  theme_classic() +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10))
ldh_norm_col
ggsave("LDH.png", ldh_norm_col, dpi = 300, units = "cm", width = 29.7, height = 21)
ggsave("LDH.eps", ldh_norm_col, dpi = 600, units = "cm", width = 29.7, height = 21)

ldh_norm_col_outl <- bind_rows(ldh_norm_iqr %>% group_by(Image_Metadata_GeneID) %>% filter(median(Abs) >= 5),
                               ldh_norm_iqr %>% group_by(Image_Metadata_GeneID) %>% filter(median(Abs) <= -5),
                               ldh_norm_iqr %>% filter(Image_Metadata_GeneID == "Neu. ctrl."),
                               ldh_norm_iqr %>% filter(Image_Metadata_GeneID == "siPLIN1"),) %>% 
  ggplot(aes(x = Image_Metadata_GeneID, y = Abs)) + geom_boxplot(fill="#004b6e") +
  labs(title = "LDH - Genes with rz >=|5|") + xlab("Count") + ylab("LDH") +
  theme_light() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
ldh_norm_col_outl
ggsave("LDH_genes_out.png", ldh_norm_col_outl, dpi = 300)

genes_ldh_out <- bind_rows(ldh_norm_iqr %>% group_by(Image_Metadata_GeneID) %>% filter(median(Abs) >= 5),
                           ldh_norm_iqr %>% group_by(Image_Metadata_GeneID) %>% filter(median(Abs) <= -5))
genes_ldh_out <- as.character(genes_ldh_out$Image_Metadata_GeneID)

##### ANOVA #####

#Norm_iqr

data_anova <- data_norm_cv[ ! data_norm_cv$Image_Metadata_GeneID %in% genes_ldh_out, ] #Change the data input
saveRDS(data_anova, "data_norm_cv_screen2.RDS")

#colnames(data_anova) <- c("Image_Metadata_GeneID", "Conc")

data_anova$Image_Metadata_GeneID <- factor(data_anova$Image_Metadata_GeneID)
data_anova$Image_Metadata_GeneID <- relevel(data_anova$Image_Metadata_GeneID, ref = "Neu. ctrl.")
levels(data_anova$Image_Metadata_GeneID)

anov_norm <- aov(Conc ~ Image_Metadata_GeneID, data=data_anova)
anov_norm
summary(anov_norm)

# Dunnett's test:
post_test <- glht(anov_norm,
                  linfct = mcp(Image_Metadata_GeneID = "Dunnett"))
post_test
sum_post_test <- summary(post_test)
sum_post_test

saveRDS(sum_post_test, "sum_post_test_normiqr_pl_screen2.RDS") # change save

##p value to dataframe

dun <- as.data.frame(sum_post_test[["test"]][["pvalues"]])
dun_rows <- levels(data_anova$Image_Metadata_GeneID) %>% as.data.frame()
dun_rows <- dun_rows[-1,] %>% as.data.frame()
colnames(dun_rows)[1] <- "Image_Metadata_GeneID"
colnames(dun)[1] <- 'p_adj'

data_anova_med <- data_anova %>% group_by(Image_Metadata_GeneID) %>%
  summarise(Conc = median(Conc))

dun <- cbind(dun_rows, dun) 
dun <- left_join(data_anova_med, dun) %>% na.omit()



##### Import ANOVA test RDS #####

#sum_dun_norm <- readRDS("sum_post_test_normiqr_pl_screen2.RDS")

dun_norm <- dun

dun_norm[dun_norm == 0] <- 0.0000000000001

dun_norm <- dun_norm %>% mutate(levels = case_when(Conc >= 2 & p_adj <= 0.05 ~ "High",
                                                 Conc <= -2 & p_adj <= 0.05 ~ "Low",
                                                 TRUE ~ "Unchanged"),
                              significance = case_when(
                                abs(Conc) >= 2 & p_adj <= 0.05 & p_adj > 0.01 ~ "p_adj 0.05", 
                                abs(Conc) >= 2 & p_adj <= 0.01 & p_adj > 0.001 ~ "p_adj 0.01",
                                abs(Conc) >= 2 & p_adj <= 0.001 ~ "p_adj 0.001", 
                                TRUE ~ "Unchanged"),
                              logp = -log10(p_adj))
head(dun_norm) 

top <- 25
top_genes_norm <- bind_rows(
  dun_norm %>% 
    filter(levels == 'High') %>% 
    arrange(p_adj, desc(abs(Conc))) %>% 
    head(top),
  dun_norm %>% 
    filter(levels == 'Low') %>% 
    arrange(p_adj, desc(abs(Conc))) %>% 
    head(top)
)
top_genes_norm
write.csv(top_genes_norm, "top_genes_norm.csv")

##### Plots ANOVA  #####

options(ggrepel.max.overlaps = Inf)
anova_plot <- ggplot(dun_norm, aes(x=Conc, y=-log(p_adj,10))) + geom_point(aes(color = levels)) +
  ylim(0,13)  + scale_color_manual(values = c("#004b6e", "#b4364a", "#ebebeb")) +
  geom_label_repel(data = top_genes_norm,
                   mapping = aes(Conc, -log(p_adj,10), label = Image_Metadata_GeneID),
                   size = 5) +
  labs(title = "Signicant alterations in glycerol levels", subtitle = "Normalization per plate", color = "Glycerol levels") +
  xlab("z-score") + ylab("-log10(adjusted p-value)") +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black") +
  geom_vline(xintercept=2, linetype="dashed", color = "black") +
  geom_vline(xintercept=-2, linetype="dashed", color = "black") +
  theme_light()  +
  theme(plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
anova_plot


options(ggrepel.max.overlaps = Inf)
anova_plot_p <- ggplot(dun_norm, aes(x=Conc, y=-log10(p_adj))) + geom_point(aes(color = significance)) +
  ylim(0,13) + scale_color_manual(values = c("#004b6e", "#5886a5", "#9dc6e0", "#c1e7ff")) +
  geom_label_repel(data = top_genes_norm,
                   mapping = aes(Conc, -log(p_adj,10), label = Image_Metadata_GeneID),
                   size = 2) +
  labs(title = "Signicant alterations in glycerol levels", subtitle = "Normalization per plate", color = "Glycerol levels") +
  xlab("z-score") + ylab("adjusted p-value") +
  geom_hline(yintercept=-log10(0.001), linetype="dashed", color = "black") +
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") +
  geom_vline(xintercept=2, linetype="dashed", color = "black") +
  geom_vline(xintercept=-2, linetype="dashed", color = "black") +
  theme_light()  +
  theme(plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
anova_plot_p


