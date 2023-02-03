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
library(data.table)
library(DataExplorer)
library(tidymodels)
library(AICcmodavg)
library(multcomp)
library(devtools)
library(rstatix)
library(DescTools)
library(ggrepel)
library(ggstatsplot)
library(GGally)


##### Import Data ####


data_file <- file.choose()
data_samples <- read.csv(data_file) 
data_samples$Signal <- as.double(as.character(data_samples$Signal))
head(data_samples)
plot_intro(data_samples)

data_std <- file.choose()
std_all <- read.csv(data_std)
head(std_all)
std_all <- na.omit(std_all)
plot_intro(std_all)


##### Linear regression #####

lr <- std_all %>% group_by(Image_Metadata_Plate, Conc) %>% summarise(Signal = median(Signal))

lm <- lr %>% nest(data = -Image_Metadata_Plate) %>%
  mutate(lm = map(data, ~lm(Conc ~ Signal, data = .)),
         results = map(lm, glance)) %>% 
  unnest(results)  

r_sq <-  ggplot(lm, aes(x = factor(Image_Metadata_Plate), y = r.squared)) +
  geom_bar(stat = "identity") +
  labs(x = "Plate", y = expression(R^{2})) + 
  labs(title = "Glycerol Standard Curve") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) 
r_sq


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
pred <- pred[,-46]
pred <- pred %>% pivot_longer(cols = 2:45,
                     names_to = "Image_Metadata_Plate",
                     values_to = "Conc")

data_fitted <- left_join(annotation, pred)


data_fitted <- right_join(data_samples, pred)
#data_fitted <- na.omit(data_fitted)

data_fitted$patient <- "P1"
data_fitted[data_fitted$ImageNumber>900,]$patient <-"P2"
data_fitted[data_fitted$ImageNumber>1860,]$patient <-"P3"
data_fitted <- data_fitted %>% select(patient, everything())


##### Descriptive Statistics #####

stats <- data_fitted %>% group_by(Image_Metadata_GeneID)  %>% 
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
            Var = var(Conc, na.rm = TRUE),
            n = n(),
            Missing = sum(is.na(Conc)))

ggplot(data_fitted, aes(x=Conc)) + geom_histogram(color="black", fill="grey") + 
  geom_vline(aes(xintercept=mean(Conc)), color="blue", linetype="dashed") +
  labs(title = "Raw data") + xlab("Count") + ylab("Concentration") +
  theme_light() +
  theme(plot.title = element_text(size = 40),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

data_fitted_iqr <- data_fitted %>% left_join(stats) %>%
  group_by(patient, Image_Metadata_GeneID) %>%
  subset(Conc > (Q1 - 1*iqr) & Conc < (Q3 + 1*iqr))
data_fitted_iqr <- data_fitted_iqr[,-c(6,8:19)]


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

high_cv <- stats_iqr %>% filter(CV > 20) %>% arrange(CV)
genes_highCV <- unique(high_cv$Image_Metadata_GeneID) 
genes_highCV

high_cv_toplot <- stats_iqr %>% group_by(Image_Metadata_GeneID) %>% summarise(CV = median(CV)) %>% arrange(CV)

genes_highcv_plot <- ggplot(high_cv_toplot, aes(x=reorder(Image_Metadata_GeneID, CV), y=CV)) + geom_col(fill = "#2a7c8e") +
  ylab("CV%") + xlab("Gene") +
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  geom_hline(yintercept=20, linetype="dashed", color = "#440154") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 18, family = "sans"),
        axis.title.y = element_text(size = 18, family = "sans"),
        axis.text.x = element_text(size = 3, angle = 90, hjust = 1, family = "sans"),
        axis.text.y = element_text(size = 18, family = "sans"))
genes_highcv_plot


data_fitted_cv <- data_fitted_iqr[ !data_fitted_iqr$Image_Metadata_GeneID %in% genes_highCV, ]

##### Batch effect plots #####

bef <- data_fitted %>% ggplot(aes(x = patient, y = Conc, fill = patient)) + geom_violin() + geom_boxplot(width = 0.1) +
  ylim(-10,100) +
  labs(title = "Raw data", subtitle = "Batch effect") +
  ylab("Concentration") + xlab("Patient") +
  scale_fill_manual(values = c("#004b6e", "#bc5090", "#ffa600")) +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
  bef

bef_iqr <- data_fitted_iqr %>% ggplot(aes(x = patient, y = Conc, fill = patient)) + geom_violin() + geom_boxplot(width = 0.1) +
  ylim(-10,100) +
  labs(title = "Raw data  | IQR outlier elimination", subtitle = "Batch effect") +
  ylab("Concentration") + xlab("Patient") +
  scale_fill_manual(values = c("#004b6e", "#bc5090", "#ffa600")) +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
bef_iqr


bef_cv <- data_fitted_cv %>% ggplot(aes(x = patient, y = Conc, fill = patient)) + geom_violin() + geom_boxplot(width = 0.1) +
  ylim(-10,100) +
  labs(title = "Raw data  | Genes with CV > 20% excluded", subtitle = "Batch effect") +
  ylab("Concentration") + xlab("Patient") +
  scale_fill_manual(values = c("#004b6e", "#bc5090", "#ffa600")) +
  theme_light() +
  theme(legend.position = "none",
        plot.title = element_text(size = 40),
        plot.subtitle =  element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
bef_cv


##### Normalization  #####

data_norm <- data_fitted  #data_fitted[,-6]
plot_intro(data_norm)

data_norm <- na.omit(data_norm)
plot_intro(data_norm)
data_norm <- data_norm[,-6]

neg <- data_norm %>% filter(Image_Metadata_GeneID == "neg")

data_norm <-normalize(population = data_norm ,variables = colnames(data_norm[6]), 
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

stats_norm_iqr<- data_norm_iqr %>% group_by(Image_Metadata_GeneID)  %>% 
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
#saveRDS(data_norm_iqr, "data_norm_iqr_screen1.RDS")

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


##### Replicate correlation #####

data_corr <- data_fitted_cv[,c(3,6)] %>% subset(Image_Metadata_GeneID!="neg" & Image_Metadata_GeneID!="Pos" & Image_Metadata_GeneID!="Ebf1")

data_neg <- data_fitted_cv[,c(3,4,6)] %>% subset(Image_Metadata_GeneID == "neg")
data_neg$Image_Metadata_GeneID <- paste(data_neg$Image_Metadata_GeneID, data_neg$Image_Metadata_Plate)
data_neg <- data_neg[,-2]
data_neg <- data_neg %>% group_by(Image_Metadata_GeneID) %>%
  mutate(Replicate = row_number(Image_Metadata_GeneID),
         Replicate2 = case_when(Replicate == 1 ~ 1,
                                Replicate == 2 ~ 1,
                                Replicate == 3 ~ 1,
                                Replicate == 4 ~ 2,
                                Replicate == 5 ~ 2,
                                Replicate == 6 ~ 2))
data_neg$Image_Metadata_GeneID <- paste(data_neg$Image_Metadata_GeneID, data_neg$Replicate2)
data_neg <- data_neg[,-c(3:4)]

data_pos <- data_fitted_cv[,c(3,4,6)] %>% subset(Image_Metadata_GeneID == "Pos")
data_pos$Image_Metadata_GeneID <- paste(data_pos$Image_Metadata_GeneID, data_pos$Image_Metadata_Plate)
data_pos <- data_pos[,-2]

data_ebf <- data_fitted_cv[,c(3,4,6)] %>% subset(Image_Metadata_GeneID == "Ebf1")
data_ebf$Image_Metadata_GeneID <- paste(data_ebf$Image_Metadata_GeneID, data_ebf$Image_Metadata_Plate)
data_ebf <- data_ebf[,-2]

data_corr <- rbind(data_corr, data_pos, data_ebf, data_neg)

data_corr <- data_corr %>% group_by(Image_Metadata_GeneID) %>%
  mutate(Replicate = row_number(Image_Metadata_GeneID)) %>% ungroup()
data_corr <- data_corr %>% select(Replicate, everything()) %>% 
  pivot_wider(names_from = Replicate, values_from = Conc) ## !!!!! CHANGE VALUES FROM
colnames(data_corr) <- c("Gene ID", "Replicate 1", "Replicate 2", "Replicate 3")

corr <- correlation(data_corr[,2:4], method = c("pearson"))

rep1vs2 <- ggscatterstats(data = data_corr[,2:4], x = "Replicate 1", y = "Replicate 2", bf.message = FALSE,
                          marginal = T, 
                          xfill = "#21918c",
                          yfill = "#5ec962",
                          smooth.line.args = list(size = 1.5, color = "#3b528b", method = "lm", formula = y ~ x,
                                                  na.rm = T)) 
rep1vs2

rep2vs3 <- ggscatterstats(data = data_corr[,2:4], x = "Replicate 2", y = "Replicate 3", bf.message = FALSE,
                          marginal = T, 
                          xfill = "#21918c",
                          yfill = "#5ec962",
                          smooth.line.args = list(size = 1.5, color = "#3b528b", method = "lm", formula = y ~ x,
                                                  na.rm = TRUE))
rep2vs3

rep1vs3 <- ggscatterstats(data = data_corr[,2:4], x = "Replicate 1", y = "Replicate 3", bf.message = FALSE,
                          marginal = T, 
                          xfill = "#21918c",
                          yfill = "#5ec962",
                          smooth.line.args = list(size = 1.5, color = "#3b528b", method = "lm", formula = y ~ x,
                                                  na.rm = TRUE))
rep1vs3


ggplot(data_corr, aes(x = `Replicate 1`, y = `Replicate 2`)) + 
  geom_point() +
  geom_bin_2d(bins=120) +
  geom_smooth(method = "lm", col = "red") +
  scale_fill_continuous(type = "viridis") +
  #scale_fill_continuous(type = "viridis") +
  theme_classic()

ggplot(data_corr, aes(x = `Replicate 2`, y = `Replicate 3`)) + 
  geom_point() +
  geom_bin2d(bins=100) +
  geom_smooth(method = "lm", col = "red") +
  scale_fill_continuous(type = "viridis") +
  theme_classic()

ggplot(data_corr, aes(x = `Replicate 1`, y = `Replicate 3`)) + 
  geom_point() +
  geom_bin2d(bins=100) +
  geom_smooth(method = "lm", col = "red") +
  scale_fill_continuous(type = "viridis") +
  theme_classic()


ggplot(data_corr) +
  geom_point(aes(x = `Replicate 1`, y = `Replicate 2`, colour = "red")) +
  geom_point(aes(x = `Replicate 2`, y = `Replicate 3`, colour = "blue")) +
  geom_point(aes(x = `Replicate 1`, y = `Replicate 3`, colour = "green")) +
  geom_smooth(aes(x = `Replicate 1`, y = `Replicate 2`, col = "red"), method = "lm", col = "red") +
  geom_smooth(aes(x = `Replicate 2`, y = `Replicate 3`, col = "blue"), method = "lm", col = "blue") +
  geom_smooth(aes(x = `Replicate 1`, y = `Replicate 3`, col = "green"), method = "lm", col = "green") +
  #ylim(0, 120) +
  theme_classic() +
  theme(legend.title = element_blank())
ggsave("rep_corr_all.eps", dpi = 300)


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
ggsave("bef_norm_iqr.png", plot = bef_norm_iqr, dpi = 300)

## Controls

c_norm <-  bind_rows(
  data_norm %>% filter(Image_Metadata_GeneID == "neg"),
  data_norm %>% filter(Image_Metadata_GeneID == "Pos"),
  data_norm %>% filter(Image_Metadata_GeneID == "Ebf1")
) 
c_norm$Image_Metadata_GeneID <- factor(c_norm$Image_Metadata_GeneID, levels = c("neg", "Pos", "Ebf1"))
c_pl_norm <-  ggplot(c_norm, aes(x = Image_Metadata_GeneID, y = Conc, fill = Image_Metadata_GeneID)) + 
  geom_violin(aes(fill = Image_Metadata_GeneID)) + geom_boxplot(aes(fill = Image_Metadata_GeneID), width = 0.1) +
  ylim(-10,40) +
  labs(title = "Normalized data | Raw", subtitle = "Normalization per plate", fill = "Gene") +
  ylab("z-score") + xlab("Gene") +
  scale_x_discrete(labels = c("Negative control", "PLIN1", "Ebf1")) +
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
ggsave("c_norm.png", c_pl_norm, dpi = 300)

c_norm_iqr <-  bind_rows(
  data_norm_iqr %>% filter(Image_Metadata_GeneID == "neg"),
  data_norm_iqr %>% filter(Image_Metadata_GeneID == "Pos"),
  data_norm_iqr %>% filter(Image_Metadata_GeneID == "Ebf1")
) 
c_norm_iqr$Image_Metadata_GeneID <- factor(c_norm_iqr$Image_Metadata_GeneID, levels = c("neg", "Pos", "Ebf1"))
c_pl_norm_iqr <-  ggplot(c_norm_iqr, aes(x = Image_Metadata_GeneID, y = Conc)) + 
  geom_violin(aes(fill = Image_Metadata_GeneID)) + geom_boxplot(aes(fill = Image_Metadata_GeneID), width = 0.1) +
  ylim(-10,40) +
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
ggsave("c_norm_iqr.png", c_pl_norm_iqr, dpi = 300)

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
ggsave("histo_norm.png")

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
ggsave("histo_normiqr.png")

##### ANOVA #####

#Norm_iqr

data_anova <- data_norm_cv #Change the data input
#colnames(data_anova) <- c("Image_Metadata_GeneID", "Conc")

data_anova$Image_Metadata_GeneID <- factor(data_anova$Image_Metadata_GeneID)
data_anova$Image_Metadata_GeneID <- relevel(data_anova$Image_Metadata_GeneID, ref = "neg")
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

saveRDS(dun, "dun_norm_cv.RDS") # change save

##### Import ANOVA test RDS #####

#sum_dun_norm <- readRDS("sum_post_test_normiqr_pl.RDS")

dun_norm <-  readRDS("Results_Screen1_ANOVA/dun_norm_cv.RDS") #readRDS("Results_Glycerol_RDS/False/dun_norm_iqr.RDS")

dun_norm[dun_norm == 0] <- 0.0000000000001
dun_norm <- dun_norm %>% mutate(Image_Metadata_GeneID = recode(Image_Metadata_GeneID, Pos = "siPLIN1", neg = "Neu. ctrl.",
                                                     HDAC7A = "HDAC7", 
                                                     P8 = "NUPR1", 	
                                                     XPMC2H = "REXO4", 	
                                                     TTRAP = "TDP2",	
                                                     ZBTB7 = "ZBTB7A")) 

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

top <- 20
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
write.csv(top_genes_norm, "top_genes_screen1.csv")


##### Plots ANOVA #####

options(ggrepel.max.overlaps = Inf)
anova_plot <- ggplot(dun_norm, aes(x=Conc, y=-log(p_adj,10))) + geom_point(aes(color = levels)) +
  ylim(0,13)  + scale_color_manual(values = c("#2a7c8e", "#fce61f", "#ebebeb")) +
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






