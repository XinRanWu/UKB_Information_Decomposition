# red & syn

install.packages(c("lmerTest", "dplyr", "tidyr", "ggplot2", "purrr", "stringr"))

library(purrr)
library(lmerTest)
library(dplyr)
library(tidyr)
library(ggplot2)

library(stringr)
library(readr)

################################################################################
#                              Aging Trajectory
################################################################################


UKB_BasicMRI <- read_csv("/public/home/zhangjie/ZJLab/UK
                         Biobank_Project/data/table/UKB_BrainMRICov.txt")
UKB_RedSynPC_2_0 = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/ukb_RedSynPCYeo7n_2_0.txt")
UKB_RedSynNC_2_0 = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/ukb_RedSynNCYeo7n_2_0.txt")
UKB_RedSyn_2_0 = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/ukb_RedSynYeo7n_2_0.txt")
UKB_AllVars <- read_csv("ZJLab/UKBiobank_Project/data/table/UKB_AllVars_NumCol.csv")

names(UKB_RedSyn_2_0)[1] <- "eid"
UKB_PID_Yeo7n = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Yeo7n.txt",na = "NaN")
data <- merge(merge(UKB_BasicMRI,UKB_PID_Yeo7n),UKB_AllVars[,c("eid","Fluid_intelligence_score_2_0")])
data$age <- data$AgeAttend_2_0
data$site <- data$Centre_2_0

data <- merge(merge(UKB_BasicMRI,UKB_RedSyn_2_0),UKB_AllVars[,c("eid","Fluid_intelligence_score_2_0")])
data$age <- data$AgeAttend_2_0
data$site <- data$Centre_2_0

# data <- read.csv("your_data_file.csv")
X_vars <- names(UKB_RedSyn_2_0)[2:ncol(UKB_RedSyn_2_0)]
dependent_var <- "Fluid_intelligence_score_2_0"
covariates <- c("Sex_0_0","BMI_2_0","Race_1","Race_2","Race_3","Race_4" ,"HeadMotion_2_0","SNR_2_0")

window_size <- 10  # 5 yrs Win Size
step_size <- 5    # 1 yrs Step 

min_age <- floor(min(data$age))
max_age <- ceiling(max(data$age))

time_windows <- seq(min_age, max_age - window_size + 1, by = step_size)
fit_models <- function(data, start_age) {
  end_age <- start_age + window_size - 1
  window_label <- paste0("time window ", start_age - min_age + 1)
  age_range <- paste0(start_age, "~", end_age)
  
  window_data <- data %>%
    filter(age >= start_age & age <= end_age)

  results <- tibble(
    time_window = window_label,
    age_range = age_range,
    variable = x_var,
    dependent_var = dependent_var,
    Estimate = NA,
    df = NA,
    t_value = NA,
    p_value = NA,
    CI_lower = NA,
    CI_upper = NA
  )
  
  n = 0
  for (x_var in X_vars) {
    n = n + 1
    formula <- as.formula(
      paste(dependent_var, "~ age +", paste(c(covariates, x_var), collapse = " + "), "+ (1|site)")
    )
    model <- lmerTest::lmer(formula, data = window_data)
    model_summary <- summary(model)
    ci <- try(confint(model, parm = x_var, level = 0.95), silent = TRUE)
    ci_lower <- ci[1]
    ci_upper <- ci[2]
    fixed_effect <- coef(model_summary)[x_var, ]
    result <- tibble(
      time_window = window_label,
      age_range = age_range,
      variable = x_var,
      dependent_var = dependent_var,
      Estimate = fixed_effect["Estimate"],
      df = fixed_effect["df"],
      t_value = fixed_effect["t value"],
      p_value = fixed_effect["Pr(>|t|)"],
      CI_lower = ci_lower,
      CI_upper = ci_upper
    )
    results[n,] <- result
  }
  return(results)
}

all_results = list()
for (i in 1:length(time_windows)){
  all_results[[i]] <- fit_models(data, time_windows[i])
}
all_results <- bind_rows(all_results)

head(all_results)

# write.csv(all_results, "model_results.csv", row.names = FALSE)

all_results <- all_results %>%
  mutate(
    start_age = as.numeric(str_extract(age_range, "^[0-9]+")),
    mid_age = start_age + (window_size - 1) / 2
  )
color_palette <- c("red", "blue", "green", "purple", "orange", "brown", "pink") 


all_results$label <- NA

all_results$label[grep("RedYeo7n_In([1-7])$", all_results$variable)] <- 'RIn'
all_results$label[grep("RedYeo7n_Btw([0-1])$", all_results$variable)] <- 'RBtw_L'
all_results$label[grep("RedYeo7n_Btw([2-9]|1[0-1])$", all_results$variable)] <- 'RBtw_LH'
all_results$label[grep("RedYeo7n_Btw(1[2-9]|2[0-1])$", all_results$variable)] <- 'RBtw_H'

all_results$label[grep("SynYeo7n_In([1-7])$", all_results$variable)] <- 'SIn'
all_results$label[grep("SynYeo7n_Btw([0-1])$", all_results$variable)] <- 'SBtw_L'
all_results$label[grep("SynYeo7n_Btw([2-9]|1[0-1])$", all_results$variable)] <- 'SBtw_LH'
all_results$label[grep("SynYeo7n_Btw(1[2-9]|2[0-1])$", all_results$variable)] <- 'SBtw_H'



ggplot(all_results[all_results$label=="RBtw_LH",], aes(x = mid_age, y = t_value, color = variable)) +
  geom_line() +
  geom_point() +
  # geom_errorbar(aes(ymin = (Estimate - CI_upper)/abs(CI_upper - CI_lower), 
  #                   ymax = (Estimate + CI_upper)/abs(CI_upper - CI_lower)), 
  #               width = 0.5) + 
  # #scale_color_manual(values = color_palette) +
  labs(
    title = "Sliding Window Analysi",
    x = "Age",
    y = "t value",
    color = "Yeo 7 network"
  ) +
  theme_minimal()

ggplot(all_results[str_detect(all_results$variable,"SynYeo7n_"),], aes(x = mid_age, y = t_value, group = variable, color = label)) +
  geom_line() +
  geom_point() +
  # geom_errorbar(aes(ymin = (Estimate - CI_upper)/abs(CI_upper - CI_lower), 
  #                   ymax = (Estimate + CI_upper)/abs(CI_upper - CI_lower)), 
  #               width = 0.5) + 
  # #scale_color_manual(values = color_palette) +
  labs(
    title = "Sliding Window Analysis",
    x = "Age",
    y = "t value",
    color = "Yeo 7 network"
  ) +
  theme_minimal()


# adni3



library(mgcv)
library(ggplot2)
library(ggpubr)

model <- gamm(RedYeo7n_Total ~ s(AgeAttend_2_0), data = data)

pred_data <- data.frame(AgeAttend_2_0 = seq(min(data$AgeAttend_2_0), max(data$AgeAttend_2_0), length.out = 100))
pred <- predict(model$gam, newdata = pred_data, se.fit = TRUE) # »ñÈ¡Ô¤²âÖµºÍ±ê×¼Îó²î
pred_data$pred <- pred$fit
pred_data$se.fit <- pred$se.fit  

ggplot(pred_data, aes(x = AgeAttend_2_0, y = pred_data$pred)) +
  # geom_point(aes(color = ..density..), stat = "density_2d", size = 2, alpha = 0.5) + 
  scale_color_viridis_c() +  
  ggplot2::geom_density_2d(colour = "gray") +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(x = AgeAttend_2_0, ymin = pred - se.fit, ymax = pred + se.fit), 
              fill = "blue", alpha = 0.2) +
  labs(title = "Trajectory of Redundary by Age", x = "Age", y = "Redundary") +
  theme_minimal()

ggplot(pred_data, aes(x = AgeAttend_2_0, y = pred_data$pred)) +
  # geom_point(aes(color = ..density..), stat = "density_2d", size = 2, alpha = 0.5) + 
  scale_color_viridis_c() +  
  geom_ribbon(aes(x = AgeAttend_2_0, ymin = pred - se.fit, ymax = pred + se.fit), 
              fill = "blue", alpha = 0.2) +
  ggplot2::geom_density_2d_filled() +
  geom_line(color = "blue", size = 1) +
  labs(title = "Trajectory of Redundary by Age", x = "Age", y = "Redundary") +
  theme_minimal()


install.packages("gratia")
library(gratia)
# 1. ÑµÁ·Ä£ÐÍ
model <- gamm(RedYeo7n_Total ~ s(AgeAttend_2_0), data = data)

# 2. ¼ÆËãÒ»½×µ¼Êý£¨±ä»¯ÂÊ£©
deriv_data <- derivatives(model$gam, select = "s(AgeAttend_2_0)", interval = "confidence")

ggplot(deriv_data, aes(x = AgeAttend_2_0, y = .derivative)) +
  geom_line(color = "blue", linewidth = 1) +  
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "blue", alpha = 0.2) + 
  labs(title = "Rate of Change in Redundary by Age", 
       x = "Age", 
       y = "Rate of Change") +
  theme_minimal()


model <- gamm(SynYeo7n_Total ~ s(AgeAttend_2_0), data = data)

pred_data <- data.frame(AgeAttend_2_0 = seq(min(data$AgeAttend_2_0), max(data$AgeAttend_2_0), length.out = 100))
pred <- predict(model$gam, newdata = pred_data, se.fit = TRUE) # »ñÈ¡Ô¤²âÖµºÍ±ê×¼Îó²î
pred_data$pred <- pred$fit
pred_data$se.fit <- pred$se.fit  

ggplot(pred_data, aes(x = AgeAttend_2_0, y = pred_data$pred)) +
  # geom_point(aes(color = ..density..), stat = "density_2d", size = 2, alpha = 0.5) + 
  # scale_color_viridis_c() +  
  geom_line(color = "red", size = 1) +
  geom_ribbon(aes(x = AgeAttend_2_0, ymin = pred - se.fit, ymax = pred + se.fit), 
              fill = "red", alpha = 0.2) +
  ggplot2::geom_density_2d(colour = "gray") +
  labs(title = "Trajectory of Synergy by Age", x = "Age", y = "Synergy") +
  theme_minimal()


model <- gamm(SynYeo7n_Total ~ s(AgeAttend_2_0), data = data)

# 2. ¼ÆËãÒ»½×µ¼Êý£¨±ä»¯ÂÊ£©
deriv_data <- derivatives(model$gam, select = "s(AgeAttend_2_0)", interval = "confidence")

ggplot(deriv_data, aes(x = AgeAttend_2_0, y = .derivative)) +
  geom_line(color = "red", linewidth = 1) +  
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "red", alpha = 0.2) + 
  labs(title = "Rate of Change in Redundary by Age", 
       x = "Age", 
       y = "Rate of Change") +
  theme_minimal()







get_col_indices <- function(col_id) {
  idx <- c()  # ´æ´¢Õâ¸öÁÐÏà¹ØµÄÏÂÈý½ÇÎ»ÖÃ
  count <- 1
  for (i in 2:7) {
    for (j in 1:(i - 1)) {
      if (j == col_id || i == col_id) {
        idx <- c(idx, count)
      }
      count <- count + 1
    }
  }
  return(idx)
}

network_names <- c("Visual", "Somatomotor", "DorsalAttention", "VentralAttention",
                   "Limbic", "Frontoparietal", "Default")

data2 <- data.frame(AgeAttend_2_0 = data$AgeAttend_2_0)

for (i in 1:7) {
  col_idx <- get_col_indices(i)
  between_names <- paste0("SynYeo7n_Btw", col_idx)
  vals <- cbind(data[[paste0("SynYeo7n_In", i)]], data[, between_names])
  data2[[network_names[i]]] <- rowMeans(vals, na.rm = TRUE)
}
data2[["Total"]] <- data[, "SynYeo7n_Total"]

for (i in 1:7) {
  col_idx <- get_col_indices(i)
  between_names <- paste0("RedYeo7n_Btw", col_idx)
  vals <- cbind(data[[paste0("RedYeo7n_In", i)]], data[, between_names])
  data2[[network_names[i]]] <- rowMeans(vals, na.rm = TRUE)
}
data2[["Total"]] <- data[, "RedYeo7n_Total"]




data2 <- merge(UKB_BasicMRI,UKB_PID_Yeo7n[,c(1:94,102:108)])#95:101
names(data2)[95:101] <- c("Visual", "Somatomotor", "DorsalAttention", "VentralAttention",
                   "Limbic", "Frontoparietal", "Default")


Yeo7Color = c(rgb(0.471,0.0710,0.522),
              rgb(0.275,0.510,0.706),
              rgb(0,0.463,0.0550),
              rgb(0.769,0.224,0.976),
              rgb(0.863,0.973,0.639),
              rgb(0.902,0.576,0.129),
              rgb(0.804,0.239,0.306))



ggplot(data2, aes(x = AgeAttend_2_0, group = Sex_0_0)) +
  geom_smooth(aes(y = Visual, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[1], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = Somatomotor, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[2], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = DorsalAttention, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[3], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = VentralAttention, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[4], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = Limbic, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[5], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = Frontoparietal, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[6], se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = Default, linetype = factor(Sex_0_0)), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[7], se = TRUE, alpha = 0.2) +
  labs(title = "GAM", x = "Age", y = " ") +
  theme_bw()



ggplot(data2, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = Visual), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[1], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Somatomotor), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[2], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = DorsalAttention), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[3], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = VentralAttention), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[4], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Limbic), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[5], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Frontoparietal), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[6], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Default), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[7], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Total), method = "gam", formula = y ~ s(x),
              color = "blue", se = TRUE, alpha = 0.5) +
  labs(title = "GAM", x = "Age", y = " ") +
  theme_bw()

ggplot(data2, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = Visual), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[1], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Somatomotor), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[2], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = DorsalAttention), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[3], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = VentralAttention), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[4], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Limbic), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[5], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Frontoparietal), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[6], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Default), method = "gam", formula = y ~ s(x),
              color = Yeo7Color[7], se = TRUE, alpha = 0.2, linetype = "33") +
  geom_smooth(aes(y = Total), method = "gam", formula = y ~ s(x),
              color = "red", se = TRUE, alpha = 0.6, size = 1.5) +
  labs(title = "GAM", x = "Age", y = " ") +
  theme_bw()






library(mgcv)
library(ggplot2)
library(gridExtra)

# ±äÁ¿Ãû
inner_vars <- paste0("RedYeo7n_In", 1:7)
btw_vars <- paste0("RedYeo7n_Btw", 1:28)

covariates <- "Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + TDI_0_0 + Vol_WB_TIV_2_0"

# ´´½¨¿ÕµÄ»æÍ¼ÁÐ±í
plots <- list()
x_min <- min(data$AgeAttend_2_0)
x_max <- max(data$AgeAttend_2_0)
# Éú³É 7x7 µÄ×ÓÍ¼
for (i in 1:7) {
  for (j in 1:7) {
    if (i == j) {
      # ¶Ô½ÇÏß: RedYeo7n_In ±äÁ¿µÄÇúÏß
      model <- gamm(as.formula(paste(inner_vars[i], "~ s(AgeAttend_2_0) +", covariates)), data = data)
      pred_data <- data.frame(
        AgeAttend_2_0 = seq(x_min, x_max, length.out = 100),
        Sex_0_0 = mean(data$Sex_0_0, na.rm = TRUE),
        BMI_2_0 = mean(data$BMI_2_0, na.rm = TRUE),
        HeadMotion_2_0 = mean(data$HeadMotion_2_0, na.rm = TRUE),
        Race_1 = mean(data$Race_1, na.rm = TRUE),
        Race_2 = mean(data$Race_2, na.rm = TRUE),
        Race_3 = mean(data$Race_3, na.rm = TRUE),
        Race_4 = mean(data$Race_4, na.rm = TRUE),
        TDI_0_0 = mean(data$TDI_0_0, na.rm = TRUE),
        Vol_WB_TIV_2_0 = mean(data$Vol_WB_TIV_2_0, na.rm = TRUE)
      )
      pred <- predict(model$gam, newdata = pred_data, se.fit = TRUE)
      pred_data$pred <- pred$fit
      pred_data$se.fit <- pred$se.fit
      
      p <- ggplot(pred_data, aes(x = AgeAttend_2_0, y = pred)) +
        geom_ribbon(aes(ymin = pred - se.fit, ymax = pred + se.fit), fill = "blue", alpha = 0.2) +
        geom_line(color = "blue", size = 1) +
        labs(title = "", x = "", y = "") +
        theme_minimal() +
        theme(plot.title = element_text(size = 8), axis.title.x = element_text(size = 8))
      # xlim(x_min, x_max) + ylim(y_min, y_max)
    } else if (i > j) {
      # ÏÂÈý½Ç: RedYeo7n_Btw ±äÁ¿µÄÇúÏß
      index <- ((j - 1) * 7) - (((j - 1) * j) / 2) + (i - j)  # ¼ÆËã±äÁ¿Ë÷Òý
      if (index <= length(btw_vars)) {
        model <- gamm(as.formula(paste(btw_vars[index], "~ s(AgeAttend_2_0) +", covariates)), data = data)
        pred_data <- data.frame(
          AgeAttend_2_0 = seq(x_min, x_max, length.out = 100),
          Sex_0_0 = mean(data$Sex_0_0, na.rm = TRUE),
          BMI_2_0 = mean(data$BMI_2_0, na.rm = TRUE),
          HeadMotion_2_0 = mean(data$HeadMotion_2_0, na.rm = TRUE),
          Race_1 = mean(data$Race_1, na.rm = TRUE),
          Race_2 = mean(data$Race_2, na.rm = TRUE),
          Race_3 = mean(data$Race_3, na.rm = TRUE),
          Race_4 = mean(data$Race_4, na.rm = TRUE),
          TDI_0_0 = mean(data$TDI_0_0, na.rm = TRUE),
          Vol_WB_TIV_2_0 = mean(data$Vol_WB_TIV_2_0, na.rm = TRUE)
        )
        pred <- predict(model$gam, newdata = pred_data, se.fit = TRUE)
        pred_data$pred <- pred$fit
        pred_data$se.fit <- pred$se.fit
        
        p <- ggplot(pred_data, aes(x = AgeAttend_2_0, y = pred)) +
          geom_ribbon(aes(ymin = pred - se.fit, ymax = pred + se.fit), fill = "blue", alpha = 0.2) +
          geom_line(color = "blue", size = 1) +
          labs(title = "", x = "", y = "") +
          theme_minimal() +
          theme(plot.title = element_text(size = 8), axis.title.x = element_text(size = 8))
        # xlim(x_min, x_max) + ylim(y_min, y_max)
      } else {
        p <- ggplot() + theme_void()  # ¿Õ°×Í¼
      }
    } else {
      p <- ggplot() + theme_void()  # ÉÏÈý½ÇÁô¿Õ
    }
    
    plots[[length(plots) + 1]] <- p
  }
}

# »æÖÆ 7x7 Íø¸ñ
grid.arrange(grobs = plots, nrow = 7, ncol = 7)



# 1500 1200




# ´´½¨Ò»¸ö7x7µÄ¾ØÕó£¬²¢ÓÃNAÌî³ä
mat <- matrix(NA, nrow = 7, ncol = 7)

# »ñÈ¡ÏÂÈý½ÇÇøÓòµÄË÷Òý
lower_tri_indices <- which(lower.tri(mat), arr.ind = TRUE)

# ÌáÈ¡ÐÐË÷ÒýºÍÁÐË÷Òý
row_indices <- lower_tri_indices[, 1]
col_indices <- lower_tri_indices[, 2]

# ´òÓ¡Ë÷Òý
row_indices
col_indices
# ±äÁ¿Ãû
inner_vars <- paste0("SynYeo7n_In", 1:7)
btw_vars <- paste0("SynYeo7n_Btw", 1:28)

# ´´½¨¿ÕµÄ»æÍ¼ÁÐ±í
covariates <- "Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + TDI_0_0 + Vol_WB_TIV_2_0"

# ´´½¨¿ÕµÄ»æÍ¼ÁÐ±í
plots <- list()
x_min <- min(data$AgeAttend_2_0)
x_max <- max(data$AgeAttend_2_0)
# Éú³É 7x7 µÄ×ÓÍ¼
for (i in 1:7) {
  for (j in 1:7) {
    if (i == j) {
      # ¶Ô½ÇÏß: RedYeo7n_In ±äÁ¿µÄÇúÏß
      model <- gamm(as.formula(paste(inner_vars[i], "~ s(AgeAttend_2_0) +", covariates)), data = data)
      pred_data <- data.frame(
        AgeAttend_2_0 = seq(x_min, x_max, length.out = 100),
        Sex_0_0 = mean(data$Sex_0_0, na.rm = TRUE),
        BMI_2_0 = mean(data$BMI_2_0, na.rm = TRUE),
        HeadMotion_2_0 = mean(data$HeadMotion_2_0, na.rm = TRUE),
        Race_1 = mean(data$Race_1, na.rm = TRUE),
        Race_2 = mean(data$Race_2, na.rm = TRUE),
        Race_3 = mean(data$Race_3, na.rm = TRUE),
        Race_4 = mean(data$Race_4, na.rm = TRUE),
        TDI_0_0 = mean(data$TDI_0_0, na.rm = TRUE),
        Vol_WB_TIV_2_0 = mean(data$Vol_WB_TIV_2_0, na.rm = TRUE)
      )
      pred <- predict(model$gam, newdata = pred_data, se.fit = TRUE)
      pred_data$pred <- pred$fit
      pred_data$se.fit <- pred$se.fit
      
      p <- ggplot(pred_data, aes(x = AgeAttend_2_0, y = pred)) +
        geom_ribbon(aes(ymin = pred - se.fit, ymax = pred + se.fit), fill = "red", alpha = 0.2) +
        geom_line(color = "red", size = 1) +
        labs(title = "", x = "", y = "") +
        theme_minimal() +
        theme(plot.title = element_text(size = 8), axis.title.x = element_text(size = 8))
      # xlim(x_min, x_max) + ylim(y_min, y_max)
    } else if (i > j) {
      # ÏÂÈý½Ç: RedYeo7n_Btw ±äÁ¿µÄÇúÏß
      index <- ((j - 1) * 7) - (((j - 1) * j) / 2) + (i - j)  # ¼ÆËã±äÁ¿Ë÷Òý
      if (index <= length(btw_vars)) {
        model <- gamm(as.formula(paste(btw_vars[index], "~ s(AgeAttend_2_0) +", covariates)), data = data)
        pred_data <- data.frame(
          AgeAttend_2_0 = seq(x_min, x_max, length.out = 100),
          Sex_0_0 = mean(data$Sex_0_0, na.rm = TRUE),
          BMI_2_0 = mean(data$BMI_2_0, na.rm = TRUE),
          HeadMotion_2_0 = mean(data$HeadMotion_2_0, na.rm = TRUE),
          Race_1 = mean(data$Race_1, na.rm = TRUE),
          Race_2 = mean(data$Race_2, na.rm = TRUE),
          Race_3 = mean(data$Race_3, na.rm = TRUE),
          Race_4 = mean(data$Race_4, na.rm = TRUE),
          TDI_0_0 = mean(data$TDI_0_0, na.rm = TRUE),
          Vol_WB_TIV_2_0 = mean(data$Vol_WB_TIV_2_0, na.rm = TRUE)
        )
        pred <- predict(model$gam, newdata = pred_data, se.fit = TRUE)
        pred_data$pred <- pred$fit
        pred_data$se.fit <- pred$se.fit
        
        p <- ggplot(pred_data, aes(x = AgeAttend_2_0, y = pred)) +
          geom_ribbon(aes(ymin = pred - se.fit, ymax = pred + se.fit), fill = "red", alpha = 0.2) +
          geom_line(color = "red", size = 1) +
          labs(title = "", x = "", y = "") +
          theme_minimal() +
          theme(plot.title = element_text(size = 8), axis.title.x = element_text(size = 8))
        # xlim(x_min, x_max) + ylim(y_min, y_max)
      } else {
        p <- ggplot() + theme_void()  # ¿Õ°×Í¼
      }
    } else {
      p <- ggplot() + theme_void()  # ÉÏÈý½ÇÁô¿Õ
    }
    
    plots[[length(plots) + 1]] <- p
  }
}

# »æÖÆ 7x7 Íø¸ñ
grid.arrange(grobs = plots, nrow = 7, ncol = 7)




################################################################################
#                                Association
################################################################################

library(circlize)

A <- matrix(c(0.0279115762153883,-0.0595012686759988,0.0539812440804207,0.00836794833398643,0.0357022072861227,0.0583145444984285,0.0614022915168292,-0.0521604963898887,0.0135408703637860,0.0162042208857329,-0.0173797568463735,
              -0.00337346360183192,-0.0118209378216270,-0.0249164371104383,0.0151367217015674,-0.0377048433350643,-0.0182063586084313,-0.0132828147635505,0.00356356583354949,-0.0358209689735020,0.00536342465528599,0.0229504026679898), nrow = 11)
B <- matrix(c(7.11946328161013e-06,9.31567747578803e-22,3.50780208585606e-18,0.178075450780731,9.05708681260064e-09,5.84591184525826e-21,4.55448853346531e-23,4.46195287635335e-17,0.0293111921139243,0.00910813113070343,0.00515489842818440,
              0.587419781205330,0.0571126106312354,6.05890354992155e-05,0.0148449027007438,1.27814705962837e-09,0.00338688646830913,0.0325378190842783,0.566307599831030,8.08702758312701e-09,0.388049840401962,0.000220781069739105), nrow = 11)  # ÏÔÖøÐÔ£º1´ú±íÏÔÖø

labels <- c('WMHyperInt_WB_Total_2_0','FA_WB_Total_2_0','MD_WB_Total_2_0','MO_WB_Total_2_0',
            'L1_WB_Total_2_0','L2_WB_Total_2_0','L3_WB_Total_2_0','ICVF_WB_Total_2_0','ISOVF_WB_Total_2_0',
            'Area_WB_Total_2_0','Thickness_WB_Total_2_0')

rownames(A) <- labels
rownames(B) <- labels

colnames(A) <- c("Synergy","Redundancy")
colnames(B) <- c("Synergy","Redundancy")

library(circlize)
library(ComplexHeatmap)

rownames(A) <- c("WMHyperInt", "FA", "MD", "MO", "L1", "L2", "L3", "ICVF", "ISOVF", "Area", "Thickness")
colnames(A) <- c("Synergy", "Redundancy")

# ·Ö×é£¨±ÈÈç°´½á¹¹·Ö²ã£©
split_var <- factor(c("WM", "WM", "WM", "WM", "WM", "WM", "WM", "WM", "WM", "GM", "GM"))

# ÑÕÉ«Ó³Éä
col_fun <- circlize::colorRamp2(c(-0.06, 0, 0.06), c("blue", "white", "red"))
circos.clear()
circos.heatmap(A,
               split = split_var,
               col = col_fun,
               cluster = FALSE,
               rownames.side = "outside",   # ÉèÖÃÐÐ±êÇ©ÏÔÊ¾ÔÚÄÚ²à
               track.height = 0.2)








# ¼ÓÔØËùÐè°ü
library(lme4)
library(ggplot2)

variables_of_interest <- c("eid","Fluid_intelligence_score_2_0", "RedYeo7n_Total", "SynYeo7n_Total", 
                           "Sex_0_0","AgeAttend_2_0","BMI_2_0", "HeadMotion_2_0", "Race_1", "Race_2", 
                           "Race_3", "Race_4","TDI_0_0", "Vol_WB_TIV_2_0","Centre_0_0")
data_tmp <- na.omit(data[,variables_of_interest])

m <- lmer(Fluid_intelligence_score_2_0 ~ Sex_0_0 + RedYeo7n_Total + AgeAttend_2_0 + 
            AgeAttend_2_0 * RedYeo7n_Total + BMI_2_0 + HeadMotion_2_0 + 
            Race_1 + Race_2 + Race_3 + Race_4 + TDI_0_0 + Vol_WB_TIV_2_0 + 
            (1 | Centre_0_0), data = data_tmp)

# ÌáÈ¡ÏµÊý
summary(m)$coefficients

# ÄâºÏ²Ð²îÄ£ÐÍ
residualModel <- lmer(RedYeo7n_Total ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + 
                        Race_1 + Race_2 + Race_3 + Race_4 + TDI_0_0 + 
                        Vol_WB_TIV_2_0 + (1 | Centre_0_0), data = data_tmp)

# ÌáÈ¡²Ð²î
data_tmp$RedYeo7n_Residual <- residuals(residualModel)

# ÄâºÏ½»»¥Ð§Ó¦Ä£ÐÍ
interactionModel <- lm(Fluid_intelligence_score_2_0 ~ RedYeo7n_Residual + AgeAttend_2_0 + 
                         AgeAttend_2_0 * RedYeo7n_Residual, data = data_tmp)

# ²é¿´½»»¥Ð§Ó¦Ä£ÐÍÏµÊý
summary(interactionModel)

# ÌáÈ¡×Ô±äÁ¿¡¢µ÷½Ú±äÁ¿ºÍÒò±äÁ¿
X <- data_tmp$RedYeo7n_Residual  # ×Ô±äÁ¿
Z <- data_tmp$AgeAttend_2_0      # µ÷½Ú±äÁ¿
Y <- data_tmp$Fluid_intelligence_score_2_0  # Òò±äÁ¿

# ¼ÆËãµ÷½Ú±äÁ¿µÄ¾ùÖµºÍ±ê×¼²î
Z_mean <- mean(Z)
Z_std <- sd(Z)

# ¶¨Òåµ÷½Ú±äÁ¿µÄÈý¸öË®Æ½
Z_low <- Z_mean - Z_std
Z_mid <- Z_mean
Z_high <- Z_mean + Z_std

# ´´½¨×Ô±äÁ¿µÄÈ¡Öµ·¶Î§
X_range <- seq(min(X), max(X), length.out = 100)

# ¼ÆËãÃ¿¸öµ÷½Ú±äÁ¿Ë®Æ½ÏÂµÄÔ¤²âÖµ
pred_low <- predict(interactionModel, newdata = data.frame(RedYeo7n_Residual = X_range, 
                                                           AgeAttend_2_0 = Z_low), 
                    interval = "confidence")

pred_mid <- predict(interactionModel, newdata = data.frame(RedYeo7n_Residual = X_range, 
                                                           AgeAttend_2_0 = Z_mid), 
                    interval = "confidence")

pred_high <- predict(interactionModel, newdata = data.frame(RedYeo7n_Residual = X_range, 
                                                            AgeAttend_2_0 = Z_high), 
                     interval = "confidence")

# ½«Ô¤²âÖµºÍÖÃÐÅÇø¼ä×ª»»ÎªÊý¾Ý¿ò
pred_data <- data.frame(
  X_range = rep(X_range, 3),
  Predicted = c(pred_low[, 1], pred_mid[, 1], pred_high[, 1]),
  Lower = c(pred_low[, 2], pred_mid[, 2], pred_high[, 2]),
  Upper = c(pred_low[, 3], pred_mid[, 3], pred_high[, 3]),
  AgeAttend = rep(c("Low", "Mid", "High"), each = length(X_range))
)

ggplot(pred_data, aes(x = X_range, y = Predicted, color = AgeAttend)) +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = AgeAttend, color = NA), alpha = 0.15) +
  labs(x = "Redundancy Residual", y = "Fluid Intelligence Score", 
       title = "Interaction Effect of Redundancy and Age") +
  scale_fill_manual(values = c("blue", "green", "red")) +
  scale_color_manual(values = c("blue", "green", "red"),  # ÉèÖÃ²»Í¬ÄêÁä¶ÎµÄÑÕÉ«
                     labels = c(sprintf('AgeAttend = %.2f (Low)', Z_low), 
                                sprintf('AgeAttend = %.2f (Mid)', Z_mid),
                                sprintf('AgeAttend = %.2f (High)', Z_high))) +  # ×Ô¶¨ÒåÍ¼Àý±êÇ©
  theme_minimal()




# ¼ÙÉèdata_tmpÒÑ¾­¼ÓÔØµ½R»·¾³ÖÐ

# ÏßÐÔ»ìºÏÐ§Ó¦Ä£ÐÍ
m <- lmer(Fluid_intelligence_score_2_0 ~ Sex_0_0 + SynYeo7n_Total + AgeAttend_2_0 + 
            AgeAttend_2_0 * SynYeo7n_Total + BMI_2_0 + HeadMotion_2_0 + 
            Race_1 + Race_2 + Race_3 + Race_4 + TDI_0_0 + Vol_WB_TIV_2_0 + 
            (1 | Centre_0_0), data = data_tmp)

# ÌáÈ¡ÏµÊý
summary(m)$coefficients

# ÄâºÏ²Ð²îÄ£ÐÍ
residualModel <- lmer(SynYeo7n_Total ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + 
                        Race_1 + Race_2 + Race_3 + Race_4 + TDI_0_0 + 
                        Vol_WB_TIV_2_0 + (1 | Centre_0_0), data = data_tmp)

# ÌáÈ¡²Ð²î
data_tmp$SynYeo7n_Residual <- residuals(residualModel)

# ÄâºÏ½»»¥Ð§Ó¦Ä£ÐÍ
interactionModel <- lm(Fluid_intelligence_score_2_0 ~ SynYeo7n_Residual + AgeAttend_2_0 + 
                         AgeAttend_2_0 * SynYeo7n_Residual, data = data_tmp)

# ²é¿´½»»¥Ð§Ó¦Ä£ÐÍÏµÊý
summary(interactionModel)

# ÌáÈ¡×Ô±äÁ¿¡¢µ÷½Ú±äÁ¿ºÍÒò±äÁ¿
X <- data_tmp$SynYeo7n_Residual  # ×Ô±äÁ¿
Z <- data_tmp$AgeAttend_2_0      # µ÷½Ú±äÁ¿
Y <- data_tmp$Fluid_intelligence_score_2_0  # Òò±äÁ¿

# ¼ÆËãµ÷½Ú±äÁ¿µÄ¾ùÖµºÍ±ê×¼²î
Z_mean <- mean(Z)
Z_std <- sd(Z)

# ¶¨Òåµ÷½Ú±äÁ¿µÄÈý¸öË®Æ½
Z_low <- Z_mean - Z_std
Z_mid <- Z_mean
Z_high <- Z_mean + Z_std

# ´´½¨×Ô±äÁ¿µÄÈ¡Öµ·¶Î§
X_range <- seq(min(X), max(X), length.out = 100)

# ¼ÆËãÃ¿¸öµ÷½Ú±äÁ¿Ë®Æ½ÏÂµÄÔ¤²âÖµ
pred_low <- predict(interactionModel, newdata = data.frame(SynYeo7n_Residual = X_range, 
                                                           AgeAttend_2_0 = Z_low), 
                    interval = "confidence")

pred_mid <- predict(interactionModel, newdata = data.frame(SynYeo7n_Residual = X_range, 
                                                           AgeAttend_2_0 = Z_mid), 
                    interval = "confidence")

pred_high <- predict(interactionModel, newdata = data.frame(SynYeo7n_Residual = X_range, 
                                                            AgeAttend_2_0 = Z_high), 
                     interval = "confidence")

# ½«Ô¤²âÖµºÍÖÃÐÅÇø¼ä×ª»»ÎªÊý¾Ý¿ò
pred_data <- data.frame(
  X_range = rep(X_range, 3),
  Predicted = c(pred_low[, 1], pred_mid[, 1], pred_high[, 1]),
  Lower = c(pred_low[, 2], pred_mid[, 2], pred_high[, 2]),
  Upper = c(pred_low[, 3], pred_mid[, 3], pred_high[, 3]),
  AgeAttend = rep(c("Low", "Mid", "High"), each = length(X_range))
)

ggplot(pred_data, aes(x = X_range, y = Predicted, color = AgeAttend)) +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = AgeAttend, color = NA), alpha = 0.15) +
  labs(x = "Synergy Residual", y = "Fluid Intelligence Score", 
       title = "Interaction Effect of Synergy and Age") +
  scale_fill_manual(values = c("blue", "green", "red")) +
  scale_color_manual(values = c("blue", "green", "red"),  # ÉèÖÃ²»Í¬ÄêÁä¶ÎµÄÑÕÉ«
                     labels = c(sprintf('Age = %.2f (Low)', Z_low), 
                                sprintf('Age = %.2f (Mid)', Z_mid),
                                sprintf('Age = %.2f (High)', Z_high))) +  # ×Ô¶¨ÒåÍ¼Àý±êÇ©
  theme_minimal()






################################################################################
#                             Interaction Curves
################################################################################









# Step 1: ·Ö×é
UKB_Fluid_intelligence = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_Fluid_intelligence.txt",na = "NaN")

data2 <- data.frame(eid = data$eid,AgeAttend_2_0 = data$AgeAttend_2_0)

for (i in 1:7) {
  col_idx <- get_col_indices(i)
  between_names <- paste0("RedYeo7n_Btw", col_idx)
  vals <- cbind(data[[paste0("RedYeo7n_In", i)]], data[, between_names])
  data2[[network_names[i]]] <- rowMeans(vals, na.rm = TRUE)
}
data2[["Total"]] <- data[, "RedYeo7n_Total"]


data2 <- merge(data2,UKB_Fluid_intelligence[,c(1,2)])
data2 <- na.omit(data2)
data2$FI_group <- cut(data2$Fluid_intelligence_score_2_0,
                      breaks = quantile(data2$Fluid_intelligence_score_2_0, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      labels = c("Low", "Mid", "High"),
                      include.lowest = TRUE)


network_names <- c("Visual", "Somatomotor", "DorsalAttention", "VentralAttention",
                   "Limbic", "Frontoparietal", "Default")

# Step 2: ¿ÉÊÓ»¯ÇúÏß
library(mgcv)
library(ggplot2)

ggplot(data2, aes(x = AgeAttend_2_0, y = Total)) +
  geom_smooth(aes(color = FI_group, fill = FI_group), method = "gam",
              formula = y ~ s(x), se = TRUE, alpha = 0.2, size = 0) +  # Ö»»­ÒõÓ°
  geom_smooth(aes(color = FI_group), method = "gam",
              formula = y ~ s(x), se = FALSE, size = 1.2,
              alpha = c(Low = 1, Mid = 0.6, High = 0.3)[data2$FI_group]) +  # ¿ØÖÆÏßÌõÍ¸Ã÷¶È
  scale_color_manual(values = rep("#08306B", 3)) +
  scale_fill_manual(values = rep("#08306B", 3)) +
  labs(title = "Separate alpha for line and shading", x = "Age", y = "Total") +
  theme_bw()

model1 <- gam(Total ~ s(AgeAttend_2_0, by = Fluid_intelligence_score_2_0) + Fluid_intelligence_score_2_0,
              data = data2)
summary(model1)

# Ä£ÐÍÎÞ½»»¥Ïî
model_null <- gam(Total ~ s(AgeAttend_2_0) + s(Fluid_intelligence_score_2_0),
                  data = data2)

# Ä£ÐÍÓÐ½»»¥Ïî
model_full <- gam(Total ~ ti(AgeAttend_2_0, Fluid_intelligence_score_2_0) + 
                    s(AgeAttend_2_0) + 
                    s(Fluid_intelligence_score_2_0),
                  data = data2)

anova(model_null, model_full, test = "Chisq")




























# Step 1: ·Ö×é
UKB_Fluid_intelligence = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_Fluid_intelligence.txt",na = "NaN")

data2 <- data.frame(eid = data$eid,AgeAttend_2_0 = data$AgeAttend_2_0)

for (i in 1:7) {
  col_idx <- get_col_indices(i)
  between_names <- paste0("RedYeo7n_Btw", col_idx)
  vals <- cbind(data[[paste0("RedYeo7n_In", i)]], data[, between_names])
  data2[[network_names[i]]] <- rowMeans(vals, na.rm = TRUE)
}
data2[["Global"]] <- data[, "RedYeo7n_Total"]


data2 <- merge(data2,UKB_Fluid_intelligence[,c(1,2)])
data2 <- na.omit(data2)
data2$FI_group <- cut(data2$Fluid_intelligence_score_2_0,
                      breaks = quantile(data2$Fluid_intelligence_score_2_0, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      labels = c("Low", "Mid", "High"),
                      include.lowest = TRUE)


library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)

# ¼ÙÉèÒÑÓÐµÄÊý¾ÝÊÇ data2£¬°üº¬ÒÔÏÂ±äÁ¿£º
# eid, AgeAttend_2_0, Fluid_intelligence_score_2_0, Visual, ..., Default, Total
# ²¢ÒÑ¾­´´½¨ºÃ data2$FI_group£¨·Ö×é±äÁ¿£©

# ¶¨ÒåÍøÂçÃû³ÆË³Ðò£¨°üÀ¨ Global£©
network_levels <- c("Visual", "Somatomotor", "DorsalAttention", "VentralAttention", 
                    "Limbic", "Frontoparietal", "Default", "Global")

# ¶ÔÓ¦ÑÕÉ«£¨Äã×Ô¶¨ÒåµÄ£©
Yeo7Color <- c(
  rgb(0.471,0.0710,0.522),
  rgb(0.275,0.510,0.706),
  rgb(0,0.463,0.0550),
  rgb(0.769,0.224,0.976),
  rgb(0.863,0.973,0.639),
  rgb(0.902,0.576,0.129),
  rgb(0.804,0.239,0.306),
  "blue"
)
names(Yeo7Color) <- network_levels

# ¶¨Òå FI ·Ö×éºÍÍ¸Ã÷¶È
FI_levels <- c("Low", "Mid", "High")
FI_alpha_levels <- c("High" = 1, "Mid" = 0.7, "Low" = 0.4)

data2$FI_group <- cut(data2$Fluid_intelligence_score_2_0,
                      breaks = quantile(data2$Fluid_intelligence_score_2_0, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      labels = c("Low", "Mid", "High"),
                      include.lowest = TRUE)


# ×ªÎª³¤¸ñÊ½
data_sub <- data2[, c("AgeAttend_2_0", "FI_group", network_levels)]

data_long <- reshape(
  data_sub,
  varying = network_levels,
  v.names = "Connectivity",
  timevar = "Network",
  times = network_levels,
  direction = "long"
)

data_long$Network <- factor(data_long$Network, levels = network_levels)
data_long$FI_group <- factor(data_long$FI_group, levels = FI_levels)
data_long$Net_FI <- interaction(data_long$Network, data_long$FI_group, sep = "_")

# ¹¹ÔìÑÕÉ«±í£ºÃ¿¸öÍøÂç¹Ì¶¨ÑÕÉ«£¬Í¨¹ýÍ¸Ã÷¶È±ä»¯³öÈýÖÖ
make_alpha <- function(hex, alpha) {
  rgb <- col2rgb(hex) / 255
  rgb(rgb[1], rgb[2], rgb[3], alpha = alpha)
}

Net_FI_levels <- as.vector(outer(network_levels, FI_levels, paste, sep = "_"))

Net_FI_colors <- setNames(
  mapply(function(net, fi) make_alpha(Yeo7Color[net], FI_alpha_levels[fi]),
         rep(network_levels, times = 3),
         rep(FI_levels, each = 8),
         SIMPLIFY = TRUE),
  Net_FI_levels
)

# ÉèÖÃ Net_FI Îª factor£¬±£³ÖË³Ðò
data_long$Net_FI <- factor(data_long$Net_FI, levels = Net_FI_levels)

# »æÍ¼
ggplot(data_long, aes(x = AgeAttend_2_0, y = Connectivity)) +
  geom_smooth(aes(group = Net_FI, color = Net_FI, linetype = FI_group),
              method = "gam", formula = y ~ s(x),
              se = TRUE, size = 1.2, alpha = 0.05) +
  scale_color_manual(values = Net_FI_colors) +
  scale_linetype_manual(values = c("High" = "solid", "Mid" = "11", "Low" = "22")) + 
  scale_y_continuous(labels = function(x) sprintf("%.4f", x)) +
  facet_wrap(~ Network, nrow = 1, scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +  # ÏßÐÍºÍÑÕÉ«¶¼ÒªÏÔÊ¾£¬·ÅÓÒ²à
  labs(title = "Connectivity vs Age by Network and FI Group",
       x = "Age at Assessment", y = "Connectivity", linetype = "FI Group")






data_plot <- subset(data_long, FI_group != "Mid")

ggplot(data_plot, aes(x = AgeAttend_2_0, y = Connectivity)) +
  geom_smooth(aes(group = interaction(Net_FI, FI_group), 
                  color = Net_FI, linetype = FI_group),
              method = "gam", formula = y ~ s(x),
              se = TRUE, size = 1.2, alpha = 0.05) +
  scale_color_manual(values = Net_FI_colors) +
  scale_linetype_manual(values = c("High" = "solid", "Low" = "22")) + 
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
  facet_wrap(~ Network, nrow = 1, scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  labs(title = "Connectivity vs Age by Network and FI Group",
       x = "Age at Assessment", y = "Connectivity",
       linetype = "FI Group", color = "Net FI")












# Step 1: ·Ö×é
UKB_Fluid_intelligence = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_Fluid_intelligence.txt",na = "NaN")

Y <- "Syn"
Yc <- "red"
data2 <- data.frame(eid = data$eid,AgeAttend_2_0 = data$AgeAttend_2_0)
for (i in 1:7) {
  col_idx <- get_col_indices(i)
  between_names <- paste0(Y,"Yeo7n_Btw", col_idx)
  vals <- cbind(data[[paste0(Y,"Yeo7n_In", i)]], data[, between_names])
  data2[[network_names[i]]] <- rowMeans(vals, na.rm = TRUE)
}
data2[["Global"]] <- data[, paste0(Y,"Yeo7n_Total")]

data2 <- merge(data2,UKB_Fluid_intelligence[,c(1,2)])
data2 <- na.omit(data2)
data2$FI_group <- cut(data2$Fluid_intelligence_score_2_0,
                      breaks = quantile(data2$Fluid_intelligence_score_2_0, probs = c(0, 1/2,1), na.rm = TRUE),
                      labels = c("Low", "High"),
                      include.lowest = TRUE)

library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)

network_levels <- c("Visual", "Somatomotor", "DorsalAttention", "VentralAttention", 
                    "Limbic", "Frontoparietal", "Default", "Global")
Yeo7Color <- c(
  rgb(0.471,0.0710,0.522),
  rgb(0.275,0.510,0.706),
  rgb(0,0.463,0.0550),
  rgb(0.769,0.224,0.976),
  rgb(0.863,0.973,0.639),
  rgb(0.902,0.576,0.129),
  rgb(0.804,0.239,0.306),
  Yc
)
names(Yeo7Color) <- network_levels

# ¶¨Òå FI ·Ö×éºÍÍ¸Ã÷¶È
FI_levels <- c("Low", "High")
FI_alpha_levels <- c("High" = 1, "Low" = 0.4)

data2$FI_group <- cut(data2$Fluid_intelligence_score_2_0,
                      breaks = quantile(data2$Fluid_intelligence_score_2_0, probs = c(0, 1/2, 1), na.rm = TRUE),
                      labels = c("Low",  "High"),
                      include.lowest = TRUE)
data_sub <- data2[, c("AgeAttend_2_0", "FI_group", network_levels)]

data_long <- reshape(
  data_sub,
  varying = network_levels,
  v.names = "Connectivity",
  timevar = "Network",
  times = network_levels,
  direction = "long"
)

data_long$Network <- factor(data_long$Network, levels = network_levels)
data_long$FI_group <- factor(data_long$FI_group, levels = FI_levels)
data_long$Net_FI <- interaction(data_long$Network, data_long$FI_group, sep = "_")

# ¹¹ÔìÑÕÉ«±í£ºÃ¿¸öÍøÂç¹Ì¶¨ÑÕÉ«£¬Í¨¹ýÍ¸Ã÷¶È±ä»¯³öÈýÖÖ
make_alpha <- function(hex, alpha) {
  rgb <- col2rgb(hex) / 255
  rgb(rgb[1], rgb[2], rgb[3], alpha = alpha)
}

Net_FI_levels <- as.vector(outer(network_levels, FI_levels, paste, sep = "_"))

Net_FI_colors <- setNames(
  mapply(function(net, fi) make_alpha(Yeo7Color[net], FI_alpha_levels[fi]),
         rep(network_levels, times = 2),
         rep(FI_levels, each = 8),
         SIMPLIFY = TRUE),
  Net_FI_levels
)

# ÉèÖÃ Net_FI Îª factor£¬±£³ÖË³Ðò
data_long$Net_FI <- factor(data_long$Net_FI, levels = Net_FI_levels)

# »æÍ¼
ggplot(data_long, aes(x = AgeAttend_2_0, y = Connectivity)) +
  geom_smooth(aes(group = Net_FI, color = Net_FI, linetype = FI_group),
              method = "gam", formula = y ~ s(x),
              se = TRUE, size = 1.2, alpha = 0.05) +
  scale_color_manual(values = Net_FI_colors) +
  scale_linetype_manual(values = c("High" = "solid",  "Low" = "22")) + 
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
  facet_wrap(~ Network, nrow = 1, scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +  # ÏßÐÍºÍÑÕÉ«¶¼ÒªÏÔÊ¾£¬·ÅÓÒ²à
  labs(title = "Connectivity vs Age by Network and FI Group",
       x = "Age at Assessment", y = "Connectivity", linetype = "FI Group")


unique(data_long$Network)
data_subset <- subset(data_long, Network == "Visual")
m1 <- gam(Connectivity ~ s(AgeAttend_2_0), data = data_subset)
m2 <- gam(Connectivity ~ s(AgeAttend_2_0, by = FI_group) + FI_group, 
          data = data_subset)

# ¼ì²é²îÒì£º¹ì¼£ÊÇ·ñ²»Í¬
anova(m1, m2, test = "Chisq")

# ²é¿´¾ßÌå½á¹û
summary(m2)

lm3 <- lm(Connectivity ~ AgeAttend_2_0 * FI_group, data = data_subset)

summary(lm3)


# 1500 260


ggplot(data_long, aes(x = AgeAttend_2_0, y = Connectivity)) +
  geom_smooth(aes(group = interaction(Net_FI, FI_group), 
                  color = Net_FI, linetype = FI_group),
              method = "gam", formula = y ~ s(x),
              se = TRUE, size = 1.2, alpha = 0.05) +
  scale_color_manual(values = Net_FI_colors) +
  scale_linetype_manual(values = c("High" = "solid", "Low" = "22")) + 
  scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
  facet_wrap(~ Network, nrow = 1, scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right") +   # ? ÏÔÊ¾ÔÚÓÒ²à
  labs(title = "Connectivity vs Age by Network and FI Group",
       x = "Age at Assessment", y = "Connectivity", 
       linetype = "FI Group", color = "Net FI")  # ? ¸øÍ¼Àý¼ÓÃû×Ö














################################################################################
#                             Moderation Analysis
################################################################################






# 1. Õë¶ÔÃ¿¸ö Network ÄâºÏ´ø½»»¥ÏîµÄÄ£ÐÍ
model_results <- data_long %>%
  group_by(Network) %>%
  group_split() %>%
  map(~ {
    gam_full <- gam(Connectivity ~ s(AgeAttend_2_0, by = FI_group) + FI_group,
                    data = ., method = "REML")
    gam_null <- gam(Connectivity ~ s(AgeAttend_2_0), 
                    data = ., method = "REML")
    anova_result <- anova(gam_null, gam_full, test = "Chisq")
    
    tibble(
      Network = unique(.$Network),
      p_value = anova_result$`Pr(>Chi)`[2]
    )
  }) %>%
  bind_rows()

# 2. FDR Ð£Õý£¨¿ÉÑ¡£©
model_results <- model_results %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

print(model_results)







library(lme4)
library(lmerTest)  # ¿ÉÌá¹© p Öµ

# ÄâºÏÄ£ÐÍ£º°üÀ¨ Age * FI_group ½»»¥Ïî
modelx <- lmerTest::lmer(RedYeo7n_Total ~ AgeAttend_2_0 * Income_0_0 +
                Sex_0_0 + BMI_2_0 + HeadMotion_2_0 +
                Race_1 + Race_2 + Race_3 + Race_4 +
                TDI_0_0 + Vol_WB_TIV_2_0 +
                (1 | Centre_0_0),
              data = data)

# Êä³ö½á¹û
summary(modelx)


library(nlme)
modelx <- lme(Redundancy_VN ~ AgeAttend_2_0 * Income_0_0 +
                Sex_0_0 + BMI_2_0 + HeadMotion_2_0 +
                Race_1 + Race_2 + Race_3 + Race_4 +
                TDI_0_0 + Vol_WB_TIV_2_0,
              random = ~ 1 | Centre_0_0,
              data = data,
              na.action = na.omit)

summary(modelx)





library(nlme)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)



UKB_PID_Yeo7n = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Yeo7n.txt",na = "NaN")

# Step 0: »ñÈ¡±äÁ¿Ãû
Y_vars <- names(UKB_RedSyn_2_0)[2:59]      # ³ýÈ¥µÚÒ»ÁÐ£¨¼ÙÉèÊÇ ID£©

UKB_XVars <- UKB_XVars_0_0
names(UKB_XVars)[2:ncol(UKB_XVars_0_0)] <- paste0("x",gsub("-0.0","",names(UKB_XVars_0_0)[2:ncol(UKB_XVars_0_0)]))
X_vars <- names(UKB_XVars)[-1]      # ³ýÈ¥µÚÒ»ÁÐ£¨¼ÙÉèÊÇ ID£©

data_merge <- merge(data,UKB_XVars)

# Step 2: ´æ´¢½á¹û
results_all <- list()

# Step 3: Ë«²ãÑ­»·
counter <- 0
total <- length(Y_vars) * length(X_vars)

for (y_var in Y_vars) {
  for (x_var in X_vars) {
    
    counter <- counter + 1
    message("Running model ", counter, " of ", total, ": Y = ", y_var, ", X = ", x_var)
    
    formula_str <- paste0(
      y_var, " ~ AgeAttend_2_0 * ", x_var, " + ",
      "Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + ",
      "Race_1 + Race_2 + Race_3 + Race_4 + ",
      "Vol_WB_TIV_2_0"
    )
    
    model_formula <- as.formula(formula_str)
    
    model <- lme(fixed = model_formula,
          random = ~1 | Centre_0_0,
          data = data_merge,
          na.action = na.omit)
    
    # »ñÈ¡¹Ì¶¨Ð§Ó¦µÄÏµÊý±í£¨°üº¬ Estimate, Std.Error, DF, t-value, p-value£©
    coef_table <- summary(model)$tTable
    
    # ÕÒµ½½»»¥Ïî£¨×¢ÒâÆ¥ÅäÃû×Ö¿ÉÄÜÊÇ AgeAttend_2_0:X_var »ò X_var:AgeAttend_2_0£©
    interaction_name <- paste0("AgeAttend_2_0:", x_var)
    alt_name <- paste0(x_var, ":AgeAttend_2_0")
    
    # ¼ì²é½»»¥ÏîÔÚÄÄÒ»ÐÐ
    row_index <- which(rownames(coef_table) %in% c(interaction_name, alt_name))
    
    # Èç¹ûÕÒµ½ÁË£¬¾ÍÌáÈ¡ÐÅÏ¢
    if (length(row_index) > 0) {
      interaction_row <- coef_table[row_index, ]
      
      results_all[[paste(y_var, x_var, sep = "_")]] <- data.frame(
        Y = y_var,
        X = x_var,
        Estimate = interaction_row["Value"],
        Std_Error = interaction_row["Std.Error"],
        t_value = interaction_row["t-value"],
        p_value = interaction_row["p-value"]
      )
      
    } else {
      cat("?? ½»»¥Ïî", interaction_name, "Î´ÕÒµ½¡£\n")
    }
    
  }
}



# Step 4: »ã×Ü½á¹û
results_interact <- bind_rows(results_all)

# Step 5: FDRÐ£Õý
results_interact <- results_interact %>%
  mutate(p_fdr = p.adjust(p_value, method = "fdr"))

# Step 6: ÅÅÐòÏÔÊ¾Ç°¼¸¸öÏÔÖøµ÷½Ú×÷ÓÃ
top_results <- results_interact %>%
  arrange(p_value) %>%
  filter(p_value < 0.05/58)

print(top_results)

# Step 7: ¿ÉÊÓ»¯ÈÈÍ¼£¨Y ¡Á X£©
heatmap_data <- results_interact %>%
  mutate(Signif = ifelse(p_value < 0.05/14, "*", "")) %>%
  select(Y, X, Estimate, t_value, p_value, Signif)

ggplot(heatmap_data, aes(x = X, y = Y, fill = t_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Signif), size = 5, color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal(base_size = 12) +
  labs(title = "Moderation Effect",
       x = "Moderator X", y = "Brain Variable Y", fill = "T Value of Interaction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Step 0: »ñÈ¡±äÁ¿Ãû
Y_vars <- names(UKB_PID_Yeo7n)[95:108]      # ³ýÈ¥µÚÒ»ÁÐ£¨¼ÙÉèÊÇ ID£©

UKB_XVars <- UKB_XVars_0_0
names(UKB_XVars)[2:ncol(UKB_XVars_0_0)] <- paste0("x",gsub("-0.0","",names(UKB_XVars_0_0)[2:ncol(UKB_XVars_0_0)]))
X_vars <- names(UKB_XVars)[-1]      # ³ýÈ¥µÚÒ»ÁÐ£¨¼ÙÉèÊÇ ID£©

data_merge <- merge(data,UKB_XVars)

# Step 2: ´æ´¢½á¹û
results_all <- list()

# Step 3: Ë«²ãÑ­»·
counter <- 0
total <- length(Y_vars) * length(X_vars)

for (y_var in Y_vars) {
  for (x_var in X_vars) {
    
    counter <- counter + 1
    message("Running model ", counter, " of ", total, ": Y = ", y_var, ", X = ", x_var)
    
    formula_str <- paste0(
      y_var, " ~ AgeAttend_2_0 + ", x_var, " + ",
      "Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + ",
      "Race_1 + Race_2 + Race_3 + Race_4 + ",
      "TDI_0_0 + Vol_WB_TIV_2_0"
    )
    
    model_formula <- as.formula(formula_str)
    
    model <- lme(fixed = model_formula,
                 random = ~1 | Centre_0_0,
                 data = data_merge,
                 na.action = na.omit)
    
    # »ñÈ¡¹Ì¶¨Ð§Ó¦µÄÏµÊý±í£¨°üº¬ Estimate, Std.Error, DF, t-value, p-value£©
    coef_table <- summary(model)$tTable
    
    # ¼ì²é½»»¥ÏîÔÚÄÄÒ»ÐÐ
    row_index <- which(rownames(coef_table) %in% x_var)
    
    # Èç¹ûÕÒµ½ÁË£¬¾ÍÌáÈ¡ÐÅÏ¢
    if (length(row_index) > 0) {
      main_row <- coef_table[row_index, ]
      
      results_all[[paste(y_var, x_var, sep = "_")]] <- data.frame(
        Y = y_var,
        X = x_var,
        Estimate = main_row["Value"],
        Std_Error = main_row["Std.Error"],
        t_value = main_row["t-value"],
        p_value = main_row["p-value"]
      )
      
    } else {
      cat("?? ½»»¥Ïî", x_var, "Î´ÕÒµ½¡£\n")
    }
    
  }
}



# Step 4: »ã×Ü½á¹û
results_main <- bind_rows(results_all)

# Step 6: ÅÅÐòÏÔÊ¾Ç°¼¸¸öÏÔÖøµ÷½Ú×÷ÓÃ
top_results <- results_main %>%
  arrange(p_value) %>%
  filter(p_value < 0.05/14)

print(top_results)

# Step 7: ¿ÉÊÓ»¯ÈÈÍ¼£¨Y ¡Á X£©
heatmap_data <- results_main %>%
  mutate(Signif = ifelse(p_value < 0.05/14, "*", "")) %>%
  select(Y, X, Estimate, t_value, p_value, Signif)

ggplot(heatmap_data, aes(x = X, y = Y, fill = t_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Signif), size = 5, color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal(base_size = 12) +
  labs(title = "Moderation Effect",
       x = "Moderator X", y = "Brain Variable Y", fill = "T Value of Interaction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





################################################################################
#                        Aging Trajectory (Polynomial)
################################################################################

library(readr)
UKB_RedSyn_2_0 <- read_csv("ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/result/ukb_RedSynYeo7n_2_0.csv")
UKB_BasicMRI <- read_csv("ZJLab/UKBiobank_Project/data/table/UKB_BrainMRICov.txt")
UKB_redsyn_tmp <- merge(UKB_BasicMRI,UKB_RedSyn_2_0)
names(UKB_redsyn_tmp)[95:152] <- paste0("Y",1:58)

values <- 8:28
n <- 7

vec2mat <- function(values,n){
  mat <- matrix(NA, nrow = n, ncol = n)
  if (length(values) > sum(lower.tri(mat))) {
    stop("Values exceed the number of elements in the lower triangle.")
  }
  mat[lower.tri(mat)] <- values
  return(mat)
}


values <- 1:28
n <- 7

vec2mat.diag <- function(values,diag_values,n){
  mat <- matrix(NA, nrow = n, ncol = n)
  if (length(values) > sum(lower.tri(mat))) {
    stop("Values exceed the number of elements in the lower triangle.")
  }
  mat[lower.tri(mat)] <- values
  if (length(diag_values) != n) {
    stop("The length of diag_values must equal the size of the matrix.")
  }
  diag(mat) <- diag_values
  return(mat)
}

vec2mat.diag(8:28,1:7,7)






library(lme4)
library(tidyr)
library(dplyr)

library(lme4)

fit_models <- function(data, Y_value) {
  # È·±£ Y_value ÊÇÁÐÃû
  Y_value <- as.name(Y_value)
  
  # ¶¨Òå¹«Ê½
  formula_base <- paste(deparse(Y_value), "~ (1|Centre_0_0) + Sex_0_0 + HeadMotion_2_0 + BMI_2_0 + ",
                        "BloodPress_2 + EchoTime_2_0 + IntensityScaling_2_0 + Race_1 + Race_2 + Race_3 + Race_4 +",
                        "Vol_WB_TIV_2_0+ Handedness_0_0 + TDI_0_0")
  
  # ¹¹½¨¸÷½×Ä£ÐÍ
  model_lin <- lmer(as.formula(paste(formula_base, "+ AgeAttend_2_0")), data = data)
  model_quad <- lmer(as.formula(paste(formula_base, "+ I(AgeAttend_2_0^2) + AgeAttend_2_0")), data = data)
  model_cubic <- lmer(as.formula(paste(formula_base, "+ I(AgeAttend_2_0^3) + I(AgeAttend_2_0^2) + AgeAttend_2_0")), data = data)
  model_quartic <- lmer(as.formula(paste(formula_base, "+ I(AgeAttend_2_0^4) + I(AgeAttend_2_0^3) + I(AgeAttend_2_0^2) + AgeAttend_2_0")), data = data)
  
  # ¼ÆËã AIC Öµ
  AIC_values <- data.frame(
    Model = c("Linear", "Quadratic", "Cubic", "Quartic"),
    AIC = c(AIC(model_lin), AIC(model_quad), AIC(model_cubic), AIC(model_quartic))
  )
  
  # Ñ¡Ôñ AIC ×îÓÅÄ£ÐÍ
  best_model_name <- AIC_values$Model[which.min(AIC_values$AIC)]
  best_model <- switch(
    best_model_name,
    Linear = model_lin,
    Quadratic = model_quad,
    Cubic = model_cubic,
    Quartic = model_quartic
  )
  
  # Èç¹û²»ÊÇÏßÐÔÄ£ÐÍ£¬¼ÆËã¹Õµã
  inflection_points <- NULL
  if (best_model_name != "Linear") {
    coefficients <- fixef(best_model)
    if (best_model_name == "Quadratic") {
      # ÌáÈ¡¶þ´ÎÏîÏµÊý
      a <- coefficients["I(AgeAttend_2_0^2)"]
      b <- coefficients["AgeAttend_2_0"]
      # ¼ÆËã¶þ´ÎÄ£ÐÍµÄ¹Õµã
      inflection_points <- -b / (2 * a)
    } else if (best_model_name == "Cubic" || best_model_name == "Quartic") {
      # ÌáÈ¡Èý´Î¼°ÒÔÉÏÏîÏµÊý
      poly_derivative <- function(x) {
        a <- ifelse("I(AgeAttend_2_0^4)" %in% names(coefficients), coefficients["I(AgeAttend_2_0^4)"], 0)
        b <- ifelse("I(AgeAttend_2_0^3)" %in% names(coefficients), coefficients["I(AgeAttend_2_0^3)"], 0)
        c <- ifelse("I(AgeAttend_2_0^2)" %in% names(coefficients), coefficients["I(AgeAttend_2_0^2)"], 0)
        d <- coefficients["AgeAttend_2_0"]
        4 * a * x^3 + 3 * b * x^2 + 2 * c * x + d
      }
      # Ñ°ÕÒµ¼ÊýÎªÁãµÄµã£¨¹Õµã£©
      inflection_points <- uniroot(poly_derivative, range(data$AgeAttend_2_0))$root
    }
  }
  
  return(list(
    Best_Model = list(Name = best_model_name, Model = best_model),
    AIC_Values = AIC_values,
    Inflection_Points = inflection_points
  ))
}

results <- list()
n = 0
for (Y_col in paste0("Y",1:58)){
  n = n + 1
  results[[n]] <- fit_models(UKB_redsyn_tmp, Y_col)
}

Inflection_Points <- list()
for (i in 1:58) {
  Inflection_Points[[i]] <- results[[i]]$Inflection_Points    
}
Inflection_Points


library(ggplot2)

plot_models_with_fits <- function(data, Y_value) {
  # È·±£ Y_value ÊÇÁÐÃû
  Y_value <- as.name(Y_value)
  
  # ¹¹½¨ÄâºÏÄ£ÐÍ
  model_lin <- lm(as.formula(paste(deparse(Y_value), "~ AgeAttend_2_0")), data = data)
  model_quad <- lm(as.formula(paste(deparse(Y_value), "~ poly(AgeAttend_2_0, 2, raw = TRUE)")), data = data)
  model_cubic <- lm(as.formula(paste(deparse(Y_value), "~ poly(AgeAttend_2_0, 3, raw = TRUE)")), data = data)
  model_quartic <- lm(as.formula(paste(deparse(Y_value), "~ poly(AgeAttend_2_0, 4, raw = TRUE)")), data = data)
  
  # Éú³ÉÊý¾Ý·¶Î§ÓÃÓÚ»æÍ¼
  age_range <- seq(min(data$AgeAttend_2_0), max(data$AgeAttend_2_0), length.out = 100)
  fit_data <- data.frame(AgeAttend_2_0 = age_range)
  
  # Ô¤²âÖµ
  fit_data$Linear <- predict(model_lin, newdata = fit_data)
  fit_data$Quadratic <- predict(model_quad, newdata = fit_data)
  fit_data$Cubic <- predict(model_cubic, newdata = fit_data)
  fit_data$Quartic <- predict(model_quartic, newdata = fit_data)
  
  # ×ª»»Îª³¤¸ñÊ½ÓÃÓÚ ggplot
  fit_long <- reshape2::melt(fit_data, id.vars = "AgeAttend_2_0", 
                             variable.name = "Model", value.name = "Fitted")
  
  # »æÍ¼
  plot <- ggplot(data, aes(x = AgeAttend_2_0, y = !!Y_value)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_line(data = fit_long, aes(x = AgeAttend_2_0, y = Fitted, color = Model), linewidth = 1) +
    scale_color_manual(values = c("Linear" = "red", "Quadratic" = "green", "Cubic" = "purple", "Quartic" = "orange")) +
    labs(title = "Model Fits for Different Degrees",
         x = "AgeAttend_2_0",
         y = deparse(Y_value),
         color = "Model") +
    theme_minimal()
  
  return(plot)
}


plot_models_with_fits(UKB_redsyn_tmp, Y_col)


library(ggplot2)
library(mgcv)

plot_models_with_fits <- function(data, Y_value) {
  # È·±£ Y_value ÊÇÁÐÃû
  Y_value <- as.name(Y_value)
  
  # ¹¹½¨ GAM Ä£ÐÍ
  model_gam <- gam(as.formula(paste(deparse(Y_value), "~ s(AgeAttend_2_0)")), data = data)
  
  # Éú³ÉÊý¾Ý·¶Î§ÓÃÓÚ»æÍ¼
  age_range <- seq(min(data$AgeAttend_2_0), max(data$AgeAttend_2_0), length.out = 100)
  fit_data <- data.frame(AgeAttend_2_0 = age_range)
  
  # Ô¤²âÖµ
  fit_data$GAM <- predict(model_gam, newdata = fit_data)
  
  # »æÍ¼
  plot <- ggplot(data, aes(x = AgeAttend_2_0, y = !!Y_value)) +
    # geom_point(alpha = 0.6, color = "blue") +
    geom_line(data = fit_data, aes(x = AgeAttend_2_0, y = GAM), color = "red", size = 1) +
    labs(title = "Scatter Plot with GAM Curve",
         x = "AgeAttend_2_0",
         y = deparse(Y_value),
         color = "Model") +
    theme_minimal()
  
  return(plot)
}

plot <- plot_models_with_fits(UKB_redsyn_tmp, Y_col)
print(plot)

# ¶¨Òå7*7µÄÍ¼¿ò²¼¾Ö
layout_matrix <- matrix(1:49, nrow = 7, ncol = 7)

# ´´½¨Ò»¸öÍ¼ÐÎÉè±¸
par(mfrow = c(7, 7), mar = c(2, 2, 2, 1)) # ÉèÖÃÃ¿¸ö×ÓÍ¼µÄ±ß¾à

# Êý¾ÝÁÐÃûºÍ¶Ô½ÇÏß»æÖÆÄ¿±ê
columns <- paste0("Y", 1:28)
index <- 1

for (row in 1:7) {
  for (col in 1:7) {
    if (row == col && index <= 7) {
      # »æÖÆ¶Ô½ÇÏßÉÏµÄÇúÏß Y1-Y7
      plot <- plot_models_with_fits(data = UKB_redsyn_tmp, Y_value = columns[index])
      print(plot)
      index <- index + 1
    } else if (row > col && index <= 28) {
      # »æÖÆÏÂÈý½ÇÉÏµÄÇúÏß Y8-Y28
      plot <- plot_models_with_fits(data = UKB_redsyn_tmp, Y_value = columns[index])
      print(plot)
      index <- index + 1
    } else {
      # ¿Õ°××ÓÍ¼
      plot.new()
    }
  }
}

# ¹Ø±ÕÍ¼ÐÎÉè±¸£¨ÈçÊ¹ÓÃ pdf µÈÊä³öÐèµ÷ÓÃ dev.off()£©








library(ggplot2)
library(lme4)

best_model <- lmer(Y_value ~ poly(AgeAttend_2_0, 4) + (1|Centre_0_0) + Sex_0_0 + HeadMotion_2_0 + BMI_2_0 + 
                     BloodPress_2 + EchoTime_2_0 + IntensityScaling_2_0 + Race_1 + Race_2 + Race_3 + Race_4 +
                     Vol_WB_TIV_2_0+ Handedness_0_0 + TDI_0_0, data = UKB_redsyn_tmp)


UKB_redsyn_tmp$Predicted <- predict(best_model, newdata = UKB_redsyn_tmp, re.form = NA) 

# »æÖÆÉ¢µãÍ¼ + ÄâºÏ¹ì¼£
ggplot(UKB_redsyn_tmp, aes(x = AgeAttend_2_0, y = Y_value)) +
  # geom_point(alpha = 0.6) +  # É¢µãÍ¼
  geom_line(data = UKB_redsyn_tmp, aes(x = AgeAttend_2_0, y = Predicted), color = "blue", size = 1) +  # ÄâºÏ¹ì¼£
  labs(title = "Scatter Plot with Fitted Curve",
       x = "AgeAttend_2_0",
       y = "Y_value") +
  theme_minimal()


UKB_redsyn_tmp_long <- UKB_redsyn_tmp[1:2000,1:101] %>%
  pivot_longer(cols = starts_with("Y"), names_to = "Y_column", values_to = "Y_value")


ggplot(UKB_redsyn_tmp_long, aes(x = AgeAttend_2_0, y = Y_value)) +
  geom_point(alpha = 0.3) +  
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # ÄâºÏÏß
  facet_wrap(~ Y_column, scales = "free_y") +  # °´ Y_variable ·ÖÃæ
  theme_minimal() +  # Ê¹ÓÃ¼ò½àÖ÷Ìâ
  labs(x = "AttendAge", y = "Y_value", title = "Y1-Y58 Scatter Plots with Fitted Lines")  # Ìí¼Ó±êÇ©




result_7n <- read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/modal_biomarker_7n.txt", sep = " ")

result_7n$xname <- result_7n$xname %>%
  gsub("_0_0", "", .) %>%            # È¥µô _0_0
  gsub("_", " ", .) %>%              # °Ñ _ Ìæ»»Îª¿Õ¸ñ
  gsub("^x(\\d+)", "\\1", .)         # È¥µôÒÔ x ¿ªÍ·ºó½ÓÊý×ÖµÄ¿ªÍ· x

names(result_7n)[1] <- "Description"
names(result_7n)[2] <- "varName"




library(readxl)
biomarker_category <- read_excel("ZJLab/UKBiobank_Project/data/table/UKB_Biomarker_Category.xlsx")
biomarker_category$`Field ID` <- as.character(biomarker_category$`Field ID`)

df.biomarker.all <- merge(result_7n,biomarker_category)

df.biomarker <- df.biomarker.all[df.biomarker.all$varName=="RedYeo7n_Total",]
df.biomarker <- df.biomarker[!is.nan(df.biomarker$t),]



df.biomarker$logpvalue <- -log10(df.biomarker$p)
df.biomarker <- df.biomarker %>% arrange(Label, varName)
df.biomarker$varRank <- 1:nrow(df.biomarker)
set1_colors <- brewer.pal(9, "Set1")
set2_colors <- brewer.pal(8, "Set2")
combined_colors <- c(set1_colors, set2_colors)
color_palette <- colorRampPalette(combined_colors)(length(unique(df.biomarker$Label)))


sig_threshold <- 0.05 / nrow(df.biomarker)
top_sig <- df.biomarker[df.biomarker$p < sig_threshold,]

category_labels <- df.biomarker %>%
  group_by(Label) %>%
  summarize(
    start = min(as.numeric(varRank)),  
    end = max(as.numeric(varRank)),    
    midpoint = (start + end) / 2,      
    .groups = 'drop'
  )

ggplot(df.biomarker, aes(x = varRank, y = logpvalue, fill = Label)) +
  geom_point(aes(size = logpvalue, color = Label, shape = beta > 0), alpha = 1) +  
  scale_shape_manual(values = c(25, 24)) + 
  geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "red") +  
  geom_tile(aes(y = -0.15, fill = Label), height = .3) +  
  scale_fill_manual(values = color_palette) +  
  scale_color_manual(values = color_palette) +  
  theme_minimal() + ylim(-2,8) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +  
  labs(y = "-log10(p)", x = element_blank(), 
       title = "Manhattan Plot") +  
  geom_text_repel(data = top_sig, aes(label = Description), size = 3, force = 20, max.overlaps = 100)  





################################################################################
#                  Manhattan Plot (Phenotypes) & Heatmap
################################################################################

result_7n <- read.csv(paste0("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/",
                             "Detecting_Compensation_Effects_of_Aging_Using_PID/modal_pheno_7n.txt"), sep = " ")

# result_7n <- read.csv("ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/modal_adult_pheno_7n.txt",sep = " ")

result_7n$xname <- result_7n$xname %>%
  gsub("_2_0", "", .)

names(result_7n)[1] <- "description"
names(result_7n)[2] <- "varName"


UKB_AllVars_Code <- read_excel("ZJLab/UKBiobank_Project/data/table/UKB_AllVarUsed.xlsx")

UKB_AllVars_Code$description <- gsub(" ","_",UKB_AllVars_Code$Description)
UKB_AllVars_Code$description <- gsub("/","_or_",UKB_AllVars_Code$description)
UKB_AllVars_Code$description <- gsub("-","_",UKB_AllVars_Code$description)
UKB_AllVars_Code$description <- gsub("'","_",UKB_AllVars_Code$description)
UKB_AllVars_Code$description <- gsub("\\(|\\)","",UKB_AllVars_Code$description)
UKB_AllVars_Code$description <- gsub(" ","",UKB_AllVars_Code$description)
UKB_AllVars_Code$description <- gsub("__","_",UKB_AllVars_Code$description)
UKB_AllVars_Code$description <- gsub(",","",UKB_AllVars_Code$description)
UKB_AllVars_Code$description <- gsub("\\+","plus",UKB_AllVars_Code$description)

df.pheno.all <- merge(result_7n,UKB_AllVars_Code,by = "description")

df.pheno.all <- df.pheno.all[c(which(df.pheno.all$Category=="Cognitive function"),
                               which(df.pheno.all$Category=="Mental health"),
                               which(df.pheno.all$Category=="Self-reported medical conditions")),]

df.pheno <- df.pheno.all[df.pheno.all$varName=="SynYeo7n_Total",]
df.pheno <- na.omit(df.pheno)
df.pheno$logpvalue <- -log10(df.pheno$p)
df.pheno <- df.pheno %>% arrange(Category, varName)
df.pheno$varRank <- 1:nrow(df.pheno)

combined_colors <- brewer.pal(3, "Set1")
color_palette <- colorRampPalette(combined_colors)(length(unique(df.pheno$Category)))

sig_threshold <- 0.05 / nrow(df.pheno)
top_sig <- df.pheno[df.pheno$p < sig_threshold,]

category_labels <- df.pheno %>%
  group_by(Category) %>%
  summarize(
    start = min(as.numeric(varRank)),  
    end = max(as.numeric(varRank)),    
    midpoint = (start + end) / 2,      
    .groups = 'drop'
  )

ggplot(df.pheno, aes(x = varRank, y = logpvalue, fill = Category)) +
  geom_point(aes(size = logpvalue, color = Category, shape = beta > 0), alpha = 1) +  
  scale_shape_manual(values = c(25, 24)) + 
  geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "red") +  
  geom_tile(aes(y = -0.15, fill = Category), height = .3) +  
  scale_fill_manual(values = color_palette) +  
  scale_color_manual(values = color_palette) +  
  theme_minimal() + ylim(-2,18) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +  
  labs(y = "-log10(p)", x = element_blank(), 
       title = "Manhattan Plot") +  
  geom_text_repel(data = top_sig, aes(label = Description), size = 3, force = 20, max.overlaps = 100)  

# Heatmap
df.pheno.all
df.pheno.all.mat <- t(matrix(df.pheno.all$t,58,4350/58))*t(matrix((df.pheno.all$p<0.05/74),58,4350/58))
colnames(df.pheno.all.mat) <- df.pheno.all$varName[1:58]
rownames(df.pheno.all.mat) <- df.pheno.all$description[seq(1,4350,58)]
df.pheno.all.mat <- df.pheno.all.mat[-1*which(rownames(df.pheno.all.mat) == "Pregnant"),]


annotation_col <- df.pheno.all[seq(1,4350,58),c(12,13)]
rownames(annotation_col) <- df.pheno.all$description[seq(1,4350,58)]
ann_colors <- list(
  Category = c("Cognitive function" = "#E41A1C",
               "Mental health" = "#377EB8",
               "Self-reported medical conditions"="#4DAF4A")
)

# °´ÕÕ annotation_col ÀïµÄ Group ±äÁ¿ÅÅÐò
new_order <- order(annotation_col$Category)  

# ÖØÐÂÅÅÁÐÊý¾ÝÁÐ
sorted_mat <- df.pheno.all.mat[new_order]

# ÖØÐÂÅÅÁÐ annotation_col ÐÐµÄË³Ðò
annotation_col <- annotation_col[new_order, , drop = FALSE]




annotation_col2 <- data.frame(Category = c(df.pheno.all[seq(1,4350,58),c(13)]))
rownames(annotation_col2) <- df.pheno.all$description[seq(1,4350,58)]


x = t(df.pheno.all.mat[,(1:7)])
rownames(x) <- c("VN","SMN","DAN","VAN","LN","FPN","DMN")
pheatmap::pheatmap(x,cluster_rows = F,cluster_col = F,
                   color = colorRampPalette(c("blue", "white", "red"))(50), 
                   breaks = seq(-9, 9, length.out = 51),
                   annotation_col = annotation_col2,
                   annotation_colors = ann_colors)
x = t(df.pheno.all.mat[,(1:7)+29])
rownames(x) <- c("VN","SMN","DAN","VAN","LN","FPN","DMN")
pheatmap::pheatmap(x,cluster_rows = F,cluster_col = F,
                   color = colorRampPalette(c("blue", "white", "red"))(50), 
                   breaks = seq(-9, 9, length.out = 51),
                   annotation_col = annotation_col2,
                   annotation_colors = ann_colors)




x = t(df.pheno.all.mat[,(8:28)])
rownames(x) <- c('VN-SMN','VN-DAN','VN-VAN','VN-LN','VN-FPN','VN-DMN','SMN-DAN','SMN-VAN',
                 'SMN-LN','SMN-FPN','SMN-DMN','DAN-VAN','DAN-LN','DAN-FPN','DAN-DMN','VAN-LN',
                 'VAN-FPN','VAN-DMN','LN-FPN','LN-DMN','FPN-DMN')
pheatmap::pheatmap(x,cluster_rows = F,cluster_col = F,
                   color = colorRampPalette(c("blue", "white", "red"))(50), 
                   breaks = seq(-9, 9, length.out = 51),
                   annotation_col = annotation_col2,
                   annotation_colors = ann_colors)

# 1300 * 550
# 1300 * 800
x = t(df.pheno.all.mat[,29+(8:28)])
rownames(x) <- c('VN-SMN','VN-DAN','VN-VAN','VN-LN','VN-FPN','VN-DMN','SMN-DAN','SMN-VAN',
                 'SMN-LN','SMN-FPN','SMN-DMN','DAN-VAN','DAN-LN','DAN-FPN','DAN-DMN','VAN-LN',
                 'VAN-FPN','VAN-DMN','LN-FPN','LN-DMN','FPN-DMN')

pheatmap::pheatmap(x,cluster_rows = F,cluster_col = F,
                   color = colorRampPalette(c("blue", "white", "red"))(50), 
                   breaks = seq(-9, 9, length.out = 51),
                   annotation_col = annotation_col2,
                   annotation_colors = ann_colors)






















################################################################################
#                      Manhattan Plot (Biomarker) & Heatmap
################################################################################


result_7n <- read.csv(paste0("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/",
                             "Detecting_Compensation_Effects_of_Aging_Using_PID/modal_biomarker_7n.txt"), sep = " ")



result_7n$xname <- result_7n$xname %>% gsub("_0_0", "", .)
result_7n$xname <- result_7n$xname %>% gsub("virus_1", "virus", .)
result_7n$xname <- result_7n$xname %>% gsub("x1gG", "1gG", .)
result_7n$xname <- result_7n$xname %>% gsub("x2mgG", "2mgG", .)
result_7n$ID = 1:nrow(result_7n)
biomarker_category <- read_excel("ZJLab/UKBiobank_Project/data/table/UKB_Biomarker_Category.xlsx")
biomarker_category$`Field ID` <- as.character(biomarker_category$`Field ID`)
biomarker_category$xname <- biomarker_category$Description %>% gsub("-", "_", .)
biomarker_category$xname <- biomarker_category$xname  %>% gsub(" ", "_", .)

df.biomarker.all <- merge(result_7n,biomarker_category)

result_7n$yname[1:58]

ROI = result_7n$yname[1:58]
Label = df.biomarker.all$Label[df.biomarker.all$yname=="SynYeo7n_Total"]
Description = df.biomarker.all$Description[df.biomarker.all$yname=="SynYeo7n_Total"]
Description[order(Label)]

df.biomarker.all <- df.biomarker.all[order(df.biomarker.all$ID),]

combined_colors <- c(set1_colors, set2_colors)
color_palette <- colorRampPalette(combined_colors)(length(unique(df.biomarker$Label)))



df.biomarker.all
df.biomarker.all.mat <- t(matrix(df.biomarker.all$t,58,nrow(df.biomarker.all)/58))*
  t(matrix((df.biomarker.all$p<0.05/119),58,nrow(df.biomarker.all)/58))
colnames(df.biomarker.all.mat) <- df.biomarker.all$yname[1:58]
rownames(df.biomarker.all.mat) <- df.biomarker.all$Description[seq(1,nrow(df.biomarker.all),58)]
# df.biomarker.all.mat <- df.biomarker.all.mat[-1*which(rownames(df.biomarker.all.mat) == "Pregnant"),]
df.biomarker.all.mat <- df.biomarker.all.mat[Description[order(Label)], ]
df.biomarker.all.mat <- df.biomarker.all.mat[which(!is.nan(df.biomarker.all.mat[,1])),]

annotation_col <- df.biomarker.all[seq(1,nrow(df.biomarker.all),58),c(13,14)]

rownames(annotation_col) <- df.biomarker.all$Description[seq(1,nrow(df.biomarker.all),58)]
new_order <- order(annotation_col$Label)  
sorted_mat <- df.biomarker.all.mat[new_order]
annotation_col <- annotation_col[new_order, , drop = FALSE]
annotation_col2 <- data.frame(Label = c(df.biomarker.all[seq(1,nrow(df.biomarker.all),58),c(14)]))
rownames(annotation_col2) <- df.biomarker.all$Description[seq(1,nrow(df.biomarker.all),58)]


unique_labels <- unique(annotation_col$Label)
label_colors <- setNames(combined_colors[1:length(unique_labels)], unique_labels)
ann_colors <- list(Label = label_colors)


# df.biomarker.all.mat <- df.biomarker.all.mat[!is.nan(df.biomarker.all.mat[,1]),]
df.biomarker.all.mat <- df.biomarker.all.mat[,ROI]

x = t(df.biomarker.all.mat[,(1:7)])
rownames(x) <- c("VN","SMN","DAN","VAN","LN","FPN","DMN")
pheatmap::pheatmap(x,cluster_rows = F,cluster_col = F,
                   color = colorRampPalette(c("blue", "white", "red"))(50), 
                   breaks = seq(-9, 9, length.out = 51),
                   annotation_col = annotation_col2,
                   annotation_colors = ann_colors)
# syn_biomark_in
# 2000 500
x = t(df.biomarker.all.mat[,(1:7)+29])
rownames(x) <- c("VN","SMN","DAN","VAN","LN","FPN","DMN")
pheatmap::pheatmap(x,cluster_rows = F,cluster_col = F,
                   color = colorRampPalette(c("blue", "white", "red"))(50), 
                   breaks = seq(-9, 9, length.out = 51),
                   annotation_col = annotation_col2,
                   annotation_colors = ann_colors)




x = t(df.biomarker.all.mat[,(8:28)])
rownames(x) <- c('VN-SMN','VN-DAN','VN-VAN','VN-LN','VN-FPN','VN-DMN','SMN-DAN','SMN-VAN',
                 'SMN-LN','SMN-FPN','SMN-DMN','DAN-VAN','DAN-LN','DAN-FPN','DAN-DMN','VAN-LN',
                 'VAN-FPN','VAN-DMN','LN-FPN','LN-DMN','FPN-DMN')
pheatmap::pheatmap(x,cluster_rows = F,cluster_col = F,
                   color = colorRampPalette(c("blue", "white", "red"))(50), 
                   breaks = seq(-9, 9, length.out = 51),
                   annotation_col = annotation_col2,
                   annotation_colors = ann_colors)

# 1300 * 550
# 1300 * 800
x = t(df.biomarker.all.mat[,29+(8:28)])
rownames(x) <- c('VN-SMN','VN-DAN','VN-VAN','VN-LN','VN-FPN','VN-DMN','SMN-DAN','SMN-VAN',
                 'SMN-LN','SMN-FPN','SMN-DMN','DAN-VAN','DAN-LN','DAN-FPN','DAN-DMN','VAN-LN',
                 'VAN-FPN','VAN-DMN','LN-FPN','LN-DMN','FPN-DMN')

pheatmap::pheatmap(x,cluster_rows = F,cluster_col = F,
                   color = colorRampPalette(c("blue", "white", "red"))(50), 
                   breaks = seq(-9, 9, length.out = 51),
                   annotation_col = annotation_col2,
                   annotation_colors = ann_colors)





################################################################################
#                                 Mediation
################################################################################


install.packages("mediation")
install.packages("Matrix")
library(mediation)   # ÖÐ½é·ÖÎö
library(lme4)        # »ìºÏÄ£ÐÍÖ§³ÖËæ»úÐ§Ó¦
library(dplyr)       # Êý¾Ý²Ù×÷
library(Matrix)

# ¼ÙÉèÄãµÄÊý¾ÝÒÑ¾­ÊÇÒ»¸ö dataframe: data_tmp
# ²¢ÇÒ°üº¬ AgeAttend_2_0, ËùÓÐµÄ M, Y ºÍ covariates ÁÐ
YList <- c("SynYeo7n_Total", "RedYeo7n_Total")
MList <- c("WMHyperInt_WB_Total_2_0", "FA_WB_Total_2_0", "MD_WB_Total_2_0",
           "MO_WB_Total_2_0", "ICVF_WB_Total_2_0", "ISOVF_WB_Total_2_0")
COVList <- c("Sex_0_0", "HeadMotion_2_0", "SNR_2_0", "TDI_0_0",
             "BMI_2_0", "Race_1", "Race_2", "Race_3", "Race_4", "Vol_WB_TIV_2_0")

# optional: convert categorical variables
data_tmp <- merge(merge(UKB_BasicMRI,UKB_RedSyn_2_0),UKB_BrainImagingVars)
data_tmp$Sex_0_0 <- factor(data_tmp$Sex_0_0)
data_tmp$Centre_0_0 <- factor(data_tmp$Centre_0_0)  # ¿ØÖÆÕ¾µã£¨½¨Òé×÷ÎªËæ»úÐ§Ó¦£©Centre_0_0
data_tmp$Vol_WB_TIV_2_0 <- data_tmp$Vol_WB_TIV_2_0/100000
data_tmp$Vol_WB_TIV_3_0 <- data_tmp$Vol_WB_TIV_3_0/100000
library(mediation)

library(fastDummies)
library(dplyr)

# Í³¼ÆÃ¿¸öÕ¾µãµÄÑù±¾Á¿
centre_counts <- data_tmp %>%
  group_by(Centre_0_0) %>%
  summarise(n = n())
# ´´½¨ dummy variables
data_tmp <- dummy_cols(data_tmp, select_columns = "Centre_0_0", remove_selected_columns = TRUE)

# ±£ÁôÑù±¾Á¿ >= 1000 µÄÕ¾µã
valid_centres <- centre_counts %>% filter(n >= 1000) %>% pull(Centre_0_0)

# ¶ÔÓ¦µÄ dummy ÁÐÃû
centre_dummy_vars <- paste0("Centre_0_0_", valid_centres)

# È·±£ÕâÐ© dummy ÁÐÈ·Êµ´æÔÚÓÚÊý¾ÝÖÐ
centre_dummy_vars <- centre_dummy_vars[centre_dummy_vars %in% colnames(data_tmp)]

YList <- c('SynYeo7n_Total', 'RedYeo7n_Total')
MList <- c('WMHyperInt_WB_Total_2_0','FA_WB_Total_2_0','MD_WB_Total_2_0','MO_WB_Total_2_0',
           'ICVF_WB_Total_2_0','ISOVF_WB_Total_2_0')
COVList <- c('Sex_0_0','HeadMotion_2_0','SNR_2_0','TDI_0_0',
             'BMI_2_0','Race_1','Race_2','Race_3','Race_4','Vol_WB_TIV_2_0', centre_dummy_vars)

results_all <- list()

for (m in MList) {
  for (y in YList) {
    vars_needed <- c("AgeAttend_2_0", y, m, COVList)
    df <- data_tmp[, vars_needed]
    df <- na.omit(df)
    colnames(df)[1:3] <- c("X", "Y", "M")
    
    model.m <- lm(M ~ X + ., data = df)
    model.y <- lm(Y ~ M + X + ., data = df)
    model.c <- lm(Y ~ X + ., data = df)
    
    med.out <- mediate(model.m, model.y, treat = "X", mediator = "M",
                       boot = TRUE, sims = 500)
    
    # ÌáÈ¡lmÄ£ÐÍÏµÊýÕªÒª
    a_sum <- summary(model.m)$coefficients["X", ]
    b_sum <- summary(model.y)$coefficients["M", ]
    c_sum <- summary(model.c)$coefficients["X", ]
    cprime_sum <- summary(model.y)$coefficients["X", ]
    
    res <- data.frame(
      X = "AgeAttend_2_0",
      M = m,
      Y = y,
      
      a_estimate = a_sum["Estimate"],
      a_std_error = a_sum["Std. Error"],
      a_t_value = a_sum["t value"],
      a_p_value = a_sum["Pr(>|t|)"],
      
      b_estimate = b_sum["Estimate"],
      b_std_error = b_sum["Std. Error"],
      b_t_value = b_sum["t value"],
      b_p_value = b_sum["Pr(>|t|)"],
      
      ab_manual = a_sum["Estimate"] * b_sum["Estimate"],
      ab_boot = med.out$d0,
      ab_diff = a_sum["Estimate"] * b_sum["Estimate"] - med.out$d0,
      ab_se = sd(med.out$d0.sims),
      ab_p = med.out$d0.p,
      ab_CI_lower = quantile(med.out$d0.sims, 0.025),
      ab_CI_upper = quantile(med.out$d0.sims, 0.975),
      
      c_estimate = c_sum["Estimate"],
      c_std_error = c_sum["Std. Error"],
      c_t_value = c_sum["t value"],
      c_p_value = c_sum["Pr(>|t|)"],
      
      cprime_estimate = cprime_sum["Estimate"],
      cprime_std_error = cprime_sum["Std. Error"],
      cprime_t_value = cprime_sum["t value"],
      cprime_p_value = cprime_sum["Pr(>|t|)"],
      
      # med.out ÖÐÕæÕýµÄ¹À¼ÆÖµ£¨·ÇzÍ³¼ÆÁ¿£©
      total_effect = med.out$tau.coef,
      total_p = med.out$tau.p,
      total_CI_lower = quantile(med.out$tau.sims, probs = 0.025),
      total_CI_upper = quantile(med.out$tau.sims, probs = 0.975),
      
      direct_effect = med.out$zeta.coef,
      direct_p = med.out$zeta.p,
      direct_CI_lower = quantile(med.out$zeta.sims, probs = 0.025),
      direct_CI_upper = quantile(med.out$zeta.sims, probs = 0.975),
      
      # ÒÔ¼° z Í³¼ÆÁ¿£¨±ê×¼»¯Í³¼ÆÁ¿£©
      total_z_stat = med.out$z0,
      total_z_p = med.out$z0.p,
      direct_z_stat = med.out$z1,
      direct_z_p = med.out$z1.p,
      
      N = nrow(df)
    )
    
    key <- paste(m, y, sep = "_")
    results_all[[key]] <- res
  }
}

results_df <- do.call(rbind, results_all)

# ²é¿´½á¹û
print(results_df)



library(lavaan)
library(dplyr)

# ÉèÖÃ±äÁ¿
YList <- c('SynYeo7n_Total', 'RedYeo7n_Total')
MList <- c('WMHyperInt_WB_Total_2_0','FA_WB_Total_2_0','MD_WB_Total_2_0',
           'MO_WB_Total_2_0','ICVF_WB_Total_2_0','ISOVF_WB_Total_2_0')
COVList <- c('Sex_0_0','HeadMotion_2_0','SNR_2_0','TDI_0_0',
             'BMI_2_0','Race_1','Race_2','Race_3','Race_4','Vol_WB_TIV_2_0', centre_dummy_vars)
Xvar <- "AgeAttend_2_0"

# Í³Ò»´¦ÀíÈ±Ê§Öµ
all_vars <- unique(c(Xvar, YList, MList, COVList))
data_tmp_complete <- data_tmp[complete.cases(data_tmp[, all_vars]), ]

results_lavaan <- list()

for (m in MList) {
  for (y in YList) {
    
    # ÌáÈ¡Êý¾Ý
    vars_needed <- c(Xvar, m, y, COVList)
    df <- data_tmp_complete[, vars_needed]
    colnames(df)[1:3] <- c("X", "M", "Y")  # ±ê×¼»¯ÃüÃû
    
    # ¹¹Ôì lavaan Ä£ÐÍ
    # ×¢Òâ£º":=" Óï·¨¶¨Òå¼ä½ÓÐ§Ó¦¡¢Ö±½ÓÐ§Ó¦¡¢×ÜÐ§Ó¦
    model_text <- paste0("
      M ~ a*X + ", paste(COVList, collapse = " + "), "
      Y ~ b*M + c_prime*X + ", paste(COVList, collapse = " + "), "
      ab := a * b
      total := c_prime + (a * b)
    ")
    
    fit <- sem(model_text, data = df, se = "bootstrap", bootstrap = 50)
    summ <- parameterEstimates(fit, ci = TRUE, level = 0.95, standardized = TRUE)
    
    # ÌáÈ¡¸ÐÐËÈ¤µÄ½á¹û£¨a, b, ab, c', total£©
    est <- summ %>% 
      filter(op %in% c("~", ":=")) %>% 
      select(lhs, rhs, est, se, pvalue, ci.lower, ci.upper)
    
    est$X <- XVar
    est$M <- m
    est$Y <- y
    
    results_lavaan[[paste(m, y, sep = "_")]] <- est
  }
}

# ºÏ²¢½á¹û
lavaan_results_df <- bind_rows(results_lavaan)
# »òÕß write.csv(results_df, "mediation_results.csv", row.names = FALSE)




################################################################################
#                                 Bubble Plot
################################################################################





library(ggplot2)
library(dplyr)
lifestyle_tmp <- read_csv("ZJLab/UKBiobank_Project/data/table/UKB_Lifestyle_List.csv")

data_tmp <- read.csv('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/modal_lifestyle_7n.txt',
                     na = "NaN",sep = " ")
library(readr)
lifestyle_tmp[66,"xname"] <- "x1160_S"
lifestyle_tmp[66,"Description"] <- "Sleep duration (<= 6 hours)"
lifestyle_tmp[67,"xname"] <- "x1160_L"
lifestyle_tmp[67,"Description"] <- "Sleep duration (>= 9 hours)"
lifestyle_tmp$xname <- paste0("x",lifestyle_tmp$`Field ID`)
data_tmp <- merge(data_tmp,lifestyle_tmp)
data_tmp$Description[data_tmp$Description=="Avoided activities or situations because of previous stressful experience in past month"] <- 
  "Avoided activities because of previous stress in past month"

# ¼ÙÉèÄãµÄÊý¾Ý¿òÎª data_tmp£¬°üÀ¨ xname, yname, t, p
# 1. ¶¨ÒåÏÔÖøÐÔãÐÖµ
alpha <- 0.05 / 14

# 2. ¹ýÂËµôÃ»ÓÐÈÎºÎÏÔÖø½á¹ûµÄ xname
sig_x <- data_tmp %>%
  group_by(xname) %>%
  summarise(any_sig = any(p < alpha), .groups = "drop") %>%
  filter(any_sig) %>%
  pull(xname)

data_plot <- data_tmp %>%
  filter(xname %in% sig_x)

data_plot$yname <- factor(data_plot$yname, levels = rev(c(paste0("Red7n_NC",1:7),paste0("Syn7n_NC",1:7))))
# °´ Category_sub µÄË³ÐòÅÅÁÐ Description
data_plot <- data_plot %>%
  arrange(Category) %>%                       # °´ Category_sub ÅÅÐò
  mutate(Description = factor(Description, levels = unique(Description)))  # ÉèÖÃ factor Ë³Ðò



data_plot <- data_plot[data_plot$p<alpha,]

p <- ggplot(data_plot, aes(x = Description, y = yname)) +
  geom_point(aes(
    size = abs(t),      # ÆøÅÝ´óÐ¡
    fill = t,           # Ìî³äÑÕÉ«
    shape = t > 0       # Èý½ÇÐÎ·½Ïò
  ),
  color = "gray80", stroke = 0) +
  scale_size_continuous(range = c(2, 10)) +
  scale_fill_distiller(palette = "RdBu", 
                       limits = c(-10, 10), 
                       direction = -1) +   # direction=1: À¶-°×-ºì£»= -1 ·´×ª
  scale_shape_manual(values = c(25, 24)) +  # FALSE=ÏÂÈý½Ç, TRUE=ÉÏÈý½Ç
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_line(color = "grey80")
  ) +
  labs(
    x = "Environmental Factor & Lifestyle",
    y = "Yeo 7 Networks",
    size = "|t|",
    fill = "t value"
  )


print(p)





data_tmp <- read.csv('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/modal_outcome_7nx.txt',
                     na = "NaN",sep = " ")
outcome_tmp <- read_csv("ZJLab/UKBiobank_Project/data/table/UKB_Outcomes_List.csv")
outcome_tmp$xname <- paste0("x",outcome_tmp$`Field ID`)
data_tmp <- merge(data_tmp,outcome_tmp)

alpha <- 0.05 / 14

# 2. ¹ýÂËµôÃ»ÓÐÈÎºÎÏÔÖø½á¹ûµÄ xname
sig_x <- data_tmp %>%
  group_by(xname) %>%
  summarise(any_sig = any(p < alpha), .groups = "drop") %>%
  filter(any_sig) %>%
  pull(xname)

data_plot <- data_tmp %>%
  filter(xname %in% sig_x)

data_plot$yname <- factor(data_plot$yname, levels = rev(c(paste0("Red7n_NC",1:7),paste0("Syn7n_NC",1:7))))
# °´ Category_sub µÄË³ÐòÅÅÁÐ Description
data_plot <- data_plot %>%
  arrange(Category) %>%                       # °´ Category_sub ÅÅÐò
  mutate(Description = factor(Description, levels = unique(Description)))  # ÉèÖÃ factor Ë³Ðò



data_plot <- data_plot[data_plot$p<alpha,]

library(RColorBrewer)

p <- ggplot(data_plot, aes(x = Description, y = yname)) +
  geom_point(aes(
    size = abs(t),      
    fill = t,           
    shape = t > 0       
  ),
  color = "gray80", stroke = 0) +
  scale_size_continuous(range = c(2, 10)) +
  scale_fill_distiller(palette = "RdBu", 
                       limits = c(-10, 10), 
                       direction = -1) +   # direction=1: À¶-°×-ºì£»= -1 ·´×ª
  scale_shape_manual(values = c(25, 24)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_line(color = "grey80")
  ) +
  labs(
    x = "Cognition & Mental / Physic Health",
    y = "Yeo 7 Networks",
    size = "|t|",
    fill = "t value"
  )

print(p)



data_tmp <- read.csv('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/mediation_results_7nx.txt',
                     na = "NaN",sep = ",")

data_tmp <- merge(data_tmp,rbind(outcome_tmp,lifestyle_tmp)[,c(2,5)],by.x = "X",by.y = "xname")
names(data_tmp)[11] <- "XDescription"
data_tmp <- merge(data_tmp,rbind(outcome_tmp,lifestyle_tmp)[,c(2,5)],by.x = "Y",by.y = "xname")
names(data_tmp)[12] <- "YDescription"

data_tmp_ab <- data_tmp[data_tmp$Path=="ab"&data_tmp$YDescription=="Overall health rating",] # Happiness;Overall health rating
data_tmp_ab <- data_tmp_ab[-1*(which(data_tmp_ab$XDescription=="Bread type")),]
data_tmp_ab$M <- factor(data_tmp_ab$M, levels = rev(c('red-VN','red-SMN','red-DAN',
                                                      'red-VAN','red-LN','red-FPN',
                                                      'red-DMN','syn-VN','syn-SMN',
                                                      'syn-DAN','syn-VAN','syn-LN',
                                                      'syn-FPN','syn-DMN')))
alpha = 0.05/(7*9)
data_tmp_ab_sig <- data_tmp_ab[data_tmp_ab$p<alpha,]



sig_x <- data_tmp_ab %>%
  group_by(XDescription) %>%
  summarise(any_sig = any(p < alpha), .groups = "drop") %>%
  filter(any_sig) %>%
  pull(XDescription)
data_tmp_ab_sig <- data_tmp_ab %>%
  filter(XDescription %in% sig_x)
data_tmp_ab_sig <- data_tmp_ab_sig[stringr::str_detect(data_tmp_ab_sig$M,"red-"),]
data_tmp_ab_sig$t[data_tmp_ab_sig$p>=alpha] <- NA


ggplot(data_tmp_ab_sig, aes(x = XDescription, y = M)) + 
  geom_point(aes(
    size = abs(t),        # ÆøÅÝ´óÐ¡±íÊ¾ |t|
    fill = t,             # ÑÕÉ«±íÊ¾ t ÖµÕý¸º
    shape = t > 0         # TRUE/1 ÉÏÈý½Ç£¬FALSE/0 ÏÂÈý½Ç
  ),
  color = "gray80", stroke = 0) +
  scale_size_continuous(range = c(2, 10)) +
  scale_fill_distiller(palette = "RdBu", 
                       limits = c(-5, 5), 
                       direction = -1) +
  scale_shape_manual(values = c(25, 24)) + # ÏÂÈý½Ç/ÉÏÈý½Ç
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_line(color = "grey80")
  ) +
  labs(
    x = "Environmental Factor & Lifestyle",
    y = "Yeo 7 Networks",
    size = "|t|",
    fill = "t value"
  )



################################################################################
#                               Cox Regression
################################################################################

UKB_PID_Yeo7n = read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Yeo7n.txt",na = "NaN")



library(survival)
library(dplyr)
library(survminer) # ¼ÓÔØ°ü

perform_cox_regression <- function(data, x_col, baseline_col, event_col, covariates) {
  data$SurvivalTime <- as.numeric(difftime(data[[event_col]], data[[baseline_col]], units = "days")) / 365.25
  data$Status = as.numeric(!is.na(data[[event_col]]))
  data$SurvivalTime[is.na(data$SurvivalTime)] <- max(data$SurvivalTime,na.rm = T)
  data_filtered <- data %>%
    filter(SurvivalTime >= 0)
  
  formula <- as.formula(
    paste0("Surv(SurvivalTime, Status) ~ ", paste(c(x_col, covariates), collapse = " + "))
  )
  
  cox_model <- coxph(formula, data = data_filtered)
  summary_res <- summary(cox_model)
  
  results <- data.frame(
    Variable = rownames(summary_res$coefficients),
    coef = summary_res$coefficients[, "coef"],
    HazardRatio = summary_res$coefficients[, "exp(coef)"],
    SE = summary_res$coefficients[, "se(coef)"],
    Z = summary_res$coefficients[, "z"],
    pValue = summary_res$coefficients[, "Pr(>|z|)"],
    lowerCI = summary_res$conf.int[, "lower .95"],
    upperCI = summary_res$conf.int[, "upper .95"]
  )
  
  x_results <- results %>% filter(Variable == x_col)
  
  return(x_results)
}

library(survminer)
plot_surv_curves <- function(data, x_col, baseline_col, event_col, covariates) {
  data$SurvivalTime <- as.numeric(difftime(data[[event_col]], data[[baseline_col]], units = "days")) / 365.25
  data$Status = as.numeric(!is.na(data[[event_col]]))
  data$SurvivalTime[is.na(data$SurvivalTime)] <- max(data$SurvivalTime,na.rm = T)
  data_filtered <- data %>%
    filter(SurvivalTime >= 0)
  
  # ¼ÆËã 0%, 20%, 40%, 60%, 80%, 100% µÄ·ÖÎ»Êý
  quantiles <- quantile(data_filtered[, x_col], probs = c(0, 0.333, 0.666, 1), na.rm = TRUE)
  
  # Ê¹ÓÃ cut() º¯Êý¸ù¾Ý·ÖÎ»Êý·Ö×é
  data_filtered$cut <- cut(data_filtered[, x_col], 
                           breaks = quantiles, 
                           labels = c("low", "mid", "high"), 
                           include.lowest = TRUE)
  
  
  
  fit <- survfit(Surv(SurvivalTime, Status) ~ cut, data = data_filtered)
  ggsurvplot(fit, data = data_filtered,
             surv.median.line = "none",
             conf.int = TRUE,
             pval = TRUE,
             palette = "aaas",
             ylim = c(0.988, 1))
}


### Data Sort ###


UKB_ICD10_Info <- read_csv("ZJLab/UKBiobank_Project/data/table/UKB_First_Occurrences.csv")
UKB_MRItmp <- read_csv("ZJLab/UKBiobank_Project/data/table/UKB_BrainMRICov.txt")
UKB_Basic <- read_delim("ZJLab/UKBiobank_Project/data/table/UKB_Basic.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)

UKB_Attend_Date = read.ukb.table(53, UKB_FieldID_Subset)
names(UKB_Attend_Date) <- c("eid",paste0("AttendDate_",0:3,"_0"))


UKB_ICD10_Diseases_Date = read.ukb.table(as.numeric(UKB_ICD10_Info$Category_ID_sub[UKB_ICD10_Info$Data_Type=="Date"]),UKB_FieldID_Subset)
UKB_ICD10_Diseases_Names = ukb.icd10.finder(UKB_ICD10_Diseases_Date,UKB_ICD10_Info[UKB_ICD10_Info$Data_Type=="Date",])
UKB_ICD10_Diseases_Names
names(UKB_ICD10_Diseases_Date) <- UKB_ICD10_Diseases_Names

nomissing_num <- colSums(!is.na(UKB_ICD10_Diseases_Date[,1:ncol(UKB_ICD10_Diseases_Date)]))
cols_to_remove <- names(nomissing_num[nomissing_num < 1000])
UKB_ICD10_Diseases_Date <- UKB_ICD10_Diseases_Date[, !(names(UKB_ICD10_Diseases_Date) %in% cols_to_remove)]

UKB_ICD10_MRItmp <- merge(merge(UKB_MRItmp,UKB_Attend_Date),UKB_ICD10_Diseases_Date)

library(dplyr)
UKB_ICD10_MRItmp <- UKB_ICD10_MRItmp %>%
  mutate(across(A04:Q82, ~ case_when(
    is.na(.) ~ 0,  
    . < AttendDate_2_0 ~ 1,  
    TRUE ~ NA_real_  
  ), .names = "{.col}_Past"))


UKB_ICD10_MRItmp <- UKB_ICD10_MRItmp[,c(1,which(str_detect(names(UKB_ICD10_MRItmp),"_Past")))]

colSums(UKB_ICD10_MRItmp[,2:467],na.rm = T)

write.csv(UKB_ICD10_MRItmp,file = "/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_ICD10_Past.csv",na = "",row.names = F)

UKB_ICD10_MRItmp <- read.csv("/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_ICD10_Past.csv")

### Cox Regression ###


UKB_ICD10_MRItmp <- merge(merge(UKB_MRItmp,UKB_Attend_Date),UKB_ICD10_Diseases_Date)
UKB_data_MRItmp <- merge(UKB_ICD10_MRItmp,UKB_PID_Yeo7n[,c(1,95:108)])

UKB_data_MRItmp[c("Sex_0_0","Race_1","Race_2","Race_3","Race_4","Centre_0_0")] <-
  lapply(UKB_data_MRItmp[c("Sex_0_0","Race_1","Race_2","Race_3","Race_4","Centre_0_0")] , factor)


plot_surv_curves(
  UKB_data_MRItmp, "Redundancy_VN", "AttendDate_2_0", "F32",
  c("AgeAttend_2_0", "BMI_2_0", "HeadMotion_2_0", "Race_1", 
    "Race_2", "Race_3", "Race_4", "TDI_0_0", "Vol_WB_TIV_2_0", "Centre_0_0")
)

model_cox_7n <- data.frame()
XList = names(UKB_PID_Yeo7n)[95:108]
YList = names(UKB_ICD10_MRItmp[,99:564])
YList <- YList[which(colSums(!is.na(UKB_ICD10_MRItmp[,99:564]))>1000)]
model_cox_7n <- data.frame() 

for (i in XList) {
  for (j in YList) {
    
    model_tmp <- try(perform_cox_regression(
      UKB_data_MRItmp, i, "AttendDate_2_0", j,
      c("AgeAttend_2_0", "BMI_2_0", "HeadMotion_2_0", "Race_1", 
        "Race_2", "Race_3", "Race_4", "Vol_WB_TIV_2_0", "Centre_0_0")
    ), silent = TRUE)
    
    if (inherits(model_tmp, "try-error")) {
      print(paste("Error in Cox Regression of", j, "~", i))
      next  
    }
    
    try({
      model_tmp$X <- i
      model_tmp$Y <- j
      model_cox_7n <- rbind(model_cox_7n, model_tmp)
      print(paste("Cox Regression of", j, "~", i, "completed"))
    }, silent = TRUE) 
  }
}

model_cox_7n[(model_cox_7n$pValue<0.05)&stringr::str_detect(model_cox_7n$Y,"F"),]




################################################################################
#                                Moderation
################################################################################


library(lme4)
library(ggplot2)
library(dplyr)

Yeo7Names = c("VN","SMN","DAN","VAN","LN","FPN","DMN")
Yeo7Colormap
Yeo7Color
# YList <- paste0("Redundancy_",Yeo7Names)
YList <- paste0("Synergy_",Yeo7Names)
UKB_data_MRItmp <- merge(merge(UKB_BasicMRI,UKB_Fluid_intelligence),UKB_PID_Yeo7n[,c(1,95:108)])
UKB_data_MRItmp <- na.omit(UKB_data_MRItmp[,c("eid","Fluid_intelligence_score_2_0","AgeAttend_2_0",
                                              "Sex_0_0","BMI_2_0","HeadMotion_2_0","Race_1",
                                              "Race_2","Race_3","Race_4","BMI_2_0","Vol_WB_TIV_2_0",
                                              "Centre_0_0",YList)])


plots <- list()

for (i in 1:7) {
  
  y_var <- YList[i]          # network-specific redundancy/synergy
  net_name <- Yeo7Names[i]
  
  ## -----------------------------
  ## Step 1: LMM ¡ú residuals
  ## -----------------------------
  lmm <- lmer(
    as.formula(paste0(
      y_var,
      " ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + ",
      "Race_1 + Race_2 + Race_3 + Race_4 + ",
      "Vol_WB_TIV_2_0 + (1|Centre_0_0)"
    )),
    data = UKB_data_MRItmp,
    REML = TRUE
  )
  
  UKB_data_MRItmp$residual <- resid(lmm)
  
  # ÄâºÏ´ø½»»¥ÏîµÄÏßÐÔÄ£ÐÍ
  interaction_model <- lm(
    Fluid_intelligence_score_2_0 ~ residual * AgeAttend_2_0,
    data = UKB_data_MRItmp
  )
  
  age_mean <- mean(UKB_data_MRItmp$AgeAttend_2_0, na.rm = TRUE)
  age_sd   <- sd(UKB_data_MRItmp$AgeAttend_2_0, na.rm = TRUE)
  
  # ´´½¨ÐÂÊý¾ÝÓÃÓÚÔ¤²â£¨Èý¸öÄêÁäË®Æ½£©
  new_data <- expand.grid(
    residual = seq(min(UKB_data_MRItmp$residual, na.rm = TRUE),
                   max(UKB_data_MRItmp$residual, na.rm = TRUE),
                   length.out = 100),
    AgeAttend_2_0 = c(age_mean - age_sd, age_mean, age_mean + age_sd)
  ) %>%
    mutate(
      Age_Group = case_when(
        AgeAttend_2_0 == age_mean - age_sd ~ "55.93 (-1SD)",
        AgeAttend_2_0 == age_mean         ~ "63.55 (Mean)",
        AgeAttend_2_0 == age_mean + age_sd ~ "71.17 (+1SD)"
      ),
      Age_Group = factor(Age_Group, levels = c("55.93 (-1SD)", "63.55 (Mean)", "71.17 (+1SD)"))
    )
  
  pred <- predict(interaction_model, newdata = new_data, interval = "confidence")
  new_data <- bind_cols(new_data, as.data.frame(pred))
  
  p <- ggplot() +
    geom_line(
      data = new_data,
      aes(x = residual, y = fit, linetype = Age_Group),
      color = Yeo7Color[i],     # ? Ö±½ÓÖ¸¶¨ hex ÑÕÉ«£¬¾ø¶ÔÓÐÐ§£¡
      size = 1,
      inherit.aes = FALSE
    ) +
    geom_ribbon(
      data = new_data,
      aes(x = residual, ymin = lwr, ymax = upr, group = Age_Group,fill = "gray"),
      alpha = 0.2, inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    scale_linetype_manual(
      values = c("55.93 (-1SD)" = "solid",
                 "63.55 (Mean)" = "dashed",
                 "71.17 (+1SD)" = "dotted"),
      labels = c("55.93 (-1SD)", "63.55 (Mean)", "71.17 (+1SD)")
    ) +
    labs(
      x = "Residual",
      y = "Fluid Intelligence",
      linetype = "Age",
      title = "Moderation Effect"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  plots[[i]] <- p
}


# library(patchwork)
# 1700 280
wrap_plots(plots, nrow = 1)
