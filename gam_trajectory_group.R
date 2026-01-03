library(circlize)

draw_chord_diagram <- function(mat, clim, colset) {
  col_fun <- colorRamp2(c(-clim, 0, clim), c("blue", "white", "red"))
  
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 0, 
             track.margin = c(0.01, 0.01),  # 设置环间距
             cell.padding = c(0.02, 0.01, 0.02, 0.01))  # 设置单元格填充，调整宽度
  
  chordDiagram(
    mat,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.1),
    grid.col = colset,  # 外围的颜色，可以自行修改
    transparency = 0.5,         # 调整透明度
    col = function(x) col_fun(x),  # 连接线的颜色
    link.lwd = abs(mat)  # 根据值的绝对值设置连接线的粗细，最大值的粗细为5
    # link.lwd = 1              # 连接线的粗细
  )
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1] + 0.1,
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5)
    )
  }, bg.border = NA)
}



library(readr)
ukb_RedSynYeo7n_2_0 <- read_csv("ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/result/ukb_RedSynYeo7n_2_0.csv")
ukb_RedSynYeo7n_2_0 <- merge(merge(UKB_Basic,UKB_rfMRIcov),ukb_RedSynYeo7n_2_0)
UKB_AllVars_NumCol <- read_csv("ZJLab/UKBiobank_Project/data/table/UKB_AllVars_NumCol.csv")
UKB_AllVars_NumCol_2_0 = UKB_AllVars_NumCol[,c(1,which(str_detect(names(UKB_AllVars_NumCol),"_2_0")))]
UKB_AllVars_NumCol_2_0 = UKB_AllVars_NumCol_2_0[,c(1,which(colSums(is.na(UKB_AllVars_NumCol_2_0))>1000))]

# UKB_AllVars_NumCol_2_0 <- UKB_AllVars_NumCol_2_0[,-1*which(names(UKB_AllVars_NumCol_2_0) %in% 
#                                                              c("Birth_weight_known_2_0",
#                                                                "Part_of_a_multiple_birth_2_0",
#                                                                "Maternal_smoking_around_birth_2_0",
#                                                                "Pregnant_2_0") )]


UKB_BrainAllVars_2_0 = merge(ukb_RedSynYeo7n_2_0,UKB_AllVars_NumCol_2_0)

xList = names(UKB_AllVars_NumCol_2_0)[3:ncol(UKB_AllVars_NumCol_2_0)]
yList = names(ukb_RedSynYeo7n_2_0)[80:ncol(ukb_RedSynYeo7n_2_0)]
s = 0;
fit.synred.allvars <- data.frame(x=NA,y=NA,Estimate=NA,StdError=NA,tvalue=NA,P=NA)
for (x in xList) {
  for (y in yList) {
    formula <- as.formula(paste(y, '~', x, '+Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+',
                                'Race_1+Race_2+Race_3+Race_4+Centre_2_0+BMI_2_0', sep = ''))
    tryCatch({
      s <- s + 1
      Fit <- lm(formula, data = UKB_BrainAllVars_2_0)
      sFit <- summary(Fit)
      fit.synred.allvars[s, 1] <- x
      fit.synred.allvars[s, 2] <- y
      fit.synred.allvars[s, 3:6] <- sFit$coefficients[x, ]
      print(paste('regress: ', x, ' ~ ', y))
    }, error = function(e) {
      message(paste("Error in regression for x:", x, "and y:", y, "-", e$message))
    })
  }
}



library(mgcv)
library(gratia)
library(nlme)

# syn-red data (network-level): 7*3 = 21

# GAM: total syn & red
gam_SynYeo7n_Total <- gam(SynYeo7n_Total ~ s(AgeAttend_2_0) + s(eid, bs = "re") + 
                            Sex_0_0 + BMI_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + 
                            BloodPress_2_0 + HeadMotion_2_0 + SNR_2_0, data = ukb_RedSynYeo7n_2_0)

ggplot(ukb_RedSynYeo7n_2_0, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = SynYeo7n_Total), method = "gam", formula = y ~ s(x), se = T) +
  labs(title = "GAM", x = "Age", y = "Syn & Red In Net") +
  theme_bw()

ggplot(ukb_RedSynYeo7n_2_0, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = RedYeo7n_Total), method = "gam", formula = y ~ s(x), se = T) +
  labs(title = "GAM", x = "Age", y = "Syn & Red In Net") +
  theme_bw()


gam_SynYeo7n_In1 <- gam(SynYeo7n_In1 ~ s(AgeAttend_2_0) + s(eid, bs = "re") + 
                          Sex_0_0 + BMI_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + 
                          BloodPress_2_0 + HeadMotion_2_0 + SNR_2_0, data = ukb_RedSynYeo7n_2_0)


gam_SynYeo7n_In7 <- gam(SynYeo7n_In7 ~ s(AgeAttend_2_0) + s(eid, bs = "re") + 
                          Sex_0_0 + BMI_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + 
                          BloodPress_2_0 + HeadMotion_2_0 + SNR_2_0, data = ukb_RedSynYeo7n_2_0)


Yeo7Color = c(rgb(0.471,0.0710,0.522),
              rgb(0.275,0.510,0.706),
              rgb(0,0.463,0.0550),
              rgb(0.769,0.224,0.976),
              rgb(0.863,0.973,0.639),
              rgb(0.902,0.576,0.129),
              rgb(0.804,0.239,0.306))

ggplot(ukb_RedSynYeo7n_2_0, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = SynYeo7n_In1), method = "gam", formula = y ~ s(x), color = Yeo7Color[1], se = T) +
  geom_smooth(aes(y = SynYeo7n_In2), method = "gam", formula = y ~ s(x), color = Yeo7Color[2], se = T) +
  geom_smooth(aes(y = SynYeo7n_In3), method = "gam", formula = y ~ s(x), color = Yeo7Color[3], se = T) +
  geom_smooth(aes(y = SynYeo7n_In4), method = "gam", formula = y ~ s(x), color = Yeo7Color[4], se = T) +
  geom_smooth(aes(y = SynYeo7n_In5), method = "gam", formula = y ~ s(x), color = Yeo7Color[5], se = T) +
  geom_smooth(aes(y = SynYeo7n_In6), method = "gam", formula = y ~ s(x), color = Yeo7Color[6], se = T) +
  geom_smooth(aes(y = SynYeo7n_In7), method = "gam", formula = y ~ s(x), color = Yeo7Color[7], se = T) +
  labs(title = "GAM", x = "Age", y = "Syn In Net") +
  theme_bw()


ggplot(ukb_RedSynYeo7n_2_0, aes(x = AgeAttend_2_0)) +
  geom_smooth(aes(y = RedYeo7n_In1), method = "gam", formula = y ~ s(x), color = Yeo7Color[1], se = T) +
  geom_smooth(aes(y = RedYeo7n_In2), method = "gam", formula = y ~ s(x), color = Yeo7Color[2], se = T) +
  geom_smooth(aes(y = RedYeo7n_In3), method = "gam", formula = y ~ s(x), color = Yeo7Color[3], se = T) +
  geom_smooth(aes(y = RedYeo7n_In4), method = "gam", formula = y ~ s(x), color = Yeo7Color[4], se = T) +
  geom_smooth(aes(y = RedYeo7n_In5), method = "gam", formula = y ~ s(x), color = Yeo7Color[5], se = T) +
  geom_smooth(aes(y = RedYeo7n_In6), method = "gam", formula = y ~ s(x), color = Yeo7Color[6], se = T) +
  geom_smooth(aes(y = RedYeo7n_In7), method = "gam", formula = y ~ s(x), color = Yeo7Color[7], se = T) +
  labs(title = "GAM", x = "Age", y = "Red In Net") +
  theme_bw()



summary(gam_SynYeo7n_In1)

# group-specific modeling
gam_group_syn <- gam(syn_total ~ cogn + s(age, by = cogn) + s(subject, bs = "re") + sex + BMI + race + bloodpress, data = data)

summary(gam_group_syn)

# 可视化差异曲线
# 使用gratia包来评估不同逆境轨迹之间的差异
# 例如：无逆境与低逆境之间的差异
diff_curve <- difference_smooths(gam_group_syn, smooth_var = "age", by = "cogn")

# 绘制差异曲线
plot(diff_curve)

# 进行LME模型以评估逆境组间SC-FC变化的显著性差异
lme_model <- lme(gam_group_syn ~ cogn * age + sex, random = ~1|subject, data = data)

# 打印模型摘要
summary(lme_model)

# 验证加速发展的假设
# 假设你的数据框名为 data_epigenetic，包含以下列：epigenetic_age_acceleration, adversity, PC1, PC2, PC3

# 回归分析比较逆境组间表观遗传年龄加速的显著性差异
regression_model <- lm(epigenetic_age_acceleration ~ adversity + PC1 + PC2 + PC3, data = data_epigenetic)

# 打印模型摘要
summary(regression_model)



# 加载所需的包
library(segmented)
library(ggplot2)

# 假设你的数据框名为data，包含以下列：x（自变量）和y（因变量）

# 首先拟合一个线性模型
linear_model <- lm(y ~ x, data = data)

# 使用segmented包找到最佳拐点并拟合分段回归模型
segmented_model <- segmented(linear_model, seg.Z = ~x)

# 打印模型摘要
summary(segmented_model)

# 获取拐点
breakpoint <- segmented_model$psi[2]

# 打印拐点
print(breakpoint)

# 可视化分段回归模型
data$fit <- fitted(segmented_model)
ggplot(data, aes(x = x, y = y)) +
  geom_point() +
  geom_line(aes(y = fit), color = "blue") +
  geom_vline(xintercept = breakpoint, linetype = "dashed", color = "red") +
  labs(
    title = "Segmented Regression with Breakpoint",
    x = "X",
    y = "Y"
  ) +
  theme_minimal()




mat <- matrix(c(
  0, 2, -3, 4, -5, 6, -1,
  2, 0, 1, -2, 3, -4, 5,
  -3, 1, 0, 2, -1, 4, -2,
  4, -2, 2, 0, 3, -1, 1,
  -5, 3, -1, 3, 0, 2, -4,
  6, -4, 4, -1, 2, 0, 3,
  -1, 5, -2, 1, -4, 3, 0
), nrow = 7, byrow = TRUE)
mat[abs(mat)<4] <-0

Yeo7Name <- c("VIS","SMN","DAN","VAN","LIM","FPN","DMN")
rownames(mat) <- Yeo7Name
colnames(mat) <- Yeo7Name
draw_chord_diagram(mat, 7, Yeo7Color)

