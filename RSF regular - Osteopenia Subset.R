library(survival)
library(survminer)
library(dplyr)
library(splines)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pacman, rio) 
library(rcompanion)
library(car)
library(skimr)
library(mice)
library(VIM)
library(devtools)
library(MAMI)
library(glmnet)
library(plotmo)
library(Hmisc)
library(rms)
library(pec)
library(riskRegression)
library(mlr)
library(randomForestSRC)
library(ggRandomForests)

# Import dataset 
MOF <- import("C:/Users/BCS de Vries/Documents/Studie/Afstuderen/Data/Cox_regression_MOF_2020_final_subset_11_02.csv")
head(MOF)
attach(MOF)
View(MOF)

# Factor categorical variables
MOF$Gender <- factor(Gender,
                     levels = c(0,1),
                     labels = c("male","female"))
MOF$pos_family_history <- factor(pos_family_history,
                                 levels = c(0,1),
                                 labels = c("no","yes"))
MOF$prior_falls <- factor(prior_falls,
                          levels = c(0,1),
                          labels = c("no","yes"))
MOF$current_vertebral_fracture<- factor(current_vertebral_fracture,
                                        levels = c(0,1),
                                        labels = c("no","yes"))
MOF$bedridden <- factor(bedridden,
                        levels = c(0,1),
                        labels = c("no","yes"))
MOF$backpain <- factor(backpain,
                       levels = c(0,1),
                       labels = c("no","yes"))
MOF$COCP <- factor(COCP,
                   levels = c(0,1),
                   labels = c("no","yes"))
MOF$Breastfeeding <- factor(Breastfeeding,
                            levels = c(0,1),
                            labels = c("no","yes"))
MOF$Duizeligheid_DBC <- factor(Duizeligheid_DBC,
                               levels = c(0,1),
                               labels = c("no","yes"))
MOF$DM <- factor(DM,
                 levels = c(0,1),
                 labels = c("no","yes"))
MOF$CVD <- factor(CVD,
                  levels = c(0,1),
                  labels = c("no","yes"))
MOF$IBD <- factor(IBD,
                  levels = c(0,1),
                  labels = c("no","yes"))
MOF$CVA <- factor(CVA,
                  levels = c(0,1),
                  labels = c("no","yes"))
MOF$Epilepsy <- factor(Epilepsy,
                       levels = c(0,1),
                       labels = c("no","yes"))
MOF$Systemic_auto_immune <- factor(Systemic_auto_immune,
                                   levels = c(0,1),
                                   labels = c("no","yes"))
MOF$rheuma_art <- factor(rheuma_art,
                         levels = c(0,1),
                         labels = c("no","yes"))
MOF$malabsorption <- factor(malabsorption,
                            levels = c(0,1),
                            labels = c("no","yes"))
MOF$renal_insufficiency <- factor(renal_insufficiency,
                                  levels = c(0,1),
                                  labels = c("no","yes"))
MOF$corticosteroids <- factor(corticosteroids,
                              levels = c(0,1),
                              labels = c("no","yes"))
MOF$Collaps_DBC <- factor(Collaps_DBC,
                          levels = c(0,1),
                          labels = c("no","yes"))
MOF$Delier_dementie_DBC <- factor(Delier_dementie_DBC,
                                  levels = c(0,1),
                                  labels = c("no","yes"))
MOF$weigth_under_67 <- factor(weigth_under_67,
                              levels = c(0,1),
                              labels = c("no","yes"))
MOF$weight_under_60 <- factor(weight_under_60,
                              levels = c(0,1),
                              labels = c("no","yes"))
MOF$Change_length <- factor(Change_length,
                            levels = c(0,1),
                            labels = c("no","yes"))
MOF$More_6_coffee <- factor(More_6_coffee,
                            levels = c(0,1),
                            labels = c("no","yes"))
MOF$Sun_exposure <- factor(Sun_exposure,
                           levels = c(0,1),
                           labels = c("no","yes"))
MOF$Fat_fish_diet <- factor(Fat_fish_diet,
                            levels = c(0,1),
                            labels = c("no","yes"))
MOF$veg_diet <- factor(veg_diet,
                       levels = c(0,1),
                       labels = c("no","yes"))
MOF$Vitam_pills <- factor(Vitam_pills,
                          levels = c(0,1),
                          labels = c("no","yes"))
MOF$Daily_margarin <- factor(Daily_margarin,
                             levels = c(0,1),
                             labels = c("no","yes"))
MOF$Lab_GFR_decreased <- factor(Lab_GFR_decreased,
                                levels = c(0,1),
                                labels = c("GFR>=90", "GFR<90"))
#Change data for imputation
length(MOF)
MOF_impute_RF <- MOF[c(3,5:50)]
View(MOF_impute_RF)

#rfsrc with imputation by random forest
RSF_NA <- rfsrc(Surv(time_w, MOF_fracture) ~ Age + Gender + prior_falls + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + duration_menopause +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption +  renal_insufficiency + corticosteroids  + Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC +  weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + Lab_calcium + Lab_albumin + Lab_TSH + Lab_vit_d + Lab_GFR_decreased + Lab_GFR_from_90 + T_score_LWK + T_score_HIP + Weight + Length,
               data=MOF_impute_RF, ntree=500, splitrule = "bs.gradient", block.size = 10, na.action = "na.impute", nimpute = 30)
RSF_NA

matplot(RSF_NA$time.interest, 100 * t(RSF_NA$survival.oob[1:500, ]),
        main = "RSF 500 Trees", xlab = "Time", ylab = "Survival",
        type = "l", lty = 1)

dev.off()

#Make imputated dataset for cross-validation and VIMP
MOF_imputed_RF <- impute.rfsrc(Surv(time_w, MOF_fracture) ~ .,
                               data = MOF_impute_RF, nimpute=30)


#Variable importance (VIMP)
ggRSF_NA <- gg_vimp(RSF_NA, newdata = MOF_imputed_RF)
ggRSF_NA_15 <- gg_vimp(RSF_NA, newdata = MOF_imputed_RF, nvar=15)
plot(ggRSF_NA_15)
plot(ggRSF_NA)
dev.off()

# Bootstrap cross-validation for c-index with confidence interval
cv_RSF_NA <-             pec::cindex(RSF_NA,
                                    formula=Surv(time_w,MOF_fracture) ~ .,
                                    data = MOF_imputed_RF, 
                                    splitMethod = "bootcv",
                                    B = 10,
                                    confInt = TRUE,
                                    keep.matrix = TRUE)
cv_RSF_NA$BootstrapCrossValCindexMat

cv_RSF_NA$BootCvCindex