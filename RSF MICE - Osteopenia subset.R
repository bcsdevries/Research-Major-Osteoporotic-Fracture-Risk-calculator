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

#rfsrc on single imputed datasets by MICE -> datasets created in 'Cox regression - Osteopenia subset'
RSF_1 <- rfsrc(Surv(time_m, MOF_fracture) ~ Age + Gender + prior_falls + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + duration_menopause +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption + renal_insufficiency + corticosteroids +  Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC + weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + Lab_calcium + Lab_albumin + Lab_TSH + Lab_vit_d + Lab_GFR_decreased + Lab_GFR_from_90 + T_score_LWK + T_score_HIP + Weight + Length,
               data=complete1, ntree=500, splitrule = "bs.gradient", block.size = 10, importance = TRUE)
RSF_2 <- rfsrc(Surv(time_m, MOF_fracture) ~ Age + Gender + prior_falls + current_vertebral_fracture + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + duration_menopause +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption + renal_insufficiency + corticosteroids +  Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC + weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + Lab_calcium + Lab_albumin + Lab_TSH + Lab_vit_d + Lab_GFR_decreased + Lab_GFR_from_90 + T_score_LWK + T_score_HIP + Weight + Length,
               data=complete2, ntree=500, splitrule = "bs.gradient", block.size = 10, importance = TRUE)
RSF_3 <- rfsrc(Surv(time_m, MOF_fracture) ~ Age + Gender + prior_falls + current_vertebral_fracture + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + duration_menopause +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption + renal_insufficiency + corticosteroids +  Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC + weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + Lab_calcium + Lab_albumin + Lab_TSH + Lab_vit_d + Lab_GFR_decreased + Lab_GFR_from_90 + T_score_LWK + T_score_HIP + Weight + Length,
               data=complete3, ntree=500, splitrule = "bs.gradient", block.size = 10, importance = TRUE)
RSF_4 <- rfsrc(Surv(time_m, MOF_fracture) ~ Age + Gender + prior_falls + current_vertebral_fracture + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + duration_menopause +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption + renal_insufficiency + corticosteroids +  Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC + weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + Lab_calcium + Lab_albumin + Lab_TSH + Lab_vit_d + Lab_GFR_decreased + Lab_GFR_from_90 + T_score_LWK + T_score_HIP + Weight + Length,
               data=complete4, ntree=500, splitrule = "bs.gradient", block.size = 10, importance = TRUE)
RSF_5 <- rfsrc(Surv(time_m, MOF_fracture) ~ Age + Gender + prior_falls + current_vertebral_fracture + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + duration_menopause +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption + renal_insufficiency + corticosteroids +  Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC + weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + Lab_calcium + Lab_albumin + Lab_TSH + Lab_vit_d + Lab_GFR_decreased + Lab_GFR_from_90 + T_score_LWK + T_score_HIP + Weight + Length,
               data=complete1, ntree=500, splitrule = "bs.gradient", block.size = 10, importance = TRUE)

# Loop with RSF on every single dataset, using 10-fold cross validation. 
matrix_c_index_total_RSF <- as.data.frame(matrix(c(1,2,3,4,5,6,7,8,9,10)))
c_index_RSF <- as.data.frame(matrix(c("training", "test" )))

for (i in 1:30){ 
  data4 <- as.data.frame(complete(MOF_imputed_good, i))
  RSF_function <- rfsrc(Surv(time_m, MOF_fracture) ~ Age + Gender + prior_falls + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + duration_menopause +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption + renal_insufficiency + corticosteroids + Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC + weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + Lab_calcium + Lab_albumin + Lab_TSH + Lab_vit_d + Lab_GFR_decreased + Lab_GFR_from_90 + T_score_LWK + T_score_HIP + Weight + Length,
                 data=data4, ntree=500, splitrule = "bs.gradient")
    cv_RSF_function <-     pec::cindex(RSF_function,
                                     formula=Surv(time_w,MOF_fracture) ~ .,
                                     data = data4, 
                                     splitMethod = "bootcv",
                                     B = 10,
                                     confInt = TRUE,
                                     keep.matrix = TRUE,)
  matrix_c_index_RSF <- as.data.frame(cv_RSF_function$BootstrapCrossValCindexMat)
  matrix_c_index_total_RSF <- cbind(matrix_c_index_total_RSF, matrix_c_index_RSF)
  c_index_data4_RSF <- as.data.frame(matrix(c(cv_RSF_function$AppCindex, cv_RSF_function$BootCvCindex)))
  c_index_RSF <- cbind(c_index_RSF, c_index_data4_RSF)
  cat("RSF nummer", i)
}

colnames(c_index_RSF) <- c("test/training", "dataframe 1", "dataframe 2", "dataframe3", "dataframe4", "dataframe5", "dataframe6", "dataframe7", "dataframe8", "dataframe9", "dataframe10", "dataframe11", "dataframe12", "dataframe13", "dataframe14", "dataframe15", "dataframe16", "dataframe17", "dataframe18", "dataframe19", "dataframe20", "dataframe21", "dataframe22", "dataframe23", "dataframe24", "dataframe25", "dataframe26", "dataframe27", "dataframe28", "dataframe29", "dataframe30")
View(c_index_RSF)

c_index_adj_RSF <- c_index_RSF[,2:31]
c_index_matrix_RSF <- data.matrix(c_index_adj_RSF)
c_index_mean_RSF <- rowMeans(c_index_matrix_RSF)

cat("concordance index training", c_index_mean_RSF[1],"\n",
    "concordance index test", c_index_mean_RSF[2])

# confidence interval using Rubin's Rules
colnames(matrix_c_index_total_RSF) <- c("#iteration", "dataframe 1", "dataframe 2", "dataframe3", "dataframe4", "dataframe5", "dataframe6", "dataframe7", "dataframe8", "dataframe9", "dataframe10", "dataframe11", "dataframe12", "dataframe13", "dataframe14", "dataframe15", "dataframe16", "dataframe17", "dataframe18", "dataframe19", "dataframe20", "dataframe21", "dataframe22", "dataframe23", "dataframe24", "dataframe25", "dataframe26", "dataframe27", "dataframe28", "dataframe29", "dataframe30")
mean <- 0
variance <- 0
for (i in 2:31){
  mean[i] <- mean(matrix_c_index_total_RSF[,i])
  variance[i] <- var(matrix_c_index_total_RSF[,i])
}
mean_good <- mean[2:31]
variance_good <- variance[2:31]
pooled<- pool.scalar(Q=mean_good, U=variance_good, n=1770)
pooled$qbar
Confidence_interval_low = pooled$qbar-1.96*sqrt(pooled$t)
Confidence_interval_high = pooled$qbar+1.96*sqrt(pooled$t)
CI <- c(Confidence_interval_low,Confidence_interval_high)
print(CI)

# Variable importance
ggRSF_1 <- gg_vimp(RSF_1, nvar=15)
plot(ggRSF_1)