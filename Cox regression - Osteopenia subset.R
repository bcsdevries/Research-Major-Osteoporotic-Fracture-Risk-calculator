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

#add Hazard function
hazard <- nelsonaalen(MOF, time_w, MOF_fracture)
hazard

#Delete variable rec_id, add hazard
ncol(MOF)
MOF_impute <- MOF[2:50]
View(MOF_impute)
MOF_impute <- cbind(MOF_impute, hazard) # add hazard
View(MOF_impute)

# Add interaction terms for imputation with MICE

  # Interaction Age (68.53) with hours per week active (25.87), mean in (). 
Meth <- make.method(MOF_impute)
expr_age.hours <- expression((Age - 68.53) * (Hours_active - 25.87))
MOF_impute$Age.hours <- with(MOF_impute, eval(expr_age.hours))
Meth["Age.hours"] <- paste("~I(", expr_age.hours, ")", sep = "")

  # Interaction Age (68.53) with T_score_HIP (-1.45), mean in (). 
expr_age.T_score_HIP <- expression((Age - 68.53) * (T_score_HIP - (-1.45)))
MOF_impute$Age.T_score_HIP <- with(MOF_impute, eval(expr_age.T_score_HIP))
Meth["Age.T_score_HIP"] <- paste("~I(", expr_age.T_score_HIP, ")", sep = "")
  
  # Interaction Age (68.53) with T_score_LWK (-1.50), mean in (). 
expr_age.T_score_LWK <- expression((Age - 68.53) * (T_score_LWK - (-1.50)))
MOF_impute$Age.T_score_LWK <- with(MOF_impute, eval(expr_age.T_score_LWK))
Meth["Age.T_score_LWK"] <- paste("~I(", expr_age.T_score_LWK, ")", sep = "")

  # Interaction Age (68.53) with Gender, gemiddelde in (). 
expr_age.Gender <- expression((Age - 68.53) * Gender_interaction)
MOF_impute$Age.Gender <- with(MOF_impute, eval(expr_age.Gender))
Meth["Age.Gender"] <- paste("~I(", expr_age.Gender, ")", sep = "")

  # Interaction Gender with hours per week active (25.87), mean in (). 
expr_Gender.hours <- expression((Gender_interaction) * (Hours_active - 25.87))
MOF_impute$Gender.hours <- with(MOF_impute, eval(expr_Gender.hours))
Meth["Gender.hours"] <- paste("~I(", expr_Gender.hours, ")", sep = "")

  # Interaction Gender with T_score_HIP (-1.45), mean in (). 
expr_Gender.T_score_HIP <- expression((Gender_interaction) * (T_score_HIP - (-1.45)))
MOF_impute$Gender.T_score_HIP <- with(MOF_impute, eval(expr_Gender.T_score_HIP))
Meth["Gender.T_score_HIP"] <- paste("~I(", expr_Gender.T_score_HIP, ")", sep = "")

  # Interaction Gender with T_score_LWK (-1.50), mean in (). 
expr_Gender.T_score_LWK <- expression((Gender_interaction) * (T_score_LWK - (-1.50)))
MOF_impute$Gender.T_score_LWK <- with(MOF_impute, eval(expr_Gender.T_score_LWK))
Meth["Gender.T_score_LWK"] <- paste("~I(", expr_Gender.T_score_LWK, ")", sep = "")

#Prepare prediction matrix of MICE
Pred <- quickpred(MOF_impute, include = c("hazard", "Gender", "Age","MOF_fracture"),
                  exclude = c( "time_m","time_w"),
                  minpuc = 0.5, mincor = 0.1)
Pred[c("duration_menopause", "Lab_GFR_from_90","children","Hours_active"), "Age.hours"] <- 0
Pred[c("duration_menopause", "Lab_GFR_from_90","children","T_score_HIP"), "Age.T_score_HIP"] <- 0
Pred[c("duration_menopause", "Lab_GFR_from_90","children","T_score_LWK"), "Age.T_score_LWK"] <- 0
Pred[c("duration_menopause", "Lab_GFR_from_90","children","Gender"), "Age.Gender"] <- 0
Pred[c("duration_menopause", "Lab_GFR_from_90","children","Hours_active"), "Gender.hours"] <- 0
Pred[c("duration_menopause", "Lab_GFR_from_90","children","T_score_HIP"), "Gender.T_score_HIP"] <- 0
Pred[c("duration_menopause", "Lab_GFR_from_90","children","T_score_LWK"), "Gender.T_score_LWK"] <- 0
length(MOF_impute)
Pred[1:57, "Gender_interaction"] <- 0

# mean number of predictors
rowSums(Pred)
select <- which(rowSums(Pred)>0)
View(select)
select_no_interaction <- select[1:29]
mean(select_no_interaction)

#Impute with conditions and correct predictormatrix
MOF_imputed_good <- mice(MOF_impute, pred = Pred, meth = Meth,
                         m=30, maxit = 50, seed=2795, print = F)
MOF_imputed_good
densityplot(MOF_imputed_good)

# complete datasets seperately
complete1 <- complete(MOF_imputed_good,1)
complete2 <- complete(MOF_imputed_good,2)
complete3 <- complete(MOF_imputed_good,3)
complete4 <- complete(MOF_imputed_good,4)
complete5 <- complete(MOF_imputed_good,5)
complete6 <- complete(MOF_imputed_good,6)
complete7 <- complete(MOF_imputed_good,7)
complete8 <- complete(MOF_imputed_good,8)
complete9 <- complete(MOF_imputed_good,9)
complete10 <- complete(MOF_imputed_good,10)
complete11 <- complete(MOF_imputed_good,11)
complete12 <- complete(MOF_imputed_good,12)
complete13 <- complete(MOF_imputed_good,13)
complete14 <- complete(MOF_imputed_good,14)
complete15 <- complete(MOF_imputed_good,15)
complete16 <- complete(MOF_imputed_good,16)
complete17 <- complete(MOF_imputed_good,17)
complete18 <- complete(MOF_imputed_good,18)
complete19 <- complete(MOF_imputed_good,19)
complete20 <- complete(MOF_imputed_good,20)
complete21 <- complete(MOF_imputed_good,21)
complete22 <- complete(MOF_imputed_good,22)
complete23 <- complete(MOF_imputed_good,23)
complete24 <- complete(MOF_imputed_good,24)
complete25 <- complete(MOF_imputed_good,25)
complete26 <- complete(MOF_imputed_good,26)
complete27 <- complete(MOF_imputed_good,27)
complete28 <- complete(MOF_imputed_good,28)
complete29 <- complete(MOF_imputed_good,29)
complete30 <- complete(MOF_imputed_good,30)

# complete a stacked dataset with weight 1/w. 
MOF_stacked <- complete1
for (i in 2:30){
  MOF_stacked <- rbind(MOF_stacked, complete(MOF_imputed_good, i))}
MOF_stacked$w <- 1/30

#check if there are no values imputed which are not possible
min(MOF_stacked$Lab_bse)
min(MOF_stacked$Lab_calcium)
min(MOF_stacked$Lab_albumin)
min(MOF_stacked$Lab_TSH)
min(MOF_stacked$Lab_vit_d)
min(MOF_stacked$Lab_GFR_from_90)
min(MOF_stacked$Hours_active)
min(MOF_stacked$Total_daily_calcium)
min(MOF_stacked$children)

# Check assumption of proportional hazard for complete model
cox_function_pha <- coxph(Surv(time_w,MOF_fracture) ~ Age + Gender + prior_falls + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + ns(duration_menopause, 3) +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption +  renal_insufficiency + corticosteroids + Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC +  weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + ns(Lab_calcium, 4) + Lab_albumin + Lab_TSH + ns(Lab_vit_d, 3) + Lab_GFR_decreased + ns(Lab_GFR_from_90,3) + T_score_LWK + T_score_HIP + Weight + Length +  Hours_active:Age +  T_score_HIP:Age + T_score_LWK:Age + Gender:Age + Hours_active:Gender + T_score_HIP:Gender + T_score_LWK:Gender,
                          data = MOF_stacked,
                          weights = w,
                          x = TRUE,
                          y = TRUE)
summary(cox_function_pha)
print(cox.zph(cox_function_pha))

# Check assumption of proportional hazard for final model
cox_function_final_model <- coxph(Surv(time_w, MOF_fracture)~ Age + Gender + prior_falls + ns(duration_menopause,3) + T_score_HIP + Epilepsy + T_score_HIP:Age,
                                  data = MOF_stacked,
                                  weights = w,
                                  x = TRUE,
                                  y = TRUE)
print(cox.zph(cox_function_final_model))

# Model for LASSO and fit with glmnet
cox_function_all <- coxph(Surv(time_w, MOF_fracture) ~ Age + Gender + prior_falls + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + ns(duration_menopause, 3) +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption  +  renal_insufficiency + corticosteroids  + Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC +  weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + ns(Lab_calcium, 4) + Lab_albumin + Lab_TSH + ns(Lab_vit_d, 3) + Lab_GFR_decreased + ns(Lab_GFR_from_90,3) + T_score_LWK + T_score_HIP + Weight + Length + Hours_active:Age + T_score_HIP:Age + T_score_LWK:Age + Gender:Age + Hours_active:Gender + T_score_HIP:Gender + T_score_LWK:Gender,                          
                          data = MOF_stacked,
                          weights = w,
                          x = TRUE,
                          y = TRUE)
fit = glmnet(
  x = cox_function_all$x,
  y = cox_function_all$y,
  alpha = 1, family = "cox")

# Compute Coxph with LASSO on every dataset to get mean penalty factor lambda (both 1 SE and min) 
coef_min_total <- 0
coef_1se_total <- 0
for (i in 1:30){ 
  data1 <- as.data.frame(complete(MOF_imputed_good, i))
  cox_function_LASSO <- coxph(Surv(time_w,MOF_fracture) ~ Age + Gender + prior_falls + pos_family_history + bedridden + backpain + children + COCP + Breastfeeding + ns(duration_menopause, 3) +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption +  renal_insufficiency + corticosteroids + Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC +  weigth_under_67 + weight_under_60 + Change_length + Hours_active + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + Lab_bse + ns(Lab_calcium, 4) + Lab_albumin + Lab_TSH + ns(Lab_vit_d, 3) + Lab_GFR_decreased + ns(Lab_GFR_from_90,3) + T_score_LWK + T_score_HIP + Weight + Length  + Hours_active:Age  + T_score_HIP:Age + T_score_LWK:Age + Gender:Age + Hours_active:Gender + T_score_HIP:Gender + T_score_LWK:Gender,
                              data = data1,
                              x = TRUE,
                              y = TRUE)
  coxlasso1 <- glmnet(
    x = cox_function_LASSO$x,
    y = cox_function_LASSO$y,
    alpha=1, family = 'cox')
  cvfit1 <- cv.glmnet(
    x = cox_function_LASSO$x,
    y = cox_function_LASSO$y, family = "cox")
  coef_1se_total <- c(coef_1se_total, cvfit1$lambda.1se)
  coef_min_total <- c(coef_min_total, cvfit1$lambda.min)
}

coef_1se_total <- coef_1se_total[2:31]
mean_LASSO_1se <- mean(coef_1se_total)
plot(cvfit1)

coef_select_1se <- coef(fit, s=mean_LASSO_1se)
select_1se <- which(or(coef_select_1se>0,coef_select_1se<0))
coef_select_1se[select_1se, 1]

coef_min_total <- coef_min_total[2:31]
mean_LASSO_min <- mean(coef_min_total)

coef_select_min <- coef(fit, s=mean_LASSO_min)
select_min <- which(or(coef_select_min>0,coef_select_min<0))
coef_select_min[select_min, 1]

#Fit final Cox model on all 30 datasets and pool for lambda 1se
fit.Cox <- with(data=MOF_imputed_good, exp=coxph(Surv(time_w,MOF_fracture) ~ Age))
Cox.pool <- summary(pool(fit.Cox))
Cox.pool

pool.HR <- exp(cbind(Cox.pool[,1], (Cox.pool[,1]-1.96*(Cox.pool[,2])), 
                     (Cox.pool[,1]+1.96*(Cox.pool[,2]))))
colnames(pool.HR) <- (c("HR", "95% LO", "95% UP"))
pool.HR

outcome_Cox <- cbind(Cox.pool, pool.HR)
outcome_Cox <- round(outcome_Cox,3)
print(outcome_Cox)

#Fit final Cox model on all 30 datasets and pool for lambda min
fit.Cox_min <- with(data=MOF_imputed_good, exp=coxph(Surv(time_w,MOF_fracture) ~ Age + Gender +  prior_falls + CVD + Delier_dementie_DBC  + Change_length + Hours_active +  Fat_fish_diet + Lab_GFR_decreased + ns(Lab_GFR_from_90, 3) + Age:T_score_HIP + T_score_HIP))
Cox.pool_min <- summary(pool(fit.Cox_min))
Cox.pool_min

pool.HR_min <- exp(cbind(Cox.pool_min[,1], (Cox.pool_min[,1]-1.96*(Cox.pool_min[,2])), 
                         (Cox.pool_min[,1]+1.96*(Cox.pool_min[,2]))))
colnames(pool.HR_min) <- (c("HR", "95% LO", "95% UP"))
pool.HR_min

outcome_Cox_min <- cbind(Cox.pool_min, pool.HR_min)
outcome_Cox_min <- round(outcome_Cox_min,3)
print(outcome_Cox_min)

# Cross validation loop for 1se
matrix_c_index_total <- as.data.frame(matrix(c(1,2,3,4,5,6,7,8,9,10)))
c_index <- as.data.frame(matrix(c("training", "test" )))

for (i in 1:30){ 
  data5 <- as.data.frame(complete(MOF_imputed_good, i))
  dd <- datadist(data5)
  options(datadist="dd")  
  cv_cox_function <- cph(Surv(time_w,MOF_fracture) ~ Age,
                         data = data5,
                         surv = TRUE,
                         x = TRUE,
                         y = TRUE)
  cv_cox_function <-     pec::cindex(cv_cox_function,
                                     formula=Surv(time_w,MOF_fracture) ~ .,
                                     data = data5, 
                                     splitMethod = "bootcv",
                                     B = 10,
                                     confInt = TRUE,
                                     keep.matrix = TRUE,)
  matrix_c_index <- as.data.frame(cv_cox_function$BootstrapCrossValCindexMat)
  matrix_c_index_total <- cbind(matrix_c_index_total, matrix_c_index)
  c_index_data2 <- as.data.frame(matrix(c(cv_cox_function$AppCindex, cv_cox_function$BootCvCindex)))
  c_index <- cbind(c_index, c_index_data2)
}

colnames(c_index) <- c("test/training", "dataframe 1", "dataframe 2", "dataframe3", "dataframe4", "dataframe5", "dataframe6", "dataframe7", "dataframe8", "dataframe9", "dataframe10", "dataframe11", "dataframe12", "dataframe13", "dataframe14", "dataframe15", "dataframe16", "dataframe17", "dataframe18", "dataframe19", "dataframe20", "dataframe21", "dataframe22", "dataframe23", "dataframe24", "dataframe25", "dataframe26", "dataframe27", "dataframe28", "dataframe29", "dataframe30")
View(c_index)

c_index_adj <- c_index[,2:31]
c_index_matrix <- data.matrix(c_index_adj)
c_index_mean <- rowMeans(c_index_matrix)
training <- "concordance index training"
test <- "concordance index test"
cat("concordance index training", c_index_mean[1],"\n",
    "concordance index test", c_index_mean[2])


# confidence interval using Rubin's Rules
colnames(matrix_c_index_total) <- c("#iteration", "dataframe 1", "dataframe 2", "dataframe3", "dataframe4", "dataframe5", "dataframe6", "dataframe7", "dataframe8", "dataframe9", "dataframe10", "dataframe11", "dataframe12", "dataframe13", "dataframe14", "dataframe15", "dataframe16", "dataframe17", "dataframe18", "dataframe19", "dataframe20", "dataframe21", "dataframe22", "dataframe23", "dataframe24", "dataframe25", "dataframe26", "dataframe27", "dataframe28", "dataframe29", "dataframe30")
mean <- 0
variance <- 0
for (i in 2:31){
  mean[i] <- mean(matrix_c_index_total[,i])
  variance[i] <- var(matrix_c_index_total[,i])
}
mean_good <- mean[2:31]
variance_good <- variance[2:31]
pooled<- pool.scalar(Q=mean_good, U=variance_good, n=1770)
pooled$qbar
Confidence_interval_low = pooled$qbar-1.96*sqrt(pooled$t)
Confidence_interval_high = pooled$qbar+1.96*sqrt(pooled$t)
CI <- c(Confidence_interval_low,Confidence_interval_high)
print(CI)

# Cross validation loop for lambda min
matrix_c_index_total_min <- as.data.frame(matrix(c(1,2,3,4,5,6,7,8,9,10)))
c_index_min <- as.data.frame(matrix(c("training", "test" )))

for (i in 1:30){ 
  data6 <- as.data.frame(complete(MOF_imputed_good, i))
  dd <- datadist(data6)
  options(datadist="dd")
  cv_cox_function_min <- cph(Surv(time_w,MOF_fracture) ~ Age + Gender +  prior_falls + CVD + Delier_dementie_DBC  + Change_length + Hours_active +  Fat_fish_diet + Lab_GFR_decreased + rcs(Lab_GFR_from_90, 3) + Age*T_score_HIP + T_score_HIP,
                             data = data6,
                             surv = TRUE,
                             x = TRUE,
                             y = TRUE)
  cv_cox_function_min <-     pec::cindex(cv_cox_function_min,
                                         formula=Surv(time_w,MOF_fracture) ~ .,
                                         data = data6, 
                                         splitMethod = "bootcv",
                                         B = 10,
                                         confInt = TRUE,
                                         keep.matrix = TRUE,)
  matrix_c_index_min <- as.data.frame(cv_cox_function_min$BootstrapCrossValCindexMat)
  matrix_c_index_total_min <- cbind(matrix_c_index_total_min, matrix_c_index_min)
  c_index_data3_min <- as.data.frame(matrix(c(cv_cox_function_min$AppCindex, cv_cox_function_min$BootCvCindex)))
  c_index_min <- cbind(c_index_min, c_index_data3_min)
}

colnames(c_index_min) <- c("test/training", "dataframe 1", "dataframe 2", "dataframe3", "dataframe4", "dataframe5", "dataframe6", "dataframe7", "dataframe8", "dataframe9", "dataframe10", "dataframe11", "dataframe12", "dataframe13", "dataframe14", "dataframe15", "dataframe16", "dataframe17", "dataframe18", "dataframe19", "dataframe20", "dataframe21", "dataframe22", "dataframe23", "dataframe24", "dataframe25", "dataframe26", "dataframe27", "dataframe28", "dataframe29", "dataframe30")
View(c_index_min)

colnames(matrix_c_index_total_min) <- c("#iteration", "dataframe 1", "dataframe 2", "dataframe3", "dataframe4", "dataframe5", "dataframe6", "dataframe7", "dataframe8", "dataframe9", "dataframe10", "dataframe11", "dataframe12", "dataframe13", "dataframe14", "dataframe15", "dataframe16", "dataframe17", "dataframe18", "dataframe19", "dataframe20", "dataframe21", "dataframe22", "dataframe23", "dataframe24", "dataframe25", "dataframe26", "dataframe27", "dataframe28", "dataframe29", "dataframe30")
View(matrix_c_index_total_min)

c_index_adj_min <- c_index_min[,2:31]
c_index_matrix_min <- data.matrix(c_index_adj_min)
c_index_mean_min <- rowMeans(c_index_matrix_min)

cat("concordance index training", c_index_mean_min[1],"\n",
    "concordance index test", c_index_mean_min[2])

#confidence interval using Rubin's Rules
mean <- 0
variance <- 0
for (i in 2:31){
  mean[i] <- mean(matrix_c_index_total_min[,i])
  variance[i] <- var(matrix_c_index_total_min[,i])
}
mean_good <- mean[2:31]
variance_good <- variance[2:31]
pooled<- pool.scalar(Q=mean_good, U=variance_good, n=1770)
pooled$qbar
Confidence_interval_low = pooled$qbar-1.96*sqrt(pooled$t)
Confidence_interval_high = pooled$qbar+1.96*sqrt(pooled$t)
CI <- c(Confidence_interval_low,Confidence_interval_high)
print(CI)

# Shiny Tool

library(rsconnect)
rsconnect::setAccountInfo(name='majorosteoporoticfracture', token='F15BD5A96C6A91A27B59E1EEE2410AEC', secret='EWSKAkJRir8ieRD2YBGYcRtkSxN2UMOuwvc2b/6e')
library(shiny)
library(scales)

ui <- fluidPage(
  titlePanel("Major Osteoporotic Fracture Risk Calculator"),
  wellPanel(fluidRow(
    column(6, numericInput("age", "Age:", 50, min = 50, max = 110)),
    column(6, selectInput("gender", "Gender:", c("male","female"), selected = "female", multiple=FALSE)),
  ),
  fluidRow(
    column(6, selectInput("prior_falls", "Prior Falls:", c("no","yes"), selected="no")),
    column(6, selectInput("cvd", "Cardiovascular disease:", c("no","yes"), selected="no")),
  ),
  fluidRow(
    column(6, selectInput("delier_dementie", "Delirium or dementia:", c("no","yes"), selected="no")),
    column(6, selectInput("change_length", "Change in length in recent years:", c("no","yes"), selected="no")),
  ),
  fluidRow(
    column(6, numericInput("hours_active", "Hours of moderate activity per week:", 0, min=0, max=168)),
    column(6, selectInput("fat_fish", "Dietary use of fat fish (at least twice a week):", c("no","yes"), selected="no")),
  ),
  fluidRow(
    column(6, numericInput("eGFR", "Last known eGFR:  (if > 90 fill in 90)", 90, min=0, max=90)),
    column(6, numericInput("T_score_hip", "T-score of the hip:",-1.5, min=-5, max = 5)),
  ),
  fluidRow(
    column(3, actionButton("go", "Predict risk"))
  )),
  fluidRow(
    column(12, verbatimTextOutput("text"))
  ),
  wellPanel(
    verbatimTextOutput("text_age"),
  )
)
server <- function(input, output) {
  newdata <- data.frame(Age=70, Gender = "female", prior_falls = "no", CVD = "yes", Delier_dementie_DBC = "no", Change_length="yes", Hours_active = 7, Fat_fish_diet = "yes", Lab_GFR_decreased = "GFR<90", Lab_GFR_from_90 = 1, T_score_HIP = -1.5)
  data <- eventReactive(input$go, {
    newdata$Gender <- factor(newdata$Gender,
                             levels = c(0,1),
                             labels = c("male","female"))
    newdata$prior_falls <- factor(newdata$prior_falls,
                                  levels = c(0,1),
                                  labels = c("no","yes"))
    newdata$Fat_fish_diet <- factor(newdata$Fat_fish_diet,
                                    levels = c(0,1),
                                    labels = c("no","yes"))
    newdata$CVD <- factor(newdata$CVD,
                          levels = c(0,1),
                          labels = c("no","yes"))
    newdata$Delier_dementie_DBC <- factor(newdata$Delier_dementie_DBC,
                                          levels = c(0,1),
                                          labels = c("no","yes"))
    newdata$Change_length <- factor(newdata$Change_length,
                                    levels = c(0,1),
                                    labels = c("no","yes"))
    newdata$Lab_GFR_decreased <- factor(newdata$Lab_GFR_decreased,
                                        levels = c(0,1),
                                        labels = c("GFR>=90", "GFR<90"))
    converted_eGFR = 90-input$eGFR
    # eGFR
    if (converted_eGFR > 0 ){
      nierfunctie = "GFR<90"
    } else {
      nierfunctie = "GFR>=90"
    }
    
    
    # Prediction on new data
    newdata <- data.frame(Age=input$age, Gender = input$gender, prior_falls = input$prior_falls, CVD = input$cvd, Delier_dementie_DBC = input$delier_dementie, Change_length=input$change_length, Hours_active = input$hours_active, Fat_fish_diet = input$fat_fish, Lab_GFR_decreased = nierfunctie, Lab_GFR_from_90 = converted_eGFR, T_score_HIP = input$T_score_hip)
    pred_survival_matrix <- as.data.frame(matrix(c("3 year", "5 year")))
    
    for (i in 1:30){
      data10 <- as.data.frame(complete(MOF_imputed_good, i))
      dd <- datadist(data10)
      options(datadist="dd")
      cox_function_final_model <- coxph(Surv(time_w, MOF_fracture) ~ Age + Gender +  prior_falls + CVD + Delier_dementie_DBC  + Change_length + Hours_active +  Fat_fish_diet + Lab_GFR_decreased + ns(Lab_GFR_from_90, 3) + Age:T_score_HIP + T_score_HIP,
                                        data = data10,
                                        x = TRUE,
                                        y = TRUE)
      fit <- predictSurvProb(cox_function_final_model, newdata = newdata, times = c(156,260))
      fit_dataframe<- as.data.frame(matrix(fit, nrow=2, ncol=1))
      pred_survival_matrix <- cbind(pred_survival_matrix, fit_dataframe)
    }
    prediction_all <- pred_survival_matrix[2:31]
    Three_year_risk = (1-rowMeans(prediction_all[1,]))
    Five_year_risk =  (1-rowMeans(prediction_all[2,]))
    Three_year_risk_text <- paste("The 3 year risk of Major Osteoporotic Fracture =", percent(Three_year_risk))
    Five_year_risk_text <- paste("The 5 year risk of Major Osteoporotic Fracture =", percent(Five_year_risk))
    paste(Three_year_risk_text, Five_year_risk_text, "" , "With Major Osteoporotic Fracture defined as a fracture from the hip, wrist, spine or humerus", sep="\n")
  })
  output$text <- renderText({
    data()
  })
  output$text_age <- renderText({paste("Disclaimer: The Major Osteoporotic Fracture Risk Calculator still needs further (external) validation, and clinical evaluation before it can support decision making.",
                                       "",
                                       "Risk factors should be used as follows:",
                                       "",
                                       "Age: Current age in years",
                                       "Gender: Male or female",
                                       "Prior falls: if prior falling occured, tick yes",
                                       "Cardiovascular disease: a history or a current comorbidity related to the cardiovascular system. Examples are arrythmia, coronary heart disease and peripherial arterial disease amongst others.",
                                       "Delirium or dementia: a history of delirium or presence of dementia.",
                                       "Change in length in recent years: fill in yes when the patients noticed diminishing length in the last few years.",
                                       "Hours of moderate activity per week: moderate activity as defined by the WHO. Examples include gardening, housework and walking domestic animals.",
                                       "Dietary use of fat fish: tick yes if a patient reports including fat fish at least twice a week in his/her diet.",
                                       "Last known eGFR: Enter the last known eGFR with no decimal. If the eGFR exceeds 90, fill in 90.",
                                       "T-score of the hip: Enter the last known T-score of the hip using a comma.",
                                       sep="\n")})
}
shinyApp(ui = ui, server = server)
