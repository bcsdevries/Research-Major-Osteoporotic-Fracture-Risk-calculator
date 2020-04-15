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

# CSV 
MOF <- import("C:/Users/BCS de Vries/Documents/Studie/Afstuderen/Data/Cox_regression_MOF_2020_final_11_02.csv")
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

#add Hazard function
hazard <- nelsonaalen(MOF, time_w, MOF_fracture)
hazard

#Delete variable rec_id, add hazard
ncol(MOF)
MOF_impute <- MOF[2:51]
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

  # Interaction Age (68.53) with Gender, mean in (). 
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
Pred[1:58, "Gender_interaction"] <- 0

# mean number of predictors
rowSums(Pred)
select <- which(rowSums(Pred)>0)
select_no_interaction <- select[1:32]
mean(select_no_interaction)

#Impute with conditions and correct predictormatrix
MOF_imputed_good <- mice(MOF_impute, pred = Pred, meth = Meth,
                         m=30, maxit = 50, seed=2795, print = F)
MOF_imputed_good
densityplot(MOF_imputed_good)

# complete databases seperately
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
cox_function_pha <- coxph(Surv(time_w,MOF_fracture) ~ Age + Gender + prior_falls + current_vertebral_fracture + pos_family_history + bedridden + backpain + sqrt(children) + COCP + Breastfeeding + ns(duration_menopause, 3) +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption +  renal_insufficiency + corticosteroids  + Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC +  weigth_under_67 + weight_under_60 + Change_length + sqrt(Hours_active) + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + I(log(Lab_bse)^2) + ns(Lab_calcium, 4) + Lab_albumin + ns(Lab_TSH, 3) + Lab_vit_d + Lab_GFR_decreased + ns(Lab_GFR_from_90,3) + T_score_LWK + T_score_HIP + Weight + Length + Hours_active:Age + T_score_HIP:Age + T_score_LWK:Age + Gender:Age + Hours_active:Gender + T_score_HIP:Gender + T_score_LWK:Gender,
                          data = MOF_stacked,
                          weights = w,
                          x = TRUE,
                          y = TRUE)
summary(cox_function_pha)
print(cox.zph(cox_function_pha))

# Check assumption of proportional hazard for final model
cox_function_final_model <- coxph(Surv(time_w, MOF_fracture)~ Age + Gender + prior_falls + current_vertebral_fracture + ns(duration_menopause,3) + T_score_HIP + Epilepsy + T_score_HIP:Age,
                          data = MOF_stacked,
                          weights = w,
                          x = TRUE,
                          y = TRUE)

print(cox.zph(cox_function_final_model))

# Model for LASSO and fit with glmnet
cox_function_all <- coxph(Surv(time_w, MOF_fracture) ~ Age + Gender + prior_falls + current_vertebral_fracture + pos_family_history + bedridden + backpain + sqrt(children) + COCP + Breastfeeding + ns(duration_menopause, 3) +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption +  renal_insufficiency + corticosteroids + Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC +  weigth_under_67 + weight_under_60 + Change_length + sqrt(Hours_active) + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + I(log(Lab_bse)^2) + ns(Lab_calcium, 4) + Lab_albumin + ns(Lab_TSH, 3) + Lab_vit_d + Lab_GFR_decreased + ns(Lab_GFR_from_90,3) + T_score_LWK + T_score_HIP + Weight + Length  + Hours_active:Age  + T_score_HIP:Age + T_score_LWK:Age + Gender:Age + Hours_active:Gender + T_score_HIP:Gender + T_score_LWK:Gender,
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
  cox_function_LASSO <- coxph(Surv(time_w,MOF_fracture) ~ Age + Gender + prior_falls + current_vertebral_fracture + pos_family_history + bedridden + backpain + sqrt(children) + COCP + Breastfeeding + ns(duration_menopause, 3) +  DM + CVD + IBD + CVA + Epilepsy + Systemic_auto_immune + rheuma_art + malabsorption +  renal_insufficiency + corticosteroids + Collaps_DBC + Delier_dementie_DBC + Duizeligheid_DBC +  weigth_under_67 + weight_under_60 + Change_length + sqrt(Hours_active) + Total_daily_calcium + More_6_coffee + Sun_exposure + Fat_fish_diet + veg_diet + Vitam_pills + Daily_margarin + I(log(Lab_bse)^2) + ns(Lab_calcium, 4) + Lab_albumin + ns(Lab_TSH, 3) + Lab_vit_d + Lab_GFR_decreased + ns(Lab_GFR_from_90,3) + T_score_LWK + T_score_HIP + Weight + Length + Hours_active:Age  + T_score_HIP:Age + T_score_LWK:Age + Gender:Age + Hours_active:Gender + T_score_HIP:Gender + T_score_LWK:Gender,
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

coef_select_1se <- coef(fit, s=mean_LASSO_1se)
select_1se <- which(or(coef_select_1se>0,coef_select_1se<0))
coef_select_1se[select_1se, 1]

coef_min_total <- coef_min_total[2:31]
mean_LASSO_min <- mean(coef_min_total)

coef_select_min <- coef(fit, s=mean_LASSO_min)
select_min <- which(or(coef_select_min>0,coef_select_min<0))
coef_select_min[select_min, 1]

plot(cvfit1)

#Fit final Cox model on all 30 datasets and pool for lambda 1se
fit.Cox <- with(data=MOF_imputed_good, exp=coxph(Surv(time_w,MOF_fracture) ~ Age + Gender +  prior_falls + current_vertebral_fracture + Epilepsy + Sun_exposure + T_score_HIP + T_score_HIP:Age + ns(duration_menopause,3)))
Cox.pool <- summary(pool(fit.Cox))
Cox.pool

pool.HR <- exp(cbind(Cox.pool[,1], (Cox.pool[,1]-1.96*(Cox.pool[,2])), 
                     (Cox.pool[,1]+1.96*(Cox.pool[,2]))))
colnames(pool.HR) <- (c("HR", "95% LO", "95% UP"))
pool.HR

outcome_Cox <- cbind(Cox.pool, pool.HR)
outcome_Cox <- round(outcome_Cox,3)
print(outcome_Cox)

# Cross validation loop for 1se
matrix_c_index_total <- as.data.frame(matrix(c(1,2,3,4,5,6,7,8,9,10)))
c_index <- as.data.frame(matrix(c("training", "test" )))

for (i in 1:30){ 
  data2 <- as.data.frame(complete(MOF_imputed_good, i))
  dd <- datadist(data2)
  options(datadist="dd")
  cv_cox_function <- cph(Surv(time_w,MOF_fracture) ~ Age + Gender +  prior_falls + current_vertebral_fracture + Epilepsy + Sun_exposure + T_score_HIP + T_score_HIP*Age + rcs(duration_menopause,3),
                          data = data2,
                          surv = TRUE,
                          x = TRUE,
                          y = TRUE)
  cv_cox_function <-     pec::cindex(cv_cox_function,
                                      formula=Surv(time_w,MOF_fracture) ~ .,
                                      data = data2, 
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
pooled<- pool.scalar(Q=mean_good, U=variance_good, n=7578)
pooled$qbar
Confidence_interval_low = pooled$qbar-1.96*sqrt(pooled$t)
Confidence_interval_high = pooled$qbar+1.96*sqrt(pooled$t)
CI <- c(Confidence_interval_low,Confidence_interval_high)
print(CI)
