### 1 - Importing packages
import numpy as np
import pandas as pd
from pysurvival.models.semi_parametric import NonLinearCoxPHModel
from pysurvival.utils.metrics import concordance_index
from sklearn.model_selection import KFold
from sklearn.preprocessing import scale
import time



#importing 30 datasets (datasets prepared using MICE imputation in R. For details: Cox-regresion complete dataset)
dataframe_result=pd.DataFrame({"iteration": [1,2,3,4,5,6,7,8,9,10]})
for number_dataset in np.arange(1,31,1):
    dataset = "Complete_total"+str(number_dataset)+".csv"
    MOF = pd.read_csv(dataset, sep = ",")
    print(MOF.head())
    N=7578
    # Defining the features
    features = ['prior_falls', 'current_vertebral_fracture', 'pos_family_history','bedridden', 'backpain', 'children',
            'COCP', 'Breastfeeding', 'duration_menopause', 'DM', 'CVD', 'IBD', 'CVA', 'Epilepsy',
            'Systemic_auto_immune', 'rheuma_art', 'malabsorption','renal_insufficiency', 'corticosteroids',
            'Collaps_DBC', 'Delier_dementie_DBC', 'Duizeligheid_DBC', 'Gender', 'Age',
            'weigth_under_67', 'weight_under_60', 'Change_length', 'Hours_active', 'Total_daily_calcium', 'More_6_coffee',
            'Sun_exposure', 'Fat_fish_diet', 'veg_diet', 'Vitam_pills', 'Daily_margarin', 'Lab_bse', 'Lab_calcium',
            'Lab_albumin', 'Lab_TSH', 'Lab_vit_d', 'Lab_GFR_decreased', 'Lab_GFR_from_90','T_score_LWK', 'T_score_HIP',
            'Weight', 'Length']
    print(features)
    #Standardize continuous variables
    MOF[['children','duration_menopause', 'Age', 'Hours_active',
     'Total_daily_calcium','Lab_bse', 'Lab_calcium','Lab_albumin', 'Lab_TSH',
     'Lab_vit_d', 'Lab_GFR_from_90','T_score_LWK','T_score_HIP', 'Weight','Length']] = scale(MOF[['children','duration_menopause', 'Age',
                                                                                                  'Hours_active','Total_daily_calcium','Lab_bse',
                                                                                                  'Lab_calcium','Lab_albumin', 'Lab_TSH',
                                                                                                  'Lab_vit_d', 'Lab_GFR_from_90','T_score_LWK',
                                                                                                  'T_score_HIP', 'Weight','Length']])
    #Drop first column
    MOF = MOF.drop(MOF.columns[0], axis=1)

# 10 finals models were constructed and averaged. 10-fold cross validation was used.
    kf = KFold(n_splits=10)
    kf.get_n_splits(MOF)
    c_index_df = pd.DataFrame({"c_index_test":[]})
    c_index_dataframe = "c-index"+str(number_dataset)
    for train_index, test_index in kf.split(MOF):
        c_index_multiple_models = pd.DataFrame({"c_index_test": [],"c_index_train": [],"duration": []})
        data_train = MOF.loc[train_index].reset_index(drop=True)
        data_test = MOF.loc[test_index].reset_index(drop=True)
        for x in np.arange(1,11,1):
            X_train, X_test = data_train[features], data_test[features]
            T_train, T_test = data_train['time'].values, data_test['time'].values
            E_train, E_test = data_train['event'].values, data_test['event'].values
            structure = [{'activation': 'RELU', 'num_units': 20},{'activation': 'SELU', 'num_units': 25}]
            nonlinear_coxph = NonLinearCoxPHModel(structure=structure)
            start = time.time()
            nonlinear_coxph.fit(X_train, T_train, E_train, lr=0.0001, init_method='xav_uniform', dropout=0.2)
            stop = time.time()
            duration = stop - start
            c_index = concordance_index(nonlinear_coxph, X_test, T_test, E_test)
            c_index_train = concordance_index(nonlinear_coxph, X_train, T_train, E_train)
            c_index_multiple_models = c_index_multiple_models.append({"c_index_test": c_index,"c_index_train": c_index_train, "duration": duration}, ignore_index=True)
            print(c_index_multiple_models)
            if len(c_index_multiple_models) == 10:
                mean_test = c_index_multiple_models["c_index_test"].mean()
                mean_train = c_index_multiple_models["c_index_train"].mean()
                mean_duration = c_index_multiple_models["duration"].mean()
                c_index_df = c_index_df.append({"c_index_test": mean_test}, ignore_index=True)
                print(c_index_df)
                c_index_multiple_models_test = pd.DataFrame({"c_index_test": []})
                c_index_multiple_models_train = pd.DataFrame({"c_index_train": []})
                duration_multiple_models = pd.DataFrame({"duration": []})
            else:
                print("")
        if len(c_index_df) == 10:
            dataframe_result[c_index_dataframe] = c_index_df
            print(dataframe_result)
        else:
            print("")

dataframe_result.to_csv('Final_result_averaged_models.csv')