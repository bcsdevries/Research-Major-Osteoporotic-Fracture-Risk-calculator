### 1 - Importing packages
import numpy as np
import pandas as pd
from pysurvival.models.semi_parametric import NonLinearCoxPHModel
from pysurvival.utils.metrics import concordance_index
from sklearn.model_selection import KFold
from sklearn.preprocessing import scale

#importing single dataset
MOF = pd.read_csv("Complete_subset1.csv", sep = ",")
print(MOF.head())
N=1770

# Defining the features
features = ['prior_falls', 'pos_family_history','bedridden', 'backpain', 'children',
            'COCP', 'Breastfeeding', 'duration_menopause', 'DM', 'CVD', 'IBD', 'CVA', 'Epilepsy',
            'Systemic_auto_immune', 'rheuma_art', 'malabsorption','renal_insufficiency', 'corticosteroids',
            'Collaps_DBC', 'Delier_dementie_DBC', 'Duizeligheid_DBC', 'Gender', 'Age',
            'weigth_under_67', 'weight_under_60', 'Change_length', 'Hours_active', 'Total_daily_calcium', 'More_6_coffee',
            'Sun_exposure', 'Fat_fish_diet', 'veg_diet', 'Vitam_pills', 'Daily_margarin', 'Lab_bse', 'Lab_calcium',
            'Lab_albumin', 'Lab_TSH', 'Lab_vit_d', 'Lab_GFR_decreased', 'Lab_GFR_from_90','T_score_LWK', 'T_score_HIP',
            'Weight', 'Length']

#Standardize continuous variables
MOF[['children','duration_menopause', 'Age', 'Hours_active',
     'Total_daily_calcium','Lab_bse', 'Lab_calcium','Lab_albumin', 'Lab_TSH',
     'Lab_vit_d', 'Lab_GFR_from_90','T_score_LWK','T_score_HIP', 'Weight','Length']] = scale(MOF[['children','duration_menopause', 'Age',
                                                                                                  'Hours_active','Total_daily_calcium','Lab_bse',
                                                                                                  'Lab_calcium','Lab_albumin', 'Lab_TSH',
                                                                                                  'Lab_vit_d', 'Lab_GFR_from_90','T_score_LWK',
                                                                                                  'T_score_HIP', 'Weight','Length']])

# Drop first column
MOF = MOF.drop(MOF.columns[0], axis=1)

# Grid search using 10-fold cross-validation (2 layers). Can be changed to 1 layer by removing num_units2 in for loop and structure.
kf = KFold(n_splits=10)
kf.get_n_splits(MOF)
dataframe_hp = pd.DataFrame(
    {"lr": [], "dropout": [], "num_units": [], "num_units2": [], "num_layers": [], "activation": [], "mean": []})
c_index_df = pd.DataFrame({"c_index": []})
print(dataframe_hp)
for lr in np.arange(0.0001, 0.01, 0.001):
    for dropout in np.arange(0.0, 0.5, 0.1):
        for num_units in np.arange(20, 50, 5):
            for num_units2 in np.arange(20,50,5):
                for train_index, test_index in kf.split(MOF):
                    data_train = MOF.loc[train_index].reset_index(drop=True)
                    data_test = MOF.loc[test_index].reset_index(drop=True)
                    X_train, X_test = data_train[features], data_test[features]
                    T_train, T_test = data_train['time'].values, data_test['time'].values
                    E_train, E_test = data_train['event'].values, data_test['event'].values
                    structure = [{'activation': 'SELU', 'num_units': num_units},{'activation': 'SELU', 'num_units': num_units2}]
                    nonlinear_coxph = NonLinearCoxPHModel(structure=structure)
                    nonlinear_coxph.fit(X_train, T_train, E_train, lr=lr, init_method='xav_uniform', dropout=dropout)
                    c_index = concordance_index(nonlinear_coxph, X_test, T_test, E_test)
                    c_index_df = c_index_df.append({"c_index": c_index}, ignore_index=True)
                    print(c_index_df)
                    if len(c_index_df) == 10:
                        mean = c_index_df["c_index"].mean()
                        dataframe_hp = dataframe_hp.append(
                            {"lr": lr, "dropout": dropout, "num_units": num_units, "num_units2": num_units2, "num_layers": 2,
                            "activation": "SELU/SELU", "mean": mean}, ignore_index=True)
                        c_index_df = pd.DataFrame({"c_index": []})
                        print(dataframe_hp)
                        maximize = dataframe_hp["mean"].idxmax()

                        print("Hyperparameters for final model:", dataframe_hp.iloc[maximize])
                    else:
                        print("")