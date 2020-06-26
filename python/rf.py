import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.model_selection import cross_val_predict
from multiprocessing import Pool

import tools

'''
Organizing the data
'''
# How to generate the RF's predictions for the whole dataset
KF = False
OOB = True

# Whether this is running on Windows
WINDOWS = True

# Whether to remove contacts who are PCR- but seropositive
OMIT_DISC = True

# Reading in the data
if WINDOWS:
    file_dir = WINDOWS_FILE_DIR
else:
    file_dir = UNIX_FILE_DIR

records = pd.read_csv(file_dir + 'records.csv')

# List of symptom names and case definitions
symptom_list = ['wheeze', 'throat', 'sob', 'nausea', 'myalgia', 'headache',
                'fatigue', 'discomf', 'diarrhea', 'cough', 'chestpain', 
                'abdpain', 'fever_chills', 'nasal_combo', 'tastesmell_combo']

# Pinning down the inputs and targets
X_bin = np.array(records[symptom_list], dtype=np.uint8)

# Casting the categorical features to one-hot arrays
cat_cols = ['ethnicity', 'sex', 'race']
X_cat = pd.concat([tools.sparsify(records[cat]) for cat in cat_cols],
                  axis=1)

# Putting all the features together
age = records.age_adult.values.reshape(-1, 1)
#cond = records.any_cond2.values.reshape(-1, 1)
X = np.concatenate([X_bin, X_cat, age], axis=1)

# Making a list of the feature names
cat_names = [name for name in X_cat.columns.values]
feature_names = symptom_list + cat_cols + ['age_adult'] #+ ['any_cond2']

# Pulling out the target
y = np.array(records.pcr_pos, dtype=np.uint8)
sero = np.array(records.sero_pos, dtype=np.uint8)

# Optionally removing sero-pos/PCR-neg contacts
if OMIT_DISC:
    # Dividing the records by sero/PCR status
    conc = np.where([not (y[i] == 0 and sero[i] == 1) 
                     for i in range(y.shape[0])])[0]
    disc = np.where([y[i] == 0 and sero[i] == 1 
                     for i in range(y.shape[0])])[0]
    
    # Removing conflicting contacts from the data
    X_disc = X[disc]
    y_disc = y[disc]
    X = X[conc]
    y = y[conc]
    
    # Giving conflicting contacts their own dataframe
    disc_df = records.iloc[disc, :]
    conc_df = records.iloc[conc, :]

'''
Running RFs in a loop to get average predictions for each person in the data
'''
if KF and not OMIT_DISC:
    print('Running the k-fold cross-validation loop.')
    n_splits = 5
    n_runs = 100
    scores = np.zeros(shape=(X.shape[0], n_runs))
    votes = np.zeros(shape=(X.shape[0], n_runs))
    np.random.seed(10221983)
    seeds = np.random.randint(1, 1e6, n_runs)
    
    for run in range(n_runs):
        print(run)
        kf = StratifiedKFold(n_splits=n_splits, 
                             shuffle=True, 
                             random_state=seeds[run])
        rf = RandomForestClassifier(n_estimators=500,
                                    class_weight='balanced',
                                    random_state=seeds[run]+1,
                                    n_jobs=-1)
        cv_score = cross_val_predict(rf, 
                                     X=X, 
                                     y=y, 
                                     cv=kf, 
                                     method='predict_proba')[:, 1]
        scores[:, run] = cv_score
        votes[:, run] = tools.threshold(cv_score)
    
    # Taking the mean of the scores and votes
    records['kf_score'] = np.mean(scores, axis=1)
    records['kf_vote'] = tools.threshold(records.kf_score)

if OOB:
    print('Running the single random forest for getting the OOB score.')
    # Running a single gigantic random forest to see how the scores change
    rf = RandomForestClassifier(n_estimators=100000,
                                oob_score=True,
                                class_weight='balanced',
                                random_state=10221983,
                                n_jobs=-1)
    rf.fit(X, y)
    oob_score = rf.oob_decision_function_[:, 1]
    
    if OMIT_DISC:
        # Getting oob scores for the concordant contacts
        conc_df['rf_score'] = oob_score
        conc_df['rf_vote'] = tools.threshold(oob_score)
        
        # Getting predicted scores for the discordant conctacts
        pos_probs = rf.predict_proba(X_disc)[:, 1]
        disc_df['rf_score'] = pos_probs
        disc_df['rf_vote'] = tools.threshold(pos_probs)
        
        # Writing to disk
        conc_df.to_csv(file_dir + 'conc_records.csv', index=False)
        disc_df.to_csv(file_dir + 'disc_records.csv', index=False)
        records = pd.concat([conc_df, disc_df], axis=0)
    
    else:
        records['oob_score'] = pd.Series(oob_score)
        records['oob_vote'] = tools.threshold(oob_score)

# Writing the results to disk
records.to_csv(file_dir + 'records.csv', index=False)
