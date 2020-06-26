'''This script runs the primary analysis. It produces 2 main tables: 

  1. rule_cis.xlsx -- the performance metrics and accompanying CIs for each rule
  2. slim_diff_cis.xlsx -- adult-child differences in sens and spec
 
We use bias-corrected and accelerated (BCa) bootstrap CIs, so each CI requires
running a jackknife to construct. When the number of bootstrap samples (N_BOOT)
is high, it will take a while to run the script, even on machines that support
multiprocessing. 
'''

import numpy as np
import pandas as pd
import itertools
import time
import pickle

from importlib import reload
from multiprocessing import Pool
from sklearn.metrics import f1_score

import tools
import multi


# Whether to adjust confidence intervals for single defs and symptoms
BONF_SINGLE = False

# Whether to omit the 13 sero+/PCR- contacts from the primary analysis
OMIT_DISC = True

# Number of bootstrap samples to take for calculating CIs
N_BOOT = 10000

# Significance level to use for statistical testing
ALPHA = 0.05

# Number of rules to consider from each group for comparison
TOP_N = 3

# What to round floating-point metrics to
ROUND = 3

# Whether to write results to Excel workbooks
EXCEL = True

# Whether to reformat confidence intervals
REFORMAT_CIS = False

# Whether this is running on Windows
WINDOWS = True

# Reading in the data
if WINDOWS:
    file_dir = WINDOWS_FILE_DIR
else:
    file_dir = UNIX_FILE_DIR

records = pd.read_csv(file_dir + 'records.csv')

# Optionally removing sero+/PCR-
if OMIT_DISC:
    concordants = np.where([not (records.pcr_pos[i] == 0
                            and records.sero_pos[i] ==1)
                            for i in range(records.shape[0])])[0]
    records = records.iloc[concordants, :]

# List of symptom names and case definitions
symptom_list = [
    'wheeze', 'throat', 'sob', 'nausea', 'myalgia', 'headache',
    'fatigue', 'discomf', 'diarrhea', 'cough', 'chestpain', 
    'abdpain', 'fever_chills', 'nasal_combo', 'tastesmell_combo'
]
case_list = [
    'ili', 'cdc', 'ari', 'cste', 'cli', 'vaccine_a1', 'vaccine_a2', 
    'vaccine_a3', 'vaccine_A_all'
]
new_list = ['mc1', 'mc2', 'mc3', 'mc4']
rule_lists = [symptom_list, case_list, new_list]
rule_names = ['symptoms', 'existing defs', 'new defs']

# Organizing by age group
kids = np.where(records.age_adult == 0)[0]
adults = np.where(records.age_adult == 1)[0]

# Pinning down the inputs and targets
X = np.array(records[symptom_list], dtype=np.uint8)
y = np.array(records.pcr_pos, dtype=np.uint8)

# Making separate inputs and targets for kids, adults, and everyone
X_list = [X, X[adults], X[kids]]
y_list = [y, y[adults], y[kids]]

# Combining the two lists into a single list of tuples for easier looping
strata = [(X_list[i], y_list[i]) for i in range(len(X_list))]
strata_idx = [np.array(list(range(X.shape[0]))), adults, kids]
strata_names = ['all', 'adults', 'kids']

# Initializing a multiprocessing.Pool to help with the calculations
p = Pool()

# Making an empty list to hold the stratum-specific CIs
ci_dfs = []

# Calculating metrics and bootstrap CIs for all rules in all strata
for i, idx in enumerate(strata_idx):
    # Getting the stratum name
    s_name = strata_names[i]
    
    # Otionally adjusting alpha for multiple comparisons
    alpha = ALPHA
    if BONF_SINGLE:
        m = len(all_rules)
        alpha = ALPHA / m
    
    # Running the loop
    single_counts = []
    single_cis = []
    all_rules = tools.flatten(rule_lists)
    
    for rule in all_rules:
        # Getting the cell counts for the 2x2 table
        stat = tools.clf_metrics(y[idx], 
                                 records[rule][idx])
        counts = stat[['tp', 'fp', 'tn', 'fn']]
        single_counts.append(counts)
        
        # And now calculating the bootstrap CI
        ci = multi.boot_cis(targets=y[idx], 
                            guesses=records[rule][idx], 
                            sample_by=records.hh_id[idx],
                            n=N_BOOT, 
                            a=alpha)
        single_cis.append(ci)
    
    #Pickling the CIs for later
    with open(file_dir + s_name + '_cis.pkl', 'wb') as f:
        pickle.dump(single_cis, f)
        f.close()
    
    # Isolating the dx and prev metrics and reshaping the CIs
    all_cis = []
    
    for j, ci in enumerate(single_cis):
        ci = ci.cis.transpose()
        dx_prev = pd.concat([
            ci.sens, ci.spec, ci.ppv,
            ci.npv, ci.j, ci.f1, 
            ci.rel_prev_diff], 
                            axis=0).round(ROUND)
        out_df = pd.DataFrame(dx_prev).transpose()
        out_df.columns = [
            'sens', 'sens.lower', 'sens.upper',
            'spec', 'spec.lower', 'spec.upper',
            'ppv', 'ppv.lower', 'ppv.upper',
            'npv', 'npv.lower', 'npv.upper',
            'j', 'j.lower', 'j.upper',
            'f1', 'f1.lower', 'f1.upper',
            'prev_diff', 'prev_diff.lower', 'prev_diff.upper' 
        ]
        out_df['rule'] = all_rules[j]
        all_cis.append(out_df)
    
    all_cis = pd.concat(all_cis, axis=0)
    
    # Optionally merging the point estimates with the CIs
    if REFORMAT_CIS:
        stat_list = [
            'sens', 'spec', 'ppv', 
            'npv', 'j', 'f1',
            'prev_diff'
        ]
        all_cis = tools.merge_cis(all_cis, stat_list, ROUND)
    
    # Adding the 2x2 counts for each rule to the df of CIs
    all_counts = pd.concat(single_counts, axis=0)
    all_cis = pd.concat([all_counts, all_cis], axis=1)
    
    # Adding the stratum-specific CIs to the overall results
    ci_dfs.append(all_cis)

# Writing to disk
if EXCEL:
    writer = pd.ExcelWriter(file_dir + 'rule_cis.xlsx')
    for i, df in enumerate(ci_dfs):
        df.to_excel(writer, sheet_name=strata_names[i])
    writer.save()

# Loading the adult and kid CIs 
with open(file_dir + 'adults_cis.pkl', 'rb') as f:
    adult_cis = pickle.load(f)
    f.close()

with open(file_dir + 'kids_cis.pkl', 'rb') as f:
    kid_cis = pickle.load(f)
    f.close()

# Getting the diff CIs
diff_cis = [tools.diff_boot_cis(kid_cis[i], adult_cis[i]) 
            for i in range(len(adult_cis))]

# Saving the full diff CIs to disk
if EXCEL:
    writer = pd.ExcelWriter(file_dir + 'diff_cis.xlsx')
    for i, df in enumerate(diff_cis):
        df.to_excel(writer, sheet_name=all_rules[i])
    writer.save()

# Making an empty df to hold diffs for sensitivity and specificity
slim_cols = [
        'rule', 'd_sens', 'd_sens_l', 'd_sens_u',
        'd_spec', 'd_spec_l', 'd_spec_u'
]
slim_diffs = np.zeros(shape=(len(all_rules), 7))
slim_diffs = pd.DataFrame(slim_diffs, columns=slim_cols)

# Recalculating the CIs for sens and spec with adjusted alpha
slim_cis = [tools.diff_boot_cis(kid_cis[i],
                                adult_cis[i],
                                a=ALPHA/2) 
            for i in range(len(adult_cis))]

# Filling in the table
for i, df in enumerate(slim_cis):
    slim_diffs.rule[i] = all_rules[i]
    slim_diffs.iloc[i, 1:] = df.iloc[4:6, 2:5].values.flatten()

# Writing the table to disk
if EXCEL:
    writer = pd.ExcelWriter(file_dir + 'slim_diff_cis.xlsx')
    slim_diffs.to_excel(writer, sheet_name='adults vs kids')
    writer.save()
