import numpy as np
import pandas as pd
import xlsxwriter
import itertools
import time

from multiprocessing import Pool
from sklearn.metrics import f1_score

import tools
import multi

'''
Setting some basic parameters for the experiment
'''
# Whether to calculate statistics for any-of-n combinations
RUN_SINGLE = False

# Whether to calculate statistics for combinations of combinations (metacombos)
RUN_META = True

# Whether to omit the 13 sero+/PCR- contacts from the primary analysis
OMIT_DISC = True

# Whether to write results to Excel workbooks
EXCEL = True

# Whether this is running on Windows
WINDOWS = True

'''
Importing and organizing the data
'''
# Reading in the data
if WINDOWS:
    file_dir = WINDOWS_FILE_DIR
else:
    file_dir = UNIX_FILE_DIR

records = pd.read_csv(file_dir + 'records.csv')

# Optionally removing folks who are sero+ but PCR-
if OMIT_DISC:
    no_disc = np.where([not (records.pcr_pos[i] == 0 and
                             records.sero_pos[i] == 1) 
                        for i in range(records.shape[0])])[0]
    records = records.iloc[no_disc, :]

# List of symptom names and case definitions
symptom_list = [
    'wheeze', 'throat', 'sob', 'nausea', 'myalgia', 'headache',
    'fatigue', 'discomf', 'diarrhea', 'cough', 'chestpain', 
    'abdpain', 'fever_chills', 'nasal_combo', 'tastesmell_combo'
    ]

# Pinning down the inputs and targets
X = np.array(records[symptom_list], dtype=np.uint8)
y = np.array(records.pcr_pos, dtype=np.uint8)

# Organizing by age group
kids = np.where(records.age_adult == 0)[0]
adults = np.where(records.age_adult == 1)[0]

# Making separate inputs and targets for kids, adults, and everyone
X_list = [X, X[adults], X[kids]]
y_list = [y, y[adults], y[kids]]

# Combining the two lists into a single list of tuples for easier looping
groups = [(X_list[i], y_list[i]) for i in range(len(X_list))]
group_idx = [np.array(list(range(X.shape[0]))), adults, kids]
group_names = ['all', 'adults', 'kids']

'''
Calculating performance for the any-of-n (single)
and m-of-n [and/or] m-of-n (meta) combinations
'''
# Starting the pool
p = Pool()

# Setting the maximum combination size
c_min = 1
c_max = 5
c_list = list(range(c_min, c_max+1))

# Generating the combos
n_symps = range(len(symptom_list))
combos = [[(list(group), k) 
           for group in itertools.combinations(n_symps, k)]
          for k in list(range(c_min, c_max+1))]
combos = tools.flatten(combos)
nums = [combo[1] for combo in combos]
col_combos = [combo[0] for combo in combos]

# Running the combo comparisons for the different age groups
cstat_list = []
cy_list = []
cnames_list = [[symptom_list[n] for n in combo] for combo in col_combos]

if RUN_SINGLE:
    for i, tup in enumerate(groups): 
        # Pulling out the features and targets for each stratum
        ftrs = tup[0]
        tgts = tup[1]
        
        # Getting the performance stats for each combo
        combo_input = [ftrs[:, cols] for cols in col_combos]
        combo_y = p.map(tools.rowsums, combo_input)
        combo_stats = pd.concat(p.starmap(tools.clf_metrics, 
                                        [(tgts, preds) for preds in combo_y]),
                                axis=0)
        combo_stats['rule'] = cnames_list
        combo_stats['combo_size'] = nums
        cstat_list.append(combo_stats)
        cy_list.append(combo_y)
    
    # Writing the combo stats to csv
    if EXCEL:
        writer = pd.ExcelWriter(file_dir + 'combo_stats.xlsx')
        for i, df in enumerate(cstat_list):
            df['se.sp'] = df.sens + df.spec
            df.to_excel(writer, sheet_name=group_names[i])
        writer.save()

# Calculating performance for the metacombinations
if RUN_META:
    meta_iter = [pair for pair in itertools.combinations(col_combos, 2)]
    metacombos = p.map(tools.unique_combo, meta_iter)
    metacombos = [c for c in metacombos if c is not None]
    
    # Combos of m for m-of-n
    mcs = [pair for pair in itertools.permutations(range(1, 6), 2)]
    mcs += [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    
    # Sums for at least 2 
    raw_pairsums = [[(X, c, mc) for c in metacombos] for mc in mcs]
    
    # Weeding out combinations that don't make sense
    pairsum_input = []
    for i, pairs in enumerate(raw_pairsums):
        good_pairs = [pair for pair in pairs if len(pair[1][0]) >= pair[2][0] 
                      and len(pair[1][1]) >= pair[2][1]]
        pairsum_input.append(good_pairs)
    
    # Empty list to hold the best indices from each run
    top_performers = []
    
    # Max number of combos to consider from 'and', 'or', and 'any'
    best_n = 10
    
    for run_num, input in enumerate(pairsum_input):
        print(run_num)
        # Getting the rowsums for each of the combos
        psums = p.starmap(tools.pairsum, input)
        
        # Converting the pair of rowsums to a single column based on 
        # the specified logical criterium ('any', 'or', or 'and')
        csums = p.map(tools.combo_sum, [ps for ps in psums])
        
        # Calculating f1 score for each of the combo sums
        and_f1s = p.starmap(f1_score, [(y, cs[:, 0]) for cs in csums])
        any_f1s = p.starmap(f1_score, [(y, cs[:, 1]) for cs in csums])
        
        # Pulling out the best from each one
        top_and = np.argsort(and_f1s)[::-1][0:best_n]
        top_any = np.argsort(any_f1s)[::-1][0:best_n]
        top_all = [top_and, top_any]
        modes = ['both', 'any']
        prefixes = ['both of ', 'any of ']
        
        # Running the full metrics on the best from each group
        for j, top_idx in enumerate(top_all):
            prefx = prefixes[j]
            mode = modes[j]
            group_best = []
            
            for i in top_idx:
                pair = input[i]
                pair_cols = pair[1]
                pair_m = pair[2]
                
                # Making the string specifying the condition
                colnames = []
                for cols in pair_cols:
                    colnames.append(' '.join([str(cnames_list[c][0]) 
                                                 for c in cols]))
                s1 = prefx + str(pair_m[0]) + ' from [' + str(colnames[0])
                s2 = '] and ' + str(pair_m[1]) + ' from [' + str(colnames[1])
                s = s1 + s2 + ']'
                mtx = tools.combo_metrics(X, y, pair_cols, pair_m, mode)
                mtx['cond'] = s
                group_best.append(mtx)
            
            # Adding the results to the list
            top_performers.append(pd.concat(group_best, axis=0))
        
        # Writing results for the top combos to disk
        top_df = pd.concat(top_performers, axis=0)
        top_df.to_csv(file_dir + 'metacombo_stats.csv', index=False)

'''
Reconstructing some of the top-performing metacombinations
'''
# Bringing in the metacombination results
mc_df = pd.read_csv(file_dir + 'top_meta.csv')

# Reconstructing some of the best combos
sm = records.sob + records.myalgia
dfc = records.discomf + records.fever_chills
smdfc = np.array(sm + dfc >= 3, dtype=np.uint8)
mc1 = np.array(records.tastesmell_combo + smdfc > 0,
               dtype=np.uint8)

dctm = records.tastesmell_combo + records.discomf
dctm = np.array(dctm >= 1, dtype=np.uint8)
wsobfc = records.wheeze + records.sob + records.fever_chills
wsobfc = np.array(wsobfc >= 2, dtype=np.uint8)
mc2 = np.array(dctm + wsobfc > 0, dtype=np.uint8)

wdsfc = records.wheeze + records.discomf + records.sob + records.fever_chills
wdsfc = np.array(wdsfc >= 2, dtype=np.uint8)
mc3 = np.array(records.tastesmell_combo + wdsfc > 0,
               dtype=np.uint8)

sob_fc = records.sob + records.fever_chills
sob_fc = np.array(sob_fc == 2, dtype=np.uint8)
mc4 = np.array(records.tastesmell_combo + sob_fc > 0, 
               dtype=np.uint8)

# Folding the proposed case defs back into the records
records['mc1'] = mc1
records['mc2'] = mc2
records['mc3'] = mc3
records['mc4'] = mc4
records.to_csv(file_dir + 'records.csv', index=False)
