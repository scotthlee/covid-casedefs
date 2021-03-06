'''Functions used for the primary analysis'''

import pandas as pd
import numpy as np

from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from scipy.stats import binom, chi2, norm
from copy import deepcopy
from multiprocessing import Pool


def threshold(probs, cutoff=.5):
    '''Converts probabilities to class guesses.
    
    Parameters
      probs: the probabilities to be cut (float in [0, 1])
      cutoff: the probability cut point (float in [0, 1])
    
    Returns
      class guesses as ints from {0, 1}
    '''
    return np.array(probs >= cutoff).astype(np.uint8)


def mcnemar_test(targets, guesses, cc=True):
    '''Runs McNemar's test for the difference in paired proportions.
    
    Parameters
      targets: the true labels (arr of {0, 1})
      guesses: the predicted labels (arr of {0, 1})
      cc: whether to perform a continuity correction (bool)
    
    Returns
      'b': number of false negatives
      'c': number of false positivies
      'stat': chi-squared statistic
      'pval': p-value from the test
    '''
    cm = confusion_matrix(true, pred)
    b = int(cm[0, 1])
    c = int(cm[1, 0])
    if cc:
        stat = (abs(b - c) - 1)**2 / (b + c)
    else:
        stat = (b - c)**2 / (b + c)
    p = 1 - chi2(df=1).cdf(stat)
    outmat = np.array([b, c, stat, p]).reshape(-1, 1)
    out = pd.DataFrame(outmat.transpose(),
                       columns=['b', 'c', 'stat', 'pval'])
    return out


def brier_score(targets, guesses):
    '''Calculates Brier score, or mean squared error.
    
    Parameters
      targets: the true labels (arr of {0, 1})
      guesses: the predicted scores (float in (0, 1) or int from {0, 1})
    
    Returns
      Brier score (float in (0, 1))
    '''
    return np.sum((guesses - targets)**2) / targets.shape[0]


def slim_metrics(df, rules, by=None):
    '''Returns number and percent positive for a set of predicted labels.
    
    Parameters
      df: a data frame holding the columns of predicted labels
      rules: column names for the predicted labels
      by: criteria to use for counting, e.g., for calculating sensitivity
    
    Returns
      a df with the rule, n positive, and percent positive
    '''
    if by is not None:
        good_idx = np.where(by == 1)[0]
        df = df.iloc[good_idx]
    N = df.shape[0]
    out = np.zeros(shape=(len(rules), 2))
    for i, rule in enumerate(rules):
        out[i, 0] = np.sum(df[rule])
        out[i, 1] = out[i, 0] / N
    out_df = pd.DataFrame(out, columns=['n', 'pct'])
    out_df['rule'] = rules
    out_df = out_df[['rule', 'n', 'pct']]
    return out_df


def clf_metrics(targets, 
                guesses,
                average_by=None,
                weighted=True,
                round=4,
                round_pval=False,
                mcnemar=False):
    '''Calculates a range of binary classification metrics for a set of class
    predictions relative to a reference standard.
    
    Keyword arugments:
      targets: the true labels (arr of {0, 1})
      guesses: the predicted labels (arr of {0, 1})
      average_by: the variable to use for macro averaging (1-d array)
      weighted: whether to weight macro averaging (bool)
      round: number of significant digits to report
      round_pval: whether to round p-values from McNemar's test (bool)
      mcnemar: whether to run McNemar's test
    
    Returns
      a one-row data frame with the following columns:
        tp: true positive count
        fp: false positive count
        tn: true negative count
        fn: false negative count
        sens: sensitivity
        spec: specificity
        ppv: positive predictive value
        npv: negative predictive value
        j: Youden's j index
        mcc: Matthews correlation coefficient
        brier: Brier score (or 1 - acc)
        f1: F1 score
        
        true_prev: true prevalence
        pred_prev: predicted prevalence
        abs_diff: absolute difference in prevalence
        rel_prev_diff: percent difference in prevalence
        mcnemar: p-value from McNemar's test (optional)   
    '''
    
    # Converting pd.Series to np.array
    stype = type(pd.Series())
    if type(guesses) == stype:
        guesses = guesses.values
    if type(targets) == stype:
        targets = targets.values
    if type(average_by) == stype:
        average_by == average_by.values
    
    # Optionally returning macro-average results
    if average_by is not None:
        return macro_clf_metrics(targets=targets,
                                 guesses=guesses,
                                 by=average_by,
                                 weighted=weighted,
                                 round=round)
    # Constructing the 2x2 table
    confmat = confusion_matrix(targets, guesses)
    tp = confmat[1, 1]
    fp = confmat[0, 1]
    tn = confmat[0, 0]
    fn = confmat[1, 0]
    
    # Calculating basic measures of diagnostic accuracy
    sens = np.round(tp / (tp + fn), round)
    spec = np.round(tn / (tn + fp), round)
    ppv = np.round(tp / (tp + fp), round)
    npv = np.round(tn / (tn + fn), round)
    f1 = np.round(2 * (sens * ppv) / (sens + ppv), round)
    j = sens + spec - 1
    mcc_num = ((tp * tn) - (fp * fn))
    mcc_denom = np.sqrt(((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    mcc = mcc_num / mcc_denom
    brier = np.round(brier_score(targets, guesses), round)
    outmat = np.array([tp, fp, tn, fn,
                       sens, spec, ppv,
                       npv, j, f1, mcc, brier]).reshape(-1, 1)
    out = pd.DataFrame(outmat.transpose(),
                       columns=['tp', 'fp', 'tn', 
                                'fn', 'sens', 'spec', 
                                'ppv', 'npv', 'j', 
                                'f1', 'mcc', 'brier'])
    
    # Calculating some additional measures based on positive calls
    true_prev = int(np.sum(targets == 1))
    pred_prev = int(np.sum(guesses == 1))
    abs_diff = (true_prev - pred_prev) * -1
    rel_diff = np.round(abs_diff / true_prev, round)
    if mcnemar:
        pval = mcnemar_test(targets, guesses).pval[0]
        if round_pval:
            pval = np.round(pval, round)
    count_outmat = np.array([true_prev, pred_prev, abs_diff, 
                             rel_diff]).reshape(-1, 1)
    count_out = pd.DataFrame(count_outmat.transpose(),
                             columns=['true_prev', 'pred_prev', 
                                      'prev_diff', 'rel_prev_diff'])
    out = pd.concat([out, count_out], axis=1)
    
    # Optionally dropping the mcnemar p-val
    if mcnemar:
        out['mcnemar'] = pval
    
    return out


def macro_clf_metrics(targets,
                      guesses,
                      by,
                      weighted=True,
                      round=4,
                      p_method='harmonic',
                      mcnemar=True):
    '''Performs weighted or unweighted macro-averaging of clf_metrics()
    by a group variable.
    
    Parameters
      targets: the true labels(arr of {0 , 1})
      guesses: the predict labels (arr of {0, 1})
      by: an array of group IDs to use for averaging (1-d array)
      weighted: whether to return a weighted average
      round: number of significant digits to return
      p_method: how to average p-values; may be 'harmonic' or 'fisher'
      7. mcnemar: whether to run McNemar's test (bool)
     
    Returns
      the df from clf_metrics() where everything has been averaged
    '''
    # Column groups for rounding later
    count_cols = ['tp', 'fp', 'tn', 'fn']
    prev_cols = ['true_prev', 'pred_prev', 'prev_diff']
    
    # Getting the indices for each group
    n = len(targets)
    group_names = np.unique(by)
    n_groups = len(group_names)
    group_idx = [np.where(by == group)[0]
                 for group in group_names]
    group_counts = np.array([len(idx) for idx in group_idx])
    
    # Calculating the groupwise statistics
    group_stats = [clf_metrics(targets[idx],
                               guesses[idx],
                               mcnemar=mcnemar) 
                   for idx in group_idx]
    
    # Casting the basic counts as proportions
    for i, df in enumerate(group_stats):
        df[count_cols] /= group_counts[i]
        df[prev_cols] /= group_counts[i]
    
    group_stats = pd.concat(group_stats, axis=0)
    
    # Calculating the weights
    if weighted:
        w = np.array(group_counts / n)
    else:
        w = np.repeat(1 / n_groups, n_groups)
    
    # Calculating the mean values
    averages = np.average(group_stats, axis=0, weights=w)
    avg_stats = pd.DataFrame(averages).transpose()
    avg_stats.columns = group_stats.columns.values
    
    # Converting the count metrics back to integers
    avg_stats[count_cols] *= n
    avg_stats[count_cols] = avg_stats[count_cols].astype(int)
    avg_stats[prev_cols] *= n
    avg_stats.rel_prev_diff = avg_stats.prev_diff / avg_stats.true_prev
    
    # Rounding off the floats
    float_cols = ['sens', 'spec', 'npv', 
                  'ppv', 'j', 'f1', 'brier']
    avg_stats[float_cols] = avg_stats[float_cols].round(round)
    avg_stats.rel_prev_diff = avg_stats.rel_prev_diff.round(round)
    
    # Getting the mean of the p-values with either Fisher's method
    # or the harmonic mean method
    if mcnemar:
        avg_stats.mcnemar = average_pvals(group_stats.mcnemar,
                                          w=w,
                                          method=p_method)
    
    return avg_stats


def average_pvals(p_vals, 
                  w=None, 
                  method='harmonic',
                  smooth=True,
                  smooth_val=1e-7):
    '''Averages p-values using either the harmonic mean or Fisher's method.
    
    Parameters
      p_vals: the p-values (arr of floats in [0, 1])
      w: the weights for averaging
      method: either 'harmonic' (default) or 'fisher' (str)
      smooth: whether to fix pvals of 0.0 (bool)
      smooth_val: the amount to use for smoothing (float)
    
    Returns
      the average p-value (single float in [0, 1])
    '''
    if smooth:
        p = p_vals + smooth_val
    else:
        p = deepcopy(p_vals)
    if method == 'harmonic':
        if w is None:
            w = np.repeat(1 / len(p), len(p))
        p_avg = 1 / np.sum(w / p)
    elif method == 'fisher':
        stat = -2 * np.sum(np.log(p))
        p_avg = 1 - chi2(df=1).cdf(stat)
    return p_avg


def boot_sample(df,
                by=None,
                size=None,
                seed=None,
                return_df=False):
    '''Returns a single bootstrap sample of rows from a data frame.
    
    Parameters
      df: the data frame holding the records (2-d array or pd.DataFrame)
      by: an array of group IDs for sampling by group instead of row (arr)
      size: the size of bootstrap samples to take, if not nrow(df) (int)
      seed: seed to use for generating the random sample (int)
      return_df: whether to return row indices (False) or the df (True)
    
    Returns
      1a. An array of bootstrap-sampled row numbers, if return_df is False; OR
      1b. A boostrap sample of the original df, if return_df is True
    '''
    # Setting the random states for the samples
    if seed is None:
        seed = np.random.randint(1, 1e6, 1)[0]
    np.random.seed(seed)
    
    # Getting the sample size
    if size is None:
        size = df.shape[0]
    
    # Sampling across groups, if group is unspecified
    if by is None:
        np.random.seed(seed)
        idx = range(size)
        boot = np.random.choice(idx,
                                size=size,
                                replace=True)
    
    # Sampling by group, if group has been specified
    else:
        levels = np.unique(by)
        level_idx = [np.where(by == level)[0]
                     for level in levels]
        boot = np.random.choice(level_idx,
                                size=len(levels),
                                replace=True) 
        boot = np.concatenate(boot).ravel()
    
    if not return_df:
        return boot
    else:
        return df.iloc[boot, :]
    

def diff_boot_cis(ref, 
                  comp, 
                  a=0.05,
                  abs_diff=False, 
                  method='bca',
                  interpolation='nearest'):
    '''Calculates boostrap confidence intervals for the difference in
    performance metrics between two competing classifiers.
    
    Parameters
      ref: the refernece multi.boot_cis object
      comp: the comparison multi.boot_cis object
      a: significance level for the intervals (float in [0, 1])
      abs_diff: whether to take the absolute value of the difference (bool)
      method: interval method; options are 'diff', 'pct', and 'bca'
      interpolation: interpolation method for np.quantile
      
    Returns
      A pd.DataFrame with the following columns: 
        ref: the reference value for the metric
        comp: the comparison value for the metric
        d: the (absolute) difference between the ref and the comp values
        lower: the lower bound for the difference
        upper: the upper bound for the difference
    '''
    # Quick check for a valid estimation method
    methods = ['pct', 'diff', 'bca']
    assert method in methods, 'Method must be pct, diff, or bca.'
    
    # Pulling out the original estiamtes
    ref_stat = pd.Series(ref.cis.stat.drop('true_prev').values)
    ref_scores = ref.scores.drop('true_prev', axis=1)
    comp_stat = pd.Series(comp.cis.stat.drop('true_prev').values)
    comp_scores = comp.scores.drop('true_prev', axis=1)
    
    # Optionally Reversing the order of comparison
    diff_scores = comp_scores - ref_scores
    diff_stat = comp_stat - ref_stat
        
    # Setting the quantiles to retrieve
    lower = (a / 2) * 100
    upper = 100 - lower
    
    # Calculating the percentiles 
    if method == 'pct':
        cis = np.nanpercentile(diff_scores,
                               q=(lower, upper),
                               interpolation=interpolation,
                               axis=0)
        cis = pd.DataFrame(cis.transpose())
    
    elif method == 'diff':
        diffs = diff_stat.values.reshape(1, -1) - diff_scores
        percents = np.nanpercentile(diffs,
                                    q=(lower, upper),
                                    interpolation=interpolation,
                                    axis=0)
        lower_bound = pd.Series(diff_stat + percents[0])
        upper_bound = pd.Series(diff_stat + percents[1])
        cis = pd.concat([lower_bound, upper_bound], axis=1)
    
    elif method == 'bca':
        # Removing true prevalence from consideration to avoid NaNs
        ref_j_means = ref.jack[1].drop('true_prev')
        ref_j_scores = ref.jack[0].drop('true_prev', axis=1)
        comp_j_means = comp.jack[1].drop('true_prev')
        comp_j_scores = comp.jack[0].drop('true_prev', axis=1)
        
        # Calculating the bias-correction factor
        n = ref.scores.shape[0]
        stat_vals = diff_stat.transpose().values.ravel()
        n_less = np.sum(diff_scores < stat_vals, axis=0)
        p_less = n_less / n
        z0 = norm.ppf(p_less)
        
        # Fixing infs in z0
        z0[np.where(np.isinf(z0))[0]] = 0.0
        
        # Estiamating the acceleration factor
        j_means = comp_j_means - ref_j_means
        j_scores = comp_j_scores - ref_j_scores
        diffs = j_means - j_scores
        numer = np.sum(np.power(diffs, 3))
        denom = 6 * np.power(np.sum(np.power(diffs, 2)), 3/2)
        
        # Getting rid of 0s in the denominator
        zeros = np.where(denom == 0)[0]
        for z in zeros:
            denom[z] += 1e-6
        
        acc = numer / denom
        
        # Calculating the bounds for the confidence intervals
        zl = norm.ppf(a / 2)
        zu = norm.ppf(1 - (a/2))
        lterm = (z0 + zl) / (1 - acc*(z0 + zl))
        uterm = (z0 + zu) / (1 - acc*(z0 + zu))
        lower_q = norm.cdf(z0 + lterm) * 100
        upper_q = norm.cdf(z0 + uterm) * 100
                                
        # Returning the CIs based on the adjusted quantiles
        cis = [np.nanpercentile(diff_scores.iloc[:, i], 
                                q=(lower_q[i], upper_q[i]),
                                interpolation=interpolation,
                                axis=0) 
               for i in range(len(lower_q))]
        cis = pd.DataFrame(cis, columns=['lower', 'upper'])
                
    cis = pd.concat([ref_stat, comp_stat, diff_stat, cis], 
                    axis=1)
    cis = cis.set_index(ref_scores.columns.values)
    cis.columns = ['ref', 'comp', 'd', 
                   'lower', 'upper']
    
    return cis


def roc_cis(rocs, alpha=0.05, round=2):
    '''Calculates upper and lower bounds for a collection of ROC curves. 
    
    Parameters
      rocs: the sklearn.metrics.roc_curve curves
      alpha: significance value for constructing the intervals
      round: number of significant digits to report
      
    Returns
      A pd.DataFrame with the following columns:
        fpr: false positive rate
        lower: lower bound of the true positive rate (tpr)
        med: median for the tpr
        upper: upper bound for the tpr
    '''
    # Getting the quantiles to make CIs
    lq = (alpha / 2) * 100
    uq = (1 - (alpha / 2)) * 100
    fprs = np.concatenate([roc[0] for roc in rocs], axis=0)
    tprs = np.concatenate([roc[1] for roc in rocs], axis=0)
    roc_arr = np.concatenate([fprs.reshape(-1, 1), 
                              tprs.reshape(-1, 1)], 
                             axis=1)
    roc_df = pd.DataFrame(roc_arr, columns=['fpr', 'tpr'])
    roc_df.fpr = roc_df.fpr.round(round)
    unique_fprs = roc_df.fpr.unique()
    fpr_idx = [np.where(roc_df.fpr == fpr)[0] for fpr in unique_fprs]
    tpr_quants = [np.percentile(roc_df.tpr[idx], q=(lq, 50, uq)) 
                  for idx in fpr_idx]
    tpr_quants = np.vstack(tpr_quants)
    quant_arr = np.concatenate([unique_fprs.reshape(-1, 1),
                                tpr_quants],
                               axis=1)
    quant_df = pd.DataFrame(quant_arr, columns=['fpr', 'lower',
                                                'med', 'upper'])
    quant_df = quant_df.sort_values('fpr')
    return quant_df


def x_at_y(x, y, yval, grid):
    '''Calculates the maximum value of x for a minimum value of y (e.g., when
     finding maximum specifiicity achievable for a given sensitivity.)
     
     Parameters
       x: column name of the metric to optimize
       y: column name of the metric for the basis of optimizing x
       yval: the minimum value of y required for x
       grid: the pd.DataFrame holding the data
     
     Returns
       The max value of df[x] that still gets you a minimum value of df[y]
    '''
    y = np.array(grid[y])
    x = np.array(grid[x])
    assert np.sum(y >= yval) > 0, 'No y vals meet the minimum'
    good_y = np.where(y >= yval)[0]
    best_x = np.max(x[good_y])
    return best_x


def merge_cis(df, stats, round=4):
    '''Merges the upper and lower columns from a boot_cis object.
    
    Parameters
      df: the boot_cis.cis object
      stats: the list of metrics to include in the reformatted table
      round: number of significant digits to report
    
    Returns
      a pd.DataFrame with the merged CI columns
    '''
    df = deepcopy(df)
    for stat in stats:
        lower = stat + '.lower'
        upper = stat + '.upper'
        new = stat + '.ci'
        l = df[lower].values.round(round)
        u = df[upper].values.round(round)
        strs = [pd.Series('(' + str(l[i]) + ', ' + str(u[i]) + ')')
                for i in range(df.shape[0])]
        df[new] = pd.concat(strs, axis=0)
        df = df.drop([lower, upper], axis=1)
    return df


def sparsify(col, 
             return_df=True,
             long_names=False):
    '''Converts a 1-d categorical array to a one-hot matrix.
    
    Parameters
      col: the 1-d array holding the categorical variable
      return_df: whether to return the matrix as a pd.DataFrame
      long_names: whether col names should include the var name
    
    Returns
      either a np.array or pd.DataFrame with the sparsified variable
    '''
    levels = np.unique(col)
    out = np.array([col == level for level in levels],
                   dtype=np.uint8).transpose()
    if long_names:
        var = col.name + '.'
        levels = [var + level for level in levels]
    columns = [col.lower() for col in levels]
    if return_df:
        out = pd.DataFrame(out, columns=columns)
    return out


def rowsums(m, min=1):
    '''Determines which rows of m have at least min 1s'''
    sums = np.sum(m, axis=1)
    return np.array(sums >= min, dtype=np.uint8)


def pairsum(X, c, min=(1, 1)):
    '''Runs rowsums() on cols c of array X.
    
    Parameters
      X: a numpy Array
      c: a tuple or list of (int) lists of column numbers
      min: a tuple of minimum counts for the two calls to rowsums()
    
    Returns
      a tuple of 1-d arrays with the row sums for each set of columns
    '''
    a = rowsums(X[:, c[0]], min=min[0])
    b = rowsums(X[:, c[1]], min=min[1])
    return (a, b)


def combo_sum(ctup):
    '''Determines whether any or both sets of columns have rowsums greater 
    than a specified minimum.
    
    Parameters
      ctup: an array of shape (n, 2) and values in {0, 1}
    
    Returns
      an array of shape (n, 2) with the combo sums
    '''
    sums = np.sum(ctup, axis=0)
    both_y = np.array(sums == 2, dtype=np.uint8)
    any_y = np.array(sums >= 1, dtype=np.uint8)
    out = np.concatenate([both_y.reshape(-1, 1),
                          any_y.reshape(-1, 1)],
                         axis=1)
    return out


def combo_metrics(X, y, 
                  cols, 
                  min=(1, 1),
                  mode='both',
                  mcnemar=True):
    '''Runs clf_metrics() on the combo_sum() of two sets of columns.
    
    Parameters
      X: the base numpy Array of columns
      y: the true labels for prediction
      cols: the 2 sets of column numbers for summing
      mode: whether to report metrics for 'both' or 'any' combo_sum
      mcnemar: whether to run McNemar's test
    
    Returns
      the output of clf_metrics() for the combined sets of columns
    '''
    ps = pairsum(X, cols, min=min)
    cs = combo_sum(ps)
    if mode == 'both':
        return clf_metrics(y, cs[:, 0], mcnemar=mcnemar)
    elif mode == 'any':
        return clf_metrics(y, cs[:, 1], mcnemar=mcnemar)


def flatten(l):
    '''Flattens a list.'''
    return [item for sublist in l for item in sublist]


def unique_combo(c):
    '''Determines if two sets of column nums have any nums in common.'''
    if len(np.intersect1d(c[0], c[1])) == 0:
        return c
    else:
        return None
