# Python code for analysis
This is the Python code we used for our analysis. In general, there are two kinds of files: modules containing functions and classes used to conduct the analysis; and scripts that reproduce certain results from the paper. Please note that it is not necessary to run each script in order to reproduce our results; with the analytic dataset we provide on this repository's home page, that can be done solely by running `primary_analysis.py`. 

## Modules and scripts
  1. `tools.py`: support functions used for the analysis
  2. `multi.py`: multiprocessing-enabled versions of functions from `tools.py`
  3. `combo_search.py`: runs the combinatorial symptom search
  4. `rf.py`: trains a random forest on the data
  5. `primary_analysis.py`: script that produces the statistics and tables in the manuscript

## Software requirements
We used Python 3.6.8 with only a few extra packages, like `scikit-learn`, `pandas`, and `numpy`. For the full list of packages and their versions, see `requirements.txt`.

## Hardware suggestions
Between the combinatorial search, random forest, and BCa bootstraps, the scripts are a bit computationally intensive. To keep runtime reasonable, we made regular use of the `multiprocessing` module to parallelize operations where possible. Because of this, for those interested in running everything from beginning to end, we recommend the use of a scientific workstation or high-performance computing cluster. We used a Dell workstation with 24 logical cores and 64GB of RAM and found things to run smoothly, although certain simulations, like the exhaustive search for compound symptom combinations, took a while (~1.5 hours). 
