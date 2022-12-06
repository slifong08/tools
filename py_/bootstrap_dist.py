import os, sys
import numpy as np
import pandas as pd

def bootstrap(data_list, size):  
    
    """
    return the discrete and relative 95% confidence intervals of a data_list 
    
    input
        list of data (list, continuous values) - any list of data values, must be continuous values. 
        size (int) - size of dataset to bootstrap from list. If None, make bootstrapped distribution from entire list
        
    method 
        1. If size is None, get the length of the list 
        2. get the median of the list
        3. set bootstrap parameters
        4. per iteration, randomly choose elements from the fold changes list w replacement
        5. append the median to the list of bootstrapped_medians
        6. turn medians into a dataframe
        7. calculate the delta distances from the population median. This centers the data.
        8. sort from largest to smallest difference
        9. get discrete 0.025 adn 0.975 quantile values of the centered median distribution. 
        10. calculate relative confidence intervals and the population median - quantile values
    """
    
    #1
    if size is None:    
        size = len(data_list) # size of distribution to bootstrap

    #2
    xbar = np.median(data_list) # get median fold-change from n shuffles (if using mean for empirical obs, use mean here)
    
    #3
    nboot = 10000 # resample 10000 times
    val = 0
    bs_medians = []
    
    #4
    while val < nboot:

        bs_dist = np.random.choice(data_list, replace = True, size = size)
        
        #5
        bsmedian = np.median(bs_dist)
        bs_medians.append(bsmedian)
        val +=1
    #6
    bs = pd.DataFrame(data = bs_medians, index = np.arange(nboot), columns = ["bs_medians"]) # make dataframe of bootstraps

    #7 center the medians distribution
    bs["deltas"] = bs.bs_medians - xbar

    #8
    bs = bs.sort_values(by = "deltas", ascending= False)
    
    #9  get discrete 95th CI
    low = bs.deltas.quantile(0.025) 
    high = bs.deltas.quantile(0.975)
    ci_discrete = [high, low]

    #10  return ci relative to data_list median 
    ci_relative = xbar - [high, low]
   
    print("CI", ci_discrete, ci_relative, xbar)
    return ci_discrete, ci_relative
