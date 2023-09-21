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
        2. get the observed stat of the list (mean, median, quantile)
        3. set bootstrap parameters
        4. per iteration, randomly choose elements from the fold changes list w replacement
        5. append the stat to the list of bootstrapped_stat
        6. turn stats into a dataframe
        7. calculate the delta distances from the population stat. This centers the data.
        8. sort from largest to smallest difference
        9. get discrete 0.025 adn 0.975 quantile values of the centered stat distribution. 
        10. calculate relative confidence intervals and actual confidence interval values (population stat - quantile values)
    """
    
    #1
    if size is None:    
        size = len(data_list) # size of distribution to bootstrap

    #2
    obs_stat = np.median(data_list) # get observed stat
    
    #3
    nboot = 10000 # resample 10000 times
    val = 0
    bs_stats = []
    
    #4
    while val < nboot:

        bs_dist = np.random.choice(data_list, replace = True, size = size)
        
        #5
        bs_stat = np.median(bs_dist)
        bs_stats.append(bs_stat)
        val +=1
    #6
    bs = pd.DataFrame(data = bs_stats, index = np.arange(nboot), columns = ["bs_stat"]) # make dataframe of bootstraps

    #7 center the stat distribution
    bs["deltas"] = bs["bs_stat"] - obs_stat

    #8
    bs = bs.sort_values(by = "deltas", ascending= False)
    
    #9  get discrete 95th CI
    low = bs.deltas.quantile(0.025) 
    high = bs.deltas.quantile(0.975)
    ci_discrete = [high, low]

    #10  return ci relative to observed stat 
    ci_relative = obs_stat - [high, low]
   
    print("CI", ci_discrete, ci_relative, obs_stat)
    return ci_discrete, ci_relative
