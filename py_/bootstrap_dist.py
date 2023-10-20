import os, sys
import numpy as np
import pandas as pd

def bootstrap(data_list, size, stat):  
    
    """
    return the discrete and 95% confidence intervals for a stat from a data_list 
    
    input
        list of data (list, continuous values) - any list of data values, must be continuous values. 
        size (int) - size of dataset to bootstrap from list. If None, make bootstrapped distribution from entire list
        stat (float | str) - quantile to bootstrap (float, 0-1) | mean (str)
        
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
    if type(stat) is float:
        obs_stat = np.quantile(data_list, stat) # get observed stat
        
    elif stat=="mean":
        obs_stat = np.mean(data_list) # get observed stat
    
    #3
    nboot = 10000 # resample 10000 times
    val = 0  # resample count
    bs_stats = []
    
    #4
    while val < nboot:

        bs_dist = np.random.choice(data_list, replace=True, size=size)
        
        #5
        if type(stat) is float:
            bs_stat = np.quantile(bs_dist, stat)
            
        elif stat=="mean":
            bs_stat = np.mean(bs_dist)
            
        bs_stats.append(bs_stat)  # add stat to list
        val +=1  # count resample
        
    #6
    bs = pd.DataFrame(data = bs_stats, 
                      index = np.arange(nboot), 
                      columns = ["bs_stat"]) # make dataframe of bootstraps

    #7 center the stat distribution
    bs["deltas"] = bs["bs_stat"] - obs_stat

    #8
    bs = bs.sort_values(by = "deltas", ascending= False)
    
    #9  get discrete 95th CI
    low = bs.deltas.quantile(0.025) 
    high = bs.deltas.quantile(0.975)
    ci_relative = [high, low]  # assume obs value is centered at zero

    #10  return ci relative to observed stat 
    ci_discrete = obs_stat - [high, low]  # assume obs value is center
   
    print(f"measure bootstrap CI of {stat} \
          quantile| mean estimate\n observed {stat} value:", 
          obs_stat,  
          "\ndiscrete diff from observed:", ci_discrete, 
          "\nrelative diff from observed:", ci_relative)
    
    return ci_discrete, ci_relative
