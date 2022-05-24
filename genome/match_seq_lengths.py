import os, sys
import pandas as pd

def custom_round(x, base=10):
    return int(base * round(float(x)/base))

def match_len(Df1, Df2, base_len, columns):
    """
    match df2 on df1 sequence lengths, rounding to the nearest length (e.g. rounding to 10 bp)
    return list of region_ids in df2 matched to length distribution in df1. 
    
    args 
    
    df1: pandas dataframe w/ columns "region id" and "len". Match df2 to df1 length distributions
    df2: pandas dataframe w/ columns "region id" and "len". Match df2 to df1 length distributions
    base_len: int, sequence length range that should be matched on in distribution.
    columns: list, two column names for df. First will be the region_id (chr1:84895-85000) and the second will be the length of that id (e.g. 105 bp)
    
    process
    
    1. keep only region_id, len in df1, df2
    2. round region lengths to nearest base_len
    3. find overlapping length range for df1, df2
    4. per length, 
    5. get a list of N sample ids in df1, df2 that match length. 
    6. determine min number of region ids w/ match length to select from df1, df2. 
    7. If there are more than N sample ids w matched length in df1 than df2 (or vice versa), 
        randomly select w/ replacement region ids from that dataframe.  
    8. append df1, df2 matching length ids.
    9 return lists of id for df1, df2
    
    notes
    
    - base_len gives better power to match sequences becuase it doesn't require that two sequence match on exact bp lengths, 
        but within a base_len range. E.g. if base_len is 10, start and stop coordinates will be rounded to nearest base10 value.
    

    """
    #columns = ["region_id", "len"]
    columns_names = ["matching_ids", "matching_len"]
    columns = list(columns)
    print(columns)

    #1 df1
    df1 = Df1[columns].drop_duplicates()  # reduce Df1 to only ids and lengths
    df1.columns = columns_names  # rename columns
    
    # same for df2
    df2 = Df2[columns].drop_duplicates()
    df2.columns = columns_names
    
    #2 round df1, df2 region lengths to the nearest base_len bps
    df1["matching_len"] = df1["matching_len"].astype(float).apply(lambda x: custom_round(x, base=base_len)) 
    df2["matching_len"] = df2["matching_len"].astype(float).apply(lambda x: custom_round(x, base=base_len))
    
    #3 find intersecting lengths
    lens = set(list(set(df1["matching_len"])) + list(set(df2["matching_len"])))

    df1_match_list, df2_match_list = [], []
    
    #4 for each length, 
    for i, length in enumerate(lens):
        
        #5 count how many regions have length of one size (e.g. 100bp) in df1, df2
        Ndf1 = df1.loc[df1["matching_len"] == length].size
        Ndf2 = df2.loc[df2["matching_len"] == length].size

        #6 get the minimum number of regions to sample from df1, df2. 
        # naturally, n regions of a length for one df will be chosen, while n regions for the other df will be randomly sampled
        # Ns = number to sample
        Ns = min(Ndf1, Ndf2)

        #7 get region ids w/ matching lengths. 
        # Sample randomly from df1, df2 
        # with replacement
        
        if length > 0 and Ns > 0:
            
            df1_ids = list(df1.loc[df1["matching_len"] == length, "matching_ids"].sample(n = Ns, replace = True)) # sample w/ replacement
            df2_ids = list(df2.loc[df2["matching_len"] == length, "matching_ids"].sample(n = Ns, replace = True)) # sample w/ replacement
            
            #8 append list of ids
            df1_match_list.extend(df1_ids), df2_match_list.extend(df2_ids)

    #9 return matched ids
    return df1_match_list, df2_match_list