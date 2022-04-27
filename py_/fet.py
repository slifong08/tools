import numpy as np
import pandas as pd

from scipy import stats
import statsmodels
import statsmodels.api as sm

collection_dict = {}

def get_2x2(a, b, c, d, comparison):


    obs = [[a,b], [c,d]]
    print(obs)
    OR, P = stats.fisher_exact(np.array(obs))
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
   
    newdf = pd.DataFrame({"a":[a], "b":[b], "c":[c], "d":[d],
                          "OR":[OR], 
                          "P":[P], 
                          "ci_lower" :[odds_ci[0]],
                          "ci_lower_diff" :[OR - odds_ci[0]],
                          "ci_upper" :[odds_ci[1]],
                          "ci_upper_diff" :[odds_ci[1]-OR],
                          "OR_log2" :[np.log2(OR)],
                          "ci_lower_log2" :[np.log2(odds_ci[0])],
                          "ci_lower_diff" :[np.log2(OR) - np.log2(odds_ci[0])],
                          "ci_upper_log2" :[np.log2(odds_ci[1])],
                          "ci_upper_diff" :[np.log2(odds_ci[1])-np.log2(OR)],
                          "comparison":[comparison]
                         })

    print(comparison, obs, OR, P)

    return newdf

def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)

    # add an asterisks column
    df["asterisks"] = None
    df.loc[df["reject_null"]== True, "asterisks"] = "*"

    return df
