def calculateEmpiricalP(obs, exp_sum_list):
    mu = np.mean(exp_sum_list)
    sigma = np.std(exp_sum_list)
    dist_from_mu = [exp - mu for exp in exp_sum_list]
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu))

    fold_change = (obs + 1.0) / (mu + 1.0)
    p_val = (p_sum + 1.0) / (len(exp_sum_list) + 1.0)

    return "%d\t%.3f\t%.3f\t%.3f\t%.3f" % (obs, mu, sigma, fold_change, p_val)
