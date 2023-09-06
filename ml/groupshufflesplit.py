from sklearn.model_selection import GroupShuffleSplit
from sklearn.linear_model import Lasso, Ridge

def cv(x_, y_, group_, model_name):
    
    """
    group cross validation strategy. 
        - report best score, alpha for crossvalidation on 
        - Lasso or Ridge regression using grouped training set.
        - if you want a different linear model, make sure to import it!
    
    
    input
        x_ (np.array) - array of training x features
        y_ (np.array) - array of training y values
        group_ (np.array) - array of group identities to split on
        model_name (str) - name of model to test. 
    
    method 
        1. instantiate best_score variable, last_test set 
        2. instantiate GroupShuffleSplit method
            10 splits
            90% training, 10% validation
        3. per split, separate training into training and validation sets
        4. check differences in test indexes, 
            update last_test set w/ current test indexes
        5. test a bunch of alpha hyperparameter values to define
            regularization penalties
        6. select model to make w training set
        7. fit model
        8. score model - pearson
        9. evaluate if model w/ alpha is the best score
        
    
    """
    # 1 variable to score the best score
    best_score = 0
    
    # evaluate difference between test sets
    last_test = set()
    
    # 2 shuffle split on group, 10 splits, 90% training, 10% validation
    gss = GroupShuffleSplit(n_splits=10, 
                            train_size=0.9, 
                            random_state=42
                           )
    
    # 3 per split
    for train_index, test_index in gss.split(x_, y_, group_):
        
        # separate x_, y_ into training and validation sets
        foldx_, valx_ = x_[train_index], x_[test_index]
        foldy_, valy_ = y_[train_index], y_[test_index]

        # 4 report the differene between test_set indexes
        print(len(set(test_index).difference(last_test)))

        # update last_test
        last_test = test_index

        # 5 test a bunch of alphas
        for alpha in [0,10,100, 1000]:
            
            # 6 call Ridge model, L2 penalization
            if model_name == "Ridge":
                model = Ridge(alpha=alpha)
                
            # call Lasso model, L1 penalization
            elif model_name == "Lasso":
                model = Lasso(alpha=alpha)
            
            # 7 fit model
            model.fit(foldx_, foldy_)
            
            # 8 evaluate score on validation set
            score = model.score(valx_, valy_)
            
            # 9 report only improved scores
            if score > best_score:
                best_score = score
                print(score, "model:", model, "alpha:" alpha)
                
                
def train_test_split_Group(X, y, group, train_size):
    """ group split for train, test

    inputs
        X (np.array) - np array of features
        y (np.array) - np array of predictor
        groups (np.array) - np array of group labels
        train_size (float) - 1 > decimal >0, fraction for training set
    
    method
        1. instantiate GroupShuffleSplit method
            1 split (train, test)
            train_size train set, (1-train_size) test set
        2. split data into train and test sets
        3. print information about training and test sets

    return
        trainx, testx (np.arrays) - arrays of X data, split into train and test
        trainy, testy (np.arrays) - arrays of Y data, split into train and test 
        train_index, test_index - arrays of training and test indexes
    
    """
    
    #1
    gss = GroupShuffleSplit(n_splits=1, train_size=train_size, random_state=42)
    
    #2
    for train_index, test_index in gss.split(X, y, group):
        
        trainx, testx = X[train_index], X[test_index]
        trainy, testy = y[train_index], y[test_index]
        
        print(f"Fold {i}:")
        print(f"  Train: index={train_index}, group={np.unique(groups[train_index])}")
        print(f"  Test:  index={test_index}, group={np.unique(groups[test_index])}")
        print(len(set(np.unique(groups[train_index])).intersection(set(np.unique(groups[test_index])))))

    return trainx, trainy, testx, testy, train_index, test_index