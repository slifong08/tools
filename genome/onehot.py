from sklearn.preprocessing import OneHotEncoder

def OneHotEncoderDNA(X):
    enc = OneHotEncoder(handle_unknown='ignore')
    return enc.fit_transform(X)
