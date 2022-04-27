import os

def split_filename(file):
    path, filename= os.path.split(file)

    sample_id = os.path.splitext(filename)[0]
    
    return path, filename, sample_id