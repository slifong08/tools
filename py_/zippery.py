import os

def rezip_file(unzipped_file_name):

    # if file is unzipped, rezip it
    
    cmd = f"gzip {unzipped_file_name}"    
    
    if os.path.exists(unzipped_file_name) is True:  # only rezip if unzipped
        os.system(cmd)

    else:
        print("zipped already. Remember to unzip next time")

    zipped_file = unzipped_file_name + ".gz"
    
    return zipped_file

def unzip_file(zipped_file_name):
    
    # if file is zipped, unzip it
    
    if os.path.exists(zipped_file_name) is True:
        cmd = f"gunzip {zipped_file_name}"
        os.system(cmd)

    else:
        print("unzipped already. Remember to rezip")

    unzipped_file = zipped_file_name.strip(".gz")
    
    return unzipped_file
    