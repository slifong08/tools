import os, sys

def make_windows(f, window_size):
    """
    return a .bed file split into n window_sizes with original .bed file information. 
    
    input - f (.bed
    """
    path, fname = os.path.split(f) # get the file path, names
    fname = fname.split(".bed")[0]
    
    outwindows = os.path.join(path, f"{fname}_windows-{window_size}.bed")  # produce a file
    
    if os.path.exists(outwindows) is False:
        # make windows
        window_cmd = f"bedtools makewindows -b {f} -n {window_size} -i winnum > {outwindows}"
        print(window_cmd)
        os.system(window_cmd)
        os.chdir(path)

        # add back original enhancer ids
        int_window_cmd = f"bedtools intersect -a {outwindows} -b {f} -wa -wb > t && mv t {outwindows}"
        print(int_window_cmd)
        os.system(int_window_cmd)
    
    return outwindows