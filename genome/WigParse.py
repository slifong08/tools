def getBigWigVal(test_bed, bw, outpath, src_bwavg):
    
    """
    write file.bed w/ wiggle track value average for coordinates in bed file
    
    input
        test_bed (str) - full path to bed file
        bw (str) - full path to bigwig file to parse
        outpath (str) - path to write outfile
        src_bwave (str) - path to bigWigAverageOverBed executable
    
    method
        1. string split the bed file id
        2. create out file
        3. check if path to bigWigAverageOverBed executable is provided. Else, provide
        4. build bigWigAverageOverBed command
        5. check that the output file does not already exist. 
            5a. if not, run command in commandline
        6. return the new outfile (str) to results
    
    return
        out (str) - full path to .bed file w/ wiggle value mean as 5th column.
        
    Note about bigWigAverageOverBed from kentutils:
    
        bigWigAverageOverBed v2 - Compute average score of big wig over each bed, which may have introns.
            
            usage:
               bigWigAverageOverBed in.bw in.bed out.tab -bedOut= <outfile>
    
    """
    
    #1
    id_name = (test_bed.split("/")[-1]).strip(".bed") # get the file id
    
    #2
    out = os.path.join(outpath, f"{id_name}-phylop100way.bed")  # out file str

    #3
    if src_bwavg is None:
        src_bwavg = "/wynton/home/ahituv/fongsl/nullomers/src/bigWigAverageOverBed"
    #4    
    cmd = f"{src_bwavg} {bw} {test_bed} {out}.tab -bedOut={out}"
    
    #5
    if os.path.exists(out) is False:
        os.system(cmd)
        print(cmd)
    #6    
    return out


def parallelBigWigVal(file_list, bw, outpath, src_bwavg):
    
    """
    run parallel jobs to parse bigwigs for mean values


    input
        file_list (list) - list of strs w/ full paths to bed files to run in parallel
        bw (str) - full path to bigwig file to parse
        outpath (str) - path to write outfile
        src_bwave (str) - path to bigWigAverageOverBed executable
        
    method
        1. run parallel jobs
        2. return list of outfiles
    
    return 
        list of .bed files intersected w/ bigwig file, average values for each coordinate set. 
        
    notes
        ncores is fixed at 16
        
    """
    ncores=16
    
    exps_phylop = Parallel(
                        n_jobs=ncores, verbose=100, prefer="threads")\
                        (delayed(getBigWigVal)\
                        (e, bw,outpath, src_bwavg) for e in file_list)
    return exps_phylop
    