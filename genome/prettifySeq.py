def prettifySeq(original, mut):
    """
    prettify sequence. 
        All matching bases will be written as "."
        All non matching bases will be written w mutated base identity. 
    """
    prettyseq = ""
    for o, m in zip(original, mut):
        if o != m:
            prettyseq += m
        else:
            prettyseq += "."

    return prettyseq