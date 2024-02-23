import os, sys

BEDFILE = sys.argv[1]

CLEAN_BED = ".".join(BEDFILE.split(".")[:-1]) + ".cleaned.bed"
print("write:", CLEAN_BED)
with open(CLEAN_BED, "w") as writer:
    with open(BEDFILE, "r") as reader:

        # per line
        for i, line in enumerate(reader):
            start, end =  int(line.split("\t")[1]), int(line.split("\t")[2])
            #print(start, end)
            # if header, write the header line
            if "start" in line.split("\t")[1]:
                writer.write(line)
                
            # test if start coordinate is larger than end coordinate
            elif end > start and start>0:
                writer.write(line)
            else:
                print("problem line", line)
writer.close(), reader.close()