import sys


mean_depth = sys.argv[1]
passing_count = sys.argv[2]
coverage = sys.argv[3]

output = {
    "passing reads": "",
    "mean coverage": "",
    "1x": "",
    "5x": "",
    "10x": "",
    "100x": "",
    "qc": ""
}

#passing read count first, see shell script for more info
with open(passing_count, 'r') as p:
    for line in p.readline():
        line = line.strip()
        output["passing reads"] = line

#Mean depth from file
with open(mean_depth, 'r') as p:
    for line in p.readline():
        line = line.strip()
        output["mean coverage"] = line


with open(coverage, 'r') as p:
    for line in p.readlines():
        line = line.strip()



