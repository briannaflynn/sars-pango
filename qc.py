import sys
import json

mean_depth = sys.argv[1]
passing_count = sys.argv[2]
coverage = sys.argv[3]
outf = sys.argv[4]

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
with open(passing_count, 'r') as c:
    for line in c.readlines():
        line = line.strip()
        print(line)
        output["passing reads"] = line

#Mean depth from file
with open(mean_depth, 'r') as x:
    for line in x.readlines():
        line = line.strip()
        print(line)
        output["mean coverage"] = line

covdata = {}
maxnum = 0
with open(coverage, 'r') as p:
    for line in p.readlines():
        line = line.strip()
        if line.startswith("all"):
            line = line.split("\t")
            maxnum = int(line[3])
            covdata[line[1]] = {
                "count": line[2],
                "percent": line[4]
            }

one = 0.0
onesum = 0
five = 0.0
fivesum = 0
ten = 0.0
tensum = 0
hundred = 0.0
hundredsum = 0

for depth, props in covdata.items():
    if int(depth) >= 1:
        onesum += int(props["count"])
    if int(depth) >= 5:
        fivesum += int(props["count"])
    if int(depth) >= 10:
        tensum += int(props["count"])
    if int(depth) >= 100:
        hundredsum += int(props["count"])

output["1x"] = "{:.1%}".format(onesum/maxnum)
output["5x"] = "{:.1%}".format(fivesum/maxnum)
output["10x"] = "{:.1%}".format(tensum/maxnum)
output["100x"] = "{:.1%}".format(hundredsum/maxnum)

#determine passing criteria
if int(output["passing reads"]) < 10000:
    output["qc"] = "Failed"
else:
    output["qc"] = "Passed"

json.dump(output, open(outf, 'w'), indent=2)






