# importing pandas
from matplotlib import pyplot as plt 
import pandas as pd
import sys

print(sys.argv)

def remove_matches(read, matches):
    new_r = ''
    lbound = 0
    for l, h in matches:
        new_r += read[lbound:l]
        lbound = h
    new_r += read[matches[-1][1]:]
    return new_r

# read text file into pandas DataFrame
df = pd.read_csv(sys.argv[1], sep='\t', header=None)

# display DataFrame
paramsArray = df.to_numpy()
paramsDict = {}

for i in paramsArray:
    if i[0] in paramsDict:
        paramsDict[i[0]].append(i[2])
    else:
        tempArray = [i[2]]
        paramsDict[i[0]] = tempArray

#print(paramsDict)

blast = pd.read_csv(sys.argv[2], sep="\t", header=None)

blastArray = blast.to_numpy()
blastDict = {}

for i in blastArray:
    qstart=int(i[1])-1
    qend=int(i[2])-1
    if i[0] in blastDict:
        blastDict[i[0]].append((qstart, qend))
    else:
        tempArray = [(qstart,qend)]
        blastDict[i[0]] = tempArray
    
for i in paramsDict:
    read = paramsDict[i][0]
    newRead = remove_matches(read, blastDict[i])
    paramsDict[i][0] = newRead

f = open(sys.argv[3], "w")
for key, value in paramsDict.items():
    f.write(">" + key + '\n')
    f.write(value[0] + '\n')
f.close()
