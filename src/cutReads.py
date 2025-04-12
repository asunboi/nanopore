# importing pandas
from matplotlib import pyplot as plt 
import pandas as pd
import sys

print(sys.argv)

def remove_matches(read, matches):
    new_r = ''
    lbound = 0
    for l, h in matches:
        new_r += read[lbound:l-1]
        lbound = h
    new_r += read[matches[-1][1]:]
    return new_r

def remaining_sections(readName, read, matches):
    readLength = len(read)
    remaining_matches = [readName]
    current_start = 1
    total_removed = 0
    
    for start, end in matches:
        if start > current_start:
            adjusted_start = current_start - total_removed
            adjusted_end = start - total_removed - 1
            adjusted_range = (adjusted_start, adjusted_end)
            remaining_matches.append([(adjusted_start, adjusted_end), (current_start, start - 1), total_removed])
        current_start = end + 1
        total_removed += end - start + 1

    if current_start <= readLength:
        adjusted_start = current_start - total_removed
        adjusted_end = readLength - total_removed
        adjusted_range = (adjusted_start, adjusted_end)
        remaining_matches.append([(adjusted_start, adjusted_end), (current_start, readLength), total_removed])

    return remaining_matches

def extract_substrings(readName, input_str, substring_ranges, strand):
    substrings = [readName]
    for start, end in substring_ranges:
        sub_str = input_str[start-1:end]
        substrings.append((start, end, sub_str, strand))
    return substrings


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
    qstart=int(i[1])
    qend=int(i[2])
    if i[0] in blastDict:
        blastDict[i[0]].append((qstart, qend))
    else:
        tempArray = [(qstart,qend)]
        blastDict[i[0]] = tempArray
    
remainingArray = []
substringArray = []

for i in paramsDict:
    read = paramsDict[i][0]
    remainingSections = remaining_sections(i, read, blastDict[i])
    substring = extract_substrings(i, read, blastDict[i], sys.argv[4])
    newRead = remove_matches(read, blastDict[i])
    remainingArray.append(remainingSections)
    substringArray.append(substring)
    paramsDict[i][0] = newRead

f = open(str(sys.argv[3]) + ".remaining.txt", "w")
for value in remainingArray:
    f.write(str(value) + '\n')
f.close()

f = open(str(sys.argv[3]) + ".substrings.txt", "w")
for value in substringArray:
    f.write(str(value) + '\n')
f.close()

f = open(str(sys.argv[3]) + ".fa", "w")
for key, value in paramsDict.items():
    if (len(value[0]) > 0):
        f.write(">" + key + '\n')
        f.write(value[0] + '\n')
f.close()
