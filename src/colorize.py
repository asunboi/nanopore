# importing pandas
from matplotlib import pyplot as plt 
import pandas as pd
import sys
from itertools import islice
import ast
import json

def color_substrings(array1, string):
    # ANSI escape codes for text color
    GREEN = "\033[92m"
    RED = "\033[91m"
    BLACK = "\033[30m"
    RESET = "\033[0m"
    
    # Process and append the substrings to a string
    colored_output = ""

    lbound = 0
    for start, end, _, operation in array1:
        colored_output += string[lbound:start-1]
        substring = string[start-1:end]
        if operation == "plus":
            colored_output += GREEN + substring + RESET
        elif operation == "minus":
            colored_output += RED + substring + RESET
        lbound = end
    colored_output += string[array1[-1][1]:]
    return colored_output

# Specify the JSON file path
json_file = "all.substrings.json"

# Read the dictionary from the JSON file
with open(json_file, "r") as file:
    substringsDict = json.load(file)
file.close()

reads = pd.read_csv("final_reads.tsv", sep='\t', header=None)

# display DataFrame
readsArray = reads.to_numpy()
readsDict = {}

for i in readsArray:
    if i[0] in readsDict:
        readsDict[i[0]].append(i[2])
    else:
        tempArray = [i[2]]
        readsDict[i[0]] = tempArray
        
f = open("coloredReads.txt", "w")
for key, value in readsDict.items():
    read = value[0]
    substringArray = substringsDict[key]
    colored_read = color_substrings(substringArray, read)
    f.write(">" + key + '\n')
    f.write(colored_read + '\n')
f.close()