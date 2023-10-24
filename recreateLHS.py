# importing pandas
from matplotlib import pyplot as plt 
import pandas as pd
import sys
from itertools import islice
import ast

def create_tuples(array1, string1):
    result_tuples = []

    for arr in array1:
        start1, end1 = arr[0]
        start2, end2 = arr[1]
        substring = string1[start1 - 1:end1]
        result_tuples.append((start2, end2, substring))

    return result_tuples

def sort_and_insert_tuples(Array1, Array2):
    result = []
    i = 0

    for a2 in Array2:
        first_int_a2 = int(a2[0])
        while i < len(Array1) and first_int_a2 > int(Array1[i][1]):
            result.append(Array1[i])
            i += 1
        
        result.append(a2)
    
    result.extend(Array1[i:])
    
    return result

def concatenate_tuples(arr):
    concatenated_string = ''.join(tup[2] for tup in arr)
    return concatenated_string

def compare_and_combine(array1, array2):
    result = []

    for a2 in array2:
        for a1 in array1:
            if a1[0][0] <= a2[0] <= a1[0][1] and a1[0][0] <= a2[1] <= a1[0][1]:
                result.append((a2[0]+a1[2],a2[1]+a1[2], a2[2]))

    return result

file_path = "all.remaining.txt"  
line_count = 0
with open(file_path, "r") as file:
    for line in file:
        line_count += 1
file.close()

substringsDict = {}

for i in range(1,line_count+1):
    if i == 1:
        with open(f'{i}_filter.substrings.txt', 'r') as file:
            for line in file:
                try:
                    tempArray = eval(line)
                    key = tempArray[0]
                    values = tempArray[1:]

                    if key in substringsDict:
                        substringsDict[key].extend(values)
                    else:
                        substringsDict[key] = values
                except Exception as e:
                    print(f"Error evaluating line: {line}")
                    print(f"Error message: {e}")
        file.close()
    else:
        tempDict = {}
        with open(f'{i-1}_filter.remaining.txt', 'r') as file:
            for line in file:
                try:
                    tempArray = ast.literal_eval(line)
                    if isinstance(tempArray, list):
                        key = tempArray[0]
                        values = tempArray[1:]
                        if key in tempDict:
                            tempDict[key].extend(values)
                        else:
                            tempDict[key] = values
                    else:
                        print("The evaluated result is not a list.")
                except (SyntaxError, ValueError):
                    print("Invalid input string.")
        file.close()
        with open(f'{i}_filter.substrings.txt', 'r') as file:
            for line in file:
                try:
                    tempArray = eval(line)
                    key = tempArray[0]
                    values = tempArray[1:]
                    remainingArray = tempDict[key]
                    newValues = compare_and_combine(remainingArray, values)
                    currentValues = substringsDict[key]
                    sortedValues = sort_and_insert_tuples(newValues, currentValues)
                    substringsDict[key] = sortedValues
                except Exception as e:
                    print(f"Error evaluating line: {line}")
                    print(f"Error message: {e}")
        file.close()

f = open("LHS" + ".fa", "w")
for key, value in substringsDict.items():
    LHS = concatenate_tuples(value)
    f.write(">" + key + '\n')
    f.write(LHS + '\n')
f.close()