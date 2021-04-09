import os

file_name = 'myData.txt'

sugar = ""

with open(file_name) as base_file:
    for lineNo, line in enumerate(base_file):
        print(line)