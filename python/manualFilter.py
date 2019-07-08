import csv;
from matplotlib import pyplot as plt

from collections import defaultdict;
import numpy as np;

def readCSV(filename, parseHeader = True):
    header = None
    data = []
    with open(filename) as csvFile:
        reader = csv.reader(csvFile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for row in reader:
            if parseHeader:
                header = row
                parseHeader = False
            else:
                data.append(row)
    return (header, data)


(header, data) = readCSV("sampleTPM.csv")
startCol = 2

MIN_THRESH = 1;  #1 TPM isn't very much...
filterX = []
filterY = []
labels = []
for row in data:
    if len(row[0]) == 0:
        #These are the transfected genes, skip.
        continue;
    maxVal = 0
    numZeros = 0
    for val in row[startCol:]:
        val = float(val)
        if val <= MIN_THRESH:
            numZeros += 1;
        if val > maxVal:
            maxVal = val;

    if numZeros >= 10 and numZeros <= len(header) - startCol - 10 and maxVal >= 100:
        filterX.append(numZeros)
        filterY.append(maxVal)
        labels.append(row[0] + " / " + row[1])


fig = plt.figure(figsize=(1,1))
plt.xlabel("Num 0's")
plt.ylabel("Max Expression")

plt.scatter(filterX, filterY);
for x,y,label in zip(filterX, filterY, labels):
    plt.annotate(label, xy=(x,y), xytext=(20,0), textcoords='offset points')
plt.show();
