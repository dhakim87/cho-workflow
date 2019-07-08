import csv;
from matplotlib import pyplot as plt;

import sys;

#Scatterplot number of core reactions VS secretion

if len(sys.argv) < 2:
    print("Usage: " + sys.argv[0] + " <inputFile>");
    print("Input File is a csv file with 4 columns - cell_line, num_core_reactions, yield_micrograms, yield_nanomole.  Each row corresponds to a cell line specific model.  This type of file can be generated via the printModels.m MATLAB script");
    exit();

FILENAME = sys.argv[1];

#input is a csv file with 4 columns- cell_line, num_core_reactions, yieldInMicrograms, yieldInNanomole
#Each row of the csv file is a cell line (which has been run through mCADRE and the rest of the workflow)
#This type of csv file is constructed from printModels in MATLAB.

cell_line = []
numCore = []
yieldMG = []
yieldNMol = []
with open("geneProteinReaction/X.txt") as csvFile:
    reader = csv.reader(csvFile, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL);
    for row in reader:
        cell_line.append(row[0]);
        numCore.append(float(row[1]));
        yieldMG.append(float(row[2]));
        yieldNMol.append(float(row[3]));

plt.scatter(numCore, yieldMG);
plt.show();
