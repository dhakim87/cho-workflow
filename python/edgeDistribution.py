#Daniel Hakim

from collections import defaultdict;
import scipy.stats
from matplotlib import pyplot as plt
import seaborn as sns;
import pandas as pd
import numpy as np

#Open all the files in modelResults
#For each file, write down the secretion micrograms of the cell line,
#and make a map from reactions to cell lines

BASE_DIR = "modelResultsKallehauge";
MIN_CELL_LINE = 101
MAX_CELL_LINE = 196

secretionMicrograms = {}
reactionToCellLine = defaultdict(set);

for i in range(MIN_CELL_LINE, MAX_CELL_LINE):
    filename = "./" + BASE_DIR + "/P6975_" + str(i) + "_core.txt";
    cellLine = "P6975_" + str(i);
    with open(filename) as f:
        data = f.read().splitlines();
        micrograms = float(data[0].split(", ")[0])
        secretionMicrograms[cellLine] = micrograms;
        for val in data[1:]:
            reaction = val.strip();
            reactionToCellLine[reaction].add(cellLine);
        
        print cellLine, secretionMicrograms[cellLine], "#Reactions:", len(data)-1;

# Each reaction splits cell lines into cell lines that have and don't have that reaction.
# This creates two sets of secretion microgram numbers.  We use a permutation test to see if
# these two sets come from the same distribution.  This informs on whether or not each reaction
# is relevant to secretion.

permutationTestResults = {}

for reaction in reactionToCellLine:
    reactionSet = reactionToCellLine[reaction];
    on = []
    off = []
    for cellLine in secretionMicrograms:
        if cellLine in reactionSet:
            on.append(secretionMicrograms[cellLine])
        else:
            off.append(secretionMicrograms[cellLine])
    on.sort();
    off.sort();

    if len(on) > 0 and len(off) > 0:
        mannWhitneyU = scipy.stats.mannwhitneyu(on, off)
        
        if (mannWhitneyU.pvalue > 0 and mannWhitneyU.pvalue < .05 / len(reactionToCellLine)):
            on_off = ['ON']*len(on) + ['OFF'] * len(off)
            df = pd.DataFrame({'secretion':on + off, reaction:on_off});
            df['log_secretion'] = np.log10(df['secretion'])
            plt.subplot(2,1,1);
            sns.swarmplot(x='secretion', y=reaction, data=df);
            plt.subplot(2,1,2);
            sns.swarmplot(x='log_secretion', y=reaction, data=df);
            plt.suptitle("P = " + str(mannWhitneyU.pvalue));
            plt.show();

        
        permutationTestResults[reaction] = (reaction, len(on), len(off), mannWhitneyU.statistic, mannWhitneyU.pvalue, (min(on), on[len(on)/2], max(on)), (min(off), off[len(off)/2], max(off)));
    else:
        permutationTestResults[reaction] = (reaction, len(on), len(off), 0, 0, -1, -1);

for reaction in reactionToCellLine:
    print reaction, permutationTestResults[reaction];

all = []
for reaction in permutationTestResults:
    all.append(permutationTestResults[reaction])

all.sort(key=lambda t:t[4], reverse=True);

print "HITS:";
for t in all:
    if t[4] > 0 and t[4] < .05:
        print t

allPVals = []
for t in all:
    if t[4] > 0:
        allPVals.append(t[4]);

plt.hist(allPVals, bins=40);
plt.title("All P Values");
plt.show();
