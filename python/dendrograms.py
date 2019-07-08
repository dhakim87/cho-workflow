import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
import pandas as pd
import seaborn as sns;
import scipy.stats;
import random;

baseDir = "../geneProteinReaction/modelResultsOld";
ubiquityScorePath = "../geneProteinReaction/ubiquityScore30.csv"

MIN_CELL_LINE = 101;
MAX_CELL_LINE = 196;

KEPT_OR_CORE = "kept";

#Part 1, Load up all of the model results, these files, one per cell line, contain secretion data, followed by all reactions kept in their cell specific models.
#Open all the files in modelResults
#For each file, write down the secretion micrograms of the cell line, and make a map from reactions to cell lines
secretionMicrograms = {}
reactionToCellLine = defaultdict(set);
cellLineToNumReactions = {}
for i in range(MIN_CELL_LINE, MAX_CELL_LINE):
    filename = baseDir + "/P6975_" + str(i) + "_"+KEPT_OR_CORE+".txt";
    cellLine = i - 101;
    with open(filename) as f:
        data = f.read().splitlines();
        micrograms = float(data[0].split(", ")[0])
        secretionMicrograms[cellLine] = micrograms;
        for val in data[1:]:
            reaction = val.strip();
            if reaction >= "r1664" and reaction <= "r2071":
                continue;
            reactionToCellLine[reaction].add(cellLine);
    
        cellLineToNumReactions[cellLine] = len(data)-1;
        print (cellLine, secretionMicrograms[cellLine], "#Reactions:", len(data)-1);

#Construct a numpy matrix that we can call linkage on with cell lines as rows and reactions as columns.  1s if the reaction is kept/core.
allReactions = []
for reaction in reactionToCellLine:
    allReactions.append(reaction);

allReactions.sort();
mat = np.zeros((95, len(allReactions)));

for reactionIndex in range(len(allReactions)):
    reaction = allReactions[reactionIndex];
    cellLines = reactionToCellLine[reaction];
    for cellLine in cellLines:
        mat[cellLine, reactionIndex] = 1;


#Part 2, load up the initial sample TPM data, filtered down to the rows that are shared across all of the cell lines (ie, not counting the TPM for the rows that are cell line specific).  Note that these tpm values won't sum to 1000000 per cell line, as they are each missing the TPM associated with their target protein.
totalTPM = [0] * 95;
rnaMat = np.zeros((95, 28476));
filename = "./sampleTPM.csv";
geneIDToRow = {}
lineNumber = 0;
with open(filename) as f:
    fileData = f.read().splitlines();
    for line in fileData[1:-95]:  #NOTICE:  We skip first line for header AND last 95 lines because this file format is stupid!
        tpms = line.split(",");
        for cellLine in range(95):
            totalTPM[cellLine] += float(tpms[cellLine+2]);
            rnaMat[cellLine, lineNumber] = float(tpms[cellLine+2])
            geneIDToRow[tpms[cellLine]] = lineNumber;
        lineNumber += 1;

#Part 2.5, load up the genes in the iCHOv1 model and filter down the TPM data to contain just these rows
filename = "./modelGenes.csv";
with open(filename) as f:
    fileData = f.read().splitlines();
    numFound = 0;
    for line in fileData:
        gene = line.strip();
        if gene in geneIDToRow:
            numFound += 1;
    rnaMat2 = np.zeros((95, numFound));
    col = 0;
    for line in fileData:
        gene = line.strip();
        print gene;
        if gene in geneIDToRow:
            foundCol = geneIDToRow[gene];
            for cl in range(95):
                rnaMat2[cl,col] = rnaMat[cl,foundCol];
            col += 1;

#Part 3, load up the ubiquity scores coming out of StanDep.  I suspect these should look (and cluster) extremely similar to the set of core reactions.

totalUbiquityScore = [0] * 95;
ubMat = np.zeros((95, 6663));
filename = ubiquityScorePath;
lineNumber = 0;
with open(filename) as f:
    fileData = f.read().splitlines();
    for line in fileData[1:]:  #Note, skip first line to remove header.
        ubScores = line.split(",");
        for cellLine in range(95):
            totalUbiquityScore[cellLine] += float(ubScores[cellLine]);
            ubMat[cellLine, lineNumber] = float(ubScores[cellLine]);
        lineNumber += 1;


from scipy.spatial import distance
matJacc = np.zeros((95,95))

for row in range(95):
    for col in range(95):
        matJacc[row,col] = distance.jaccard(mat[row,:], mat[col, :]);

linked = linkage(mat, metric='jaccard');
linked2 = linkage(rnaMat, metric='euclidean');
linked3 = linkage(rnaMat2, metric='euclidean');
linkedEuclideanOfJaccard = linkage(matJacc, metric='euclidean');
#linkedReactionByCellLine = linkage(mat.transpose(), metric='jaccard');

plt.subplot(2,2,1);
dendrogram(linked, orientation='top', distance_sort='descending', show_leaf_counts=True)
plt.title(KEPT_OR_CORE);
plt.subplot(2,2,2);
dendrogram(linked2, orientation='top', distance_sort='descending', show_leaf_counts=True)
plt.title("RNA")
plt.subplot(2,2,3);
dendrogram(linked3, orientation='top', distance_sort='descending', show_leaf_counts=True)
plt.title("RNA (Genes In Model Only)")
plt.subplot(2,2,4);
dendrogram(linkedEuclideanOfJaccard, orientation='top', distance_sort='descending', show_leaf_counts=True)
plt.title("Euclidean Distance Of Jaccard")
#plt.subplot(2,2,4);
#dendrogram(linkedReactionByCellLine, orientation='top', distance_sort='descending', show_leaf_counts=True)
#plt.title("Clustered Reactions")

plt.show()

cluster = AgglomerativeClustering(n_clusters=2, affinity='euclidean')
cluster.fit_predict(mat)

for cellLine in range(95):
    print cellLine, cluster.labels_[cellLine];

secretionList = []
reactionCountList = []
g1ReactionCounts = []
g2ReactionCounts = []
g1Secretion = []
g2Secretion = []
g1TPMNonTarget = []
g2TPMNonTarget = []
g1TotalUB = []
g2TotalUB = []

for cellLine in range(95):
    secretionList.append(secretionMicrograms[cellLine]);
    reactionCountList.append(cellLineToNumReactions[cellLine]);
    if cluster.labels_[cellLine] == 0:
        g1ReactionCounts.append(cellLineToNumReactions[cellLine]);
        g1Secretion.append(secretionMicrograms[cellLine]);
        g1TPMNonTarget.append(totalTPM[cellLine])
        g1TotalUB.append(totalUbiquityScore[cellLine]);
    else:
        g2ReactionCounts.append(cellLineToNumReactions[cellLine]);
        g2Secretion.append(secretionMicrograms[cellLine]);
        g2TPMNonTarget.append(totalTPM[cellLine]);
        g2TotalUB.append(totalUbiquityScore[cellLine]);

secretionMWU = scipy.stats.mannwhitneyu(g1Secretion, g2Secretion)
reactionCountMWU = scipy.stats.mannwhitneyu(g1ReactionCounts, g2ReactionCounts)
tpmNonTargetMWU = scipy.stats.mannwhitneyu(g1TPMNonTarget, g2TPMNonTarget)
totalUbiquityScoreMWU = scipy.stats.mannwhitneyu(g1TotalUB, g2TotalUB)

data = pd.DataFrame({'secretion':secretionList, 'numReactions':reactionCountList, 'dendroGroup':cluster.labels_, 'totalTPM': totalTPM, 'totalUB':totalUbiquityScore});

plt.subplot(4,1,1);
sns.swarmplot(x='dendroGroup', y='secretion', data=data);
plt.title("P = " + str(secretionMWU.pvalue));
plt.subplot(4,1,2);
sns.swarmplot(x='dendroGroup', y='totalTPM', data=data);
plt.title("P = " + str(tpmNonTargetMWU.pvalue));
plt.subplot(4,1,3);
sns.swarmplot(x='dendroGroup', y='numReactions', data=data);
plt.title("P = " + str(reactionCountMWU.pvalue));
plt.subplot(4,1,4);
sns.swarmplot(x='dendroGroup', y='totalUB', data=data);
plt.title("P = " + str(totalUbiquityScoreMWU.pvalue));

plt.show();

#Identify reactions that tend towards one group or the other
X = [];
Y = [];
C = [];

funkyReactions = []

for reaction in range(len(allReactions)):
    aCount = 0;
    bCount = 0;
    for cellLine in range(95):
        if mat[cellLine, reaction] == 1:
            if cluster.labels_[cellLine] == 0:
                aCount += 1;
            else:
                bCount += 1;
    X.append(aCount);
    Y.append(bCount);
    if aCount == 0 and bCount > 40:
        C.append('r')
        funkyReactions.append(reaction);
        print(allReactions[reaction]);
    else:
        C.append('k')

for i in range(len(X)):
    X[i] += random.random()-.5;
    Y[i] += random.random()-.5;
    
plt.scatter(X,Y, color=C);
plt.show();

#Looks like the culprit in the bimodal distribution is 200+ reactions that are all based on the same 2-3 genes.  So we're going to make a swarm plot of the tpm values of these genes

col1 = geneIDToRow["100768667"]
col2 = geneIDToRow["100760403"]
col3 = geneIDToRow["100689315"]

tpm1 = rnaMat[:, col1]
tpm2 = rnaMat[:, col2]
tpm3 = rnaMat[:, col3]

#100768667 is the basis of reactions 5886 - 5911
#100760403 is the basis of reactions 5912 - 6286
#100689315 is the basis of reactions 5886 - 6286
#Suspect two subgroups of ubiquity scores based on these ranges.

ub1 = ubMat[:, 5889]
ub2 = ubMat[:, 5900] #Should be identical to ub1
ub3 = ubMat[:, 5913]
ub4 = ubMat[:, 6000] #Should be identical to ub3

targetGeneTPM = pd.DataFrame({
    '100768667 (tpm)':tpm1, '5889 (ub)':ub1,
    '100760403 (tpm)':tpm2, '5900 (ub)':ub2,
    '100689315 (tpm)':tpm3, '5913 (ub)':ub3, '6000 (ub)':ub4,
    "cluster":cluster.labels_
});

plt.subplot(3,3,1);
sns.swarmplot(x='cluster', y='100768667 (tpm)', data=targetGeneTPM);
plt.subplot(3,3,2);
plt.hist(tpm1, bins=15);
plt.subplot(3,3,3);
sns.swarmplot(x='cluster', y='5900 (ub)', data=targetGeneTPM);

plt.subplot(3,3,4);
sns.swarmplot(x='cluster', y='100760403 (tpm)', data=targetGeneTPM);
plt.subplot(3,3,5);
plt.hist(tpm2, bins=15);
plt.subplot(3,3,6);
sns.swarmplot(x='cluster', y='6000 (ub)', data=targetGeneTPM);

plt.subplot(3,3,7);
sns.swarmplot(x='cluster', y='100689315 (tpm)', data=targetGeneTPM);
plt.subplot(3,3,8);
plt.hist(tpm3, bins=15);
plt.show();

plt.show();
