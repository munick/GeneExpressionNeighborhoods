"""
GENController: Gene Expression Neighborhood Controller
======================================================
An evaluation of gene expression clustering via proximity graphs.

GENController is called with the following arguments (in order):
similarityMeasure - any measure defined in RelationalMethod.py (e.g. Manhattan, Euclidean, Pearson, etc.)
geneExpressionFileName - path to the input file containing gene expression levels, note the basepath is assumed to be data
geneNameFileName - path to the input file containing gene names, note the basepath is assumed to be data
"""

from RelationalMethod import defaultRelationMethods as relationMethods
import ProximityGraph as pg
import csv
import sys
import numpy as np
from scipy import stats


def calculateRelationMatrix(dataset, relationMethod):
    """
    Given a dataset and a distance method, computes a 2d matrix of distances between each of the rows in the dataset

    Parameters
    ----------
    :param dataset: array_like
        a 2d numpy array of normalized gene expressions.
    :param relationMethod: method
        any method that will calculate the distance or correlation between two 1d arrays.

    Returns
    ----------
    :returns: 2darray
        2d matrix of distances or correlations between each of the rows in the dataset
    """
    relationMatrix = np.zeros((len(dataset), len(dataset)))
    # loop through each row and then each row below it
    enumerated = list(enumerate(dataset))
    for i, row in enumerated:
        for j, comparisonRow in enumerated[i+1:]:
            distance = relationMethod(row, comparisonRow)
            if isinstance(distance, tuple):
                # from what I can tell, distance is always first and p value second
                distance = distance[0]
            relationMatrix[i, j] = relationMatrix[j, i] = distance
    return relationMatrix


# export edges to csv
def exportEdgesToCSV(graphType, relationMethod, edgeData):
    with open('output/yeast/' + relationMethod + graphType + 'Edges.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(["source","target"])
        for edge in edgeData:
            writer.writerow(list(edge))


# exports gene names and groups to csv
def exportGeneGroupToCSV(graphType, relationMethod, geneGroupData, geneNames):
    with open('output/yeast/' + relationMethod + graphType + 'Groups.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(["name","group"])
        for i, gene in enumerate(geneNames):
            writer.writerow([gene, geneGroupData.get(i)])


# Get distance method from command line arguments
if len(sys.argv) >= 4:
    distanceMethod = sys.argv[1]
    geneExpressionFileName = 'data/' + sys.argv[2]
    geneNameFileName = 'data/' + sys.argv[3]
else:
    # default to Euclidean distance
    if len(sys.argv) == 1:
        distanceMethod = "Euclidean"
    else:
        distanceMethod = sys.argv[1]
    # default to yeast gene expression dataset
    geneExpressionFileName = "data/yeast/yeastEx.txt"
    geneNameFileName = "data/yeast/yeastNames.txt"


# Read gene expression data and gene names files
data  = np.genfromtxt(geneExpressionFileName, dtype=np.uint16)
genes = np.genfromtxt(geneNameFileName, dtype=("|S10"))

# normalize the data using zscore
normalizedData = stats.zscore(data, 1)

# replace nans with infinite values
whereAreNaNs = np.isnan(normalizedData)
normalizedData[whereAreNaNs] = np.Inf

# np.savetxt("normalized.txt", normalizedData, '%.2f')
print "Calculating relationship with %s method ..." % distanceMethod
relationMatrix = calculateRelationMatrix(normalizedData, relationMethods[distanceMethod].method);
print "... %s calculation complete" % distanceMethod

# replace infinite values with NaNs
relationMatrix[np.isinf(relationMatrix)] = np.NaN

# replace diagonal with Nans
np.fill_diagonal(relationMatrix, np.NaN)

## Nearest Neighbor
nnPointGroup, nnEdges = pg.getNearestNeighborGraph(relationMatrix, relationMethods[distanceMethod].bestScore)
## Relative Neighbor
# nnPointGroup, nnEdges = getRelativeNeighborGraph(relationMatrix, relationMethods[distanceMethod].bestScore, relationMethods[distanceMethod].worstScore)

# export to csv
exportEdgesToCSV("NN", distanceMethod, nnEdges)
exportGeneGroupToCSV("NN", distanceMethod, nnPointGroup, genes)


