"""
GENController: Gene Expression Neighborhood Controller
======================================================
An evaluation of gene expression clustering via proximity graphs.

GENController is called with the following arguments (in order):
similarityMeasure - any measure defined in RelationalMethod.py (e.g. Manhattan, Euclidean, Pearson, etc.)
geneExpressionFileName - path to the input file containing gene expression levels, note the basepath is assumed to be data
geneNameFileName - path to the input file containing gene names, note the basepath is assumed to be data
"""

import sys
sys.path.insert(0,'/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
from RelationalMethod import defaultRelationMethods as relationMethods
import ProximityGraph as pg
import csv
import os
import itertools
import time
import numpy as np
from scipy import stats
from scipy import spatial

baseInputDir = "data/"
baseOutputDir = "output/"


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
def exportEdgesToCSV(subdirectory, graphType, relationMethod, edgeData):
    with open(subdirectory + '/' + relationMethod + graphType + 'Edges.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(["source","target"])
        for edge in edgeData:
            writer.writerow(list(edge))


# exports gene names and groups to csv
def exportGeneGroupToCSV(subdirectory, graphType, relationMethod, geneGroupData, geneNames):
    with open(subdirectory + '/' + relationMethod + graphType + 'Groups.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(["name","group"])
        for i, gene in enumerate(geneNames):
            writer.writerow([gene, geneGroupData.get(i)])


# export gene names, groups and edges to csv
def exportToCSV(subdirectory, graphType, relationMethod, edgeData, geneGroupData, geneNames):
    # create output dir if necessary
    outputDirectory = subdirectory.replace(baseInputDir, baseOutputDir, 1)
    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)
    # export to csv
    exportEdgesToCSV(outputDirectory, graphType, relationMethod, edgeData)
    exportGeneGroupToCSV(outputDirectory, graphType, relationMethod, geneGroupData, geneNames)


if __name__ == "__main__":
    # As of now this only runs on posix systems (mainly, osx and *nix operating systems)
    if os.name != "posix":
        print "GeneExpressionNeighborhoods can only be run on posix systems. Try osx or linux."
        sys.exit()

    # Get distance method from command line arguments
    if len(sys.argv) >= 4:
        distanceMethod = sys.argv[1]
        geneExpressionFileName = baseInputDir + sys.argv[2]
        geneNameFileName = baseInputDir + sys.argv[3]
        if len(sys.argv) >= 5:
            skipNormalization = sys.argv[4] == "skip"
        else:
            skipNormalization = False
    else:
        # default to Euclidean distance
        if len(sys.argv) == 1:
            distanceMethod = "Euclidean"
        else:
            distanceMethod = sys.argv[1]
        # default to yeast gene expression dataset
        geneExpressionFileName = baseInputDir + "yeast/yeastEx.txt"
        geneNameFileName = baseInputDir + "yeast/yeastNames.txt"
        # normalize by default
        skipNormalization = False

    start_time = time.time()
    print start_time, " (start time)"
    # get subdirectory
    subdirectory = os.path.dirname(geneExpressionFileName)

    # INPUT
    # Read gene expression data and gene names files
    data  = np.genfromtxt(geneExpressionFileName, dtype=np.float32)
    genes = np.genfromtxt(geneNameFileName, dtype=("|S10"))

    #PREPROCESSING
    if skipNormalization:
        print "Skipping Normalization/Pre-Processing"
        normalizedData = data
    else:
        print "Pre-Processing data ..."
        # normalize the data using zscore
        normalizedData = stats.zscore(data, 1)

        # replace nans with infinite values
        whereAreNaNs = np.isnan(normalizedData)
        normalizedData[whereAreNaNs] = np.Inf

        #multiply by 100
        normalizedData = normalizedData * 100
        print "... Pre-Processing complete"

    # np.savetxt("normalized.txt", normalizedData, '%.2f')
    print "Calculating relationship with %s method ..." % distanceMethod
    relationMatrix = calculateRelationMatrix(normalizedData, relationMethods[distanceMethod].method)
    print "... %s calculation complete" % distanceMethod
    done_with_distance_time = time.time()
    print done_with_distance_time, " (end distance matrix time)"
    print done_with_distance_time - start_time, " seconds (total runtime for distance matrix)"

    # replace infinite values with NaNs
    relationMatrix[np.isinf(relationMatrix)] = np.NaN

    # replace diagonal with Nans
    np.fill_diagonal(relationMatrix, np.NaN)

    # Nearest Neighbor
    print "Calculating NN ..."
    nnPointGroupMapping, nnEdges = pg.getNearestNeighborGraph(relationMatrix, relationMethods[distanceMethod].bestScore)
    print "... NN calculation complete"
    # OUTPUT
    exportToCSV(subdirectory,"NN", distanceMethod, nnEdges, nnPointGroupMapping, genes)
    done_with_nn_time = time.time()
    print done_with_nn_time, " (end NN time)"
    print done_with_nn_time - done_with_distance_time, " seconds (runtime for NN)"

    ## Relative Neighbor Multi Processor
    allEdges = itertools.combinations(range(len(data)), 2)
    print "Calculating RN_MP ..."
    rnPointGroupMapping, rnEdges = pg.getRelativeNeighborGraphMP(allEdges, relationMatrix, relationMethods[distanceMethod].bestScore, relationMethods[distanceMethod].worstScore)
    print "... RN_MP calculation complete"
    # OUTPUT
    exportToCSV(subdirectory, "RN_MP", distanceMethod, rnEdges, rnPointGroupMapping, genes)
    done_with_rn_time = time.time()
    print done_with_rn_time, " (end RN time)"
    print done_with_rn_time - done_with_nn_time, " seconds (runtime for RN)"

    # Relative Neighbor
    print "Calculating GG ..."
    ggPointGroupMapping, ggEdges = pg.getGabrielNeighborGraph(rnEdges, relationMatrix, relationMethods[distanceMethod].bestScore, relationMethods[distanceMethod].worstScore)
    print "... GG calculation complete"
    # OUTPUT
    exportToCSV(subdirectory, "GG", distanceMethod, ggEdges, ggPointGroupMapping, genes)
    end_time = time.time()
    print end_time, " (end time)"
    print end_time - done_with_rn_time, " seconds (runtime for GG)"
    print end_time - start_time, " seconds (total runtime)"
