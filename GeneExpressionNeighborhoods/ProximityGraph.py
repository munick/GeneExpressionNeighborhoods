"""
RelativeNighborGraph
=====================
Classes and methods defining graph types and converting similarity or distance matrices to sets of points and edges.
"""

import sys
sys.path.insert(0,'/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
import numpy as np
import math
import itertools
from scipy.spatial import Delaunay


# Delaunay triangulation
def getDelaunayTriangulationEdges(points):
    tri = Delaunay(points)
    edges = set()
    for triangle in tri.simplices:
        for edge in itertools.combinations(triangle, 2):
            edges.add(edge)
    return edges


# nearest neighbor graph
def getNearestNeighborGraph(relationMatrix, getBestScore):
    edgesWithWeights = set()
    edges = set()
    groups = []
    pointGroup = dict()
    for i, vector in enumerate(relationMatrix):
        nn = getBestScore(vector)

        if not np.isnan(nn):
            edges.add(frozenset((i, nn)))
            # print "%d" % nn
            # print "%d, %d:  %s" % (i, nn, type(relationMatrix[i, nn]))
            edgesWithWeights.add((frozenset((i, nn)), relationMatrix[i, nn]))
        else:
            nn = i

        groupFound = False
        for gIndex, group in enumerate(groups):
            if not group.isdisjoint({i, nn}):
                group.update({i, nn})
                groupFound = True
                pointGroup.update({i: gIndex})
                break
        if not groupFound:
            pointGroup.update({i: len(groups)})
            groups.append({i, nn})

    return pointGroup, edges


# relative neighbor graph
def getRelativeNeighborGraph(inputEdges, relationMatrix, getBestScore, getWorstScore):
    """
    Returns the relative neighborhood graph of the given relational matrix.

    :param relationMatrix: pairwise distance matrix
    :param getBestScore: returns highest similarity or smallest distance
    :param getWorstScore: returns lowest similarity or largest distance
    """
    edgesWithWeights = set()
    edges = set()
    groups = []
    pointGroup = dict()
    # loop through rows of distance/similarity matrix ************************************************* N
    for p, q in inputEdges:
        relationPQ = relationMatrix[p, q]
        row = relationMatrix[p]
        # maxJRow = getBestScore(relationMatrix[q])
        # non-numerical distances/similarities will not be counted as edges
        if np.isnan(relationPQ):
            isEdge = False
        # if there is a numeric value
        else:
            isEdge = True   # assume edge until proven wrong
            # loop through all columns in the ith row
            # relationPR is weight of edge p,r ***************************************************** (N^3)/2
            for r, relationPR in enumerate(row):
                # skip rows p and q and any points for which there is no distance value
                if p != r != q and (not np.isnan(relationPR)) and (not np.isnan(relationMatrix[q, r])):
                    # for triangle prq, if pq is the longest distance, then p and q are not neighbors
                    lengths = [relationPR, relationMatrix[q, r]]
                    if lengths[getWorstScore(lengths)] < relationPQ:
                        isEdge = False  # not an edge!
                        break           # break to next q
        # if p and q are neighbors
        if isEdge:
            edges.add(frozenset((p, q)))                            # add (p,q) tuple to edges set
            edgesWithWeights.add((frozenset((p, q)), relationPQ))   # add ((p,q), weight) to weighted edges set

    return pointGroup, edges


# Gabriel graph
def getGabrielNeighborGraph(inputEdges, relationMatrix, getBestScore, getWorstScore):
    """
    Returns the gabriel graph of the given relational matrix.

    :param relationMatrix: pairwise distance matrix
    :param getBestScore: returns highest similarity or smallest distance
    :param getWorstScore: returns lowest similarity or largest distance
    """
    edgesWithWeights = set()
    edges = set()
    groups = []
    pointGroup = dict()
    # loop through rows of distance/similarity matrix ************************************************* N
    for p, q in inputEdges:
        relationPQ = relationMatrix[p, q]
        row = relationMatrix[p]
        # maxJRow = getBestScore(relationMatrix[q])
        # non-numerical distances/similarities will not be counted as edges
        if np.isnan(relationPQ):
            isEdge = False
        # if there is a numeric value
        else:
            isEdge = True # assume edge until proven wrong
            # loop through all columns in the ith row
            # relationPR is weight of edge p,r ***************************************************** (N^3)/2
            for r, relationPR in enumerate(row):
                # skip rows p and q and any where there are no distance values
                if p != r != q and (not np.isnan(relationPR)) and (not np.isnan(relationMatrix[q, r])):
                    # if angle prq is > pi/2, then it's not an edge
                    # another calculation: if d(p,q) > sqrt of (d(p,r)^2 + d(r,q)^2)) then not an edge
                    lengths = math.sqrt(relationPR**2 + relationMatrix[q, r]**2)
                    if relationPQ > lengths:
                        isEdge = False  # not an edge!
                        break           # break to next q
        # if p and q are neighbors
        if isEdge:
            edges.add(frozenset((p, q)))                            # add (p,q) tuple to edges set
            edgesWithWeights.add((frozenset((p, q)), relationPQ))   # add ((p,q), weight) to weighted edges set

    return pointGroup, edges