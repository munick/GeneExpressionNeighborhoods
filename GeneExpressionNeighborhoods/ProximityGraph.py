"""
RelativeNighborGraph
=====================
Classes and methods defining graph types and converting similarity or distance matrices to sets of points and edges.
"""

import sys
sys.path.insert(0,'/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
import numpy as np
import math
import multiprocessing

_globalRelationalMatrix = np.zeros((1, 2))


# nearest neighbor graph
def getNearestNeighborGraph(relationMatrix, getBestScore):
    edges = set()
    for i, relationArray in enumerate(relationMatrix):
        try:
            nn = getBestScore(relationArray)
            if not np.isnan(nn):
                edges.add(frozenset((i, nn)))
        except ValueError:  # numpy 1.8 raises exception instead of warn & null value!
            continue

    return getPointGroupMapping(edges), edges


# relative neighbor graph
def getRelativeNeighborGraph(inputEdges, relationMatrix, getBestScore, getWorstScore):
    """
    Returns the relative neighborhood graph of the given relational matrix.

    :param relationMatrix: pairwise distance matrix
    :param getBestScore: returns highest similarity or smallest distance
    :param getWorstScore: returns lowest similarity or largest distance
    """
    edges = set()
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

    return getPointGroupMapping(edges), edges


# relative neighbor graph
def getRelativeNeighborGraphMP(inputEdges, relationMatrix, getBestScore, getWorstScore):
    """
    Returns the relative neighborhood graph of the given relational matrix.

    :param relationMatrix: pairwise distance matrix
    :param getBestScore: returns highest similarity or smallest distance
    :param getWorstScore: returns lowest similarity or largest distance
    """
    global _globalRelationalMatrix
    _globalRelationalMatrix = relationMatrix                                # put relationMatrix in global mem
    inputEdges = list(inputEdges)
    edges = set()                                                           # set of edges found
    queue = multiprocessing.Queue()                                         # results queue
    numberOfProcesses = min(len(inputEdges), multiprocessing.cpu_count())   # number of processors to use
    edgesPerProcess = len(inputEdges)/numberOfProcesses                     # number of edges per processor
    inputEdgeSubsets = []                                                   # contains edge set for each process
    subsetIndex = -1
    print "Dividing Edges"
    # divide inputEdges into even subsets
    for i, edge in enumerate(inputEdges):
        if i % edgesPerProcess == 0:
            inputEdgeSubsets.append(set())
            subsetIndex += 1
        inputEdgeSubsets[subsetIndex].add(edge)

    print "Running Processes: ", numberOfProcesses
    # run processes
    for i in range(numberOfProcesses):
        p = multiprocessing.Process(target=rngWorker, args=(inputEdgeSubsets[i], queue))
        p.start()

    # collect results
    for i in range(numberOfProcesses):
        edges = edges.union(queue.get())
        print i, " process done"

    return getPointGroupMapping(edges), edges


def rngWorker(inputEdges, queue):
    """
    Work done by each processor
    :param inputEdges: set of edges (p, q)
    :param queue: shared Queue to place results
    """
    edges = set()
    for p, q in inputEdges:
        relationPQ = _globalRelationalMatrix[p, q]
        row = _globalRelationalMatrix[p]
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
                if p != r != q and (not np.isnan(relationPR)) and (not np.isnan(_globalRelationalMatrix[q, r])):
                    # for triangle prq, if pq is the longest distance, then p and q are not neighbors
                    lengths = [relationPR, _globalRelationalMatrix[q, r]]
                    if lengths[np.nanargmax(lengths)] < relationPQ:
                        isEdge = False      # not an edge!
                        break               # break to next q
        # if p and q are neighbors
        if isEdge:
            edges.add(frozenset((p, q)))    # add (p,q) tuple to edges set
    queue.put(edges)


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
    # loop through rows of distance/similarity matrix ************************************************* N
    for p, q in inputEdges:
        relationPQ = relationMatrix[p, q]
        row = relationMatrix[p]
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

    return getPointGroupMapping(edges), edges


# Identify Clusters/Neighborhoods in edges
def getPointGroupMapping(edges):
    pointGroupMapping = dict()
    pointGroups = dict()
    groupIndex = 0
    # loop through all edges
    for edge in edges:
        matchingGroupIndexes = set()
        # find groups they belong to
        for index, group in pointGroups.iteritems():
            if not edge.isdisjoint(group):
                # they have points in common!
                matchingGroupIndexes.add(index)
        # if they belong to one or more groups
        if len(matchingGroupIndexes) > 0:
            # take one index out of matchingGroupIndexes
            indexToKeep = matchingGroupIndexes.pop()
            # if matchingGroupIndexes has any more indexes
            for indexToRemove in matchingGroupIndexes:
                # merge them into indexToKeep and remove them from pointGroups
                pointGroups[indexToKeep].update(pointGroups.pop(indexToRemove))
            # finally, add edges to indexToKeep pointGroup
            pointGroups[indexToKeep].update(edge)

        # if they don't belong to an existing group
        if len(matchingGroupIndexes) == 0:
            # add edges to new group
            pointGroups[groupIndex] = set(edge)
            groupIndex += 1

        # create point -> groupId mapping
        for groupId, points in enumerate(pointGroups.values()):
            for point in points:
                pointGroupMapping[point] = groupId
    return pointGroupMapping


