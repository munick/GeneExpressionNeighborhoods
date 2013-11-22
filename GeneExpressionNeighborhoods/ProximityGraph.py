"""
RelativeNighborGraph
=====================
Classes and methods defining graph types and converting similarity or distance matrices to sets of points and edges.
"""

import RelationalMethod
import numpy as np
import math


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
def getRelativeNeighborGraph(relationMatrix, getBestScore, getWorstScore):
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
    for p, row in enumerate(relationMatrix):
        # maxIRow = getBestScore(row)
        print (p+1, row[p+1:])
        # loop from column p+1 through to the last column of ith row ********************************** (N^2)/2
        # essentially, we're looping through all weighted edges, let relationPQ = weight of edge p,q
        for q, relationPQ in enumerate(row[p+1:]):
            q += p + 1
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
                    # skip rows p and q
                    if p != r != q:
                        # for triangle prq, if pq is the longest distance, then p and q are not neighbors
                        lengths = [relationPR, relationMatrix[q, r]]
                        if lengths[getWorstScore(lengths)] < relationPQ:
                            isEdge = False  # not an edge!
                            break           # break to next q
            # if p and q are neighbors
            if isEdge:
                edges.add(frozenset((p, q)))                            # add (p,q) tuple to edges set
                edgesWithWeights.add((frozenset((p, q)), relationPQ))   # add ((p,q), weight) to weighted edges set
                # now, deal with groups ...
                groupFound = False                      # flag to keep track of weather a group was found or not
                # loop through list of groups
                for gIndex, group in enumerate(groups):
                    # if either p or q are in the group
                    if not group.isdisjoint({p, q}):
                        group.update({p, q})            # add p and q to the group
                        groupFound = True               # set groupFound flag
                        pointGroup.update({p: gIndex})  # mark gIndex to be the group for p
                        break
                # if no group is found for p and q
                if not groupFound:
                    pointGroup.update({p: len(groups)}) # mark new group index as index for p
                    groups.append({p, q})               # make new group with p, q

    return pointGroup, edges


# Gabriel graph
def getGabrielNeighborGraph(relationMatrix, getBestScore, getWorstScore):
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
    for p, row in enumerate(relationMatrix):
        # maxIRow = getBestScore(row)
        print (p+1, row[p+1:])
        # loop from column p+1 through to the last column of ith row ********************************** (N^2)/2
        # essentially, we're looping through all weighted edges, let relationPQ = weight of edge p,q
        for q, relationPQ in enumerate(row[p+1:]):
            q += p + 1
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
                    # skip rows p and q
                    if p != r != q:
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
                # now, deal with groups ...
                groupFound = False                      # flag to keep track of weather a group was found or not
                # loop through list of groups
                for gIndex, group in enumerate(groups):
                    # if either p or q are in the group
                    if not group.isdisjoint({p, q}):
                        group.update({p, q})            # add p and q to the group
                        groupFound = True               # set groupFound flag
                        pointGroup.update({p: gIndex})  # mark gIndex to be the group for p
                        break
                # if no group is found for p and q
                if not groupFound:
                    pointGroup.update({p: len(groups)}) # mark new group index as index for p
                    groups.append({p, q})               # make new group with p, q

    return pointGroup, edges