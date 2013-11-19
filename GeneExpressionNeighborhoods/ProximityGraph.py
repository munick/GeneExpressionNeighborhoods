"""
RelativeNighborGraph
=====================
Classes and methods defining graph types and converting similarity or distance matrices to sets of points and edges.
"""

import RelationalMethod
import numpy as np


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
    edgesWithWeights = set()
    edges = set()
    groups = []
    pointGroup = dict()
    for i, row in enumerate(relationMatrix):
        # maxIRow = getBestScore(row)
        print (i+1, row[i+1:])
        for j, relation in enumerate(row[i+1:]):
            # maxJRow = getBestScore(relationMatrix[j])
            if np.isnan(relation):
                isEdge = False

            else:
                isEdge = True # assume edge until proven wrong
                for r, iDistance in enumerate(row):
                    if getWorstScore([iDistance, relationMatrix[j, r], relation]) == 2:
                        isEdge = False # not an edge!
                        break

            if isEdge:
                edges.add(frozenset((i, j)))
                edgesWithWeights.add((frozenset((i, j)), relation))

                groupFound = False
                for gIndex, group in enumerate(groups):
                    if not group.isdisjoint({i, j}):
                        group.update({i, j})
                        groupFound = True
                        pointGroup.update({i: gIndex})
                        break
                if not groupFound:
                    pointGroup.update({i: len(groups)})
                    groups.append({i, j})

    return pointGroup, edges


    # for r in points:
    #     if max(dist2(p,r),dist2(q,r)) < dist2(p,q): return False
