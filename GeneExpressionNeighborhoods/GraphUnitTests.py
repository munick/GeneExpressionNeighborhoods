import sys
sys.path.insert(0,'/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
import unittest
import time
import itertools
import ProximityGraph as pg
import numpy as np
from scipy import spatial


class TestRelativeNeighborhoodGraph(unittest.TestCase):
    def setUp(self):
        # import data files and edge files
        self.points = np.genfromtxt("data/testRNGWiki/points.csv")
        self.genes  = np.genfromtxt("data/testRNGWiki/yeastNames.txt")
        edges = np.genfromtxt("data/testRNGWiki/edges.csv", delimiter=",")
        self.expectedEdges = set()
        for edge in edges:
            self.expectedEdges.add(frozenset(edge))
        condensedDistanceMatrix = spatial.distance.pdist(self.points,metric='euclidean')
        distanceMatrix = spatial.distance.squareform(condensedDistanceMatrix)
        self.allEdges = itertools.combinations(range(len(self.points)), 2)
        self.pointGroups, self.edges = pg.getRelativeNeighborGraph(self.allEdges, distanceMatrix, np.nanargmin, np.nanargmax)

    def testEdges(self):
        self.assertIsNotNone(self.edges, "No edges were returned")
        unexpectedEdges = self.edges.difference(self.expectedEdges)
        self.assertTrue(len(unexpectedEdges) == 0, "" + str(len(unexpectedEdges)) +
                                                   " edges were found in the RNG that were not in the expected edge set")
        missingEdges = self.expectedEdges.difference(self.edges)
        self.assertTrue(len(missingEdges) == 0, "" + str(len(missingEdges)) +
                                                " edges were expected in the RNG but not found")

    def testGroups(self):
        self.assertIsNotNone(self.pointGroups, "No pointGroups were returned")
        # todo: add unit tests and expected results for pointGroups


class TestGabrielGraph(unittest.TestCase):
    def setUp(self):
        # import data files and edge files
        self.points = np.genfromtxt("data/testGGWiki/GGpoints.txt")
        self.genes  = np.genfromtxt("data/testGGWiki/yeastNames.txt")
        edges = np.genfromtxt("data/testGGWiki/GGedges.csv", delimiter=",")
        self.expectedEdges = set()
        for edge in edges:
            self.expectedEdges.add(frozenset(edge))
        condensedDistanceMatrix = spatial.distance.pdist(self.points,metric='euclidean')
        self.distanceMatrix = spatial.distance.squareform(condensedDistanceMatrix)
        self.allEdges = itertools.combinations(range(len(self.points)), 2)
        self.pointGroups, self.edges = pg.getGabrielNeighborGraph(self.allEdges, self.distanceMatrix, np.nanargmin, np.nanargmax)

    def testEdges(self):
        self.assertIsNotNone(self.edges, "No edges were returned in the GG")
        unexpectedEdges = self.edges.difference(self.expectedEdges)
        self.assertTrue(len(unexpectedEdges) == 0, "" + str(len(unexpectedEdges)) +
                                                   " edges were found in the GG that were not in the expected edge set")
        missingEdges = self.expectedEdges.difference(self.edges)
        self.assertTrue(len(missingEdges) == 0, "" + str(len(missingEdges)) +
                                                " edges were expected in the GG but not found")


class TestRNGGGSubsetness(unittest.TestCase):
    def setUp(self):
        # import data files and edge files
        # ggEdges = np.genfromtxt("data/testGGWiki/GGedges.csv", delimiter=",")
        ggEdges = np.genfromtxt("output/yeast/EuclideanGGEdges.csv", delimiter=",")
        self.expectedGGEdges = set()
        for edge in ggEdges:
            self.expectedGGEdges.add(frozenset(edge))

        # rngEdges = np.genfromtxt("data/testRNGWiki/edges.csv", delimiter=",")
        rngEdges = np.genfromtxt("output/yeast/EuclideanRNEdges.csv", delimiter=",")
        self.expectedRNGEdges = set()
        for edge in rngEdges:
            self.expectedRNGEdges.add(frozenset(edge))

    def testSubsetness(self):
        # there shouldnt be any edges inGG that are not in RNG
        edgesInRNGButNotGG = self.expectedGGEdges.difference(self.expectedRNGEdges)
        edgesInGGButNotRNG = self.expectedRNGEdges.difference(self.expectedGGEdges)
        self.assertTrue(len(edgesInGGButNotRNG) == 0,
                        "There shouldn't be any edges in GG that are not in RNG, but there are " +
                        str(len(edgesInGGButNotRNG)) + "!")
        # self.assertTrue(len(edgesInRNGButNotGG) == 0,
        #                 " " +
        #                 str(len(edgesInRNGButNotGG)) + "!")

class TestDTRNGGGSubsetness(unittest.TestCase):
    def setUp(self):
        # import data files and edge files
        points = np.genfromtxt("data/yeast/yeastEx.txt")
        # ggEdges = np.genfromtxt("data/testGGWiki/GGedges.csv", delimiter=",")
        ggEdges = np.genfromtxt("output/yeast_11-23/EuclideanGGEdges.csv", delimiter=",")
        self.expectedGGEdges = set()
        for edge in ggEdges:
            self.expectedGGEdges.add(frozenset(edge))

        # rngEdges = np.genfromtxt("data/testRNGWiki/edges.csv", delimiter=",")
        rngEdges = np.genfromtxt("output/yeast_11-23/EuclideanRNEdges.csv", delimiter=",")
        self.expectedRNGEdges = set()
        for edge in rngEdges:
            self.expectedRNGEdges.add(frozenset(edge))
        start_time = time.time()
        print start_time, "start time"
        dtEdges = pg.getDelaunayTriangulationEdges(points)
        run_time = time.time() - start_time
        print run_time, "seconds"
        self.expectedDTEdges = dtEdges

    def testSubsetness(self):
        # there shouldnt be any edges inGG that are not in RNG
        self.assertTrue(self.expectedDTEdges >= self.expectedRNGEdges >= self.expectedGGEdges, "DT !>= RNG !>= GG")
        edgesInRNGButNotGG = self.expectedGGEdges.difference(self.expectedRNGEdges)
        edgesInGGButNotRNG = self.expectedRNGEdges.difference(self.expectedGGEdges)
        self.assertTrue(len(edgesInGGButNotRNG) == 0,
                        "There shouldn't be any edges in GG that are not in RNG, but there are " +
                        str(len(edgesInGGButNotRNG)) + "!")
        # self.assertTrue(len(edgesInRNGButNotGG) == 0,
        #                 " " +
        #                 str(len(edgesInRNGButNotGG)) + "!")



if __name__ == '__main__':
    unittest.main()