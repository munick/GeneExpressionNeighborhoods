import unittest
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
        self.pointGroups, self.edges = pg.getRelativeNeighborGraph(distanceMatrix, np.nanargmin, np.nanargmax)

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
        self.pointGroups, self.edges = pg.getGabrielNeighborGraph(self.distanceMatrix, np.nanargmin, np.nanargmax)

    def testEdges(self):
        self.assertIsNotNone(self.edges, "No edges were returned in the GG")
        unexpectedEdges = self.edges.difference(self.expectedEdges)
        self.assertTrue(len(unexpectedEdges) == 0, "" + str(len(unexpectedEdges)) +
                                                   " edges were found in the GG that were not in the expected edge set")
        missingEdges = self.expectedEdges.difference(self.edges)
        self.assertTrue(len(missingEdges) == 0, "" + str(len(missingEdges)) +
                                                " edges were expected in the GG but not found")


if __name__ == '__main__':
    unittest.main()