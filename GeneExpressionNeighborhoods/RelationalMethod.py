"""
RelationalMethod
================
A wrapper for numpy distance metric methods and similarity measure methods. The wrapper adds bestScore and worstScore
methods that return the nearest/most similar and the furthest/least similar values (respectively) in an array.
"""
import sys
sys.path.insert(0,'/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
from abc import ABCMeta
from abc import abstractproperty
import numpy as np
from scipy import stats
from scipy import spatial

class RelationalMethod:
    """
    Abstract Base Class for concrete RelationalMethod implementations that defines method, bestScore and worstScore.
    method is the distance or similarity method. bestScore is the method that returns the index of the nearest
    or most similar values in an array. worstScore is the method that returns the index of the furthest
    or least similar values in an array.
    """
    __metaclass__ = ABCMeta
    @abstractproperty
    def method(self):
        pass
    @abstractproperty
    def bestScore(self):
        pass
    @abstractproperty
    def worstScore(self):
        pass


class DistanceMethod:
    """
    A concrete class implementing RelationalMethod for Distance Methods. The bestScore is numpy.nanargmin and worstScore
    is numpy.nanargmax.

    Parameters:
    method: the distance method to wrap around
    """
    def __init__(self, method):
        self.method = method
        self.bestScore = np.nanargmin
        self.worstScore = np.nanargmax


class CorrelationMethod:
    """
    A concrete class implementing RelationalMethod for Correlation Methods. The bestScore is numpy.nanargmax of the
    absolute values of the array and worstScore is numpy.nanargmin, also of absolute values of the array.

    Parameters:
    method: the correlation method to wrap around
    """
    def __init__(self, method):
        self.method = method

    def bestScore(array):
        return np.nanargmax(np.absolute(array))

    def worstScore(array):
        return np.nanargmin(np.absolute(array))

# Register DistanceMethod and CorrelationMethod as implementations of RelationalMethod
RelationalMethod.register(DistanceMethod)
RelationalMethod.register(CorrelationMethod)

# Special Case Distance Methods
def customMahalonbis(vectorA, vectorB):
#     get inverted cov to pass to scipy mahalonbis
    covariance = np.cov(vectorA, vectorB)
    inverseCov = np.linalg.inv(covariance)
    return spatial.distance.mahalanobis(vectorA, vectorB, inverseCov)

# default map of relational methods
defaultRelationMethods = {
    "Manhattan": DistanceMethod(spatial.distance.cityblock),
    "Euclidean": DistanceMethod(spatial.distance.euclidean),
    # "Mahalonbis": DistanceMethod(customMahalonbis),
    "Cosine": DistanceMethod(spatial.distance.cosine),
    "Pearson": CorrelationMethod(stats.spearmanr),
    "Spearman": CorrelationMethod(stats.spearmanr)
}
