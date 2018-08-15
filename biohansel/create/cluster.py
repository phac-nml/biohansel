import attr
import numpy as np

from typing import Dict

@attr.s
class Cluster(object):
    distance_matrix = attr.ib(default=None)
    clustering_array = attr.ib(default=None)
    flat_clusters= attr.ib(default=None)

  