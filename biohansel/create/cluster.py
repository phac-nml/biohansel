import attr


@attr.s
class Cluster(object):
    distance_matrix = attr.ib(default=None)
    clustering_array = attr.ib(default=None)
    hierarchical_clusters= attr.ib(default=None)

  