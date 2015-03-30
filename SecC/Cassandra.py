# cython: boundscheck=False, c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, cdivision_warnings=True, profile=True

from cassandra.cluster import Cluster

import cPickle

"""
Contains utilities relating to working with a kvp databases.
"""


def StartCassandra(keyspace="SecC", IP="127.0.0.1"):
    """
    Starts a Cassandra session for a given database.
    """
    cluster = Cluster(IP)
    session = cluster.connect(keyspace)
    return session

