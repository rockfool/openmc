"""
openmc.fuel
==============

A one-dimension fuel performance feedback front-end tool.
"""

from .dummy_comm import DummyCommunicator
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    have_mpi = True
except ImportError:
    comm = DummyCommunicator()
    have_mpi = False

from .variables import *  
from .materials import * 
from .operator  import * 
 
 
