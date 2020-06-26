"""
openmc.thsolver
==============

A thermal/hydraulic feedback front-end tool.
"""

from .dummy_comm import DummyCommunicator
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    have_mpi = True
except ImportError:
    comm = DummyCommunicator()
    have_mpi = False

from .constants import *
from .variables import * 
from .init      import * 
from .materials import * 
from .solver    import * 
 
