#
#
#	solver.py
#  -----------
#
#	This is the solver module. It generates a FEmodel object, by reading
# 	data from the InputFile object, unless the FEmodel already exists.
#   It also has the functions necessary to build stiffness matrices, 
#   mass matrices and load vectors. After building the FEModel object, 
#   the requested solutions are calculated using the solutions objects.
#   Once solutions are finished calculated, results are stored in
#   *.res, *.out and *.csv files.
#

import os
import pickle
import numpy as np
import scipy.sparse as sp
import sys

sys.path.insert(1, '../Objects')

from timeit import time
from reader import *

from nodes import *
from materials import *
from sections import *
from elements import *
from meshes import *
from loads import *
from boundaries import *
from constraints import *
from dampings import *
from tables import *
from solutions import *






class Solver(object):
    '''
Base class for finite element model. It first instantiates
all nodes, materials, sections, elements, meshes, loads,
boundaries and solutions listed in the InputData object
if necessary. It assembles the global stiffness matrix and
mass matrix, and then runs the solutions as they are set up
by the user.
'''
    def __init__(self):
        self.name = ''
