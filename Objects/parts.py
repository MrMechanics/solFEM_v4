#
#
#   parts.py
#  --------------
#   Part object, holding geometry and/-or mesh,
#   linking them together. Either from inported
#   step-file geometry, geometric shapes made 
#   from scratch with user input, imported mesh
#   from mesh files, or mesh built from scratch.
#

from geometries import *
from meshes import *


class Part:
    '''
Part object. Holds geometry and/or mesh for the
user to interact with.

    variables:
    --------------
    part.name
    
    part.mesh
    part.nodesets
    part.elementsets
    
    part.coordSys
    part.lines
    part.edges
    part.faces
    part.shell
    part.solid
    
    part.displayLists


'''
    def __init__(self, name, mesh=None):
        self.name = name

        self.mesh = mesh
        self.nodesets = {}
        self.elementsets = {}

        self.coordSys = CoordSys3D(Point3D(0.,0.,0.),Vector3D(1.,0.,0.),Vector3D(0.,1.,0.))
        self.lines = {}
        self.edges = {}
        self.faces = {}
        self.shell = {}
        self.solid = {}

        self.displayLists = {'nodes':       None,
                             'wireframe':   None,
                             'shaded':      None,
                             'lines':       None,
                             'faces':       None,
                             'seeds':       None}


    def readStepFile(self,filename):
        '''
    Import part geometry from step file.
    '''
        pass
    
    
    def importMesh(self,filename):
        '''
    Import mesh from *.sol file or from file
    generated in some other program. Supported
    file types are *.bdf (Nastran), 
    *.inp (ABAQUS) and *.dat (FreeCAD).
    '''
        pass