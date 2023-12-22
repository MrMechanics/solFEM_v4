#
#
#   modeler.py
#  -----------------
#   This is the modeler module. This is where all
#   individual parts are combined into one large
#   Finite Element model. It is also where the
#   Finite Element model is manipulated to apply
#   constraints, loads, boundary conditions, etc.,
#   for setting up an analysis.
#   
#   It lets you send the model directly to the
#   solver, or to generate a *.sol file which can
#   also be run in the solver. It also imports
#   mesh and some other data from *.sol, *.dat, 
#   *.mdl, *.inp, and *.bdf files, as well as
#   geometry from *.step files, and results from
#   *.out files.
#

from mesher import *
from parts import *
from meshes import *
from elements import *
from nodes import *
from loads import *
from boundaries import *
from constraints import *




class FEmodel(object):
    '''
Class describing the Finite Element model.
Combining all parts with materials, sections,
constraints, loads, boundary conditions, etc.,
for the Finite Element Analysis.
'''
    def __init__(self, interface):

        self.interface = interface
        self.mesher = Mesher()
        
        self.materials = {'6061-T6_aluminum (m kg N)':          {'Elasticity': 689e8, 'Poisson ratio': 0.35, 'Density':   2700.},
                          '1010_carbon_steel (m kg N)':         {'Elasticity': 205e9, 'Poisson ratio': 0.29, 'Density':   7870.},
                          '316_stainless_steel (m kg N)':       {'Elasticity': 193e9, 'Poisson ratio': 0.27, 'Density':   7870.},
                          'Grade2_titanium (m kg N)':           {'Elasticity': 105e9, 'Poisson ratio': 0.37, 'Density':   4510.},
                          '6061-T6_aluminum (mm kg mN kPa)':    {'Elasticity': 689e5, 'Poisson ratio': 0.35, 'Density':  2.7e-6},
                          '1010_carbon_steel (mm kg mN kPa)':   {'Elasticity': 205e6, 'Poisson ratio': 0.29, 'Density': 7.87e-6},
                          '316_stainless_steel (mm kg mN kPa)': {'Elasticity': 193e6, 'Poisson ratio': 0.27, 'Density': 7.87e-6},
                          'Grade2_titanium (mm kg mN kPa)':     {'Elasticity': 105e6, 'Poisson ratio': 0.37, 'Density': 4.51e-6}}
        self.sections = {}
        self.parts = {}

        self.currentPart = None
        self.currentSolution = None
        self.currentResults = {'solution': None, 'result': None, 'subresult': None}

        self.nodesSelected = False
        self.selected_nodes = {}
        self.elementsSelected = False
        self.selected_elements = {}
        self.linesSelected = False
        self.selected_lines = {}
        self.facesSelected = False
        self.selected_faces = {}
        self.selectOption = 'lines'
        
        self.results = {}
        self.scaleShearBendDiagram = 1.
        self.scale_factor = 20.
        self.colors = {'selected':       (1.00, 0.00, 0.00, 1.0),
                       'loads':          (0.6875, 0.3984375, 0.375, 1.0),
                       'displacements':  (0.4140625, 0.48828125, 0.5546875, 1.0),
                       'constraints':    (0.7890625, 0.55859375, 0.2578125, 1.0)}
        self.displayLists = {'selected_nodes': None,
                             'selected_elements': None,
                             'selected_lines': None,
                             'selected_faces': None,
                             'solutions': {}}




    def clearModel(self):
        '''
    Clears out all data in model from
    current session.
    '''
        pass
    
    

    def setCurrentPart(self,part_name):
        '''
    Sets the current Part to part_name so the FEmodel knows what
    part is currently active, and the viewer knows what part to
    render if not viewing assembly.
    '''
        self.currentPart = self.parts[part_name]
        self.interface.viewer.updateDisplayList()
        
           

    def selectedFeaturesDisplayList(self):
        '''
    Create a displaylist for the currently selected
    lines, faces, nodes or elements.
    '''
        if self.selectOption == 'lines':
            if len(self.selected_lines) != 0:
                self.linesSelected = True
                self.displayLists['selected_lines'] = glGenLists(1)
                glNewList(self.displayLists['selected_lines'], GL_COMPILE)
                glLineWidth(5.0)
                glColor3f(self.colors['selected'][0],
                          self.colors['selected'][1],
                          self.colors['selected'][2])
                for line in self.selected_lines:
                    for point in range(len(self.selected_lines[line].points)-1):
                        glBegin(GL_LINES)
                        coord = self.selected_lines[line].points[point]
                        glVertex3f(coord.x(),coord.y(),coord.z())
                        coord = self.selected_lines[line].points[point+1]
                        glVertex3f(coord.x(),coord.y(),coord.z())
                        glEnd()
                glEndList()        
            
        if self.selectOption == 'faces':
            pass
        if self.selectOption == 'nodes':
            glPointSize(10.0)
        else:
            pass



    