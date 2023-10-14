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
        
        self.materials = {'6061-T6_aluminum (m kg N)': 	 		{'Elasticity': 689e8, 'Poisson ratio': 0.35, 'Density':   2700.},
						  '1010_carbon_steel (m kg N)':	 		{'Elasticity': 205e9, 'Poisson ratio': 0.29, 'Density':   7870.},
						  '316_stainless_steel (m kg N)':   	{'Elasticity': 193e9, 'Poisson ratio': 0.27, 'Density':   7870.},
						  'Grade2_titanium (m kg N)':			{'Elasticity': 105e9, 'Poisson ratio': 0.37, 'Density':   4510.},
						  '6061-T6_aluminum (mm kg mN kPa)':	{'Elasticity': 689e5, 'Poisson ratio': 0.35, 'Density':  2.7e-6},
						  '1010_carbon_steel (mm kg mN kPa)':	{'Elasticity': 205e6, 'Poisson ratio': 0.29, 'Density': 7.87e-6},
						  '316_stainless_steel (mm kg mN kPa)':	{'Elasticity': 193e6, 'Poisson ratio': 0.27, 'Density': 7.87e-6},
						  'Grade2_titanium (mm kg mN kPa)':	 	{'Elasticity': 105e6, 'Poisson ratio': 0.37, 'Density': 4.51e-6}}
        self.sections = {}
        self.parts = { 'part-1': Part('part-1') }

        self.currentPart = self.parts['part-1']
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
        self.selectOption = 'Nodes'
		
        self.results = {}
        self.scaleShearBendDiagram = 1.
        self.scale_factor = 20.
        self.displayLists = {}


    def clearModel(self):
        '''
    Clears out all data in model from
    current session.
    '''
        pass
    
    
    
    