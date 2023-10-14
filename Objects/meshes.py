#
#
#	meshes.py
#  ---------
#	Mesh object holding all nodes and elements 
#	used in finite element model. Also used to store results.
#	





class Mesh(object):
    '''
Class for mesh. Holds information on node coordinates
element nodes, element type, and analysis results.
It is stored to binary file for easy access of results
to the viewer.
'''
    def __init__(self,nodes,elements):
        self.nodes = nodes
        self.elements = elements
        self.solutions = []

        self.externalNodes()


    def externalNodes(self):
        self.external = {}
        



