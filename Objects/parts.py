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

import numpy as np

from geometries import *
from meshes import *
from reader import *
from OpenGL.GL import *
from OpenGL.GLU import *




class Part(object):
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

        self.view_scope = {'max': [ 1., 1., 1.],
                           'min': [-1.,-1.,-1.]}
        self.view_radius = 2.
        self.colors = {'part_lines':         (0.05, 0.10, 0.05, 1.0),
                       'part_faces':         (0.32, 0.43, 0.49, 1.0),
                       'element_lines_pre':  (0.50, 0.50, 0.50, 1.0),
                       'element_lines_post': (0.50, 0.50, 0.50, 1.0),
                       'element_faces_pre':  (0.50, 0.50, 0.50, 1.0),
                       'element_faces_post': (0.50, 0.50, 0.50, 1.0)}
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
        self.step_geometry = StepFileData(filename)
        geom = self.step_geometry
        for f in geom.advanced_face:

            # type of face
            if geom.advanced_face[f][-2] in geom.plane:
                face_type = 'plane'
            elif geom.advanced_face[f][-2] in geom.cylindrical_surface:
                face_type = 'cylindrical'
            elif geom.advanced_face[f][-2] in geom.toroidal_surface:
                face_type = 'toroidal'
            elif geom.advanced_face[f][-2] in geom.b_spline_surface:
                face_type = 'b_spline'
            else:
                face_type = 'UNKNOWN'
                print('\n  UNKNOWN face type for face number:', f)
            self.faces[f] = Face(f,face_type)

            # face edges
            for e in geom.advanced_face[f][:-2]:
                if e in geom.face_outer_bound:
                    bound = geom.face_outer_bound
                else: 
                    bound = geom.face_bound
                if e not in self.edges:
                    self.edges[e] = Edge(e)
                self.faces[f].edges[e] = self.edges[e]

                # edge lines
                for edge_loop in geom.edge_loop[bound[e][0]]:
                    l = geom.oriented_edge[edge_loop][2]
                    if l not in self.lines:
                        # line type
                        if geom.edge_curve[l][2] in geom.line:
                            self.lines[l] = Line(l)
                            p1 = geom.cartesian_point[geom.vertex_point[geom.edge_curve[l][0]]]
                            p2 = geom.cartesian_point[geom.vertex_point[geom.edge_curve[l][1]]]
                            self.lines[l].newPoints([p1,p2])
                        elif geom.edge_curve[l][2] in geom.circle:
                            self.lines[l] = Arc(l)
                            p1 = geom.cartesian_point[geom.vertex_point[geom.edge_curve[l][0]]]
                            p2 = geom.cartesian_point[geom.vertex_point[geom.edge_curve[l][1]]]
                            r = geom.circle[geom.edge_curve[l][2]][1]
                            c = geom.cartesian_point[geom.axis2_placement_3D[geom.circle[geom.edge_curve[l][2]][0]][0]]
                            z_vec = geom.direction[geom.axis2_placement_3D[geom.circle[geom.edge_curve[l][2]][0]][1]]
                            x_vec = geom.direction[geom.axis2_placement_3D[geom.circle[geom.edge_curve[l][2]][0]][2]]
                            y_vec = np.cross(z_vec,x_vec)
                            self.lines[l].setRadius(r)
                            self.lines[l].setCenter(c)
                            self.lines[l].setAxis(x_vec,y_vec)
                            self.lines[l].newPoints([p1,p2])
                        elif geom.edge_curve[l][2] in geom.ellipse:
                            self.lines[l] = Ellipse(l)
                            p1 = geom.cartesian_point[geom.vertex_point[geom.edge_curve[l][0]]]
                            p2 = geom.cartesian_point[geom.vertex_point[geom.edge_curve[l][1]]]
                            r1 = geom.ellipse[geom.edge_curve[l][2]][1]
                            r2 = geom.ellipse[geom.edge_curve[l][2]][2]
                            c = geom.cartesian_point[geom.axis2_placement_3D[geom.ellipse[geom.edge_curve[l][2]][0]][0]]
                            z_vec = geom.direction[geom.axis2_placement_3D[geom.ellipse[geom.edge_curve[l][2]][0]][1]]
                            x_vec = geom.direction[geom.axis2_placement_3D[geom.ellipse[geom.edge_curve[l][2]][0]][2]]
                            y_vec = np.cross(z_vec,x_vec)
                            self.lines[l].setRadius(r1,r2)
                            self.lines[l].setCenter(c)
                            self.lines[l].setAxis(x_vec,y_vec)
                            self.lines[l].newPoints([p1,p2])
                        elif geom.edge_curve[l][2] in geom.b_spline_curve:
                            self.lines[l] = Spline(l)
                        else:
                            print('\nUNKNOWN line type for line number:', l)

                    self.edges[e].lines[l] = self.lines[l]

        print('len(self.lines):', len(self.lines))

        self.view_radius = max((geom.x_max - geom.x_min, geom.y_max - geom.y_min, geom.z_max-geom.z_min))*0.5
        self.view_scope = {'max': [geom.x_max, geom.y_max, geom.z_max],
                           'min': [geom.x_min, geom.y_min, geom.z_min]}
        self.generateDisplayLists('geometry')



    
    def importMesh(self,filename):
        '''
    Import mesh from *.sol file or from file
    generated in some other program. Supported
    file types are *.bdf (Nastran), 
    *.inp (ABAQUS) and *.dat (FreeCAD).
    '''
        pass
    
    
    
    
    def generateDisplayLists(self,displ_type='geometry'):
        '''
    Generates an updated displaylist for the part to
    be rendered in the viewer.
    '''
        if displ_type == 'geometry':
            self.displayLists['lines'] = glGenLists(1)
#            self.displayLists['faces'] = glGenLists(1)
#            self.displayLists['seeds'] = glGenLists(1)

            glNewList(self.displayLists['lines'], GL_COMPILE)
            glLineWidth(3.0)
            glColor3f(self.colors['part_lines'][0], 
                      self.colors['part_lines'][1],
                      self.colors['part_lines'][2])
            for l in self.lines:
                if self.lines[l].type == 'line':
                    glBegin(GL_LINES)
                    glVertex3f(self.lines[l].points[0].x(),self.lines[l].points[0].y(),self.lines[l].points[0].z())
                    glVertex3f(self.lines[l].points[1].x(),self.lines[l].points[1].y(),self.lines[l].points[1].z())
                    glEnd()
                elif self.lines[l].type == 'arc':
                    for p in range(len(self.lines[l].points)-1):
                        glBegin(GL_LINES)
                        glVertex3f(self.lines[l].points[p].x(),self.lines[l].points[p].y(),self.lines[l].points[p].z())
                        glVertex3f(self.lines[l].points[p+1].x(),self.lines[l].points[p+1].y(),self.lines[l].points[p+1].z())
                        glEnd()
                else:
                    pass

            glEndList()

            
            
        elif displ_type == 'mesh':
            pass
        else:
            pass