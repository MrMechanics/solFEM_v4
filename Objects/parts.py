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
                       'part_faces':         (0.66796875, 0.6171875, 0.44921875, 1.0),
                       'part_faces_inside':  (0.62, 0.63, 0.19, 1.0),
                       'element_lines_pre':  (0.50, 0.50, 0.50, 1.0),
                       'element_lines_post': (0.50, 0.50, 0.50, 1.0),
                       'element_faces_pre':  (0.336, 0.447, 0.588, 1.0),
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
            elif geom.advanced_face[f][-2] in geom.conical_surface:
                face_type = 'conical'
            elif geom.advanced_face[f][-2] in geom.surface_of_revolution:
                face_type = 'revolution'
            elif geom.advanced_face[f][-2] in geom.surface_of_linear_extrusion:
                face_type = 'extrusion'
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
                            self.lines[l].setRadiuses(r2,r1)
                            self.lines[l].setCenter(c)
                            self.lines[l].setAxis(x_vec,y_vec)
                            self.lines[l].newPoints([p1,p2])
                        elif geom.edge_curve[l][2] in geom.b_spline_curve:
                            self.lines[l] = Spline(l)
                            self.lines[l].newPoints(geom,geom.edge_curve[l][2])
                        else:
                            print('\nUNKNOWN line type for line number:', l)

                    self.edges[e].lines[l] = self.lines[l]
                # get edge points
                self.edges[e].getEdgePoints()

            # generate surface mesh for rendering
            self.faces[f].surfaceNormal(geom)
            self.faces[f].getFacePoints()
            self.mesher.meshFace(self.faces[f])


#        print('number of faces:', len(self.faces))
#        print(self.faces.keys())
#        for f in self.faces:
#            print('\n\nface:', f)
#            print('face type:', self.faces[f].type)
#            print('face points:', self.faces[f].points)
#            print('face centroid:', self.faces[f].centroid)
#            if self.faces[f].type == 'plane':
#                print('face_normal:', self.faces[f].normal_v)
#            print('face_edges:', self.faces[f].edges.keys(), end='')
#            for e in self.faces[f].edges:
#                print('edge_points:\n', self.faces[f].edges[e].points)
#                print('centroid:', self.faces[f].edges[e].centroid)
#                print('\nedge', e, 'lines:')
#                for l in self.faces[f].edges[e].lines:
#                    print(l, end=', ')

        print('\n\t New Part: '+self.name)
        print('/----- ---------- ------  --------- ------\\')
        if len(geom.advanced_brep_shape_representation) == 1:
            print('  Solid model')
        elif len(geom.manifold_surface_shape_representation) == 1:
            print('  Shell model')
#        elif 2D planar model?
        else:
            print('Unknown model type')
        print('\n  Number of faces:    ', len(geom.advanced_face))
        print('  Number of vertices: ', len(geom.vertex_point))
		
        self.x_max = max(self.lines[l].points[j].x() for l in self.lines for j in range(len(self.lines[l].points)))
        self.x_min = min(self.lines[l].points[j].x() for l in self.lines for j in range(len(self.lines[l].points)))
        self.y_max = max(self.lines[l].points[j].y() for l in self.lines for j in range(len(self.lines[l].points)))
        self.y_min = min(self.lines[l].points[j].y() for l in self.lines for j in range(len(self.lines[l].points)))
        self.z_max = max(self.lines[l].points[j].z() for l in self.lines for j in range(len(self.lines[l].points)))
        self.z_min = min(self.lines[l].points[j].z() for l in self.lines for j in range(len(self.lines[l].points)))
		
        print(f'\n  x-min: {self.x_min:.2f} \tx-max: {self.x_max:.2f}')
        print(f'  y-min: {self.y_min:.2f} \ty-max: {self.y_max:.2f}')
        print(f'  z-min: {self.z_min:.2f} \tz-max: {self.z_max:.2f}')
#        self.center_of_mass = (0., 0., 0.)
#        print('\n  Center of mass:', self.center_of_mass)
#        self.volume = 0.
#        print('  Volume:', self.volume)
        print('\\----- ---------- ------  --------- ------/')

        self.view_radius = max((self.x_max - self.x_min, self.y_max - self.y_min, self.z_max-self.z_min))*0.5
        self.view_scope = {'max': [self.x_max, self.y_max, self.z_max],
                           'min': [self.x_min, self.y_min, self.z_min]}
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
            self.displayLists['faces'] = glGenLists(1)
            self.displayLists['seeds'] = glGenLists(1)



            # -----------
            # DRAW FACES
            # -------------------
            glNewList(self.displayLists['faces'], GL_COMPILE)
            for f in self.faces:
                if self.faces[f].type in ['conical', 'toroidal', 'plane', 'cylindrical']:
                    glColor3f(self.colors['part_faces'][0], 
                              self.colors['part_faces'][1],
                              self.colors['part_faces'][2])
                    for e in self.faces[f].g_mesh['elements']:
                        glBegin(GL_TRIANGLES)
                        elm = self.faces[f].g_mesh['elements'][e]
                        glVertex3f(self.faces[f].g_mesh['nodes'][elm[0]][0],
                                   self.faces[f].g_mesh['nodes'][elm[0]][1],
                                   self.faces[f].g_mesh['nodes'][elm[0]][2])
                        glVertex3f(self.faces[f].g_mesh['nodes'][elm[1]][0],
                                   self.faces[f].g_mesh['nodes'][elm[1]][1],
                                   self.faces[f].g_mesh['nodes'][elm[1]][2])
                        glVertex3f(self.faces[f].g_mesh['nodes'][elm[2]][0],
                                   self.faces[f].g_mesh['nodes'][elm[2]][1],
                                   self.faces[f].g_mesh['nodes'][elm[2]][2])
                        glEnd()
                        glBegin(GL_TRIANGLES)
                        glVertex3f(self.faces[f].g_mesh['nodes'][elm[2]][0],
                                   self.faces[f].g_mesh['nodes'][elm[2]][1],
                                   self.faces[f].g_mesh['nodes'][elm[2]][2])
                        glVertex3f(self.faces[f].g_mesh['nodes'][elm[1]][0],
                                   self.faces[f].g_mesh['nodes'][elm[1]][1],
                                   self.faces[f].g_mesh['nodes'][elm[1]][2])
                        glVertex3f(self.faces[f].g_mesh['nodes'][elm[0]][0],
                                   self.faces[f].g_mesh['nodes'][elm[0]][1],
                                   self.faces[f].g_mesh['nodes'][elm[0]][2])
                        glEnd()
                else:
                    pass
            glEndList()



            # -----------
            # DRAW LINES
            # -------------------
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
                elif self.lines[l].type == 'ellipse':
                    for p in range(len(self.lines[l].points)-1):
                        glBegin(GL_LINES)
                        glVertex3f(self.lines[l].points[p].x(),self.lines[l].points[p].y(),self.lines[l].points[p].z())
                        glVertex3f(self.lines[l].points[p+1].x(),self.lines[l].points[p+1].y(),self.lines[l].points[p+1].z())
                        glEnd()
                elif self.lines[l].type == 'spline':
                    for p in range(len(self.lines[l].points)-1):
                        glBegin(GL_LINES)
                        glVertex3f(self.lines[l].points[p].x(),self.lines[l].points[p].y(),self.lines[l].points[p].z())
                        glVertex3f(self.lines[l].points[p+1].x(),self.lines[l].points[p+1].y(),self.lines[l].points[p+1].z())
                        glEnd()
                else:
                    pass

            # draw surface normals for debugging
            for f in self.faces:
                if self.faces[f].type == 'planeiii':
                    glColor3f(0.9,0,0)
                    glBegin(GL_LINES)
                    glVertex3f(self.faces[f].normal[0].x(),self.faces[f].normal[0].y(),self.faces[f].normal[0].z())
                    glVertex3f(self.faces[f].normal[1].x(),self.faces[f].normal[1].y(),self.faces[f].normal[1].z())
                    glEnd()
                if self.faces[f].type == 'cylindricaliii':
                    glColor3f(0.9,0,0)
                    if self.faces[f].inwards == True:
                        glColor3f(0,0.9,0)
                    glBegin(GL_LINES)
                    glVertex3f(self.faces[f].normal[0].x(),self.faces[f].normal[0].y(),self.faces[f].normal[0].z())
                    glVertex3f(self.faces[f].normal[1].x(),self.faces[f].normal[1].y(),self.faces[f].normal[1].z())
                    glEnd()
                if self.faces[f].type == 'toroidaliii':
                    glColor3f(0.9,0,0)
                    if self.faces[f].inwards == True:
                        glColor3f(0,0.9,0)
                    glBegin(GL_LINES)
                    glVertex3f(self.faces[f].normal[0].x(),self.faces[f].normal[0].y(),self.faces[f].normal[0].z())
                    glVertex3f(self.faces[f].normal[1].x(),self.faces[f].normal[1].y(),self.faces[f].normal[1].z())
                    glEnd()
                if self.faces[f].type == 'conical':
                    glColor3f(0.9,0,0)
                    if self.faces[f].inwards == True:
                        glColor3f(0,0.9,0)
                    glBegin(GL_LINES)
                    glVertex3f(self.faces[f].normal[0].x(),self.faces[f].normal[0].y(),self.faces[f].normal[0].z())
                    glVertex3f(self.faces[f].normal[1].x(),self.faces[f].normal[1].y(),self.faces[f].normal[1].z())
                    glEnd()
                    
            glEndList()


            
            # -----------
            # DRAW SEEDS
            # -------------------
            glNewList(self.displayLists['seeds'], GL_COMPILE)

            glPointSize(8.0)
            glColor3f(0.4, 0.65, 0.4)
            glBegin(GL_POINTS)
            for e in self.edges:
                for p in self.edges[e].points:
                    glVertex3f(p[0],p[1],p[2])
            glEnd()

            glEndList()


        elif displ_type == 'mesh':
            pass
        else:
            pass