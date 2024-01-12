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
        self.colors = {'part_lines':         (0.37109375, 0.43359375, 0.3671875,  1.0),
                       'part_faces':         (0.66796875, 0.6171875, 0.44921875,  1.0),
                       'part_faces_inside':  (0.62, 0.63, 0.19, 1.0),
                       'element_lines_pre':  (0.37109375, 0.43359375, 0.3671875,  1.0),
                       'element_lines_post': (0.328125,   0.3984375,  0.41796875, 1.0),
                       'element_faces_pre':  (0.4609375,  0.703125,   0.46484375, 1.0),
                       'element_faces_post': (0.46484375, 0.640625,   0.6875,     1.0),
                       'nodes_pre':          (0.37109375, 0.43359375, 0.3671875,  1.0),
                       'nodes_post':         (0.328125,   0.3984375,  0.41796875, 1.0),
                       'seeds':              (0.6484375,  0.8046875,  0.65234375, 1.0)}
        self.displayLists = {'nodes':       None,
                             'wireframe':   None,
                             'elements':    None,
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
        for l in self.lines:
            self.mesher.seedLine(self.lines[l],5.)
        for e in self.edges:
            self.edges[e].updateSeeds()
        self.generateDisplayLists('geometry')



    
    def importMesh(self,filename):
        '''
    Import mesh from *.sol file or from file
    generated in some other program. Supported
    file types are *.bdf (Nastran), 
    *.inp (ABAQUS) and *.dat (FreeCAD).
    '''
        pass
    
    
    
    
    def generateDisplayLists(self,displ_type='geometry',draw_lines=True,draw_faces=True,draw_seeds=True):
        '''
    Generates an updated displaylist for the part to
    be rendered in the viewer.
    '''
#        print(f'self.displayLists (before): {self.displayLists}')
        if displ_type == 'geometry':
#            print('INSIDE GEOMETRY PART OF generateDisplayLists()')
            # -----------
            # DRAW FACES
            # -------------------
            if draw_faces:
                self.displayLists['faces'] = glGenLists(1)
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
            if draw_lines:
                self.displayLists['lines'] = glGenLists(1)
                glNewList(self.displayLists['lines'], GL_COMPILE)
                glLineWidth(3.0)
                glColor3f(self.colors['part_lines'][0], 
                          self.colors['part_lines'][1],
                          self.colors['part_lines'][2])
                glBegin(GL_LINES)
                for l in self.lines:
                    if self.lines[l].type == 'line':
                        glVertex3f(self.lines[l].points[0].x(),self.lines[l].points[0].y(),self.lines[l].points[0].z())
                        glVertex3f(self.lines[l].points[1].x(),self.lines[l].points[1].y(),self.lines[l].points[1].z())
                    elif self.lines[l].type == 'arc':
                        for p in range(len(self.lines[l].points)-1):
                            glVertex3f(self.lines[l].points[p].x(),self.lines[l].points[p].y(),self.lines[l].points[p].z())
                            glVertex3f(self.lines[l].points[p+1].x(),self.lines[l].points[p+1].y(),self.lines[l].points[p+1].z())
                    elif self.lines[l].type == 'ellipse':
                        for p in range(len(self.lines[l].points)-1):
                            glVertex3f(self.lines[l].points[p].x(),self.lines[l].points[p].y(),self.lines[l].points[p].z())
                            glVertex3f(self.lines[l].points[p+1].x(),self.lines[l].points[p+1].y(),self.lines[l].points[p+1].z())
                    elif self.lines[l].type == 'spline':
                        for p in range(len(self.lines[l].points)-1):
                            glVertex3f(self.lines[l].points[p].x(),self.lines[l].points[p].y(),self.lines[l].points[p].z())
                            glVertex3f(self.lines[l].points[p+1].x(),self.lines[l].points[p+1].y(),self.lines[l].points[p+1].z())
                    else:
                        pass
                glEnd()
    
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
                    if self.faces[f].type == 'conicaliii':
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
            if draw_seeds:
                self.displayLists['seeds'] = glGenLists(1)
                glNewList(self.displayLists['seeds'], GL_COMPILE)
    
                glPointSize(8.0)
                glColor3f(self.colors['seeds'][0], 
                          self.colors['seeds'][1],
                          self.colors['seeds'][2])
                glBegin(GL_POINTS)
                for e in self.edges:
                    for s in self.edges[e].seeds:
                        glVertex3f(s[0],s[1],s[2])
                glEnd()
    
                glEndList()


        elif displ_type == 'mesh':
#            print('INSIDE MESH PART OF generateDisplayLists()')
            # -----------
            # DRAW ELEMENT FACES
            # -------------------
            self.displayLists['elements'] = glGenLists(1)
            glNewList(self.displayLists['elements'], GL_COMPILE)
            for f in self.faces:
                glColor3f(self.colors['element_faces_pre'][0], 
                          self.colors['element_faces_pre'][1],
                          self.colors['element_faces_pre'][2])
                if len(self.faces[f].mesh['elements']) != 0:
                    for e in self.faces[f].mesh['elements']:
                        elm = self.faces[f].mesh['elements'][e]
                        glBegin(GL_TRIANGLES)
                        glVertex3f(self.faces[f].mesh['nodes'][elm[0]][0],
                                   self.faces[f].mesh['nodes'][elm[0]][1],
                                   self.faces[f].mesh['nodes'][elm[0]][2])
                        glVertex3f(self.faces[f].mesh['nodes'][elm[1]][0],
                                   self.faces[f].mesh['nodes'][elm[1]][1],
                                   self.faces[f].mesh['nodes'][elm[1]][2])
                        glVertex3f(self.faces[f].mesh['nodes'][elm[2]][0],
                                   self.faces[f].mesh['nodes'][elm[2]][1],
                                   self.faces[f].mesh['nodes'][elm[2]][2])
                        glEnd()
                        glBegin(GL_TRIANGLES)
                        glVertex3f(self.faces[f].mesh['nodes'][elm[2]][0],
                                   self.faces[f].mesh['nodes'][elm[2]][1],
                                   self.faces[f].mesh['nodes'][elm[2]][2])
                        glVertex3f(self.faces[f].mesh['nodes'][elm[1]][0],
                                   self.faces[f].mesh['nodes'][elm[1]][1],
                                   self.faces[f].mesh['nodes'][elm[1]][2])
                        glVertex3f(self.faces[f].mesh['nodes'][elm[0]][0],
                                   self.faces[f].mesh['nodes'][elm[0]][1],
                                   self.faces[f].mesh['nodes'][elm[0]][2])
                        glEnd()
            glEndList()

            # -----------
            # DRAW ELEMENT LINES
            # -------------------
            self.displayLists['wireframe'] = glGenLists(1)
            glNewList(self.displayLists['wireframe'], GL_COMPILE)
            glLineWidth(3.0)
            glColor3f(self.colors['element_lines_pre'][0], 
                      self.colors['element_lines_pre'][1],
                      self.colors['element_lines_pre'][2])
            for f in self.faces:
                if len(self.faces[f].mesh['elements']) != 0:
                    glBegin(GL_LINES)
                    for e in self.faces[f].mesh['elements']:
                        elm = self.faces[f].mesh['elements'][e]
                        glVertex3f(self.faces[f].mesh['nodes'][elm[0]][0],
                                   self.faces[f].mesh['nodes'][elm[0]][1],
                                   self.faces[f].mesh['nodes'][elm[0]][2])
                        glVertex3f(self.faces[f].mesh['nodes'][elm[1]][0],
                                   self.faces[f].mesh['nodes'][elm[1]][1],
                                   self.faces[f].mesh['nodes'][elm[1]][2])
            
                        glVertex3f(self.faces[f].mesh['nodes'][elm[1]][0],
                                   self.faces[f].mesh['nodes'][elm[1]][1],
                                   self.faces[f].mesh['nodes'][elm[1]][2])
                        glVertex3f(self.faces[f].mesh['nodes'][elm[2]][0],
                                   self.faces[f].mesh['nodes'][elm[2]][1],
                                   self.faces[f].mesh['nodes'][elm[2]][2])
            
                        glVertex3f(self.faces[f].mesh['nodes'][elm[2]][0],
                                   self.faces[f].mesh['nodes'][elm[2]][1],
                                   self.faces[f].mesh['nodes'][elm[2]][2])
                        glVertex3f(self.faces[f].mesh['nodes'][elm[0]][0],
                                   self.faces[f].mesh['nodes'][elm[0]][1],
                                   self.faces[f].mesh['nodes'][elm[0]][2])
                    glEnd()
            glEndList()
        else:
            pass
#        print(f'self.displayLists (after): {self.displayLists}')
