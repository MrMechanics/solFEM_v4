#
#
#   viewer.py
#  ---------------
#   This is the viewer module. It holds the viewer object,
#   which renders geometry, nodes, elements, constraints,
#   loads, boundary conditions, text, etc...
#   The viewer is interactive, letting the user select,
#   unselect, rotate, zoom and change views using keyboard
#   shortcuts.




from cameras import *
from geometries import *

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
except:
    print(' Error PyOpenGL not installed properly!!')
    sys.exit()

try:
    from PyQt5 import QtGui, QtWidgets, QtCore
    from PyQt5.QtOpenGL import *
except:
    print(' Error PyQt5 not installed properly!!')
    sys.exit()


import numpy as np





class Viewer(QGLWidget):
    '''
3D viewer for rendering the geometry
and mesh with or without contour plots,
animations, loads, etc.
'''
    def __init__(self, interface):
        QGLWidget.__init__(self, interface)
        
        self.interface = interface
        self.model = interface.model
        
        self.setMouseTracking(True)
        self.setMinimumSize(500, 500)
        self.width = 100
        self.height = 100

        self.viewAssembly = False
        self.viewGeometry = True
        self.viewMesh = False
        self.viewBoundaries = False
        self.viewConstraints = False
        self.viewLoads = False
        self.viewSolutions = False
        self.viewResults = False
        self.viewOrigin = False
        self.viewNodes = False
        self.viewShaded = True
        self.viewWireframe = False
        self.viewAveraged = False
        self.viewMeshTree = True
        self.viewAnimate = True
        self.viewFrame = 0
        self.veiwFrameRising = True
        self.viewAnimationSpeed = [0.05+(0.005**2)*math.sin((math.pi*i)/12) for i in range(7)]
        self.viewAnimationSpeed = self.viewAnimationSpeed[::-1]
        self.viewAnimationSpeed += self.viewAnimationSpeed[:0:-1]
        self.viewLoadingMessage = False

        self.activeCTRL = False
        self.activeSHIFT = False
        self.activeALT = False
        self.oldx = self.oldy = 0
        self.mouseButtonPressed = False
        self.activeSelection = False
        self.selectionRectangleStart = [0,0]
        self.selectionRectangleEnd = [0,0]

        self.colors = {'background_pre':    (0.60546875, 0.68359375, 0.5546875, 1.0),
                       'background_post':   (0.6796875, 0.73046875, 0.77734375, 1.0),
                       'boundaries':        (0.4140625, 0.48828125, 0.5546875),
                       'boundaries_rot':    (0.336, 0.447, 0.588, 1.0),
                       'loads':             (0.6875, 0.3984375, 0.375, 1.0),
                       'constraints':       (0.7890625, 0.55859375, 0.2578125, 1.0)}
        self.currentDisplayList = { 'part':          'None',
                                    'solution':      'None',
                                    'result':        'None',
                                    'subresult':     'None',
                                    'info':          'None',
                                    'avg_info':      'None',
                                    'max_val':         None,
                                    'min_val':         None,
                                    'avg_max_val':     None,
                                    'avg_min_val':     None,
                                    'view_radius':       2,
                                    'view_scope':   {'max': [ 1., 1., 1.],
                                                     'min': [-1.,-1.,-1.] },
                                    'displayLists': {'orientation':  None,
                                                     'nodes':        None,
                                                     'wireframe':    None,
                                                     'shaded':       None,
                                                     'average':      None }}

        self.camera = Camera()
        self.camera.setSceneRadius( self.currentDisplayList['view_radius'] )
        self.camera.reset()
        self.modelCentered = False

        self.coordSys0 = CoordSys3D(Point3D(0.,0.,0.),Vector3D(1.,0.,0.),Vector3D(0.,1.,0.))
        self.coordSys0_centered = self.coordSys0




    def paintGL(self):
        '''
    Draw all the graphics in the viewer.
    First the triad to show the current view angle,
    then the current pre-built 3D displaylists,
    and finally the 2D informational text overlay.
    '''
        # Render all 3D objects in viewer
        # -------------------------------
        glMatrixMode( GL_PROJECTION )
        glLoadIdentity()
        self.camera.transform()
        glMatrixMode( GL_MODELVIEW )
        glLoadIdentity()

        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )

        glDepthFunc( GL_LEQUAL )
        glEnable( GL_DEPTH_TEST )
        glEnable( GL_CULL_FACE )
        glFrontFace( GL_CCW )
        glDisable( GL_LIGHTING )
        glShadeModel( GL_SMOOTH )
        
        # for transparent shear and bending moment diagrams
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        self.viewport = glGetIntegerv( GL_VIEWPORT )
        self.projection = glGetDoublev( GL_PROJECTION_MATRIX )
        self.view_matrix = glGetDoublev( GL_MODELVIEW_MATRIX )
        
        if self.modelCentered:
            glTranslatef(-(self.currentDisplayList['view_scope']['max'][0]+self.currentDisplayList['view_scope']['min'][0])/2.,
                         -(self.currentDisplayList['view_scope']['max'][1]+self.currentDisplayList['view_scope']['min'][1])/2.,
                         -(self.currentDisplayList['view_scope']['max'][2]+self.currentDisplayList['view_scope']['min'][2])/2.)

        if self.viewOrigin:
            glLineWidth(5.0)
            glColor(1.0, 0.0, 0.0)
            glBegin(GL_LINES)
            glVertex(0.,0.,0.)
            glVertex(0.15*self.currentDisplayList['view_radius'], 0., 0.)
            glEnd()
            glColor(0.0, 1.0, 0.0)
            glBegin(GL_LINES)
            glVertex(0.,0.,0.)
            glVertex(0.,0.15*self.currentDisplayList['view_radius'], 0.)
            glEnd()
            glColor(0.0, 0.0, 1.0)
            glBegin(GL_LINES)
            glVertex(0.,0.,0.)
            glVertex(0., 0., 0.15*self.currentDisplayList['view_radius'])
            glEnd()



        # Render all of the 3D model
        # --------------------------
        cDL = self.currentDisplayList

        if self.viewAssembly:
            pass
        
        if self.viewGeometry:
            pass

        if self.viewMesh:
            pass
        
        if self.viewBoundaries:
            pass
        
        if self.viewConstraints:
            pass
        
        if self.viewLoads:
            pass
        
        if self.viewSolutions:
            pass
        
        if self.viewResults:
            pass
        
        if self.viewNodes:
            if cDL['displayLists']['nodes'] != None:
                glCallList(cDL['displayLists']['nodes'])

        if self.viewWireframe:
            if cDL['displayLists']['wireframe'] != None:
                glCallList(cDL['displayLists']['wireframe'])
        else:
            if cDL['displayLists']['wireframe'] != None:
                glCallList(cDL['displayLists']['wireframe'])
            if self.viewAveraged:
                if cDL['displayLists']['average'] != None:
                    glCallList(cDL['displayLists']['average'])
            else:
                if cDL['displayLists']['shaded'] != None:
                    glCallList(cDL['displayLists']['shaded'])



        # Render what is selected
        # --------------------------
        if self.model.nodesSelected:
            pass
        elif self.model.elementsSelected:
            pass

        elif self.model.linesSelected:
            if self.model.displayLists['selected_lines'] != None:
                glCallList(self.model.displayLists['selected_lines'])
            else:
                self.model.selected_lines.clear()
                self.model.linesSelected = False
        
        elif self.model.facesSelected:
            if self.model.displayLists['selected_lines'] != None:
                glCallList(self.model.displayLists['selected_lines'])
            if self.model.displayLists['selected_faces'] != None:
                glCallList(self.model.displayLists['selected_faces'])
            else:
                self.model.selected_lines.clear()
                self.model.linesSelected = False
                
        else:
            pass



        # Draw RGB triad in left corner
        # -----------------------------
        glViewport(0, 0, self.width//3, self.height//3)
        if self.modelCentered:
            glTranslatef((self.currentDisplayList['view_scope']['max'][0]+self.currentDisplayList['view_scope']['min'][0])/2.,
                         (self.currentDisplayList['view_scope']['max'][1]+self.currentDisplayList['view_scope']['min'][1])/2.,
                         (self.currentDisplayList['view_scope']['max'][2]+self.currentDisplayList['view_scope']['min'][2])/2.)

        view_length = -(self.camera.position - self.camera.target).length()

        glLineWidth(5.0)
        glColor(1.0, 0.0, 0.0)
        # X-dir
        glBegin(GL_LINES)
        glVertex( self.camera.target.x(), self.camera.target.y(), self.camera.target.z())
        glVertex( self.camera.target.x() - 1.*view_length*0.2, self.camera.target.y(), self.camera.target.z())
        glEnd()
        glColor(0.0, 1.0, 0.0)
        # Y-dir
        glBegin(GL_LINES)
        glVertex( self.camera.target.x(), self.camera.target.y(), self.camera.target.z())
        glVertex( self.camera.target.x(), self.camera.target.y() - 1.*view_length*0.2, self.camera.target.z())
        glEnd()
        glColor(0.0, 0.0, 1.0)
        # Z-dir
        glBegin(GL_LINES)
        glVertex( self.camera.target.x(), self.camera.target.y(), self.camera.target.z())
        glVertex( self.camera.target.x(), self.camera.target.y(), self.camera.target.z() - 1.*view_length*0.2)
        glEnd()

        glLineWidth(2.0)
        glColor(1.0, 0.0, 0.0)
        # 3D-X
        glBegin(GL_LINES)
        glVertex( self.camera.target.x() - view_length*0.25, self.camera.target.y() - view_length*0.02, self.camera.target.z() + view_length*0.01)
        glVertex( self.camera.target.x() - view_length*0.21, self.camera.target.y() + view_length*0.02, self.camera.target.z() - view_length*0.01)
        glEnd()
        glBegin(GL_LINES)
        glVertex( self.camera.target.x() - view_length*0.21, self.camera.target.y() - view_length*0.02, self.camera.target.z() - view_length*0.01)
        glVertex( self.camera.target.x() - view_length*0.25, self.camera.target.y() + view_length*0.02, self.camera.target.z() + view_length*0.01)
        glEnd()
        glColor(0.0, 1.0, 0.0)
        # 3D-Y
        glBegin(GL_LINES)
        glVertex( self.camera.target.x(), self.camera.target.y() - view_length*0.23, self.camera.target.z())
        glVertex( self.camera.target.x() + view_length*0.02, self.camera.target.y() - view_length*0.25, self.camera.target.z() - view_length*0.01)
        glEnd()
        glBegin(GL_LINES)
        glVertex( self.camera.target.x(), self.camera.target.y() - view_length*0.23, self.camera.target.z())
        glVertex( self.camera.target.x() - view_length*0.02, self.camera.target.y() - view_length*0.25, self.camera.target.z() + view_length*0.01)
        glEnd()
        glBegin(GL_LINES)
        glVertex( self.camera.target.x(), self.camera.target.y() - view_length*0.23, self.camera.target.z())
        glVertex( self.camera.target.x(), self.camera.target.y() - view_length*0.21, self.camera.target.z())
        glEnd()
        glColor(0.0, 0.0, 1.0)
        # 3D-Z
        glBegin(GL_LINES)
        glVertex( self.camera.target.x() - view_length*0.02, self.camera.target.y() - view_length*0.02, self.camera.target.z() - view_length*0.235)
        glVertex( self.camera.target.x() + view_length*0.02, self.camera.target.y() + view_length*0.02, self.camera.target.z() - view_length*0.245)
        glEnd()
        glBegin(GL_LINES)
        glVertex( self.camera.target.x() - view_length*0.02, self.camera.target.y() - view_length*0.02, self.camera.target.z() - view_length*0.235)
        glVertex( self.camera.target.x() + view_length*0.02, self.camera.target.y() - view_length*0.02, self.camera.target.z() - view_length*0.245)
        glEnd()
        glBegin(GL_LINES)
        glVertex( self.camera.target.x() + view_length*0.02, self.camera.target.y() + view_length*0.02, self.camera.target.z() - view_length*0.245)
        glVertex( self.camera.target.x() - view_length*0.02, self.camera.target.y() + view_length*0.02, self.camera.target.z() - view_length*0.235)
        glEnd()

        glViewport(0, 0, self.width, self.height)
        
        
        
        # Render all 2D overlay
        # ---------------------
        glClear(GL_DEPTH_BUFFER_BIT)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluOrtho2D(0, self.width, self.height, 0)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()


        # change views using left side buttons
        if self.viewAssembly:
            pass
        if self.viewGeometry:
            glColor3f(1., 1., 1.)
            if self.model.currentPart != None:
                self.renderText(40, self.height-40, self.model.currentPart.name, QtGui.QFont( 'helvetica', 18 ) )
#            if self.viewMeshTree:
#                self.drawGeometryTree()
        elif self.viewMesh:
            pass
#            if self.viewMeshTree:
#                self.drawMeshTree()
        elif self.viewBoundaries:
            pass
        elif self.viewConstraints:
            pass
        elif self.viewLoads:
            pass            
        elif self.viewSolutions:
            pass
#            if self.viewMeshTree:
#                self.drawMeshTree()
#            if self.model.currentSolution != None:
#                self.writeInfo(self.currentSolution)
        elif self.viewResults:
            pass
        else:
            pass


        # draw rectangle of selection
        if self.mouseButtonPressed == True:
            if self.activeSelection == True:
                self.drawRectangle()


        # animate eigenmodes
        if self.viewAnimate == True and self.model.currentResults['result'] == 'Eigenmodes':
            if self.viewFrame == 12:
                self.veiwFrameRising = False
                self.viewFrame = 11
            elif self.viewFrame == 0:
                self.veiwFrameRising = True
                self.viewFrame = 1
            else:
                if self.veiwFrameRising == True:
                    self.viewFrame += 1
                else:
                    self.viewFrame -= 1
            time.sleep(self.viewAnimationSpeed[self.viewFrame])
            self.update()


        # display loading message
        if self.viewLoadingMessage == True:
            glColor3f(1., 1., 1.)
            self.renderText(int(self.width/2.)-50, int(self.height/2.), 'Loading...', QtGui.QFont( 'helvetica', 16 ) )


        glFlush()
        
        


        
    def resizeGL(self, widthInPixels, heightInPixels):
        '''
    Resize the viewer.
    '''
        self.camera.setViewportDimensions(widthInPixels, heightInPixels)
        self.width = widthInPixels
        self.height = heightInPixels
        glViewport(0, 0, widthInPixels, heightInPixels)


    def initializeGL(self):
        '''
    Initialize the viewer.
    '''
        glClearColor(self.colors['background_pre'][0], 
                     self.colors['background_pre'][1], 
                     self.colors['background_pre'][2], 
                     self.colors['background_pre'][3])
        glClearDepth(1.0)




    def mouseMoveEvent(self, mouseEvent):
        '''
    Update the camera view when holding down the
    different buttons while dragging the mouse
    around on the screen.
    '''
        if int(mouseEvent.buttons()) != QtCore.Qt.NoButton:
            # user is dragging
            delta_x = mouseEvent.x() - self.oldx
            delta_y = self.oldy - mouseEvent.y()
            if int(mouseEvent.buttons()) & QtCore.Qt.LeftButton:
                if (self.activeCTRL and self.activeALT):
                    self.camera.orbit(self.oldx,self.oldy,mouseEvent.x(),mouseEvent.y())
                else:
                    self.activeSelection = True
                    self.selectionRectangleEnd = [mouseEvent.x(), mouseEvent.y()]
            elif int(mouseEvent.buttons()) & QtCore.Qt.RightButton :
                self.camera.dollyCameraForward( 3*(delta_x+delta_y), False )
            elif int(mouseEvent.buttons()) & QtCore.Qt.MidButton :
                self.camera.translateSceneRightAndUp( delta_x, delta_y )
            self.update()
        self.oldx = mouseEvent.x()
        self.oldy = mouseEvent.y()


    def mouseDoubleClickEvent(self, mouseEvent):
        '''
    Clear the model selection variables when the
    left mousebutton is doubble-clicked.
    '''
        self.model.selected_nodes.clear()
        self.model.nodesSelected = False
        self.model.displayLists['selected_nodes'] = None
        self.model.selected_elements.clear()
        self.model.elementsSelected = False
        self.model.displayLists['selected_elements'] = None
        self.model.selected_lines.clear()
        self.model.linesSelected = False
        self.model.displayLists['selected_lines'] = None
        self.model.selected_faces.clear()
        self.model.facesSelected = False
        self.model.displayLists['selected_faces'] = None
        self.update()


    def mousePressEvent(self, e):
        '''
    Update the mouseButtonPressed and
    selectionRectangleStart variables to say
    that the left mouse button has been pressed.
    '''
        if self.mouseButtonPressed == False:
            self.selectionRectangleStart = [e.x(), e.y()]
        self.mouseButtonPressed = True


    def mouseReleaseEvent(self, e):
        '''
    Update the mouseButtonPressed and
    activeSelection variables to say that
    the left mouse button has been released.
    '''
        self.mouseButtonPressed = False
        if self.activeSelection:
            if self.model.selectOption in ['nodes', 'elements']:
                self.model.selected_nodes = self.select()
                if len(self.model.selected_lines) == 0:
                    self.model.displayLists['selected_nodes'] = None
                if len(self.model.selected_faces) == 0:
                    self.model.displayLists['selected_elements'] = None
            elif self.model.selectOption in ['lines', 'faces']:
                self.model.selected_lines = self.select()
                if len(self.model.selected_lines) == 0:
                    self.model.displayLists['selected_lines'] = None
                if len(self.model.selected_faces) == 0:
                    self.model.displayLists['selected_faces'] = None
            else:
                pass
            self.model.selectedFeaturesDisplayList()
            self.activeSelection = False
            self.update()
        
        
        

    def select(self):
        '''
    Select nodes, elements, lines or faces by
    clicking the left mouse button and dragging
    a selection window encompassing what is to
    be selected.
    '''
        if self.selectionRectangleEnd[0] < self.selectionRectangleStart[0]:
            x = [self.selectionRectangleEnd[0], self.selectionRectangleStart[0],
                 self.selectionRectangleEnd[0], self.selectionRectangleStart[0]]
        else:
            x = [self.selectionRectangleStart[0], self.selectionRectangleEnd[0],
                 self.selectionRectangleStart[0], self.selectionRectangleEnd[0]]
        if self.selectionRectangleEnd[1] < self.selectionRectangleStart[1]:
            y = [self.height-self.selectionRectangleEnd[1], self.height-self.selectionRectangleEnd[1],
                 self.height-self.selectionRectangleStart[1], self.height-self.selectionRectangleStart[1]]
        else:
            y = [self.height-self.selectionRectangleStart[1], self.height-self.selectionRectangleStart[1],
                 self.height-self.selectionRectangleEnd[1], self.height-self.selectionRectangleEnd[1]]
        startpoint = [(0.,0.,0.),(0.,0.,0.),(0.,0.,0.),(0.,0.,0.)]
        endpoint = [(0.,0.,0.),(0.,0.,0.),(0.,0.,0.),(0.,0.,0.)]
        for i in range(4):
            try:
                startpoint[i] = gluUnProject(x[i], y[i], 0., self.view_matrix, self.projection, self.viewport)
                endpoint[i] = gluUnProject(x[i], y[i], 1., self.view_matrix, self.projection, self.viewport)
            except ValueError:
                pass
            else:
                ray_direction = (endpoint[i][0] - startpoint[i][0],
                                 endpoint[i][1] - startpoint[i][1],
                                 endpoint[i][2] - startpoint[i][2])
                endpoint[i] = (startpoint[i][0] + ray_direction[0],
                               startpoint[i][1] + ray_direction[1],
                               startpoint[i][2] + ray_direction[2])
        P = np.zeros((8,3))
        P[0] = np.array([startpoint[0][0], startpoint[0][1], startpoint[0][2]]) # front top left
        P[1] = np.array([startpoint[1][0], startpoint[1][1], startpoint[1][2]]) # front top right
        P[2] = np.array([startpoint[2][0], startpoint[2][1], startpoint[2][2]]) # front bottom left
        P[3] = np.array([startpoint[3][0], startpoint[3][1], startpoint[3][2]]) # front bottom right
        P[4] = np.array([  endpoint[0][0],   endpoint[0][1],   endpoint[0][2]]) # back top left
        P[5] = np.array([  endpoint[1][0],   endpoint[1][1],   endpoint[1][2]]) # back top right
        P[6] = np.array([  endpoint[2][0],   endpoint[2][1],   endpoint[2][2]]) # back bottom left
        P[7] = np.array([  endpoint[3][0],   endpoint[3][1],   endpoint[3][2]]) # back bottom right
        Frustum = []
        Frustum.append(np.cross((P[5]-P[1]),(P[4]-P[0]))) # top plane normal vector
        Frustum.append(np.cross((P[7]-P[3]),(P[5]-P[1]))) # right plane normal vector
        Frustum.append(np.cross((P[6]-P[2]),(P[7]-P[3]))) # bottom plane normal vector
        Frustum.append(np.cross((P[4]-P[0]),(P[6]-P[2]))) # left plane normal vector

        if self.model.selectOption in ['nodes', 'elements']:
            if len(self.model.currentPart.mesh.nodes) != 0:
                selected_nodes = {}
                meshnodes = self.currentPart.mesh.nodes
            for node in meshnodes:
                if self.modelCentered:
                    point_to_check = np.array([meshnodes[node].coord[0] + self.coordSys0_centered.origin.x(),
                                               meshnodes[node].coord[1] + self.coordSys0_centered.origin.y(),
                                               meshnodes[node].coord[2] + self.coordSys0_centered.origin.z()])
                else:
                    point_to_check = np.array([meshnodes[node].coord[0],
                                               meshnodes[node].coord[1],
                                               meshnodes[node].coord[2]])
                if np.dot(P[0]-point_to_check,Frustum[0]) < 0:
                    pass
                elif np.dot(P[1]-point_to_check,Frustum[1]) < 0:
                    pass
                elif np.dot(P[3]-point_to_check,Frustum[2]) < 0:
                    pass
                elif np.dot(P[2]-point_to_check,Frustum[3]) < 0:
                    pass
                else:
                    selected_nodes[node] = meshnodes[node]
            self.model.selectedFeaturesDisplayList()
            return selected_nodes

        elif self.model.selectOption in ['lines', 'faces']: 
            selected_lines = {}
            if len(self.model.currentPart.lines) != 0:
                partlines = self.model.currentPart.lines
                for line in partlines:
                    all_line_points_selected = True
                    for point in range(len(partlines[line].points)):
                        if self.modelCentered:
                            point_to_check = np.array([partlines[line].points[point].x() + self.coordSys0_centered.origin.x(),
                                                       partlines[line].points[point].y() + self.coordSys0_centered.origin.y(),
                                                       partlines[line].points[point].z() + self.coordSys0_centered.origin.z()])
                        else:
                            point_to_check = np.array([partlines[line].points[point].x(), 
                                                       partlines[line].points[point].y(), 
                                                       partlines[line].points[point].z()])
                        if np.dot(P[0]-point_to_check,Frustum[0]) < 0:
                            all_line_points_selected = False
                            break
                        elif np.dot(P[1]-point_to_check,Frustum[1]) < 0:
                            all_line_points_selected = False
                            break
                        elif np.dot(P[3]-point_to_check,Frustum[2]) < 0:
                            all_line_points_selected = False
                            break
                        elif np.dot(P[2]-point_to_check,Frustum[3]) < 0:
                            all_line_points_selected = False
                            break
                        else:
                            pass
                    if all_line_points_selected:
                        selected_lines[line] = partlines[line]
            return selected_lines




    def drawRectangle(self):
        '''
    Draw a rectangle from where in the viewer the 
    mouse pointer is clicked to where it is released.
    '''
        glLineWidth(2.0)
        glColor3f(0., 0., 0.)

        glEnable(GL_LINE_STIPPLE)
        factor = 3
        pattern = 0x5555
        glLineStipple(factor, pattern)
 
        glBegin(GL_LINE_LOOP)
        glVertex2f(self.selectionRectangleStart[0], self.selectionRectangleStart[1])
        glVertex2f(self.selectionRectangleEnd[0],   self.selectionRectangleStart[1])
        glVertex2f(self.selectionRectangleEnd[0],   self.selectionRectangleEnd[1])
        glVertex2f(self.selectionRectangleStart[0], self.selectionRectangleEnd[1])
        glEnd()

        glDisable(GL_LINE_STIPPLE)
        
        
        
        
    def updateDisplayList(self):
        '''
    Update the currentDisplayList variable to show the 
    current geometry or mesh with the selected loads and 
    boundary conditions, or the selected results.
    '''
        cDL = self.currentDisplayList
        if self.viewGeometry:
            if self.model.currentPart != None:
                cDL['part'] = self.model.currentPart.name
                cDL['view_radius'] = self.model.currentPart.view_radius
                cDL['view_scope'] = self.model.currentPart.view_scope
                cDL['displayLists']['orientation'] = None
                cDL['displayLists']['nodes'] = self.model.currentPart.displayLists['seeds']
                cDL['displayLists']['wireframe'] = self.model.currentPart.displayLists['lines']
                cDL['displayLists']['shaded'] = self.model.currentPart.displayLists['faces']
                cDL['displayLists']['average'] = None

        elif self.viewResults:
            pass
                
        elif self.viewMesh:
            pass
            
        else:
            pass

        self.camera.setSceneRadius( cDL['view_radius'] )
        
        
        