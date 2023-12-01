#
#
#   interfacer.py
#  -----------------
#   This is the interface module. It uses pyqt5
#   to create a graphical user interface for the
#   user to interact with the finite element model
#
#

import sys
import re

sys.path.insert(1, '../Objects')

from viewer import *
from modeler import *
from mesher import *
from helper import *
from selector import *
from converter import *

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
    
    


class Interface(QtWidgets.QMainWindow):
    '''
The main user interface (GUI). This holds all
toolbars and menus. The FE-viewer (Viewer object)
is an OpenGL widget running inside this framework.
'''
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)

        self.setWindowIcon(QtGui.QIcon('../Icons/icon_view_result.png'))
        self.setWindowTitle('Finite Element Viewer')
        self.statusBar().showMessage('  ready  ')

        self.model = FEmodel(self)
        self.viewer = Viewer(self)

        exit = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_exit.png'),'Exit', self)
        exit.setShortcut('Ctrl+Q')
        exit.setStatusTip('Exit application')
        exit.triggered.connect(self.close)
        
        newsession = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_new_file.png'),'New', self)
        newsession.setStatusTip('Clear out current session')
        newsession.triggered.connect(self.clearModel)

        openfile = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_open.png'),'Open', self)
        openfile.setStatusTip('Open *.out or *.mdl file')
        openfile.triggered.connect(self.openFile)

        importfrom = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_import.png'),'Import', self)
        importfrom.setStatusTip('Import mesh from file')
        importfrom.triggered.connect(self.importFrom)

        export = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_export.png'),'Export', self)
        export.setStatusTip('Export mesh to *.sol file')
        export.triggered.connect(self.exportMesh)

        savefile = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_save.png'),'Save', self)
        savefile.setStatusTip('Save current session to *.mdl file')
        savefile.triggered.connect(self.saveFile)

        quicksave = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_quick_save.png'),'Quick Save', self)
        quicksave.setShortcut('Ctrl+S')
        quicksave.setStatusTip('Save current session to *.mdl file')
        quicksave.triggered.connect(self.quickSaveFile)

        selectnodes = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_select_nodes.png'),'Select Nodes', self)
        selectnodes.setStatusTip('Select nodes')
        selectnodes.triggered.connect(self.selectNodes)

        selectelements = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_select_elements.png'),'Select Elements', self)
        selectelements.setStatusTip('Select elements')
        selectelements.triggered.connect(self.selectElements)

        selectlines = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_select_lines.png'),'Select Lines', self)
        selectlines.setStatusTip('Select lines')
        selectlines.triggered.connect(self.selectLines)

        selectfaces = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_select_faces.png'),'Select Faces', self)
        selectfaces.setStatusTip('Select faces')
        selectfaces.triggered.connect(self.selectFaces)

        seedlines = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_seed_lines.png'),'Seed Lines', self)
        seedlines.setStatusTip('Seed lines')
        seedlines.triggered.connect(self.seedLines)

        seedpart = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_seed_geometry.png'),'Seed Part', self)
        seedpart.setStatusTip('Seed geometry')
        seedpart.triggered.connect(self.seedPart)

        meshfaces = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_mesh_face.png'),'Mesh Faces', self)
        meshfaces.setStatusTip('Mesh faces')
        meshfaces.triggered.connect(self.meshFaces)

        meshpart = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_mesh_geometry.png'),'Mesh Part', self)
        meshpart.setStatusTip('Mesh part')
        meshpart.triggered.connect(self.meshPart)

        delete = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_delete.png'),'Delete', self)
        delete.setStatusTip('Delete sets, mesh, material, solutions, boundaries, loads...')
        delete.triggered.connect(self.deleteItem)

        resetview = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_reset_view.png'),'Reset view', self)
        resetview.setShortcut('R')
        resetview.setStatusTip('Reset view to origin')
        resetview.triggered.connect(self.viewer.camera.reset)
        resetview.triggered.connect(self.viewer.update)

        centerview = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_center_view.png'),'Center view', self)
        centerview.setShortcut('C')
        centerview.setStatusTip('Center/uncenter view on model')
        centerview.triggered.connect(self.centerModel)
        centerview.triggered.connect(self.viewer.camera.reset)
        centerview.triggered.connect(self.viewer.update)
        
        viewleft = QtWidgets.QAction(QtGui.QIcon(''), 'View Left', self)
        viewleft.setShortcut('1')
        viewleft.setStatusTip('View Left')
        viewleft.triggered.connect(self.viewLeft)
        viewleft.triggered.connect(self.viewer.update)

        viewfront = QtWidgets.QAction(QtGui.QIcon(''), 'View Front', self)
        viewfront.setShortcut('2')
        viewfront.setStatusTip('View Front')
        viewfront.triggered.connect(self.viewFront)
        viewfront.triggered.connect(self.viewer.update)

        viewright = QtWidgets.QAction(QtGui.QIcon(''), 'View Right', self)
        viewright.setShortcut('3')
        viewright.setStatusTip('View Right')
        viewright.triggered.connect(self.viewRight)
        viewright.triggered.connect(self.viewer.update)

        viewback = QtWidgets.QAction(QtGui.QIcon(''), 'View Back', self)
        viewback.setShortcut('4')
        viewback.setStatusTip('View Back')
        viewback.triggered.connect(self.viewBack)
        viewback.triggered.connect(self.viewer.update)

        viewover = QtWidgets.QAction(QtGui.QIcon(''), 'View Over', self)
        viewover.setShortcut('5')
        viewover.setStatusTip('View Over')
        viewover.triggered.connect(self.viewOver)
        viewover.triggered.connect(self.viewer.update)

        viewunder = QtWidgets.QAction(QtGui.QIcon(''), 'View Under', self)
        viewunder.setShortcut('6')
        viewunder.setStatusTip('View Under')
        viewunder.triggered.connect(self.viewUnder)
        viewunder.triggered.connect(self.viewer.update)

        viewangled = QtWidgets.QAction(QtGui.QIcon(''), 'View XYZ', self)
        viewangled.setShortcut('7')
        viewangled.setStatusTip('View XYZ')
        viewangled.triggered.connect(self.viewAngled)
        viewangled.triggered.connect(self.viewer.update)

        nodesview = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_nodes.png'),'Nodes', self)
        nodesview.setShortcut('N')
        nodesview.setStatusTip('Toggle view of nodes on/off')
        nodesview.triggered.connect(self.nodesView)

        hide = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_hide_elements.png'),'Hide Elements', self)
        hide.setStatusTip('Hide selected elements')
        hide.triggered.connect(self.hideElements)

        show = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_show_elements.png'),'Show Elements', self)
        show.setStatusTip('Show all elements in mesh')
        show.triggered.connect(self.showElements)

        wireframe = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_wireframe.png'),'Wireframe', self)
        wireframe.setShortcut('W')
        wireframe.setStatusTip('Change view to wireframe mode')
        wireframe.triggered.connect(self.wireframeView)

        shaded = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_shaded.png'),'Shaded', self)
        shaded.setShortcut('S')
        shaded.setStatusTip('Change view to shaded mode')
        shaded.triggered.connect(self.shadedView)

        origin = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_origin_triad.png'),'Origin', self)
        origin.setShortcut('O')
        origin.setStatusTip('Show origin coordinate system')
        origin.triggered.connect(self.showOrigin)
        
        meshtree = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),'Mesh Tree', self)
        meshtree.setShortcut('H')
        meshtree.setStatusTip('Toggle Mesh Tree On/Off')
        meshtree.triggered.connect(self.showMeshTree)

        assemblyview = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_current_assembly.png'),'View Assembly...', self)
        assemblyview.setStatusTip('View assembly')
        assemblyview.triggered.connect(self.viewAssembly)

        partview = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_current_geometry.png'),'Current Part...', self)
        partview.setStatusTip('Select what part to view')
        partview.triggered.connect(self.selectPart)

        solutionview = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_current_solution.png'),'Current Solution...', self)
        solutionview.setStatusTip('Select what solution to view')
        solutionview.triggered.connect(self.selectSolution)

        resultview = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_current_results.png'),'Current Result...', self)
        resultview.setStatusTip('Select what result to view')
        resultview.triggered.connect(self.selectResult)

        highlightnode = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_highlight_node.png'),'Node', self)
        highlightnode.setStatusTip('Highlight node')
        highlightnode.triggered.connect(self.highlightNode)

        highlightelement = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_highlight_element.png'),'Element', self)
        highlightelement.setStatusTip('Highlight element')
        highlightelement.triggered.connect(self.highlightElement)

        highlightnodeset = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_highlight_nodeset.png'),'Highlight Nodeset', self)
        highlightnodeset.setStatusTip('Highlight nodeset')
        highlightnodeset.triggered.connect(self.highlightNodeSet)

        highlightelementset = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_highlight_elementset.png'),'Highlight Elementset', self)
        highlightelementset.setStatusTip('Highlight elementset')
        highlightelementset.triggered.connect(self.highlightElementSet)

        node = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_new_node.png'),'Create Node', self)
        node.setStatusTip('Create a new node')
        node.triggered.connect(self.createNode)

        movenodes = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_move_nodes.png'),'Move Nodes', self)
        movenodes.setStatusTip('Move the selected nodes')
        movenodes.triggered.connect(self.moveNodes)

        fusenodes = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_fuse_nodes.png'),'Fuse Nodes', self)
        fusenodes.setStatusTip('Fuse selected nodes within specified tolerance')
        fusenodes.triggered.connect(self.fuseNodes)

        pointmass = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_new_point_mass.png'),'Create Point Mass', self)
        pointmass.setStatusTip('Create a new point mass element')
        pointmass.triggered.connect(self.createPointMass)

        element = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_new_elements.png'),'Create Elements', self)
        element.setStatusTip('Create new elements')
        element.triggered.connect(self.createElements)

        extrude = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_extrude_elements.png'),'Extrude Elements', self)
        extrude.setStatusTip('Create new elements by extruding from selected elements')
        extrude.triggered.connect(self.extrudeElements)

        convert = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_convert_elements.png'),'Convert Elements', self)
        convert.setStatusTip('Convert elements from one type to another')
        convert.triggered.connect(self.convertElements)

        insert = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_insert_elements.png'),'Insert Elements', self)
        insert.setShortcut('I')
        insert.setStatusTip('Insert elements between selected nodes')
        insert.triggered.connect(self.insertElements)

        splitbeams = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_split_beam.png'),'Split Beams', self)
        splitbeams.setStatusTip('Split beam elements into smaller elements')
        splitbeams.triggered.connect(self.splitBeams)

        beamorient = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_beam_orient.png'),'Beam Orientation', self)
        beamorient.setStatusTip('Set orientation on beam elements')
        beamorient.triggered.connect(self.beamOrientation)

        newmesh = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_create_mesh.png'),'Create Mesh', self)
        newmesh.setStatusTip('Create a new mesh from scratch')
        newmesh.triggered.connect(self.createNewMesh)

        resizeelements = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_resize_elements.png'),'Resize Elements', self)
        resizeelements.setStatusTip('Resize selected elements')
        resizeelements.triggered.connect(self.resizeElements)

        getinfo = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_get_info.png'),'Node/Element Info', self)
        getinfo.setStatusTip('Print out info about selected node/element')
        getinfo.triggered.connect(self.getNodeElementInfo)

        copynodes = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_copy_nodes.png'),'Copy Nodes', self)
        copynodes.setStatusTip('Copy nodes and offset')
        copynodes.triggered.connect(self.copyNodes)

        copy = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_copy.png'),'Copy Elements', self)
        copy.setStatusTip('Copy elements and offset')
        copy.triggered.connect(self.copyElements)

        mirrorcopy = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_mirror_copy.png'),'Mirror Elements', self)
        mirrorcopy.setStatusTip('Mirror copy selected elements about x-, y- or z-plane')
        mirrorcopy.triggered.connect(self.mirrorCopyElements)

        move = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_move_elements.png'),'Move Elements', self)
        move.setStatusTip('Move all elements in elementset by specified coordinates')
        move.triggered.connect(self.moveElements)

        rotate = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_rotate_elements.png'),'Rotate Elements', self)
        rotate.setStatusTip('Rotate all elements in elementset around specified axis')
        rotate.triggered.connect(self.rotateElements)

        distance = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_distance.png'),'Measure Distance', self)
        distance.setStatusTip('Measure distance between two selected nodes')
        distance.triggered.connect(self.measureDistance)

        renumberNodes = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_link.png'),'Renumber Nodes', self)
        renumberNodes.setStatusTip('Renumber nodes from input')
        renumberNodes.triggered.connect(self.renumberNodes)

        renumberElements = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_link.png'),'Renumber Elements', self)
        renumberElements.setStatusTip('Renumber elements from input')
        renumberElements.triggered.connect(self.renumberElements)

        nodeset = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_create_nodeset.png'),'Create Nodeset', self)
        nodeset.setStatusTip('Create a nodeset')
        nodeset.triggered.connect(self.createNodeset)

        elementset = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_create_elementset.png'),'Create Elementset', self)
        elementset.setStatusTip('Create an elementset')
        elementset.triggered.connect(self.createElementset)

        material = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_new_material.png'),'Create Material', self)
        material.setStatusTip('Create a material')
        material.triggered.connect(self.newMaterial)

        section = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_create_section.png'),'Create Section', self)
        section.setStatusTip('Create a Section')
        section.triggered.connect(self.newSection)

        beamsection = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_beam_section.png'),'Modify Beam Section', self)
        beamsection.setStatusTip('Modify Beam Section')
        beamsection.triggered.connect(self.newBeamSection)

        applySection = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_apply_section.png'),'Apply Section', self)
        applySection.setStatusTip('Apply section to elementset')
        applySection.triggered.connect(self.applySection)

        static = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_static.png'),'Static', self)
        static.setStatusTip('Static solver')
        static.triggered.connect(self.solveStatic)

        eigenmodes = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_eigenmodes.png'),'Eigenmodes', self)
        eigenmodes.setStatusTip('Eigenmodes solver')
        eigenmodes.triggered.connect(self.solveEigenmodes)

        modal = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_modal_dynamic.png'),'Modal Dynamic', self)
        modal.setStatusTip('Modal Dynamic solver')
        modal.triggered.connect(self.solveModalDynamic)

        plastic = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_static_plastic.png'),'Static Plastic', self)
        plastic.setStatusTip('Static Plastic solver')
        plastic.triggered.connect(self.solveStaticPlastic)

        heattransfer = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_heat_transfer.png'),'Heat Transfer', self)
        heattransfer.setStatusTip('Heat Transfer solver')
        heattransfer.triggered.connect(self.solveHeatTransfer)

        touchlock = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_touch_lock.png'),'Touch Lock', self)
        touchlock.setStatusTip('Apply touchlock constraints between two nodesets')
        touchlock.triggered.connect(self.applyTouchLockConstraint)

        spiderlock = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_spider.png'),'Spider', self)
        spiderlock.setStatusTip('Create a spider between one node and a set of nodes')
        spiderlock.triggered.connect(self.createSpider)
        
        uniformload = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_uniform_load.png'),'Uniform Load', self)
        uniformload.setStatusTip('Apply uniform load to nodeset')
        uniformload.triggered.connect(self.applyUniformLoad)

        concentratedload = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_concentrated_load.png'),'Concentrated Load', self)
        concentratedload.setStatusTip('Apply concentrated load to nodeset')
        concentratedload.triggered.connect(self.applyConcentratedLoad)
        
        distributedload = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_distributed_load.png'),'Distributed Load', self)
        distributedload.setStatusTip('Apply distributed load to elementset (beam elements only)')
        distributedload.triggered.connect(self.applyDistributedLoad)

        torqueload = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_torque_load.png'),'Torque Load', self)
        torqueload.setStatusTip('Apply torque load to node (beam elements only)')
        torqueload.triggered.connect(self.applyTorqueLoad)

        gravityload = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_gravity_load.png'),'Gravity Load', self)
        gravityload.setStatusTip('Apply gravity load to elementset')
        gravityload.triggered.connect(self.applyGravityLoad)

        dynamicload = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_dynamic_load.png'),'Dynamic Load', self)
        dynamicload.setStatusTip('Apply dynamic load to nodeset')
        dynamicload.triggered.connect(self.applyDynamicLoad)

        displacement = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_boundary.png'),'Displacement', self)
        displacement.setStatusTip('Apply displacement boundary to nodeset')
        displacement.triggered.connect(self.applyDisplacement)

        newsolfile = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_new_file.png'),'Write sol-file...', self)
        newsolfile.setStatusTip('Write a sol-file that can be run directly in solver')
        newsolfile.triggered.connect(self.newSolFile)

        runsolver = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_settings.png'),'Run sol-file in solver...', self)
        runsolver.setStatusTip('Run sol-file in the solver to generate results')
        runsolver.triggered.connect(self.runSolver)

        scalefactor = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_scale_factor.png'),'Scale Factor', self)
        scalefactor.setStatusTip('Set the scale factor for view of displacements')
        scalefactor.triggered.connect(self.setScaleFactor)

        scalediagram = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_scale_diagram.png'),'Scale Diagram', self)
        scalediagram.setStatusTip('Adjust the scale for shear and bending moment diagrams')
        scalediagram.triggered.connect(self.setShearBendingDiagram)

        average = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_average.png'),'Average Stress/Strain', self)
        average.setStatusTip('View average stresses and strains')
        average.triggered.connect(self.viewAverage)

        animationspeed = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_animation_speed.png'),'Speed', self)
        animationspeed.setStatusTip('Set the frame to frame speed for animations')
        animationspeed.triggered.connect(self.setAnimationSpeed)

        animationonoff = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_start_pause.png'),'Play/Pause', self)
        animationonoff.setShortcut('P')
        animationonoff.setStatusTip('Turn animation on or off')
        animationonoff.triggered.connect(self.setAnimationOnOff)

        previousmode = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_previous.png'),'Previous Mode', self)
        previousmode.setStatusTip('Change to previous eigenmode')
        previousmode.triggered.connect(self.previousEigenmode)

        nextmode = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_next.png'),'Next Mode', self)
        nextmode.setStatusTip('Change to next eigenmode')
        nextmode.triggered.connect(self.nextEigenmode)

        empty1 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty1.setEnabled(False)
        empty2 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty2.setEnabled(False)
        empty3 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty3.setEnabled(False)
        empty4 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty4.setEnabled(False)
        empty5 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty5.setEnabled(False)
        empty6 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty6.setEnabled(False)
        empty7 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty7.setEnabled(False)
        empty8 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty8.setEnabled(False)
        empty9 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty9.setEnabled(False)
        empty10 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty10.setEnabled(False)
        empty11 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty11.setEnabled(False)
        empty12 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty12.setEnabled(False)
        empty13 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty13.setEnabled(False)
        empty14 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty14.setEnabled(False)
        empty15 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty15.setEnabled(False)
        empty16 = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_empty.png'),' ', self)
        empty16.setEnabled(False)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(newsession)
        fileMenu.addAction(importfrom)
        fileMenu.addAction(export)
        fileMenu.addAction(openfile)
        fileMenu.addAction(savefile)
        fileMenu.addAction(quicksave)
        fileMenu.addAction(exit)
        editMenu = menubar.addMenu('&Edit')
        editMenu.addAction(selectnodes)
        editMenu.addAction(selectelements)
        editMenu.addAction(selectlines)
        editMenu.addAction(selectfaces)
        editMenu.addAction(delete)
        viewMenu = menubar.addMenu('&View')
        viewMenu.addAction(resetview)
        viewMenu.addAction(centerview)
        viewMenu.addAction(wireframe)
        viewMenu.addAction(shaded)
        viewMenu.addAction(nodesview)
        viewMenu.addAction(origin)
        viewMenu.addAction(meshtree)
        highlightMenu = viewMenu.addMenu('&Highlight')
        highlightMenu.addAction(highlightnode)
        highlightMenu.addAction(highlightelement)
        highlightMenu.addAction(highlightnodeset)
        highlightMenu.addAction(highlightelementset)
        presetViewMenu = viewMenu.addMenu('&Views')
        presetViewMenu.addAction(viewleft)
        presetViewMenu.addAction(viewfront)
        presetViewMenu.addAction(viewright)
        presetViewMenu.addAction(viewback)
        presetViewMenu.addAction(viewover)
        presetViewMenu.addAction(viewunder)
        presetViewMenu.addAction(viewangled)
        geometryMenu = menubar.addMenu('&Geometry')
        geometryMenu.addAction(partview)
        geometryMenu.addAction(seedlines)
        geometryMenu.addAction(seedpart)
        geometryMenu.addAction(meshfaces)
        geometryMenu.addAction(meshpart)
        meshMenu = menubar.addMenu('&Mesh')
        meshMenu.addAction(newmesh)
        meshMenu.addAction(getinfo)
        nodeMenu = meshMenu.addMenu('&Nodes')
        nodeMenu.addAction(node)
        nodeMenu.addAction(copynodes)
        nodeMenu.addAction(movenodes)
        nodeMenu.addAction(fusenodes)
        nodeMenu.addAction(distance)
        nodeMenu.addAction(nodeset)
        nodeMenu.addAction(renumberNodes)
        elementMenu = meshMenu.addMenu('&Elements')
        elementMenu.addAction(element)
        elementMenu.addAction(extrude)
        elementMenu.addAction(pointmass)
        elementMenu.addAction(insert)
        elementMenu.addAction(convert)
        elementMenu.addAction(copy)
        elementMenu.addAction(mirrorcopy)
        elementMenu.addAction(move)
        elementMenu.addAction(rotate)
        elementMenu.addAction(resizeelements)
        elementMenu.addAction(splitbeams)
        elementMenu.addAction(beamorient)
        elementMenu.addAction(elementset)
        elementMenu.addAction(renumberElements)
        materialMenu = meshMenu.addMenu('&Materials')
        materialMenu.addAction(material)
        sectionMenu = meshMenu.addMenu('&Sections')
        sectionMenu.addAction(section)
        sectionMenu.addAction(applySection)
        sectionMenu.addAction(beamsection)
        solverMenu = menubar.addMenu('&Solution')
        solverMenu.addAction(solutionview)
        solutionMenu = solverMenu.addMenu('&Solutions')
        solutionMenu.addAction(static)
        solutionMenu.addAction(eigenmodes)
        solutionMenu.addAction(modal)
        solutionMenu.addAction(plastic)
        solutionMenu.addAction(heattransfer)
        constraintsMenu = solverMenu.addMenu('&Constraints')
        constraintsMenu.addAction(touchlock)
        constraintsMenu.addAction(spiderlock)
        boundMenu = solverMenu.addMenu('&Boundaries')
        boundMenu.addAction(displacement)
        loadsMenu = solverMenu.addMenu('&Loads')
        loadsMenu.addAction(uniformload)
        loadsMenu.addAction(concentratedload)
        loadsMenu.addAction(distributedload)
        loadsMenu.addAction(torqueload)
        loadsMenu.addAction(gravityload)
        loadsMenu.addAction(dynamicload)
        solverMenu.addAction(newsolfile)
        solverMenu.addAction(runsolver)
        resultMenu = menubar.addMenu('&Result')
        resultMenu.addAction(resultview)
        resultMenu.addAction(scalefactor)
        resultMenu.addAction(scalediagram)
        resultMenu.addAction(average)
        aniMenu = resultMenu.addMenu('&Animation')
        aniMenu.addAction(animationonoff)
        aniMenu.addAction(previousmode)
        aniMenu.addAction(nextmode)
        aniMenu.addAction(animationspeed)

        main_toolbar1 = QtWidgets.QToolBar('Main Toolbar Upper')
        main_toolbar1.setIconSize(QtCore.QSize(24,24))
        main_toolbar1.setMovable(False)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(empty1)
        main_toolbar1.addAction(partview)
        main_toolbar1.addAction(solutionview)
        main_toolbar1.addAction(resultview)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(resetview)
        main_toolbar1.addAction(shaded)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(highlightnode)
        main_toolbar1.addAction(highlightelement)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(newmesh)
        main_toolbar1.addAction(element)
        main_toolbar1.addAction(extrude)
        main_toolbar1.addAction(insert)
        main_toolbar1.addAction(convert)
        main_toolbar1.addAction(resizeelements)
        main_toolbar1.addAction(rotate)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(nodeset)
        main_toolbar1.addAction(elementset)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(material)
        main_toolbar1.addAction(beamsection)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(static)
        main_toolbar1.addAction(plastic)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(touchlock)
        main_toolbar1.addAction(spiderlock)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(uniformload)
        main_toolbar1.addAction(concentratedload)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(average)
        main_toolbar1.addAction(animationonoff)
        main_toolbar1.addAction(animationspeed)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()
        main_toolbar1.addAction(newsolfile)
        main_toolbar1.addSeparator()
        main_toolbar1.addSeparator()

        main_toolbar2 = QtWidgets.QToolBar('Main Toolbar Middle')
        main_toolbar2.setIconSize(QtCore.QSize(24,24))
        main_toolbar2.setMovable(False)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(selectlines)
        main_toolbar2.addAction(selectfaces)
        main_toolbar2.addAction(selectnodes)
        main_toolbar2.addAction(selectelements)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(centerview)
        main_toolbar2.addAction(wireframe)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(highlightnodeset)
        main_toolbar2.addAction(highlightelementset)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(node)
        main_toolbar2.addAction(fusenodes)
        main_toolbar2.addAction(copynodes)
        main_toolbar2.addAction(movenodes)
        main_toolbar2.addAction(copy)
        main_toolbar2.addAction(mirrorcopy)
        main_toolbar2.addAction(move)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(hide)
        main_toolbar2.addAction(show)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(section)
        main_toolbar2.addAction(applySection)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(eigenmodes)
        main_toolbar2.addAction(modal)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(displacement)
        main_toolbar2.addAction(empty2)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(distributedload)
        main_toolbar2.addAction(torqueload)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(scalefactor)
        main_toolbar2.addAction(previousmode)
        main_toolbar2.addAction(nextmode)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()
        main_toolbar2.addAction(runsolver)
        main_toolbar2.addSeparator()
        main_toolbar2.addSeparator()

        main_toolbar3 = QtWidgets.QToolBar('Main Toolbar Lower')
        main_toolbar3.setIconSize(QtCore.QSize(24,24))
        main_toolbar3.setMovable(False)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(delete)
        main_toolbar3.addAction(empty3)
        main_toolbar3.addAction(empty4)
        main_toolbar3.addAction(empty5)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(origin)
        main_toolbar3.addAction(nodesview)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(getinfo)
        main_toolbar3.addAction(distance)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(seedlines)
        main_toolbar3.addAction(seedpart)
        main_toolbar3.addAction(meshfaces)
        main_toolbar3.addAction(meshpart)
        main_toolbar3.addAction(pointmass)
        main_toolbar3.addAction(beamorient)
        main_toolbar3.addAction(splitbeams)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(empty6)
        main_toolbar3.addAction(empty7)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(empty8)
        main_toolbar3.addAction(empty9)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(heattransfer)
        main_toolbar3.addAction(empty10)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(empty11)
        main_toolbar3.addAction(empty12)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(gravityload)
        main_toolbar3.addAction(dynamicload)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(scalediagram)
        main_toolbar3.addAction(empty13)
        main_toolbar3.addAction(empty14)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()
        main_toolbar3.addAction(empty15)
        main_toolbar3.addSeparator()
        main_toolbar3.addSeparator()

        btnViewAssembly = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_view_assembly.png'),'Switch view between Assembly and Part', self)
        btnViewAssembly.triggered.connect(self.btnViewAssemblyAction)
        btnViewAssembly.setStatusTip('Switch view between Assembly and Part')
        btnViewGeometry = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_view_geometry.png'),'Change view to Geometry', self)
        btnViewGeometry.triggered.connect(self.btnViewGeometryAction)
        btnViewGeometry.setStatusTip('Change view to Geometry')
        btnViewMesh = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_view_mesh.png'),'Change view to Mesh', self)
        btnViewMesh.triggered.connect(self.btnViewMeshAction)
        btnViewMesh.setStatusTip('Change view to Mesh')
        btnViewConstraint = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_view_constraint.png'),'Change view to Constraints', self)
        btnViewConstraint.triggered.connect(self.btnViewConstraintAction)
        btnViewConstraint.setStatusTip('Change view to Constraints')
        btnViewBoundary = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_view_boundary.png'),'Change view to Boundaries', self)
        btnViewBoundary.triggered.connect(self.btnViewBoundaryAction)
        btnViewBoundary.setStatusTip('Change view to Boundaries')
        btnViewLoad = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_view_load.png'),'Change view to Loads', self)
        btnViewLoad.triggered.connect(self.btnViewLoadAction)
        btnViewLoad.setStatusTip('Change view to Loads')
        btnViewSolution = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_view_solution.png'),'Change view to Solutions', self)
        btnViewSolution.triggered.connect(self.btnViewSolutionAction)
        btnViewSolution.setStatusTip('Change view to Solution')
        btnViewResult = QtWidgets.QAction(QtGui.QIcon('../Icons/icon_view_result.png'),'Change view to Results', self)
        btnViewResult.triggered.connect(self.btnViewResultAction)
        btnViewResult.setStatusTip('Change view to Result')

        view_toolbar = QtWidgets.QToolBar('View Toolbar')
        view_toolbar.setIconSize(QtCore.QSize(48,48))
        view_toolbar.setMovable(False)
        view_toolbar.addAction(empty1)
        view_toolbar.addAction(empty2)
        view_toolbar.addSeparator()
        view_toolbar.addAction(btnViewAssembly)
        view_toolbar.addSeparator()
        view_toolbar.addAction(btnViewGeometry)
        view_toolbar.addSeparator()
        view_toolbar.addAction(btnViewMesh)
        view_toolbar.addSeparator()
        view_toolbar.addAction(btnViewConstraint)
        view_toolbar.addSeparator()
        view_toolbar.addAction(btnViewBoundary)
        view_toolbar.addSeparator()
        view_toolbar.addAction(btnViewLoad)
        view_toolbar.addSeparator()
        view_toolbar.addAction(btnViewSolution)
        view_toolbar.addSeparator()
        view_toolbar.addAction(btnViewResult)
        view_toolbar.addSeparator()

        self.addToolBar(QtCore.Qt.LeftToolBarArea, view_toolbar)

        parentWidget = QtWidgets.QWidget()
        tool1 = QtWidgets.QToolBar()
        tool1.addWidget(main_toolbar1)
        tool2 = QtWidgets.QToolBar()
        tool2.addWidget(main_toolbar2)
        tool3 = QtWidgets.QToolBar()
        tool3.addWidget(main_toolbar3)

        vbox = QtWidgets.QVBoxLayout()
        vbox.setContentsMargins(0,5,5,5)
        vbox.setSpacing(0)

        vbox.addWidget(tool1)
        vbox.addWidget(tool2)
        vbox.addWidget(tool3)

        self.viewer.setSizePolicy( QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding )
        vbox.addWidget(self.viewer)

        parentWidget.setLayout(vbox)
        self.setCentralWidget(parentWidget)

        self.centerWindow()
        self.resize(1200,800)
        self.viewAngled()


    def centerWindow(self):
        '''
    Center Window on the Current Screen,
    with Multi-Monitor support
    '''
        self.showNormal()
        window_geometry = self.frameGeometry()
        self.resize(int(QtWidgets.QDesktopWidget().screenGeometry().width() // 1.3),
                    int(QtWidgets.QDesktopWidget().screenGeometry().height() // 1.5))
        mousepointer_position = QtWidgets.QApplication.desktop().cursor().pos()
        screen = QtWidgets.QApplication.desktop().screenNumber(mousepointer_position)
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        window_geometry.moveCenter(centerPoint)
        return bool(not self.move(window_geometry.topLeft()))

    

    # File operations
    # ---------------
    def clearModel(self):
        '''
    Clears the entire model and starts
    with a new one.
    '''
        print('NOT READY!')
    def openFile(self):
        print('NOT READY!')
    def importFrom(self):
        '''
    Import meshes, solutions, boundaries, loads, 
    elementsets and nodesets from a *.sol file,
    *.out file or *.mesh file with the help of PyQt's
    built-in QFileDialog class. 
    '''
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', './')[0]
        part_name = 'part-'+str(len(self.model.parts)+1)
        self.model.parts[part_name] = Part(part_name)
        self.model.parts[part_name].mesher = self.model.mesher
        
        if filename[-4:] == '.out':
            print('\n\tImport mesh from *.out file not ready')
        elif filename[-4:] == '.sol':
            self.model.parts[part_name].importMesh(filename)
        elif filename[-4:] == '.bdf':
            self.model.parts[part_name].importMesh(filename)
        elif filename[-4:] == '.inp':
            self.model.parts[part_name].importMesh(filename)
        elif filename[-4:] == '.dat':
            self.model.parts[part_name].importMesh(filename)
        elif filename[-4:] in ['.stp', '.STP']:
            self.model.parts[part_name].readStepFile(filename)
        elif filename[-5:] in ['.step', '.STEP']:
            self.model.parts[part_name].readStepFile(filename)
        elif filename[-4:] == '.mdl':
            print('\n\tImport from *.mdl file not ready.')
        else:
            print('\n\tUnknown file type. Accepted files are:')
            print('\t*.sol, *.bdf, *.inp, *.dat, *.STP, *.STEP, *.mdl')

        self.model.setCurrentPart(part_name)
        self.centerModel()
        self.viewAngled()

    def exportMesh(self):
        print('NOT READY!')
    def saveFile(self):
        print('NOT READY!')
    def quickSaveFile(self):
        print('NOT READY!')



    # Specify what can be selected
    # using the mouse pointer
    # -----------------------------
    def selectNodes(self):
        '''
    Set mouse pointer to select nodes.
    '''
        if len(self.model.selected_elements) != 0:
            self.model.selected_elements.clear()
            self.viewer.update()
        self.model.selectOption = 'Nodes'
        self.statusBar().showMessage('  Selecting... nodes' )

    def selectElements(self):
        '''
    Set mouse pointer to select elements.
    '''
        if len(self.model.selected_nodes) != 0:
            self.model.selected_nodes.clear()
            self.viewer.update()
        self.model.selectOption = 'Elements'
        self.statusBar().showMessage('  Selecting... elements' )        

    def selectLines(self):
        '''
    Set mouse pointer to select lines.
    '''
        if len(self.model.selected_lines) != 0:
            self.model.selected_lines.clear()
            self.viewer.update()
        self.model.selectOption = 'Lines'
        self.statusBar().showMessage('  Selecting... lines' )

    def selectFaces(self):
        '''
    Set mouse pointer to select faces.
    '''
        if len(self.model.selected_faces) != 0:
            self.model.selected_faces.clear()
            self.viewer.update()
        self.model.selectOption = 'Faces'
        self.statusBar().showMessage('  Selecting... faces' )



    def getNodeElementInfo(self):
        for l in self.model.selected_lines:
            line = self.model.selected_lines[l]
            print('\nline number:', l)
            print('type:', line.type)
            print('length:', line.length)
            if line.type == 'arc':
                print('radius:', line.radius)
                print('axis:', line.axis)
            print('points:', line.points)


    def deleteItem(self):
        print('NOT READY!')




    # Change how the mesh and geometry is viewed
    # ------------------------------------------
    def nodesView(self):
        '''
    Toggle visibility of nodes on and off
    '''
        if self.viewer.viewNodes == True:
            self.viewer.viewNodes = False
        else:
            self.viewer.viewNodes = True
        self.viewer.update()
    def wireframeView(self):
        '''
    Change to wireframe view.
    '''
        self.viewer.viewShaded = False
        self.viewer.viewWireframe = True
        self.viewer.update()
    def shadedView(self):
        '''
    Change to shaded view.
    '''
        self.viewer.viewShaded = True
        self.viewer.viewWireframe = False
        self.viewer.update()
    def showOrigin(self):
        '''
    Wether or not to show an RGB triad at the 
    origin (coordinates 0., 0., 0.).
    '''
        if self.viewer.viewOrigin == False:
            self.viewer.viewOrigin = True
        else:
            self.viewer.viewOrigin = False
        self.viewer.update()

    def selectPart(self):
        print('NOT READY!')
    def selectSolution(self):
        print('NOT READY!')
    def selectResult(self):
        print('NOT READY!')

    def showElements(self):
        print('NOT READY!')
    def hideElements(self):
        print('NOT READY!')
    def showMeshTree(self):
        print('NOT READY!')
    def viewAssembly(self):
        print('NOT READY!')



    # Camera manipulations
    # --------------------
    def keyPressEvent(self, e):
        '''
    Capture the user holding down CTRL, SHIFT
    or ALT, to be used while clicking the mouse
    for manipulation of the camera.
    '''
        if e.key() == QtCore.Qt.Key_Control:
            self.viewer.activeCTRL = True
        if e.key() == QtCore.Qt.Key_Shift:
            self.viewer.activeSHIFT = True
        if e.key() == QtCore.Qt.Key_Alt:
            self.viewer.activeALT = True
    def keyReleaseEvent(self, e):
        '''
    Capture the user releasing CTRL, SHIFT
    or ALT, to be used while clicking the mouse
    for manipulation of the camera.
    '''
        if e.key() == QtCore.Qt.Key_Control:
            self.viewer.activeCTRL = False
        if e.key() == QtCore.Qt.Key_Shift:
            self.viewer.activeSHIFT = False
        if e.key() == QtCore.Qt.Key_Alt:
            self.viewer.activeALT = False
    def centerModel(self):
        '''
    Center the model in the viewer.
    '''
        if self.viewer.modelCentered == True:
            self.viewer.modelCentered = False
            self.viewer.coordSys0_centered = self.viewer.coordSys0
        else:
            self.viewer.modelCentered = True
            self.viewer.coordSys0_centered = CoordSys3D( \
                 Point3D(-(self.viewer.currentDisplayList['view_scope']['max'][0]+self.viewer.currentDisplayList['view_scope']['min'][0])/2.,
                         -(self.viewer.currentDisplayList['view_scope']['max'][1]+self.viewer.currentDisplayList['view_scope']['min'][1])/2.,
                         -(self.viewer.currentDisplayList['view_scope']['max'][2]+self.viewer.currentDisplayList['view_scope']['min'][2])/2.),
                     Vector3D(1.,0.,0.),Vector3D(0.,1.,0.))
    def viewLeft(self):
        '''
    View model from the left.
    '''
        tangent = math.tan( self.viewer.camera.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        distanceFromTarget = self.viewer.camera.sceneRadius / tangent
        print('\n\tView Left')
        self.viewer.camera.position = Point3D(-distanceFromTarget*0.707,0,0)
        self.viewer.camera.target = Point3D(0,0,0)
        self.viewer.camera.up = self.viewer.camera.ground.returnCopy()
    def viewFront(self):
        '''
    View model from the front.
    '''
        tangent = math.tan( self.viewer.camera.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        distanceFromTarget = self.viewer.camera.sceneRadius / tangent
        print('\n\tView Front')
        self.viewer.camera.position = Point3D(0,0,distanceFromTarget*0.707)
        self.viewer.camera.target = Point3D(0,0,0)
        self.viewer.camera.up = self.viewer.camera.ground.returnCopy()
    def viewRight(self):
        '''
    View model from the right.
    '''
        tangent = math.tan( self.viewer.camera.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        distanceFromTarget = self.viewer.camera.sceneRadius / tangent
        print('\n\tView Right')
        self.viewer.camera.position = Point3D(distanceFromTarget*0.707,0,0)
        self.viewer.camera.target = Point3D(0,0,0)
        self.viewer.camera.up = self.viewer.camera.ground.returnCopy()
    def viewBack(self):
        '''
    View model from the back.
    '''
        tangent = math.tan( self.viewer.camera.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        distanceFromTarget = self.viewer.camera.sceneRadius / tangent
        print('\n\tView Back')
        self.viewer.camera.position = Point3D(0,0,-distanceFromTarget*0.707)
        self.viewer.camera.target = Point3D(0,0,0)
        self.viewer.camera.up = self.viewer.camera.ground.returnCopy()
    def viewOver(self):
        '''
    View model from above.
    '''
        tangent = math.tan( self.viewer.camera.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        distanceFromTarget = self.viewer.camera.sceneRadius / tangent
        print('\n\tView Over')
        self.viewer.camera.position = Point3D(0,distanceFromTarget*0.707,0)
        self.viewer.camera.target = Point3D(0,0,0)
        self.viewer.camera.up = Vector3D(0,0,-1)
    def viewUnder(self):
        '''
    View model from under.
    '''
        tangent = math.tan( self.viewer.camera.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        distanceFromTarget = self.viewer.camera.sceneRadius / tangent
        print('\n\tView Under')
        self.viewer.camera.position = Point3D(0,-distanceFromTarget*0.707,0)
        self.viewer.camera.target = Point3D(0,0,0)
        self.viewer.camera.up = Vector3D(0,0,1)
    def viewAngled(self):
        '''
    View model at an angle.
    '''
        tangent = math.tan( self.viewer.camera.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        distanceFromTarget = self.viewer.camera.sceneRadius / tangent
        print('\n\tView Angled')
        self.viewer.camera.position = Point3D(0.505*distanceFromTarget,0.505*distanceFromTarget,0.505*distanceFromTarget)
        self.viewer.camera.target = Point3D(0,0,0)
        self.viewer.camera.up = self.viewer.camera.ground.returnCopy()








    def seedLines(self):
        print('NOT READY!')
    def seedPart(self):
        print('NOT READY!')
    def meshFaces(self):
        print('NOT READY!')
    def meshPart(self):
        print('NOT READY!')

    def highlightNode(self):
        print('NOT READY!')
    def highlightElement(self):
        print('NOT READY!')
    def highlightNodeSet(self):
        print('NOT READY!')
    def highlightElementSet(self):
        print('NOT READY!')
    def createNode(self):
        print('NOT READY!')
    def moveNodes(self):
        print('NOT READY!')
    def fuseNodes(self):
        print('NOT READY!')
    def createPointMass(self):
        print('NOT READY!')
    def createElements(self):
        print('NOT READY!')
    def extrudeElements(self):
        print('NOT READY!')
    def convertElements(self):
        print('NOT READY!')
    def insertElements(self):
        print('NOT READY!')
    def splitBeams(self):
        print('NOT READY!')
    def beamOrientation(self):
        print('NOT READY!')
    def createNewMesh(self):
        print('NOT READY!')
    def resizeElements(self):
        print('NOT READY!')
    def copyNodes(self):
        print('NOT READY!')
    def copyElements(self):
        print('NOT READY!')
    def mirrorCopyElements(self):
        print('NOT READY!')
    def moveElements(self):
        print('NOT READY!')
    def rotateElements(self):
        print('NOT READY!')
    def measureDistance(self):
        print('NOT READY!')
    def renumberNodes(self):
        print('NOT READY!')
    def renumberElements(self):
        print('NOT READY!')
    def createNodeset(self):
        print('NOT READY!')
    def createElementset(self):
        print('NOT READY!')

        
    def newMaterial(self):
        print('NOT READY!')
    def newSection(self):
        print('NOT READY!')
    def newBeamSection(self):
        print('NOT READY!')
    def applySection(self):
        print('NOT READY!')


    def solveStatic(self):
        print('NOT READY!')
    def solveEigenmodes(self):
        print('NOT READY!')
    def solveModalDynamic(self):
        print('NOT READY!')
    def solveStaticPlastic(self):
        print('NOT READY!')
    def solveHeatTransfer(self):
        print('NOT READY!')
    def applyTouchLockConstraint(self):
        print('NOT READY!')
    def createSpider(self):
        print('NOT READY!')
    def applyUniformLoad(self):
        print('NOT READY!')
    def applyConcentratedLoad(self):
        print('NOT READY!')
    def applyDistributedLoad(self):
        print('NOT READY!')
    def applyTorqueLoad(self):
        print('NOT READY!')
    def applyGravityLoad(self):
        print('NOT READY!')
    def applyDynamicLoad(self):
        print('NOT READY!')
    def applyDisplacement(self):
        print('NOT READY!')
    def newSolFile(self):
        print('NOT READY!')
    def runSolver(self):
        print('NOT READY!')
    def setScaleFactor(self):
        print('NOT READY!')
    def setShearBendingDiagram(self):
        print('NOT READY!')
    def viewAverage(self):
        print('NOT READY!')
    def setAnimationSpeed(self):
        print('NOT READY!')
    def setAnimationOnOff(self):
        print('NOT READY!')
    def previousEigenmode(self):
        print('NOT READY!')
    def nextEigenmode(self):
        print('NOT READY!')

    def btnViewAssemblyAction(self):
        if self.viewer.viewAssembly == True:
            self.viewer.viewAssembly = False
        else:
            self.viewer.viewAssembly = True
        self.viewer.updateDisplayList()
        self.statusBar().showMessage('  ASSEMBLY  ')
        glClearColor(self.viewer.colors['background_pre'][0], 
                     self.viewer.colors['background_pre'][1], 
                     self.viewer.colors['background_pre'][2], 
                     self.viewer.colors['background_pre'][3])
        glClearDepth(1.0)
        self.viewer.update()

    def btnViewGeometryAction(self):
        self.viewer.viewGeometry = True
        self.viewer.viewMesh = False
        self.viewer.viewResults = False
        self.viewer.updateDisplayList()
        self.statusBar().showMessage('  GEOMETRY  ')
        glClearColor(self.viewer.colors['background_pre'][0], 
                     self.viewer.colors['background_pre'][1], 
                     self.viewer.colors['background_pre'][2], 
                     self.viewer.colors['background_pre'][3])
        glClearDepth(1.0)
        self.viewer.update()
        
    def btnViewMeshAction(self):
        self.viewer.viewGeometry = False
        self.viewer.viewMesh = True
        self.viewer.viewResults = False
        self.viewer.updateDisplayList()
        self.statusBar().showMessage('  MESH  ')
        glClearColor(self.viewer.colors['background_pre'][0], 
                     self.viewer.colors['background_pre'][1], 
                     self.viewer.colors['background_pre'][2], 
                     self.viewer.colors['background_pre'][3])
        glClearDepth(1.0)
        self.viewer.update()
        
    def btnViewConstraintAction(self):
        self.viewer.viewResults = False
        self.viewer.viewBoundaries = False
        self.viewer.viewLoads = False
        self.viewer.viewConstraints = True
        self.viewer.viewSolutions = False
        self.viewer.updateDisplayList()
        self.statusBar().showMessage('  CONSTRAINTS  ')
        glClearColor(self.viewer.colors['background_pre'][0], 
                     self.viewer.colors['background_pre'][1], 
                     self.viewer.colors['background_pre'][2], 
                     self.viewer.colors['background_pre'][3])
        glClearDepth(1.0)
        self.viewer.update()
        
    def btnViewBoundaryAction(self):
        self.viewer.viewResults = False
        self.viewer.viewBoundaries = True
        self.viewer.viewLoads = False
        self.viewer.viewConstraints = False
        self.viewer.viewSolutions = False
        self.viewer.updateDisplayList()
        self.statusBar().showMessage('  BOUNDARIES  ')
        glClearColor(self.viewer.colors['background_pre'][0], 
                     self.viewer.colors['background_pre'][1], 
                     self.viewer.colors['background_pre'][2], 
                     self.viewer.colors['background_pre'][3])
        glClearDepth(1.0)
        self.viewer.update()
        
    def btnViewLoadAction(self):
        self.viewer.viewResults = False
        self.viewer.viewBoundaries = False
        self.viewer.viewLoads = True
        self.viewer.viewConstraints = False
        self.viewer.viewSolutions = False
        self.viewer.updateDisplayList()
        self.statusBar().showMessage('  LOADS  ')
        glClearColor(self.viewer.colors['background_pre'][0], 
                     self.viewer.colors['background_pre'][1], 
                     self.viewer.colors['background_pre'][2], 
                     self.viewer.colors['background_pre'][3])
        glClearDepth(1.0)
        self.viewer.update()
        
    def btnViewSolutionAction(self):
        self.viewer.viewResults = False
        self.viewer.viewBoundaries = False
        self.viewer.viewLoads = False
        self.viewer.viewConstraints = False
        self.viewer.viewSolutions = True
        self.viewer.updateDisplayList()
        self.statusBar().showMessage('  SOLUTIONS  ')
        glClearColor(self.viewer.colors['background_pre'][0], 
                     self.viewer.colors['background_pre'][1], 
                     self.viewer.colors['background_pre'][2], 
                     self.viewer.colors['background_pre'][3])
        glClearDepth(1.0)
        self.viewer.update()
        
    def btnViewResultAction(self):
        self.viewer.viewGeometry = False
        self.viewer.viewMesh = False
        self.viewer.viewResults = True
        self.viewer.viewBoundaries = False
        self.viewer.viewLoads = False
        self.viewer.viewConstraints = False
        self.viewer.viewSolutions = False
        self.viewer.updateDisplayList()
        self.statusBar().showMessage('  RESULTS  ')
        glClearColor(self.viewer.colors['background_post'][0], 
                     self.viewer.colors['background_post'][1], 
                     self.viewer.colors['background_post'][2], 
                     self.viewer.colors['background_post'][3])
        glClearDepth(1.0)
        self.viewer.update()



























if __name__ == '__main__':

    app = QtWidgets.QApplication(['Finite Element Viewer'])
    window = Interface()
    window.show()
    app.exec_()

