#
#
#   mesher.py
#  ------------------
#   This is the mesher module. It has the ability to
#   generate meshes for individual parts, either from
#   scratch or from geometry. The generated meshes are
#   stored in part objects.
#

import sys
sys.path.insert(1, '../Objects')
from timeit import time
import random

import numpy as np
from scipy.spatial import Delaunay, KDTree
from scipy.interpolate import interp1d
from geometries import *

import matplotlib.pyplot as plt # for debugging





class Mesher(object):
    '''
Mesher module. Used for all meshing operations.
Takes lines, faces or parts as input and generates
seeds/meshes for them. Can also be used to manipulate
nodes and elements in existing meshes, or build new
meshes from scratch.
'''
    def __init__(self):
        pass



    def seedLine(self,line,element_size,weight=None):
        '''
    Creates seeds at coordinates between start point
    and end point of line. At least 2 seeds (start and
    stop), and the rest spaced out between those.
    '''
        points = [p.coordinates for p in line.points]
#        print(f'element_size: {element_size}')
        n_seeds = int(round(line.length/element_size))
#        print(f'points: {points}')
#        print(f'n_seeds: {n_seeds}')

        line.seeds = []
        if n_seeds <= 2:
            n_seeds = 2
            line.seeds.append(points[0])
            line.seeds.append(points[-1])
        else:
            line.seeds = self.interpolateCurve(points, n_seeds)
#        print(f'seeds: {seeds}')


    def meshFace(self,face,seeded=False):
        '''
    Creates a surface mesh for a given Face based on
    either the seeded edges or all the points on the
    edge lines of the face. It first converts the
    edge point coordinates to a 2D parameter space,
    and then uses bridson sampling and delaunay to
    generate the surface triangles.
    '''

        
        # -----------
        # mesh PLANE faces
        # -------------------------
        if face.type == 'plane':
            # move all points so centroid is at origin
            centered = np.empty((0, 3))
            centered_e = {}
            for e in face.edges:
                centered_e[e] = face.edges[e].points - face.centroid_np
                centered = np.vstack([centered,centered_e[e]])
            # flip points to x-y plane
            need_to_flip = False
            if face.normal_v == Vector3D(0,0,1) or face.normal_v == Vector3D(0,0,-1):
                flipped = centered
                flipped_e = centered_e
            else:
                need_to_flip = True
                flip_angle = self.angle_between(np.array(face.normal_v.coordinates),np.array([0.,0.,1.]))
                axis = np.cross(np.array(face.normal_v.coordinates),np.array([0.,0.,1.]))
                axis = self.unit_vector(axis)
                a0 = Point3D(0,0,0)
                a1 = Vector3D(axis[0],axis[1],axis[2])
                flipped = np.empty((0, 3))
                flipped_e = {}
                for e in face.edges:
                    flipped_e[e] = []
                    for p in centered_e[e]:
                        p1 = Point3D(p[0],p[1],p[2])
                        flipped_e[e].append(self.rotate_point_about_axis(p1, a0, a1, flip_angle))
                    flipped_e[e] = np.array(flipped_e[e])
                    flipped = np.vstack([flipped,flipped_e[e]])
            # generate points inside face with bridson sampling
            max_min = [np.min(flipped[:, 0]), np.min(flipped[:, 1]), np.max(flipped[:, 0]), np.max(flipped[:, 1])]
            cell_size = max_min[2] - max_min[0]
            for e in face.edges:
                for l in face.edges[e].lines:
                    if face.edges[e].lines[l].length < cell_size:
                        cell_size = face.edges[e].lines[l].length
            cell_size = cell_size*0.7
            outside_point = np.min(flipped, axis=0) - 23
            outside_point[0] -= 0.666*np.pi
            internal = self.bridson_sampling_2D(max_min, flipped_e, flipped, cell_size, outside_point)
            # genrate triangles from face points and remove those that are outside edges
            tri = Delaunay(internal[:,:2], qhull_options='Qt QJ')
            nodes, elements = self.keep_internal_triangles_2D(flipped_e, internal, tri.simplices, outside_point)
            # flip back from x-y plane
            for node in nodes:
                nodes[node] = Point3D(nodes[node][0], nodes[node][1], 0.)
                if need_to_flip:
                    nodes[node] = np.array(self.rotate_point_about_axis(nodes[node], a0, a1, -flip_angle))
                else:
                    nodes[node] = np.array(nodes[node].coordinates)
                # move all points back so the centroid is back where it was
                nodes[node] = nodes[node] + face.centroid_np
            if face.normal_v == Vector3D(0,0,-1):
                for elm in elements:
                    elements[elm] = elements[elm][::-1]
            face.g_mesh['nodes'] = nodes
            face.g_mesh['elements'] = elements


        # -----------
        # mesh CYLINDRICAL faces
        # -------------------------
        elif face.type == 'cylindrical':
            # add extra points on edge of cylinder so there are an
            # equal amount of points on that as there is on the arcs
            smallest_length = 10e5
            for e in face.edges:
                for p in range(len(face.edges[e].points)-1):
                    if np.linalg.norm(face.edges[e].points[p]-face.edges[e].points[p+1]) < smallest_length:
                        smallest_length = np.linalg.norm(face.edges[e].points[p]-face.edges[e].points[p+1])
            new_points_e = {}
            new_points = np.empty((0, 3))
            for e in face.edges:
                new_points_e[e] = []
                for p in range(len(face.edges[e].points)):
                    if abs(np.linalg.norm(face.edges[e].points[p]-face.edges[e].points[p-1]) - smallest_length) < 0.1:
                        new_points_e[e].append(face.edges[e].points[p])
                    else:
                        inserts = int(np.linalg.norm(face.edges[e].points[p]-face.edges[e].points[p-1])/smallest_length) + 1 
                        points = face.edges[e].points
                        for i in range(inserts-1):
                            new_points_e[e].append([points[p-1][0]+(i+1)*(points[p][0]-points[p-1][0])/inserts,
                                                    points[p-1][1]+(i+1)*(points[p][1]-points[p-1][1])/inserts,
                                                    points[p-1][2]+(i+1)*(points[p][2]-points[p-1][2])/inserts])
                        new_points_e[e].append(points[p])
                new_points_e[e] = np.array(new_points_e[e])
                new_points = np.vstack([new_points,new_points_e[e]])
            # move all points so centroid is at origin
            centered = np.empty((0, 3))
            centered_e = {}
            for e in face.edges:
                centered_e[e] = face.edges[e].points - np.array(face.normal[0].coordinates)
                centered_e[e] = new_points_e[e] - np.array(face.normal[0].coordinates)
                centered = np.vstack([centered,centered_e[e]])
            # flip base to x-y plane
            need_to_flip = False
            if face.normal_v == Vector3D(0,0,1) or face.normal_v == Vector3D(0,0,-1):
                flipped = centered
                flipped_e = centered_e
            else:
                need_to_flip = True
                flip_angle = self.angle_between(np.array(face.normal_v.coordinates),np.array([0.,0.,1.])) - np.pi
                flip_axis = np.cross(np.array(face.normal_v.coordinates),np.array([0.,0.,1.]))
                flip_axis = self.unit_vector(flip_axis)
                a0 = Point3D(0,0,0)
                a1 = Vector3D(flip_axis[0],flip_axis[1],flip_axis[2])
                flipped = np.empty((0, 3))
                flipped_e = {}
                for e in face.edges:
                    flipped_e[e] = []
                    for p in centered_e[e]:
                        p1 = Point3D(p[0],p[1],p[2])
                        flipped_e[e].append(self.rotate_point_about_axis(p1, a0, a1, flip_angle))
                    flipped_e[e] = np.array(flipped_e[e])
                    flipped = np.vstack([flipped,flipped_e[e]])
            # convert all points to parameter space (z, theta)
            converted = np.empty((0, 2))
            converted_e = {}
            axis = np.array([0,0,1])
            base_point = np.array([0.,0.,0.]) #np.array(face.normal[0].coordinates)
            radius = face.radius
            for e in flipped_e:
                converted_e[e] = []
                for p in flipped_e[e]:
                    converted_e[e].append(self.cylindrical_coordinates(p, base_point, radius, axis))
                converted_e[e] = np.array(converted_e[e])
                converted_e[e] = self.correct_angular_discontinuity(converted_e[e])
                converted_e[e][:,1:] = converted_e[e][:,1:]*10.
                converted = np.vstack([converted,converted_e[e]])
            # generate points inside face with bridson sampling
            max_min = [np.min(converted[:, 0]), np.min(converted[:, 1]), np.max(converted[:, 0]), np.max(converted[:, 1])]
            cell_size = max_min[2] - max_min[0]
            for e in face.edges:
                for l in face.edges[e].lines:
                    if face.edges[e].lines[l].length < cell_size:
                        cell_size = face.edges[e].lines[l].length
            cell_size = cell_size*0.1
            outside_point = np.min(converted, axis=0) - 23
            outside_point[0] -= 0.666*np.pi
            internal = self.bridson_sampling_2D(max_min, converted_e, converted, cell_size, outside_point)
            # genrate triangles from face points and remove those that are outside edges
            tri = Delaunay(internal[:,:2], qhull_options='Qt QJ')
            nodes, elements = self.keep_internal_triangles_2D(converted_e, internal, tri.simplices, outside_point)
            # convert back from cylindrical coordinates
            for node in nodes:
                nodes[node][1] = nodes[node][1]*0.1
                nodes[node] = np.array(self.cylindrical_coordinates(nodes[node], base_point, radius, axis, True))
            # flip back from x-y plane
            for node in nodes:
                nodes[node] = Point3D(nodes[node][0], nodes[node][1], nodes[node][2])
                if need_to_flip:
                    nodes[node] = np.array(self.rotate_point_about_axis(nodes[node], a0, a1, -flip_angle))
                else:
                    nodes[node] = np.array(nodes[node].coordinates)
                # move all points back so the centroid is back where it was
                nodes[node] = nodes[node] + np.array(face.normal[0].coordinates)
            face.g_mesh['nodes'] = nodes
            face.g_mesh['elements'] = elements


        # -----------
        # mesh TOROIDAL faces
        # -------------------------
        elif face.type == 'toroidal':
            # move all points so centroid is at origin
            centered = np.empty((0, 3))
            centered_e = {}
            for e in face.edges:
                centered_e[e] = face.edges[e].points - np.array(face.normal[0].coordinates)
                centered = np.vstack([centered,centered_e[e]])
            # flip base to x-y plane
            need_to_flip = False
            if face.normal_v == Vector3D(0,0,1) or face.normal_v == Vector3D(0,0,-1):
                flipped = centered
                flipped_e = centered_e
            else:
                need_to_flip = True
                flip_angle = self.angle_between(np.array(face.normal_v.coordinates),np.array([0.,0.,1.])) - np.pi
                flip_axis = np.cross(np.array(face.normal_v.coordinates),np.array([0.,0.,1.]))
                flip_axis = self.unit_vector(flip_axis)
                a0 = Point3D(0,0,0)
                a1 = Vector3D(flip_axis[0],flip_axis[1],flip_axis[2])
                flipped = np.empty((0, 3))
                flipped_e = {}
                for e in face.edges:
                    flipped_e[e] = []
                    for p in centered_e[e]:
                        p1 = Point3D(p[0],p[1],p[2])
                        flipped_e[e].append(self.rotate_point_about_axis(p1, a0, a1, flip_angle))
                    flipped_e[e] = np.array(flipped_e[e])
                    flipped = np.vstack([flipped,flipped_e[e]])
            # convert all points to parameter space (phi, theta)
            converted = np.empty((0, 2))
            converted_e = {}
            axis = np.array([0,0,1])
            base_point = np.array([0.,0.,0.]) #np.array(face.normal[0].coordinates)
            major_radius = face.major_radius
            minor_radius = face.minor_radius
            for e in flipped_e:
                converted_e[e] = []
                for p in flipped_e[e]:
                    converted_e[e].append(self.toroidal_coordinates(p, base_point, major_radius, minor_radius, axis))
                converted_e[e] = np.array(converted_e[e])
                converted_e[e] = self.correct_angular_discontinuity(converted_e[e])
                converted_e[e] = self.correct_angular_discontinuity(converted_e[e],0)
                converted = np.vstack([converted,converted_e[e]])
            # generate points inside face with bridson sampling
            max_min = [np.min(converted[:, 0]), np.min(converted[:, 1]), np.max(converted[:, 0]), np.max(converted[:, 1])]
            cell_size = max_min[2] - max_min[0]
            for e in face.edges:
                for l in face.edges[e].lines:
                    if face.edges[e].lines[l].length < cell_size:
                        cell_size = face.edges[e].lines[l].length
            cell_size = cell_size*0.1
            outside_point = np.min(converted, axis=0) - 23
            outside_point[0] -= 0.666*np.pi
            internal = self.bridson_sampling_2D(max_min, converted_e, converted, cell_size, outside_point)
            # genrate triangles from face points and remove those that are outside edges
            tri = Delaunay(internal[:,:2], qhull_options='Qt QJ')
            nodes, elements = self.keep_internal_triangles_2D(converted_e, internal, tri.simplices, outside_point)
#            self.plot_2D_face_mesh(nodes,elements,max_min[0],max_min[1],max_min[2],max_min[3])            
            # convert back from toroidal coordinates
            for node in nodes:
                nodes[node] = np.array(self.toroidal_coordinates(nodes[node], base_point, major_radius, minor_radius, axis, True))
            # flip back from x-y plane
            for node in nodes:
                nodes[node] = Point3D(nodes[node][0], nodes[node][1], nodes[node][2])
                if need_to_flip:
                    nodes[node] = np.array(self.rotate_point_about_axis(nodes[node], a0, a1, -flip_angle))
                else:
                    nodes[node] = np.array(nodes[node].coordinates)
                # move all points back so the centroid is back where it was
                nodes[node] = nodes[node] + np.array(face.normal[0].coordinates)
            face.g_mesh['nodes'] = nodes
            face.g_mesh['elements'] = elements


        # -----------
        # mesh CONICAL faces
        # -------------------------
        elif face.type == 'conical':
            # move all points so centroid is at origin
            centered = np.empty((0, 3))
            centered_e = {}
            for e in face.edges:
                centered_e[e] = face.edges[e].points - np.array(face.normal[0].coordinates)
                centered = np.vstack([centered,centered_e[e]])
            # flip base to x-y plane
            need_to_flip = False
            if face.normal_v == Vector3D(0,0,1) or face.normal_v == Vector3D(0,0,-1):
                flipped = centered
                flipped_e = centered_e
            else:
                need_to_flip = True
                flip_angle = self.angle_between(np.array(face.normal_v.coordinates),np.array([0.,0.,1.])) - np.pi
                flip_axis = np.cross(np.array(face.normal_v.coordinates),np.array([0.,0.,1.]))
                flip_axis = self.unit_vector(flip_axis)
                a0 = Point3D(0,0,0)
                a1 = Vector3D(flip_axis[0],flip_axis[1],flip_axis[2])
                flipped = np.empty((0, 3))
                flipped_e = {}
                for e in face.edges:
                    flipped_e[e] = []
                    for p in centered_e[e]:
                        p1 = Point3D(p[0],p[1],p[2])
                        flipped_e[e].append(self.rotate_point_about_axis(p1, a0, a1, flip_angle))
                    flipped_e[e] = np.array(flipped_e[e])
                    flipped = np.vstack([flipped,flipped_e[e]])
            # convert all points to parameter space (r, theta)
            converted = np.empty((0, 2))
            converted_e = {}
            axis = np.array([0,0,1])
            base_point = np.array([0.,0.,0.])
            base_radius = face.base_radius
            semi_angle = face.semi_angle
            flip_z = True
            if face.normal_v == Vector3D(0,0,1):
                flip_z = False
            for e in flipped_e:
                converted_e[e] = []
                for p in flipped_e[e]:
                    converted_e[e].append(self.conical_coordinates(p, base_radius, semi_angle, flip_z))
                converted_e[e] = np.array(converted_e[e])
                converted_e[e] = self.correct_angular_discontinuity(converted_e[e])
                converted = np.vstack([converted,converted_e[e]])
            # generate points inside face with bridson sampling
            max_min = [np.min(converted[:, 0]), np.min(converted[:, 1]), np.max(converted[:, 0]), np.max(converted[:, 1])]
            cell_size = max_min[2] - max_min[0]
            for e in face.edges:
                for l in face.edges[e].lines:
                    if face.edges[e].lines[l].length < cell_size:
                        cell_size = face.edges[e].lines[l].length
            cell_size = cell_size*0.1
            outside_point = np.min(converted, axis=0) - 23
            outside_point[0] -= 0.666*np.pi
            internal = self.bridson_sampling_2D(max_min, converted_e, converted, cell_size, outside_point)
            # genrate triangles from face points and remove those that are outside edges
            tri = Delaunay(internal[:,:2], qhull_options='Qt QJ')
            nodes, elements = self.keep_internal_triangles_2D(converted_e, internal, tri.simplices, outside_point)
#            self.plot_2D_face_mesh(nodes,elements,max_min[0],max_min[1],max_min[2],max_min[3])            
            # convert back from conical coordinates
            for node in nodes:
                nodes[node] = np.array(self.conical_coordinates(nodes[node], base_radius, semi_angle, flip_z, True))
            # flip back from x-y plane
            for node in nodes:
                nodes[node] = Point3D(nodes[node][0], nodes[node][1], nodes[node][2])
                if need_to_flip:
                    nodes[node] = np.array(self.rotate_point_about_axis(nodes[node], a0, a1, -flip_angle))
                else:
                    nodes[node] = np.array(nodes[node].coordinates)
                # move all points back so the centroid is back where it was
                nodes[node] = nodes[node] + np.array(face.normal[0].coordinates)
            face.g_mesh['nodes'] = nodes
            face.g_mesh['elements'] = elements


        # -----------
        # mesh REVOLUTION faces
        # -------------------------
        elif face.type == 'revolution':
            pass


        # -----------
        # mesh EXTRUSION faces
        # -------------------------
        elif face.type == 'extrusion':
            pass


        # -----------
        # mesh B_SPLINE faces
        # -------------------------
        elif face.type == 'b_spline_surface':
            pass





    def meshFace(self,face,seeded=False,element_size=5.):
        '''
    Creates a surface mesh for a given Face based on
    either the seeded edges or all the points on the
    edge lines of the face. It first converts the
    edge point coordinates to a 2D parameter space,
    and then uses bridson sampling and delaunay to
    generate the surface triangles.
    '''
        if face.type in ['plane', 'cylindrical', 'toroidal', 'conical']:
            
            # add extra points on edge of cylinder so there are an
            # equal amount of points on that as there is on the arcs
            if face.type == 'cylindrical':
                if not seeded:
                    smallest_length = 10e5
                    for e in face.edges:
                        for p in range(len(face.edges[e].points)-1):
                            if np.linalg.norm(face.edges[e].points[p]-face.edges[e].points[p+1]) < smallest_length:
                                smallest_length = np.linalg.norm(face.edges[e].points[p]-face.edges[e].points[p+1])
                    new_points_e = {}
                    new_points = np.empty((0, 3))
                    for e in face.edges:
                        new_points_e[e] = []
                        for p in range(len(face.edges[e].points)):
                            if abs(np.linalg.norm(face.edges[e].points[p]-face.edges[e].points[p-1]) - smallest_length) < 0.1:
                                new_points_e[e].append(face.edges[e].points[p])
                            else:
                                inserts = int(np.linalg.norm(face.edges[e].points[p]-face.edges[e].points[p-1])/smallest_length) + 1 
                                points = face.edges[e].points
                                for i in range(inserts-1):
                                    new_points_e[e].append([points[p-1][0]+(i+1)*(points[p][0]-points[p-1][0])/inserts,
                                                            points[p-1][1]+(i+1)*(points[p][1]-points[p-1][1])/inserts,
                                                            points[p-1][2]+(i+1)*(points[p][2]-points[p-1][2])/inserts])
                                new_points_e[e].append(points[p])
                        new_points_e[e] = np.array(new_points_e[e])
                        new_points = np.vstack([new_points,new_points_e[e]])

            # move all points so centroid is at origin
            centered = np.empty((0, 3))
            centered_e = {}
            for e in face.edges:
                if face.type == 'plane':
                    if not seeded:
                        centered_e[e] = face.edges[e].points - face.centroid_np
                    else:
                        centered_e[e] = face.edges[e].seeds - face.centroid_np
                elif face.type in ['cylindrical', 'toroidal', 'conical']:
                    if not seeded:
                        centered_e[e] = face.edges[e].points - np.array(face.normal[0].coordinates)
                        if face.type == 'cylindrical':
                            centered_e[e] = new_points_e[e] - np.array(face.normal[0].coordinates)
                    else:
                        centered_e[e] = face.edges[e].seeds - np.array(face.normal[0].coordinates)
                centered = np.vstack([centered,centered_e[e]])

            # flip to x-y plane
            need_to_flip = False
            if face.normal_v == Vector3D(0,0,1) or face.normal_v == Vector3D(0,0,-1):
                flipped = centered
                flipped_e = centered_e
            else:
                need_to_flip = True
                flip_angle = self.angle_between(np.array(face.normal_v.coordinates),np.array([0.,0.,1.]))
                if face.type in ['cylindrical', 'toroidal', 'conical']:
                    flip_angle -= np.pi
                flip_axis = np.cross(np.array(face.normal_v.coordinates),np.array([0.,0.,1.]))
                flip_axis = self.unit_vector(flip_axis)
                a0 = Point3D(0,0,0)
                a1 = Vector3D(flip_axis[0],flip_axis[1],flip_axis[2])
                flipped = np.empty((0, 3))
                flipped_e = {}
                for e in face.edges:
                    flipped_e[e] = []
                    for p in centered_e[e]:
                        p1 = Point3D(p[0],p[1],p[2])
                        flipped_e[e].append(self.rotate_point_about_axis(p1, a0, a1, flip_angle))
                    flipped_e[e] = np.array(flipped_e[e])
                    flipped = np.vstack([flipped,flipped_e[e]])

            # convert all points to parameter space
            converted = np.empty((0, 2))
            converted_e = {}
            axis = np.array([0,0,1])
            base_point = np.array([0.,0.,0.])
            if face.type == 'cylindrical':
                radius = face.radius
            elif face.type == 'toroidal':
                major_radius = face.major_radius
                minor_radius = face.minor_radius
            elif face.type == 'plane':
                converted = flipped[:,:2]
            elif face.type == 'conical':
                base_radius = face.base_radius
                semi_angle = face.semi_angle
                flip_z = True
                if face.normal_v == Vector3D(0,0,1):
                    flip_z = False
            for e in flipped_e:
                converted_e[e] = []
                if face.type == 'plane':
                    converted_e[e] = flipped_e[e][:,:2]
                for p in flipped_e[e]:
                    if face.type == 'cylindrical':
                        converted_e[e].append(self.cylindrical_coordinates(p, base_point, radius, axis))
                    elif face.type == 'toroidal':
                        converted_e[e].append(self.toroidal_coordinates(p, base_point, major_radius, minor_radius, axis))
                    elif face.type == 'conical':
                        converted_e[e].append(self.conical_coordinates(p, base_radius, semi_angle, flip_z))
                if face.type in ['cylindrical', 'toroidal', 'conical']:
                    converted_e[e] = np.array(converted_e[e])
                    converted_e[e] = self.correct_angular_discontinuity(converted_e[e])
                if face.type == 'cylindrical':
                    converted_e[e][:,1:] = converted_e[e][:,1:]*10.
                if face.type == 'toroidal':
                    converted_e[e] = self.correct_angular_discontinuity(converted_e[e],0)
                converted = np.vstack([converted,converted_e[e]])
                
            # generate points inside face with bridson sampling
            max_min = [np.min(converted[:, 0]), np.min(converted[:, 1]), np.max(converted[:, 0]), np.max(converted[:, 1])]
            cell_size = max_min[2] - max_min[0]
            for e in face.edges:
                for l in face.edges[e].lines:
                    if face.edges[e].lines[l].length < cell_size:
                        cell_size = face.edges[e].lines[l].length
            if not seeded:
                if face.type == 'plane':
                    cell_size = cell_size*0.7
                elif face.type in ['toroidal', 'cylindrical', 'conical']:
                    cell_size = cell_size*0.1
            else:
                cell_size = element_size
            outside_point = np.min(converted, axis=0) - 23
            outside_point[0] -= 0.666*np.pi
            internal = self.bridson_sampling_2D(max_min, converted_e, converted, cell_size, outside_point)
#            if seeded:
#                fig = plt.gcf()
#                fig.set_size_inches(8., 8.)
#                plt.axis([max_min[0]-1,max_min[2]+1,max_min[1]-1,max_min[3]+1])
#                for p in internal:
#                    plt.plot(p[0], p[1], 'o', color='green')
#                for p in converted:
#                    plt.plot(p[0], p[1], 'o', color='blue')
#                plt.show() 

            # genrate triangles from face points and remove those that are outside edges
            tri = Delaunay(internal[:,:2], qhull_options='Qt QJ')
            nodes, elements = self.keep_internal_triangles_2D(converted_e, internal, tri.simplices, outside_point)
            
#            if seeded:
#                self.plot_2D_face_mesh(nodes,elements,max_min[0],max_min[1],max_min[2],max_min[3])       

            # convert back from parameter space
            for node in nodes:
                if face.type == 'cylindrical':
                    nodes[node][1] = nodes[node][1]*0.1
                    nodes[node] = np.array(self.cylindrical_coordinates(nodes[node], base_point, radius, axis, True))
                elif face.type == 'toroidal':
                    nodes[node] = np.array(self.toroidal_coordinates(nodes[node], base_point, major_radius, minor_radius, axis, True))
                elif face.type == 'conical':
                    nodes[node] = np.array(self.conical_coordinates(nodes[node], base_radius, semi_angle, flip_z, True))
            # flip back from x-y plane
            for node in nodes:
                if face.type == 'plane':
                    nodes[node] = Point3D(nodes[node][0], nodes[node][1], 0.)
                elif face.type in ['cylindrical', 'toroidal', 'conical']:
                    nodes[node] = Point3D(nodes[node][0], nodes[node][1], nodes[node][2])
                if need_to_flip:
                    nodes[node] = np.array(self.rotate_point_about_axis(nodes[node], a0, a1, -flip_angle))
                else:
                    nodes[node] = np.array(nodes[node].coordinates)
                # move all points back to original position
                if face.type == 'plane':
                    nodes[node] = nodes[node] + face.centroid_np
                elif face.type in ['cylindrical', 'toroidal', 'conical']:
                    nodes[node] = nodes[node] + np.array(face.normal[0].coordinates)    
            # reverse face direction if plane face with face normal [0,0,-1]
            if face.type == 'plane':
                if face.normal_v == Vector3D(0,0,-1):
                    for elm in elements:
                        elements[elm] = elements[elm][::-1]
            if not seeded:
                face.g_mesh['nodes'] = nodes
                face.g_mesh['elements'] = elements
            else:
                face.mesh['nodes'] = nodes
                face.mesh['elements'] = elements



        
    
    
    def meshSolid(self,part,element_size):
        '''
    Creates a tetrahedral mesh for a given solid part
    using the seeded edges of that part. It first uses
    the meshFace function to mesh all the faces, and
    then fill in the volume.
    '''
        nodes = {}
        elements = {}

        # first mesh any faces that have not already been meshed
        for f in part.faces:
            if len(part.faces[f].mesh['elements']) == 0:
                self.meshFace(part.faces[f],True,element_size)
        
        # next use the face meshes to generate nodes inside the volume
#        internal = self.bridson_sampling_3D(max_min, converted_e, converted, cell_size, outside_point)
        
        # next use face meshes and nodes inside volume to generate tetrahedra
        
        
        # next remove tetrahedra that are located outside the volume
        
        
        # next merge all nodes that are in the same location
        
        
        # finally ensure any nodes that are located on a face are attached to that face

        
        return nodes, elements
    
    
    




    # ------------
    # various helper functions 
    # --------------------------

    def interpolateCurve(self, points, num_points):
        '''
     Creates a given number of points evenly spaced out among the
     points that are already there.
     '''
        # Convert the list of points to a NumPy array
        points = np.array(points)
        # Calculate the cumulative distance along the curve
        distances = np.cumsum(np.sqrt(np.sum(np.diff(points, axis=0)**2, axis=1)))
        distances = np.insert(distances, 0, 0)  # Add a 0 at the beginning for the start
        # Normalize the distances to the range [0, 1]
        distances /= distances[-1]
        # Interpolation functions for each dimension
        interpolate_x = interp1d(distances, points[:, 0], kind='linear')
        interpolate_y = interp1d(distances, points[:, 1], kind='linear')
        interpolate_z = interp1d(distances, points[:, 2], kind='linear')
        # Generate the evenly spaced points along the curve
        new_distances = np.linspace(0, 1, num_points)
        new_points = np.zeros((num_points, 3))
        new_points[:, 0] = interpolate_x(new_distances)
        new_points[:, 1] = interpolate_y(new_distances)
        new_points[:, 2] = interpolate_z(new_distances)
        return new_points


    def unit_vector(self, vector):
        '''
    Returns the unit vector of a vector.
    '''
        return vector / np.linalg.norm(vector)
    
    
    def angle_between(self, v1, v2):
        '''
    Returns the angle between two vectors.
    '''
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


    def rotate_point_about_axis(self,point,axis_p0,axis_p1,angle):
        '''
    First moves point, axis_p0 and axis_p1 together so that 
    axis0 is at the origin. Then creates the two quaternions 
    and rotates the point about the arbitrary axis defined
    by axis_p0 and axis_p1. Finally moves the point back so that 
    axis_p0 is at its original position again. 
    '''
        axis1_0 = axis_p1 - axis_p0
        point_0 = point - axis_p0
        qPnt0 = Quaternion()
        qPnt0.vectorToQuat(point_0)
        qRot = Quaternion()
        qRot.axisAngleToQuat(axis1_0.normalized(),angle)
        qPnt1 = qRot.multi(qPnt0.multi(qRot.conj()))
        return [qPnt1.i+axis_p0.x(), qPnt1.j+axis_p0.y(), qPnt1.k+axis_p0.z()]


    # Assume cylinder_axis is a unit vector along the cylinder's axis
    # base_point is a point on the cylinder's axis
    # point is the Cartesian coordinate of the edge point
    def cylindrical_coordinates(self, point, base_point, radius, cylinder_axis, reverse=False):
        '''
     Convert Cartesian coordinates to 2D cylindrical coordinates [z, theta], 
     or back to cartesian coordinates from toroidal coordinates [x,y,z].
     '''
        if reverse:
            z = point[0]
            theta = point[1]
            # Calculate x and y using radius and theta
            x = radius * np.cos(theta-0.5*np.pi)
            y = radius * np.sin(theta-0.5*np.pi)
            # Transform this point in the plane (x, y) back into 3D space
            # Find two orthogonal vectors in the plane perpendicular to cylinder_axis
            v1 = np.cross(cylinder_axis, [1, 0, 0])
            if np.linalg.norm(v1) == 0:  # If cylinder_axis is parallel to [1, 0, 0], choose another vector
                v1 = np.cross(cylinder_axis, [0, 1, 0])
            v1 = v1 / np.linalg.norm(v1)
            v2 = np.cross(cylinder_axis, v1)
            # Combine x, y with the base point and the height along the cylinder axis
            new_point = base_point + z * cylinder_axis + x * v1 + y * v2
            return new_point
        else:
            # Project point onto the plane perpendicular to the cylinder's axis
            projected_point = point - np.dot(point - base_point, cylinder_axis) * cylinder_axis
            # Calculate angle
            theta = np.arctan2(projected_point[1], projected_point[0])  # Assuming the base point is at the origin
            # Calculate height
            z = np.dot(point - base_point, cylinder_axis)
            return [z,theta]


    def correct_angular_discontinuity(self,points,xy=1):
        '''
     Corrects angle from -pi to plus pi or the other way around
     depending on neighbouring points.
     '''
        corrected_points = points.copy()
        for i in range(1, len(points)):
            delta_theta = corrected_points[i][xy] - corrected_points[i - 1][xy]
            if delta_theta > np.pi:
                corrected_points[i][xy] -= 2 * np.pi
            elif delta_theta < -np.pi:
                corrected_points[i][xy] += 2 * np.pi
        return corrected_points


    def toroidal_coordinates(self, point, base_point, R, r, torus_normal, reverse=False):
        '''
     Convert Cartesian coordinates to 2D toroidal coordinates [phi, theta], 
     or back to cartesian coordinates from toroidal coordinates [x,y,z].
     '''
        if reverse:
            phi = point[0]
            theta = point[1]
            # Calculate the position on the larger circle (centerline of the torus tube)
            x_major = (R + r * np.cos(phi)) * np.cos(theta-0.5*np.pi)
            y_major = (R + r * np.cos(phi)) * np.sin(theta-0.5*np.pi)
            z_major = r * np.sin(phi)
            # Transform this point in 3D space
            # Note: torus_normal is a unit vector normal to the plane of the major circle of the torus
            # Find two orthogonal vectors in the plane of the major circle
            v1 = np.cross(torus_normal, [1, 0, 0])
            if np.linalg.norm(v1) == 0:
                v1 = np.cross(torus_normal, [0, 1, 0])
            v1 = v1 / np.linalg.norm(v1)
            v2 = np.cross(torus_normal, v1)
            # Combine the coordinates with the base point
            new_point = base_point + x_major * v1 + y_major * v2 + z_major * torus_normal
            return new_point
        else:
            # Convert from Cartesian to toroidal coordinates
            # Adjust point relative to base_point
            adjusted_point = point - base_point
            # Calculate theta (angle in the major circle plane)
            theta = np.arctan2(adjusted_point[1], adjusted_point[0])
            # Project point onto the plane of the major circle to find the minor radius
            projected_point = adjusted_point - np.dot(adjusted_point, torus_normal) * torus_normal
            minor_radius_component = np.linalg.norm(projected_point) - R
            # Calculate the vector pointing towards the point in the minor circle's plane
            minor_circle_vector = adjusted_point - projected_point
            # Compute phi using the vector in the plane of the minor circle
            # Assuming that the torus is centered at the origin and the major circle is in the xy-plane
            phi = np.arctan2(minor_circle_vector[2], minor_radius_component)
            return [phi, theta]


    def conical_coordinates(self, point, base_radius, semi_angle, flip_z, reverse=False):
        '''
     Convert Cartesian coordinates to 2D conical coordinates [r, theta] (with semi_angle to
     find the height), or back to cartesian coordinates from conical coordinates [x,y,z].
     '''
        if reverse:
            r = point[0]
            theta = point[1]
            # Calculate the Cartesian coordinates
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            z = (r - base_radius) * np.tan(0.5*np.pi - semi_angle)
            if flip_z:
                z = -z
            return [x, y, z]
        else:
            # Calculate the radial distance r in the xy-plane
            r = np.sqrt(point[0]**2 + point[1]**2)
            # Calculate the polar angle theta
            theta = np.arctan2(point[1], point[0])
            # Adjust theta to be in the range [0, 2pi]
            if theta < 0:
                theta += 2 * np.pi
            return [r, theta]


    def bridson_sampling_2D(self, max_min, edges, points, radius, outside_point, k=30):
        '''
     Generates points evenly spaced out inside the area given by a set of edges.
     '''        
        start_time = time.time()
#        print('\n\nStarting bridson_sampling_2D() ...')
        width = max_min[2]-max_min[0]
        height = max_min[3]-max_min[1]
        # Step 0: Initialize the necessary variables
        cell_size = radius/np.sqrt(3)
        grid_width = int(np.ceil(width / cell_size))
        grid_height = int(np.ceil(height / cell_size))
        grid = np.empty((grid_width, grid_height), dtype=object)
        process_list = []
        sample_points = []
        internal_points = []
        for p in points:
            sample_points.append((p[0],p[1]))
            grid[self.grid_index_2D(p, cell_size, max_min[0], max_min[1])] = (p[0], p[1])
    
        # Step 1: Randomly select the initial point and add it to the grid and process list
        initial_point = random.choice(sample_points)
        process_list.append(initial_point)
        grid[self.grid_index_2D(initial_point, cell_size, max_min[0], max_min[1])] = initial_point
    
        # Step 2: While the process list is not empty, continue sampling
        while process_list:
            point = random.choice(process_list)  # Select a random point from process list
            process_list.remove(point)
            for _ in range(k):  # Try to find k valid points in the annulus
                new_point = self.get_point_in_annulus_2D(point, radius)
                if max_min[0] <= new_point[0] < max_min[2] and max_min[1] <= new_point[1] < max_min[3]:  # Check if the point is within the domain
                    if self.is_valid_point_2D(new_point, grid, cell_size, radius, max_min[0], max_min[1]):
                        if self.is_within_boundary_2D(new_point, edges, outside_point):
                            process_list.append(new_point)
                            internal_points.append(new_point)
                            sample_points.append(new_point)
                            grid[self.grid_index_2D(new_point, cell_size, max_min[0], max_min[1])] = new_point
        sample_points = np.array(sample_points)
    
        end_time = time.time()
#        print('... finished in {:.3f} seconds'.format(end_time-start_time))
    
        return sample_points
    
    
    def grid_index_2D(self, point, cell_size, x_min, y_min):
        '''
     Helper function for bridson_sampling_2D.
     '''
        return int((point[0] - x_min) / cell_size), int((point[1] - y_min) / cell_size)
    
    
    def is_valid_point_2D(self, new_point, grid, cell_size, radius, x_min, y_min):
        '''
     Helper function for bridson_sampling_2D.
     '''
        x, y = self.grid_index_2D(new_point, cell_size, x_min, y_min)
        for i in range(max(0, x - 2), min(x + 3, grid.shape[0])):
            for j in range(max(0, y - 2), min(y + 3, grid.shape[1])):
                point = grid[i, j]
                if point is not None:
                    distance = np.linalg.norm(np.array(point) - np.array(new_point))
                    if distance < radius:
                        return False
        return True
    
    
    def get_point_in_annulus_2D(self, center_point, radius):
        '''
     Helper function for bridson_sampling_2D.
     '''
        random_angle_phi = 2 * np.pi * random.random()
        random_radius = random.uniform(radius, 2 * radius)
        random_x = center_point[0] + random_radius * np.cos(random_angle_phi)
        random_y = center_point[1] + random_radius * np.sin(random_angle_phi)
        return random_x, random_y
    
    
    def is_on_boundary_2D(self, point, edges):
        '''
     Helper function for bridson_sampling_2D.
     '''
        for e in edges:
            for p in edges[e]:
                if np.linalg.norm(point - p[:2]) < 0.2:
                    return True
        
        
    def is_within_boundary_2D(self, point, edges, outside_point):
        '''
     Helper function for bridson_sampling_2D.
     '''
        intersections = 0
        for e in edges:
            for p in range(len(edges[e])):
                line_segment = np.array((edges[e][p-1],edges[e][p]))
                if self.does_intersect_2D(point, outside_point, line_segment):
                    intersections += 1
        return intersections % 2 == 1

       
    def does_intersect_2D(self, p1, p2, line_segment):
        '''
     Helper function for bridson_sampling_2D.
     '''
        if self.intersect_2D(p1, p2, line_segment[0], line_segment[1]):
            return True
        return False


    def intersect_2D(self, p1, p2, p3, p4):
        '''
     Helper function for bridson_sampling_2D.
     '''
        d = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
        if d == 0:
            return False
        u = ((p3[0] - p1[0]) * (p4[1] - p3[1]) - (p3[1] - p1[1]) * (p4[0] - p3[0])) / d
        v = ((p3[0] - p1[0]) * (p2[1] - p1[1]) - (p3[1] - p1[1]) * (p2[0] - p1[0])) / d
        if u >= 0 and u <= 1 and v >= 0 and v <= 1:
            return True
        return False
    

    def is_flat_triangle_2D(self, p1, p2, p3):
        '''
     Helper function for keep_internal_triangles_2D.
     '''
        AB = p2 - p1
        AC = p3 - p1
        normal = np.cross(AB, AC)
        if np.linalg.norm(normal) < 0.01:
            return True
        return False
    
    
    def keep_internal_triangles_2D(self, edges, points, triangles, outside_point):
        '''
     Remove all triangles that are outside the edges boundary.
     '''
        internal_triangles = []
        for triangle in triangles:
            keep_triangle = True
            p1 = points[triangle[0]]
            p2 = points[triangle[1]]
            p3 = points[triangle[2]]
            triangle_center_point = np.array([np.mean(np.array([p1[0],p2[0],p3[0]])),
                                              np.mean(np.array([p1[1],p2[1],p3[1]]))])
            if not self.is_within_boundary_2D(triangle_center_point, edges, outside_point):
#                if not self.is_on_boundary_2D(triangle_center_point, edges):
                keep_triangle = False
            if self.is_flat_triangle_2D(p1,p2,p3):
                keep_triangle = False
            if keep_triangle:
                internal_triangles.append(triangle)
        nodes = {}
        elements = {}
        element_number = 1
        for triangle in internal_triangles:
            if triangle[0] not in nodes:
                nodes[triangle[0]] = points[triangle[0]]
            if triangle[1] not in nodes:
                nodes[triangle[1]] = points[triangle[1]]
            if triangle[2] not in nodes:
                nodes[triangle[2]] = points[triangle[2]]
            elements[element_number] = [triangle[0], triangle[1], triangle[2]]
            element_number += 1
        return nodes, elements


    def plot_2D_face_mesh(self, nodes,elements,x_min,y_min,x_max,y_max):
        # for debugging purposes only
        fig = plt.gcf()
        fig.set_size_inches(8., 8.)
        plt.axis([x_min-1,x_max+1,y_min-1,y_max+1])
        for element in elements:
            plt.plot([nodes[elements[element][0]][0], nodes[elements[element][1]][0]],
                     [nodes[elements[element][0]][1], nodes[elements[element][1]][1]],color='blue')
            plt.plot([nodes[elements[element][0]][0], nodes[elements[element][2]][0]],
                     [nodes[elements[element][0]][1], nodes[elements[element][2]][1]],color='blue')
            plt.plot([nodes[elements[element][2]][0], nodes[elements[element][1]][0]],
                     [nodes[elements[element][2]][1], nodes[elements[element][1]][1]],color='blue')
        for node in nodes:
            plt.plot(nodes[node][0], nodes[node][1], 'o', color='green')
        plt.show()   