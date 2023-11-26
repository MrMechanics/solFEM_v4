# -*- coding: utf-8 -*-
"""
Created on Mon May  8 08:27:48 2023

@author: MATMATH
"""

import sys
sys.path.insert(1, '../Objects')
sys.path.insert(1, '../Modules')

import numpy as np
from scipy.spatial import Delaunay, KDTree
import matplotlib.pyplot as plt
#from skimage.measure import marching_cubes
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import random
import matplotlib.tri as mtri
from quaternion import *
from timeit import time





def is_on_boundary(point, all_edge_seeds):
    for s in all_edge_seeds:
        if np.linalg.norm(point - s) < 0.2:
            return True
        
        
        
def is_within_boundary(point, edges2D, lines2D, all_edge_seeds):
    outside_point = np.min(all_edge_seeds, axis=0) - 1
    outside_point[0] -= 0.666*np.pi
    intersections = 0
    for edge in edges2D:
        for line in edges2D[edge]:
            for p in range(len(lines2D[line])-1):
                line_segment = np.vstack((lines2D[line][p],lines2D[line][p+1]))
                if does_intersect(point, outside_point, line_segment):
                    intersections += 1
    return intersections % 2 == 1


   
def does_intersect(p1, p2, line_segment):
    if intersect(p1, p2, line_segment[0], line_segment[1]):
        return True
    return False



def intersect(p1, p2, p3, p4):
    d = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
    if d == 0:
        return False
    u = ((p3[0] - p1[0]) * (p4[1] - p3[1]) - (p3[1] - p1[1]) * (p4[0] - p3[0])) / d
    v = ((p3[0] - p1[0]) * (p2[1] - p1[1]) - (p3[1] - p1[1]) * (p2[0] - p1[0])) / d
    if u >= 0 and u <= 1 and v >= 0 and v <= 1:
        return True
    return False



def get_rad(index, samples, r_max):
    r1 = 1.1*np.linalg.norm(samples[index-1] - samples[index])
    r2 = 1.1*np.linalg.norm(samples[index] - samples[index+1])
    return min(r_max, min(r1,r2))



def get_new_points(point, radius, samples):
    new_points = []
    for t in range(14):
        theta = (t/14)*2*np.pi
        r = radius*1.1
        x = point[0] + r*np.cos(theta)
        y = point[1] + r*np.sin(theta)
        new_points.append([x,y])
    return new_points



def keep_point(point, radius, samples, edges2D, lines2D, all_edge_seeds, outside=False):
    if not outside:
        if not is_within_boundary(point, edges2D, lines2D, all_edge_seeds):
            return False
    else:
        if is_within_boundary(point, edges2D, lines2D, all_edge_seeds):
            return False
    # Check if the point intersects any other points
    for s in samples:
        if np.linalg.norm(point - s) < 0.75*radius:
            return False
    return True



def is_valid(point, radius, seeds2D):
    # Check if the point intersects any other points
    for s in seeds2D:
        if np.linalg.norm(point - s) < radius:
            return False
    return True



# bridson_sampling(max_min, face_nodes_kdtree, all_face_nodes, all_face_elements, max_element_size)
# max_min = [x_min, y_min, x_max, y_max]
def bridson_sampling_2D(max_min, edges2D, lines2D, all_edge_seeds, radius, k=30):
    start_time = time.time()
    print('\n\nStarting bridson_sampling_2D() ...')
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
    for node in all_edge_seeds:
        sample_points.append((node[0],node[1]))
        grid[grid_index_2D(node, cell_size, max_min[0], max_min[1])] = (node[0], node[1])

    # Step 1: Randomly select the initial point and add it to the grid and process list
#    initial_point = (random.uniform(0, width), random.uniform(0, height), random.uniform(0, depth))
    initial_point = random.choice(sample_points)
    process_list.append(initial_point)
#    sample_points.append(initial_point)
    grid[grid_index_2D(initial_point, cell_size, max_min[0], max_min[1])] = initial_point

    # Step 2: While the process list is not empty, continue sampling
    while process_list:
        point = random.choice(process_list)  # Select a random point from process list
        process_list.remove(point)
        for _ in range(k):  # Try to find k valid points in the annulus
            new_point = get_point_in_annulus_2D(point, radius)
            if max_min[0] <= new_point[0] < max_min[2] and max_min[1] <= new_point[1] < max_min[3]:  # Check if the point is within the domain
                if is_valid_point_2D(new_point, grid, cell_size, radius, max_min[0], max_min[1]):
                    if is_within_boundary(new_point, edges2D, lines2D, all_edge_seeds):
                        process_list.append(new_point)
                        internal_points.append(new_point)
                        sample_points.append(new_point)
                        grid[grid_index_2D(new_point, cell_size, max_min[0], max_min[1])] = new_point
    sample_points = np.array(sample_points)

    end_time = time.time()
    print('all_edge_seeds:', len(all_edge_seeds))
    print('sample_points:', len(sample_points))
    print('... finished in {:.3f} seconds'.format(end_time-start_time))

    return sample_points

def grid_index_2D(point, cell_size):
    return int(point[0] / cell_size), int(point[1] / cell_size)

def grid_index_2D(point, cell_size, x_min, y_min):
    return int((point[0] - x_min) / cell_size), int((point[1] - y_min) / cell_size)

def is_valid_point_2D(new_point, grid, cell_size, radius, x_min, y_min):
    x, y = grid_index_2D(new_point, cell_size, x_min, y_min)
    for i in range(max(0, x - 2), min(x + 3, grid.shape[0])):
        for j in range(max(0, y - 2), min(y + 3, grid.shape[1])):
            point = grid[i, j]
            if point is not None:
                distance = np.linalg.norm(np.array(point) - np.array(new_point))
                if distance < radius:
                    return False
    return True

def get_point_in_annulus_2D(center_point, radius):
    random_angle_phi = 2 * np.pi * random.random()
    random_radius = random.uniform(radius, 2 * radius)
    random_x = center_point[0] + random_radius * np.cos(random_angle_phi)
    random_y = center_point[1] + random_radius * np.sin(random_angle_phi)
    return random_x, random_y








def poisson_disk_sampling(edges2D, lines2D, all_edge_seeds, r_max):
    samples = all_edge_seeds
    boundary_samples = samples
    previous_samples = []
    
    # first layer of points inside boundary
    for p in range(len(samples)-1):
        current_point = samples[p]
        current_rad = get_rad(p, samples, r_max)
        new_points = get_new_points(current_point, current_rad, samples)
        for n in new_points:
            if keep_point(n, current_rad, samples, edges2D, lines2D, all_edge_seeds):
                previous_samples.append(n)
                samples = np.vstack((samples,np.array(n)))
    previous_samples = np.array(previous_samples)

    # add more layers of points inside until no more space    
    for n in range(10):
        new_samples = []
        for p in range(len(previous_samples)-1):
            current_point = previous_samples[p]
            current_rad = get_rad(p, previous_samples, r_max)
            new_points = get_new_points(current_point, current_rad, previous_samples)
            for n in new_points:
                if keep_point(n, current_rad, samples, edges2D, lines2D, all_edge_seeds):
                    new_samples.append(n)
                    samples = np.vstack((samples,np.array(n)))
        if len(new_samples) == 0:
            break
        previous_samples = np.array(new_samples)

    # add a layer of points outside boundary
    for p in range(len(boundary_samples)-1):
        current_point = boundary_samples[p]
        current_rad = get_rad(p, boundary_samples, r_max)
        new_points = get_new_points(current_point, current_rad, boundary_samples)
        for n in new_points:
            if keep_point(n, current_rad, samples, edges2D, lines2D, all_edge_seeds, outside=True):
                samples = np.vstack((samples,np.array(n)))
#    fig = plt.gcf()
#    fig.set_size_inches(8., 8.)
#    plt.plot(samples[:,0], samples[:,1], '.', color='green')
#    plt.plot(all_edge_seeds[:,0], all_edge_seeds[:,1], 'o', color='blue')
#    plt.show()

    return samples



def is_flat_triangle(p1,p2,p3):
    AB = p2 - p1
    AC = p3 - p1
    normal = np.cross(AB, AC)
    if np.linalg.norm(normal) == 0.:
        return True
    return False


def keep_internal_triangles(edges2D, lines2D, all_points, all_edge_seeds, triangles):
    # remove all triangles that are outside the boundary
    internal_triangles = []
#    print('triangles', triangles)
    for triangle in triangles:
        keep_triangle = True
#        print('triangle', triangle)
        p1 = all_points[triangle[0]]
        p2 = all_points[triangle[1]]
        p3 = all_points[triangle[2]]
#        print('p1, p2, p3:\n', p1, p2, p3)
        triangle_center_point = np.array([np.mean(np.array([p1[0],p2[0],p3[0]])),
                                          np.mean(np.array([p1[1],p2[1],p3[1]]))])
#        print('triangle_center_point:', triangle_center_point)
        if not is_within_boundary(triangle_center_point, edges2D, lines2D, all_edge_seeds):
#            print('triangle is outside boundary!', triangle_center_point)
            if not is_on_boundary(triangle_center_point, all_edge_seeds):
#                print('triangle is not on boundary!', triangle_center_point)
                keep_triangle = False
        if is_flat_triangle(p1,p2,p3):
            keep_triangle = False
        if keep_triangle:
            internal_triangles.append(triangle)
    # remove all points that are outside the boundary
    nodes = {}
    elements = {}
    element_number = 1
    for triangle in internal_triangles:
        if triangle[0] not in nodes:
            nodes[triangle[0]] = {'coord': all_points[triangle[0]], 'number': triangle[0]}
        if triangle[1] not in nodes:
            nodes[triangle[1]] = {'coord': all_points[triangle[1]], 'number': triangle[1]}
        if triangle[2] not in nodes:
            nodes[triangle[2]] = {'coord': all_points[triangle[2]], 'number': triangle[2]}
        elements[element_number] = [nodes[triangle[0]]['number'], nodes[triangle[1]]['number'], nodes[triangle[2]]['number'],]
#        elements[element_number] = triangle.tolist()
        element_number += 1

    return nodes, elements, internal_triangles



def create_2D_face_mesh(edges2D, lines2D, all_edge_seeds, max_element_size):

    # generate points to fill the inside of the edges
    # and one layer of points on the outside of all edges
#    poisson_sample_points = poisson_disk_sampling(edges2D, lines2D, all_edge_seeds, max_element_size)
    print('all_edge_seeds:', all_edge_seeds)
    x_min = min(all_edge_seeds[:,0])
    y_min = min(all_edge_seeds[:,1])
    x_max = max(all_edge_seeds[:,0])
    y_max = max(all_edge_seeds[:,1])
    max_min = [x_min,y_min,x_max,y_max]
    poisson_sample_points = bridson_sampling_2D(max_min, edges2D, lines2D, all_edge_seeds, max_element_size, k=30)

    # group sample points and edge points for meshing
    all_points = np.vstack((poisson_sample_points, all_edge_seeds))
    all_points = np.unique(all_points, axis=0)
    
    tolerance = 0.1
    arr_rounded = np.round(all_points / tolerance) * tolerance
    unique_rounded, indices = np.unique(arr_rounded, axis=0, return_index=True)
    unique_coords = all_points[indices]
    all_points = unique_coords
    
    tri = Delaunay(all_points, qhull_options='Qt QJ')
    nodes, elements, internal_triangles = keep_internal_triangles(edges2D, lines2D, all_points, all_edge_seeds, tri.simplices)
    x_min, y_min = np.min(all_points, axis=0) -4
    x_max, y_max = np.max(all_points, axis=0) +4    
    plot_2D_face_mesh(nodes,elements,x_min,y_min,x_max,y_max)
    
    return nodes, elements
    

    
def plot_2D_face_mesh(nodes,elements,x_min,y_min,x_max,y_max):
    fig = plt.gcf()
    fig.set_size_inches(8., 8.)
    plt.axis([x_min,x_max,y_min,y_max])
    for element in elements:
        plt.plot([nodes[elements[element][0]]['coord'][0], nodes[elements[element][1]]['coord'][0]],
                 [nodes[elements[element][0]]['coord'][1], nodes[elements[element][1]]['coord'][1]],color='blue')
        plt.plot([nodes[elements[element][0]]['coord'][0], nodes[elements[element][2]]['coord'][0]],
                 [nodes[elements[element][0]]['coord'][1], nodes[elements[element][2]]['coord'][1]],color='blue')
        plt.plot([nodes[elements[element][2]]['coord'][0], nodes[elements[element][1]]['coord'][0]],
                 [nodes[elements[element][2]]['coord'][1], nodes[elements[element][1]]['coord'][1]],color='blue')
#    for node in nodes:
#        plt.plot(nodes[node][0], nodes[node][1], 'o', color='green')
    plt.show()    



def unit_vector(vector):
    return vector / np.linalg.norm(vector)



def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))



def check_if_plane_face(face_num, all_edge_seeds):
    normals = []
    origin = None
    print('\n\nface', face_num)
    vectors = []
    v1 = np.array(all_edge_seeds[-3])-np.array(all_edge_seeds[0])
    v2 = np.array(all_edge_seeds[3])-np.array(all_edge_seeds[0])
    cp = np.cross(v1,v2)
    for seed in all_edge_seeds:
        v = np.array(seed)-np.array(all_edge_seeds[0])
        dot_product = cp.dot(v)
        if abs(dot_product) > 0.01:
            print('Is NOT a plane face!')
            return cp, False

    print('Is a plane face!')
    return cp, True



# Function to transform a set of 3D points onto the XY-plane
def flip_to_xy_plane(face_num, part, all_edge_seeds, cp):
    # Calculate the centroid of the points
    centroid = np.mean(all_edge_seeds, axis=0)
    print('centroid:', centroid)

    # Subtract the centroid from the points to center them at the origin
    centered_seeds = all_edge_seeds - centroid
    centered_lines = {}
    for edge in part['faces'][face_num]['edges']:
        for line in part['faces'][face_num]['edges'][edge]['lines']:
            centered_lines[line] = part['lines'][line] - centroid

    normal = cp
    normal /= np.linalg.norm(normal)
    print('normal:', cp)
    zv = np.array([0,0,1])
    axis = np.cross(cp,zv)
    print('axis:', axis)
    # Calculate the angle between the normal vector and the z-axis
    angle = angle_between(cp,zv)
    angle_in_deg = 180*angle/np.pi
    print('angle:', angle)
    print('angle_in_deg:', 180*angle/np.pi)

    rotated_seeds = []
    rotated_lines = {}
    if angle_in_deg in [0.,180.]:
        zcoord = centered_seeds[0][2]
        rotated_seeds = centered_seeds[:,:2]
        print('face', face_num, 'already in xy-plane')
        for line in centered_lines:
            rotated_lines[line] = centered_lines[line][:,:2]
    else:
        for seed in centered_seeds:
            rotated_seeds.append(rotatePointAboutAxis(seed,centroid,centroid+axis,angle))
        zcoord = rotated_seeds[0][2]
        rotated_seeds = np.array(rotated_seeds)[:,:2]
        for line in centered_lines:
            rotated_lines[line] = []
            for line_point in centered_lines[line]:
                rotated_lines[line].append(rotatePointAboutAxis(line_point,centroid,centroid+axis,angle))
            rotated_lines[line] = np.array(rotated_lines[line])[:,:2]
    return rotated_seeds, rotated_lines, axis, angle, centroid, zcoord



def flip_back_to_3D(nodes, elements, axis, angle, centroid, zcoord, node_num1, element_num1):
    # add zcoord and renumber nodes and elements
    flipped_nodes = {}
    flipped_elements = {}
    for node in nodes:
        flipped_nodes[node_num1] = {'coord': np.hstack((nodes[node]['coord'],np.array([zcoord]))),
                                    'number': node_num1}
        nodes[node]['number'] = node_num1
        node_num1 += 1
    for element in elements:
        flipped_elements[element_num1] = [nodes[elements[element][0]]['number'],
                                          nodes[elements[element][1]]['number'],
                                          nodes[elements[element][2]]['number']]
        element_num1 += 1

    # flip all nodes using axis and angle
    # and then move them back to their original
    # position using centroid
    angle_in_deg = 180*angle/np.pi
    for node in flipped_nodes:
        if angle_in_deg in [0.,180.]:
            flipped_nodes[node]['coord'] = flipped_nodes[node]['coord']
        else:
            flipped_nodes[node]['coord'] = np.array(rotatePointAboutAxis(flipped_nodes[node]['coord'],centroid,centroid+axis,-angle))
        flipped_nodes[node]['coord'] += centroid
    
#    print('flipped_nodes:', flipped_nodes.keys())
#    print('flipped_elements:', flipped_elements)

    return flipped_nodes, flipped_elements



def plot_3D_face_mesh(nodes,elements):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for element in elements:
        l = np.vstack((nodes[elements[element][0]]['coord'], nodes[elements[element][1]]['coord']))
        ax.plot(l[:,0],l[:,1],l[:,2],color='b')
        l = np.vstack((nodes[elements[element][0]]['coord'], nodes[elements[element][2]]['coord']))
        ax.plot(l[:,0],l[:,1],l[:,2],color='b')
        l = np.vstack((nodes[elements[element][1]]['coord'], nodes[elements[element][2]]['coord']))
        ax.plot(l[:,0],l[:,1],l[:,2],color='b')
    plt.show()    








# bridson_sampling(max_min, face_nodes_kdtree, all_face_nodes, all_face_elements, max_element_size)
# max_min = [x_min, y_min, z_min, x_max, y_max, z_max]
def bridson_sampling(max_min, face_nodes_kdtree, face_nodes, all_face_nodes, face_elements, outside_point, radius, k=30):
    start_time = time.time()
    print('\n\nStarting bridson_sampling() ...')

    print('max_min:', max_min)
    width = max_min[3]-max_min[0]
    height = max_min[4]-max_min[1]
    depth = max_min[5]-max_min[2]
    # Step 0: Initialize the necessary variables
    cell_size = radius/np.sqrt(3)
    grid_width = int(np.ceil(width / cell_size))
    grid_height = int(np.ceil(height / cell_size))
    grid_depth = int(np.ceil(depth / cell_size))
    grid = np.empty((grid_width, grid_height, grid_depth), dtype=object)
    process_list = []
    sample_points = []
    internal_points = []
    for node in face_nodes:
        sample_points.append(node)
        grid[grid_index(node, cell_size)] = node

    # Step 1: Randomly select the initial point and add it to the grid and process list
#    initial_point = (random.uniform(0, width), random.uniform(0, height), random.uniform(0, depth))
    initial_point = random.choice(sample_points)
    process_list.append(initial_point)
#    sample_points.append(initial_point)
    grid[grid_index(initial_point, cell_size)] = initial_point

    # Step 2: While the process list is not empty, continue sampling
    while process_list:
        point = random.choice(process_list)  # Select a random point from process list
        process_list.remove(point)
        for _ in range(k):  # Try to find k valid points in the annulus
            new_point = get_point_in_annulus(point, radius)
            if 0 <= new_point[0] < width and 0 <= new_point[1] < height and 0 <= new_point[2] < depth:  # Check if the point is within the domain
                if is_valid_point(new_point, grid, cell_size, radius):
                    if is_within_volume(new_point, all_face_nodes, face_elements, max_min, outside_point):
                        process_list.append(new_point)
                        internal_points.append(new_point)
                        sample_points.append(new_point)
                        grid[grid_index(new_point, cell_size)] = new_point
    kdtrees = []
    kdtrees.append(face_nodes_kdtree)
    if len(internal_points) != 0:
        kdtrees.append(KDTree(np.array(internal_points)))
    sample_points = np.array(sample_points)

    end_time = time.time()
    print('... finished in {:.3f} seconds'.format(end_time-start_time))

    return sample_points, kdtrees

def grid_index(point, cell_size):
    return int(point[0] / cell_size), int(point[1] / cell_size), int(point[2] / cell_size)

def is_valid_point(new_point, grid, cell_size, radius):
    x, y, z = grid_index(new_point, cell_size)
    for i in range(max(0, x - 2), min(x + 3, grid.shape[0])):
        for j in range(max(0, y - 2), min(y + 3, grid.shape[1])):
            for k in range(max(0, z - 2), min(z + 3, grid.shape[2])):
                point = grid[i, j, k]
                if point is not None and np.linalg.norm(np.array(point) - np.array(new_point)) < radius:
                    return False
    return True

def get_point_in_annulus(center_point, radius):
    random_angle_phi = 2 * np.pi * random.random()
    random_angle_theta = np.arccos(2 * random.random() - 1)
    random_radius = random.uniform(radius, 2 * radius)
    random_x = center_point[0] + random_radius * np.sin(random_angle_theta) * np.cos(random_angle_phi)
    random_y = center_point[1] + random_radius * np.sin(random_angle_theta) * np.sin(random_angle_phi)
    random_z = center_point[2] + random_radius * np.cos(random_angle_theta)
    return random_x, random_y, random_z









def poisson_sphere_sampling(face_nodes, face_elements, r_max):
    start_time = time.time()
    print('\n\nStarting poisson_sphere_sampling() ...')

    x_min = min(face_nodes[i]['coord'][0] for i in face_nodes)
    y_min = min(face_nodes[i]['coord'][1] for i in face_nodes)
    z_min = min(face_nodes[i]['coord'][2] for i in face_nodes)
    x_max = max(face_nodes[i]['coord'][0] for i in face_nodes)
    y_max = max(face_nodes[i]['coord'][1] for i in face_nodes)
    z_max = max(face_nodes[i]['coord'][2] for i in face_nodes)
    max_min = [x_min,y_min,z_min,x_max,y_max,z_max]
    outside_point = np.array([x_min-1.,y_min-1.,z_min-1.])
    outside_point[0] -= 0.666*np.pi
    outside_point[1] -= 0.232*np.pi
    outside_point[2] -= 0.499*np.pi
    print('outside_point:', outside_point)

    boundary_samples = []
    for node in face_nodes:
        boundary_samples.append(face_nodes[node]['coord'].tolist())
    boundary_samples = np.array(boundary_samples)
    boundary_samples = np.unique(boundary_samples, axis=0)
    # create a kdtree for each layer
    kdtrees = []
    kdtrees.append(KDTree(boundary_samples))
    samples = boundary_samples.tolist()
#    internal_samples = []
    previous_samples = []
    
    # first layer of points inside boundary
    print('sampling first layer inside volume...')
    print('len(samples):', len(samples))
    for p in range(len(samples)-3):
        current_point = samples[p]
        current_rad = get_rad_3D(p, samples, r_max, kdtrees, 0)
#        print('current_rad:', current_rad)
        new_points = get_new_points_3D(current_point, current_rad, 6)
#        print('new_points:', new_points)
        for n in new_points:
            if keep_point_3D(n, 0.95*current_rad, samples, face_nodes, face_elements, max_min, outside_point, previous_samples, kdtrees, 0):
                previous_samples.append(n)
#                internal_samples.append(n)
                samples.append(n)
#                samples = np.vstack((samples,np.array(n)))
#    previous_samples = np.array(previous_samples)
    kdtrees.append(KDTree(np.array(previous_samples)))

    # add more layers of points inside until no more space    
    for l in range(4):
        print('sampling next layer inside volume...')
        print('len(previous_samples):', len(previous_samples))
        new_samples = []
        for p in range(len(previous_samples)-3):
            current_point = previous_samples[p]
            current_rad = get_rad_3D(p, previous_samples, r_max, kdtrees, l+1)
            new_points = get_new_points_3D(current_point, current_rad, 6)
            for n in new_points:
                if keep_point_3D(n, 0.95*current_rad, samples, face_nodes, face_elements, max_min, outside_point, new_samples, kdtrees, l+1):
                    new_samples.append(n)
#                    internal_samples.append(n)
                    samples.append(n)
        if len(new_samples) == 0:
            break
        previous_samples = new_samples
        kdtrees.append(KDTree(np.array(new_samples)))

    samples = np.array(samples)
#    internal_samples = np.array(internal_samples)
#    fig = plt.figure()
#    ax = plt.axes(projection='3d')
#    ax.scatter(internal_samples[:,0], internal_samples[:,1], internal_samples[:,2], color='g')
#    ax.scatter(boundary_samples[:,0], boundary_samples[:,1], boundary_samples[:,2], color='b')
#    plt.show()

    end_time = time.time()
    print('... finished in {:.3f} seconds'.format(end_time-start_time))

    return samples, kdtrees



def check_line_triangle_intersection(line_start, line_end, triangle):
    # Step 3: Calculate the direction vector of the line segment
    direction = line_end - line_start
    
    # Convert direction vector to floating-point type
    direction = direction.astype(float)
    
    # Step 4: Calculate the normal vector of the triangle
    AB = triangle[1] - triangle[0]
    AC = triangle[2] - triangle[0]
    normal = np.cross(AB, AC)
    
    # Step 5: Normalize the direction vector and the normal vector
#    direction /= np.linalg.norm(direction)
    normal /= np.linalg.norm(normal)
    
    # Step 6: Calculate the dot product
    dot = np.dot(direction, normal)
    
    # Check if line is parallel to triangle plane
    if abs(dot) < 1e-9:
        return False
    
    # Step 7: Calculate the vector T
    T = line_start - triangle[0]
    
    # Step 8: Calculate the parameter t
    t = -np.dot(T, normal) / dot
    
    # Step 9: Check if t is within [0, 1]
    if t < 0 or t > 1:
        return False
    
    # Step 10: Calculate the intersection point
    intersection_point = line_start + t * direction
    
    # Step 11: Check if the intersection point lies inside the triangle
    v0 = triangle[2] - triangle[0]
    v1 = triangle[1] - triangle[0]
    v2 = intersection_point - triangle[0]
    
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)
    
    denom = dot00 * dot11 - dot01**2
    
    u = (dot11 * dot02 - dot01 * dot12) / denom
    v = (dot00 * dot12 - dot01 * dot02) / denom
    
    return (u >= 0) and (v >= 0) and (u + v <= 1)
    

        
def is_within_volume(point, nodes, triangles, max_min, outside_point):
    intersections = 0
    for triangle in triangles:
        triangle_coords = np.array([[nodes[triangles[triangle][0]]['coord'][0],nodes[triangles[triangle][0]]['coord'][1],nodes[triangles[triangle][0]]['coord'][2]],
                                    [nodes[triangles[triangle][1]]['coord'][0],nodes[triangles[triangle][1]]['coord'][1],nodes[triangles[triangle][1]]['coord'][2]],
                                    [nodes[triangles[triangle][2]]['coord'][0],nodes[triangles[triangle][2]]['coord'][1],nodes[triangles[triangle][2]]['coord'][2]]])
        if point[0] < max_min[0] or \
           point[1] < max_min[1] or \
           point[2] < max_min[2] or \
           point[0] > max_min[3] or \
           point[1] > max_min[4] or \
           point[2] > max_min[5]:
               pass
        else:
            if check_line_triangle_intersection(outside_point, point, triangle_coords):
                intersections += 1
                    
    if intersections == 0:
        return False
    else:
        return intersections % 2 == 1


   
def get_rad_3D(index, samples, r_max, kdtrees, layer_num):
    nearest = kdtrees[layer_num].query_ball_point(samples[index],r_max)
    if len(nearest) == 0:
        return r_max
    else:
        r = [r_max]
        for p in nearest:
            r_tmp = 1.1*np.linalg.norm(kdtrees[layer_num].data[p] - samples[index])
            if r_tmp > 0.001:
                r.append(r_tmp)
        return min(r)



def get_new_points_3D(center, radius, num_points=14):

    points = [[center[0]-radius, center[1], center[2]],
              [center[0]+radius, center[1], center[2]],
              [center[0], center[1]-radius, center[2]],
              [center[0], center[1]+radius, center[2]],
              [center[0], center[1], center[2]-radius],
              [center[0], center[1], center[2]+radius]]
    if num_points == 14:
        pi_4 = 0.7071067811865475
        points.append([center[0]+pi_4*radius, center[1]+pi_4*radius, center[2]+pi_4*radius])
        points.append([center[0]+pi_4*radius, center[1]+pi_4*radius, center[2]-pi_4*radius])
        points.append([center[0]+pi_4*radius, center[1]-pi_4*radius, center[2]+pi_4*radius])
        points.append([center[0]+pi_4*radius, center[1]-pi_4*radius, center[2]-pi_4*radius])
        points.append([center[0]-pi_4*radius, center[1]+pi_4*radius, center[2]+pi_4*radius])
        points.append([center[0]-pi_4*radius, center[1]+pi_4*radius, center[2]-pi_4*radius])
        points.append([center[0]-pi_4*radius, center[1]-pi_4*radius, center[2]+pi_4*radius])
        points.append([center[0]-pi_4*radius, center[1]-pi_4*radius, center[2]-pi_4*radius])
    
    return points



def keep_point_3D(point, radius, samples, nodes, triangles, max_min, outside_point, new_samples, kdtrees, layer_num):
    # Check if the point intersects any other points
    for layer in kdtrees:
        if len(layer.query_ball_point(point,radius)) != 0:
            return False
    for n in new_samples:
        if np.linalg.norm([point[0] - n[0], point[1] - n[1], point[2] - n[2]]) <= radius:
            return False
    # check if the point is inside the volume
    if layer_num <= 1:
        if not is_within_volume(point, nodes, triangles, max_min, outside_point):
            return False

    return True



def is_valid_3D(point, radius, samples):
    # Check if the point intersects any other points
    for s in samples:
        if np.linalg.norm(point - s) < radius:
            return False
    return True



def are_points_on_plane(points, tolerance=1e-6):
    # Convert the points into a NumPy array
#    points = np.array(points)

    # Subtract the first point from all other points to get vectors
    vectors = points[1:] - points[0]

    # Calculate the cross product of the first two vectors
    cross_product = np.cross(vectors[0], vectors[1])

    # Check if the cross product is within the tolerance for all other vectors
    for vector in vectors[2:]:
        if np.abs(np.dot(cross_product, vector)) > tolerance:
            return False

    return True








#number of tets total: 922
#flat: 51
#internal_nodes: 425
#centerpoint_outside: 97


#number of tets total: 956
#flat_tets sorted:
#[False  True]
#[ 50 906]
#internal_tets sorted:
#[False  True]
#[485 471]
#len(center_points): 471
#tets_within_volume sorted:
#[False  True]
#[485 471]



   
def keep_internal_tetrahedrons_vectorized(all_nodes, face_nodes, face_elements, max_min, outside_point, kdtrees):
    start_time = time.time()
    print('\n\nStarting keep_internal_tetrahedrons_vectorized() ...')

    keep_tetrahedrons = []
    removed = 0
    count = 0
    tri = Delaunay(all_nodes,qhull_options='Qt')
    print('\n\nnumber of tets total:', len(tri.simplices))
    # Vectorize the volume calculation for all tetrahedrons
    V = np.einsum('ij,ij->i', all_nodes[tri.simplices[:, 0]] - all_nodes[tri.simplices[:, 3]], 
                     np.cross(all_nodes[tri.simplices[:, 1]] - all_nodes[tri.simplices[:, 3]], 
                              all_nodes[tri.simplices[:, 2]] - all_nodes[tri.simplices[:, 3]]))
    # Swap points if points are counter-clockwise
    mask = V > 0
    tri.simplices[mask] = tri.simplices[mask][:, [1, 0, 2, 3]]


    # Create an initial mask of all True values
    keep = np.ones(len(tri.simplices), dtype=bool)
    # Create an initial mask of all False values
    checked = np.zeros(len(tri.simplices), dtype=bool)
    
    # Check if tetrahedrons are flat
    flat_tets = are_points_on_plane_vectorized(all_nodes[tri.simplices[:, :4]], 0.1)

    # Update the mask for flat tets
    keep &= ~flat_tets
    checked = flat_tets
    
    # Check if any of the tetrahedron nodes are not on volume boundary and keep those
    internal_tets_full = np.zeros(len(flat_tets), dtype=bool)
    non_flat_tets = tri.simplices[~checked]
    internal_tets_full[~checked] = any_tet_nodes_not_on_boundary(all_nodes[non_flat_tets], kdtrees)
    
    # Update the mask for internal tets
    checked |= internal_tets_full
    
    # Calculate center points for non-flat and internal tetrahedrons
    center_points = np.mean(all_nodes[tri.simplices[~checked]], axis=1)
    
    # Prepare line start points and line end points
    line_start = outside_point
    line_ends = center_points
    
    # Prepare triangles
    triangles = face_nodes[face_elements]
    
    # Create a new mask of False values for tets within volume
    tets_within_volume_full = np.zeros(len(tri.simplices), dtype=bool)

    # Check if center points are within volume only for remaining tets
    tets_within_volume = is_within_volume_vectorized(line_start, line_ends, triangles)

    # Update the full mask only for remaining tets
    tets_within_volume_full[~checked] = tets_within_volume

    # Update the mask for tets within volume
    keep[~checked] &= tets_within_volume

    # Now you can use the mask to index tri.simplices
    keep_tetrahedrons = tri.simplices[keep]

    end_time = time.time()
    print('... finished in {:.3f} seconds'.format(end_time-start_time))
    return keep_tetrahedrons



def is_within_volume_vectorized(line_start, line_ends, triangles):
    # Initialization
    intersections = np.zeros(len(line_ends), dtype=bool)
    
    # Iterate over each center point (line end)
    for i, line_end in enumerate(line_ends):
        # Calculate the direction vector for this line segment
        direction = line_end - line_start
        
        # Iterate over each triangle
        for triangle in triangles:
            # Calculate the normal vector for this triangle
            AB = triangle[1] - triangle[0]
            AC = triangle[2] - triangle[0]
            normal = np.cross(AB, AC)
            
            # Normalize the normal vector
            normal /= np.linalg.norm(normal)
            
            # Calculate the dot product of the direction vector and the normal vector
            dot = np.dot(direction, normal)
            
            # If the line is not parallel to the triangle plane
            if np.abs(dot) > 1e-9:
                # Calculate the vector T for this line segment
                T = line_start - triangle[0]
                
                # Calculate the parameter t for this line segment
                t = -np.dot(T, normal) / dot
                
                # If t is within [0, 1], the line intersects the plane defined by the triangle
                if 0 <= t <= 1:
                    # Calculate the intersection point
                    intersection_point = line_start + t * direction
                    
                    # Check if the intersection point lies inside the triangle
                    v0 = triangle[2] - triangle[0]
                    v1 = triangle[1] - triangle[0]
                    v2 = intersection_point - triangle[0]
                    
                    dot00 = np.dot(v0, v0)
                    dot01 = np.dot(v0, v1)
                    dot02 = np.dot(v0, v2)
                    dot11 = np.dot(v1, v1)
                    dot12 = np.dot(v1, v2)
                    
                    denom = dot00 * dot11 - dot01**2
                    u = (dot11 * dot02 - dot01 * dot12) / denom
                    v = (dot00 * dot12 - dot01 * dot02) / denom
                    
                    # If the intersection point is inside the triangle
                    if u >= 0 and v >= 0 and u + v <= 1:
                        # Update the intersections array for this line end
                        intersections[i] = ~intersections[i]

    return intersections




def is_within_volume_vectorized(line_start, line_ends, triangles):
    # Initialization
    intersections = np.zeros(len(line_ends), dtype=bool)

    # Calculate the direction vector for each line segment
    directions = line_ends - line_start

    # Iterate over each triangle
    for triangle in triangles:
        # Calculate the normal vector for this triangle
        AB = triangle[1] - triangle[0]
        AC = triangle[2] - triangle[0]
        normal = np.cross(AB, AC)

        # Normalize the normal vector
        normal /= np.linalg.norm(normal)

        # Calculate the dot product of the direction vectors and the normal vector
        dots = np.dot(directions, normal)

        # Calculate the vector T for each line segment
        Ts = line_start - triangle[0]

        # Calculate the parameters t for each line segment
        ts = -np.dot(Ts, normal) / dots

        # The lines that intersect the plane defined by the triangle have t within [0, 1]
        t_mask = np.logical_and(0 <= ts, ts <= 1)

        # Calculate the intersection points
        intersection_points = line_start + np.expand_dims(ts[t_mask], axis=1) * directions[t_mask]

        # Check if the intersection points lie inside the triangle
        v0 = triangle[2] - triangle[0]
        v1 = triangle[1] - triangle[0]
        v2 = intersection_points - triangle[0]

        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot02 = np.dot(v0, v2.T)
        dot11 = np.dot(v1, v1)
        dot12 = np.dot(v1, v2.T)

        denom = dot00 * dot11 - dot01**2
        us = (dot11 * dot02 - dot01 * dot12) / denom
        vs = (dot00 * dot12 - dot01 * dot02) / denom

        # The intersection points that are inside the triangle
        uv_mask = np.logical_and.reduce([us >= 0, vs >= 0, us + vs <= 1])

        # Update the intersections array for these line ends
        intersections[t_mask] ^= uv_mask

    return intersections

#def is_within_volume_vectorized(points, nodes, triangles, max_min, outside_point):
    # First, we create an array of all triangle coordinates
#    triangle_coords = np.array([[
#        nodes[triangles[triangle][i]]['coord'] 
#        for i in range(3)
#    ] for triangle in triangles])

    # Now, we can use numpy's functions to perform computations on the entire array
#    within_volume = np.zeros(len(points), dtype=bool)
#    for i, point in enumerate(points):
#        intersections = np.array([
#            check_line_triangle_intersection(outside_point, point, coords)
#            for coords in triangle_coords
#        ])

        # Count intersections
#        intersection_count = np.count_nonzero(intersections)

#       within_volume[i] = intersection_count % 2 == 1

#    return within_volume


def check_line_triangle_intersection(line_start, line_end, triangle):
    # Step 3: Calculate the direction vector of the line segment
    direction = line_end - line_start
    
    # Convert direction vector to floating-point type
    direction = direction.astype(float)
    
    # Step 4: Calculate the normal vector of the triangle
    AB = triangle[1] - triangle[0]
    AC = triangle[2] - triangle[0]
    normal = np.cross(AB, AC)
    
    # Step 5: Normalize the direction vector and the normal vector
    normal /= np.linalg.norm(normal)
    
    # Step 6: Calculate the dot product
    dot = np.dot(direction, normal)
    
    # Check if line is parallel to triangle plane
    if abs(dot) < 1e-9:
        return False
    
    # Step 7: Calculate the vector T
    T = line_start - triangle[0]
    
    # Step 8: Calculate the parameter t
    t = -np.dot(T, normal) / dot
    
    # Step 9: Check if t is within [0, 1]
    if t < 0 or t > 1:
        return False
    
    # Step 10: Calculate the intersection point
    intersection_point = line_start + t * direction
    
    # Step 11: Check if the intersection point lies inside the triangle
    v0 = triangle[2] - triangle[0]
    v1 = triangle[1] - triangle[0]
    v2 = intersection_point - triangle[0]
    
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)
    
    denom = dot00 * dot11 - dot01**2
    
    u = (dot11 * dot02 - dot01 * dot12) / denom
    v = (dot00 * dot12 - dot01 * dot02) / denom
    
    return (u >= 0) and (v >= 0) and (u + v <= 1)


def any_tet_nodes_not_on_boundary(tet_nodes, kdtrees):
    nodes_boundary_check = np.array([
        [len(kdtrees[0].query_ball_point(node, 0.01)) for node in tet]
        for tet in tet_nodes
    ])
    return np.any(nodes_boundary_check == 0, axis=1)


def are_points_on_plane_vectorized(points, tol):
    AB = points[:, 1] - points[:, 0]
    AC = points[:, 2] - points[:, 0]
    AD = points[:, 3] - points[:, 0]
    
    cross_product = np.cross(AB, AC)
    dot_product = np.einsum('ij,ij->i', cross_product, AD)
    volume = np.abs(dot_product) / 6  # volume of tetrahedron

    return volume <= tol

#def is_within_volume_vectorized(point, nodes, triangles, max_min, outside_point):
    # First, we create an array of all triangle coordinates
#    triangle_coords = np.array([[
#        nodes[triangles[triangle][i]]['coord'] 
#        for i in range(3)
#    ] for triangle in triangles])

    # Now, we can use numpy's functions to perform computations on the entire array
#    intersections = np.array([
#        check_line_triangle_intersection(outside_point, point, coords)
#        for coords in triangle_coords
#    ])

    # Count intersections
#    intersection_count = np.count_nonzero(intersections)

#    return intersection_count % 2 == 1










def keep_internal_tetrahedrons(all_nodes, face_nodes, face_elements, max_min, outside_point, kdtrees):
    start_time = time.time()
    print('\n\nStarting keep_internal_tetrahedrons() ...')
    keep_tetrahedrons = []
    removed = 0
    count = 0
    print('len(all_nodes):', len(all_nodes))
    print('all_nodes:', all_nodes)
    tri = Delaunay(all_nodes,qhull_options='Qt')
    print('number of tets total:', len(tri.simplices))
    for i, simplex in enumerate(tri.simplices):
        A, B, C, D = all_nodes[simplex]
    
        V = np.dot(A-D, np.cross(B-D, C-D))
    
        if V > 0:
            # If the points are counter-clockwise, we swap two points to make it clockwise.
            tri.simplices[i][1], tri.simplices[i][2] = tri.simplices[i][2], tri.simplices[i][1]
    flat = 0
    internal_nodes = 0
    centerpoint_outside = 0
    for tet in tri.simplices:
        count += 1
        if count%500 == 0:
            print('count:', count)
            print('removed:', removed)
        p1 = all_nodes[tet[0]]
        p2 = all_nodes[tet[1]]
        p3 = all_nodes[tet[2]]
        p4 = all_nodes[tet[3]]
        points = np.array([[all_nodes[tet[0]][0], all_nodes[tet[0]][1], all_nodes[tet[0]][2]],
                           [all_nodes[tet[1]][0], all_nodes[tet[1]][1], all_nodes[tet[1]][2]],
                           [all_nodes[tet[2]][0], all_nodes[tet[2]][1], all_nodes[tet[2]][2]],
                           [all_nodes[tet[3]][0], all_nodes[tet[3]][1], all_nodes[tet[3]][2]]])
        # check if tetrahedron is flat
        if are_points_on_plane(points,0.1):
            removed += 1
            flat += 1
            continue
        # check if any nodes not on volume boundary
        if len(kdtrees[0].query_ball_point(p1,0.01)) == 0 or len(kdtrees[0].query_ball_point(p2,0.01)) == 0 or \
           len(kdtrees[0].query_ball_point(p3,0.01)) == 0 or len(kdtrees[0].query_ball_point(p4,0.01)) == 0:
            keep_tetrahedrons.append(tet)
            internal_nodes += 1
            continue
        # check if center point of tetrahedron is outside volume
        center_point = np.mean(points,axis=0)
        if not is_within_volume(center_point, face_nodes, face_elements, max_min, outside_point):
            centerpoint_outside += 1
            removed += 1
        else:
            keep_tetrahedrons.append(tet)
    print('flat:', flat)
    print('internal_nodes:', internal_nodes)
    print('centerpoint_outside:', centerpoint_outside)
    end_time = time.time()
    print('... finished in {:.3f} seconds'.format(end_time-start_time))
    return keep_tetrahedrons     




def fuseNodes(nodes, elements):
    new_nodes = {}
    new_elements = {}
    for node in range(len(nodes)):
        new_nodes[node+1] = nodes[node]
    for element in range(len(elements)):
        new_elements[element+1] = {'nodes': [elements[element][0]+1, elements[element][1]+1, elements[element][2]+1, elements[element][3]+1],
                                   'number': element+1}
#    print('new_nodes:', new_nodes)
#    print('new_elements:', new_elements)
    tolerance = 0.1
    fused_nodes = []
    for node1 in new_nodes:
        for node2 in new_nodes:
            if node1 != node2:
                distance = np.sqrt((new_nodes[node2][0] - new_nodes[node1][0])**2 + \
    					   		      (new_nodes[node2][1] - new_nodes[node1][1])**2 + \
    					   		      (new_nodes[node2][2] - new_nodes[node1][2])**2)
                if distance <= tolerance:
                    fused_nodes.append(tuple(sorted([node1,node2])))
    fused_nodes = sorted(set(fused_nodes))
    to_be_deleted = sorted([i[1] for i in fused_nodes])
    confirm_to_be_deleted = []
    for element in new_elements:
        node_in_element = False
        elm_nodes = []
        for node in range(len(new_elements[element]['nodes'])):
            if new_elements[element]['nodes'][node] in to_be_deleted:
                # also check that no nodes that are in
                # the same element are fused together
                if len(elm_nodes) == 0:
                    for elm_node in new_elements[element]['nodes']:
                        elm_nodes.append(new_elements[element]['nodes'][node])
                for node_combo in range(len(fused_nodes))[::-1]:
                    if new_elements[element]['nodes'][node] in fused_nodes[node_combo]:
                        if fused_nodes[node_combo][0] in elm_nodes and fused_nodes[node_combo][1] in elm_nodes:
                            print('\tCannot fuse together two nodes that are in the same element:')
                            print('\tElement:', mesh.elements[element].number)
                            print('\nNodes:', fused_nodes[node_combo][0], 'and', fused_nodes[node_combo][1])
                            pass
                        else:
                            new_elements[element]['nodes'][node] = fused_nodes[node_combo][0]
                            confirm_to_be_deleted.append(fused_nodes[node_combo][1])
    print('Deleting nodes ...')
    deleted = []
    for node in to_be_deleted:
        if node in confirm_to_be_deleted:
            if node in new_nodes:
                del new_nodes[node]
                deleted.append(node)
    return new_nodes, new_elements











def plot_tetrahedrons(points, tets):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    count = 0
    for tet in tets:
        pts = points[tet, :]
        c_l = 'g'
        l_w = '0.1'
        if count == 1:
            c_l = 'r'
            l_w = '1.'
        ax.plot3D(pts[[0,1],0], pts[[0,1],1], pts[[0,1],2], color=c_l, lw=l_w)
        ax.plot3D(pts[[0,2],0], pts[[0,2],1], pts[[0,2],2], color=c_l, lw=l_w)
        ax.plot3D(pts[[0,3],0], pts[[0,3],1], pts[[0,3],2], color=c_l, lw=l_w)
        ax.plot3D(pts[[1,2],0], pts[[1,2],1], pts[[1,2],2], color=c_l, lw=l_w)
        ax.plot3D(pts[[1,3],0], pts[[1,3],1], pts[[1,3],2], color=c_l, lw=l_w)
        ax.plot3D(pts[[2,3],0], pts[[2,3],1], pts[[2,3],2], color=c_l, lw=l_w)
        count += 1

#    ax.scatter(points[:,0], points[:,1], points[:,2], color='b')
    plt.show()    





# function to create 3D mesh from part (edges, lines and seeds)
def mesh_part(part, max_element_size):
    nodes = {}
    elements = {}

    # first mesh faces
    all_face_nodes = {}
    all_face_elements = {}
    for face_num in part['faces']:
        first = True
        for edge in part['faces'][face_num]['edges']:
            for line in range(len(part['faces'][face_num]['edges'][edge]['lines'])):
                if line == 0:
                    edge_seeds = part['seeds'][part['faces'][face_num]['edges'][edge]['lines'][line]]
                else:
                    edge_seeds = np.vstack((edge_seeds, part['seeds'][part['faces'][face_num]['edges'][edge]['lines'][line]]))
            if first:
                all_edge_seeds = edge_seeds
                first = False
            else:
                all_edge_seeds = np.vstack((all_edge_seeds,edge_seeds))
        all_edge_seeds = np.unique(all_edge_seeds, axis=0)

        # check if plane face
        cp, is_plane = check_if_plane_face(face_num, all_edge_seeds)

        if is_plane:
            # flip face to XY-plane
            all_edge_seeds, lines2D, axis, angle, centroid, zcoord = flip_to_xy_plane(face_num,part,all_edge_seeds,cp)
            all_edge_seeds = np.unique(all_edge_seeds, axis=0)
            edges2D = {}
            for edge in part['faces'][face_num]['edges']:
                edges2D[edge] = part['faces'][face_num]['edges'][edge]['lines']
            nodes, elements = create_2D_face_mesh(edges2D, lines2D, all_edge_seeds, max_element_size)
#            print('nodes:', nodes)
#            print('elements:', elements)
        
        # flip 2D mesh back to face plane
        if is_plane:
            print('all_face_nodes:', len(all_face_nodes))
            print('all_face_elements:', len(all_face_elements))
            if len(all_face_nodes) == 0:
                node_num1 = 1
                element_num1 = 1
            else:
                node_num1 = max(all_face_nodes.keys())+1
                element_num1 = max(all_face_elements.keys())+1
            print('node_num1:', node_num1)
            print('element_num1:', element_num1)
            nodes, elements = flip_back_to_3D(nodes, elements, axis, angle, centroid, zcoord, node_num1, element_num1)
            for node in nodes:
                all_face_nodes[node] = nodes[node]
            for element in elements:
                all_face_elements[element] = elements[element]


#    print('all_face_nodes:', all_face_nodes.keys())        
#    print('all_face_element:', all_face_elements)        
#    plot_3D_face_mesh(all_face_nodes, all_face_elements)
    all_face_nodes_numpy = np.zeros((len(all_face_nodes),3))
#    all_face_triangles_numpy = np.zeros((len(all_face_elements),3))
    all_face_triangles_numpy = []
    for node in all_face_nodes:
        all_face_nodes_numpy[node-1] = all_face_nodes[node]['coord']
#    print('all_face_nodes_numpy:', all_face_nodes_numpy)
    for element in all_face_elements:
        all_face_triangles_numpy.append([x-1 for x in all_face_elements[element]])
#        all_face_triangles_numpy[element-1] = [int(x) for x in all_face_elements[element]]
#    print('all_face_triangles_numpy:', all_face_triangles_numpy)

    # second mesh interior
    test_points = np.array([[ 40., 10., 10.],
                            [ 14., 10., 33.],
                            [ 30., 20., 20.],
                            [  5., 49., 15.]])
#    test_points = np.random.rand(10,3)
#    test_points[:,0] *= 60
#    test_points[:,0] -= 5
#    test_points[:,1] *= 60
#    test_points[:,1] -= 5
#    test_points[:,2] *= 30
#    test_points[:,2] -= 5

    x_min = min(all_face_nodes[i]['coord'][0] for i in all_face_nodes)
    y_min = min(all_face_nodes[i]['coord'][1] for i in all_face_nodes)
    z_min = min(all_face_nodes[i]['coord'][2] for i in all_face_nodes)
    x_max = max(all_face_nodes[i]['coord'][0] for i in all_face_nodes)
    y_max = max(all_face_nodes[i]['coord'][1] for i in all_face_nodes)
    z_max = max(all_face_nodes[i]['coord'][2] for i in all_face_nodes)
    outside_point = np.array([x_min-1.,y_min-1.,z_min-1.])
    outside_point[0] -= 0.666*np.pi
    outside_point[1] -= 0.232*np.pi
    outside_point[2] -= 0.499*np.pi
    print('\n\n\noutside_point:', outside_point)

    max_min = [x_min, y_min, z_min, x_max, y_max, z_max]

    count = 1
    for n in all_face_nodes:
        if count == 1:
            face_nodes = all_face_nodes[n]['coord']
        else:
            face_nodes = np.vstack((face_nodes,all_face_nodes[n]['coord']))
        count += 1
    face_nodes_kdtree = KDTree(face_nodes)

    for p in test_points:
        print('point', p, 'within volume:', is_within_volume(p, all_face_nodes, all_face_elements, max_min, outside_point))
#    all_nodes, kdtrees = poisson_sphere_sampling(all_face_nodes, all_face_elements, max_element_size)
    all_nodes, kdtrees = bridson_sampling(max_min, face_nodes_kdtree, face_nodes, all_face_nodes, all_face_elements, outside_point, max_element_size)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter(test_points[:,0], test_points[:,1], test_points[:,2], color='r')
    ax.scatter(all_nodes[:,0], all_nodes[:,1], all_nodes[:,2], color='b')
    plt.show()

    elements = keep_internal_tetrahedrons_vectorized(all_nodes, all_face_nodes_numpy, all_face_triangles_numpy, max_min, outside_point, kdtrees)



    nodes, elements = fuseNodes(all_nodes,elements)

    # third plot the tetrahedrons
#    plot_tetrahedrons(all_nodes, elements)



#    nodes = np.array([[0,0,0],[1,1,1],[1,0,0],[0,1,0]])
#    elements = {1: {'nodes': [0,1,2,3]}}

    return nodes, elements













# import part
from test_part7 import *






# TEST THE 2D MESHER WITH THIS CODE BLOCK
#--------------------------------------------------------------------------#
part = plate_with_hole
# get 2D lines and seeds
#edges2D = {}
#lines2D = {}
#facenum = 17
#for edge in part['faces'][facenum]['edges']:
#    edges2D[edge] = part['faces'][facenum]['edges'][edge]['lines']
#    for line in part['faces'][facenum]['edges'][edge]['lines']:
#        lines2D[line] = np.array(part['lines'][line])[:,:2]
#seeds2D = {}
#for edge in part['faces'][facenum]['edges']:
#    for line in part['faces'][facenum]['edges'][edge]['lines']:
#        seeds2D[line] = np.array(part['seeds'][line])[:,:2]
#for line in range(len(seeds2D)):
#    if line == 0:
#        all_edge_seeds = seeds2D[list(seeds2D.keys())[line]]
#    else:
#        all_edge_seeds = np.vstack((all_edge_seeds,seeds2D[list(seeds2D.keys())[line]]))
#all_edge_seeds = np.unique(all_edge_seeds, axis=0)
#print('edges2D', edges2D)
#print('lines2D', lines2D)
#print('seeds2D', seeds2D)
#print('all_edge_seeds', all_edge_seeds)
#max_element_size = 7
#create_2D_face_mesh(edges2D, lines2D, all_edge_seeds, max_element_size)
#--------------------------------------------------------------------------#




# TEST THE 3D MESHER WITH THIS CODE BLOCK
#--------------------------------------------------------------------------#
part = plate_with_hole
max_element_size = 4
nodes, elements = mesh_part(part, max_element_size)

try:
    fobj_out1 = open('tet_mesh1.sol','w+')

except OSError as e:
    print('\n\n  *** ERROR!!!', e)

else:
    fobj_out1.write('#\n#\n#\n')
    for n in nodes:
#        if is_within_boundary(all_points_poisson[n],part):
        fobj_out1.write('NODE, '+str(n)+', '+str(nodes[n][0])+', '+str(nodes[n][1])+', '+str(nodes[n][2])+'\n')
    for elem in elements:
        fobj_out1.write('ELEMENT, TET4N, '+str(elem)+', 0, '+str(elements[elem]['nodes'][0])+', '+str(elements[elem]['nodes'][1])+', '+str(elements[elem]['nodes'][2])+', '+str(elements[elem]['nodes'][3])+'\n')
        elem += 1
    fobj_out1.close()



#--------------------------------------------------------------------------#


