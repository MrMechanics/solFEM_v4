#
#
#    reader.py
#  -----------
#
#    This is the reader module. It has objects used for reading
#   mesh files (*.sol, *.dat, *.inp, *.bdf), model files (*.mdl)
#   and step-files (*.step, *.STEP, *.stp, *.STP).
#


import os
import sys
sys.path.insert(1, '../work_directory')
import re
from timeit import time





class StepFileData(object):
    '''
Class for geometry input read from step-files. This
object first reads a step-file, possibly does some
checks, and then stores the data from that step-file
for easy access by the Part object.
'''
    def __init__(self,filename):

        self.filename = filename
        for i in range(len(filename)):
            if filename[-i-1] == '/' or filename[-i-1] == '\\':
                if filename[-4:] in ['STEP', 'step']:
                    self.name = str(filename[-i:-4])
                if filename[-3:] in ['STP', 'stp']:
                    self.name = str(filename[-i:-4])
                break
            else:
                if filename[-4:] in ['STEP', 'step']:
                    self.name = str(filename[:-5])
                if filename[-3:] in ['STP', 'stp']:
                    self.name = str(filename[:-5])
        
        print('\nReading step file... ', end='')
        self.parse_step_file()
        print('Finished!')

        
        
    def parse_step_file(self):
        self.round_to_decimal_place = 3
        
        self.cartesian_point = {}
        self.direction = {}
        self.vertex_point = {}
        self.vector = {}
        self.axis1_placement = {}
        self.axis2_placement_3D = {}
        self.line = {}
        self.circle = {}
        self.ellipse = {}
        self.b_spline_curve = {}
        self.edge_curve = {}
        self.oriented_edge = {}
        self.edge_loop = {}
        self.plane = {}
        self.cylindrical_surface = {}
        self.toroidal_surface = {}
        self.conical_surface = {}
        self.surface_of_revolution = {}
        self.surface_of_linear_extrusion = {}
        self.b_spline_surface = {}
        self.face_bound = {}
        self.face_outer_bound = {}
        self.advanced_face = {}
        self.closed_shell = {}
        self.open_shell = {}
        self.manifold_solid_brep = {}
        self.shell_based_surface_model = {}
        self.advanced_brep_shape_representation = {}
        self.manifold_surface_shape_representation = {}
    
        with open(self.filename, 'r') as f:
            entry = ""
            for line in f:
                entry += line.strip()
                if ';' in line:  # Indicates the end of an entry
                    parsed_data = self.parse_line(entry)
                    if parsed_data:
                        data_type, index, data = parsed_data
    
                        if data_type == "CARTESIAN_POINT":
                            self.cartesian_point[index] = data
                        elif data_type == "DIRECTION":
                            self.direction[index] = data
                        elif data_type == "VERTEX_POINT":
                            self.vertex_point[index] = data
                        elif data_type == "VECTOR":
                            self.vector[index] = data
                        elif data_type == "AXIS1_PLACEMENT":
                            self.axis1_placement[index] = data 
                        elif data_type == "AXIS2_PLACEMENT_3D":
                            self.axis2_placement_3D[index] = data 
                        elif data_type == "LINE":
                            self.line[index] = data
                        elif data_type == "CIRCLE":
                            self.circle[index] = data 
                        elif data_type == "ELLIPSE":
                            self.ellipse[index] = data 
                        elif data_type == "B_SPLINE_CURVE_WITH_KNOTS":
                            self.b_spline_curve[index] = data
                        elif data_type == "EDGE_CURVE":
                            self.edge_curve[index] = data 
                        elif data_type == "ORIENTED_EDGE":
                            self.oriented_edge[index] = data 
                        elif data_type == "EDGE_LOOP":
                            self.edge_loop[index] = data 
                        elif data_type == "PLANE":
                            self.plane[index] = data 
                        elif data_type == "CYLINDRICAL_SURFACE":
                            self.cylindrical_surface[index] = data 
                        elif data_type == "TOROIDAL_SURFACE":
                            self.toroidal_surface[index] = data 
                        elif data_type == "CONICAL_SURFACE":
                            self.conical_surface[index] = data 
                        elif data_type == "SURFACE_OF_REVOLUTION":
                            self.surface_of_revolution[index] = data 
                        elif data_type == "SURFACE_OF_LINEAR_EXTRUSION":
                            self.surface_of_linear_extrusion[index] = data 
                        elif data_type == "B_SPLINE_SURFACE":
                            self.b_spline_surface[index] = data 
                        elif data_type == "FACE_BOUND":
                            self.face_bound[index] = data 
                        elif data_type == "FACE_OUTER_BOUND":
                            self.face_outer_bound[index] = data 
                        elif data_type == "ADVANCED_FACE":
                            self.advanced_face[index] = data 
                        elif data_type == "CLOSED_SHELL":
                            self.closed_shell[index] = data 
                        elif data_type == "OPEN_SHELL":
                            self.open_shell[index] = data 
                        elif data_type == "MANIFOLD_SOLID_BREP":
                            self.manifold_solid_brep[index] = data 
                        elif data_type == "SHELL_BASED_SURFACE_MODEL":
                            self.shell_based_surface_model[index] = data 
                        elif data_type == "ADVANCED_BREP_SHAPE_REPRESENTATION":
                            self.advanced_brep_shape_representation[index] = data 
                        elif data_type == "MANIFOLD_SURFACE_SHAPE_REPRESENTATION":
                            self.manifold_surface_shape_representation[index] = data 
                    
                    entry = ""  # Reset the entry for the next round
                    
                    
    def parse_line(self,line):
        match = re.match(r'#(\d+)\s*=\s*(\w+)\s*\(\s*\'[^\']*\'\s*,(.+)\)\s*;', line)
        if not match:
            if 'B_SPLINE_SURFACE' in line:
                index_match = re.match(r'#(\d+)', line)
                if not index_match:
                    return None
                index = int(index_match.group(1))
                return 'B_SPLINE_SURFACE', index, self.parse_b_spline_surface(line)
            if 'BOUNDED_CURVE' in line:
                index_match = re.match(r'#(\d+)', line)
                if not index_match:
                    return None
                index = int(index_match.group(1))
                return 'B_SPLINE_CURVE_WITH_KNOTS', index, self.parse_b_spline_curve(line)
            return None
    
        index = int(match.group(1))
        data_type = match.group(2)
        data = match.group(3)
    
        if data_type in ["CARTESIAN_POINT", "DIRECTION"]:
            coords = tuple(round(float(val),self.round_to_decimal_place) for val in re.findall(r'(-?\d+\.?\d*[eE]?[-+]?\d*)', data))
            return data_type, index, coords
    
        elif data_type == "VERTEX_POINT":
            vertex = int(re.search(r'#(\d+)', data).group(1))
            return data_type, index, vertex
    
        elif data_type == "VECTOR":
            components = [comp.strip() for comp in data.split(',')]
            ref_point = int(re.search(r'#(\d+)', components[0]).group(1))
            magnitude = round(float(components[1]),self.round_to_decimal_place)
            return data_type, index, (ref_point, magnitude)
    
        elif data_type == "AXIS1_PLACEMENT":
            components = [comp.strip() for comp in data.split(',')]
            refs = tuple(int(val) for val in re.findall(r'#(\d+)', data))
            return data_type, index, refs
        
        elif data_type == "AXIS2_PLACEMENT_3D":
            components = [comp.strip() for comp in data.split(',')]
            refs = tuple(int(val) for val in re.findall(r'#(\d+)', data))
            return data_type, index, refs
        
        elif data_type == "LINE":
            points = tuple(int(val) for val in re.findall(r'#(\d+)', data))
            return data_type, index, points
    
        elif data_type == "CIRCLE":
            components = [comp.strip() for comp in data.split(',')]
            ref_point = int(re.search(r'#(\d+)', components[0]).group(1))
            radius = round(float(components[1]),self.round_to_decimal_place)
            return data_type, index, (ref_point, radius)
    
        elif data_type == "ELLIPSE":
            components = [comp.strip() for comp in data.split(',')]
            ref_point = int(re.search(r'#(\d+)', components[0]).group(1))
            major_radius = round(float(components[1]),self.round_to_decimal_place)
            minor_radius = round(float(components[2]),self.round_to_decimal_place)
            return data_type, index, (ref_point, major_radius, minor_radius)
    
        elif data_type == "B_SPLINE_CURVE_WITH_KNOTS":
            # Extract degree
            degree = int(re.search(r'\d+', data).group())
            
            # Use the regex pattern to capture the control points, booleans, knot multiplicities, and knot values
            pattern = r"\(([^)]*)\)[\s,]*\.UNSPECIFIED\.[\s,]*\.(F|T)\.[\s,]*\.(F|T)\.[\s,]*\(([^)]*)\)[\s,]*\(([^)]*)\)"
            matches = re.search(pattern, data)
            
            # If not matched, raise error
            if not matches:
                raise ValueError("Unexpected data format for B_SPLINE_CURVE_WITH_KNOTS")
            
            # Extract control points
            control_points_matches = re.findall(r'#(\d+)', matches.group(1))
            control_points = tuple(map(int, control_points_matches))
            
            # Extract planar and closed booleans
            planar = True if matches.group(2) == "T" else False
            closed = True if matches.group(3) == "T" else False
            
            # Extract knot multiplicities
            knot_multiplicities_matches = re.findall(r'(\d+)', matches.group(4))
            knot_multiplicities = tuple(map(int, knot_multiplicities_matches))
            
            # Extract knot values
            knot_value_matches = re.findall(r'(-?\d+\.\d+e?-?\d+)', matches.group(5))
            knot_values = tuple(map(float, knot_value_matches))
            
            return data_type, index, {
                'degree': degree, 
                'control_points': control_points, 
                'planar': planar, 
                'closed': closed, 
                'knot_multiplicities': knot_multiplicities, 
                'knot_values': knot_values
            }
    
        elif data_type == "EDGE_CURVE":
            components = [comp.strip() for comp in data.split(',')]
        
            # Extracting refs
            try:
                ref1 = int(components[0][1:]) if "#" in components[0] else None
                ref2 = int(components[1][1:]) if "#" in components[1] else None
                ref3 = int(components[2][1:]) if "#" in components[2] else None
            except ValueError:
                return None
    
            bool_val = components[3] == ".T."
            result = (ref1, ref2, ref3, bool_val)
            return data_type, index, result
    
        elif data_type == "ORIENTED_EDGE":
            components = [comp.strip() for comp in data.split(',')]
        
            # Extracting refs
            try:
                ref1 = int(components[0][1:]) if "#" in components[0] else None
                ref2 = int(components[1][1:]) if "#" in components[1] else None
                ref3 = int(components[2][1:]) if "#" in components[2] else None
            except ValueError:
                return None
        
            bool_val = components[3] == ".T."
            result = (ref1, ref2, ref3, bool_val)
            return data_type, index, result
        
        elif data_type == "EDGE_LOOP":
            # Remove newline characters for consistent parsing
            cleaned_data = data.replace('\n', '').replace('\r', '')
            
            # Extract all references between the parentheses
            refs_matches = re.search(r'\(([^)]+)\)', cleaned_data)
            if not refs_matches:
                return None  # or handle this case as needed
            refs_str = refs_matches.group(1)
            refs = tuple(int(val) for val in re.findall(r'#(\d+)', refs_str))
            
            return data_type, index, refs
        
        elif data_type == "PLANE":
            plane = int(re.search(r'#(\d+)', data).group(1))
            return data_type, index, plane
    
        elif data_type == "CYLINDRICAL_SURFACE":
            components = [comp.strip() for comp in data.split(',')]
            ref_point = int(re.search(r'#(\d+)', components[0]).group(1))
            radius = round(float(components[1]),self.round_to_decimal_place)
            return data_type, index, (ref_point, radius)
    
        elif data_type == "TOROIDAL_SURFACE":
            components = [comp.strip() for comp in data.split(',')]
            ref_point = int(re.search(r'#(\d+)', components[0]).group(1))
            major_radius = round(float(components[1]),self.round_to_decimal_place)
            minor_radius = round(float(components[2]),self.round_to_decimal_place)
            return data_type, index, (ref_point, major_radius, minor_radius)
    
        elif data_type == "CONICAL_SURFACE":
            components = [comp.strip() for comp in data.split(',')]
            ref_point = int(re.search(r'#(\d+)', components[0]).group(1))
            radius = round(float(components[1]),self.round_to_decimal_place)
            semi_angle = round(float(components[2]),self.round_to_decimal_place)
            return data_type, index, (ref_point, radius, semi_angle)
    
        elif data_type == "SURFACE_OF_REVOLUTION":
            components = [comp.strip() for comp in data.split(',')]
            refs = tuple(int(val) for val in re.findall(r'#(\d+)', data))
            return data_type, index, refs
        
        elif data_type == "SURFACE_OF_LINEAR_EXTRUSION":
            components = [comp.strip() for comp in data.split(',')]
            refs = tuple(int(val) for val in re.findall(r'#(\d+)', data))
            return data_type, index, refs
        
        elif data_type == "FACE_BOUND":
            # Remove all newline characters
            cleaned_data = data.replace('\n', '').replace('\r', '')
            
            # Extract the reference using regex
            ref_match = re.search(r'#(\d+)', cleaned_data)
            if not ref_match:
                return None  # or handle this case as needed
            ref = int(ref_match.group(1))
            
            bool_val = ".T." in cleaned_data
            return data_type, index, (ref, bool_val)
        
        elif data_type == "FACE_OUTER_BOUND":
            # Remove all newline characters
            cleaned_data = data.replace('\n', '').replace('\r', '')
            
            # Extract the reference using regex
            ref_match = re.search(r'#(\d+)', cleaned_data)
            if not ref_match:
                return None  # or handle this case as needed
            ref = int(ref_match.group(1))
            
            bool_val = ".T." in cleaned_data
            return data_type, index, (ref, bool_val)
    
        elif data_type == "ADVANCED_FACE":
            # Remove all newline characters
            cleaned_data = data.replace('\n', '').replace('\r', '')
            
            # Find the indices of the opening and closing parentheses for the ref values
            start_idx = cleaned_data.find('(')
            end_idx = cleaned_data.find(')', start_idx)
            
            if start_idx != -1 and end_idx != -1:
                # Extract the content between the parentheses
                refs_str = cleaned_data[start_idx+1:end_idx]
                
                # Extract the individual references
                refs = tuple(map(int, re.findall(r'#(\d+)', refs_str)))
            else:
                return None  # or handle this case as needed
            
            last_ref_match = re.search(r'#(\d+)', cleaned_data[end_idx:])
            if not last_ref_match:
                return None  # or handle this case as needed
            last_ref = int(last_ref_match.group(1))
            
            bool_val = ".T." in cleaned_data[end_idx:]
            return data_type, index, (*refs, last_ref, bool_val)
    
        elif data_type == "CLOSED_SHELL":
            # Remove all newline characters
            cleaned_data = data.replace('\n', '').replace('\r', '')
            
            # Find the indices of the opening and closing parentheses
            start_idx = cleaned_data.find('(')
            end_idx = cleaned_data.rfind(')')
            
            if start_idx != -1 and end_idx != -1:
                # Extract the content between the parentheses
                refs_str = cleaned_data[start_idx+1:end_idx]
                
                # Extract the individual references
                refs = tuple(map(int, re.findall(r'#(\d+)', refs_str)))
                return data_type, index, refs
            else:
                return None
        
        elif data_type == "OPEN_SHELL":
            # Remove all newline characters
            cleaned_data = data.replace('\n', '').replace('\r', '')
            
            # Find the indices of the opening and closing parentheses
            start_idx = cleaned_data.find('(')
            end_idx = cleaned_data.rfind(')')
            
            if start_idx != -1 and end_idx != -1:
                # Extract the content between the parentheses
                refs_str = cleaned_data[start_idx+1:end_idx]
                
                # Extract the individual references
                refs = tuple(map(int, re.findall(r'#(\d+)', refs_str)))
                return data_type, index, refs
            else:
                return None
        
        elif data_type == "MANIFOLD_SOLID_BREP":
            ref_match = re.search(r'#(\d+)', data)
            if ref_match:
                ref = int(ref_match.group(1))
                return data_type, index, ref
            else:
                return None
    
        elif data_type == "SHELL_BASED_SURFACE_MODEL":
            ref_match = re.search(r'#(\d+)', data)
            if ref_match:
                ref = int(ref_match.group(1))
                return data_type, index, ref
            else:
                return None
    
        elif data_type == "ADVANCED_BREP_SHAPE_REPRESENTATION":
            components = [comp.strip() for comp in data.split(',')]
            
            # Extract all references from the components
            all_refs = re.findall(r'#(\d+)', data)
            
            if not all_refs or len(all_refs) < 2:
                return None  # Safety check in case we don't get expected data format
            
            refs = tuple(int(val) for val in all_refs[:-1])
            last_ref = int(all_refs[-1])
            
            return data_type, index, (refs, last_ref)
        
        elif data_type == "MANIFOLD_SURFACE_SHAPE_REPRESENTATION":
            components = [comp.strip() for comp in data.split(',')]
            
            # Extract all references from the components
            all_refs = re.findall(r'#(\d+)', data)
            
            if not all_refs or len(all_refs) < 2:
                return None  # Safety check in case we don't get expected data format
            
            refs = tuple(int(val) for val in all_refs[:-1])
            last_ref = int(all_refs[-1])
            
            return data_type, index, (refs, last_ref)
        
        else:   
            return None
    
    
    def parse_b_spline_surface(self,line):
        # B_SPLINE_SURFACE extraction
        pattern_bspline = r"B_SPLINE_SURFACE\s*\(\s*(\d+),\s*(\d+),"
        match_bspline = re.search(pattern_bspline, line)
        u_degree = int(match_bspline.group(1))
        v_degree = int(match_bspline.group(2))
        
        # Pattern to capture entire tuple structure
        pattern_control_points_tuple = r"\(\s*(?:#\d+\s*,\s*)*#\d+\s*\)"
            
        # Find all tuples in the line
        tuples_found = re.findall(pattern_control_points_tuple, line)
        control_points = []
        for t in tuples_found:
            # Extract all individual control points from the tuple
            control_points.append([int(x) for x in re.findall(r'#(\d+)', t)])
    
        # Pattern to capture boolean structure
        pattern_boolean_structure = r"\.UNSPECIFIED\.\s*,\s*\.(F|T)\.\s*,\s*\.(F|T)\.\s*,\s*\.(F|T)\."
        
        # Find all booleans in the line
        booleans_found = re.search(pattern_boolean_structure, line)
        
        closed_u = booleans_found.group(1) == 'T'
        closed_v = booleans_found.group(2) == 'T'
        polynomial = booleans_found.group(3) == 'T'
    
        # B_SPLINE_SURFACE_WITH_KNOTS extraction
        pattern_knots = r"B_SPLINE_SURFACE_WITH_KNOTS\s*\(\s*\(([^)]*)\s*\),\s*\(([^)]*)\s*\),\s*\(([^)]*)\s*\),\s*\(([^)]*)\s*\),\s*\.UNSPECIFIED\.\s*\)"
        match_knots = re.search(pattern_knots, line)
        if not match_knots:
            raise ValueError(f"Unexpected format for B_SPLINE_SURFACE_WITH_KNOTS in line: {line}")
    
        u_multiplicities = tuple(map(int, match_knots.group(1).split(',')))
        v_multiplicities = tuple(map(int, match_knots.group(2).split(',')))
        u_knots = tuple(map(float, match_knots.group(3).split(',')))
        v_knots = tuple(map(float, match_knots.group(4).split(',')))
    
        # RATIONAL_B_SPLINE_SURFACE extraction (optional)
        pattern_rational = r"RATIONAL_B_SPLINE_SURFACE\s*\(\s*\(\s*((\(\s*[+-]?\d+\.\d+(?:e[+-]?\d+)?\s*,\s*[+-]?\d+\.\d+(?:e[+-]?\d+)?\s*,\s*[+-]?\d+\.\d+(?:e[+-]?\d+)?\s*\)\s*,?\s*)+)\)\s*\)"
        match_rational = re.search(pattern_rational, line)
        weights = None
        if match_rational:
            # Here, the goal is to get all numbers from the nested tuples
            raw = re.findall(r'(-?\d+\.\d+e?-?\d+)', match_rational.group(1))
            n = len(control_points)
            m = len(control_points[0])
            weights = []
            for u in range(n):
                weights.append([])
                for v in range(m):
                    weights[u].append(raw[v])
    
        return {
            'type': 'BOUNDED_SURFACE B_SPLINE_SURFACE',
            'u_degree': u_degree,
            'v_degree': v_degree,
            'control_points': control_points,
            'closed_u': closed_u,
            'closed_v': closed_v,
            'polynomial': polynomial,
            'u_multiplicities': u_multiplicities,
            'v_multiplicities': v_multiplicities,
            'u_knots': u_knots,
            'v_knots': v_knots,
            'weights': weights
        }


    def parse_b_spline_curve(self,line):
        # B_SPLINE_CURVE extraction
        pattern_bspline = r"B_SPLINE_CURVE\s*\(\s*(\d+),\s*\(\s*(#\d+(?:\s*,\s*#\d+)*)\s*\)\s*,"
        match_bspline = re.search(pattern_bspline, line)
        if not match_bspline:
            raise ValueError(f"Unexpected format for B_SPLINE_CURVE in line: {line}")
    
        degree = int(match_bspline.group(1))
        control_points = tuple(map(int, re.findall(r'#(\d+)', match_bspline.group(2))))
    
        # B_SPLINE_CURVE_WITH_KNOTS extraction
        pattern_knots = r"B_SPLINE_CURVE_WITH_KNOTS\s*\(\s*\(([^)]*)\)\s*,\s*\(([^)]*)\)\s*,"
        match_knots = re.search(pattern_knots, line)
        if not match_knots:
            raise ValueError(f"Unexpected format for B_SPLINE_CURVE_WITH_KNOTS in line: {line}")
    
        multiplicities = tuple(map(int, match_knots.group(1).split(',')))
        knots = tuple(map(float, match_knots.group(2).split(',')))
    
        # RATIONAL_B_SPLINE_CURVE extraction (optional)
        pattern_rational = r"RATIONAL_B_SPLINE_CURVE\s*\(\s*\(([^)]*)\)\s*\)"
        match_rational = re.search(pattern_rational, line)
        weights = None
        if match_rational:
            weights = tuple(map(float, re.findall(r'(-?\d+\.\d+e?-?\d+)', match_rational.group(1))))
    
        return {
            'type': 'B_SPLINE_CURVE',
            'degree': degree,
            'control_points': control_points,
            'knot_multiplicities': multiplicities,
            'knot_values': knots,
            'weights': weights
        }




        



class SolFileData(object):
    '''
Class for input read from *.sol-files. This object
reads the data from the *.sol-file into a format that
is accessible to the FEModel object.
'''
    def __init__(self,inputfile):

        self.name = inputfile[:-4]
        for i in range(len(inputfile)):
            if inputfile[-i-1] == '/' or inputfile[-i-1] == '\\':
                self.name = str(inputfile[-i:-4])
                break
        self.nodes = {}
        self.nodesets = {}
        self.elements = {}
        self.elementsets = {}
        self.materials = {}
        self.meshes = {}
        self.sections = {}
        self.beamOrients = {}
        self.loads = {}
        self.boundaries = {}
        self.constraints = {}
        self.dampings = {}
        self.tables = {}
        self.solutions = {}

        read_file_start = time.time()
        print('\n\tReading input file:', self.name, '...', end=' ')
        self.input_error = self.readInput(inputfile)
        if self.input_error == False:
            read_file_stop = time.time()
            self.read_time = read_file_stop - read_file_start
            print('%.3f seconds\n' % (self.read_time))


    def readInput(self,inputfile):
        '''
    Read nodes, materials, element sections,
    elements, boundary conditions, loads and
    solutions into the FEModel object
    '''
        try:
            fobj = open(inputfile, 'r')

        except OSError as e:
            print('\n\n  *** ERROR!!!', e)

        else:
            current_solution = ''
            current_results = ''
            input_error = False
            line_number = 1

            for eachLine in fobj:

                line = [x.strip() for x in eachLine.split(',')]
                if line[0] == '':
                    pass
                elif(line[0][0] == '#'):
                    pass
                elif(line[0] == 'NODE'):
                    self.nodes[int(line[1])] = {'coord':[float(x) for x in line[2:]]}

                elif(line[0] == 'SET_NODES'):
                    tmpNodeSet = []
                    for n in line[2:]:
                        if '-' not in n:
                            tmpNodeSet.append(int(n))
                        else:
                            r_n = n.split('-')
                            for n_i in range(int(r_n[0]),int(r_n[1])):
                                tmpNodeSet.append(n_i)
                            tmpNodeSet.append(int(r_n[1]))
                    self.nodesets[int(line[1])] = tmpNodeSet

                elif(line[0] == 'MATERIAL'):
                    self.materials[line[2]] = {'type':line[1],'properties': {}}
                    if self.materials[line[2]]['type'] == 'Isotropic':
                        if len(line[3:]) >= 2:
                            self.materials[line[2]]['properties']['E-modulus'] = float(line[3])
                            self.materials[line[2]]['properties']['poisson ratio'] = float(line[4])
                        else:
                            print('\n\tERROR: (line number '+str(line_number)+')')
                            print('\tNot enough material properties specified for ', line[2])
                            input_error = True
                            break
                        if len(line[3:]) > 2:
                            self.materials[line[2]]['properties']['density'] = float(line[5])
                        if len(line[3:]) > 3:
                            self.materials[line[2]]['properties']['thermal expansion coefficient'] = float(line[6])
                        if len(line[3:]) > 4:
                            self.materials[line[2]]['properties']['conductivity'] = float(line[7])
                        if len(line[3:]) > 5:
                            self.materials[line[2]]['properties']['specific heat'] = float(line[8])
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown material type: ', line[1])
                        input_error = True
                        break

                elif(line[0] == 'SECTION'):
                    self.sections[line[2]] = {'type': line[1],
                                              'material': line[3],
                                              'properties': {}}
                    if self.sections[line[2]]['type'] == 'PlaneSect':
                        self.sections[line[2]]['properties']['thickness'] = float(line[4])
                        if len(line) > 5:
                            if line[5] == 'planestrain':
                                self.sections[line[2]]['properties']['planestrain'] = True
                            else:
                                print('\n\tERROR: (line number '+str(line_number)+')')
                                print('\tUnknown input for SECTION ', line[2], ': ', line[5])
                                input_error = True
                                break
                        else:
                            self.sections[line[2]]['properties']['planestrain'] = False
                    elif self.sections[line[2]]['type'] == 'PlateSect':
                        self.sections[line[2]]['properties']['thickness'] = float(line[4])
                    elif self.sections[line[2]]['type'] == 'RodSect':
                        self.sections[line[2]]['properties']['area'] = float(line[4])
                    elif self.sections[line[2]]['type'] == 'BeamSect':
                        self.sections[line[2]]['properties']['area'] = float(line[4])
                        self.sections[line[2]]['properties']['Izz'] = float(line[5])
                        if len(line) > 6 and line[6] != 'CrossSection':
                            self.sections[line[2]]['properties']['Iyy']  = float(line[6])
                    elif self.sections[line[2]]['type'] == 'SolidSect':
                        pass
                    elif self.sections[line[2]]['type'] == 'MassSect':
                        self.sections[line[2]]['material'] = None
                        self.sections[line[2]]['properties']['m'] = float(line[3])
                        self.sections[line[2]]['properties']['m_rx'] = float(line[4])
                        self.sections[line[2]]['properties']['m_ry'] = float(line[5])
                        self.sections[line[2]]['properties']['m_rz'] = float(line[6])
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown SECTION type: ', line[1])
                        input_error = True
                        break
                    if 'CrossSection' in line:
                        c_n = line.index('CrossSection')
                        if line[c_n+1] == 'Rectangle':
                            self.sections[line[2]]['CrossSection'] = {'Type': line[c_n+1],
                                                                      'width, w': float(line[c_n+2]),
                                                                      'height, h': float(line[c_n+3]),
                                                                      'inner width, iw': float(line[c_n+4]),
                                                                      'inner height, ih': float(line[c_n+5])}
                        elif line[c_n+1] == 'Circle':
                            self.sections[line[2]]['CrossSection'] = {'Type': line[c_n+1],
                                                                      'radius, r': float(line[c_n+2]),
                                                                      'inner radius, ir': float(line[c_n+3])}
                        elif line[c_n+1] == 'L-Beam':
                            self.sections[line[2]]['CrossSection'] = {'Type': line[c_n+1],
                                                                      'side thickness, st': float(line[c_n+2]),
                                                                      'bottom width, bw': float(line[c_n+3]),
                                                                      'bottom thickness, bt': float(line[c_n+4]),
                                                                      'height, h': float(line[c_n+5])}
                        elif line[c_n+1] == 'I-Beam':
                            self.sections[line[2]]['CrossSection'] = {'Type': line[c_n+1],
                                                                      'top width, tw': float(line[c_n+2]),
                                                                      'top thickness, tt': float(line[c_n+3]),
                                                                      'middle thickness, mt': float(line[c_n+4]),
                                                                      'bottom width, bw': float(line[c_n+5]),
                                                                      'bottom thickness, bt': float(line[c_n+6]),
                                                                      'height, h': float(line[c_n+7])}
                        elif line[c_n+1] == 'C-Beam':
                            self.sections[line[2]]['CrossSection'] = {'Type': line[c_n+1],
                                                                      'top width, tw': float(line[c_n+2]),
                                                                      'top thickness, tt': float(line[c_n+3]),
                                                                      'middle thickness, mt': float(line[c_n+4]),
                                                                      'bottom width, bw': float(line[c_n+5]),
                                                                      'bottom thickness, bt': float(line[c_n+6]),
                                                                      'height, h': float(line[c_n+7])}
                        elif line[c_n+1] == 'T-Beam':
                            self.sections[line[2]]['CrossSection'] = {'Type': line[c_n+1],
                                                                      'top width, tw': float(line[c_n+2]),
                                                                      'top thickness, tt': float(line[c_n+3]),
                                                                      'middle thickness, mt': float(line[c_n+4]),
                                                                      'height, h': float(line[c_n+5])}
                        else:
                            print('\n\tERROR: (line number '+str(line_number)+')')
                            print('\tUnknown SECTION CrossSection type: ', line[c_n+1])
                            input_error = True
                            break

                elif(line[0] == 'BEAMORIENT'):
                    self.beamOrients[line[2]] = {'type': line[1],
                                                 'elementset': int(line[3])}
                    if self.beamOrients[line[2]]['type'] == 'BEAM2N2D':
                        self.beamOrients[line[2]]['x-vec'] = [float(line[4]), float(line[5])]
                    elif self.beamOrients[line[2]]['type'] == 'BEAM2N':
                        self.beamOrients[line[2]]['x-vec'] = [float(line[4]), float(line[5]), float(line[6])]
                        self.beamOrients[line[2]]['y-vec'] = [float(line[7]), float(line[8]), float(line[9])]
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown BEAMSECTION type: ', line[1])
                        input_error = True
                        break

                elif(line[0] == 'ELEMENT'):
                    self.elements[int(line[2])] = {'type': line[1],
                                                   'section': line[3],
                                                   'nodes': [int(x) for x in line[4:]]}
                    if line[1] in ['MASS1N', 'MASS1N2D', 'ROD2N', 'ROD2N2D', 'BEAM2N', 'BEAM2N2D', 'TRI3N',
                                   'TRI6N', 'QUAD4N', 'QUAD8N', 'TET4N', 'TET10N', 'HEX8N', 'HEX20N']:
                        pass
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown ELEMENT type: ', line[1])
                        input_error = True
                        break

                elif(line[0] == 'SET_ELEMENTS'):
                    tmpElmSet = []
                    for n in line[2:]:
                        if '-' not in n:
                            tmpElmSet.append(int(n))
                        else:
                            r_n = n.split('-')
                            for n_i in range(int(r_n[0]),int(r_n[1])):
                                tmpElmSet.append(n_i)
                            tmpElmSet.append(int(r_n[1]))
                    self.elementsets[int(line[1])] = tmpElmSet

                elif(line[0] == 'LOAD'):
                    self.loads[line[2]] = {'type': line[1],
                                           'vector': [float(x) for x in line[5:]]}
                    if self.loads[line[2]]['type'] == 'Gravity':
                        self.loads[line[2]]['acceleration'] = float(line[4])
                        self.loads[line[2]]['elementset'] = int(line[3])
                    elif self.loads[line[2]]['type'] in ['ForceConcentrated', 'Force', 'ForceDynamic', 'Acceleration']:
                        self.loads[line[2]]['force'] = float(line[4])
                        self.loads[line[2]]['nodeset'] = int(line[3])
                    elif self.loads[line[2]]['type'] == 'Torque':
                        self.loads[line[2]]['torque'] = float(line[4])
                        self.loads[line[2]]['nodeset'] = int(line[3])
                    elif self.loads[line[2]]['type'] == 'Pressure':
                        self.loads[line[2]]['pressure'] = float(line[4])
                        self.loads[line[2]]['nodeset'] = int(line[3])
                    elif self.loads[line[2]]['type'] == 'ForceDistributed':
                        self.loads[line[2]]['force'] = float(line[4])
                        self.loads[line[2]]['elementset'] = int(line[3])
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown LOAD type: ', line[1])
                        input_error = True
                        break

                elif(line[0] == 'BOUNDARY'):
                    self.boundaries[line[2]] = {'type': line[1],
                                                'nodeset': int(line[3]),
                                                'value': float(line[4]),
                                                'DOFs': [int(x) for x in line[5:]]}
                    if line[1] in ['Displacement']:
                        pass
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown BOUNDARY type: ', line[1])
                        input_error = True
                        break

                elif(line[0] == 'CONSTRAINT'):
                    self.constraints[line[2]] = {'type': line[1],
                                                 'nodeset1': int(line[3]),
                                                 'nodeset2': int(line[4])}
                    if line[1] == 'TouchLock':
                        self.constraints[line[2]]['tolerance'] = float(line[5])
                        self.constraints[line[2]]['DOFs'] = [int(x) for x in line[6:]]
                    else:
                        self.constraints[line[2]]['DOFs'] = [int(x) for x in line[6:]]
                    if line[1] in ['NodeLock', 'TouchLock']:
                        pass
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown CONSTRAINT type: ', line[0])
                        input_error = True
                        break

                elif(line[0] == 'DAMPING'):
                    self.dampings[line[2]] = {'type': line[1]}
                    if line[1] == 'Viscous':
                        self.dampings[line[2]]['damping_ratio'] = float(line[3])
                    elif line[1] == 'Frequency':
                        self.dampings[line[2]]['damping_ratio'] = 1.
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown DAMPING type: ', line[1])
                        input_error = True
                        break

                elif(line[0] == 'TABLE'):
                    self.tables[line[2]] = {'type': line[1],
                                            'filename': line[4]}
                    if line[1] == 'Acceleration':
                        self.tables[line[2]]['type'] = 'AccelTable'
                        self.tables[line[2]]['boundary'] = line[3]
                    elif line[1] == 'StressStrain':
                        self.tables[line[2]]['type'] = 'StressStrainTable'
                        self.tables[line[2]]['material'] = line[3]
                    elif line[1] == 'ForceDynamic':
                        self.tables[line[2]]['type'] = 'ForceTable'
                        self.tables[line[2]]['load'] = line[3]
                    elif line[1] == 'DampingRatio':
                        self.tables[line[2]]['type'] = 'DampingTable'
                        self.tables[line[2]]['damping'] = line[3]
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown TABLE type: ', line[1])
                        input_error = True
                        break

                elif(line[0] == 'SOLUTION'):
                    self.solutions[line[1]] = {'type':line[2],
                                               'meshes': {},
                                               'constraints': [],
                                               'loads': [],
                                               'boundaries': [],
                                               'dampings': [],
                                               'results': {} }
                    current_solution = line[1]
                    self.meshes[current_solution] = 'all'
                    if line[2] in ['Static', 'StaticPlastic', 'Eigenmodes', 'ModalDynamic']:
                        pass
                    else:
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tUnknown SOLUTION type: ', line[2])
                        input_error = True
                        break

                elif(line[0] == 'MESHES'):
                    if current_solution == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tMESHES need to be specified AFTER the SOLUTION they are used in!')
                        input_error = True
                        break
                    else:
                        self.solutions[current_solution]['meshes'] = [int(x) for x in line[1:]]
                        self.meshes[current_solution] = [int(x) for x in line[1:]]

                elif(line[0] == 'CONSTRAINTS'):
                    if current_solution == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tCONSTRAINTS need to be specified AFTER the SOLUTION they are used in!')
                        input_error = True
                        break
                    else:
                        self.solutions[current_solution]['constraints'].append(line[1])

                elif(line[0] == 'LOADS'):
                    if current_solution == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tLOADS need to be specified AFTER the SOLUTION they are used in!')
                        input_error = True
                        break
                    else:
                        self.solutions[current_solution]['loads'].append(line[1])

                elif(line[0] == 'BOUNDARIES'):
                    if current_solution == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tBOUNDARIES need to be specified AFTER the SOLUTION they are used in!')
                        input_error = True
                        break
                    else:
                        self.solutions[current_solution]['boundaries'].append(line[1])

                elif(line[0] == 'DAMPINGS'):
                    if current_solution == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tDAMPINGS need to be specified AFTER the SOLUTION they are used in!')
                        input_error = True
                        break
                    else:
                        self.solutions[current_solution]['dampings'].append(line[1])

                elif(line[0] == 'RESULTS'):
                    if current_solution == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tRESULTS need to be specified AFTER the SOLUTION they are used in!')
                        input_error = True
                        break
                    else:
                        self.solutions[current_solution]['results'] = {}
                        current_results = line[1]

                elif(line[0] == 'DISPLACEMENT'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tDISPLACEMENT needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['displacement'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['displacement']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['displacement'] = {}
                            self.solutions[current_solution]['results']['displacement']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['displacement']['result DOF'] = int(line[line.index('text')+2])

                elif(line[0] == 'ACCELERATION'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tACCELERATION needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['acceleration'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['acceleration']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['acceleration'] = {}
                            self.solutions[current_solution]['results']['acceleration']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['acceleration']['result DOF'] = int(line[line.index('text')+2])

                elif(line[0] == 'FRF_ACCEL'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tFRF_ACCEL needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['frf_accel'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['frf_accel']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['frf_accel'] = {}
                            self.solutions[current_solution]['results']['frf_accel']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['frf_accel']['result DOF'] = int(line[line.index('text')+2])

                elif(line[0] == 'SRS_ACCEL'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tFRF_ACCEL needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['srs_accel'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['srs_accel']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['srs_accel'] = {}
                            self.solutions[current_solution]['results']['srs_accel']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['srs_accel']['result DOF'] = int(line[line.index('text')+2])

                elif(line[0] == 'VELOCITY'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tVELOCITY needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['velocity'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['velocity']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['velocity'] = {}
                            self.solutions[current_solution]['results']['velocity']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['velocity']['result DOF'] = int(line[line.index('text')+2])

                elif(line[0] == 'NODEFORCE'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tNODEFORCE needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['nodeforce'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['nodeforce']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['nodeforce'] = {}
                            self.solutions[current_solution]['results']['nodeforce']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['nodeforce']['result DOF'] = int(line[line.index('text')+2])

                elif(line[0] == 'ELEMENTFORCE'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tELEMENTFORCE needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['elementforce'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['elementforce']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['elementforce'] = {}
                            self.solutions[current_solution]['results']['elementforce']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['elementforce']['result DOF'] = int(line[line.index('text')+2])


                elif(line[0] == 'STRESS'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tSTRESS needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['stress'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['stress']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['stress'] = {}
                            self.solutions[current_solution]['results']['stress']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['stress']['result DOF'] = int(line[line.index('text')+2])

                elif(line[0] == 'STRAIN'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tSTRAIN needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        line += ['pass', 'pass']
                        if 'plot' in line:
                            self.solutions[current_solution]['results']['strain'] = \
                                        {'plot': int(line[line.index('plot')+1])}
                            if line[line.index('plot')+2].isdigit():
                                self.solutions[current_solution]['results']['strain']['result DOF'] = int(line[line.index('plot')+2])
                        if 'text' in line:
                            if 'plot' not in line:
                                self.solutions[current_solution]['results']['strain'] = {}
                            self.solutions[current_solution]['results']['strain']['text'] = int(line[line.index('text')+1])
                            if line[line.index('text')+2].isdigit():
                                self.solutions[current_solution]['results']['strain']['result DOF'] = int(line[line.index('text')+2])

                elif(line[0] == 'MODESHAPES'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tMODESHAPES need to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        self.solutions[current_solution]['results']['modeshapes'] = int(line[1])

                elif(eachLine[0:14] == '\tENERGYDENSITY'):
                    if current_results == '':
                        print('\n\tERROR: (line number '+str(line_number)+')')
                        print('\tENERGYDENSITY needs to be specified AFTER the RESULTS they are from!')
                        input_error = True
                        break
                    else:
                        self.solutions[current_solution]['results']['energydensity'] = True

                else:
                    print('\n\tINPUT WARNING: (line number '+str(line_number)+')')
                    print('\tUnknown input '+line[0]+'...   Ignored!')

                line_number +=1





            # reality check for sol-file
        # ---------------------------------


            # check that there are nodes
            if len(self.nodes) == 0:
                print('\n\tERROR:\n\tNo nodes have been defined.')
                input_error = True
        
            # check that there are elements
            if len(self.elements) == 0:
                print('\n\tERROR:\n\tNo elements have been defined.')
                input_error = True

            # check if all elements and nodes specified in
            # element- and nodesets are actually defined
            for nodeset in self.nodesets:
                if not all(node in self.nodes for node in self.nodesets[nodeset]):
                    print('\n\tERROR:\n\tNodeset '+str(nodeset)+' contains nodes that have not been defined.')
                    input_error = True

            for elementset in self.elementsets:
                if not all(element in self.elements for element in self.elementsets[elementset]):
                    print('\n\tERROR:\n\tElementset '+str(elementset)+' contains elements that have not been defined.')
                    input_error = True

            for element in self.elements:
                if self.elements[element]['section'] not in self.sections:
                    print('\n\tERROR:\n\tElement '+str(element)+' has section which has not been defined.')
                    input_error = True

                # check that elements have the right
                # number of nodes
                if self.elements[element]['type'] in ['ROD2N2D', 'ROD2N', 'BEAM2N2D', 'BEAM2N']:
                    if len(self.elements[element]['nodes']) != 2:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] == 'TRI3N':
                    if len(self.elements[element]['nodes']) != 3:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] == 'TRI6N':
                    if len(self.elements[element]['nodes']) != 6:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] == 'QUAD4N':
                    if len(self.elements[element]['nodes']) != 4:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] == 'QUAD8N':
                    if len(self.elements[element]['nodes']) != 8:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] == 'TET4N':
                    if len(self.elements[element]['nodes']) != 4:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] == 'TET10N':
                    if len(self.elements[element]['nodes']) != 10:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] == 'HEX8N':
                    if len(self.elements[element]['nodes']) != 8:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] == 'HEX20N':
                    if len(self.elements[element]['nodes']) != 20:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                elif self.elements[element]['type'] in ['MASS1N', 'MASS1N2D']:
                    if len(self.elements[element]['nodes']) != 1:
                        print('\n\tERROR:\n\tElement '+str(element)+' does not have the right number of nodes.')
                        input_error = True
                else:
                    pass

            # check if specified section
            # material is actually defined
            for sect in self.sections:
                if self.sections[sect]['type'] != 'MassSect':
                    if self.sections[sect]['material'] not in self.materials:
                        print('\n\tERROR:\n\tSection '+sect+' uses material that has not been defined.')
                        input_error = True

            # check if specifed tables
            # can be accessed
            for table in self.tables:
                if not os.path.isfile(self.tables[table]['filename']):
                    print('\n\tERROR:\n\tTable '+table+' uses file '+self.tables[table]['filename']+' which does not exist.')
                    input_error = True
            
            # check that nodelock constraints have
            # nodeset with only one node in it
            for constraint in self.constraints:
                if self.constraints[constraint]['nodeset1'] not in self.nodesets or \
                        self.constraints[constraint]['nodeset2'] not in self.nodesets:
                    print('\n\tERROR:\n\tConstraint '+constraint+' has nodeset that is not defined.')
                    input_error = True
                    break
                if self.constraints[constraint]['type'] == 'NodeLock':
                    not_ok = True
                    if len(self.nodesets[self.constraints[constraint]['nodeset1']]) == 1:
                        not_ok = False
                    if len(self.nodesets[self.constraints[constraint]['nodeset2']]) == 1:
                        not_ok = False
                    if not_ok:
                        print('\n\tERROR:\n\tConstraint '+constraint+' must have at least one nodeset with only one node in it.')
                        input_error = True

            for sol in self.solutions:
                # check that results are specified for
                # every solution
                if len(self.solutions[sol]['results']) == 0:
                    print('\n\tERROR:\n\tSolution '+sol+' has no results requested.')
                    input_error = True
                    break

                # check if constraints in solution
                # have been defined
                for constraint in self.solutions[sol]['constraints']:
                    if constraint not in self.constraints:
                        print('\n\tERROR:\n\tSolution '+sol+' has constraints that are not defined.')
                        input_error = True
                    else:
                        pass

                # check if specified loads are defined
                for load in self.solutions[sol]['loads']:
                    if load not in self.loads:
                        print('\n\tERROR:\n\tSolution '+sol+' has loads that are not defined.')
                        input_error = True
                        break
                    else:
                        if 'nodeset' in self.loads[load]:
                            if self.loads[load]['nodeset'] not in self.nodesets:
                                print('\n\tERROR:\n\tLoad '+load+' has nodeset that is not defined.')
                                input_error = True
                        if 'elementset' in self.loads[load]:
                            if self.loads[load]['elementset'] not in self.elementsets:
                                print('\n\tERROR:\n\tLoad '+load+' has elementset that is not defined.')
                                input_error = True

                if self.solutions[sol]['type'] in ['Static', 'StaticPlastic']:
                    # check that results requested are
                    # supported for solution type
                    for result in self.solutions[sol]['results']:
                        if result not in ['displacement', 'nodeforce', 'elementforce', 'stress', 'strain']:
                            print('\n\tERROR:\n\t'+result+' not supported for solution type '+self.solutions[sol]['type'])
                            input_error = True

                    # check if static solutions are
                    # constrained with boundary conditions
                    if len(self.solutions[sol]['boundaries']) == 0:
                        print('\n\tERROR:\n\tSolution '+sol+' has no boundary conditions applied.')
                        input_error = True

                    # check if boundary conditions in solution
                    # have been defined
                    else:
                        for boundary in self.solutions[sol]['boundaries']:
                            if boundary not in self.boundaries:
                                print('\n\tERROR:\n\tSolution '+sol+' has boundary conditions that are not defined.')
                                input_error = True
                            elif self.boundaries[boundary]['nodeset'] not in self.nodesets:
                                print('\n\tERROR:\n\tBoundary '+boundary+' has nodeset that is not defined.')
                                input_error = True
                            else:
                                pass

                    # check if specified loads are
                    # applicable to solution type
                    for load in self.solutions[sol]['loads']:
                        if load in self.loads:
                            if self.loads[load]['type'] not in ['Force', 'ForceConcentrated', 'Torque', 'ForceDistributed', 'Gravity', 'Pressure']:
                                print('\n\tERROR:\n\tSolution '+sol+' has loads that can not be used in Static solution.')
                                input_error = True

                elif self.solutions[sol]['type'] in ['Eigenmodes']:
                    if len(self.solutions[sol]['results']) == 0:
                        print('\n\tERROR:\n\tSolution '+sol+' has no results requested.')
                        input_error = True

                elif self.solutions[sol]['type'] in ['ModalDynamic']:
                    # check that results requested are
                    # supported for solution type
                    for result in self.solutions[sol]['results']:
                        if result not in ['displacement', 'velocity', 'acceleration', 'frf_accel', 'srs_accel', 'modeshapes']:
                            print('\n\tERROR:\n\t'+result+' not supported for solution type '+self.solutions[sol]['type'])
                            input_error = True

                    # check if boundary conditions in solution
                    # have been defined
                    else:
                        for boundary in self.solutions[sol]['boundaries']:
                            if boundary not in self.boundaries:
                                print('\n\tERROR:\n\tSolution '+sol+' has boundary conditions that are not defined.')
                                input_error = True
                            elif self.boundaries[boundary]['nodeset'] not in self.nodesets:
                                print('\n\tERROR:\n\tBoundary '+boundary+' has nodeset that is not defined.')
                                input_error = True
                            else:
                                pass

                    # check if specified loads are
                    # applicable to solution type
                    for load in self.solutions[sol]['loads']:
                        if self.loads[load]['type'] not in ['ForceDynamic', 'Acceleration']:
                            print('\n\tERROR:\n\tSolution '+sol+' has loads that can not be used in ModalDynamic solution.')
                            input_error = True

                else:
                    pass


                # check that damping is applied to
                # modal dynamics solution
                

                
            fobj.close()
            return input_error





if __name__ == '__main__':
    stpf = StepFileData('test_part3_AP214.step')
        