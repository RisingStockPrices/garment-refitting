import bpy
import numpy as np
from mathutils import Vector
import os


###############################################
#       Correspondence Computation for        #
# Characters with Different Mesh Connectivity #
###############################################


# We assume Mesh Segmentation Results are given in advance
# load mesh segmentation results for source and target models
labels_female = np.load("C:\\Users\\Spock_the_Wizard\\Desktop\\researchSpring21\\clothes\\labels.npy")
labels_male = np.load('C:\\Users\\Spock_the_Wizard\\Desktop\\researchSpring21\\clothes\\labels_male.npy')

female = bpy.data.objects['F']
male = bpy.data.objects['M']

vertices_f = female.data.vertices
vertices_m = male.data.vertices

segments_f = []
segments_m = []

# This is what we're here for
correspondences = []

for i in range(16):
    segments_f.append([])
    segments_m.append([])
for v in vertices_m:
    seg_idx = labels_male[v.index]
    segments_m[seg_idx].append(v)
for v in vertices_f:
    seg_idx = labels_female[v.index]
    segments_f[seg_idx].append(v)
    

def center_of_mass(vertices):
    center= Vector.Fill(3)
    for v in vertices:
        center+= v.co
    center /= len(vertices)
    return center

def translated_coords(vertices, trnsl_vector):
    coords = []
    for v in vertices:
        coords.append(v.co+trnsl_vector)
    return coords

def get_coords (vertices):
    coords = []
    for v in vertices:
        coords.append(v.co)
    return coords

# Creates temporary object for visualization
def make_temp_mesh(vertices):
    obj_name = "temp_visualization"
    mesh_data = bpy.data.meshes.new(obj_name + "_data")
    obj = bpy.data.objects.new(obj_name, mesh_data)
    bpy.context.collection.objects.link(obj)
    mesh_data.from_pydata(vertices, [], [])

# Checks if ray cast result is in current segment
def is_current_segment(face_idx, target_obj, seg_vertices):
    for idx in target_obj.data.polygons[face_idx].vertices:
        for v in seg_vertices:
            if v.index==idx:
                return True
    return False
    
# takes in segments_f or segments_m as source and target,
def get_correspondence(seg_idx, _source, _target):
    
    coords = []
    post_process_coords = []
    source = _source[seg_idx]
    target = _target[seg_idx]
    
    center_source = center_of_mass(source)
    center_target = center_of_mass(target)
    
    # translate source to target's center of mass
    translated_source = translated_coords(source, (center_target-center_source))
    sum = 0
    # attempt ray cast for all source points to target object
    # get an average of the distance ratio between original vector and resulting vector
    # for ones that don't have a hit -> apply scaling afterwards (post-processing)
    for v in source:
        direction = v.co-center_source
        res, loc, nor, idx = male.ray_cast(center_target,direction)
        if res and is_current_segment(idx, male, target):
            coords.append((v.index,loc))
            sum += (loc-center_target).length/direction.length
        else:
            post_process_coords.append((v.index,direction))
        
    average_scale = sum / len(coords)
    print(average_scale)
    for idx,vec in post_process_coords:
        coords.append((idx,average_scale*vec+center_target))
    
    if len(source)!=len(coords):
        print('Error! different length in GET_CORRESPONDENCE')
        exit(2)
        
    return coords

# Idx 0 is the HEAD (takes up over 3,000 vertices, which is pretty wasteful)
# Our code works for garments that go below the head
for i in range(1,16):
    res = get_correspondence(i,segments_f, segments_m)
    for v in res:
        correspondences.append(v)
        
print(correspondences)