import bpy
from mathutils.interpolate import poly_3d_calc
from mathutils import Vector
import bmesh
import numpy as np


###############################################
#       Correspondence Computation for        #
# Characters with Different Mesh Connectivity #
###############################################

# We assume Mesh Segmentation Results are given in advance
# load mesh segmentation results for source and target models
labels_female = np.load("C:\\Users\\Spock_the_Wizard\\Desktop\\researchSpring21\\clothes\\labels.npy")
labels_male = np.load('C:\\Users\\Spock_the_Wizard\\Desktop\\researchSpring21\\clothes\\labels_male.npy')

female = bpy.data.objects['source_f']
male = bpy.data.objects['target_m']

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
        #exit(2)
        
    return coords

# Idx 0 is the HEAD (takes up over 3,000 vertices, which is pretty wasteful)
# Our code works for garments that go below the head
for i in range(1,16):
    res = get_correspondence(i,segments_f, segments_m)
    for v in res:
        correspondences.append(v)
        

###############################################
#               Deformation for               #
# Characters with Different Mesh Connectivity #
###############################################

source = bpy.data.objects['source_f']
garment = bpy.data.objects['garment_skirt']
target = bpy.data.objects['target_m']

garment_vertices = garment.data.vertices
source_vertices = source.data.vertices
target_vertices = target.data.vertices

garment_polygons = garment.data.polygons
source_polygons = source.data.polygons

deformed_garment_vertices = []

def get_corresponding_vertex(source_idx):
    for idx,coords in correspondences:
        if source_idx==idx:
            return coords
    print(source_idx)
    return None

for v in garment_vertices:
    # get projection from garment vertex to source model
    res, loc, nor, idx = source.closest_point_on_mesh(v.co)
    
    face_vertex_index = source_polygons[idx].vertices
    # interpolate projected point from source onto target
    source_face_vertices = [source_vertices[vx].co for vx in face_vertex_index]
    bary_weights = poly_3d_calc(source_face_vertices, loc)
    target_projection_point = Vector((0,0,0))

    # get correspondence from source to target    
    for i in range(3):
        co = get_corresponding_vertex(face_vertex_index[i])
        if co==None:
            print('correspondence doesnt exist')
            #exit(-1)
            #return
        else:
            target_projection_point += bary_weights[i] * co
        
    scale = 1.0
    deformed_garment_vertices.append(scale * (target_projection_point + (v.co - loc)))
  

mesh = bpy.data.meshes.new("mesh")
new_garment = bpy.data.objects.new("deformed_garment", mesh)

col = bpy.data.collections.get("Collection")
col.objects.link(new_garment)
bpy.context.view_layer.objects.active = new_garment

mesh = bpy.context.object.data
bm = bmesh.new()

face_list = []
for p in garment_polygons:
    idx_list = [ v for v in p.vertices]
    face_list.append(idx_list)

mesh.from_pydata(deformed_garment_vertices, [], face_list)