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
#labels_female = np.load("C:\\Users\\Spock_the_Wizard\\Desktop\\researchSpring21\\clothes\\labels.npy")
#labels_male = np.load('C:\\Users\\Spock_the_Wizard\\Desktop\\researchSpring21\\clothes\\labels_male.npy')
labels_female = np.load("C:\\Users\\SpocktheWizard\\Desktop\\garment-refitting\\clothes\\labels.npy")
labels_male = np.load('C:\\Users\\SpocktheWizard\\Desktop\\garment-refitting\\clothes\\labels_male.npy')

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
    """
    if 1<= seg_idx <=2:
        seg_idx = 1
    elif 4<=seg_idx<=5:
        seg_idx = 4
    elif 7<=seg_idx<=8:
        seg_idx = 7
    elif 10<=seg_idx<=11:
        seg_idx = 10
    elif 13<=seg_idx<=14:
        seg_idx = 13
    """
    segments_m[seg_idx].append(v)
    
for v in vertices_f:
    seg_idx = labels_female[v.index]
    """
    if 1<= seg_idx <=2:
        seg_idx = 1
    elif 4<=seg_idx<=5:
        seg_idx = 4
    elif 7<=seg_idx<=8:
        seg_idx = 7
    elif 10<=seg_idx<=11:
        seg_idx = 10
    elif 13<=seg_idx<=14:
        seg_idx = 13
    """
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
def is_current_segment(face_idx, target_obj, seg_vertices, strict=True):
    if strict:
        # for the three vertices making up the polygon,
        # check if each one is a valid member of the segment
        for idx in target_obj.data.polygons[face_idx].vertices:
            okay = False
            for v in seg_vertices:
                if v.index==idx:
                    # found a match
                    okay = True
                    break
            if not okay:
                return False
        return True
    else:
        # looser restrictions
        # at least one vertex should be a valid member
        for idx in target_obj.data.polygons[face_idx].vertices:
            for v in seg_vertices:
                if v.index==idx:
                    return True
        return False

def bbox(vertices):
    x_min = float("inf")
    x_max = float("-inf")
    y_min = float("inf")
    y_max = float("-inf")
    z_min = float("inf")
    z_max = float("-inf")
    
    for v in vertices:
        if v[0]>x_max:
            x_max = v[0]
        if v[0]<x_min:
            x_min = v[0]
        if v[1]>y_max:
            y_max = v[1]
        if v[1]<y_min:
            y_min = v[1]
        if v[2]>z_max:
            z_max = v[2]
        if v[2]<z_min:
            z_min = v[2]
    
    center = Vector((x_min+x_max,y_min+y_max,z_min+z_max))
    center /=2
    
    return center, x_max-x_min, y_max-y_min, z_max-z_min
        
def get_nearest_point(pt, points):
    dist = []
    for v in points:
        distance = (pt-v.co).length
        dist.append((v, distance))
    
    dist.sort(key=lambda x:x[1])
    return dist[0][0]

# takes in segments_f or segments_m as source and target,
def get_correspondence(seg_idx, _source, _target):
    
    coords = []
    post_process_coords = []
    source = _source[seg_idx]
    target = _target[seg_idx]
    
    if len(source)==0:
        return coords
    
    # compute bounding box for source and target segments
    center_source, x_s, y_s, z_s = bbox(get_coords(source))
    center_target, x_t, y_t, z_t = bbox(get_coords(target))
    #print(bbox_source)
    
    for v in source:
        direction = v.co - center_source
        
        # scale up the direction
        direction[0] *= x_t/x_s
        direction[1] *= y_t/y_s
        direction[2] *= z_t/z_s

        
        res, loc, nor, idx = male.ray_cast(center_target,direction)
        if res is True and is_current_segment(idx,male,target):    
            coords.append((v.index,loc))
        else:
            coords.append((v.index,center_target+direction))
                # update direction to better fit the segment
        
        #coords.append((v.index,loc))
    
    """
    
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
        if res and is_current_segment(idx, male, target,strict=False):
            coords.append((v.index,loc))
            sum += (loc-center_target).length/direction.length
        else:
            clo = get_nearest_point(center_target+direction,target)
            coords.append((v.index,clo.co))
            #print('closest point',clo)
            #r, l, n, i = male.closest_point_on_mesh(center_target+direction)
            #coords.append((v.index, l))
            #print('undealt idx:',v.index)
    
    # make_temp_mesh(list(map(lambda x:x[1],coords)))
    """
    if len(source)!=len(coords):
        print('Segment [',seg_idx,']: Successfully recovered',len(coords),'vertices out of',len(source))
    
    return coords

#get_correspondence(15, segments_f, segments_m)


# Idx 0 is the HEAD (takes up over 3,000 vertices, which is pretty wasteful)
# Our code works for garments that go below the head
for i in range(16):
    res = get_correspondence(i,segments_f, segments_m)
    for v in res:
        correspondences.append(v)
        

###############################################
#               Deformation for               #
# Characters with Different Mesh Connectivity #
###############################################

source = bpy.data.objects['source_f']
garment = bpy.data.objects['garment_pants']
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
    #print(source_idx)
    return None

valid_idx = []
for v in garment_vertices:
    # bind garment vertex to closets point in source mesh
    res, loc, nor, idx = source.closest_point_on_mesh(v.co)
    
    # compute barycentric coordinates of the closest point
    face_vertex_index = source_polygons[idx].vertices
    source_face_vertices = [source_vertices[vx].co for vx in face_vertex_index]
    bary_weights = poly_3d_calc(source_face_vertices, loc)
    target_projection_point = Vector((0,0,0))

    # for each vertex in the entailing triangle, 
    # get correspondending point in target    
    check = False
    for i in range(3):
        co = get_corresponding_vertex(face_vertex_index[i])
        if co==None:
            #print('correspondence doesnt exist for segment')
            check = True
            break
        else:
            target_projection_point += bary_weights[i] * co
    if check:
        continue
    scale = 1.5
    deformed_garment_vertices.append( (target_projection_point + scale*(v.co - loc)))
    valid_idx.append(v.index)

mesh = bpy.data.meshes.new("mesh")
new_garment = bpy.data.objects.new("deformed_garment", mesh)

col = bpy.data.collections.get("Collection")
col.objects.link(new_garment)
bpy.context.view_layer.objects.active = new_garment

mesh = bpy.context.object.data
bm = bmesh.new()

face_list = []
for p in garment_polygons:
    skip = False
    idx_list = []
    for v in p.vertices:
        try:
            # get new index of the vertex in deformed vertices list
            new_idx = valid_idx.index(v)
            idx_list.append(new_idx)
        except ValueError:
            # move on to the next polygon
            break
    # valid polygon in the deformed garment
    if len(idx_list) == len(p.vertices):
        face_list.append(idx_list)
        
    #idx_list = [ v for v in p.vertices]
    #face_list.append(idx_list)
print('successfully recovered ',len(face_list),'out of',len(garment_polygons))
mesh.from_pydata(deformed_garment_vertices, [], face_list)
