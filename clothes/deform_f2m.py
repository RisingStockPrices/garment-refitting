import bpy
from mathutils.interpolate import poly_3d_calc
from mathutils import Vector
import mathutils
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

# evaluation criteria
unmatched_correspondences = 0

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
def is_current_segment(face_idx, target_obj, seg_vertices, strict=False):
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

def visualize_bbox(center, x, y, z):
    obj_name = "bbox"
    mesh_data = bpy.data.meshes.new(obj_name + "_data")
    obj = bpy.data.objects.new(obj_name, mesh_data)
    bpy.context.collection.objects.link(obj)
    
    vertices = []
    vertices.append(center+Vector((x,y,z))*0.5)
    vertices.append(center+Vector((x,y,-z))*0.5)
    vertices.append(center+Vector((x,-y,z))*0.5)
    vertices.append(center+Vector((x,-y,-z))*0.5)
    vertices.append(center+Vector((-x,y,z))*0.5)
    vertices.append(center+Vector((-x,y,-z))*0.5)
    vertices.append(center+Vector((-x,-y,z))*0.5)
    vertices.append(center+Vector((-x,-y,-z))*0.5)
    
    face_list = []
    face_list.append([0,1,3,2])
    face_list.append([0,1,5,4])
    face_list.append([1,3,7,5])
    face_list.append([4,5,7,6])
    face_list.append([2,3,7,6])
    face_list.append([0,2,6,4])
    mesh_data.from_pydata(vertices, [], face_list)
    
def compute_human_grid(human):
    x_rightarm_torso = -0.2
    x_leftarm_torso = 0.2
    z_torso_legs = 0.69
    x_legs = 0
    
    # order : torso-leftarm-rightarm-leftleg-rightleg
    bboxs = []
    i = 1
    while i<16:
        vertices = []
        for idx in range(i, i+3):
            vertices = vertices+human[idx]
        bboxs.append(list(bbox(get_coords(vertices))))
        i+=3
    
    # adjust bounding boxes so that they don't overlap
    [center_torso, x_torso, y_torso, z_torso] = bboxs[0]
    for i in range(1,5):
        [center,x,y,z] = bboxs[i]
        # left arm
        if i==1:
            bound_r = center_torso[0]+0.5*x_torso
            bound_l = center[0]+0.5*x
            bboxs[i][0][0] = 0.5*(bound_l+bound_r)
            bboxs[i][1] = abs(bound_l-bound_r)
        # right arm
        elif i==2:
            bound_l = center_torso[0]-0.5*x_torso
            bound_r = center[0]-0.5*x
            bboxs[i][0][0] = 0.5*(bound_l+bound_r)
            bboxs[i][1] = abs(bound_l-bound_r)
        else:
            bound_u = center_torso[2]-0.5*z_torso
            bound_d = center[2]-0.5*z
            bboxs[i][0][2] = 0.5*(bound_u+bound_d)
            bboxs[i][3] = abs(bound_u-bound_d)
        
    for bx in bboxs:
        visualize_bbox(bx[0],bx[1],bx[2],bx[3])
    
    return bboxs

def create_plane(value, point, axis='x'):
    size = 1
    vertices = []
    obj_name = "plane"
    mesh_data = bpy.data.meshes.new(obj_name + "_data")
    obj = bpy.data.objects.new(obj_name, mesh_data)
    bpy.context.collection.objects.link(obj)
    
    if axis=='x':
        ref = Vector((value,point[1],point[2]))
        vertices.append([ref[0],ref[1]+size,ref[2]+size])
        vertices.append([ref[0],ref[1]+size,ref[2]-size])
        vertices.append([ref[0],ref[1]-size,ref[2]+size])
        vertices.append([ref[0],ref[1]-size,ref[2]-size])
        mesh_data.from_pydata(vertices, [], [[0,1,3,2]])
    elif axis=='z':
        ref = Vector((point[0],point[1],value))
        vertices.append([ref[0]+size,ref[1]+size,ref[2]])
        vertices.append([ref[0]+size,ref[1]-size,ref[2]])
        vertices.append([ref[0]-size,ref[1]+size,ref[2]])
        vertices.append([ref[0]-size,ref[1]-size,ref[2]])
        mesh_data.from_pydata(vertices, [], [[0,1,3,2]])
        
    
TORSO = 0
LEFTARM = 1
RIGHTARM = 2
LEFTLEG = 3
RIGHTLEG = 4
ETC = -1

CENTER = 0
X_LEN = 1
Y_LEN = 2
Z_LEN = 3

bboxs_target = compute_human_grid(segments_m)
bboxs_source = compute_human_grid(segments_f)

source_torso_right = bboxs_source[0][0][0]-0.5*bboxs_source[0][1]
source_torso_left = bboxs_source[0][0][0]+0.5*bboxs_source[0][1]
source_torso_bottom = bboxs_source[0][0][2]-0.5*bboxs_source[0][3]
source_torso_up = bboxs_source[0][0][2]+0.5*bboxs_source[0][3]

target_torso_right = bboxs_target[0][0][0]-0.5*bboxs_target[0][1]
target_torso_left = bboxs_target[0][0][0]+0.5*bboxs_target[0][1]
target_chest_bottom = bboxs_target[LEFTARM][CENTER][2]-0.5*bboxs_target[LEFTARM][Z_LEN]
# forgot to de-separate the legs.... f...

def get_source_cross_section_vector(v, horizontal=True, side='None'):
    # get cross section in source
    for f in bm_source.faces:
        f.select = True
    edges = [e for e in bm_source.edges]
    faces = [f for f in bm_source.faces]
    geom = []
    geom.extend(edges)
    geom.extend(faces)
    
    if horizontal:
        no = [0,0,1]
    else:
        no = [1,0,0]
        
    # get cross section curve
    result = bmesh.ops.bisect_plane(bm_source,
                              dist=0.001,
                              geom=geom,
                              plane_co=v.co,
                              plane_no=no)
                              
    geom_cut = result['geom_cut']
    count=0
    center_mass = Vector((0,0,0))
    if side=='None':
        for item in geom_cut:
            if isinstance(item, bmesh.types.BMVert):
                if source_torso_left>=item.co[0]>=source_torso_right:
                    center_mass = center_mass + item.co
                    count+=1
        center_mass /= count
    else:
        for item in geom_cut:
            if isinstance(item, bmesh.types.BMVert):
                if (side=='left' and item.co[0]>0) or (side=='right' and item.co[0]<0):
                    center_mass = center_mass + item.co
                    count+=1
        center_mass /= count
    
    #make_temp_mesh([center_mass, v.co])
    #print(geom_cut)
    
    return v.co-center_mass

def get_target_chest_correspondence(z_value, x_value,v_src):
    
    for f in bm_target.faces:
        f.select = True
    edges = [e for e in bm_target.edges]
    faces = [f for f in bm_target.faces]
    geom = []
    geom.extend(edges)
    geom.extend(faces)
        
    # get cross section curve
    result = bmesh.ops.bisect_plane(bm_target,
                              dist=0.001,
                              geom=geom,
                              plane_co=[0,0,z_value],
                              plane_no=[0,0,1])
                              
    geom_cut = result['geom_cut']
    candidates = []
    
    #verts = []
    count = 0
    center_mass = Vector((0,0,0))
    # get the curve's center of mass

    for item in geom_cut:
        if isinstance(item, bmesh.types.BMVert):
            if target_torso_left>=item.co[0]>=target_torso_right:
                #verts.append(item)
                center_mass = center_mass + item.co
                count+=1
        elif isinstance(item,bmesh.types.BMEdge):
            vert1 = item.verts[0].co[0]
            vert2 = item.verts[1].co[0]
            if (vert1>vert2 and vert1>=x_value>=vert2) or (vert1<vert2 and vert2>=x_value>=vert1):
                """print(item.verts[0].co, item.verts[1].co, x_value,'check')
                intersection = mathutils.geometry.intersect_line_plane(item.verts[0].co,item.verts[1].co,[x_value,0,0],[-1,0,0])
                if intersection != None:
                    candidates.append(intersection)
                    """
                candidates.append([x_value,0.5*(item.verts[0].co[1]+item.verts[1].co[1]),z_value])
                    
    center_mass /= count   
      
    if len(candidates)!=2:
        print('length of candidates is not 2 but: ',len(candidates))       
    for pt in candidates:
        if v_src[1]>0:
            if pt[1]>center_mass[1]:
                return center_mass, Vector(pt)
        else:
            if pt[1]<center_mass[1]:
                return center_mass, Vector(pt)
    
    return center_mass, None


def get_target_cross_section_center_mass(value, horizontal=True,side='None',v_src=[]):
    for f in bm_target.faces:
        f.select = True
    edges = [e for e in bm_target.edges]
    faces = [f for f in bm_target.faces]
    geom = []
    geom.extend(edges)
    geom.extend(faces)
    
    if horizontal:
        co = [0,0,value]
        no = [0,0,1]
    else:
        co = [value,0,0]
        no = [1,0,0]
        
    # get cross section curve
    result = bmesh.ops.bisect_plane(bm_target,
                              dist=0.001,
                              geom=geom,
                              plane_co=co,
                              plane_no=no)
                              
    geom_cut = result['geom_cut']
    verts = []
    edges = []
    count = 0
    center_mass = Vector((0,0,0))
    # get the curve's center of mass
    if side=='None':
        if value < 0.70:
            if v_src[0]>0:
                for item in geom_cut:
                    if isinstance(item, bmesh.types.BMVert) and item.co[0]>0:
                        center_mass = center_mass + item.co
                        count+=1
            else:
                for item in geom_cut:
                    if isinstance(item, bmesh.types.BMVert) and item.co[0]<0:
                        center_mass = center_mass + item.co
                        count+=1
        else:
            for item in geom_cut:
                if isinstance(item, bmesh.types.BMVert):
                    if target_torso_left>=item.co[0]>=target_torso_right:
                        center_mass = center_mass + item.co
                        count+=1
        center_mass /= count                    
    else:
        for item in geom_cut:
            if isinstance(item, bmesh.types.BMVert):
                if (side=='left' and item.co[0]>0) or (side=='right' and item.co[0]<0):
                    center_mass = center_mass + item.co
                    count+=1
        center_mass /= count
    
    #make_temp_mesh([center_mass, v.co])
    #print(geom_cut)
    
    return center_mass


def compute_correspondence(v, source, target):
    
    # get gundam-grid segment that v is involved in
    seg_idx = ETC
    if v.co[2]<source_torso_bottom:
        if v.co[0]>0:
            seg_idx = LEFTLEG
        else:
            seg_idx = RIGHTLEG
    elif v.co[2]<source_torso_up:
        if v.co[0]<source_torso_right:
            seg_idx = RIGHTARM
        elif v.co[0]>source_torso_left:
            seg_idx = LEFTARM
        else:
            seg_idx = TORSO
    
    if seg_idx==ETC:
        print('invalid segment!')
        return (v.index,None)
    
    horizontal = not (LEFTARM<=seg_idx<=RIGHTARM)
    
    if seg_idx==LEFTLEG or seg_idx==LEFTARM:
        side='left'
    elif seg_idx==TORSO:
        side='None'
    else:
        side='right'
    source_vec = get_source_cross_section_vector(v,horizontal=horizontal,side=side)
    
    
    # compute main axis scale
    if not horizontal:
        scale = abs(v.co[0]-(bboxs_source[seg_idx][CENTER][0]-0.5*bboxs_source[seg_idx][X_LEN]))/bboxs_source[seg_idx][X_LEN]
        cross_section = bboxs_target[seg_idx][CENTER][0]-0.5*bboxs_target[seg_idx][X_LEN]+scale*bboxs_target[seg_idx][X_LEN]
    else:
        # get z scale
        scale = abs(v.co[2]-(bboxs_source[seg_idx][CENTER][2]-0.5*bboxs_source[seg_idx][Z_LEN]))/bboxs_source[seg_idx][Z_LEN]
        cross_section = bboxs_target[seg_idx][CENTER][2]-0.5*bboxs_target[seg_idx][Z_LEN]+scale*bboxs_target[seg_idx][Z_LEN]
    """
    if seg_idx==TORSO and cross_section>=target_chest_bottom:
        # get additional x scale....
        x_scale = abs(v.co[0]-(bboxs_source[seg_idx][CENTER][0]-0.5*bboxs_source[seg_idx][X_LEN]))/bboxs_source[seg_idx][X_LEN]
        x_value = bboxs_target[seg_idx][CENTER][0]-0.5*bboxs_target[seg_idx][X_LEN]+scale*bboxs_target[seg_idx][X_LEN]
        target_center_mass, loc= get_target_chest_correspondence(cross_section, x_value,v.co)
        if loc!=None:
            res, loc, nor, idx = target.ray_cast(target_center_mass, loc-target_center_mass)
        else:
            res = False
        
    else:
        target_center_mass = get_target_cross_section_center_mass(cross_section, horizontal=horizontal, side=side,v_src=v.co)
        res, loc, nor, idx = target.ray_cast(target_center_mass, source_vec)
    """
    target_center_mass = get_target_cross_section_center_mass(cross_section, horizontal=horizontal, side=side,v_src=v.co)
    res, loc, nor, idx = target.ray_cast(target_center_mass, source_vec)
    
    if res:
        if seg_idx==TORSO and not(target_torso_left>=loc[0]>=target_torso_right):
            # find closest point
            if loc[0]>target_torso_left:
                loc[0] = target_torso_left
            else:
                loc[0] = target_torso_right
        
        return (v.index,loc,nor)
    else:
        print('no ray cast result for idx: ',v.index)
        return (v.index,None)


###############################################
#               Deformation for               #
# Characters with Different Mesh Connectivity #
###############################################

source = bpy.data.objects['source_f']
garment = bpy.data.objects['garment_pants']
target = bpy.data.objects['target_m']

bm_source = bmesh.new()
bm_source.from_mesh(bpy.data.objects['source_f'].data)
bm_target = bmesh.new()
bm_target.from_mesh(bpy.data.objects['target_m'].data)

garment_vertices = garment.data.vertices
source_vertices = source.data.vertices
target_vertices = target.data.vertices

garment_polygons = garment.data.polygons
source_polygons = source.data.polygons

deformed_garment_vertices = []


def get_corresponding_vertex(source_idx):
    global unmatched_correspondences
    for v in correspondences:
        if source_idx==v[0]:
            if v[1]==None:
                return None
            else:
                return v
    # not computed yet
    res = compute_correspondence(source_vertices[face_vertex_index[i]],source,target)
    correspondences.append(res)
    
    
    if res[1]!=None:
        return res
    else:
        unmatched_correspondences = unmatched_correspondences+1
        return None

unmatched_correspondences = 0
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
        #co = compute_correspondence(source_vertices[face_vertex_index[i]],source,target)
        
        co = get_corresponding_vertex(face_vertex_index[i])
        
        if co==None:
            #print('correspondence doesnt exist for segment')
            check = True
            break
        else:
            target_projection_point += bary_weights[i] * co[1]
        
    
    if check:
        continue
    scale = 0.02#1.5*(v.co-loc).length/co[2].length#0.03
    # may want to customize per vertex here (intersection check, post-processing)
    deformed_garment_vertices.append(target_projection_point + scale*co[2])
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
        

print('Correspondence Evaluation:',unmatched_correspondences, '/',len(correspondences))
print('Garment Polygon Evaluation:',len(face_list),'/',len(garment_polygons))
mesh.from_pydata(deformed_garment_vertices, [], face_list)
