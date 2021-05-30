import bpy
from mathutils.interpolate import poly_3d_calc
from mathutils import Vector
import bmesh

default = bpy.data.objects['default']
garment = bpy.data.objects['garment']
input = bpy.data.objects['input']

garment_vertices = garment.data.vertices
default_vertices = default.data.vertices
input_vertices = input.data.vertices

garment_polygons = garment.data.polygons
default_polygons = default.data.polygons

deformed_garment_vertices = []


for v in garment_vertices:
    res, loc, nor, idx = default.closest_point_on_mesh(v.co)
    
    face_vertex_index = default_polygons[idx].vertices
    # interpolate projected point from default onto input
    default_face_vertices = [default_vertices[vx].co for vx in face_vertex_index]
    bary_weights = poly_3d_calc(default_face_vertices, loc)
    input_projection_point = Vector((0,0,0))

    for i in range(3):#vx in face_vertex_index:
        input_projection_point += bary_weights[i] * input_vertices[face_vertex_index[i]].co
        
    scale = 1.0
    deformed_garment_vertices.append(scale * (input_projection_point + (v.co - loc)))
    
    """
    fac_vertices = [ default_vertices[v].co for v in default_face_verticesidx].vertices]
    bary_weights = poly_3d_calc(fac_vertices, loc)
    projection_input = Vector((0,0,0))
    i=0
    for vx in body.data.polygons[idx].vertices:
        projection_input += bary_weights[i] * input_vertices[vx].co
        i = i+1
    
    scale = 1.0
    deformed_garment_vertices.append(scale*(projection_input+(v.co-loc)))
    
    """
    #print(projection_input)
    #input_fac_vertices = [input_vertices[v].co for v in body.data.polygons[idx].vertices]
    #print(input_fac_vertices)
    
    # get barycentric coordinates
    # corresponding point in other mesh 
    # warp to align normal between faces
    
    # how about scaling...
    


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