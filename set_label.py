import bpy
import numpy as np
obj = bpy.data.objects['F0']
len_obj = len(obj.data.vertices)

"""
bpy.ops.object.mode_set(mode = 'OBJECT')
obj = bpy.context.active_object
bpy.ops.object.mode_set(mode = 'EDIT') 
bpy.ops.mesh.select_mode(type="VERT")
bpy.ops.mesh.select_all(action = 'DESELECT')
bpy.ops.object.mode_set(mode = 'OBJECT')
obj.data.vertices[0].select = True
bpy.ops.object.mode_set(mode = 'EDIT') 
"""

vg_name_list = ["Head", "Chest", "Abdomen", "Pelvis","LeftArm", "LeftForearm","LeftHand", "RightArm", "RightForearm", "RightHand","LeftUpperLeg", "LeftLowerLeg", "LeftFoot", "RightUpperLeg", "RightLowerLeg", "RightFoot"]

check = False
for v in obj.data.vertices:
    if len(v.groups) != 1:
        check = True
        print(v)

if check==False:
    print('everything looks fine')

vgroups = obj.vertex_groups
labels = []
for v in obj.data.vertices:
    group_idx = v.groups[0].group
    
    group_name = vgroups[group_idx].name
    # find index of name in vg_name_list
    idx = vg_name_list.index(group_name)
    if idx<0 or idx> 15:
        print('couldnt find index for '+v)
        break
    # push in labels list
    labels.append(idx)

if len(labels)==len(obj.data.vertices):
    print('label created for each vertex')

np.save('saved_labels',labels)