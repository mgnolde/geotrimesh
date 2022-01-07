import bpy
import os
import sys
import math

input_filename = sys.argv[5]
output_filename = sys.argv[6]

try:
	xrot = float(sys.argv[7])
	yrot = float(sys.argv[8])
	zrot = float(sys.argv[9])
except:
	xrot = 0
	yrot = 0
	zrot = 0


print(input_filename, output_filename)

bpy.ops.object.select_all(action="SELECT")
bpy.ops.object.delete() 


if input_filename.split('.')[-1].lower() == "obj":
    bpy.ops.import_scene.obj(filepath=input_filename)
elif input_filename.split('.')[-1].lower() == "fbx":
    bpy.ops.import_scene.fbx(filepath=input_filename)
elif input_filename.split('.')[-1].lower() == "gltf":
    bpy.ops.import_scene.gltf(filepath=input_filename)
elif input_filename.split('.')[-1].lower() == "x3d":
    bpy.ops.import_scene.x3d(filepath=input_filename)
elif input_filename.split('.')[-1].lower() == "dxf":
    bpy.ops.import_scene.dxf(filepath=input_filename)
else:
    print("Input file format not supported")
    sys.exit()


objects = bpy.context.scene.objects

for obj in objects:
    obj.select_set(obj.type == "EMPTY")
    obj.select_set(True)    
    bpy.context.view_layer.objects.active = obj


#bpy.ops.object.select_all(action="SELECT")

#obj = bpy.context.active_object

#ob = bpy.context.object
obj.rotation_euler = (xrot,yrot,zrot) 
#obj.rotation_euler = (0,0,math.radians(-90)) 
#obj.rotation_euler = (0,math.radians(-90),0) 
#obj.rotation_euler = (math.radians(-90),0,0) 

#obj.rotation_euler[0] = math.radians(-90)
#obj.rotation_euler[1] = math.radians(-90)
#obj.rotation_euler[2] = math.radians(-90)
 

if output_filename.split('.')[-1].lower() == "gltf":
    bpy.ops.export_scene.gltf(filepath=output_filename, export_yup=False)
elif output_filename.split('.')[-1].lower() == "obj":
    bpy.ops.export_scene.obj(filepath=output_filename, export_yup=False)
    #bpy.ops.export_scene.obj(filepath=output_filename[:-5] + '_test.obj')
else:
    print("Output file format not supported")
    sys.exit()

bpy.ops.wm.quit_blender()
sys.exit()
