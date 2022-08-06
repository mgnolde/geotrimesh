import os
import sys
from pathlib import Path
import bpy,os

argv = sys.argv
dash_pos = argv.index("--")
out_dirpath = argv[dash_pos+1:][0]
component = argv[dash_pos+1:][1]


for filename in os.listdir(out_dirpath):
    filename_stem = str(Path(filename).stem)
    filename_suffix = str(Path(filename).suffix)

    in_filepath = str(Path(out_dirpath,filename))
    out_filepath = str(Path(out_dirpath,filename_stem + "_zup" + ".obj"))

    if (filename_suffix == ".glb") and (component in filename_stem) and (not "zup" in filename_stem):

        print(filename)

        if os.path.isfile(out_filepath):
            os.remove(out_filepath)

        bpy.ops.import_scene.gltf(filepath=in_filepath)
        #bpy.ops.export_scene.x3d(filepath=out_filepath,axis_forward="-Y",axis_up="X")
        bpy.ops.export_scene.obj(filepath=out_filepath,axis_forward="-Y",axis_up="X",use_materials=True)

        #newmtl Material
        #Ns 323.999994
        #Ka 1.000000 1.000000 1.000000
        #Kd 0.800000 0.800000 0.800000
        #Ks 0.500000 0.500000 0.500000
        #Ke 0.000000 0.000000 0.000000
        #Tr 0.000000
        #Ni 10.000000
        #d 1.000000
        #illum 2
        #map_Kd ortho__1_0.png

