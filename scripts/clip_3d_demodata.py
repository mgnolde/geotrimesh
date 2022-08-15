import os
import sys
from pathlib import Path
import trimesh
import geopandas as gpd
from shapely.geometry import Point, box, MultiPolygon, MultiPoint
from shapely.validation import make_valid



def _check_face_inside_box(vertices, face, box, offset):

    x_offset, y_offset = offset

    points = [
        [vertices[face[0]][0], vertices[face[0]][1]],
        [vertices[face[1]][0], vertices[face[1]][1]],
        [vertices[face[2]][0], vertices[face[2]][1]],
    ]


    inside = False
    for point in points:
        p = Point(point[0] + x_offset, point[1] + y_offset)
        if p.intersects(box):
            inside = True

    return inside

## Get these values by opening the original file in the FZKViewer (appear on load)
x_min_orig = 2677116.375000
y_min_orig = 1241839.025000
x_max_orig = 2689381.985000
y_max_orig = 1254150.950000



data_dirpath = Path(os.getcwd(), "demodata")
boundary_filepath = Path(data_dirpath, "bbox.gpkg")
buildings_filepath = Path(os.sep,"home","mic","devenv","geotrimesh","data","orig","zh","feat","buildings","ganze_stadt.glb")
buildings_out_filepath = Path(os.sep,"home","mic","devenv","geotrimesh","data","orig","zh","feat","buildings","ganze_stadt_demo.glb")

boundary = gpd.read_file(boundary_filepath).dissolve().explode(index_parts=True)
print(boundary.total_bounds)

margin = 100
boundary_box = box(
    boundary.total_bounds[0]-margin, 
    boundary.total_bounds[1]-margin, 
    boundary.total_bounds[2]+margin, 
    boundary.total_bounds[3]+margin
    )

scene = trimesh.load(buildings_filepath)
scene_out = trimesh.Scene()

xyz_min, xyz_max = scene.bounds
x_min, y_min, z_min = xyz_min
x_max, y_max, z_max = xyz_max

x_offset = x_min_orig - x_min 
y_offset = y_min_orig - y_min 


for key_id, key in enumerate(scene.geometry):
    if type(scene.geometry[key]).__name__ == "Trimesh":
        mesh = scene.geometry[key]
        

        print(key_id, "faces:", mesh.faces.shape[0], "vertices:", mesh.vertices.shape[0], "bounds:", mesh.bounds)

        face_mask = []
        for face in mesh.faces:
            is_inside = _check_face_inside_box(mesh.vertices, face, boundary_box, (x_offset, y_offset))
            face_mask.append(is_inside)
            #print(is_inside)

        mesh.update_faces(face_mask)
        print(key_id, "faces:", mesh.faces.shape[0], "vertices:", mesh.vertices.shape[0], "bounds:", mesh.bounds)



        vertices_transformed = []
        for vertex_id, vertex in enumerate(mesh.vertices):
            vertices_transformed.append(
                [
                    vertex[0] + x_offset,
                    vertex[1] + y_offset,
                    vertex[2],
                ]
            )


        mesh_new = trimesh.Trimesh(
            vertices=vertices_transformed, 
            faces=mesh.faces,
            face_normals=mesh.face_normals, 
            vertex_normals=mesh.vertex_normals
        )

        scene_out.add_geometry(
            mesh_new, node_name=str(key_id), geom_name=str(key_id)
        )



with open(buildings_out_filepath,"wb") as f:
    f.write(
        trimesh.exchange.gltf.export_glb(
            scene_out, include_normals=True
        )
    )