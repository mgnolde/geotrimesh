import os
import sys
from shapely.validation import make_valid
import vtk
from shapely.geometry import (
    Point,
    Polygon,
    mapping,
    shape,
    box,
    MultiPolygon,
    MultiPoint,
)
import trimesh
import copy
import numpy as np
import geopandas as gpd
from shapely.ops import unary_union
from shapely.geometry import MultiPolygon
import logging
import tempfile
from pathlib import Path
import rasterio
from rasterio.windows import from_bounds
import rasterio.crs
import math
from collections import namedtuple
import subprocess
from PIL import Image, ImageOps
from rasterio.transform import Affine
import open3d



def _get_mesh_from_scene(scene):

    meshes = []

    for key in scene.geometry:
        if type(scene.geometry[key]).__name__ == "Trimesh":
            mesh_orig = scene.geometry[key]
            mesh_orig_new = trimesh.Trimesh(vertices=mesh_orig.vertices, faces=mesh_orig.faces)
            meshes.append(mesh_orig_new)

    return trimesh.util.concatenate(meshes)


class GeoSceneSet():

    def __init__(self):
        self.terrain = None
        self.features = None
        self.ortho = None
        self.dem = None

    class TilingScheme():

        def __init__(self, boundary, filepaths, height=None, width=None):

            self.id = None
            self.scheme_dims = None
            self.dimns = None
            self.res = None
            self.gdf = None

            self.tiles = []

            gdf_dict = {}
            gdf_dict["geometry"] = []


            upper_left = Point(boundary.bounds.minx, boundary.bounds.maxy)
            lower_right = Point(boundary.bounds.maxx, boundary.bounds.miny)

            center = ((float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0, (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0)
            scale = (100000.0, 100000.0, 100000.0)

            upper_left_aoi, res, nodata = self.define_aoi_origin(upper_left, filepaths)

            tile = namedtuple("tile", ["id", "scheme_dims", "dims", "res", "nodata", "geom"])


            if not (width and height):

                height = math.ceil((upper_left.y - lower_right.y) / res[1]) + 1
                width = math.ceil((lower_right.x - upper_left.x) / res[0]) + 1

                geom = box(
                    upper_left_aoi.x, 
                    upper_left_aoi.y - (height * res[1]), 
                    upper_left_aoi.x + (width * res[0]), 
                    upper_left_aoi.y
                )

                geom_target = box(
                    (upper_left_aoi.x - center[0]) * scale[0], 
                    ((upper_left_aoi.y - (height * res[1])) - center[1]) * scale[1], 
                    ((upper_left_aoi.x + (width * res[0])) - center[0]) * scale[0], 
                    (upper_left_aoi.y - center[1]) * scale[1]
                )




                gdf_dict["geometry"].append(geom)

                tile = namedtuple("tile", ["id", "scheme_dims", "dims", "res", "nodata", "geom"])

                tile.id = (0,0)
                tile.scheme_dims = (1, 1)
                tile.dims = (height, width)
                tile.res = res
                tile.nodata = nodata
                tile.geom = geom
                tile.geom_target = geom_target
                tile.total_bounds = None


                self.tiles.append(tile)

            else:

                count_x = math.ceil(((lower_right.x - upper_left_aoi.x) / res[0]) / width)
                count_y = math.ceil(((upper_left_aoi.y - lower_right.y) / abs(res[1])) / height)

                for y in range(0,count_y):
                    for x in range(0,count_x):

                        geom = box(
                            upper_left_aoi.x + (x * width * res[0]), 
                            upper_left_aoi.y - ((y+1) * height * res[1]), 
                            upper_left_aoi.x + ((x+1) * width * res[0]), 
                            upper_left_aoi.y - (y * height * res[1])
                        )


                        geom_target = box(
                            (upper_left_aoi.x - center[0]) * scale[0], 
                            ((upper_left_aoi.y - (height * res[1])) - center[1]) * scale[1], 
                            ((upper_left_aoi.x + (width * res[0])) - center[0]) * scale[0], 
                            (upper_left_aoi.y - center[1]) * scale[1]
                        )


                        gdf_dict["geometry"].append(geom)

                        tile = namedtuple("tile", ["id", "scheme_dims", "dims", "res", "nodata", "geom"])

                        tile.id = (y,x)
                        tile.scheme_dims = (count_y, count_x)
                        tile.dims = (height, width)
                        tile.res = res
                        tile.nodata = nodata
                        tile.geom = geom
                        tile.geom_target = geom_target
                        tile.total_bounds = None



                        self.tiles.append(tile)


            self.gdf = gpd.GeoDataFrame.from_dict(gdf_dict)
            self.gdf = self.gdf.set_crs(boundary.crs, allow_override=True)


        def define_aoi_origin(self, point, filepaths):

            for filepath in filepaths:
                src = rasterio.open(filepath)
                src_box = box(*src.bounds)

                if point.intersects(src_box):
                    dist_left_pix = math.floor((point.x - src.bounds.left) / src.res[0])
                    dist_top_pix = math.floor((src.bounds.top - point.y) / src.res[1])

                    origin = Point(src.bounds.left + (dist_left_pix * src.res[0]), src.bounds.top - (dist_top_pix * src.res[1]))

                    if src.nodatavals[0]:
                        nodata = src.nodatavals[0]
                    else:
                        nodata = -99999

                    return(origin, src.res, nodata)


    class Ortho():

        def __init__(
            self,
            tilingscheme=None,
            out_dirpath=tempfile.gettempdir(),
            filepaths=None,
            tiles=None,
            boundary=None,
        ):

            for tile in tiles:
                for glb in [["terrain",str(Path(out_dirpath, "terrain" + "__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".glb"))],
                    ["buildings",str(Path(out_dirpath, "buildings" + "__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".glb"))]]:

                    mosaic, width, height = self.generate_ortho_mosaic(filepaths, tile.total_bounds)

                    mosaic_rearranged = np.transpose(mosaic, axes=[1, 2, 0])

                    newImg1 = Image.fromarray(mosaic_rearranged.astype('uint8'), 'RGB')
                    newImg1.save(Path(out_dirpath, "ortho__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".png"),"PNG")

                    im = Image.open(Path(out_dirpath, "ortho__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".png"))
                    im = ImageOps.flip(im)
                    scene = trimesh.load(str(glb[1]))
                    scene_out = trimesh.Scene()

                    meshes = []
                    center = ((float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0, (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0)
                    scale = (100000.0, 100000.0, 100000.0)


                    for key in scene.geometry:
                        if type(scene.geometry[key]).__name__ == "Trimesh":
                            mesh_orig = scene.geometry[key]


                            uv = []
                            for vert in mesh_orig.vertices:
                                print(vert)
                                x_local, y_local, z_local = vert

                                x = (float(x_local) / scale[0]) + center[0]
                                y = (float(y_local) / scale[1]) +  center[1]

                                x_offset = x - tile.total_bounds[0]
                                y_offset = tile.total_bounds[3] - y

                                dem_x_dist = tile.total_bounds[2] - tile.total_bounds[0]
                                dem_y_dist = tile.total_bounds[3] - tile.total_bounds[1]

                                x_ratio = round(x_offset / float(dem_x_dist), 2)
                                y_ratio = round(y_offset / float(dem_y_dist), 2)
                                uv.append(np.array([x_ratio,y_ratio]))

                            material = trimesh.visual.texture.SimpleMaterial(image=im)
                            color_visuals = trimesh.visual.TextureVisuals(uv=uv, image=im, material=material)
                            mesh_orig_new = trimesh.Trimesh(vertices=mesh_orig.vertices, faces=mesh_orig.faces, visual=color_visuals, validate=True, process=False)

                            scene_out.add_geometry(mesh_orig_new, node_name="mesh", geom_name="mesh")


                    with open(Path(out_dirpath, "result_" + glb[0] + "__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".glb"), 'wb') as f:
                        f.write(trimesh.exchange.gltf.export_glb(scene_out, include_normals=True))



        def generate_ortho_mosaic(self, filepaths, bounds):

            left, bottom, right, top = bounds
            with rasterio.open(str(filepaths[0])) as src:
                res_x, res_y = src.res

            with rasterio.open(filepaths[0]) as src:
                target_array = src.read(1, window=from_bounds(left, bottom, right, top, src.transform), boundless=True, fill_value=0)
                rast = np.zeros((3, target_array.shape[0], target_array.shape[1]), dtype=np.uint8)

            for filepath in filepaths:
                loc = np.amax(rast,axis=0)==0
                with rasterio.open(filepath) as src:
                    for i in range(0,3):
                        rast[i][loc] = src.read(i+1, window=from_bounds(left, bottom, right, top, src.transform), boundless=True, fill_value=0)[loc]

            return rast, target_array.shape[0], target_array.shape[1]


    class Dem():

        def __init__(self, filepaths, bounds, tile=None):
            return None


    class Terrain():
    
        def __init__(self, out_dirpath=tempfile.gettempdir(), tiles=None, filepaths=[], tile=None, boundary=None, openscad_bin_filepath=Path("openscad")):

            for tile in tiles:

                self.mosaic = self.generate_terrain_mosaic(filepaths, tile)
                self.mesh = self.generate_terrain_mesh(self.mosaic, boundary, tile, out_dirpath, openscad_bin_filepath)

                open3d_mesh = open3d.geometry.TriangleMesh()
                open3d_mesh.vertices = open3d.utility.Vector3dVector(self.mesh.vertices)
                open3d_mesh.triangles = open3d.utility.Vector3iVector(self.mesh.faces)

                open3d_mesh = open3d_mesh.simplify_quadric_decimation(target_number_of_triangles=int(round(len(self.mesh.vertices)/100,2)))

                self.mesh = trimesh.Trimesh(vertices=open3d_mesh.vertices, faces=open3d_mesh.triangles)


                self.scene = trimesh.Scene()
                self.scene.add_geometry(self.mesh, node_name="mesh", geom_name="mesh")
                tile.total_bounds = tile.geom.bounds


                with open(Path(out_dirpath, "terrain__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".glb"), 'wb') as f:
                    f.write(trimesh.exchange.gltf.export_glb(self.scene, include_normals=True))



        def generate_terrain_mosaic(self, filepaths, tile):

            rast = np.ones(tile.dims, dtype=np.float32) * tile.nodata

            for filepath in filepaths:
                with rasterio.open(filepath) as src:
                    left, bottom, right, top = tile.geom.bounds
                    rast_tmp = src.read(1, window=from_bounds(left, bottom, right, top, src.transform), boundless=True, fill_value=src.nodatavals[0])

                    rast[rast<5000] = rast_tmp[rast<5000]

            return rast


        def _get_polyhedron_faces(self, rast, tile, boundary):

            left, bottom, right, top = tile.geom.bounds
            center = ((float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0, (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0)
            scale = (100000.0, 100000.0, 100000.0)

            polyhedron_faces_array = np.zeros(
                (tile.dims[0], tile.dims[1], 16, 3), dtype=np.int32
            )
            polyhedron_points_floor_array = np.zeros(
                (tile.dims[0] * tile.dims[1], 3), dtype=np.float32
            )

            cnt = 0
            dem_x_min = None
            dem_x_max = None
            dem_y_min = None
            dem_y_max = None

            for i in range(0, tile.dims[0]):

                logging.info(f'Row: {i} of {tile.dims[0]}')

                for j in range(0, tile.dims[1]):

                    i0_coord = top - (
                        tile.res[1] * i
                    ) 
                    j0_coord = left + (tile.res[0] * j) 

                    z_a = rast[i][j] * scale[2]

                    dem_x = j0_coord
                    dem_y = i0_coord

                    polyhedron_points_floor_array[(i * tile.dims[1]) + j][0] = (dem_x - center[0]) * scale[0]
                    polyhedron_points_floor_array[(i * tile.dims[1]) + j][1] = (dem_y - center[1]) * scale[1]
                    polyhedron_points_floor_array[(i * tile.dims[1]) + j][2] = z_a


                    if not dem_x_min or dem_x < dem_x_min:
                        dem_x_min = dem_x
                    if not dem_x_max or dem_x > dem_x_max:
                        dem_x_max = dem_x

                    if not dem_y_min or dem_y < dem_y_min:
                        dem_y_min = dem_y
                    if not dem_y_max or dem_y > dem_y_max:
                        dem_y_max = dem_y

                    if i < tile.dims[0] - 1 and j < tile.dims[1] - 1:

                        z_b = rast[i + 1][j] * scale[2]
                        z_c = rast[i][j + 1] * scale[2]
                        z_d = rast[i + 1][j + 1] * scale[2]

                        point_a_ceil = (i * tile.dims[1]) + j
                        point_b_ceil = ((i + 1) * tile.dims[1]) + j
                        point_c_ceil = (i * tile.dims[1]) + j + 1
                        point_d_ceil = ((i + 1) * tile.dims[1]) + j + 1

                        point_a_floor = (tile.dims[0] * tile.dims[1]) + (i * tile.dims[1]) + j
                        point_b_floor = (tile.dims[0] * tile.dims[1]) + ((i + 1) * tile.dims[1]) + j
                        point_c_floor = (tile.dims[0] * tile.dims[1]) + (i * tile.dims[1]) + j + 1
                        point_d_floor = (
                            (tile.dims[0] * tile.dims[1]) + ((i + 1) * tile.dims[1]) + j + 1
                        )

                        if z_a > -5000 and z_b > -5000 and z_c > -5000:

                            polyhedron_faces_array[i][j][0][:] = np.array(
                                [point_c_ceil, point_b_ceil, point_a_ceil]
                            )  ## ceiling
                            polyhedron_faces_array[i][j][1][:] = np.array(
                                [point_b_floor, point_a_floor, point_a_ceil]
                            )  ## left sidev
                            polyhedron_faces_array[i][j][2][:] = np.array(
                                [point_b_ceil, point_b_floor, point_a_ceil]
                            )
                            polyhedron_faces_array[i][j][3][:] = np.array(
                                [point_c_floor, point_b_floor, point_b_ceil]
                            )  ## right side (diagonal)
                            polyhedron_faces_array[i][j][4][:] = np.array(
                                [point_c_ceil, point_c_floor, point_b_ceil]
                            )
                            polyhedron_faces_array[i][j][5][:] = np.array(
                                [point_a_floor, point_c_floor, point_c_ceil]
                            )  ## top side
                            polyhedron_faces_array[i][j][6][:] = np.array(
                                [point_a_ceil, point_a_floor, point_c_ceil]
                            )
                            polyhedron_faces_array[i][j][7][:] = np.array(
                                [point_a_floor, point_b_floor, point_c_floor]
                            )  ## floor

                        if z_b > -5000 and z_d > -5000 and z_c > -5000:

                            polyhedron_faces_array[i][j][8][:] = np.array(
                                [point_c_ceil, point_d_ceil, point_b_ceil]
                            )  ## ceiling
                            polyhedron_faces_array[i][j][9][:] = np.array(
                                [point_d_ceil, point_d_floor, point_b_floor]
                            )  ## bottom side
                            polyhedron_faces_array[i][j][10][:] = np.array(
                                [point_d_ceil, point_b_floor, point_b_ceil]
                            )
                            polyhedron_faces_array[i][j][11][:] = np.array(
                                [point_c_ceil, point_c_floor, point_d_floor]
                            )  ## right side
                            polyhedron_faces_array[i][j][12][:] = np.array(
                                [point_c_ceil, point_d_floor, point_d_ceil]
                            )
                            polyhedron_faces_array[i][j][13][:] = np.array(
                                [point_b_ceil, point_b_floor, point_c_floor]
                            )  ## left side (diagonal)
                            polyhedron_faces_array[i][j][14][:] = np.array(
                                [point_c_ceil, point_b_ceil, point_c_floor]
                            )
                            polyhedron_faces_array[i][j][15][:] = np.array(
                                [point_b_floor, point_d_floor, point_c_floor]
                            )  ## floor

                    cnt += 1

            return polyhedron_faces_array, polyhedron_points_floor_array


        def _get_polyhedron_faces_clean(self, polyhedron_faces_array, polyhedron_points_floor_array, tile):

            polyhedron_faces_clean = []

            for i in range(0, tile.dims[0]):

                for j in range(0, tile.dims[1]):

                    for polyhedron_face_id in range(0, 16):

                        polyhedron_face = polyhedron_faces_array[i][j][polyhedron_face_id][
                            :
                        ].tolist()
                        polyhedron_face_cnt = 0

                        if not (
                            polyhedron_face[0] == 0
                            and polyhedron_face[1] == 0
                            and polyhedron_face[2] == 0
                        ):

                            i_bottom = i - 1 if i > 0 else 0
                            i_top = i + 1 if i < tile.dims[0] - 1 else tile.dims[0] - 1
                            j_left = j - 1 if i > 0 else 0
                            j_right = j + 1 if j < tile.dims[1] - 1 else tile.dims[1]

                            for i_neighbour in range(i_bottom, i_top + 1):
                                for j_neighbour in range(j_left, j_right + 1):

                                    array1 = np.array(
                                        [
                                            polyhedron_face[0],
                                            polyhedron_face[1],
                                            polyhedron_face[2],
                                        ]
                                    )
                                    array2 = np.array(
                                        [
                                            polyhedron_face[0],
                                            polyhedron_face[2],
                                            polyhedron_face[1],
                                        ]
                                    )
                                    array3 = np.array(
                                        [
                                            polyhedron_face[1],
                                            polyhedron_face[0],
                                            polyhedron_face[2],
                                        ]
                                    )
                                    array4 = np.array(
                                        [
                                            polyhedron_face[1],
                                            polyhedron_face[2],
                                            polyhedron_face[0],
                                        ]
                                    )
                                    array5 = np.array(
                                        [
                                            polyhedron_face[2],
                                            polyhedron_face[1],
                                            polyhedron_face[0],
                                        ]
                                    )
                                    array6 = np.array(
                                        [
                                            polyhedron_face[2],
                                            polyhedron_face[0],
                                            polyhedron_face[1],
                                        ]
                                    )

                                    loc1 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][j_neighbour][:]
                                            == array1,
                                            axis=1,
                                        )
                                    )
                                    loc2 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][j_neighbour][:]
                                            == array2,
                                            axis=1,
                                        )
                                    )
                                    loc3 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][j_neighbour][:]
                                            == array3,
                                            axis=1,
                                        )
                                    )
                                    loc4 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][j_neighbour][:]
                                            == array4,
                                            axis=1,
                                        )
                                    )
                                    loc5 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][j_neighbour][:]
                                            == array5,
                                            axis=1,
                                        )
                                    )
                                    loc6 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][j_neighbour][:]
                                            == array6,
                                            axis=1,
                                        )
                                    )

                                    polyhedron_face_neighbour_cnt = (
                                        len(loc1[0])
                                        + len(loc2[0])
                                        + len(loc3[0])
                                        + len(loc4[0])
                                        + len(loc5[0])
                                        + len(loc6[0])
                                    )
                                    polyhedron_face_cnt += polyhedron_face_neighbour_cnt

                            if polyhedron_face_cnt == 1:

                                polyhedron_faces_clean.append(polyhedron_face)

                            else:
                                pass

            return polyhedron_faces_clean



        def _get_polyhedron_points(self, polyhedron_points_floor_array, dem_extrude_height_up, dem_extrude_height_down):
            polyhedron_points = []

            for l in range(1, -1, -1):
                for polyhedron_point_id in range(0, polyhedron_points_floor_array.shape[0]):

                    if l == 1:
                        z_offset = dem_extrude_height_up
                    if l == 0:
                        z_offset = abs(dem_extrude_height_down) * -1


                    polyhedron_point = [
                        polyhedron_points_floor_array[polyhedron_point_id][0],
                        polyhedron_points_floor_array[polyhedron_point_id][1],
                        polyhedron_points_floor_array[polyhedron_point_id][2] - z_offset,
                    ]

                    polyhedron_points.append(polyhedron_point)

            return polyhedron_points


        def _get_polygon_intersect(self, boundary=None, tile=None, polygon_extrude_height_down=None, polygon_extrude_height_up=None, center=None, scale=None):

            center = ((float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0, (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0)

            left, bottom, right, top = tile.geom.bounds
            scale = (100000.0, 100000.0, 100000.0)

            boundary_tile = boundary.cx[left:right, bottom:top]
            polys = []


            if boundary_tile.shape[0] > 0:

                for index,row in boundary_tile.iterrows():
                    geom = row["geometry"]

                    polygon_points = []
                    polygon_paths = []
                    coord_id = 0
                    point_cnt = 0

                    polygon_path = []
                    for coord_x, coord_y in zip(geom.exterior.coords.xy[0],geom.exterior.coords.xy[1]):
                        polygon_path.append(coord_id)
                        polygon_points.append([(coord_x - center[0]) * scale[0], (coord_y - center[1]) * scale[1]])
                        coord_id +=1 
                    polygon_paths.append(polygon_path)

                    if len(geom.interiors) > 0:
                        for interior in geom.interiors:
                            polygon_path = []
                            for coord_x, coord_y in zip(interior.coords.xy[0],interior.coords.xy[1]):
                                polygon_path.append(coord_id)
                                polygon_points.append([(coord_x - center[0]) * scale[0], (coord_y - center[1]) * scale[1]])
                                coord_id +=1 
                            polygon_paths.append(polygon_path)

                    polys.append([polygon_paths, polygon_points])

            return polys


        def _write_scad(self, polyhedron_points, polyhedron_faces_clean, polys, polygon_extrude_height_down, polygon_extrude_height_up, scad_filepath):

            c = []
            c.append("module dem() {")
            c.append(
                "  polyhedron(points={}, faces={});".format(
                    polyhedron_points, polyhedron_faces_clean
                )
            )
            c.append("}")


            for poly in polys:
                polygon_paths, polygon_points = poly

                c.append("intersection() {")
                c.append(
                    "translate([0,0,{}]) linear_extrude(height={}) polygon({},{});".format(
                        abs(polygon_extrude_height_down) * -1,
                        abs(polygon_extrude_height_down) + abs(polygon_extrude_height_up),
                        polygon_points,
                        polygon_paths,
                    )
                )

                c.append("dem();")
                c.append("}")

            with open(scad_filepath, "w") as scad_file:
                for command_line in c:
                    scad_file.write(command_line + "\n")



        def _generate_stl(self, openscad_bin_filepath, scad_filepath, stl_filepath):
            subprocess_commands = [str(openscad_bin_filepath), scad_filepath, '-o', stl_filepath]
            output = subprocess.check_output(subprocess_commands, shell=False)

            stl = trimesh.load(str(stl_filepath))
            return stl


        def generate_terrain_mesh(self, rast=None, boundary=None, tile=None, out_dirpath=tempfile.gettempdir(), openscad_bin_filepath=Path("openscad")):

            scale = (100000.0, 100000.0, 100000.0)

            import logging
            polygon_extrude_height_down = 1000 * scale[2]
            polygon_extrude_height_up = 1000 * scale[2]
            dem_extrude_height_down = .5 * scale[2]
            dem_extrude_height_up = .5 * scale[2]



            polyhedron_faces_array, polyhedron_points_floor_array = self._get_polyhedron_faces(rast, tile, boundary)
            polyhedron_faces_clean = self._get_polyhedron_faces_clean(polyhedron_faces_array, polyhedron_points_floor_array, tile)

            polyhedron_points = self._get_polyhedron_points(polyhedron_points_floor_array, dem_extrude_height_up, dem_extrude_height_down)

            polys = self._get_polygon_intersect(boundary, tile, polygon_extrude_height_down, polygon_extrude_height_up)

            scad_filepath = str(Path(out_dirpath, "terrain" + "__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".scad"))
            stl_filepath = str(Path(out_dirpath, "terrain" + "__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".stl"))



            self._write_scad(polyhedron_points, polyhedron_faces_clean, polys, polygon_extrude_height_down, polygon_extrude_height_up, scad_filepath)

            mesh = self._generate_stl(openscad_bin_filepath, scad_filepath, stl_filepath)

            return mesh


    class Features():
    
        def __init__(
            self,
            description="feature",
            tilingscheme=None,
            out_dirpath=tempfile.gettempdir(),
            filepaths=None,
            recombine_bodies=False,
            simplify_factor=None,
            tiles=None,
            boundary=None,
        ):

            #tiles_total_bounds = {}
            for tile in tiles:

                self.process(
                    tilingscheme=tilingscheme,
                    out_dirpath=out_dirpath,
                    filepaths=filepaths,
                    recombine_bodies=recombine_bodies,
                    simplify_factor=simplify_factor,
                    tile=tile,
                    description=description,
                    boundary=boundary,
                )

                tile.geom.bounds

                tile.total_bounds = (
                    min(self.map2d.total_bounds[0], tile.total_bounds[0] if tile.total_bounds else tile.geom.bounds[0]),
                    min(self.map2d.total_bounds[1], tile.total_bounds[1] if tile.total_bounds else tile.geom.bounds[1]),
                    max(self.map2d.total_bounds[2], tile.total_bounds[2] if tile.total_bounds else tile.geom.bounds[2]),
                    max(self.map2d.total_bounds[3], tile.total_bounds[3] if tile.total_bounds else tile.geom.bounds[3])
                )

                with open(Path(out_dirpath, description + "__" + '_'.join(str(x) for x in tile.id) + ".glb"), 'wb') as f:
                    f.write(trimesh.exchange.gltf.export_glb(self.scene, include_normals=True))


        def process(self, tilingscheme=None,
            out_dirpath=tempfile.gettempdir(),
            filepaths=None,
            recombine_bodies=False,
            simplify_factor=None,
            tile=None,
            description="feat",
            boundary=None):

            if recombine_bodies:

                ## combine bodies which overlap in 2d space (xy) into common nodes

                bodies_multipolygons, bodies = self._get_map2d(filepaths=filepaths, tile=tile, boundary=boundary)

                map2d_unfiltered = gpd.GeoDataFrame.from_dict({"geometry": unary_union(bodies_multipolygons)})
                map2d_unfiltered = map2d_unfiltered.set_crs(boundary.crs, allow_override=True)

                features2d, features3d = self._do_map2d_binning(tilingscheme=tilingscheme, out_dirpath=out_dirpath, tile=tile, boundary=boundary, map2d_unfiltered=map2d_unfiltered, bodies_multipolygons=bodies_multipolygons, bodies=bodies)

                if len(features2d) > 0:
                    features2d_dict = {}
                    features2d_dict["geometry"] = []
                    for feature2d in features2d:
                        features2d_dict["geometry"].append(feature2d)

                    self.map2d = gpd.GeoDataFrame.from_dict(features2d_dict)
                    self.map2d = self.map2d.set_crs(boundary.crs, allow_override=True)


                    scene = trimesh.Scene()
                    for node_id, node in enumerate(features3d):
                        scene.add_geometry(node, node_name=str(node_id), geom_name=str(node_id))

                    self.scene = scene
                    self.mesh = _get_mesh_from_scene(self.scene)


                else:
                    self.scene = trimesh.Scene()
                    self.mesh = trimesh.Trimesh()
                    self.map2d = gpd.GeoDataFrame()

            else:

                scene = trimesh.Scene()
                for filepath in filepaths:

                    feature = trimesh.load(filepath)
                    mesh = _get_mesh_from_scene(feature)
                    for node_id, node in enumerate(mesh.geometry):
                        if node.centroid.intersection(tile.geom):
                            scene.add_geometry(node, node_name=str(node_id), geom_name=str(node_id))

                self.scene = scene
                self.mesh = _get_mesh_from_scene(self.scene)


        def _get_map2d(self, filepaths, tile=None, boundary=None):

            logging.info("Affiliation")

            bodies_total = []
            multipolygons_total = []

            for filepath in filepaths:

                feature = trimesh.load(filepath)
                mesh = _get_mesh_from_scene(feature)

                mesh_components = mesh.split(only_watertight=False)

                for body in mesh_components:

                    polygons_body = []

                    vertices = body.vertices
                    faces = body.faces

                    for face in faces:

                        face_poly = self._poly_from_triangle(vertices, face)

                        if face_poly and tile.geom.buffer((tile.geom.bounds[2]-tile.geom.bounds[0]), cap_style=3).intersects(face_poly):
                            polygons_body.append(face_poly)


                    if len(polygons_body) > 0:
                        multipolygon = MultiPolygon(polygons_body)                      
                        bodies_total.append(body)
                        multipolygons_total.append(multipolygon)
                        logging.info(f"Body added {len(bodies_total)}")


            return multipolygons_total, bodies_total


        def _do_map2d_binning(self, tilingscheme, out_dirpath=None, tile=None, boundary=None, map2d_unfiltered=None, bodies_multipolygons=None, bodies=None, scale=(1.0,1.0,1.0)):

            logging.info("Binning")
            center = ((float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0, (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0, 0)

            scale = (100000.0, 100000.0, 100000.0)

            features2d = []
            features3d = []


            ## loop over 2d geometries touching current tile

            for index, row in map2d_unfiltered[:].iterrows():

                if row["geometry"].centroid.intersects(tile.geom):

                    print(index, map2d_unfiltered.shape[0])
                    bodies_feature = []
                    bodies_feature_transformed = []              
                    multipolygons_feature = []

                    lowestpoints = []
                    lowestpoints_geo = []
                    lowestpoints_target = []
                    z_dist_max = 0


                    ## loop over 3d bodies to find the ones touching current 2d geometry

                    for body_id,(multipolygon,body) in enumerate(zip(bodies_multipolygons,bodies)):

                        if row["geometry"].buffer(1.0).intersects(unary_union(multipolygon)):

                            vertices = body.vertices

                            for vertex_id, vertex in enumerate(vertices):

                                if lowestpoints == [] or round(vertex[2],10) <= round(lowestpoints[0][2],10):
                                    if lowestpoints == [] or round(vertex[2],10) < round(lowestpoints[0][2],10):
                                        lowestpoints_geo = []
                                        lowestpoints = []
                                        lowestpoints_target = []

                                    lowestpoints_geo.append(Point(vertex[0], vertex[1]))
                                    lowestpoints.append((vertex[0], vertex[1], vertex[2]))
                                    lowestpoints_target.append(((vertex[0] - center[0]) * scale[0], (vertex[1] - center[1]) * scale[1], (vertex[2] - center[2]) * scale[2]))


                            lowestpoints_geo_mp = MultiPoint(lowestpoints_geo)
                            bodies_feature.append(body)
                            multipolygons_feature.append(multipolygon)


                    ## loop over adjectent tile to find terrain z-coordinate corresponding to lowest feature coordinate

                    for tile_adjacent in tilingscheme.tiles:

                        glb_filepath = str(Path(out_dirpath, "terrain" + "__" + str(tile_adjacent.id[0]) + "_" + str(tile_adjacent.id[1]) + ".glb"))

                        if os.path.isfile(glb_filepath) and lowestpoints_geo_mp.intersects(tile_adjacent.geom):
                            glb_reader = vtk.vtkGLTFReader()
                            glb_reader.SetFileName(glb_filepath)
                            glb_reader.Update()
                            polydata = vtk.vtkCompositeDataGeometryFilter()
                            polydata.SetInputConnection(glb_reader.GetOutputPort())
                            polydata.Update()
                            stl = polydata.GetOutput()

                            for lowestpoint, lowestpoint_target in zip(lowestpoints, lowestpoints_target):

                                z_values = self._get_mesh_elevation_from_xy(stl, (lowestpoint_target[0], lowestpoint_target[1]))
                                print(z_values, lowestpoint_target, stl.GetBounds(), glb_filepath, stl.GetBounds()[0] < lowestpoint_target[0] < stl.GetBounds()[1], stl.GetBounds()[2] < lowestpoint_target[1] < stl.GetBounds()[3])

                                if z_values:
                                    z_value_min = min(z_values)
                                    z_dist = lowestpoint_target[2] - z_value_min
                                    if not z_dist_max or z_dist > z_dist_max:
                                        z_dist_max = z_dist


                    ## loop over 3d bodies touching current 2d geometry and transform vertices in x,y, and z

                    for body in bodies_feature:
                        body_transformed = copy.deepcopy(body)
                        vertices_transformed = []
                        vertices = body.vertices

                        for vertex_id, vertex in enumerate(vertices):
                            vertices_transformed.append([(vertex[0] - center[0]) * scale[0], (vertex[1] - center[1]) * scale[1], ((vertex[2] - center[2]) * scale[2]) - z_dist_max])

                        body_transformed.vertices = vertices_transformed       
                        bodies_feature_transformed.append(body_transformed)



                    ## create proper features from body 2d representation (multipolygon) and 3d bodies

                    feature2d = unary_union(multipolygons_feature)

                    features2d.append(unary_union(multipolygons_feature))
                    features3d.append(trimesh.util.concatenate(bodies_feature_transformed))
                    logging.info(f"Feature3d added {len(features3d)}")

            return features2d, features3d


        def _get_mesh_elevation_from_xy(self, mesh, coord):

            x_coord, y_coord = coord
            pSource = [x_coord, y_coord, -9999999999]
            pTarget = [x_coord, y_coord, 9999999999]

            obbTree = vtk.vtkOBBTree()
            obbTree.SetDataSet(mesh)
            obbTree.BuildLocator()

            pointsVTKintersection = vtk.vtkPoints()
            code = obbTree.IntersectWithLine(pSource, pTarget, pointsVTKintersection, None)

            pointsVTKIntersectionData = pointsVTKintersection.GetData()
            noPointsVTKIntersection = pointsVTKIntersectionData.GetNumberOfTuples()
            pointsIntersection = []
            for idx in range(noPointsVTKIntersection):
                _tup = pointsVTKIntersectionData.GetTuple3(idx)
                pointsIntersection.append(_tup)

            if not len(pointsIntersection) == 0:
                return [pointsIntersection[0][2], pointsIntersection[1][2]]
            else:
                return None


        def _poly_from_triangle(self, vertices, face):

            points = [
                [vertices[face[0]][0], vertices[face[0]][1]],
                [vertices[face[1]][0], vertices[face[1]][1]],
                [vertices[face[2]][0], vertices[face[2]][1]]
            ]

            face_poly_mp = MultiPoint(points)
            face_poly = face_poly_mp.convex_hull
            face_poly = make_valid(face_poly)

            if face_poly.geom_type == "Polygon":
                face_poly_out = face_poly
            elif face_poly.geom_type in ["LineString", "Point"]:
                face_poly_out = face_poly.buffer(0.1)
            else:
                face_poly_out = None

            return face_poly_out


    class Buildings(Features):
        def __init__(
            self,
            description="buildings",
            tilingscheme=None,
            out_dirpath=tempfile.gettempdir(),
            filepaths=None,
            recombine_bodies=False,
            simplify_factor=None,
            tiles=None,
            boundary=None,
        ):

            for tile in tiles:

                self.process(
                    tilingscheme=tilingscheme,
                    out_dirpath=out_dirpath,
                    filepaths=filepaths,
                    recombine_bodies=recombine_bodies,
                    simplify_factor=simplify_factor,
                    tile=tile,
                    description=description,
                    boundary=boundary,
                )

                tile.total_bounds = self.map2d.total_bounds

                with open(Path(out_dirpath, description + "__" + '_'.join(str(x) for x in tile.id) + ".glb"), 'wb') as f:
                    f.write(trimesh.exchange.gltf.export_glb(self.scene, include_normals=True))



    def get_dem_information(self, dirpath):
        return [None, None, None, None]


    def define_tiles(self, boundary, origin, height=10, width=10):
        pass




