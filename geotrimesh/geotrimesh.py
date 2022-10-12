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
import warnings
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

logging.basicConfig(level=logging.INFO)
warnings.filterwarnings("ignore")


def _get_mesh_from_scene(scene):

    meshes = []

    for key in scene.geometry:
        if type(scene.geometry[key]).__name__ == "Trimesh":
            mesh_orig = scene.geometry[key]
            mesh_orig_new = trimesh.Trimesh(
                vertices=mesh_orig.vertices,
                faces=mesh_orig.faces,
                face_normals=mesh_orig.face_normals,
                vertex_normals=mesh_orig.vertex_normals,
            )
            meshes.append(mesh_orig_new)

    return trimesh.util.concatenate(meshes)


class GeoSceneSet:
    def __init__(self):
        self.terrain = None
        self.features = None
        self.ortho = None
        self.dem = None

    class TilingScheme:
        def __init__(self, boundary, filepaths, height=None, width=None):

            logging.info("Generating tiling scheme")

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

            center = (
                (float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0,
                (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0,
            )

            # scale = (100000.0, 100000.0, 100000.0)
            scale = (1.0, 1.0, 1.0)

            print(upper_left, filepaths)
            upper_left_aoi, res, nodata = self.define_aoi_origin(upper_left, filepaths)

            tile = namedtuple(
                "tile", ["id", "scheme_dims", "dims", "res", "nodata", "geom"]
            )

            if not (width and height):

                height = math.ceil((upper_left.y - lower_right.y) / res[1]) + 1
                width = math.ceil((lower_right.x - upper_left.x) / res[0]) + 1

                geom = box(
                    upper_left_aoi.x,
                    upper_left_aoi.y - (height * res[1]),
                    upper_left_aoi.x + (width * res[0]),
                    upper_left_aoi.y,
                )

                geom_target = box(
                    (upper_left_aoi.x - center[0]) * scale[0],
                    ((upper_left_aoi.y - (height * res[1])) - center[1]) * scale[1],
                    ((upper_left_aoi.x + (width * res[0])) - center[0]) * scale[0],
                    (upper_left_aoi.y - center[1]) * scale[1],
                )

                gdf_dict["geometry"].append(geom)

                tile = namedtuple(
                    "tile", ["id", "scheme_dims", "dims", "res", "nodata", "geom"]
                )

                tile.id = (0, 0)
                tile.scheme_dims = (1, 1)
                tile.dims = (height, width)
                tile.res = res
                tile.nodata = nodata
                tile.geom = geom
                tile.geom_target = geom_target
                tile.total_bounds = tile.geom.bounds
                tile.valid = True

                self.tiles.append(tile)

            else:

                count_x = math.ceil(
                    ((lower_right.x - upper_left_aoi.x) / res[0]) / width
                )
                count_y = math.ceil(
                    ((upper_left_aoi.y - lower_right.y) / abs(res[1])) / height
                )

                for y in range(0, count_y):
                    for x in range(0, count_x):

                        geom = box(
                            upper_left_aoi.x + (x * width * res[0]),
                            upper_left_aoi.y - ((y + 1) * height * res[1]),
                            upper_left_aoi.x + ((x + 1) * width * res[0]),
                            upper_left_aoi.y - (y * height * res[1]),
                        )

                        geom_target = box(
                            (upper_left_aoi.x - center[0]) * scale[0],
                            ((upper_left_aoi.y - (height * res[1])) - center[1])
                            * scale[1],
                            ((upper_left_aoi.x + (width * res[0])) - center[0])
                            * scale[0],
                            (upper_left_aoi.y - center[1]) * scale[1],
                        )

                        gdf_dict["geometry"].append(geom)

                        tile = namedtuple(
                            "tile",
                            ["id", "scheme_dims", "dims", "res", "nodata", "geom"],
                        )

                        tile.id = (y, x)
                        tile.scheme_dims = (count_y, count_x)
                        tile.dims = (height, width)
                        tile.res = res
                        tile.nodata = nodata
                        tile.geom = geom
                        tile.geom_target = geom_target
                        tile.total_bounds = tile.geom.bounds
                        tile.valid = True


                        self.tiles.append(tile)

            self.gdf = gpd.GeoDataFrame.from_dict(gdf_dict)
            self.gdf = self.gdf.set_crs(boundary.crs, allow_override=True)

        def define_aoi_origin(self, point, filepaths):

            for filepath in filepaths:
                src = rasterio.open(str(filepath))
                src_box = box(*src.bounds)

                if point.intersects(src_box):
                    dist_left_pix = math.floor((point.x - src.bounds.left) / src.res[0])
                    dist_top_pix = math.floor((src.bounds.top - point.y) / src.res[1])

                    origin = Point(
                        src.bounds.left + (dist_left_pix * src.res[0]),
                        src.bounds.top - (dist_top_pix * src.res[1]),
                    )

                    if src.nodatavals[0]:
                        nodata = src.nodatavals[0]
                    else:
                        nodata = -99999

                    return (origin, src.res, nodata)

    class Ortho:
        def __init__(
            self,
            tilingscheme=None,
            out_dirpath=tempfile.gettempdir(),
            filepaths=None,
            tiles=None,
            boundary=None,
            y_up=True,
            include_texture=False
        ):

            logging.info("Setting texture coordinates")

            for tile in tiles:

                if tile.valid:
                    for glb_id, glb in enumerate([
                        [
                            "terrain",
                            str(
                                Path(
                                    out_dirpath,
                                    "terrain"
                                    + "__"
                                    + str(tile.id[0])
                                    + "_"
                                    + str(tile.id[1])
                                    + ".glb",
                                )
                            ),
                        ],
                        [
                            "buildings",
                            str(
                                Path(
                                    out_dirpath,
                                    "buildings"
                                    + "__"
                                    + str(tile.id[0])
                                    + "_"
                                    + str(tile.id[1])
                                    + ".glb",
                                )
                            ),
                        ],
                    ]):

                        scene = trimesh.load(str(glb[1]))




                        scene_out = trimesh.Scene()



                        mosaic, width, height = self.generate_ortho_mosaic(
                            filepaths, tile
                        )

                        mosaic_rearranged = np.transpose(mosaic, axes=[1, 2, 0])

                        newImg1 = Image.fromarray(mosaic_rearranged.astype("uint8"), "RGB")
                        newImg1 = ImageOps.flip(newImg1)

                        newImg1.save(
                            Path(
                                out_dirpath,
                                "ortho__"
                                + str(tile.id[0])
                                + "_"
                                + str(tile.id[1])
                                + ".png",
                            ),
                            "PNG",
                        )


                        if include_texture:
                            im = Image.open(
                                Path(
                                    out_dirpath,
                                    "ortho__"
                                    + str(tile.id[0])
                                    + "_"
                                    + str(tile.id[1])
                                    + ".png",
                                )
                            )                   
                        else:
                            im = None




                        meshes = []
                        center = (
                            (float(boundary.bounds.minx) + float(boundary.bounds.maxx))
                            / 2.0,
                            (float(boundary.bounds.miny) + float(boundary.bounds.maxy))
                            / 2.0,
                        )

                        # scale = (100000.0, 100000.0, 100000.0)
                        scale = (1.0, 1.0, 1.0)

                        left, bottom, right, top = tile.total_bounds
                        #left -= tile.res[1] * 1.0
                        #bottom -= tile.res[0] * 0.5
                        #right += tile.res[1] * 0.5
                        #top += tile.res[0] * 1.0
                        print("[", tile.id, "]", "total_bounds", left, bottom, right, top)


                        #print((tile.total_bounds[2] - tile.total_bounds[0]) / tile.res[0])
                        #print((tile.total_bounds[3] - tile.total_bounds[1]) / tile.res[1])
                        #sys.exit()


                        x_ratio_min = None
                        x_ratio_max = None
                        y_ratio_min = None
                        y_ratio_max = None

                        vertices_nr = []
                        for key_id,key in enumerate(scene.geometry):

                            if type(scene.geometry[key]).__name__ == "Trimesh":
                                mesh_orig = scene.geometry[key]
                                vertices_nr.append(len(mesh_orig.vertices))

                                uv = []
                                vertices_switched_axes = []
                                vertex_normals_switched_axes = []
                                for vert, vertex_normal in zip(mesh_orig.vertices,mesh_orig.vertex_normals):
                                    vert_x_local, vert_y_local, vert_z_local = vert
                                    vertices_switched_axes.append(
                                        (vert_x_local, vert_z_local, vert_y_local*-1)
                                    )

                                    x = (float(vert_x_local) / scale[0]) + center[0]
                                    y = (float(vert_y_local) / scale[1]) + center[1]



                                    x_offset = x - left# + tile.res[0]
                                    y_offset = top - y

                                    dem_x_dist = right - left
                                    dem_y_dist = top - bottom

                                    x_ratio = round(x_offset / float(dem_x_dist), 10)
                                    y_ratio = round(y_offset / float(dem_y_dist), 10)

                                    if not x_ratio_min or x_ratio < x_ratio_min:
                                        x_ratio_min = x_ratio
                                    if not x_ratio_max or x_ratio > x_ratio_max:
                                        x_ratio_max = x_ratio

                                    if not y_ratio_min or y_ratio < y_ratio_min:
                                        y_ratio_min = y_ratio
                                    if not y_ratio_max or y_ratio > y_ratio_max:
                                        y_ratio_max = y_ratio



                                    uv.append(np.array([x_ratio, y_ratio]))

                                material = trimesh.visual.texture.SimpleMaterial(image=im)


                                color_visuals = trimesh.visual.TextureVisuals(
                                    uv=uv, image=im, material=material
                                )


                                mesh_orig.visual=color_visuals
                                mesh_orig.vertices=vertices_switched_axes




                                clean_output = False
                                if clean_output:

                                    open3d_mesh = open3d.geometry.TriangleMesh()
                                    open3d_mesh.vertices = open3d.utility.Vector3dVector(
                                        mesh_orig.vertices
                                    )
                                    open3d_mesh.triangles = open3d.utility.Vector3iVector(
                                        mesh_orig.faces
                                    )

                                    open3d_mesh.orient_triangles()
                                    open3d_mesh.compute_triangle_normals(normalized=True)
                                    open3d_mesh.compute_vertex_normals(normalized=True)
                                    open3d_mesh.remove_degenerate_triangles()
                                    open3d_mesh.remove_duplicated_triangles()
                                    open3d_mesh.remove_duplicated_vertices()
                                    open3d_mesh.remove_unreferenced_vertices()

                                    mesh_orig.vertices=open3d_mesh.vertices
                                    mesh_orig.faces=open3d_mesh.triangles
                                    mesh_orig.face_normals=open3d_mesh.triangle_normals
                                    mesh_orig.vertex_normals=mesh_orig.face_normals



                                ## Force viewers to use polygon normals (by removing vertex normals)
                                vertex_normals_new = []
                                for normal in mesh_orig.vertex_normals:
                                    vertex_normals_new.append([0.0, 0.0, 0.0])
                                mesh_orig.vertex_normals = vertex_normals_new


                                scene_out.add_geometry(
                                    mesh_orig, node_name=str(key_id), geom_name=str(key_id)
                                )


                        with open(
                            Path(
                                out_dirpath,
                                "result_"
                                + glb[0]
                                + "__"
                                + str(tile.id[0])
                                + "_"
                                + str(tile.id[1])
                                + ".glb",
                            ),
                            "wb",
                        ) as f:
                            f.write(
                                trimesh.exchange.gltf.export_glb(
                                    scene_out, include_normals=True
                                )
                            )



        def generate_ortho_mosaic(self, filepaths, tile):

            left, bottom, right, top = tile.total_bounds

            with rasterio.open(str(filepaths[0])) as src:
                res_x, res_y = src.res

            with rasterio.open(filepaths[0]) as src:
                target_array = src.read(
                    1,
                    window=from_bounds(left, bottom, right, top, src.transform),
                    boundless=True,
                    fill_value=0,
                )
                rast = np.zeros(
                    (3, target_array.shape[0], target_array.shape[1]), dtype=np.uint8
                )

            for filepath in filepaths:
                loc = np.amax(rast, axis=0) == 0
                with rasterio.open(filepath) as src:
                    for i in range(0, 3):
                        rast[i][loc] = src.read(
                            i + 1,
                            window=from_bounds(left, bottom, right, top, src.transform),
                            boundless=True,
                            fill_value=0,
                        )[loc]

            return rast, target_array.shape[0], target_array.shape[1]

    class Dem:
        def __init__(self, filepaths, bounds, tile=None):
            return None

    class Terrain:
        def __init__(
            self,
            out_dirpath=tempfile.gettempdir(),
            tiles=None,
            filepaths=[],
            tile=None,
            boundary=None,
            openscad_bin_filepath=Path("openscad"),
        ):

            logging.info("Processing Terrain")

            for tile in tiles:

                self.mosaic, self.z_bounds_total = self.generate_terrain_mosaic(
                    filepaths, tile
                )

                self.mesh = self.generate_terrain_mesh(
                    self.mosaic,
                    boundary,
                    self.z_bounds_total,
                    tile,
                    out_dirpath,
                    openscad_bin_filepath,
                )


                self.scene = trimesh.Scene()
                self.scene.add_geometry(self.mesh, node_name="mesh", geom_name="mesh")

                with open(
                    Path(
                        out_dirpath,
                        "terrain__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".glb",
                    ),
                    "wb",
                ) as f:
                    f.write(
                        trimesh.exchange.gltf.export_glb(
                            self.scene, include_normals=True
                        )
                    )

        def generate_terrain_mosaic(self, filepaths, tile):

            #print(tile.id)
            rast = np.ones((tile.dims[0]+1, tile.dims[1]+1), dtype=np.float32) * tile.nodata
            z_min_total = None
            z_max_total = None

            for filepath in filepaths:
                print(filepath, tile.geom.bounds, tile.res)
                with rasterio.open(filepath) as src:
                    left, bottom, right, top = tile.geom.bounds
                    #left -= tile.res[1]
                    bottom -= tile.res[0]
                    right += tile.res[1]
                    #top += tile.res[0]


                    meta = src.read(1)

                    rast_tmp = src.read(
                        1,
                        window=from_bounds(left, bottom, right, top, src.transform),
                        boundless=True,
                        fill_value=src.nodatavals[0],
                    )

                    #print(left, bottom, right, top)
                    #print(rast_tmp)


                    rast[rast==tile.nodata] = rast_tmp[rast==tile.nodata]


                    z_min = np.nanmin(meta[meta != src.nodatavals[0]])
                    z_max = np.nanmax(meta[meta != src.nodatavals[0]])

                    if not z_min_total or z_min < z_min_total:
                        z_min_total = z_min
                    if not z_max_total or z_max > z_max_total:
                        z_max_total = z_max

            return rast, (z_min_total, z_max_total)

        def _get_polyhedron_faces(self, rast, tile, boundary, z_bounds):

            left, bottom, right, top = tile.geom.bounds
            left -= tile.res[1]
            bottom -= tile.res[0]
            right += tile.res[1]
            top += tile.res[0]


            center = (
                (float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0,
                (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0,
                (float(z_bounds[0]) + float(z_bounds[1])) / 2.0,
            )

            scale = (1.0, 1.0, 1.0)
            # scale = (100000.0, 100000.0, 100000.0)

            polyhedron_faces_array = np.zeros(
                ((tile.dims[0]+1), (tile.dims[1]+1), 16, 3), dtype=np.int32
            )
            polyhedron_points_floor_array = np.zeros(
                ((tile.dims[0]+1) * (tile.dims[1]+1), 3), dtype=np.float32
            )

            cnt = 0
            dem_x_min = None
            dem_x_max = None
            dem_y_min = None
            dem_y_max = None

            for i in range(0, tile.dims[0]+1):

                #logging.info(f"Row: {i} of {tile.dims[0]}")

                for j in range(0, tile.dims[1]+1):

                    i0_coord = top - (tile.res[1] * i)
                    j0_coord = left + (tile.res[0] * j)

                    z_a = rast[i][j] * scale[2]

                    dem_x = j0_coord
                    dem_y = i0_coord

                    polyhedron_points_floor_array[(i * (tile.dims[1]+1)) + j][0] = (
                        dem_x - center[0]
                    ) * scale[0]
                    polyhedron_points_floor_array[(i * (tile.dims[1]+1)) + j][1] = (
                        dem_y - center[1]
                    ) * scale[1]
                    polyhedron_points_floor_array[(i * (tile.dims[1]+1)) + j][2] = (
                        z_a - center[2]
                    )

                    if not dem_x_min or dem_x < dem_x_min:
                        dem_x_min = dem_x
                    if not dem_x_max or dem_x > dem_x_max:
                        dem_x_max = dem_x

                    if not dem_y_min or dem_y < dem_y_min:
                        dem_y_min = dem_y
                    if not dem_y_max or dem_y > dem_y_max:
                        dem_y_max = dem_y

                    if i < (tile.dims[0]+1) - 1 and j < (tile.dims[1]+1) - 1:

                        z_b = rast[i + 1][j] * scale[2]
                        z_c = rast[i][j + 1] * scale[2]
                        z_d = rast[i + 1][j + 1] * scale[2]

                        point_a_ceil = (i * (tile.dims[1]+1)) + j
                        point_b_ceil = ((i + 1) * (tile.dims[1]+1)) + j
                        point_c_ceil = (i * (tile.dims[1]+1)) + j + 1
                        point_d_ceil = ((i + 1) * (tile.dims[1]+1)) + j + 1

                        point_a_floor = (
                            ((tile.dims[0]+1) * (tile.dims[1]+1)) + (i * (tile.dims[1]+1)) + j
                        )
                        point_b_floor = (
                            ((tile.dims[0]+1) * (tile.dims[1]+1)) + ((i + 1) * (tile.dims[1]+1)) + j
                        )
                        point_c_floor = (
                            ((tile.dims[0]+1) * (tile.dims[1]+1)) + (i * (tile.dims[1]+1)) + j + 1
                        )
                        point_d_floor = (
                            ((tile.dims[0]+1) * (tile.dims[1]+1))
                            + ((i + 1) * (tile.dims[1]+1))
                            + j
                            + 1
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

        def _get_polyhedron_faces_clean(
            self, polyhedron_faces_array, polyhedron_points_floor_array, tile
        ):

            polyhedron_faces_clean = []

            for i in range(0, tile.dims[0]+1):

                for j in range(0, tile.dims[1]+1):

                    for polyhedron_face_id in range(0, 16):

                        polyhedron_face = polyhedron_faces_array[i][j][
                            polyhedron_face_id
                        ][:].tolist()
                        polyhedron_face_cnt = 0

                        if not (
                            polyhedron_face[0] == 0
                            and polyhedron_face[1] == 0
                            and polyhedron_face[2] == 0
                        ):

                            i_bottom = i - 1 if i > 0 else 0
                            i_top = i + 1 if i < (tile.dims[0]+1) - 1 else (tile.dims[0]+1) - 1
                            j_left = j - 1 if i > 0 else 0
                            j_right = j + 1 if j < (tile.dims[1]+1) - 1 else (tile.dims[1]+1)

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
                                            polyhedron_faces_array[i_neighbour][
                                                j_neighbour
                                            ][:]
                                            == array1,
                                            axis=1,
                                        )
                                    )
                                    loc2 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][
                                                j_neighbour
                                            ][:]
                                            == array2,
                                            axis=1,
                                        )
                                    )
                                    loc3 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][
                                                j_neighbour
                                            ][:]
                                            == array3,
                                            axis=1,
                                        )
                                    )
                                    loc4 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][
                                                j_neighbour
                                            ][:]
                                            == array4,
                                            axis=1,
                                        )
                                    )
                                    loc5 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][
                                                j_neighbour
                                            ][:]
                                            == array5,
                                            axis=1,
                                        )
                                    )
                                    loc6 = np.where(
                                        np.all(
                                            polyhedron_faces_array[i_neighbour][
                                                j_neighbour
                                            ][:]
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

        def _get_polyhedron_points(
            self,
            polyhedron_points_floor_array,
            dem_extrude_height_up,
            dem_extrude_height_down,
        ):
            polyhedron_points = []

            for l in range(1, -1, -1):
                for polyhedron_point_id in range(
                    0, polyhedron_points_floor_array.shape[0]
                ):

                    if l == 1:
                        z_offset = dem_extrude_height_up
                    if l == 0:
                        z_offset = abs(dem_extrude_height_down) * -1

                    polyhedron_point = [
                        polyhedron_points_floor_array[polyhedron_point_id][0],
                        polyhedron_points_floor_array[polyhedron_point_id][1],
                        polyhedron_points_floor_array[polyhedron_point_id][2]
                        - z_offset,
                    ]

                    polyhedron_points.append(polyhedron_point)

            return polyhedron_points

        def _get_polygon_intersect(
            self,
            boundary=None,
            tile=None,
            polygon_extrude_height_down=None,
            polygon_extrude_height_up=None,
            center=None,
            scale=None,
        ):

            center = (
                (float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0,
                (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0,
            )

            left, bottom, right, top = tile.geom.bounds

            scale = (1.0, 1.0, 1.0)
            # scale = (100000.0, 100000.0, 100000.0)

            boundary_tile = boundary.cx[left:right, bottom:top]
            polys = []

            if boundary_tile.shape[0] > 0:

                for index, row in boundary_tile.iterrows():
                    geom = row["geometry"]

                    polygon_points = []
                    polygon_paths = []
                    coord_id = 0
                    point_cnt = 0

                    polygon_path = []
                    for coord_x, coord_y in zip(
                        geom.exterior.coords.xy[0], geom.exterior.coords.xy[1]
                    ):
                        polygon_path.append(coord_id)
                        polygon_points.append(
                            [
                                (coord_x - center[0]) * scale[0],
                                (coord_y - center[1]) * scale[1],
                            ]
                        )
                        coord_id += 1
                    polygon_paths.append(polygon_path)

                    if len(geom.interiors) > 0:
                        for interior in geom.interiors:
                            polygon_path = []
                            for coord_x, coord_y in zip(
                                interior.coords.xy[0], interior.coords.xy[1]
                            ):
                                polygon_path.append(coord_id)
                                polygon_points.append(
                                    [
                                        (coord_x - center[0]) * scale[0],
                                        (coord_y - center[1]) * scale[1],
                                    ]
                                )
                                coord_id += 1
                            polygon_paths.append(polygon_path)

                    polys.append([polygon_paths, polygon_points])

            return polys

        def _write_scad(
            self,
            polyhedron_points,
            polyhedron_faces_clean,
            polys,
            polygon_extrude_height_down,
            polygon_extrude_height_up,
            scad_filepath,
        ):

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
                        abs(polygon_extrude_height_down)
                        + abs(polygon_extrude_height_up),
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
            subprocess_commands = [
                str(openscad_bin_filepath),
                scad_filepath,
                "-o",
                stl_filepath,
            ]
            output = subprocess.check_output(subprocess_commands, shell=False)

            stl = trimesh.load(str(stl_filepath))
            return stl

        def generate_terrain_mesh(
            self,
            rast=None,
            boundary=None,
            z_bounds=None,
            tile=None,
            out_dirpath=tempfile.gettempdir(),
            openscad_bin_filepath=Path("openscad"),
        ):

            # scale = (100000.0, 100000.0, 100000.0)
            scale = (1.0, 1.0, 1.0)

            import logging

            polygon_extrude_height_down = 1000 * scale[2]
            polygon_extrude_height_up = 1000 * scale[2]
            dem_extrude_height_down = 0.5 * scale[2]
            dem_extrude_height_up = 0.5 * scale[2]

            polyhedron_faces_array, polyhedron_points_floor_array = self._get_polyhedron_faces(
                rast, tile, boundary, z_bounds
            )
            polyhedron_faces_clean = self._get_polyhedron_faces_clean(
                polyhedron_faces_array, polyhedron_points_floor_array, tile
            )

            polyhedron_points = self._get_polyhedron_points(
                polyhedron_points_floor_array,
                dem_extrude_height_up,
                dem_extrude_height_down,
            )

            polys = self._get_polygon_intersect(
                boundary, tile, polygon_extrude_height_down, polygon_extrude_height_up
            )

            if len(polys) > 0:

                tile.valid = True
                scad_filepath = str(
                    Path(
                        out_dirpath,
                        "terrain"
                        + "__"
                        + str(tile.id[0])
                        + "_"
                        + str(tile.id[1])
                        + ".scad",
                    )
                )
                stl_filepath = str(
                    Path(
                        out_dirpath,
                        "terrain" + "__" + str(tile.id[0]) + "_" + str(tile.id[1]) + ".stl",
                    )
                )

                self._write_scad(
                    polyhedron_points,
                    polyhedron_faces_clean,
                    polys,
                    polygon_extrude_height_down,
                    polygon_extrude_height_up,
                    scad_filepath,
                )

                mesh = self._generate_stl(
                    openscad_bin_filepath, scad_filepath, stl_filepath
                )

                return mesh

            else:
                tile.valid = False
                return None



    class Features:
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
            offset=None,
            margin = 150,
            number_of_neighbours=3
        ):

            logging.info("Processing Features")

            def _get_mesh_elevation_from_xy(mesh, coord):

                x_coord, y_coord = coord
                pSource = [x_coord, y_coord, -9999999999]
                pTarget = [x_coord, y_coord, 9999999999]

                obbTree = vtk.vtkOBBTree()
                obbTree.SetDataSet(mesh)
                obbTree.BuildLocator()

                pointsVTKintersection = vtk.vtkPoints()
                code = obbTree.IntersectWithLine(
                    pSource, pTarget, pointsVTKintersection, None
                )

                pointsVTKIntersectionData = pointsVTKintersection.GetData()
                noPointsVTKIntersection = pointsVTKIntersectionData.GetNumberOfTuples()
                pointsIntersection = []
                for idx in range(noPointsVTKIntersection):
                    _tup = pointsVTKIntersectionData.GetTuple3(idx)
                    pointsIntersection.append(_tup)

                if not len(pointsIntersection) == 0:
                    return max([pointsIntersection[0][2], pointsIntersection[1][2]])
                else:
                    return None

            def translate_vertices(mesh, offset):

                x_offset, y_offset, z_offset = offset

                vertices = mesh.vertices

                vertices_translated = []

                for vertex in vertices:
                    x, y, z = vertex
                    vertices_translated.append(
                        [x + x_offset, y + y_offset, z + z_offset]
                    )

                    #print((x, y, z), (x+x_offset, y+y_offset, z+z_offset))


                mesh_new = trimesh.Trimesh(
                    vertices=vertices_translated,
                    faces=mesh.faces,
                    face_normals=mesh.face_normals,
                    vertex_normals=mesh.vertex_normals,
                )

                return mesh_new

            def get_stats(mesh):

                z_min = None
                z_max = None
                x_min = None
                x_max = None
                y_min = None
                y_max = None

                faces = mesh.faces
                vertices = mesh.vertices

                for face in faces:

                    points = [
                        [
                            vertices[face[0]][0],
                            vertices[face[0]][1],
                            vertices[face[0]][2],
                        ],
                        [
                            vertices[face[1]][0],
                            vertices[face[1]][1],
                            vertices[face[1]][2],
                        ],
                        [
                            vertices[face[2]][0],
                            vertices[face[2]][1],
                            vertices[face[2]][2],
                        ],
                    ]

                    for point in points:

                        z_max = max(z_max, point[2]) if z_max else point[2]

                        if not z_min or (z_min and point[2] < z_min):
                            z_min = point[2]
                            z_min_x = point[0]
                            z_min_y = point[1]

                return (z_min_x, z_min_y, z_min)

            def filter_faces_by_geom(mesh, geom):

                vertices = mesh.vertices
                faces = mesh.faces

                face_mask = []
                for face in faces:

                    points = [
                        [vertices[face[0]][0], vertices[face[0]][1]],
                        [vertices[face[1]][0], vertices[face[1]][1]],
                        [vertices[face[2]][0], vertices[face[2]][1]],
                    ]

                    inside = True
                    for point in points:

                        p = Point(point[0], point[1])
                        if not p.intersects(geom):
                            inside = False

                    #if inside:
                    #    print(inside, points, geom)
                    #sys.exit()
                    face_mask.append(inside)

                mesh_out = copy.deepcopy(mesh)
                mesh_out.update_faces(face_mask)

                return mesh_out

            def get_2d_representation(meshes, boundary):

                geoms = []
                for mesh in meshes:

                    faces = mesh.faces
                    vertices = mesh.vertices

                    for face in faces:

                        points = [
                            (vertices[face[0]][0], vertices[face[0]][1]),
                            (vertices[face[1]][0], vertices[face[1]][1]),
                            (vertices[face[2]][0], vertices[face[2]][1]),
                        ]

                        poly = make_valid(MultiPoint(points).convex_hull)
                        if poly.geom_type == "Polygon":
                            geoms.append(poly)

                gdf = gpd.GeoDataFrame.from_dict({"geometry": unary_union(geoms)})

                return gdf

            center = (
                (float(boundary.bounds.minx) + float(boundary.bounds.maxx)) / 2.0,
                (float(boundary.bounds.miny) + float(boundary.bounds.maxy)) / 2.0,
                0,
            )

            for tile in tiles:

                print(tile.id)
                if True: #tile.valid and not (str(tile.total_bounds[0]) == "nan" ):

                    features_out_filepath = str(
                        Path(
                            out_dirpath,
                            description
                            + "__"
                            + str(tile.id[0])
                            + "_"
                            + str(tile.id[1])
                            + ".glb",
                        )
                    )

                    map2d_out_filepath = str(
                        Path(
                            out_dirpath,
                            description
                            + "__"
                            + str(tile.id[0])
                            + "_"
                            + str(tile.id[1])
                            + ".gpkg",
                        )
                    )
                    map2d_total_bounds_out_filepath = str(
                        Path(
                            out_dirpath,
                            description + "_total_bounds"
                            + "__"
                            + str(tile.id[0])
                            + "_"
                            + str(tile.id[1])
                            + ".gpkg",
                        )
                    )

                    map2d_adjacent_out_filepath = str(
                        Path(
                            out_dirpath,
                            description
                            + "_adjacent"
                            + "__"
                            + str(tile.id[0])
                            + "_"
                            + str(tile.id[1])
                            + ".gpkg",
                        )
                    )

                    adjacent_tiles = []
                    adjacent_tiles_x_min = tile.geom.bounds[0] - margin
                    adjacent_tiles_y_min = tile.geom.bounds[1] - margin
                    adjacent_tiles_x_max = tile.geom.bounds[2] + margin
                    adjacent_tiles_y_max = tile.geom.bounds[3] + margin


                    adjacent_tiles = []

                    for adjacent_tile_id_y in range(
                        tile.id[0] - number_of_neighbours,
                        tile.id[0] + number_of_neighbours + 1,
                    ):
                        for adjacent_tile_id_x in range(
                            tile.id[1] - number_of_neighbours,
                            tile.id[1] + number_of_neighbours + 1,
                        ):

                            for adjacent_tile in tiles:

                                if adjacent_tile.id == (
                                    adjacent_tile_id_y,
                                    adjacent_tile_id_x
                                ):

                                    adjacent_tiles.append(adjacent_tile)

                    adjacent_box = box(
                        adjacent_tiles_x_min,
                        adjacent_tiles_y_min,
                        adjacent_tiles_x_max,
                        adjacent_tiles_y_max,
                    )

                    for filepath in filepaths:
                        
                        scene = trimesh.load(filepath)

                        #x_offset = center[0]
                        #y_offset = center[1]
                        #z_offset = 0

                        x_min_orig = 2677116.375000
                        y_min_orig = 1241839.025000
                        x_max_orig = 2689381.985000
                        y_max_orig = 1254150.950000

                        #x_offset = (x_min_orig + x_max_orig) / 2.0
                        #y_offset = (y_min_orig + y_max_orig) / 2.0
                        #z_offset = 0

                        xyz_min, xyz_max = scene.bounds
                        x_min, y_min, z_min = xyz_min
                        x_max, y_max, z_max = xyz_max

                        x_offset = x_min_orig - x_min
                        y_offset = y_min_orig - y_min
                        z_offset = 0


                        meshes = []
                        for key_id, key in enumerate(scene.geometry):
                            if type(scene.geometry[key]).__name__ == "Trimesh":
                                mesh = scene.geometry[key]
                                mesh_new = translate_vertices(mesh, (x_offset, y_offset, z_offset))
                                mesh_new = filter_faces_by_geom(mesh_new, adjacent_box)
                                meshes.append(mesh_new)

                        #print(len(meshes))

                    map2d_adjacent = get_2d_representation(meshes, boundary)
                    if map2d_adjacent.shape[0] > 0:
                        map2d_adjacent = map2d_adjacent.set_crs(boundary.crs, allow_override=True)
                        map2d_adjacent.to_file(map2d_adjacent_out_filepath)


                    map2d_dict = {}
                    map2d_dict["geometry"] = []
                    map2d_total_bounds_dict = {}
                    map2d_total_bounds_dict["geometry"] = []


                    scene_out = trimesh.Scene()
                    for index, row in map2d_adjacent.iterrows():
                        for mesh_id, mesh in enumerate(meshes):
                            if len(mesh.faces) > 0:

                                print(row["geometry"])
                                if row["geometry"].centroid.intersects(tile.geom):

                                    mesh_component = filter_faces_by_geom(
                                        mesh, row["geometry"]
                                    )
                                    z_lowest = get_stats(mesh_component)

                                    z_lowest_terrain = None

                                    for adjacent_tile in adjacent_tiles:

                                        if Point(z_lowest[0], z_lowest[1]).intersects(
                                            adjacent_tile.geom
                                        ):

                                            print("adjacent", adjacent_tile.id)
                                            glb_filepath = str(
                                                Path(
                                                    out_dirpath,
                                                    "terrain"
                                                    + "__"
                                                    + str(adjacent_tile.id[0])
                                                    + "_"
                                                    + str(adjacent_tile.id[1])
                                                    + ".glb",
                                                )
                                            )

                                            glb_reader = vtk.vtkGLTFReader()
                                            glb_reader.SetFileName(glb_filepath)
                                            glb_reader.Update()
                                            polydata = vtk.vtkCompositeDataGeometryFilter()
                                            polydata.SetInputConnection(
                                                glb_reader.GetOutputPort()
                                            )
                                            polydata.Update()
                                            stl = polydata.GetOutput()

                                            ## check if source actually contains a mesh
                                            if stl.GetPoints():

                                                scale = (1.0, 1.0, 1.0)
                                                z_lowest_terrain = _get_mesh_elevation_from_xy(
                                                    stl,
                                                    (
                                                        (z_lowest[0] - center[0]) * scale[0],
                                                        (z_lowest[1] - center[1]) * scale[1],
                                                    ),
                                                )

                                                if not z_lowest_terrain:
                                                    print(adjacent_tile.id[0], adjacent_tile.id[1], "error at", z_lowest[0], z_lowest[1], adjacent_tile.geom)




                                    if z_lowest_terrain:

                                        x_offset = center[0] * -1
                                        y_offset = center[1] * -1
                                        z_offset = (z_lowest[2] - z_lowest_terrain) * -1

                                        mesh_component = translate_vertices(
                                            mesh_component, (x_offset, y_offset, z_offset)
                                        )

                                        scene_out.add_geometry(
                                            mesh_component,
                                            node_name=str(index),
                                            geom_name=str(index),
                                        )

                                        map2d_dict["geometry"].append(row["geometry"])




                    map2d = gpd.GeoDataFrame.from_dict(map2d_dict)

                    if map2d.shape[0] > 0:
                        map2d = map2d.set_crs(boundary.crs, allow_override=True)
                        map2d.to_file(map2d_out_filepath)

                    #print(tile.total_bounds)
                    #print(map2d.total_bounds)

                    map2d_left, map2d_bottom, map2d_right, map2d_top = map2d.total_bounds

                    if not np.isnan(map2d_left):

                        if map2d_left < tile.total_bounds[0]:
                            map2d_left = tile.total_bounds[0] - (math.ceil((tile.total_bounds[0] - map2d_left) / tile.res[0]) * tile.res[0])

                        if map2d_bottom < tile.total_bounds[1]:
                            map2d_bottom = tile.total_bounds[1] - (math.ceil((tile.total_bounds[1] - map2d_bottom) / tile.res[1]) * tile.res[1])

                        if map2d_right > tile.total_bounds[2]:
                            map2d_right = tile.total_bounds[2] + (math.ceil((map2d_right - tile.total_bounds[2]) / tile.res[0]) * tile.res[0])

                        if map2d_top > tile.total_bounds[3]:
                            map2d_top = tile.total_bounds[3] + (math.ceil((map2d_top - tile.total_bounds[3]) / tile.res[1]) * tile.res[1])

                        #tile.total_bounds = (
                        #    min(map2d_left - (2*tile.res[0]), tile.total_bounds[0] - (2*tile.res[0])),
                        #    min(map2d_bottom - (2*tile.res[1]), tile.total_bounds[1] - (2*tile.res[1])),
                        #    max(map2d_right + (2*tile.res[0]), tile.total_bounds[2] + (2*tile.res[0])),
                        #    max(map2d_top + (2*tile.res[1]), tile.total_bounds[3] + (2*tile.res[1])),
                        #)

                        tile.total_bounds = (
                            min(map2d_left - tile.res[0], tile.total_bounds[0] - tile.res[0]),
                            min(map2d_bottom - tile.res[1], tile.total_bounds[1] - tile.res[1]),
                            max(map2d_right + tile.res[0], tile.total_bounds[2] + tile.res[0]),
                            max(map2d_top + tile.res[1], tile.total_bounds[3] + tile.res[1]),
                        )

                    print(tile.total_bounds[0], tile.total_bounds[1], tile.total_bounds[2], tile.total_bounds[3])
                    map2d_total_bounds_dict["geometry"].append(box(tile.total_bounds[0], tile.total_bounds[1], tile.total_bounds[2], tile.total_bounds[3]))
                    map2d_total_bounds = gpd.GeoDataFrame.from_dict(map2d_total_bounds_dict)

                    if map2d_total_bounds.shape[0] > 0:
                        map2d_total_bounds = map2d_total_bounds.set_crs(boundary.crs, allow_override=True)
                        map2d_total_bounds.to_file(map2d_total_bounds_out_filepath)



                    with open(features_out_filepath, "wb") as f:
                        f.write(
                            trimesh.exchange.gltf.export_glb(
                                scene_out, include_normals=True
                            )
                        )

