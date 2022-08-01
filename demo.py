from pathlib import Path
import copy
import trimesh
import fiona
import geopandas as gpd
from geotrimesh import GeoSceneSet
import sys
import os
import shutil
import logging
import tempfile

logging.basicConfig(level=logging.INFO)

data_dirpath = Path(os.getcwd(), "demodata")
boundary_filepath = Path(data_dirpath, "bbox.gpkg")
buildings_filepaths = [Path(data_dirpath, "data.gltf")]
dem_filepaths = [Path(data_dirpath, "dtm_26830_12470_clip.tif")]
ortho_filepaths = [Path(data_dirpath, "2507_clip.tif")]
trees_filepaths = []

boundary = gpd.read_file(boundary_filepath).dissolve().explode(index_parts=True)
out_dirpath = Path(os.getcwd())

zurich = GeoSceneSet()

tilingscheme = GeoSceneSet.TilingScheme(boundary, dem_filepaths, height=256, width=256)
tilingscheme.gdf.to_file(Path(out_dirpath, "tiles.gpkg"))

zurich.terrain = GeoSceneSet.Terrain(
    out_dirpath=out_dirpath,
    filepaths=dem_filepaths, 
    tiles=tilingscheme.tiles[2:3],
    boundary=boundary
    )

zurich.buildings = GeoSceneSet.Features("buildings",
    tilingscheme=tilingscheme,
    out_dirpath=out_dirpath,
    filepaths=buildings_filepaths,
    recombine_bodies=True,
    boundary=boundary, 
    tiles=tilingscheme.tiles[2:3],
    )

zurich.ortho = GeoSceneSet.Ortho(
    tilingscheme=tilingscheme,
    out_dirpath=out_dirpath,
    filepaths=ortho_filepaths,
    boundary=boundary, 
    tiles=tilingscheme.tiles[2:3],
    )
