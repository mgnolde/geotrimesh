from pathlib import Path
import geopandas as gpd
from geotrimesh import GeoSceneSet
import tempfile
import os
import logging
import argparse

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument("out_dir", help="", nargs="?", default=Path(tempfile.gettempdir()))
args = parser.parse_args()
out_dirpath = args.out_dir

data_dirpath = Path(os.getcwd(), "demodata")
boundary_filepath = Path(data_dirpath, "bbox.gpkg")
buildings_filepaths = [Path(data_dirpath, "zurich_lod2_clip.glb")]

dem_filepaths = [Path(data_dirpath, "dtm_26830_12470_clip.tif")]
ortho_filepaths = [Path(data_dirpath, "2507_clip.tif")]
trees_filepaths = []

boundary = gpd.read_file(boundary_filepath).dissolve().explode(index_parts=True)
zurich = GeoSceneSet()

tilingscheme = GeoSceneSet.TilingScheme(boundary, dem_filepaths, height=256, width=256)
tilingscheme.gdf.to_file(Path(out_dirpath, "tiles.gpkg"))

#zurich.terrain = GeoSceneSet.Terrain(
#    out_dirpath=out_dirpath,
#    filepaths=dem_filepaths,
#    tiles=tilingscheme.tiles[0:3],
#    boundary=boundary
#)

zurich.buildings = GeoSceneSet.Features(
    "buildings",
    tilingscheme=tilingscheme,
    out_dirpath=out_dirpath,
    filepaths=buildings_filepaths,
    recombine_bodies=True,
    boundary=boundary,
    tiles=tilingscheme.tiles[0:3],
)

zurich.ortho = GeoSceneSet.Ortho(
    tilingscheme=tilingscheme,
    out_dirpath=out_dirpath,
    filepaths=ortho_filepaths,
    boundary=boundary,
    tiles=tilingscheme.tiles[0:3],
)
