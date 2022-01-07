
DEM_FILENAME_BASE="gtopo30_1996_clip_pottenstein_25832"
#DEM_FILENAME_BASE="srtm3_2000_clip_pottenstein_25832"

mkdir -p ./out

python3 ./geotrimesh/generate_terrain.py \
  --dem ./demodata/${DEM_FILENAME_BASE}.tif \
  --clippoly ./demodata/osm_2017_pottenstein_simple_25832.shp \
  --outdir ./out \

openscad \
  ./out/${DEM_FILENAME_BASE}.scad \
  -o ./out/${DEM_FILENAME_BASE}.stl \
