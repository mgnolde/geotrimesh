DEM_FILEPATH="./../demodata/gtopo30_1996_clip_pottenstein_25832.tif"
ADM_FILEPATH="./../demodata/osm_2017_pottenstein_simple_25832.gpkg"

python3 ./../geotrimesh/export.py \
  --dem ${DEM_FILEPATH} \
  --clippoly ${ADM_FILEPATH} \
  --out ./out.scad

openscad \
  ./out.scad \
  -o ./out.stl
