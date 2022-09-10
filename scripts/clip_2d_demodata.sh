data_dirpath="./data"
demodata_dirpath="./demodata"

margin=200
x_min=2683090.75
x_max=2683311.75
y_min=1247705.25
y_max=1247968.875


x_min_fin=$(echo "$x_min - $margin"|bc)
y_min_fin=$(echo "$y_min - $margin"|bc)
x_max_fin=$(echo "$x_max + $margin"|bc)
y_max_fin=$(echo "$y_max + $margin"|bc)

rm -f ${demodata_dirpath}/2507_clip.tif
gdal_translate \
  -projwin ${x_min_fin} ${y_max_fin} ${x_max_fin} ${y_min_fin} \
  -of GTiff \
  -ot Byte \
  -a_srs EPSG:2056 \
  -co COMPRESS=JPEG \
  -co PHOTOMETRIC=YCBCR \
  -b 1 \
  -b 2 \
  -b 3 \
  ${data_dirpath}/2507.tif \
  ${demodata_dirpath}/2507_clip.tif


rm -f ${demodata_dirpath}/dtm_26830_12470_clip.tif
gdal_translate \
  -projwin ${x_min_fin} ${y_max_fin} ${x_max_fin} ${y_min_fin} \
  -of GTiff \
  -ot Float32 \
  -a_srs EPSG:2056 \
  -co COMPRESS=LZW \
  ${data_dirpath}/dtm_26830_12470.tif \
  ${demodata_dirpath}/dtm_26830_12470_clip.tif

rm -f ${demodata_dirpath}/dtm_26830_12470_clip_lq.tif
gdal_translate \
  -projwin ${x_min_fin} ${y_max_fin} ${x_max_fin} ${y_min_fin} \
  -of GTiff \
  -ot Float32 \
  -a_srs EPSG:2056 \
  -co COMPRESS=LZW \
  -tr 5.0 5.0 \
  ${data_dirpath}/dtm_26830_12470.tif \
  ${demodata_dirpath}/dtm_26830_12470_clip_lq.tif

