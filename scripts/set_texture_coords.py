#!/usr/bin/env python

import sys
import argparse
import gdal, ogr
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--dem", action="store", type=str, required=True)
parser.add_argument("--clippoly", action="store", type=str, required=True)
parser.add_argument("--infile", action="store", type=str, required=True)
parser.add_argument("--outfile", action="store", type=str, required=True)

args = parser.parse_args()

dem_filepath = args.dem
clippoly_filepath = args.clippoly
in_filepath = args.infile
out_filepath = args.outfile


clippoly_driver = ogr.GetDriverByName("ESRI Shapefile")
clippoly_datasource = clippoly_driver.Open(clippoly_filepath, 0)
clippoly_layer = clippoly_datasource.GetLayer()
clippoly_layer_extent = clippoly_layer.GetExtent()

clippoly_layer_extent_xmin = clippoly_layer_extent[0]
clippoly_layer_extent_xmax = clippoly_layer_extent[1]
clippoly_layer_extent_ymin = clippoly_layer_extent[2]
clippoly_layer_extent_ymax = clippoly_layer_extent[3]

clippoly_layer_extent_xcent = (
    clippoly_layer_extent_xmax + clippoly_layer_extent_xmin
) / 2.0
clippoly_layer_extent_ycent = (
    clippoly_layer_extent_ymax + clippoly_layer_extent_ymin
) / 2.0

print(clippoly_layer_extent_xcent, clippoly_layer_extent_ycent)


dem_dataset = gdal.Open(dem_filepath)
dem_tmp_cols = dem_dataset.RasterXSize
dem_tmp_rows = dem_dataset.RasterYSize
dem_geotransform = dem_dataset.GetGeoTransform()
dem_xres = dem_geotransform[1]
dem_yres = abs(dem_geotransform[5])
dem_xmin = dem_geotransform[0]
dem_ymax = dem_geotransform[3]
dem_xmax = dem_xmin + dem_tmp_cols * dem_xres
dem_ymin = dem_ymax - dem_tmp_rows * dem_yres
dem_band = dem_dataset.GetRasterBand(1)
dem_tmp_array = dem_band.ReadAsArray(0, 0, dem_tmp_cols, dem_tmp_rows).astype(
    np.float32
)
dem_nodata = dem_band.GetNoDataValue()
dem_xcent = (dem_xmin + dem_xmax) / 2.0
dem_ycent = (dem_ymin + dem_ymax) / 2.0
dem_ydist = dem_ymax - dem_ymin
dem_xdist = dem_xmax - dem_xmin


print("rows, cols:", dem_tmp_rows, dem_tmp_cols)


x_ratio_min = None
x_ratio_max = None
y_ratio_min = None
y_ratio_max = None

x_min = None
x_max = None
y_min = None
y_max = None


for i in range(0, 2):

    if i == 1:
        suffix = ".mtl"
    else:
        suffix = ""

    with open(out_filepath + suffix, "w") as out_file:
        with open(in_filepath + suffix) as in_file:
            for in_line in in_file:

                if not in_line.startswith("vt "):

                    out_file.write(in_line.strip() + "\n")

                    if in_line.startswith("v "):
                        [ftype, x_local, z_local, y_local] = in_line.strip().split(" ")

                        x = float(x_local) + clippoly_layer_extent_xcent
                        y = float(y_local) + clippoly_layer_extent_ycent

                        # print(x, y, dem_xmin, dem_xmax, dem_ymin, dem_ymax)
                        # sys.exit()

                        x_offset = x - dem_xmin
                        y_offset = dem_ymax - y

                        # print(y_dist)

                        dem_x_dist = dem_xmax - dem_xmin
                        dem_y_dist = dem_ymax - dem_ymin

                        x_ratio = round(x_offset / float(dem_x_dist), 6)
                        y_ratio = round(y_offset / float(dem_y_dist), 6)

                        # if y_dist < 0:
                        #    print(y_expansion, y_dist, y_ratio)

                        """
                        if x_ratio < 0.0:
                            x_ratio = 0.0
                        if x_ratio >= 0.999999:
                            x_ratio = 0.999999
                        if y_ratio < 0.0:
                            y_ratio = 0.0
                        if y_ratio >= 0.999999:
                            y_ratio = 0.999999
                        """

                        if not x_min or x < x_min:
                            x_min = x
                        if not x_max or x > x_max:
                            x_max = x

                        if not y_min or y < y_min:
                            y_min = y
                        if not y_max or y > y_max:
                            y_max = y

                        if not x_ratio_min or x_ratio < x_ratio_min:
                            x_ratio_min = x_ratio
                        if not x_ratio_max or x_ratio > x_ratio_max:
                            x_ratio_max = x_ratio

                        if not y_ratio_min or y_ratio < y_ratio_min:
                            y_ratio_min = y_ratio
                        if not y_ratio_max or y_ratio > y_ratio_max:
                            y_ratio_max = y_ratio

                        out_file.write("vt " + str(x_ratio) + " " + str(y_ratio) + "\n")


print(round(dem_xmin, 1), round(dem_ymin, 1), round(dem_xmax, 1), round(dem_ymax, 1))
print(round(x_min, 1), round(y_min, 1), round(x_max, 1), round(y_max, 1))
print(
    round(x_ratio_min, 4),
    round(y_ratio_min, 4),
    round(x_ratio_max, 4),
    round(y_ratio_max, 4),
)