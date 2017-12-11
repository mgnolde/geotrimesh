from osgeo import gdal, ogr, osr
import sys
import os
import math
import numpy as np
import logging
from scipy.spatial import Delaunay
import argparse
import json
import time


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



class ElevationMesh(object):
    
    def __init__(self):
        gdal.UseExceptions()
        pass


    def generate_mesh(self, dem=None, orthophoto=None, boundaries=None, dem_nodata=None, orthophoto_nodata=None, tiles_size=None, tiles_bbox=None, mesh_prefix='out', 
                        mesh_path=os.path.join(os.getcwd(),'out'), mesh_shapefile=False, scale_xy=None, z_exaggeration=1.0, projection='orig', centering=True, indexed_colors=False, 
                        coloring_mode='orthophoto', mesh_format='x3d', orthophoto_bitdepth=8, tiles_naming_convention='xy', force_out_mesh_overwrite=True):

        in_dem_filename = dem
        in_orthophoto_filename = orthophoto
        in_boundaries_filename = boundaries
        in_dem_nodata_ext=dem_nodata
        in_orthophoto_nodata_ext=orthophoto_nodata
        out_mesh_filename_prefix=mesh_prefix
        out_mesh_path=mesh_path



        in_dem = gdal.Open(in_dem_filename)
        in_dem_res_x = float(in_dem.GetGeoTransform()[1])
        in_dem_res_y = float(abs(in_dem.GetGeoTransform()[5]))
        in_dem_cols = in_dem.RasterXSize
        in_dem_rows = in_dem.RasterYSize
        in_dem_prj=in_dem.GetProjection()
        in_dem_srs=osr.SpatialReference(wkt=in_dem_prj)

        in_orthophoto = gdal.Open(in_orthophoto_filename)
        in_boundaries = gdal.Open(in_boundaries_filename)

        in_dem_extent_x_min = float(in_dem.GetGeoTransform()[0])
        in_dem_extent_y_max = float(in_dem.GetGeoTransform()[3])
        in_dem_extent_x_max = float(in_dem_extent_x_min + (in_dem_cols * in_dem_res_x))
        in_dem_extent_y_min = float(in_dem_extent_y_max - (in_dem_rows * in_dem_res_y))

        in_dem_band = in_dem.GetRasterBand(1)
        in_boundaries_band = in_boundaries.GetRasterBand(1)
        in_orthophoto_red_band = in_orthophoto.GetRasterBand(1)
        in_orthophoto_green_band = in_orthophoto.GetRasterBand(2)
        in_orthophoto_blue_band = in_orthophoto.GetRasterBand(3)

        z_grid = in_dem_band.ReadAsArray(0, 0, in_dem_cols, in_dem_rows).astype(np.int16)

        red_grid = in_orthophoto_red_band.ReadAsArray(0, 0, in_dem_cols, in_dem_rows).astype(np.uint8)
        green_grid = in_orthophoto_green_band.ReadAsArray(0, 0, in_dem_cols, in_dem_rows).astype(np.uint8)
        blue_grid = in_orthophoto_blue_band.ReadAsArray(0, 0, in_dem_cols, in_dem_rows).astype(np.uint8)
        
        boundaries_grid = in_boundaries_band.ReadAsArray(0, 0, in_dem_cols, in_dem_rows).astype(np.uint8)

        #red_grid_min, red_grid_max = (307, 1769)
        #green_grid_min, green_grid_max = (540,1559)
        #blue_grid_min, blue_grid_max = (711,1557)

        #red_grid_norm = 2 * ((red_grid-red_grid_min) * 100) / (red_grid_max-red_grid_min)
        #green_grid_norm = 2 * ((green_grid-green_grid_min) * 100) / (green_grid_max-green_grid_min)
        #blue_grid_norm = 2 * ((blue_grid-blue_grid_min) * 100) / (blue_grid_max-blue_grid_min)

       
        x = np.linspace(in_dem_extent_x_min + (in_dem_res_x/2.0), in_dem_extent_x_max - (in_dem_res_x/2.0), in_dem_cols, endpoint=True)
        y = np.linspace(in_dem_extent_y_min + (in_dem_res_y/2.0), in_dem_extent_y_max - (in_dem_res_y/2.0), in_dem_rows, endpoint=True)

        #x = np.arange(in_dem_extent_x_min + (in_dem_res_x/2.0), in_dem_extent_x_max - (in_dem_res_x/2.0), in_dem_res_x)
        #y = np.arange(in_dem_extent_y_min + (in_dem_res_y/2.0), in_dem_extent_y_max - + (in_dem_res_y/2.0), in_dem_res_y)

        print(in_dem_extent_x_max - (in_dem_res_x/2.0))
        print(x[-1])

        x_grid, y_grid = np.meshgrid(x,y)

        print(in_dem_rows, in_dem_cols)

        triangles = np.zeros(shape=( ((in_dem_rows*in_dem_cols)*2), 3), dtype=np.int64)


        triangle_id = 0
        #for in_dem_y in range(0, 5):
        #    for in_dem_x in range(0, 5):


        for in_dem_y in range(0, in_dem_rows-2):
            for in_dem_x in range(0, in_dem_cols-2):

                #if not(boundaries_grid[in_dem_y][in_dem_x] == boundaries_grid[in_dem_y+1][in_dem_x] == boundaries_grid[in_dem_y][in_dem_x+1] == 1):
                if (boundaries_grid[in_dem_y][in_dem_x] == boundaries_grid[in_dem_y+1][in_dem_x] == boundaries_grid[in_dem_y][in_dem_x+1] == 1):
                    triangles[triangle_id][0] = ((in_dem_y+1) * in_dem_cols) + in_dem_x
                    triangles[triangle_id][1] = (in_dem_y * in_dem_cols) + in_dem_x+1
                    triangles[triangle_id][2] = (in_dem_y * in_dem_cols) + in_dem_x
                    triangle_id += 1

                #if not(boundaries_grid[in_dem_y+1][in_dem_x] == boundaries_grid[in_dem_y+1][in_dem_x+1] == boundaries_grid[in_dem_y][in_dem_x+1] == 1):
                if (boundaries_grid[in_dem_y+1][in_dem_x] == boundaries_grid[in_dem_y+1][in_dem_x+1] == boundaries_grid[in_dem_y][in_dem_x+1] == 1):
                    triangles[triangle_id][0] = ((in_dem_y+1) * in_dem_cols) + in_dem_x
                    triangles[triangle_id][1] = ((in_dem_y+1) * in_dem_cols) + in_dem_x+1
                    triangles[triangle_id][2] = (in_dem_y * in_dem_cols) + in_dem_x+1
                    triangle_id += 1

        triangles_count = triangle_id



        print(x_grid.flatten().shape)
        print(y_grid.flatten().shape)
        #print(z_grid.flatten().shape)

        self.write_vtu(x_grid, y_grid, z_grid, red_grid, green_grid, blue_grid, triangles)




    def write_vtu(self, x_grid, y_grid, z_grid, red_grid, green_grid, blue_grid, triangles):


        with open('out_mesh.vtu', 'w') as out_mesh:     
    
            out_mesh.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">' + '\n')
            out_mesh.write('  <UnstructuredGrid>' + '\n')
            out_mesh.write('    <Piece NumberOfPoints="' + str(len(x_grid.flatten())) + '" NumberOfCells="' + str(len(triangles[:,0])) + '">' + '\n')


            out_mesh.write('      <PointData Scalars="Colors">' + '\n')
            out_mesh.write('        <DataArray type="UInt8" Name="Colors" format="ascii" RangeMin="0" RangeMax="255" NumberOfComponents="3">' + '\n')
            out_mesh.write('          ')

            for color_id, (color_red, color_green, color_blue) in enumerate(zip(red_grid.flatten(), green_grid.flatten(), blue_grid.flatten())):
                out_mesh.write(str(color_red) + ' ' + str(color_green) + ' ' + str(color_blue) + ' ')

            out_mesh.write('' + '\n')
            out_mesh.write('        </DataArray>' + '\n')
            out_mesh.write('      </PointData>' + '\n')

            out_mesh.write('      <Points>' + '\n')
            out_mesh.write('        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="' + str(min(np.nanmin(x_grid), np.nanmin(y_grid))) + '" RangeMax="' + str(min(np.nanmax(x_grid), np.nanmax(y_grid))) + '">' + '\n')
    
    
            out_mesh.write('          ')

            for coord_id, (coord_x, coord_y, coord_z) in enumerate(zip(x_grid.flatten(), y_grid.flatten(), z_grid.flatten())):
                out_mesh.write(str(coord_x) + ' ' + str(coord_y) + ' ' + str(coord_z) + ' ')
    
    
            out_mesh.write('' + '\n')
    
            out_mesh.write('        </DataArray>' + '\n')
            out_mesh.write('      </Points>' + '\n')
            out_mesh.write('      <Cells>' + '\n')
            out_mesh.write('        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="' +  str(0) + '" RangeMax="' +  str(len(triangles[:,0])) + '">' + '\n')
    
            out_mesh.write('          ')
            for triangle_id, triangle in enumerate(triangles):
                point_a, point_b, point_c = triangle
                out_mesh.write(str(point_a) + ' ' + str(point_b) + ' ' + str(point_c) + ' ')
    
   
            out_mesh.write('' + '\n')
    
            out_mesh.write('        </DataArray>' + '\n')
    
    
            out_mesh.write('        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="' + '3' + '" RangeMax="' + str(len(triangles[:,0])*3) + '">' + '\n')
    
            out_mesh.write('          ')
            for triangle_id, triangle in enumerate(triangles):
                out_mesh.write(str((triangle_id+1)*3) + ' ')
    
   
            out_mesh.write('' + '\n')
    
            out_mesh.write('        </DataArray>' + '\n')
            out_mesh.write('        <DataArray type="UInt8" Name="types" format="ascii" RangeMin="5" RangeMax="5">' + '\n')
            
    
            out_mesh.write('          ')
            for triangle_id, triangle in enumerate(triangles):
                out_mesh.write('5' + ' ')
    
    
            out_mesh.write('' + '\n')
                   
            out_mesh.write('        </DataArray>' + '\n')
            out_mesh.write('      </Cells>' + '\n')
            out_mesh.write('    </Piece>' + '\n')
            out_mesh.write('  </UnstructuredGrid>' + '\n')
            out_mesh.write('</VTKFile>' + '\n')
    
