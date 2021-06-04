#import meshio
import argparse
import vtk
import os
import json
import sys
from osgeo import ogr, osr, gdal
import math
import numpy as np


def get_mesh_elevation_from_xy(mesh, x_coord, y_coord):

    pSource = [x_coord, y_coord, -9999]
    pTarget = [x_coord, y_coord, 9999]

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


    #print(pointsIntersection)

    if not len(pointsIntersection) == 0:
        z_values = [pointsIntersection[0][2], pointsIntersection[1][2]]
    else:
        z_values = [None, None]

    return (z_values)


def coord_to_pix(x, y, rast_xmin, rast_ymax, rast_xres, rast_yres, cols, rows, offset_x, offset_y):
    #print((x + offset_x), (y + offset_y), rast_xmin, rast_ymax, rast_xres, rast_yres)

    col = math.ceil(((x + offset_x) - rast_xmin) / rast_xres)
    row = math.ceil((rast_ymax - (y + offset_y)) / rast_yres)

    if row < 0:
        row = 0
    if row > rows-1:
        row = rows-1

    if col < 0:
        col = 0
    if col > cols-1:
        col = cols-1


    return row, col


def generate_features(in_dtm_filepath, in_dom_filepath, in_scad_dtm_filepath, in_scad_dom_filepath, in_clippoly_filepath, in_mesh_dtm_filepath, in_mesh_dom_filepath, in_feat_filepath, asset_dirpath, zmean_total, task):


    print(in_dtm_filepath, in_dom_filepath, in_scad_dtm_filepath, in_scad_dom_filepath, in_clippoly_filepath, in_mesh_dtm_filepath, in_mesh_dom_filepath, in_feat_filepath)


    dtm_dataset = gdal.Open(in_dtm_filepath)
    dtm_cols = dtm_dataset.RasterXSize
    dtm_rows = dtm_dataset.RasterYSize
    dtm_geotransform = dtm_dataset.GetGeoTransform()
    dtm_xres = dtm_geotransform[1]
    dtm_yres = abs(dtm_geotransform[5])
    dtm_xmin = dtm_geotransform[0]
    dtm_ymax = dtm_geotransform[3]
    dtm_xmax = dtm_xmin + dtm_cols * dtm_xres
    dtm_ymin = dtm_ymax - dtm_rows * dtm_yres
    dtm_band = dtm_dataset.GetRasterBand(1)
    dtm_array = dtm_band.ReadAsArray(0, 0, dtm_cols, dtm_rows).astype(np.float32)
    dtm_nodata = dtm_band.GetNoDataValue()

    dom_dataset = gdal.Open(in_dom_filepath)
    dom_cols = dom_dataset.RasterXSize
    dom_rows = dom_dataset.RasterYSize
    dom_geotransform = dom_dataset.GetGeoTransform()
    dom_xres = dom_geotransform[1]
    dom_yres = abs(dom_geotransform[5])
    dom_xmin = dom_geotransform[0]
    dom_ymax = dom_geotransform[3]
    dom_xmax = dom_xmin + dom_cols * dom_xres
    dom_ymin = dom_ymax - dom_rows * dom_yres
    dom_band = dom_dataset.GetRasterBand(1)
    dom_array = dom_band.ReadAsArray(0, 0, dom_cols, dom_rows).astype(np.float32)
    dom_nodata = dom_band.GetNoDataValue()

    print("dtm:", dtm_rows, dtm_cols)


    in_mesh_dtm_dirpath = os.path.dirname(in_mesh_dtm_filepath)
    in_mesh_dtm_filename = os.path.basename(in_mesh_dtm_filepath)
    in_mesh_dtm_filename_base, in_mesh_dtm_filename_ext = os.path.splitext(in_mesh_dtm_filename)

    in_mesh_dom_dirpath = os.path.dirname(in_mesh_dom_filepath)
    in_mesh_dom_filename = os.path.basename(in_mesh_dom_filepath)
    in_mesh_dom_filename_base, in_mesh_dom_filename_ext = os.path.splitext(in_mesh_dom_filename)

    try:
        in_mesh_dtm_filename_base_prefix = '_'.join(in_mesh_dtm_filename_base.split('_')[0:-2])
        in_mesh_dom_filename_base_prefix = '_'.join(in_mesh_dom_filename_base.split('_')[0:-2])
        in_mesh_tiley = int(in_mesh_dtm_filename_base.split('_')[-2])
        in_mesh_tilex = int(in_mesh_dtm_filename_base.split('_')[-1])
    except:
        in_mesh_dtm_filename_base_prefix = None
        in_mesh_dom_filename_base_prefix = None
        in_mesh_tiley = None
        in_mesh_tilex = None


    in_mesh_neighbours_filenames_list = []
    in_mesh_neighbours_filenames_list.append([in_mesh_dtm_filename_base_prefix + '_' + str(in_mesh_tiley-1) + '_' + str(in_mesh_tilex-1), in_mesh_dom_filename_base_prefix + '_' + str(in_mesh_tiley-1) + '_' + str(in_mesh_tilex-1)])
    in_mesh_neighbours_filenames_list.append([in_mesh_dtm_filename_base_prefix + '_' + str(in_mesh_tiley-1) + '_' + str(in_mesh_tilex), in_mesh_dom_filename_base_prefix + '_' + str(in_mesh_tiley-1) + '_' + str(in_mesh_tilex)])
    in_mesh_neighbours_filenames_list.append([in_mesh_dtm_filename_base_prefix + '_' + str(in_mesh_tiley-1) + '_' + str(in_mesh_tilex+1), in_mesh_dom_filename_base_prefix + '_' + str(in_mesh_tiley-1) + '_' + str(in_mesh_tilex+1)])
    in_mesh_neighbours_filenames_list.append([in_mesh_dtm_filename_base_prefix + '_' + str(in_mesh_tiley) + '_' + str(in_mesh_tilex-1), in_mesh_dom_filename_base_prefix + '_' + str(in_mesh_tiley) + '_' + str(in_mesh_tilex-1)])
    in_mesh_neighbours_filenames_list.append([in_mesh_dtm_filename_base_prefix + '_' + str(in_mesh_tiley) + '_' + str(in_mesh_tilex+1), in_mesh_dom_filename_base_prefix + '_' + str(in_mesh_tiley) + '_' + str(in_mesh_tilex+1)])
    in_mesh_neighbours_filenames_list.append([in_mesh_dtm_filename_base_prefix + '_' + str(in_mesh_tiley+1) + '_' + str(in_mesh_tilex-1), in_mesh_dom_filename_base_prefix + '_' + str(in_mesh_tiley+1) + '_' + str(in_mesh_tilex-1)])
    in_mesh_neighbours_filenames_list.append([in_mesh_dtm_filename_base_prefix + '_' + str(in_mesh_tiley+1) + '_' + str(in_mesh_tilex), in_mesh_dom_filename_base_prefix + '_' + str(in_mesh_tiley+1) + '_' + str(in_mesh_tilex)])
    in_mesh_neighbours_filenames_list.append([in_mesh_dtm_filename_base_prefix + '_' + str(in_mesh_tiley+1) + '_' + str(in_mesh_tilex+1), in_mesh_dom_filename_base_prefix + '_' + str(in_mesh_tiley+1) + '_' + str(in_mesh_tilex+1)])




    c = []

    c2 = []
 

    if task.lower() == "extrude":
    
        current_block_is_dem = False
    
        with open(in_scad_dom_filepath, 'r') as scad_file:
            for scad_line in scad_file:
    
                if 'dem' in scad_line.lower():
                    current_block_is_dem = True 
    
                if '}' in scad_line.lower() and current_block_is_dem:
                    current_block_is_dem = False
    
                if 'polyhedron' in scad_line.lower() and current_block_is_dem:
                    c2.append('module dom() {')
                    c2.append('  ' + scad_line.strip())
                    c2.append('}')
     
    
        current_block_is_dem = False
    
        with open(in_scad_dtm_filepath, 'r') as scad_file:
            for scad_line in scad_file:
    
                if 'dem' in scad_line.lower():
                    current_block_is_dem = True 
    
                if '}' in scad_line.lower() and current_block_is_dem:
                    current_block_is_dem = False
    
                if 'polyhedron' in scad_line.lower() and current_block_is_dem:
                    c.append('module dem() {')
                    c.append('  ' + scad_line.strip())
                    c.append('}')
                    c.append('')
                    c.append('module feat() {')
                    c.append('  union() {')
                    
                    c2.append('module dem() {')
                    c2.append('  ' + scad_line.strip())
                    c2.append('}')
    
    
    
        
        #c2.append('rotate([-90,0,0])')
        c2.append('union() {')
    



    polygon_extrude_height = 10000.0
    polyhedron_extrude_height = 100.0

    source_srs = osr.SpatialReference()
    source_srs.ImportFromEPSG(4326)
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(2056)

    transform_srs = osr.CoordinateTransformation(source_srs, target_srs)
    transform_reverse_srs = osr.CoordinateTransformation(target_srs, source_srs)



    in_clippoly_driver = ogr.GetDriverByName("ESRI Shapefile")
    in_clippoly_datasource = in_clippoly_driver.Open(in_clippoly_filepath, 0)
    in_clippoly_layer = in_clippoly_datasource.GetLayer()

    in_clippoly_layer_extent = in_clippoly_layer.GetExtent()
    in_clippoly_layer_extent_xmin = in_clippoly_layer_extent[0]
    in_clippoly_layer_extent_xmax = in_clippoly_layer_extent[1]
    in_clippoly_layer_extent_ymin = in_clippoly_layer_extent[2]
    in_clippoly_layer_extent_ymax = in_clippoly_layer_extent[3]

    in_clippoly_layer_extent_xcent = (in_clippoly_layer_extent_xmax + in_clippoly_layer_extent_xmin) / 2.0
    in_clippoly_layer_extent_ycent = (in_clippoly_layer_extent_ymax + in_clippoly_layer_extent_ymin) / 2.0




    in_mesh_dtm_reader = vtk.vtkSTLReader()
    in_mesh_dtm_reader.SetFileName(in_mesh_dtm_filepath)
    in_mesh_dtm_reader.Update()
    in_mesh_dtm = in_mesh_dtm_reader.GetOutput()
    in_mesh_dtm_bounds = in_mesh_dtm.GetBounds()
    in_mesh_dtm_xmin, in_mesh_dtm_xmax, in_mesh_dtm_ymin, in_mesh_dtm_ymax, in_mesh_dtm_zmin, in_mesh_dtm_zmax = in_mesh_dtm_bounds

    in_mesh_dom_reader = vtk.vtkSTLReader()
    in_mesh_dom_reader.SetFileName(in_mesh_dom_filepath)
    in_mesh_dom_reader.Update()
    in_mesh_dom = in_mesh_dom_reader.GetOutput()
    in_mesh_dom_bounds = in_mesh_dom.GetBounds()
    in_mesh_dom_xmin, in_mesh_dom_xmax, in_mesh_dom_ymin, in_mesh_dom_ymax, in_mesh_dom_zmin, in_mesh_dom_zmax = in_mesh_dom_bounds



    in_mesh_dtm_ring = ogr.Geometry(ogr.wkbLinearRing)
    in_mesh_dtm_ring.AddPoint(in_mesh_dtm_xmin + in_clippoly_layer_extent_xcent, in_mesh_dtm_ymin + in_clippoly_layer_extent_ycent)
    in_mesh_dtm_ring.AddPoint(in_mesh_dtm_xmax + in_clippoly_layer_extent_xcent, in_mesh_dtm_ymin + in_clippoly_layer_extent_ycent)
    in_mesh_dtm_ring.AddPoint(in_mesh_dtm_xmax + in_clippoly_layer_extent_xcent, in_mesh_dtm_ymax + in_clippoly_layer_extent_ycent)
    in_mesh_dtm_ring.AddPoint(in_mesh_dtm_xmin + in_clippoly_layer_extent_xcent, in_mesh_dtm_ymax + in_clippoly_layer_extent_ycent)
    in_mesh_dtm_ring.AddPoint(in_mesh_dtm_xmin + in_clippoly_layer_extent_xcent, in_mesh_dtm_ymin + in_clippoly_layer_extent_ycent)
    in_mesh_dtm_bbox = ogr.Geometry(ogr.wkbPolygon)
    in_mesh_dtm_bbox.AddGeometry(in_mesh_dtm_ring)
    #in_mesh_dtm_bbox_latlon = in_mesh_dtm_bbox.Clone()
    #in_mesh_dtm_bbox_latlon.Transform(transform_reverse_srs)

    #print(in_mesh_dtm_xmin, in_mesh_dtm_xmax)
    print(in_mesh_dtm_ymin, in_mesh_dtm_ymax)


    #print(in_mesh_dtm_bounds)







    in_feat_driver = ogr.GetDriverByName("GPKG")
    in_feat_datasource = in_feat_driver.Open(in_feat_filepath, 0)
    in_feat_layer = in_feat_datasource.GetLayer()
    #in_feat_layer.SetSpatialFilter(in_mesh_dtm_bbox)
    in_feat_featcount = in_feat_layer.GetFeatureCount()


    print("File:", in_feat_filepath)
    print("Features:", in_feat_featcount)

    for feature_id, feature in enumerate(in_feat_layer):


        #if feature_id == 0:
        if True:
            #print(feature_id)

            geom_sub_dtm_minz_field = float(feature.GetField("DTM_MINZ"))
            geom_sub_dom_maxz_field = float(feature.GetField("DOM_MAXZ"))


            geom = feature.GetGeometryRef()

            #geom = in_mesh_dtm_bbox.Intersection(geom_orig)

            #in_buildings_geom_proj = in_buildings_geom.Clone()
            #in_buildings_geom_proj.Transform(transform_srs)

            #x_coord = in_buildings_geom.Centroid().GetX() - in_clippoly_layer_extent_xcent
            #y_coord = in_buildings_geom.Centroid().GetY() - in_clippoly_layer_extent_ycent
            #print(in_buildings_feature_id, x_coord, y_coord, get_mesh_elevation_from_xy(in_mesh_dtm, x_coord, y_coord))

            #if feature.GetField("name") != None:
            #    feature_name = feature.GetField("name")
            #else:
            #    feature_name = ""
            
            feature_name = ""       
            geom_name = str(geom.GetGeometryName())




            if geom_name.lower() == 'multipolygon':

                print('multipolygon', geom.GetGeometryCount())

                for sub_id in range(0, geom.GetGeometryCount()):


                    geom_sub = geom.GetGeometryRef(sub_id)

                    geom_sub_minx = None
                    geom_sub_maxx = None
                    geom_sub_miny = None
                    geom_sub_maxy = None
                    geom_sub_dtm_minz = None
                    geom_sub_dtm_maxz = None
                    geom_sub_dom_minz = None
                    geom_sub_dom_maxz = None
                    geom_sub_z_information_complete = True

                    polygon_points = []
                    polygon_paths = []
                    point_cnt = 0


                    centroid_x_tmp = geom_sub.Centroid().GetX()
                    centroid_y_tmp = geom_sub.Centroid().GetY()
                    centroid_x = centroid_x_tmp - in_clippoly_layer_extent_xcent
                    centroid_y = centroid_y_tmp - in_clippoly_layer_extent_ycent


                    for bound_id in range(0, geom_sub.GetGeometryCount()):

                        geom_sub_bound = geom_sub.GetGeometryRef(bound_id)
                        geom_sub_bound.Simplify(0.1)


                        geom_sub_bound_json = json.loads(geom_sub_bound.ExportToJson())
                        polygon_path = geom_sub_bound_json['coordinates']
                        polygon_points_geom = []

                        for point_id in range(0, geom_sub_bound.GetPointCount()):

                            geom_sub_bound_point_x_tmp, geom_sub_bound_point_y_tmp, dummy = geom_sub_bound.GetPoint(point_id)
                            #geom_sub_bound_point_x_tmp = geom_sub.Centroid().GetX()
                            #geom_sub_bound_point_y_tmp = geom_sub.Centroid().GetY()


                            geom_sub_bound_point_x = geom_sub_bound_point_x_tmp - in_clippoly_layer_extent_xcent
                            geom_sub_bound_point_y = geom_sub_bound_point_y_tmp - in_clippoly_layer_extent_ycent


                            #print(in_buildings_geom_sub_bound_point_x, in_buildings_geom_sub_bound_point_y, dummy)

                            #print(geom_sub_bound_point_z)

                            geom_sub_minx = geom_sub_bound_point_x if (not geom_sub_minx or geom_sub_bound_point_x < geom_sub_minx) else geom_sub_minx 
                            geom_sub_maxx = geom_sub_bound_point_x if (not geom_sub_maxx or geom_sub_bound_point_x > geom_sub_maxx) else geom_sub_maxx 
                            geom_sub_miny = geom_sub_bound_point_y if (not geom_sub_miny or geom_sub_bound_point_y < geom_sub_miny) else geom_sub_miny
                            geom_sub_maxy = geom_sub_bound_point_y if (not geom_sub_maxy or geom_sub_bound_point_y > geom_sub_maxy) else geom_sub_maxy 

                            #geom_sub_bound_dtm_point_z = get_mesh_elevation_from_xy(in_mesh_dtm, geom_sub_bound_point_x, geom_sub_bound_point_y)[1]
                            #geom_sub_bound_dom_point_z = get_mesh_elevation_from_xy(in_mesh_dom, geom_sub_bound_point_x, geom_sub_bound_point_y)[1]

                            row, col = coord_to_pix(centroid_x, centroid_y, dtm_xmin, dtm_ymax, dtm_xres, dtm_yres, dtm_cols, dtm_rows, in_clippoly_layer_extent_xcent, in_clippoly_layer_extent_ycent)
                            #print(row, col)
                            geom_sub_bound_dtm_point_z = dtm_array[row][col]
                            geom_sub_bound_dom_point_z = dom_array[row][col]



                            if geom_sub_bound_dtm_point_z:

                                geom_sub_dtm_minz = geom_sub_bound_dtm_point_z if (not geom_sub_dtm_minz or geom_sub_bound_dtm_point_z < geom_sub_dtm_minz) else geom_sub_dtm_minz
                                geom_sub_dtm_maxz = geom_sub_bound_dtm_point_z if (not geom_sub_dtm_maxz or geom_sub_bound_dtm_point_z > geom_sub_dtm_maxz) else geom_sub_dtm_maxz
                                geom_sub_dom_minz = geom_sub_bound_dom_point_z if (not geom_sub_dom_minz or geom_sub_bound_dom_point_z < geom_sub_dom_minz) else geom_sub_dom_minz
                                geom_sub_dom_maxz = geom_sub_bound_dom_point_z if (not geom_sub_dom_maxz or geom_sub_bound_dom_point_z > geom_sub_dom_maxz) else geom_sub_dom_maxz

                            else:

                                geom_sub_z_information_complete = False
                                geom_sub_dtm_minz = None
                                geom_sub_dtm_maxz = None
                                geom_sub_dom_minz = None
                                geom_sub_dom_maxz = None



                            polygon_path.append(point_cnt)
                            polygon_points_geom.append([geom_sub_bound_point_x, geom_sub_bound_point_y])
                            point_cnt +=1


                        if geom_sub_z_information_complete:
                            polygon_paths.append(polygon_path)
                            polygon_points += polygon_points_geom


                    if len(polygon_paths) > 0:
                        #print("extruded:", geom_sub_dom_maxz-geom_sub_dtm_minz)
                        feat_line = '    translate([0,0,{}]) linear_extrude(height={}) polygon({},{});'.format(geom_sub_dtm_minz_field, geom_sub_dom_maxz_field-geom_sub_dtm_minz_field, polygon_points, polygon_paths)
                        feat_line2 = '    intersection() { ' + 'translate([0,0,{}]) linear_extrude(height={}) polygon({},{});'.format(-9999, 20000, polygon_points, polygon_paths) + ' dem(); }'
                        feat_line3 = '    intersection() { ' + 'translate([0,0,{}]) linear_extrude(height={}) polygon({},{});'.format(-9999, 20000, polygon_points, polygon_paths) + ' dom(); }'

                        c.append(feat_line)
                        c2.append(feat_line)   ## output extruded polygon
                        #c2.append(feat_line2) ## output section cropped from original mesh
                        #c2.append(feat_line3) ## output section cropped from DOM mesh
 


                #print('[', feature_id, ']', '[', sub_id, ']', geom_sub_dtm_minz_field, geom_sub_dom_maxz_field, geom_sub_z_information_complete)




            elif geom_name.lower() == 'polygon':
            
                print('polygon', geom.GetGeometryCount())

                geom_dtm_minz_field = feature.GetField("DTM_MINZ") 
                geom_dom_maxz_field = feature.GetField("DOM_MAXZ") 


                geom_minx = None
                geom_maxx = None
                geom_miny = None
                geom_maxy = None
                geom_dtm_minz = None
                geom_dtm_maxz = None
                geom_dom_minz = None
                geom_dom_maxz = None
                geom_z_information_complete = True

                polygon_points = []
                polygon_paths = []
                point_cnt = 0

                centroid_x_tmp = geom.Centroid().GetX()
                centroid_y_tmp = geom.Centroid().GetY()
                centroid_x = centroid_x_tmp - in_clippoly_layer_extent_xcent
                centroid_y = centroid_y_tmp - in_clippoly_layer_extent_ycent



                for bound_id in range(0, geom.GetGeometryCount()):

                    geom_bound = geom.GetGeometryRef(bound_id)
                    geom_bound.Simplify(0.1)

                    geom_bound_json = json.loads(geom_bound.ExportToJson())
                    polygon_path = geom_bound_json['coordinates']
                    polygon_points_geom = []

                    for point_id in range(0, geom_bound.GetPointCount()):

                        geom_bound_point_x_tmp, geom_bound_point_y_tmp, dummy = geom_bound.GetPoint(point_id)
                        #print("Coord", geom_bound_point_x_tmp, geom_bound_point_y_tmp)


 

                        #print("Centr", geom_bound_point_x_tmp, geom_bound_point_y_tmp)

                        geom_bound_point_x = geom_bound_point_x_tmp - in_clippoly_layer_extent_xcent
                        geom_bound_point_y = geom_bound_point_y_tmp - in_clippoly_layer_extent_ycent


                        #print(in_buildings_geom_bound_point_x, in_buildings_geom_bound_point_y, dummy)

                        geom_minx = geom_bound_point_x if (not geom_minx or geom_bound_point_x < geom_minx) else geom_minx 
                        geom_maxx = geom_bound_point_x if (not geom_maxx or geom_bound_point_x > geom_maxx) else geom_maxx 
                        geom_miny = geom_bound_point_y if (not geom_miny or geom_bound_point_y < geom_miny) else geom_miny 
                        geom_maxy = geom_bound_point_y if (not geom_maxy or geom_bound_point_y > geom_maxy) else geom_maxy


                        #geom_bound_dtm_point_z = get_mesh_elevation_from_xy(in_mesh_dtm, geom_bound_point_x, geom_bound_point_y)[1]
                        #geom_bound_dom_point_z = get_mesh_elevation_from_xy(in_mesh_dom, geom_bound_point_x, geom_bound_point_y)[1]

    

                        row, col = coord_to_pix(centroid_x, centroid_y, dtm_xmin, dtm_ymax, dtm_xres, dtm_yres, dtm_cols, dtm_rows, in_clippoly_layer_extent_xcent, in_clippoly_layer_extent_ycent)
                        #print(row, col)
                        geom_bound_dtm_point_z = dtm_array[row][col]
                        geom_bound_dom_point_z = dom_array[row][col]



                        #print("point_z", geom_bound_point_x, geom_bound_point_y)
                        #print("val_z", geom_bound_dtm_point_z, geom_bound_dom_point_z) 
                        #print()

                        if geom_bound_dtm_point_z:

                            geom_dtm_minz = geom_bound_dtm_point_z if (not geom_dtm_minz or geom_bound_dtm_point_z < geom_dtm_minz) else geom_dtm_minz 
                            geom_dtm_maxz = geom_bound_dtm_point_z if (not geom_dtm_maxz or geom_bound_dtm_point_z > geom_dtm_maxz) else geom_dtm_maxz
                            geom_dom_minz = geom_bound_dom_point_z if (not geom_dom_minz or geom_bound_dom_point_z < geom_dom_minz) else geom_dom_minz 
                            geom_dom_maxz = geom_bound_dom_point_z if (not geom_dom_maxz or geom_bound_dom_point_z > geom_dom_maxz) else geom_dom_maxz

                        else:

                            geom_z_information_complete = False
                            geom_dtm_minz = None
                            geom_dtm_maxz = None
                            geom_dom_minz = None
                            geom_dom_maxz = None



                        polygon_path.append(point_cnt)
                        polygon_points_geom.append([geom_bound_point_x, geom_bound_point_y])
                        point_cnt +=1


                    if geom_z_information_complete:
                        polygon_paths.append(polygon_path)
                        polygon_points += polygon_points_geom



                if task.lower() == "extrude":

                    if len(polygon_paths) > 0:
                        feat_line = '    translate([0,0,{}]) linear_extrude(height={}) polygon({},{});'.format(geom_dtm_minz_field, geom_dom_maxz_field-geom_dtm_minz_field, polygon_points, polygon_paths)
                        feat_line2 = '    intersection() { ' + 'translate([0,0,{}]) linear_extrude(height={}) polygon({},{});'.format(-9999, 20000, polygon_points, polygon_paths) + ' dem(); }'
                        feat_line3 = '    intersection() { ' + 'translate([0,0,{}]) linear_extrude(height={}) polygon({},{});'.format(-9999, 20000, polygon_points, polygon_paths) + ' dom(); }'
                        c.append(feat_line)
                        c2.append(feat_line)   ## output extruded polygon
                        #c2.append(feat_line2) ## output section cropped from original mesh
                        #c2.append(feat_line3) ## output section cropped from DOM mesh


                if task.lower() == "insert":

                    asset_filepath = os.path.join(asset_dirpath, "quaking_aspen.stl")

                    if len(polygon_paths) > 0:
                        feat_line = '     translate([{},{},{}]) rotate([90,0,0]) import("{}", convexity=10);'.format(centroid_x, centroid_y, geom_dtm_minz_field, asset_filepath)
                        #c.append(feat_line)
                        c2.append(feat_line)   ## output extruded polygon


                    #print('[', feature_id, ']', geom_dom_maxz-geom_dtm_minz, geom_dtm_minz, geom_dom_maxz, feature_name)


    if task.lower() == "extrude":


        c.append('  }')
        c.append('}')
        c.append('')
        #c.append('rotate([-90,0,0])')
        c.append('difference() {')
        c.append('  dem();')
        c.append('  feat();')
        c.append('}')
    
    
        c2.append('}')
 

    return c, c2





in_clippoly_filepath = os.path.join(os.sep, 'mnt', 'e', 'zh', 'wambachers_osm__boundaries__adm', 'data', 'district_zurich_al6_al6_2056.shp')
in_mesh_dtm_filepath = os.path.join(os.sep, 'mnt', 'e', 'zh', 'gis_zh__dom_dtm__lidar', 'dtm_mosaic_tiled_offset_500x500_scad', 'dtm_mosaic_tiled_offset_500x500_11_18.stl')
in_mesh_dom_filepath = os.path.join(os.sep, 'mnt', 'e', 'zh', 'gis_zh__dom_dtm__lidar', 'dom_mosaic_tiled_offset_500x500_scad', 'dom_mosaic_tiled_offset_500x500_11_18.stl')



#"""

scad_filename_base = 'test'
proc_target_os = 'win'

openscad_bin_filepath = 'openscad'
scad_dirpath = os.path.join(os.sep, 'mnt', 'c', 'Users', 'mic', 'dev')
proc_dirpath = os.path.join(os.sep, 'mnt', 'c', 'Users', 'mic', 'dev')




#dem_filepath = os.path.join(os.sep, 'mnt', 'e', 'zh', 'gis_zh__dom_dtm__lidar', 'dtm_mosaic_tiled_offset_500x500', 'dtm_mosaic_tiled_offset_500x500_34_33.tif')
#clippoly_filepath = os.path.join(os.sep, 'mnt', 'e', 'zh', 'wambachers_osm__boundaries__adm', 'data', 'district_zurich_al6_al6_2056.shp')



parser = argparse.ArgumentParser()
parser.add_argument('--dtm', action='store', type=str, required=True)
parser.add_argument('--dom', action='store', type=str, required=True)
parser.add_argument('--mesh_dtm', action='store', type=str, required=True)
parser.add_argument('--mesh_dom', action='store', type=str, required=True)
parser.add_argument('--scad_dtm', action='store', type=str, required=True)
parser.add_argument('--scad_dom', action='store', type=str, required=True)
parser.add_argument('--zmin', action='store', type=str, required=False)
parser.add_argument('--zmax', action='store', type=str, required=False)
parser.add_argument('--feat', action='store', type=str, required=True)
parser.add_argument('--task', action='store', type=str, required=True)
parser.add_argument('--suffix', action='store', type=str, required=True)
parser.add_argument('--clippoly', action='store', type=str, required=True)
parser.add_argument('--assetdir', action='store', type=str, required=True)
parser.add_argument('--outdir', action='store', type=str, required=True)
#parser.add_argument('--id', action='store', type=int)

args = parser.parse_args()

in_dtm_filepath = args.dtm
in_dom_filepath = args.dom
in_scad_dtm_filepath = args.scad_dtm
in_scad_dom_filepath = args.scad_dom
in_mesh_dtm_filepath = args.mesh_dtm
in_mesh_dom_filepath = args.mesh_dom
in_clippoly_filepath = args.clippoly
filename_suffix = args.suffix
task = args.task
proc_dirpath = args.outdir
asset_dirpath = args.assetdir


try:
    zmin_total = float(args.zmin)
    zmax_total = float(args.zmax)
    zmean_total = zmin_total + ((zmax_total - zmin_total) / 2.0)
except:
    zmean_total = 0
    
    
#in_buildings_filepath = os.path.join(os.sep, 'mnt', 'e', 'zh', 'geofabrik_osm__features__adm', 'osm', 'gis_osm_buildings_a_free_1_clip_2056.gpkg')
in_feat_filepath = args.feat


#sys.exit()
in_mesh_dtm_dirpath, in_mesh_dtm_filename = os.path.split(in_mesh_dtm_filepath)
in_mesh_dtm_filename_base = in_mesh_dtm_filename.split('.')[0]
scad_dirpath = args.outdir

scad_dtm_filename_base = in_mesh_dtm_filename_base + '_mold'
scad_feat_filename_base = in_mesh_dtm_filename_base + '_feat_' + filename_suffix


mold_command_lines, feat_command_lines = generate_features(in_dtm_filepath, in_dom_filepath, in_scad_dtm_filepath, in_scad_dom_filepath, in_clippoly_filepath, in_mesh_dtm_filepath, in_mesh_dom_filepath, in_feat_filepath, asset_dirpath, zmean_total, task)


if task == "extrude":
    if len(mold_command_lines) > 0:
        with open(os.path.join(scad_dirpath, scad_dtm_filename_base + '.scad'), 'w') as scad_file:
            for mold_command_line in mold_command_lines:
                scad_file.write(mold_command_line + '\n')


if task == "extrude":
    with open(os.path.join(scad_dirpath, scad_feat_filename_base + '.scad'), 'w') as scad_file:
        for feat_command_line in feat_command_lines:
            scad_file.write(feat_command_line + '\n')
if task == "insert":
    for feat_command_line_id, feat_command_line in enumerate(feat_command_lines):
        with open(os.path.join(scad_dirpath, scad_feat_filename_base + str(feat_command_line_id) + '.scad'), 'w') as scad_file:
            scad_file.write(feat_command_line + '\n')






if not proc_target_os == 'win':

    subprocess_bin = openscad_bin_filepath
    subprocess_commands = [subprocess_bin, '-o', os.path.join(proc_dirpath, scad_filename_base + '.stl'), os.path.join(scad_dirpath, scad_filename_base + '.scad')]
    output = subprocess.check_output(subprocess_commands, shell=False)

else:

    openscad_win_bin_filepath = '\\'.join(['C:', '"Program Files"', 'OpenSCAD', 'openscad.com'])
    openscad_linux_bin_filepath = os.path.join(os.sep, 'home', 'mic', 'prog', 'OpenSCAD-2019.05-x86_64', 'squashfs-root', 'usr', 'bin', 'openscad')

    scad_win_dirpath = '\\'.join(['C:', 'Users', 'mic', 'dev'])
    proc_win_dirpath = '\\'.join(['C:', 'Users', 'mic', 'dev'])

    subprocess_dtm_args = ['-o', scad_dtm_filename_base + '.stl', scad_dtm_filename_base + '.scad']
    subprocess_feat_args = ['-o', scad_feat_filename_base + '.stl', scad_feat_filename_base + '.scad']

    subprocess_win_dtm_command = [openscad_win_bin_filepath] + subprocess_dtm_args
    subprocess_win_feat_command = [openscad_win_bin_filepath] + subprocess_feat_args
    subprocess_linux_dtm_command = [openscad_linux_bin_filepath] + subprocess_dtm_args
    subprocess_linux_feat_command = [openscad_linux_bin_filepath] + subprocess_feat_args

    with open(os.path.join(scad_dirpath, scad_dtm_filename_base + '.bat'), 'w') as batch_file:
        batch_file.write(' '.join(subprocess_win_dtm_command) + '\n')
        batch_file.write(' '.join(subprocess_win_feat_command) + '\n')

    with open(os.path.join(scad_dirpath, scad_dtm_filename_base + '.sh'), 'w') as shell_file:
        shell_file.write(' '.join(subprocess_linux_dtm_command) + '\n')
        shell_file.write(' '.join(subprocess_linux_feat_command) + '\n')


#"""
