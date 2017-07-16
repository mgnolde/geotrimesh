## *************************************************************
## Licensed to the Apache Software Foundation (ASF) under one
## or more contributor license agreements.  See the NOTICE file
## distributed with this work for additional information
## regarding copyright ownership.  The ASF licenses this file
## to you under the Apache License, Version 2.0 (the
## "License"); you may not use this file except in compliance
## with the License.  You may obtain a copy of the License at
##
##   http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing,
## software distributed under the License is distributed on an
## "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
## KIND, either express or implied.  See the License for the
## specific language governing permissions and limitations
## under the License.
## *************************************************************

## Author: Michael Nolde
## URL: http://www.flatpolar.org
## 2017/06

import os
from geotrimesh import mesh


natearth_full=os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'boundaries', 'naturalearth_global_201611', 'largescale', 'physical', 'orig', 'ne_10m_land.shp')
gtopo30_full=os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'dem', 'gtopo30_global_1996', 'gtopo30_merge_1996.tif')
bluemarble_full=os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'ortho', 'bluemarbleng_200407', 'bmng_merge_200407.tif')
mesh_path_full = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'tocomee', 'wgs84_full')
mesh_prefix_full = 'tocomee_wgs84_full'

#natearth_intermediate = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'boundaries', 'naturalearth_global_201611', 'largescale', 'physical', 'orig', 'ne_10m_land.shp')
#gtopo30_intermediate = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'dem', 'gtopo30_global_1996', 'gtopo30_merge_1996.tif')
#bluemarble_intermediate = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'ortho', 'bluemarbleng_200407', 'bmng_merge_200407.tif')
#mesh_path_intermediate = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'tocomee', 'wgs84_intermediate')
#mesh_prefix_intermediate = 'tocomee_wgs84_intermediate'

#natearth_crude = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'boundaries', 'naturalearth_global_201611', 'largescale', 'physical', 'orig', 'ne_10m_land.shp')
#gtopo30_crude = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'dem', 'gtopo30_global_1996', 'gtopo30_merge_1996.tif')
#bluemarble_crude = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'ortho', 'bluemarbleng_200407', 'bmng_merge_200407.tif')
#mesh_path_crude = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'tocomee', 'wgs84_crude')
#mesh_prefix_crude = 'tocomee_wgs84_crude'


elevation = mesh.ElevationMesh()

elevation.generate_mesh( dem=gtopo30_full, 
						 orthophoto=bluemarble_full, 
						 boundaries=natearth_full, 
						 dem_nodata=-9999,
						 orthophoto_nodata=-9999,
						 tiles_size=1,
						 mesh_path=mesh_path_full, 
						 mesh_prefix=mesh_prefix_full,
						 projection='orig', 
						 centering=True, 
						 indexed_colors=True, 
						 coloring_mode='orthophoto', 
						 mesh_format='x3d' )
