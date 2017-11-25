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

osm = os.path.join('demodata','osm_2017_pottenstein_simple_25832.shp')
srtm3 = os.path.join('demodata','srtm3_2000_clip_pottenstein_25832.tif')
gtopo30 = os.path.join('demodata','gtopo30_1996_clip_pottenstein_25832.tif')
bluemarble = os.path.join('demodata','bmng_2004_clip_pottenstein_25832.tif')

elevation = mesh.ElevationMesh()

## Py/Matplotlib
elevation.generate_mesh( dem=gtopo30, orthophoto=bluemarble, boundaries=osm, 
					     mesh_prefix='pottenstein_gtopo30', mesh_format='py', centering = False, projection='orig', scale_xy = 1.0, indexed_colors=False)
#elevation.generate_mesh( dem=srtm3, orthophoto=bluemarble, boundaries=osm, 
#					     mesh_prefix='pottenstein_srtm3', mesh_format='py', centering = False, projection='orig', scale_xy = 1.0, indexed_colors=False)

## VTK
#elevation.generate_mesh( dem=gtopo30, orthophoto=bluemarble, boundaries=osm, 
#					     mesh_prefix='pottenstein_gtopo30', mesh_format='vtu', centering=False, projection='orig', scale_xy=1.0, indexed_colors=False)
#elevation.generate_mesh( dem=srtm3, orthophoto=bluemarble, boundaries=osm, 
#					     mesh_prefix='pottenstein_srtm3', mesh_format='vtu', centering=False, projection='orig', scale_xy=1.0, indexed_colors=False)

## X3D
#elevation.generate_mesh( dem=gtopo30, orthophoto=bluemarble, boundaries=osm, 
#					     mesh_prefix='pottenstein_gtopo30', mesh_format='x3d', centering = False, projection='orig', scale_xy=1.0, z_exaggeration=1.0, indexed_colors=False)
#elevation.generate_mesh( dem=srtm3, orthophoto=bluemarble, boundaries=osm, 
#					     mesh_prefix='pottenstein_srtm3', mesh_format='x3d', centering = False, projection='orig', scale_xy=1.0, z_exaggeration=1.0, indexed_colors=False)
