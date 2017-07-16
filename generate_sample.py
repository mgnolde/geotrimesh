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

osm = os.path.join('demodata','osm_2017_pottenstein_simple.shp')
srtm3 = os.path.join('demodata','srtm3_2000_clip_pottenstein.tif')
gtopo30 = os.path.join('demodata','gtopo30_1996_clip_pottenstein.tif')
bluemarble = os.path.join('demodata','bmng_2004_clip_pottenstein.tif')

elevation = mesh.ElevationMesh()

elevation.generate_mesh( dem=gtopo30, orthophoto=bluemarble, boundaries=osm, 
					     mesh_prefix='pottenstein', z_exaggeration=5)
