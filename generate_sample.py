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

osm = os.path.join('demodata','osm201704_bound_adm8_pottenstein.shp')
srtm3 = os.path.join('demodata','srtm3_n049_e011_3arc_v2_clip.tif')
bluemarble = os.path.join('demodata','bmng_200407_21600x21600_c1_clip.tif')

elevation = mesh.ElevationMesh()
elevation.generate_mesh(dem=srtm3, orthophoto=bluemarble, boundaries=osm, mesh_prefix='pottenstein')
