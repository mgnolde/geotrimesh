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


## South Dakota

bbox_sd = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'dem', 'srtm_southdakota_1arcsec_2014', 'n43_w103_1arc_v3_bbox_32613_innerbox_aoi1.shp')
tandem_sd = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'dem', 'srtm_southdakota_1arcsec_2014', 'n43_w103_1arc_v3_32613_innerbox.tif')
landsat8_sd = os.path.join(os.sep, 'media', 'sf_D_DRIVE', 'data', 'ortho', 'landsat8', 'southdakota_merged', 'southdakota_merged_clipped_composite_innerbox.tif')
                                                                                                                                                                                                                              
srtm1_sd_centroid_x = 701774.336819788
srtm1_sd_centroid_y = 4819497.48369646

#srtm1_sd_centroid_x = 728526
#srtm1_sd_centroid_y = 4811606

#tiles_bbox=[srtm1_sd_centroid_x - 1000, srtm1_sd_centroid_y - 1000, srtm1_sd_centroid_x + 1000, srtm1_sd_centroid_y + 1000],
#tiles_size=1000,



elevation = mesh.ElevationMesh()
elevation.generate_mesh( dem=tandem_sd,
                         orthophoto=landsat8_sd,
                         orthophoto_bitdepth=16,
                         boundaries=bbox_sd,
                         mesh_prefix='southdakota_srtm1_landsat8_32613',
                         mesh_format='x3d',
                         centering = True,
						 scale_xy=1.0,
                         projection='orig',
                         indexed_colors=False)

