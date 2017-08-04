# geoTriMesh


Converts GIS data to colored 3D meshes (X3D)

Usage:

	from geotrimesh import mesh

	elevation = mesh.ElevationMesh()
	elevation.generate_mesh(dem='srtm3.tif', orthophoto='bluemarble.tif', boundaries='osm.shp', mesh_prefix='sample')

![alt text](./demodata/sample_lq.png "Himalaya")
*Tile of the TOCOMEE dataset (visualized in Blender), created with geoTriMesh*

![alt text](./demodata/sample3_lq.png "Pottenstein mesh")
*Mesh of city of Pottenstein with 5x z-exaggeration  (visualized in Meshlab)*

![alt text](./demodata/sample4_lq.png "Pottenstein mesh detail")
*Detail mesh view*

![alt text](./demodata/sample5_lq.png "Globe")
*Complete TOCOMEE dataset (crude resolution), created with geoTriMesh



Options:

dem
| 	dem_nodata
| 	orthophoto
| 	orthophoto_nodata
| 	boundaries
| 	tiles_size
| 	tiles_bbox
| 	mesh_prefix
| 	mesh_path
| 	mesh_shapefile
| 	scale_xy
| 	z_exaggeration
| 	projection
| 	centering
| 	indexed_colors
| 	coloring_mode
| 	mesh_format
