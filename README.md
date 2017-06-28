# geoTriMesh

Converts GIS data to colored 3D meshes (X3D)

Usage:

	from geotrimesh import mesh

	elevation = mesh.ElevationMesh()
	elevation.generate_mesh(dem='srtm3.tif', ortho='bluemarble.tif', bound='osm.shp', mesh_prefix='sample')


![alt text](./demodata/sample.png "Himalaya")
*Tile of the TOCOMEE dataset (imported into Blender), created with geoTriMesh*



Options:
	dem
	dem_nodata=-9999
	orthophoto
	orthophoto_nodata=-9999,
	boundaries
	tiles_size=360
	tiles_bbox=(-180,-90,180,90)
	mesh_prefix='out'
	mesh_path=os.getcwd()
	mesh_shapefile=False
	scale_xy=0.000001
	z_exaggeration=1.0
	projection='orig'
	centering=True
	indexed_colors=True
	coloring_mode='orthophoto'
	mesh_format='x3d'):
