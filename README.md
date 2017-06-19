# geoTriMesh

Converts GIS data to colored 3D meshes (X3D)

Usage:

	import geotrimesh

	elevation_mesh = geotrimesh.ElevationMesh()
	elevation_mesh.generate_mesh(dem='srtm3.tif', ortho='bluemarble.tif', bound='osm.shp', mesh_prefix='sample')


![alt text](./demodata/sample.png "Himalaya")
*Tile of the TOCOMEE dataset (imported into Blender), created with geoTriMesh*
