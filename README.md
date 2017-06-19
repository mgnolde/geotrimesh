# geoTriMesh
Generator for 3D meshes in X3D format from GIS data

Usage:

	import geotrimesh

	elevation_mesh = geotrimesh.ElevationMesh()
	elevation_mesh.generate_mesh(dem='srtm3.tif', ortho='bluemarble.tif', bound='osm.shp', mesh_prefix='output')





