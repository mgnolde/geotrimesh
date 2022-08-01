Generate sceneries from georeferenced 2D and 3D sources.

```python
from geotrimesh import GeoSceneSet

zurich = GeoSceneSet()

tilingscheme = GeoSceneSet.TilingScheme(boundary, dem_filepaths, height=256, width=256)

zurich.terrain = GeoSceneSet.Terrain(
    out_dirpath=out_dirpath,
    filepaths=dem_filepaths, 
    tiles=tilingscheme.tiles[2:3],
    boundary=boundary
    )

zurich.buildings = GeoSceneSet.Features("buildings",
    tilingscheme=tilingscheme,
    out_dirpath=out_dirpath,
    filepaths=buildings_filepaths,
    recombine_bodies=True,
    boundary=boundary, 
    tiles=tilingscheme.tiles[2:3],
    )

zurich.ortho = GeoSceneSet.Ortho(
    tilingscheme=tilingscheme,
    out_dirpath=out_dirpath,
    filepaths=ortho_filepaths,
    boundary=boundary, 
    tiles=tilingscheme.tiles[2:3],
    )
```

![zurich](docs/zurich_1_0.png)
![zurich](docs/tile_1_0_edges_2.png)
