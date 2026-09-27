# Third-party data

## Natural Earth coastline

`src/plasic_app/resources/coastlines_110m.json` is derived from the
Natural Earth 1:110m Coastline dataset.

- Project: [Natural Earth](https://www.naturalearthdata.com/)
- Dataset: 1:110m physical vectors, coastline
- License: public domain
- Packaged form: 134 line segments and 5,128 longitude/latitude points

The coordinates are bundled with PlaSiC Studio so map and globe outlines work
offline. Cartopy, Shapely, and a network connection are not runtime
dependencies.
