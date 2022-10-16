## Feature Surface

This program extracts inlet/outlet surfaces of 3D grids. Grids with any type of cell and various file extentions such as `*.vtk` and `*.vtu` can be used.

![AAA mesh with 11 featured surface detected.](https://raw.github.com/ameli/feature-surface/master/doc/figures/AAA-mesh.jpg "AAA mesh with 11 featured surface detected.")

## User Guide
Please refer to [user guide page](http://ameli.github.io/feature-surface/), or PDF documentation at [doc/UserGuide.pdf](https://github.com/ameli/feature-surface/raw/master/doc/UserGuide.pdf), or the [wiki](https://github.com/ameli/feature-surface/wiki/Feature-Surface) of the project.

## ParaView Plugin
In ParaView, from _tools_ menu open _Manage Plugins_. Then load [bin/libFeatureSurface.so]https://github.com/ameli/feature-surface/blob/master/bin/libFeatureSurfacePlugin.so?raw=true) file. Next, in _Filters_ menu go to _Extensions_ and apply _Feature Surface_ filter to your pipeline.

## License
Its free and open source under GNU/zlib license. Please see [License.txt](https://raw.github.com/ameli/feature-surface/master/License.txt) for terms.

## Author
Siavash Ameli  
[Shadden Research Group](http://shaddenlab.berkeley.edu/)  
University of California, Berkeley
