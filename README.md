

This is the data repository of the article *Generation of 3D Representative Volume Elements (RVEs) of Nacre* submitted to Software Impacts.

### Contents

Two scripts are provided:
* `a01_Nacre_2D_RVE.m`: MatLab code for generating randomized layers of nacre (polygonal mineral
  platelets separated by a protein matrix) \
  geometry (size of RVE, layer thickness, matrix thickness, average platelet size) is parameterized
* `a02_Nacre_FULL_3D_RVE.py`: Python script to read in layer data and generate nacre representative volume element
    in <span style="font-variant:small-caps;">Abaqus</span> 


Additional functions for MatLab script are defined in .m files.


### External functions
   - Matlab AddOn "Mapping Toolbox"
   - inpoly2.m | INPOLY: A fast points-in-polygon test by Darren Engwirda \
     https://github.com/dengwirda/inpoly \
     *tested last with Version 3.0.0.0* 
   - VoronoiLimit.m | VoronoiLimit(vararg​in) by Jakob Sievers \
     https://de.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit-varargin \
     *tested last with Version 3.0.2.2* 
   - InterX.m | Curve intersections by NS \
     https://de.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections \
     *tested last with Version 1.5.0.0*


### Further information

For more information, the reader is referred to the article "Python codes to generate skeletal muscle models on each hierarchical level".

### Licence

This project is licensed under the terms of the MIT license.
