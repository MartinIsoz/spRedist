# Solid phase redistribution algorithm for OpenFOAM

Custom implementation of a bounded scalar redistribution over a FVM mesh

Martin Isoz         (https://github.com/MartinIsoz)
Marie Plachá        (http://vscht.cz//monolith/)

Code capabilities
-----------------
* let alpha be a volume fraction (of a solid phase), alpha in <0,1>
* redistributes alpha from cells with alpha > 1 into the closest AVAILABLE neigbors
* redistribution is performed over processor-processor boundaries in parallel runs
* arbitrary stencil is available for alpha redistribution (not limited to face-neighbors
  of the cell with alpha > 1
* does not redistribute into cells with alpha >= 1 already
* suitable for coupling Eulerian CFD codes with Lagrangian particle trapping
* suitable for packing limiting in Eulerian-Eulerian sedimentation codes such as
  the sedFoam (https://github.com/SedFoam/)
  
Compatibility
-------------
The code is prepared for compilation with OpenFOAMv8 (https://openfoam.org/version/8/)

Compilation
-----------
Note: compilation has to be performed in a terminal with OpenFOAMv8 sourced

* navigate to ./src and run wclean && wmake

Cite this work as
-----------------

Isoz, M.; Placha, M.: A parallel algorithm for flux-based bounded scalar redistribution. In the proceedings of the conference Topical Problems of Fluid Mechanics 2022, TPFM 2022.
