# Fusion-Heat-Exchangers
Multiphysics modeling of heat exchangers for fusion systems

NOTE: to use CadQuery you may need to use a separate environment. In particular install
cadquery (master)
jupyter-cadquery

OpenFOAM version: 2506

OpenFOAM workflow:
gmshToFoam partioned_cylinder_with_festim_interfaces.msh
splitMeshRegions -cellZones -overwrite
chtMultiRegionFoam
5. ./Allclean

