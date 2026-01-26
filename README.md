# Fusion-Heat-Exchangers
Multiphysics modeling of heat exchangers for fusion systems

NOTE: to use CadQuery you need to use another environment:

cq_environment.yml

To use OpenFOAM and FESTIM, use 

hx_environment.yml

To generate your mesh, run:

python meshing/shellTube/meshing_shell_tube.py

This should give an output .msh file, which name will be reflected in the Python meshing script. To load this mesh in OpenFOAM, change directory to go to:

openfoam/shellTube

Then run a command to copy the mesh:

cp ../../meshing/shellTube/[mesh_name].msh
gmshToFoam [mesh_name].msh

Usually the meshes from GMSH are in millimeters, so we can scale our mesh to meters using the following OpenFOAM command:

transformPoints -scale 0.001

Check your mesh using:

checkMesh

Then run the following to split your mesh into regions:
splitMeshRegions -cellZones -overwrite

Finally run this command to run your case:
chtMultiRegionFoam

To view results, use this command and view in Paraview:

touch case.foam


OpenFOAM details:
OpenFOAM version: 2506

