# Fusion-Heat-Exchangers
Multiphysics modeling of heat exchangers for fusion systems

## Creating your geometry
To generate your mesh, run:

`python meshing/shellTube/meshing_shell_tube.py
`

This should give an output .msh file, whose name is specified in the Python meshing script. 

> [!NOTE]
> To use CadQuery (for geometry creation) you need to use another environment:
>
> `cq_environment.yml`
>
> To use OpenFOAM and FESTIM, use 
>
> `hx_environment.yml`

## Running your OpenFOAM simulation
To load this mesh in OpenFOAM, change directory to go to:

`openfoam/shellTube`

Then run a command to copy the mesh:

`cp ../../meshing/shellTube/[mesh_name].msh hx.msh`

`gmshToFoam [mesh_name].msh gmshToFoam.log 2>&1`

Check your mesh using:

`checkMesh`

Then run the following to split your mesh into regions:

`splitMeshRegions -cellZones -overwrite`

Finally run this command to run your case:

`chtMultiRegionFoam`

To view results, use this command and view in Paraview:

`touch case.foam`

To clean your case (remove data for all timsteps except initial), run:

`./Allclean`

> [!NOTE]
> If changing meshes, you need to delete the `polyMesh` folders for each domain:
>
> `rm -rf constant/breeder/polyMesh constant/coolant/polyMesh constant/walls/polyMesh`
>
> and then run the following to copy your new mesh and split the mesh regions:
>
> `cp ../../meshing/shellTube/[mesh_name].msh hx.msh`
>
> `gmshToFoam [mesh_name].msh gmshToFoam.log 2>&1`
>
> `splitMeshRegions -cellZones -overwrite`

> [!IMPORTANT]
> Using OpenFOAM version: 2506