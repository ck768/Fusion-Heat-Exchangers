# Fusion-Heat-Exchangers
Multiphysics modeling of heat exchangers for fusion systems

## How To Run 

Clone the repository: 

 ```
git clone https://github.com/ck768/Fusion-Heat-Exchangers
cd Fusion-Heat-Exchangers
```

If you do not have `mamba` installed, first download the [Miniforge installer for your system](https://conda-forge.org/download/). Then, run the following in your terminal, replacing the script name with the correct bash script for your system. Example uses Apple Silicon installer.

```
bash Miniforge3-Darwin-arm64.sh
```

Verify the installation by running:

```
mamba --version
```

This package requires two environments to run, one for CADQuery and one for OpenFOAM/FESTIM. To set up a new environment with the right dependencies for CADQuery, use: 

Then, set up a new environment with the right dependencies (e.g. dolfinx, FESTIM) using:

```
mamba create -f cadquery_env.yml 
```

Then, activate the environment to use CADQuery for geometry creation: 

```
mamba activate cadquery-env
```

Deactive the environment using 
```
conda deactivate
```

and set up a second environment with the right meshing, OpenFOAM, and FESTIM dependencies with: 
```
conda env create -f environment.yml
```

Activate this environment using 
```
conda activate fusion-hx-env
```
