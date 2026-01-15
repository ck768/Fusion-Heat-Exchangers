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

Then, set up a new environment with the right dependencies (e.g. dolfinx, FESTIM) using:

```
mamba create -f environment.yml 
```

Then, activate the environment: 

```
mamba activate fusion-hx-env
```

Now `cadquery ` and `jupter-cadquery` can be installed with: 

```
mamba install -c conda-forge -c cadquery cadquery=master
pip install jupyter-cadquery
pip install cadquery-ocp==7.7.2.0
pip install path.py
```