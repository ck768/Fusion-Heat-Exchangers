import festim as F
from foam2dolfinx import OpenFOAMReader, find_closest_value
from dolfinx.io import gmsh
from mpi4py import MPI
import numpy as np

## Reading mesh

mesh_path = "/Users/ckhurana/FESTIM/FESTIM-dev/openfoam/OpenFOAM/shellTubeHX/Fusion-Heat-Exchangers/openfoam/shellTube/hx.msh"

mesh_data = gmsh.read_from_msh(
    mesh_path, MPI.COMM_WORLD, 0, gdim=3
)

mesh = mesh_data.mesh

ft = mesh_data.facet_tags
ct = mesh_data.cell_tags

from dolfinx import plot
import pyvista

fdim = mesh.topology.dim - 1
tdim = mesh.topology.dim
mesh.topology.create_connectivity(fdim, tdim)
topology, cell_types, x = plot.vtk_mesh(mesh, fdim, ft.indices)



breeder_mat = F.Material(D_0=1e-3, E_D=0, K_S_0=10, E_K_S=0)
coolant_mat = F.Material(D_0=1e-3, E_D=0, K_S_0=10, E_K_S=0)
wall_mat = F.Material(D_0=1e-4, E_D=0, K_S_0=5, E_K_S=0)

breeder_vol = F.VolumeSubdomain(id=3, material=breeder_mat)
coolant_vol = F.VolumeSubdomain(id=2, material=coolant_mat)
walls_vol = F.VolumeSubdomain(id=1, material=wall_mat)

coolant_inlet_marker = 4
coolant_outlet_marker = 5
coolant_walls_marker = 6
breeder_inlet_marker = 7
breeder_outlet_marker = 8
breeder_walls_marker = 9
walls_surf_marker = 10
walls_breeder_interfaces_marker = 11
walls_coolant_interfaces_marker = 12

walls = F.SurfaceSubdomain(id=walls_surf_marker)
coolant_inlet = F.SurfaceSubdomain(id=coolant_inlet_marker)
coolant_outlet = F.SurfaceSubdomain(id=coolant_outlet_marker)
breeder_inlet = F.SurfaceSubdomain(id=breeder_inlet_marker)
breeder_outlet = F.SurfaceSubdomain(id=breeder_outlet_marker)
walls_breeder_interfaces = F.SurfaceSubdomain(id=walls_breeder_interfaces_marker)
walls_coolant_interfaces = F.SurfaceSubdomain(id=walls_coolant_interfaces_marker)

my_model = F.HydrogenTransportProblemDiscontinuous()

my_model.mesh = F.Mesh(mesh)

# we need to pass the meshtags to the model directly
my_model.facet_meshtags = ft
my_model.volume_meshtags = ct

my_model.subdomains = [
    walls_coolant_interfaces, 
    walls_breeder_interfaces, 
    coolant_outlet, 
    coolant_inlet, 
    breeder_outlet, 
    breeder_inlet, 
    walls, 
    walls_vol, 
    coolant_vol, 
    breeder_vol
]

my_model.surface_to_volume = {
    breeder_inlet: breeder_vol,
    breeder_outlet: breeder_vol,
    coolant_inlet: coolant_vol,
    coolant_outlet: coolant_vol,
    walls: walls_vol,
}

my_model.interfaces = [
    F.Interface(id=walls_breeder_interfaces_marker, subdomains=[breeder_vol, walls_vol], penalty_term=1e10),
    F.Interface(id=walls_coolant_interfaces_marker, subdomains=[coolant_vol, walls_vol], penalty_term=1e10),
]

H = F.Species("H", subdomains=[breeder_vol, coolant_vol, walls_vol])
my_model.species = [H]

# my_model.temperature = lambda t: get_my_T_field(t)
my_model.temperature = 400

my_model.boundary_conditions = [
    # F.FixedConcentrationBC(subdomain=outlet, value=0, species=H),
    # F.FixedConcentrationBC(subdomain=top_surface, value=1, species=H),
    F.FixedConcentrationBC(subdomain=breeder_inlet, value=1, species=H),
    F.FixedConcentrationBC(subdomain=coolant_inlet, value=2, species=H),
    # F.SievertsBC(subdomain=top_slab_interface, S_0=slab.K_S_0, E_S=slab.E_K_S, pressure=100, species=H),
    # F.SievertsBC(subdomain=bottom_slab_interface, S_0=slab.K_S_0, E_S=slab.E_K_S, pressure=100, species=H)
]

my_reader = OpenFOAMReader(filename="/Users/ckhurana/FESTIM/FESTIM-dev/openfoam/OpenFOAM/shellTubeHX/Fusion-Heat-Exchangers/openfoam/shellTube/case.foam", cell_type=10)
vel1 = my_reader.create_dolfinx_function(t=5, name="U", subdomain="breeder")
vel2 = my_reader.create_dolfinx_function(t=5, name="U", subdomain="coolant")

advection_term_breeder = F.AdvectionTerm(
    velocity=vel1,
    subdomain=breeder_vol,
    species=H,
)

advection_term_coolant = F.AdvectionTerm(
    velocity=vel2,
    subdomain=coolant_vol,
    species=H,
)

my_model.exports = [
    F.VTXSpeciesExport(filename="breeder.bp", field=H, subdomain=breeder_vol),
    F.VTXSpeciesExport(filename="coolant.bp", field=H, subdomain=coolant_vol),
    F.VTXSpeciesExport(filename="walls.bp", field=H, subdomain=walls_vol),
]

my_model.advection_terms = [advection_term_breeder, advection_term_coolant]
my_model.settings = F.Settings(atol=1e-10, rtol=1e-10, transient=False)

print("Initialising FESTIM model...")
my_model.initialise()

print("Running FESTIM solve...")
my_model.run()

print("FESTIM run completed.")