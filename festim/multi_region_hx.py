import festim as F
from foam2dolfinx import OpenFOAMReader, find_closest_value
from dolfinx.io import gmsh, VTXWriter
from mpi4py import MPI
from dolfinx import plot
import os
import sys
import dolfinx
import numpy as np
from dolfinx import fem
from festim.helpers import nmm_interpolate
from dolfinx.mesh import meshtags, exterior_facet_indices
from scipy.spatial import cKDTree
from scifem import assemble_scalar
import ufl
from dolfinx.log import set_log_level, LogLevel

# add openfoam/ to path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, ".."))
sys.path.append(parent_dir)

class SurfaceAdvectionFlux(F.SurfaceFlux):
    """Computes the advection flux of a field on a given surface

    Args:
        field (festim.Species): species for which the surface flux is computed
        surface (festim.SurfaceSubdomain1D): surface subdomain
        filename (str, optional): name of the file to which the surface flux is exported

    Attributes:
        see `festim.SurfaceFlux`
    """

    def __init__(self, field, surface, filename, velocity_field):

        super().__init__(field=field, surface=surface, filename=filename)
        self.velocity_field = velocity_field

    @property
    def title(self):
        return f"{self.field.name} advection flux surface {self.surface.id}"

    def compute(self, u, ds: ufl.Measure, entity_maps=None):
        if isinstance(u, ufl.indexed.Indexed):
            mesh = self.field.sub_function_space.mesh
        else:
            mesh = u.function_space.mesh

        n = ufl.FacetNormal(mesh)

        # Diffusive flux
        surface_flux = assemble_scalar(
            fem.form(
                -self.D * ufl.dot(ufl.grad(u), n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )

        # Advective flux — interpolate velocity onto the submesh first
        from dolfinx.fem import Function, functionspace
        vel_space = functionspace(mesh, self.velocity_field.function_space.ufl_element())
        vel_local = Function(vel_space)
        vel_local.interpolate(self.velocity_field)

        advective_flux = assemble_scalar(
            fem.form(
                u * ufl.inner(vel_local, n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )

        self.value = surface_flux + advective_flux
        self.data.append(self.value)

def read_openfoam_data(file_name, subdomain):
    """
    Read OpenFOAM data from a file and return the pressure, velocity, temperature, and viscosity (if it exists) fields.
    """
    print("Reading OpenFOAM data...")
    openfoam_reader = OpenFOAMReader(filename=file_name, cell_type=10)

    final_time = max(openfoam_reader.times)

    T = openfoam_reader.create_dolfinx_function_with_point_data(t=final_time, name="T", subdomain=subdomain)

    try: # read u, p fields if exist
        u = openfoam_reader.create_dolfinx_function_with_point_data(t=final_time, name="U", subdomain=subdomain)
        p = openfoam_reader.create_dolfinx_function_with_point_data(t=final_time, name="p", subdomain=subdomain)
    except Exception: 
        u = None
        p = None

    try:
        # read turbulent viscosity if it exists
        nut = openfoam_reader.create_dolfinx_function_with_point_data(t=final_time, name="nut")
    except Exception:
        nut = None

    facet_meshtags = openfoam_reader.create_facet_meshtags()
    volume_meshtags = openfoam_reader.create_cell_meshtags()

    mesh = openfoam_reader.dolfinx_meshes_dict["_global"]

    return p, u, T, mesh, nut, facet_meshtags, volume_meshtags


def export_openfoam_data(p, u, T, results_folder):
    """
    Export OpenFOAM data to VTX files.
    """

    print("Exporting OpenFOAM data")
    if p != None and u != None:
        writer_p = VTXWriter(
            MPI.COMM_WORLD,
            results_folder+"/pressure.bp",
            p,
            "BP5",
        )
        writer_u = VTXWriter(
            MPI.COMM_WORLD,
            results_folder+"/velocity.bp",
            u,
            "BP5",
        )
        writer_p.write(t=0)
        writer_u.write(t=0)
    
    writer_T = VTXWriter(
        MPI.COMM_WORLD,
        results_folder+"/temp.bp",
        T,
        "BP5",
    )
    
    writer_T.write(t=0)


def build_festim_model(results_folder):

    ## coolant fields
    p_coolant, u_coolant, T_coolant, openfoam_mesh, nut_coolant, facet_meshtags, volume_meshtags = (
        read_openfoam_data(
            parent_dir+"/openfoam/steadyThickHX/hx.foam", subdomain="fluid_1_tubeside"
        )
    )

    ## breeder fields
    p_breeder, u_breeder, T_breeder, _, nut_breeder, _, _ = (
        read_openfoam_data(
            parent_dir+"/openfoam/steadyThickHX/hx.foam", subdomain="fluid_2_shellside"
        )
    )

    ## solid fields 
    p, u, T_baffles, _, _, _, _ = (
        read_openfoam_data(
            parent_dir+"/openfoam/steadyThickHX/hx.foam", subdomain="solid_baffles"
        )
    )
    p, u, T_soiid, _, _, _, _ = (
        read_openfoam_data(
            parent_dir+"/openfoam/steadyThickHX/hx.foam", subdomain="solid_pipes"
        )
    )
    p, u, T_shell, _, _, _, _ = (
        read_openfoam_data(
            parent_dir+"/openfoam/steadyThickHX/hx.foam", subdomain="solid_shell"
        )
    )

    # export_openfoam_data(p, u, T_baffles, results_folder=results_folder)

    my_model = F.HydrogenTransportProblemDiscontinuous()
    my_model.mesh = F.Mesh(openfoam_mesh)
    my_model.facet_meshtags = facet_meshtags
    my_model.volume_meshtags = volume_meshtags

    ## Defining material properties
    breeder_mat = F.Material(D_0=1e-3, E_D=0, K_S_0=10, E_K_S=0)
    coolant_mat = F.Material(D_0=1e-3, E_D=0, K_S_0=10, E_K_S=0)
    wall_mat = F.Material(D_0=1e-4, E_D=0, K_S_0=5, E_K_S=0)

    ## Defining subdomains and boundary markers from boundary summary

    # breeder
    breeder_vol = F.VolumeSubdomain(id=1, material=breeder_mat) # fluid_1_tubeside
    breeder_inlet = F.SurfaceSubdomain(id=2) # bc_inner
    breeder_outlet = F.SurfaceSubdomain(id=1)

    # coolant 
    coolant_vol = F.VolumeSubdomain(id=2, material=coolant_mat) # fluid_2_shellside
    coolant_inlet = F.SurfaceSubdomain(id=7)
    coolant_outlet = F.SurfaceSubdomain(id=6)
    
    # walls
    shell =  F.VolumeSubdomain(id=5, material=wall_mat)
    shell_wall = F.SurfaceSubdomain(id=18)

    pipes =  F.VolumeSubdomain(id=4, material=wall_mat)
    baffles =  F.VolumeSubdomain(id=3, material=wall_mat)

    my_model.subdomains = [ 
        breeder_vol, 
        breeder_inlet,
        breeder_outlet,
        coolant_vol, 
        coolant_inlet,
        coolant_outlet, 
        shell,
        shell_wall,
        pipes,
        baffles,
    ]

    my_model.surface_to_volume = {
        breeder_inlet: breeder_vol,
        breeder_outlet: breeder_vol,
        coolant_inlet: coolant_vol,
        coolant_outlet: coolant_vol,
        shell_wall: shell,
    }

    # interfaces
    my_model.interfaces = [
        F.Interface(id=22, subdomains=[coolant_vol, shell], penalty_term=1e5),
        F.Interface(id=23, subdomains=[coolant_vol, pipes], penalty_term=1e5),
        F.Interface(id=24, subdomains=[coolant_vol, baffles], penalty_term=1e5),
        F.Interface(id=25, subdomains=[baffles, shell], penalty_term=1e5),
        F.Interface(id=26, subdomains=[breeder_vol, pipes], penalty_term=1e5),
        F.Interface(id=27, subdomains=[breeder_vol, shell], penalty_term=1e5),
        F.Interface(id=28, subdomains=[breeder_vol, baffles], penalty_term=1e5),
        F.Interface(id=29, subdomains=[baffles, pipes], penalty_term=1e5),
    ]

    T = F.Species("T", subdomains=my_model.volume_subdomains)
    my_model.species = [T]

    my_model.temperature = 400

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=coolant_inlet, value=0, species=T),
        F.FixedConcentrationBC(subdomain=breeder_inlet, value=2, species=T),
    ]

    advection_term_breeder = F.AdvectionTerm(
        velocity=u_breeder,
        subdomain=breeder_vol,
        species=T,
    )

    advection_term_coolant = F.AdvectionTerm(
        velocity=u_coolant,
        subdomain=coolant_vol,
        species=T,
    )
    my_model.advection_terms = [advection_term_breeder, advection_term_coolant]

    my_model.exports = [
    F.VTXSpeciesExport(filename=results_folder+"/inner_fluid.bp",  field=T, subdomain=breeder_vol),
    F.VTXSpeciesExport(filename=results_folder+"/outer_fluid.bp",  field=T, subdomain=coolant_vol),
    F.VTXSpeciesExport(filename=results_folder+"/shell.bp",        field=T, subdomain=shell),
    F.VTXSpeciesExport(filename=results_folder+"/pipes.bp",        field=T, subdomain=pipes),
    F.VTXSpeciesExport(filename=results_folder+"/baffles.bp",      field=T, subdomain=baffles),
]

    my_model.settings = F.Settings(atol=1e-8, rtol=1e-10, transient=False)
    my_model.petsc_options = {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }

    return my_model

if __name__ == "__main__":

    my_model = build_festim_model(results_folder="hx_results")

    dolfinx.log.set_log_level(dolfinx.log.LogLevel.INFO)

    my_model.initialise()
    my_model.run() 