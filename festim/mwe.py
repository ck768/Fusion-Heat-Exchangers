import festim as F
import io4dolfinx
from mpi4py import MPI
import os
import sys
import dolfinx
from dolfinx import fem
from festim.helpers import nmm_interpolate
from dolfinx.mesh import (
    create_submesh,
)
from scifem import assemble_scalar
from basix.ufl import element
import ufl
from dolfinx.log import set_log_level, LogLevel
from Nb_recombination.Nb_recombination import nb_recomb
import h_transport_materials as htm
from openfoam_to_festim import read_openfoam_data, save_openfoam_data_to_checkpoint
from pathlib import Path


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

        vel_space = functionspace(
            mesh, self.velocity_field.function_space.ufl_element()
        )
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


def build_festim_model(openfoam_folder, results_folder):

    ### READ OPENFOAM DATA ###

    ## solid fields
    p, u, T_baffles, _, _, _, _ = read_openfoam_data(
        parent_dir + openfoam_folder, subdomain="solid_baffles"
    )
    p, u, T_pipes, _, _, _, _ = read_openfoam_data(
        parent_dir + openfoam_folder, subdomain="solid_pipes"
    )
    p, u, T_shell, _, _, _, _ = read_openfoam_data(
        parent_dir + openfoam_folder, subdomain="solid_shell"
    )

    ## coolant fields
    (
        p_coolant,
        u_coolant,
        T_coolant,
        openfoam_mesh,
        nut_coolant,
        facet_meshtags,
        volume_meshtags,
    ) = read_openfoam_data(parent_dir + openfoam_folder, subdomain="fluid_2_shellside")

    ## breeder fields
    p_breeder, u_breeder, T_breeder, _, nut_breeder, _, _ = read_openfoam_data(
        parent_dir + openfoam_folder, subdomain="fluid_1_tubeside"
    )

    my_model = F.HydrogenTransportProblemDiscontinuous()
    my_model.mesh = F.Mesh(openfoam_mesh)
    mesh = my_model.mesh.mesh
    my_model.facet_meshtags = facet_meshtags
    my_model.volume_meshtags = volume_meshtags

    ### VELOCITY INTERPOLATION ###
    
    # breeder
    entities_breeder = volume_meshtags.find(1) # breeder cells to interpolate to 
    breeder_submesh, _, _, _ = create_submesh(
        my_model.mesh.mesh, dim=mesh.topology.dim, entities=entities_breeder
    )
    V_festim = fem.functionspace(
        breeder_submesh, ("DG", 0, (mesh.geometry.dim,)) # max velocity matches with DG, is a bit off with CG
    )
    festim_breeder_velocity = fem.Function(V_festim)
    
    nmm_interpolate(f_in=u_breeder, f_out=festim_breeder_velocity)

    # coolant
    entities_coolant = volume_meshtags.find(2) # coolant cells to interpolate to 
    coolant_submesh, _, _, _ = create_submesh(
        my_model.mesh.mesh, dim=mesh.topology.dim, entities=entities_coolant
    )
    V_festim = fem.functionspace(
        coolant_submesh, ("DG", 0, (mesh.geometry.dim,)) # max velocity matches with DG, is a bit off with CG
    )
    festim_coolant_velocity = fem.Function(V_festim)
    
    nmm_interpolate(f_in=u_coolant, f_out=festim_coolant_velocity)

    ### NUT INTERPOLATION ###

    # breeder
    V_festim = fem.functionspace(breeder_submesh, ("DG", 0))
    festim_breeder_nut = fem.Function(V_festim)
    nmm_interpolate(f_in=nut_breeder, f_out=festim_breeder_nut)

    # coolant
    V_festim = fem.functionspace(coolant_submesh, ("DG", 0))
    festim_coolant_nut = fem.Function(V_festim)
    nmm_interpolate(f_in=nut_coolant, f_out=festim_coolant_nut)

    ## MATERIALS ##

    # flibe
    flibe_diffusivity = (
        htm.diffusivities.filter(material=htm.FLIBE)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )

    D_0_flibe = flibe_diffusivity.pre_exp.magnitude  # m2/s,
    E_D_flibe = flibe_diffusivity.act_energy.magnitude  # eV
    Sc = 0.7

    # use inlet temps
    def D_fluid(T):
        return D_0_flibe * ufl.exp(-E_D_flibe / (F.k_B * T))

    # breeder
    inlet_breeder_temp = 908
    D_turb_breeder = festim_breeder_nut / Sc

    D_expr = D_fluid(inlet_breeder_temp) + D_turb_breeder
    V = fem.functionspace(mesh, ("DG", 0))
    D_breeder = fem.Function(V)
    D_breeder.interpolate(fem.Expression(D_expr, V.element.interpolation_points))

    flibe_solubility = (
        htm.solubilities.filter(material=htm.FLIBE)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")
        .mean()
    )

    breeder_mat = F.Material(
        D=D_breeder,
        K_S_0=flibe_solubility.pre_exp.magnitude,
        E_K_S=flibe_solubility.act_energy.magnitude,
    )

    # breeder_mat = F.Material(
    #     D_0=1,
    #     E_D=0,
    #     K_S_0=2,
    #     E_K_S=0,
    # )

    coolant_mat = F.Material(
        D_0=2,
        E_D=0,
        K_S_0=1,
        E_K_S=0,
    )

    wall_mat = F.Material(
        D_0=3,
        E_D=0,
        K_S_0=1,
        E_K_S=0,
    )

    ## Defining subdomains and boundary markers from boundary summary
    # think the best way to do this is to define two different materials with the different corresponding D_nuts

    # breeder
    breeder_vol = F.VolumeSubdomain(id=1, material=breeder_mat)  # fluid_1_tubeside
    breeder_inlet = F.SurfaceSubdomain(id=2)  # bc_inner
    breeder_outlet = F.SurfaceSubdomain(id=1)

    # coolant
    coolant_vol = F.VolumeSubdomain(id=2, material=coolant_mat)  # fluid_2_shellside
    coolant_inlet = F.SurfaceSubdomain(id=7)
    coolant_outlet = F.SurfaceSubdomain(id=6)

    # walls
    shell = F.VolumeSubdomain(id=5, material=wall_mat)
    shell_wall = F.SurfaceSubdomain(id=18)

    pipes = F.VolumeSubdomain(id=4, material=wall_mat)
    baffles = F.VolumeSubdomain(id=3, material=wall_mat)

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
    penalty_term = 1e8
    my_model.interfaces = [
        F.Interface(id=22, subdomains=[coolant_vol, shell], penalty_term=penalty_term),
        F.Interface(id=23, subdomains=[coolant_vol, pipes], penalty_term=penalty_term),
        F.Interface(
            id=24, subdomains=[coolant_vol, baffles], penalty_term=penalty_term
        ),
        F.Interface(id=25, subdomains=[baffles, shell], penalty_term=penalty_term),
        F.Interface(id=26, subdomains=[breeder_vol, pipes], penalty_term=penalty_term),
        F.Interface(id=27, subdomains=[breeder_vol, shell], penalty_term=penalty_term),
        F.Interface(
            id=28, subdomains=[breeder_vol, baffles], penalty_term=penalty_term
        ),
        F.Interface(id=29, subdomains=[baffles, pipes], penalty_term=penalty_term),
    ]

    T = F.Species("T", subdomains=my_model.volume_subdomains)
    my_model.species = [T]

    my_model.temperature = lambda x: (
        400.0 + x[0]
    )  # placeholder, overridden by openfoam fields

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=coolant_inlet, value=0, species=T),
        F.FixedConcentrationBC(subdomain=breeder_inlet, value=2, species=T),
    ]

    advection_term_breeder = F.AdvectionTerm(
        velocity=festim_breeder_velocity,
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
        F.VTXSpeciesExport(
            filename=results_folder + "/inner_fluid.bp", field=T, subdomain=breeder_vol
        ),
        F.VTXSpeciesExport(
            filename=results_folder + "/outer_fluid.bp", field=T, subdomain=coolant_vol
        ),
        F.VTXSpeciesExport(
            filename=results_folder + "/shell.bp", field=T, subdomain=shell
        ),
        F.VTXSpeciesExport(
            filename=results_folder + "/pipes.bp", field=T, subdomain=pipes
        ),
        F.VTXSpeciesExport(
            filename=results_folder + "/baffles.bp", field=T, subdomain=baffles
        ),
    ]

    my_model.settings = F.Settings(atol=1e-8, rtol=1e-10, transient=False)
    my_model.petsc_options = {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }

    my_model.initialise()

    ### INTERPOLATE OPENFOAM TEMPERATURE FIELDS ###
    tdim = 3

    ## breeder
    T_breeder_sub = breeder_vol.sub_T
    entities_breeder = volume_meshtags.find(breeder_vol.id)
    new_breeder_mesh, _, _, _ = create_submesh(
        my_model.mesh.mesh, dim=tdim, entities=entities_breeder
    )
    V_breeder = fem.functionspace(new_breeder_mesh, ("CG", 1))
    T_breeder_from_openfoam = dolfinx.fem.Function(V_breeder)
    T_breeder_from_openfoam.interpolate(T_breeder)
    nmm_interpolate(f_out=T_breeder_sub, f_in=T_breeder_from_openfoam)

    # coolant
    T_coolant_sub = coolant_vol.sub_T
    entities_coolant = volume_meshtags.find(coolant_vol.id)
    new_coolant_mesh, _, _, _ = create_submesh(
        my_model.mesh.mesh, dim=tdim, entities=entities_coolant
    )
    V_coolant = fem.functionspace(new_coolant_mesh, ("CG", 1))
    T_coolant_from_openfoam = dolfinx.fem.Function(V_coolant)
    T_coolant_from_openfoam.interpolate(T_coolant)
    nmm_interpolate(f_out=T_coolant_sub, f_in=T_coolant_from_openfoam)

    # shell
    T_shell_sub = shell.sub_T
    entities_shell = volume_meshtags.find(shell.id)
    new_shell_mesh, _, _, _ = create_submesh(
        my_model.mesh.mesh, dim=tdim, entities=entities_shell
    )
    V_shell = fem.functionspace(new_shell_mesh, ("CG", 1))
    T_shell_from_openfoam = dolfinx.fem.Function(V_shell)
    T_shell_from_openfoam.interpolate(T_shell)
    nmm_interpolate(f_out=T_shell_sub, f_in=T_shell_from_openfoam)

    # pipes
    T_pipes_sub = pipes.sub_T
    entities_pipes = volume_meshtags.find(pipes.id)
    new_pipes_mesh, _, _, _ = create_submesh(
        my_model.mesh.mesh, dim=tdim, entities=entities_pipes
    )
    V_pipes = fem.functionspace(new_pipes_mesh, ("CG", 1)) # DG better
    T_pipes_from_openfoam = dolfinx.fem.Function(V_pipes)
    T_pipes_from_openfoam.interpolate(T_pipes)
    nmm_interpolate(f_out=T_pipes_sub, f_in=T_pipes_from_openfoam)

    # baffles
    T_baffles_sub = baffles.sub_T
    entities_baffles = volume_meshtags.find(baffles.id)
    new_baffles_mesh, _, _, _ = create_submesh(
        my_model.mesh.mesh, dim=tdim, entities=entities_baffles
    )
    V_baffles = fem.functionspace(new_baffles_mesh, ("CG", 1))
    T_baffles_from_openfoam = dolfinx.fem.Function(V_baffles)
    T_baffles_from_openfoam.interpolate(T_baffles)
    nmm_interpolate(f_out=T_baffles_sub, f_in=T_baffles_from_openfoam)

    return my_model


if __name__ == "__main__":
    my_model = build_festim_model(
        openfoam_folder="/openfoam/turbulentHX_nonconstantrho/hx.foam",
        results_folder="test",
    )

    dolfinx.log.set_log_level(dolfinx.log.LogLevel.INFO)

    my_model.run()
