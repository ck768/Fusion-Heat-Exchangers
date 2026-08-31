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

    # READ OPENFOAM MESH

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

    # fluid fields
    if os.path.isdir(
        f"{results_folder}/checkpoints"
    ):  # check if checkpointing file exists
        coolant_checkpoint_file = Path(f"{results_folder}/checkpoints/coolant_checkpoint.bp")
        breeder_checkpoint_file = Path(f"{results_folder}/checkpoints/breeder_checkpoint.bp")
        openfoam_mesh = io4dolfinx.read_mesh(coolant_checkpoint_file, MPI.COMM_WORLD)
    else:
        os.mkdir(f"{results_folder}/checkpoints")

        ## coolant fields
        (
            p_coolant,
            u_coolant,
            T_coolant,
            openfoam_mesh,
            nut_coolant,
            facet_meshtags,
            volume_meshtags,
        ) = read_openfoam_data(
            parent_dir + openfoam_folder, subdomain="fluid_2_shellside"
        )
        # save openfoam data to checkpoint file
        save_openfoam_data_to_checkpoint(
            p_coolant,
            u_coolant,
            T_coolant,
            openfoam_mesh,
            nut_coolant,
            facet_meshtags,
            volume_meshtags,
            checkpoint_file=f"{results_folder}/checkpoints/coolant_checkpoint.bp",
        )

        ## breeder fields
        p_breeder, u_breeder, T_breeder, _, nut_breeder, _, _ = read_openfoam_data(
            parent_dir + openfoam_folder, subdomain="fluid_1_tubeside"
        )
        save_openfoam_data_to_checkpoint(
            p_breeder,
            u_breeder,
            T_breeder,
            openfoam_mesh,
            nut_breeder,
            facet_meshtags,
            volume_meshtags,
            checkpoint_file=f"{results_folder}/checkpoints/breeder_checkpoint.bp",
        )

    my_model = F.HydrogenTransportProblemDiscontinuous()
    my_model.mesh = F.Mesh(openfoam_mesh)
    mesh = my_model.mesh.mesh

    facet_meshtags = io4dolfinx.read_meshtags(breeder_checkpoint_file, openfoam_mesh, meshtag_name="facet_tags")
    volume_meshtags = io4dolfinx.read_meshtags(breeder_checkpoint_file, openfoam_mesh, meshtag_name="cell_tags")
    my_model.facet_meshtags = facet_meshtags
    my_model.volume_meshtags = volume_meshtags

    # VELOCITY INTERPOLATION
    el = element(
        "Lagrange",
        mesh.topology.cell_name(),
        1,
        shape=(mesh.geometry.dim,),
    )
    V_festim = fem.functionspace(mesh, el)

    V_CG1_vec = fem.functionspace(openfoam_mesh, ("DG", 0, (3,)))
    V_CG1 = fem.functionspace(openfoam_mesh, ("DG", 0))

    breeder_velocity = fem.Function(V_CG1_vec, name="U")
    io4dolfinx.read_function(breeder_checkpoint_file, breeder_velocity, name="U")

    u_openfoam = breeder_velocity
    V_openfoam = u_openfoam.function_space

    breeder_festim_velocity = fem.Function(V_festim)

    breeder_cells = my_model.volume_meshtags.find(1)  # breeder cells to interpolate to
    interpolation_data = fem.create_interpolation_data(
        V_to=V_festim, V_from=V_openfoam, cells=breeder_cells
    )
    breeder_festim_velocity.interpolate_nonmatching(
        u_openfoam, cells=breeder_cells, interpolation_data=interpolation_data
    )
    exit()

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

    # use inlet temps
    def D_fluid(T):
        return D_0_flibe * ufl.exp(-E_D_flibe / (F.k_B * T))

    # TURBULENT VISOCOSITY INTERPOLATION
    # V_festim = fem.functionspace(mesh, el)

    # V_CG1 = fem.functionspace(openfoam_mesh, ("DG", 0))

    # # coolant nut 
    # nut = fem.Function(V_CG1, name="nut")
    # io4dolfinx.read_function(coolant_checkpoint_file, nut, name="nut")
    # nut.x.array[nut.x.array < 0.0] = 0.0  # ensure no negative eddy viscosity

    # coolant_cells = my_model.volume_meshtags.find(2)  # coolant cells to interpolate to
    # interpolation_data = fem.create_interpolation_data(
    #     V_to=V_festim, V_from=V_openfoam, cells=coolant_cells
    # )

    # coolant_festim_nut = fem.Function(V_festim)
    # coolant_festim_nut.interpolate_nonmatching(
    #     nut, cells=coolant_cells, interpolation_data=interpolation_data
    # )

    # breeder nut 
    nut = fem.Function(V_CG1, name="nut")
    io4dolfinx.read_function(breeder_checkpoint_file, nut, name="nut")
    nut.x.array[nut.x.array < 0.0] = 0.0  # ensure no negative eddy viscosity

    interpolation_data = fem.create_interpolation_data(
        V_to=V_festim, V_from=V_openfoam, cells=breeder_cells
    )

    breeder_festim_nut = fem.Function(V_festim)
    breeder_festim_nut.interpolate_nonmatching(
        nut, cells=breeder_cells, interpolation_data=interpolation_data
    )

    exit()

    # coolant diffusivity 
    inlet_coolant_temp = 800

    Sc = 0.7
    D_turb_coolant = coolant_festim_nut / Sc

    D_expr = D_fluid(inlet_coolant_temp) + D_turb_coolant
    V = fem.functionspace(mesh, ("CG", 1))
    D_coolant = fem.Function(V)
    D_coolant.interpolate(fem.Expression(D_expr, V.element.interpolation_points))

    # breeder diffusivity 
    inlet_breeder_temp = 908

    Sc = 0.7
    D_turb_breeder = breeder_festim_nut / Sc

    D_expr = D_fluid(inlet_breeder_temp) + D_turb_coolant
    V = fem.functionspace(mesh, ("CG", 1))
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

    coolant_mat = F.Material(
        D=D_coolant,
        K_S_0=flibe_solubility.pre_exp.magnitude,
        E_K_S=flibe_solubility.act_energy.magnitude,
    )

    inconel_diffusivity = htm.diffusivities.filter(material=htm.Inconel)

    membrane_diffusivity = (
        htm.diffusivities.filter(material=htm.NIOBIUM)
        .filter(exclude=True, isotope="H")
        .filter(exclude=True, isotope="D")[0]
    )

    membrane_solubility = htm.solubilities.filter(material=htm.NIOBIUM)[0]
    membrane_solubility_law = "SIEVERT"

    wall_mat = F.Material(
        D_0=membrane_diffusivity.pre_exp.magnitude,
        E_D=membrane_diffusivity.act_energy.magnitude,
        K_S_0=membrane_solubility.pre_exp.to("mol/m**3/(Pa**0.5)").magnitude,
        E_K_S=membrane_solubility.act_energy.magnitude,
        solubility_law=membrane_solubility_law,
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

    # u = htm.ureg
    # membrane_recombo = nb_recomb

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

    ### interpolate openfoam temperature fields ###
    tdim = 3

    # breeder
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
    V_pipes = fem.functionspace(new_pipes_mesh, ("CG", 1))
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
        results_folder="hx_results",
    )

    dolfinx.log.set_log_level(dolfinx.log.LogLevel.INFO)

    my_model.run()
