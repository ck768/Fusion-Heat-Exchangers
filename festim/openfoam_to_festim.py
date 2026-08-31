from foam2dolfinx import OpenFOAMReader
from dolfinx.io import VTXWriter
from mpi4py import MPI
import io4dolfinx


def read_openfoam_data(file_name, subdomain):
    """
    Read OpenFOAM data from a file and return the pressure, velocity, temperature, and viscosity (if it exists) fields.
    """
    print("Reading OpenFOAM data...")
    openfoam_reader = OpenFOAMReader(filename=file_name, cell_type=10)

    final_time = max(openfoam_reader.times)

    T = openfoam_reader.create_dolfinx_function_with_cell_data(
        t=final_time, name="T", subdomain=subdomain
    )

    try:  # read u, p fields if exist
        u = openfoam_reader.create_dolfinx_function_with_cell_data(
            t=final_time, name="U", subdomain=subdomain
        )
        p = openfoam_reader.create_dolfinx_function_with_cell_data(
            t=final_time, name="p", subdomain=subdomain
        )
    except Exception:
        u = None
        p = None

    try:
        # read turbulent viscosity if it exists
        nut = openfoam_reader.create_dolfinx_function_with_cell_data(
            t=final_time, name="nut", subdomain=subdomain
        )
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
            results_folder + "/pressure.bp",
            p,
            "BP5",
        )
        writer_u = VTXWriter(
            MPI.COMM_WORLD,
            results_folder + "/velocity.bp",
            u,
            "BP5",
        )
        writer_p.write(t=0)
        writer_u.write(t=0)

    writer_T = VTXWriter(
        MPI.COMM_WORLD,
        results_folder + "/temp.bp",
        T,
        "BP5",
    )

    writer_T.write(t=0)


def save_openfoam_data_to_checkpoint(
    T, checkpoint_file, p=None, u=None, mesh=None, nut=None, facet_meshtags=None, volume_meshtags=None, 
):

    if mesh:
        io4dolfinx.write_mesh(checkpoint_file, mesh)
    if u:
        io4dolfinx.write_function(checkpoint_file, u, time=0.0, name="U")
    if p:
        io4dolfinx.write_function(checkpoint_file, p, time=0.0, name="p")
    if nut:
        io4dolfinx.write_function(checkpoint_file, nut, time=0.0, name="nut")
    if facet_meshtags:
        io4dolfinx.write_meshtags(
            checkpoint_file, mesh, facet_meshtags, meshtag_name="facet_tags"
        )
        io4dolfinx.write_meshtags(
            checkpoint_file, mesh, volume_meshtags, meshtag_name="cell_tags"
        )

    io4dolfinx.write_function(checkpoint_file, T, time=0.0, name="T")


if __name__ == "__main__":
    pass
