import gmsh

###############################################################################
# MESH SIZE PARAMETERS
###############################################################################
max_element_size = 25.0

gmsh.initialize()
gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
gmsh.model.add("hx_openfoam")

cad_file_path = "/home/ckhurana/Fusion-Heat-Exchangers/meshing/shellTube/hx_fixed.brep"
gmsh.merge(cad_file_path)
gmsh.model.occ.synchronize()

###############################################################################
# EXTRACT AND FRAGMENT VOLUMES
###############################################################################
volumes = [e for e in gmsh.model.occ.getEntities() if e[0] == 3]
print(f"Extracted {len(volumes)} raw volumes from CAD.")
print(f"Volume tags: {[v[1] for v in volumes]}")

gmsh.model.occ.fragment(volumes, [])
gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

final_volumes = gmsh.model.getEntities(dim=3)
print(f"Final number of volumes: {len(final_volumes)}")
print(f"Final volume tags: {[v[1] for v in final_volumes]}")

surfaces = gmsh.model.getEntities(dim=2)
print(f"Number of surfaces: {len(surfaces)}")

###############################################################################
# VOLUME GROUPS  — adjust tags if fragmentation renumbers them
###############################################################################
shell_vols        = [1]
pipe_vols         = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
baffle_vols       = [12, 13, 14, 15, 16, 17]
inner_fluid_vols  = [18]
outer_fluid_vols  = [19]

###############################################################################
# BOUNDARY CONDITION SURFACES  — flat inlet/outlet faces
###############################################################################
inner_inlet  = [265]
inner_outlet = [264]
outer_inlet  = [267]
outer_outlet = [266]

###############################################################################
# HELPER: find shared surfaces between two lists of volumes
###############################################################################
def get_interface_surfaces(vol_tags_a, vol_tags_b):
    """
    Returns surface tags that bound BOTH a volume in group A
    and a volume in group B — i.e. the shared interface.
    """
    def surfaces_of(vol_tags):
        surfs = set()
        for tag in vol_tags:
            # getBoundary returns (dim, tag) pairs; take absolute tag (sign = orientation)
            bounds = gmsh.model.getBoundary([(3, tag)], oriented=False)
            surfs.update(abs(b[1]) for b in bounds)
        return surfs

    surfs_a = surfaces_of(vol_tags_a)
    surfs_b = surfaces_of(vol_tags_b)
    return list(surfs_a & surfs_b)   # intersection = shared faces

def get_outer_surfaces(vol_tags, exclude_surfaces):
    """
    Returns surfaces that bound the given volumes but are NOT
    in exclude_surfaces — i.e. the external wall faces.
    """
    all_surfs = set()
    for tag in vol_tags:
        bounds = gmsh.model.getBoundary([(3, tag)], oriented=False)
        all_surfs.update(abs(b[1]) for b in bounds)
    return list(all_surfs - set(exclude_surfaces))

###############################################################################
# COMPUTE INTERFACES
###############################################################################
print("Computing interfaces...")

# fluid-solid interfaces
iface_shell_outer_fluid  = get_interface_surfaces(shell_vols,   outer_fluid_vols)
iface_shell_inner_fluid  = get_interface_surfaces(shell_vols,   inner_fluid_vols)
iface_pipes_inner_fluid  = get_interface_surfaces(pipe_vols,    inner_fluid_vols)
iface_pipes_outer_fluid  = get_interface_surfaces(pipe_vols,    outer_fluid_vols)
iface_baffles_outer_fluid = get_interface_surfaces(baffle_vols, outer_fluid_vols)

# solid-solid interfaces
iface_shell_pipes        = get_interface_surfaces(shell_vols,   pipe_vols)
iface_shell_baffles      = get_interface_surfaces(shell_vols,   baffle_vols)
iface_pipes_baffles      = get_interface_surfaces(pipe_vols,    baffle_vols)

print(f"  shell  <-> outer_fluid : {len(iface_shell_outer_fluid)} surfaces")
print(f"  shell  <-> inner_fluid : {len(iface_shell_inner_fluid)} surfaces")
print(f"  pipes  <-> inner_fluid : {len(iface_pipes_inner_fluid)} surfaces")
print(f"  pipes  <-> outer_fluid : {len(iface_pipes_outer_fluid)} surfaces")
print(f"  baffles <-> outer_fluid: {len(iface_baffles_outer_fluid)} surfaces")
print(f"  shell  <-> pipes       : {len(iface_shell_pipes)} surfaces")
print(f"  shell  <-> baffles     : {len(iface_shell_baffles)} surfaces")
print(f"  pipes  <-> baffles     : {len(iface_pipes_baffles)} surfaces")

# outer shell wall (external surface — not shared with any fluid)
all_interface_surfs = (
    iface_shell_outer_fluid + iface_shell_inner_fluid +
    iface_pipes_inner_fluid + iface_pipes_outer_fluid +
    iface_baffles_outer_fluid +
    iface_shell_pipes + iface_shell_baffles + iface_pipes_baffles +
    inner_inlet + inner_outlet + outer_inlet + outer_outlet
)
shell_outer_wall = get_outer_surfaces(shell_vols, all_interface_surfs)
print(f"  shell outer wall       : {len(shell_outer_wall)} surfaces")

##############################################################################
# PHYSICAL GROUPS — VOLUMES
##############################################################################
gmsh.model.addPhysicalGroup(3, shell_vols,       tag=1,  name="solid_shell")
gmsh.model.addPhysicalGroup(3, pipe_vols,        tag=2,  name="solid_pipes")
gmsh.model.addPhysicalGroup(3, baffle_vols,      tag=3,  name="solid_baffles")
gmsh.model.addPhysicalGroup(3, inner_fluid_vols, tag=4,  name="fluid_1_tubeside")
gmsh.model.addPhysicalGroup(3, outer_fluid_vols, tag=5,  name="fluid_2_shellside")

###############################################################################
# PHYSICAL GROUPS — BOUNDARY CONDITIONS
###############################################################################
gmsh.model.addPhysicalGroup(2, inner_inlet,  tag=10, name="bc_inner_inlet")
gmsh.model.addPhysicalGroup(2, inner_outlet, tag=11, name="bc_inner_outlet")
gmsh.model.addPhysicalGroup(2, outer_inlet,  tag=12, name="bc_outer_inlet")
gmsh.model.addPhysicalGroup(2, outer_outlet, tag=13, name="bc_outer_outlet")

###############################################################################
# PHYSICAL GROUPS — INTERFACES
###############################################################################
if iface_shell_outer_fluid:
    gmsh.model.addPhysicalGroup(2, iface_shell_outer_fluid,   tag=20, name="iface_shell_outerfluid")
if iface_shell_inner_fluid:
    gmsh.model.addPhysicalGroup(2, iface_shell_inner_fluid,   tag=21, name="iface_shell_innerfluid")
if iface_pipes_inner_fluid:
    gmsh.model.addPhysicalGroup(2, iface_pipes_inner_fluid,   tag=22, name="iface_pipes_innerfluid")
if iface_pipes_outer_fluid:
    gmsh.model.addPhysicalGroup(2, iface_pipes_outer_fluid,   tag=23, name="iface_pipes_outerfluid")
if iface_baffles_outer_fluid:
    gmsh.model.addPhysicalGroup(2, iface_baffles_outer_fluid, tag=24, name="iface_baffles_outerfluid")
if iface_shell_pipes:
    gmsh.model.addPhysicalGroup(2, iface_shell_pipes,         tag=25, name="iface_shell_pipes")
if iface_shell_baffles:
    gmsh.model.addPhysicalGroup(2, iface_shell_baffles,       tag=26, name="iface_shell_baffles")
if iface_pipes_baffles:
    gmsh.model.addPhysicalGroup(2, iface_pipes_baffles,       tag=27, name="iface_pipes_baffles")
if shell_outer_wall:
    gmsh.model.addPhysicalGroup(2, shell_outer_wall,          tag=28, name="shell_outer_wall")

###############################################################################
# MESH SETTINGS
###############################################################################
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
# gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1)
# gmsh.option.setNumber("Mesh.Algorithm3D", 1)   # Delaunay tet
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 0)
###############################################################################
# GENERATE MESH
###############################################################################
try:
    gmsh.model.mesh.generate(3)
    print("3D mesh generation successful!")

    num_tets = len(gmsh.model.mesh.getElementsByType(4)[0])
    num_tris = len(gmsh.model.mesh.getElementsByType(2)[0])
    print(f"Generated {num_tets} tetrahedra")
    print(f"Generated {num_tris} triangles")

    if num_tets == 0:
        print("WARNING: No tetrahedra — check volume tags and physical groups")
    else:
        gmsh.model.mesh.optimize("Netgen")
        print("Mesh optimization complete!")

except Exception as e:
    print(f"Mesh generation failed: {e}")

###############################################################################
# SAVE
###############################################################################
output_file = "hx_fixed.msh"
gmsh.write(output_file)
print(f"Mesh saved to {output_file}")

num_tets = len(gmsh.model.mesh.getElementsByType(4)[0])
print(f"Total tetrahedra: {num_tets}")

# gmsh.fltk.run()
gmsh.finalize()