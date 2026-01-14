import gmsh
from collections import defaultdict

print("GMSH Version:", gmsh.__version__)

gmsh.initialize()
gmsh.option.setString("Geometry.OCCTargetUnit", "M")
gmsh.model.add("partitioned_cylinder")

# -----------------------------
# USER PARAMETERS
# -----------------------------
cad_file_path = "partitioned_cylinder.step"
H = 200.0
tol = 1e-6

# Physical IDs
top = 1
slab = 2
bottom = 3

top_inlet = 4
top_outlet = 2
bottom_inlet = 14
bottom_outlet = 11
top_slab_wall = 3
bottom_slab_wall = 13

# -----------------------------
# IMPORT CAD
# -----------------------------
gmsh.model.occ.importShapes(cad_file_path)
gmsh.model.occ.synchronize()

# -----------------------------
# FRAGMENT (CRITICAL)
# -----------------------------
volumes = gmsh.model.occ.getEntities(3)
print(f"Raw volumes: {len(volumes)}")

gmsh.model.occ.fragment(volumes, [])
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

final_volumes = gmsh.model.getEntities(3)
surfaces = gmsh.model.getEntities(2)

print(f"Final volumes: {len(final_volumes)}")
print(f"Total surfaces: {len(surfaces)}")

# -----------------------------
# BUILD volume → surface map
# -----------------------------
volume_surfaces = {}

for dim, vtag in final_volumes:
    bnd = gmsh.model.getBoundary([(3, vtag)], oriented=False, combined=False)
    volume_surfaces[vtag] = [s[1] for s in bnd if s[0] == 2]

# -----------------------------
# FIND SHARED SURFACES
# -----------------------------
surface_count = defaultdict(int)
for slist in volume_surfaces.values():
    for s in slist:
        surface_count[s] += 1

shared_surfaces = [s for s, c in surface_count.items() if c > 1]
print("Shared surfaces:", shared_surfaces)

# -----------------------------
# CLASSIFY VOLUMES
# -----------------------------
centroids = {
    v: gmsh.model.occ.getCenterOfMass(3, v)
    for _, v in final_volumes
}

top_vol = max(centroids, key=lambda v: centroids[v][1])
bottom_vol = min(centroids, key=lambda v: centroids[v][1])
slab_vol = list(set(centroids) - {top_vol, bottom_vol})[0]

gmsh.model.addPhysicalGroup(3, [top_vol], top, name="top_fluid")
gmsh.model.addPhysicalGroup(3, [slab_vol], slab, name="slab")
gmsh.model.addPhysicalGroup(3, [bottom_vol], bottom, name="bottom_fluid")

# -----------------------------
# CLASSIFY BOUNDARY SURFACES
# -----------------------------
gmsh.model.addPhysicalGroup(2, [top_inlet], name="top_inlet")
gmsh.model.addPhysicalGroup(2, [top_outlet], name="top_outlet")
gmsh.model.addPhysicalGroup(2, [bottom_inlet], name="bottom_inlet")
gmsh.model.addPhysicalGroup(2, [bottom_outlet], name="bottom_outlet")
gmsh.model.addPhysicalGroup(2, [top_slab_wall], name="top_slab_wall")
gmsh.model.addPhysicalGroup(2, [bottom_slab_wall], name="bottom_slab_wall")

# -----------------------------
# OPTIONAL: WALLS
# -----------------------------
wall_surfaces = [
    s for _, s in surfaces
    if s not in shared_surfaces
]

gmsh.model.addPhysicalGroup(2, wall_surfaces, 20, name="walls")

# -----------------------------
# MESH REFINEMENT PARAMETERS
# -----------------------------
h_bulk   = 5.0/100     # default mesh size
h_inlet  = 0.5/100     # near inlet/outlet
h_slab   = 0.3/100     # near slab interface
dist_min = 1.0/100
dist_max = 10.0/100

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h_slab)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h_bulk)

field_id = 1
fields = []

# -----------------------------
# INLET / OUTLET REFINEMENT
# -----------------------------
inlet_surfaces = [
    top_inlet,
    top_outlet,
    bottom_inlet,
    bottom_outlet,
]

gmsh.model.mesh.field.add("Distance", field_id)
gmsh.model.mesh.field.setNumbers(field_id, "FacesList", inlet_surfaces)
gmsh.model.mesh.field.setNumber(field_id, "Sampling", 100)

gmsh.model.mesh.field.add("Threshold", field_id + 1)
gmsh.model.mesh.field.setNumber(field_id + 1, "InField", field_id)
gmsh.model.mesh.field.setNumber(field_id + 1, "SizeMin", h_inlet)
gmsh.model.mesh.field.setNumber(field_id + 1, "SizeMax", h_bulk)
gmsh.model.mesh.field.setNumber(field_id + 1, "DistMin", dist_min)
gmsh.model.mesh.field.setNumber(field_id + 1, "DistMax", dist_max)

fields.append(field_id + 1)
field_id += 2

# -----------------------------
# SLAB INTERFACE REFINEMENT
# -----------------------------
gmsh.model.mesh.field.add("Distance", field_id)
gmsh.model.mesh.field.setNumbers(field_id, "FacesList", shared_surfaces)
gmsh.model.mesh.field.setNumber(field_id, "Sampling", 100)

gmsh.model.mesh.field.add("Threshold", field_id + 1)
gmsh.model.mesh.field.setNumber(field_id + 1, "InField", field_id)
gmsh.model.mesh.field.setNumber(field_id + 1, "SizeMin", h_slab)
gmsh.model.mesh.field.setNumber(field_id + 1, "SizeMax", h_bulk)
gmsh.model.mesh.field.setNumber(field_id + 1, "DistMin", dist_min)
gmsh.model.mesh.field.setNumber(field_id + 1, "DistMax", dist_max)

fields.append(field_id + 1)
field_id += 2

# -----------------------------
# VOLUME-BASED SIZING
# -----------------------------
gmsh.model.mesh.field.add("Constant", field_id)
gmsh.model.mesh.field.setNumber(field_id, "VIn", h_bulk)
gmsh.model.mesh.field.setNumbers(field_id, "VolumesList", [top_vol, bottom_vol])
fields.append(field_id)
field_id += 1

gmsh.model.mesh.field.add("Constant", field_id)
gmsh.model.mesh.field.setNumber(field_id, "VIn", h_slab)
gmsh.model.mesh.field.setNumbers(field_id, "VolumesList", [slab_vol])
fields.append(field_id)
field_id += 1

# -----------------------------
# COMBINE ALL FIELDS
# -----------------------------
gmsh.model.mesh.field.add("Min", field_id)
gmsh.model.mesh.field.setNumbers(field_id, "FieldsList", fields)
gmsh.model.mesh.field.setAsBackgroundMesh(field_id)


# -----------------------------
# MESH
# -----------------------------
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.mesh.generate(3)
gmsh.write("partitioned_cylinder_with_festim_interfaces.msh")

gmsh.fltk.run()
gmsh.finalize()
