import gmsh

###############################################
###### CREATE .STL FILES FROM CAD MODEL ######
###############################################

# LOAD CAD AND INITIALIZE MESH

gmsh.initialize()
gmsh.option.setString(
    "Geometry.OCCTargetUnit", "MM"
)  # make sure gmsh reads .step file in mm, what CADQuery exports in
gmsh.model.add("hx_openfoam")

cad_file_path = "hx.step"

# entities = gmsh.model.occ.importShapes(cad_file_path)
gmsh.merge(cad_file_path)
gmsh.model.occ.synchronize()

# EXTRACT ALL VOLUMES
volumes = [e for e in gmsh.model.occ.getEntities() if e[0] == 3]

print(f"Extracted {len(volumes)} raw volumes from CAD.")

#### FRAGMENT VOLUMES & GENERATE SHARED SURFACES #####
print("Fragmenting volumes to define interfaces...")
gmsh.model.occ.fragment(volumes, [])
gmsh.model.occ.synchronize()

# FINAL VOLUMES AFTER FRAGMENT
final_volumes = gmsh.model.getEntities(dim=3)
print(f"Final number of volumes: {len(final_volumes)}")

# interface surfaces
surfaces = gmsh.model.getEntities(dim=2)

walls_coolant_interfaces = [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 24]
walls_breeder_interfaces = [
    3,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    33,
    34,
    35,
    36,
    37,
    38,
    39,
    41,
    42,
    43,
    44,
    46,
    47,
    48,
    49,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    79,
    80,
    81,
    82,
]
coolant_inlet = 85
coolant_outlet = 92
coolant_walls = [83, 84, 86, 87, 88, 89, 90, 91]
breeder_inlet = 97
breeder_outlet = 96
breeder_walls = [93, 94, 95]
walls = [1, 32, 40, 45, 70, 78]

walls_coolant_interface_marker = 4
walls_breeder_interface_marker = 5
coolant_inlet_marker = 6
coolant_outlet_marker = 7
coolant_walls_marker = 8
breeder_inlet_marker = 9
breeder_outlet_marker = 10
breeder_walls_marker = 11
walls_marker = 12

# .STL FILE GENERATION FOR OPENFOAM
gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 2)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # needed for openfoam
gmsh.option.setNumber(
    "Mesh.Binary", 0
)  # 1 for binary, 0 for ASCII, snappyHexMesh needs ASCII
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1)
gmsh.option.setNumber("Mesh.MeshSizeMax", 1)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(
    2
)  # need this line in order for the refinement to work on the surfaces

# wall coolant interfaces
gmsh.model.addPhysicalGroup(
    2,
    walls_coolant_interfaces,
    walls_coolant_interface_marker,
    name="walls_coolant_interface",
)

filename = "coolant_to_walls.stl"
gmsh.write(filename)
print(f"Exported {filename}")

surface_physical_groups = gmsh.model.getPhysicalGroups(2)
for group in surface_physical_groups:
    rm_tag = group[1]
    gmsh.model.removePhysicalGroups([(2, rm_tag)])


# wall breeder interfaces
gmsh.model.addPhysicalGroup(
    2,
    walls_breeder_interfaces,
    walls_breeder_interface_marker,
    name="walls_breeder_interface",
)

filename = "breeder_to_walls.stl"
gmsh.write(filename)
print(f"Exported {filename}")

surface_physical_groups = gmsh.model.getPhysicalGroups(2)
for group in surface_physical_groups:
    rm_tag = group[1]
    gmsh.model.removePhysicalGroups([(2, rm_tag)])

# coolant
gmsh.model.addPhysicalGroup(
    2, [coolant_inlet], coolant_inlet_marker, name="coolant_inlet"
)
gmsh.model.addPhysicalGroup(
    2, [coolant_outlet], coolant_outlet_marker, name="coolant_outlet"
)
gmsh.model.addPhysicalGroup(
    2, coolant_walls, coolant_walls_marker, name="coolant_walls"
)

filename = "coolant.stl"
gmsh.write(filename)
print(f"Exported {filename}")

surface_physical_groups = gmsh.model.getPhysicalGroups(2)
for group in surface_physical_groups:
    rm_tag = group[1]
    gmsh.model.removePhysicalGroups([(2, rm_tag)])


# breeder
gmsh.model.addPhysicalGroup(
    2, [breeder_inlet], breeder_inlet_marker, name="breeder_inlet"
)
gmsh.model.addPhysicalGroup(
    2, [breeder_outlet], breeder_outlet_marker, name="breeder_outlet"
)
gmsh.model.addPhysicalGroup(
    2, breeder_walls, breeder_walls_marker, name="breeder_walls"
)

filename = "breeder.stl"
gmsh.write(filename)
print(f"Exported {filename}")

surface_physical_groups = gmsh.model.getPhysicalGroups(2)
for group in surface_physical_groups:
    rm_tag = group[1]
    gmsh.model.removePhysicalGroups([(2, rm_tag)])

# walls
gmsh.model.addPhysicalGroup(2, walls, walls_marker, name="walls")
filename = "walls.stl"

gmsh.write(filename)
print(f"Exported {filename}")

surface_physical_groups = gmsh.model.getPhysicalGroups(2)
for group in surface_physical_groups:
    rm_tag = group[1]
    gmsh.model.removePhysicalGroups([(2, rm_tag)])

##### TAG & NAME PHYSICAL GROUPS #####
# need to open mesh in gmsh gui to determine the proper tagging as below

# volumes
walls_marker = 1
coolant_marker = 2
breeder_marker = 3

gmsh.model.addPhysicalGroup(3, [1], walls_marker, name=f"walls")
gmsh.model.addPhysicalGroup(3, [2], coolant_marker, name=f"coolant")
gmsh.model.addPhysicalGroup(3, [3], breeder_marker, name=f"breeder")

# surfaces between walls and coolant
gmsh.model.addPhysicalGroup(
    2,
    walls_coolant_interfaces,
    walls_coolant_interface_marker,
    name="walls_coolant_interface",
)

# surfaces between walls and breeder
gmsh.model.addPhysicalGroup(
    2,
    walls_breeder_interfaces,
    walls_breeder_interface_marker,
    name="walls_breeder_interface",
)


# other surfaces
coolant_inlet_marker = 6
coolant_outlet_marker = 7
breeder_inlet_marker = 8
breeder_outlet_marker = 9
walls_marker = 10

gmsh.model.addPhysicalGroup(2, walls, walls_marker, name="walls")
gmsh.model.addPhysicalGroup(
    2, [coolant_inlet], coolant_inlet_marker, name="coolant_inlet"
)
gmsh.model.addPhysicalGroup(
    2, [coolant_outlet], coolant_outlet_marker, name="coolant_outlet"
)
gmsh.model.addPhysicalGroup(
    2, [breeder_inlet], breeder_inlet_marker, name="breeder_inlet"
)
gmsh.model.addPhysicalGroup(
    2, [breeder_outlet], breeder_outlet_marker, name="breeder_outlet"
)


# ##### MESH SIZE & REFINEMENT #####
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)

##### SYNC & GENERATE MESH #####
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)

##### SAVE MESH #####
output_file = "hx.msh"
gmsh.write(output_file)
gmsh.finalize()
