import gmsh

## Definining mesh sizes
bl_min = 3.0
bl_max = 15.0
bl_dist_min = 1.0
bl_dist_max = 5.0
curvature_factor = 25.0
max_element_size = 20.0

gmsh.initialize()
gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
gmsh.model.add("hx_openfoam")

cad_file_path = "hx.step"
gmsh.merge(cad_file_path)
gmsh.model.occ.synchronize()

# EXTRACT ALL VOLUMES
volumes = [e for e in gmsh.model.occ.getEntities() if e[0] == 3]
print(f"Extracted {len(volumes)} raw volumes from CAD.")
print(f"Volume tags: {[v[1] for v in volumes]}")

# FRAGMENT VOLUMES
print("Fragmenting volumes to define interfaces...")
gmsh.model.occ.fragment(volumes, [])
gmsh.model.occ.synchronize()

# CLEAN GEOMETRY
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

final_volumes = gmsh.model.getEntities(dim=3)
print(f"Final number of volumes: {len(final_volumes)}")
print(f"Final volume tags: {[v[1] for v in final_volumes]}")

surfaces = gmsh.model.getEntities(dim=2)
print(f"Number of surfaces: {len(surfaces)}")

# Surface IDs
walls_coolant_interfaces = [2,4,5,6,7,8,9,10,11,12,13,24]
walls_breeder_interfaces = [3,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,33,34,35,36,37,38,39,41,42,43,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,71,72,73,74,75,76,77,79,80,81,82]
coolant_inlet = 85
coolant_outlet = 92
coolant_walls = [83,84,86,87,88,89,90,91]
breeder_inlet = 97
breeder_outlet = 96
breeder_walls = [93,94,95]
walls = [1,32,40,45,70,78]

##### MESH SETTINGS FIRST - BEFORE PHYSICAL GROUPS #####

# Force deterministic meshing
gmsh.option.setNumber("General.NumThreads", 1)
gmsh.option.setNumber("Mesh.RandomFactor", 1e-9)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)

# Geometry tolerances
gmsh.option.setNumber("Geometry.Tolerance", 1e-5)
gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-5)

# Mesh size controls
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", curvature_factor)
# gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.01)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", max_element_size)

# Quality settings
gmsh.option.setNumber("Mesh.ElementOrder", 1)
gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)

# OpenFOAM compatibility
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 0)

##### MESH REFINEMENT #####

all_interface_surfaces = (
    walls_coolant_interfaces + 
    walls_breeder_interfaces + 
    walls + 
    coolant_walls +
    breeder_walls
)

distance_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance_field, "SurfacesList", all_interface_surfaces)

threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "InField", distance_field)
gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", bl_min)
gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", bl_max)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", bl_dist_min)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", bl_dist_max)

gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

##### SYNC BEFORE PHYSICAL GROUPS #####
gmsh.model.occ.synchronize()

##### TAG PHYSICAL GROUPS - CRITICAL: GET ACTUAL VOLUME TAGS #####

# Get actual volume tags after fragmentation
actual_volume_tags = [v[1] for v in final_volumes]
print(f"Creating physical groups for volumes: {actual_volume_tags}")

# IMPORTANT: Check if you have exactly 3 volumes, otherwise adjust
if len(actual_volume_tags) != 3:
    print(f"WARNING: Expected 3 volumes, got {len(actual_volume_tags)}")
    print("You may need to inspect the geometry in gmsh GUI to identify correct volume tags")

# Volumes - USE ACTUAL TAGS
walls_vol_marker = 1
coolant_vol_marker = 2
breeder_vol_marker = 3

# Assign physical groups to all volumes (CRITICAL!)
if len(actual_volume_tags) >= 3:
    gmsh.model.addPhysicalGroup(3, [actual_volume_tags[0]], walls_vol_marker, name="walls")
    gmsh.model.addPhysicalGroup(3, [actual_volume_tags[1]], coolant_vol_marker, name="coolant")
    gmsh.model.addPhysicalGroup(3, [actual_volume_tags[2]], breeder_vol_marker, name="breeder")
else:
    print("ERROR: Not enough volumes found after fragmentation!")
    # Create physical group with all volumes as fallback
    gmsh.model.addPhysicalGroup(3, actual_volume_tags, 1, name="all_volumes")

# External surfaces only
coolant_inlet_marker = 4
coolant_outlet_marker = 5
coolant_walls_marker = 6
breeder_inlet_marker = 7
breeder_outlet_marker = 8
breeder_walls_marker = 9
walls_surf_marker = 10
walls_breeder_interfaces_marker = 11
walls_coolant_interfaces_marker = 12

gmsh.model.addPhysicalGroup(2, [coolant_inlet], coolant_inlet_marker, name="coolant_inlet")
gmsh.model.addPhysicalGroup(2, [coolant_outlet], coolant_outlet_marker, name="coolant_outlet")
gmsh.model.addPhysicalGroup(2, coolant_walls, coolant_walls_marker, name="coolant_walls")
gmsh.model.addPhysicalGroup(2, [breeder_inlet], breeder_inlet_marker, name="breeder_inlet")
gmsh.model.addPhysicalGroup(2, [breeder_outlet], breeder_outlet_marker, name="breeder_outlet")
gmsh.model.addPhysicalGroup(2, breeder_walls, breeder_walls_marker, name="breeder_walls")
gmsh.model.addPhysicalGroup(2, walls, walls_surf_marker, name="walls_surface")
# gmsh.model.addPhysicalGroup(2, walls_breeder_interfaces, walls_breeder_interfaces_marker, name="walls_breeder_interfaces")
# gmsh.model.addPhysicalGroup(2, walls_coolant_interfaces, walls_coolant_interfaces_marker, name="walls_coolant_interfaces")

##### GENERATE MESH #####
gmsh.model.occ.synchronize()

try:
    gmsh.model.mesh.generate(3)
    print("3D mesh generation successful!")
    
    # Check what was actually meshed
    num_tets = len(gmsh.model.mesh.getElementsByType(4)[0])  # Type 4 = tetrahedra
    num_tris = len(gmsh.model.mesh.getElementsByType(2)[0])  # Type 2 = triangles
    print(f"Generated {num_tets} tetrahedra")
    print(f"Generated {num_tris} triangles")
    
    if num_tets == 0:
        print("\nWARNING: No tetrahedra generated!")
        print("This usually means volumes are not in physical groups or mesh failed")
        
        # Launch GUI to inspect
        print("Launching gmsh GUI for inspection...")
        gmsh.fltk.run()
    else:
        # Optimize mesh quality
        gmsh.model.mesh.optimize("Netgen")
        print("Mesh optimization complete!")
    
except Exception as e:
    print(f"Mesh generation failed: {e}")
    gmsh.fltk.run()

##### SAVE MESH #####
output_file = "hx_bl_MM.msh"
gmsh.write(output_file)
print(f"Mesh saved to {output_file}")

# Optional: Check mesh statistics
num_tets = len(gmsh.model.mesh.getElementsByType(4)[0])
print(f"Generated {num_tets} tetrahedra")
gmsh.fltk.run()

gmsh.finalize()