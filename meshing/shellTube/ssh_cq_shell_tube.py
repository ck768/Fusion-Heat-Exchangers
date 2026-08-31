import cadquery as cq

###############################################################################
# PARAMETERS
###############################################################################
R = 50
H = 500
fillet_r = 30
shell_thickness = 5

inner_fluid_r = 30
inner_fluid_depth = 50

outer_fluid_r = 15
outer_fluid_depth = 100
outer_fluid_inlet_z = H * 0.8
outer_fluid_outlet_z = H * 0.2

baffle_radius = R
baffle_thickness = 10
baffle_start = fillet_r + H / 20
baffle_end = H - H / 20 - fillet_r

hole_diameter = 15
col_spacing = baffle_radius / 2
row_spacing = baffle_radius / 2
holes_per_row = [3, 4, 3]

num_semi_baffles = 4
pipe_wall = 1

full_side_nozzle_length = R + shell_thickness + outer_fluid_depth
full_axial_nozzle_length = H + 2 * inner_fluid_depth  # spans entire body + both nozzles

###############################################################################
# HELPERS
###############################################################################


def generate_hole_positions(holes_per_row, col_spacing, row_spacing):
    holes = []
    n_rows = len(holes_per_row)
    row_offsets = [(i - (n_rows - 1) / 2) * row_spacing for i in range(n_rows)]
    for row_idx, n_holes in enumerate(holes_per_row):
        y = row_offsets[row_idx]
        col_offsets = [(i - (n_holes - 1) / 2) * col_spacing for i in range(n_holes)]
        for x in col_offsets:
            holes.append((x, y))
    return holes


###############################################################################
# 1.  MAIN SHELL BODY  (filleted outer minus filleted inner)
###############################################################################
outer_cyl = (
    cq.Workplane("front")
    .circle(R + shell_thickness)
    .extrude(H)
    .edges("front")
    .fillet(fillet_r)
    .edges("back")
    .fillet(fillet_r)
)

inner_bore = (
    cq.Workplane("front")
    .circle(R)
    .extrude(H)
    .edges("front")
    .fillet(fillet_r)
    .edges("back")
    .fillet(fillet_r)
)

shell_main = outer_cyl.cut(inner_bore)

###############################################################################
# 2.  AXIAL NOZZLES — conform-and-connect (same logic as side nozzles)
#
#  a) bore cylinder: full length from z=-inner_fluid_depth to z=H+inner_fluid_depth
#  b) shell around bore: same full length, hollowed
#  c) trim bore:  cut away portion INSIDE main inner bore  → only nozzle stub remains
#  d) trim shell: cut away portion INSIDE main outer cyl   → only nozzle wall remains
#     the trimmed shell wall has a curved face flush with the filleted end cap
#  e) fuse trimmed nozzle shells into shell_main
###############################################################################

# --- BOTTOM NOZZLE (-Z direction, starts below z=0) ---

# a) full bore from nozzle tip to well into main body
axial_bot_bore_full = (
    cq.Workplane("front")
    .workplane(offset=-inner_fluid_depth)
    .circle(inner_fluid_r)
    .extrude(inner_fluid_depth + fillet_r + shell_thickness)  # past fillet into body
)

# b) full shell wall
axial_bot_shell_full = (
    cq.Workplane("front")
    .workplane(offset=-inner_fluid_depth)
    .circle(inner_fluid_r + shell_thickness)
    .extrude(inner_fluid_depth + fillet_r + shell_thickness)
    .cut(axial_bot_bore_full)
)

# c) trim bore: remove portion inside main inner bore (inside the main cylinder)
axial_bot_bore_trimmed = axial_bot_bore_full.cut(inner_bore)

# d) trim shell: remove portion inside main outer cylinder
axial_bot_shell_trimmed = axial_bot_shell_full.cut(outer_cyl)

# --- TOP NOZZLE (+Z direction, starts above z=H) ---

axial_top_bore_full = (
    cq.Workplane("front")
    .workplane(offset=H - fillet_r - shell_thickness)
    .circle(inner_fluid_r)
    .extrude(inner_fluid_depth + fillet_r + shell_thickness)
)

axial_top_shell_full = (
    cq.Workplane("front")
    .workplane(offset=H - fillet_r - shell_thickness)
    .circle(inner_fluid_r + shell_thickness)
    .extrude(inner_fluid_depth + fillet_r + shell_thickness)
    .cut(axial_top_bore_full)
)

axial_top_bore_trimmed = axial_top_bore_full.cut(inner_bore)
axial_top_shell_trimmed = axial_top_shell_full.cut(outer_cyl)

# e) fuse into main shell
shell_with_axial_nozzles = shell_main.union(axial_bot_shell_trimmed).union(
    axial_top_shell_trimmed
)

###############################################################################
# 3.  SIDE NOZZLES — conform-and-connect
###############################################################################
full_side_nozzle_length = R + shell_thickness + outer_fluid_depth

# inlet (+Y)
inlet_bore_full = (
    cq.Workplane("XZ")
    .transformed(offset=(0, outer_fluid_inlet_z, 0))
    .circle(outer_fluid_r)
    .extrude(full_side_nozzle_length)
)
inlet_shell_full = (
    cq.Workplane("XZ")
    .transformed(offset=(0, outer_fluid_inlet_z, 0))
    .circle(outer_fluid_r + shell_thickness)
    .extrude(full_side_nozzle_length)
    .cut(inlet_bore_full)
)
inlet_bore_trimmed = inlet_bore_full.cut(inner_bore)
inlet_shell_trimmed = inlet_shell_full.cut(outer_cyl)

# outlet (-Y)
outlet_bore_full = (
    cq.Workplane("XZ")
    .transformed(offset=(0, outer_fluid_outlet_z, 0))
    .circle(outer_fluid_r)
    .extrude(-full_side_nozzle_length)
)
outlet_shell_full = (
    cq.Workplane("XZ")
    .transformed(offset=(0, outer_fluid_outlet_z, 0))
    .circle(outer_fluid_r + shell_thickness)
    .extrude(-full_side_nozzle_length)
    .cut(outlet_bore_full)
)
outlet_bore_trimmed = outlet_bore_full.cut(inner_bore)
outlet_shell_trimmed = outlet_shell_full.cut(outer_cyl)

# fuse all nozzle shells into one solid shell
solid_shell = shell_with_axial_nozzles.union(inlet_shell_trimmed).union(
    outlet_shell_trimmed
)

###############################################################################
# 4.  TUBE BUNDLE
###############################################################################
holes = generate_hole_positions(holes_per_row, col_spacing, row_spacing)
dz = baffle_end - baffle_start + baffle_thickness

tube_wall_parts = []
tube_bore_parts = []
tube_solid_parts = []

for x, y in holes:
    r_out = hole_diameter / 2
    r_in = r_out - pipe_wall
    tube_solid = (
        cq.Workplane("front")
        .workplane(offset=baffle_start)
        .center(x, y)
        .circle(r_out)
        .extrude(dz)
    )
    tube_bore = (
        cq.Workplane("front")
        .workplane(offset=baffle_start)
        .center(x, y)
        .circle(r_in)
        .extrude(dz)
    )
    tube_wall_parts.append(tube_solid.cut(tube_bore))
    tube_bore_parts.append(tube_bore)
    tube_solid_parts.append(tube_solid)

solid_pipes = tube_wall_parts[0]
for t in tube_wall_parts[1:]:
    solid_pipes = solid_pipes.union(t)

pipes_for_cut = tube_solid_parts[0]
for t in tube_solid_parts[1:]:
    pipes_for_cut = pipes_for_cut.union(t)


###############################################################################
# 5.  BAFFLES
###############################################################################
def make_baffle_disc(z_pos):
    return (
        cq.Workplane("front")
        .workplane(offset=z_pos)
        .circle(baffle_radius)
        .extrude(baffle_thickness)
        .cut(pipes_for_cut)
    )


def make_semi_baffle(index, z_pos):
    disc = (
        cq.Workplane("front")
        .workplane(offset=z_pos)
        .circle(baffle_radius)
        .extrude(baffle_thickness)
    )
    if index % 2 == 0:
        cut_box = (
            cq.Workplane("front")
            .workplane(offset=z_pos)
            .rect(2 * baffle_radius, 2 * baffle_radius)
            .extrude(baffle_thickness)
            .translate((0, baffle_radius / 2, 0))
        )
    else:
        cut_box = (
            cq.Workplane("front")
            .workplane(offset=z_pos)
            .rect(2 * baffle_radius, 2 * baffle_radius)
            .extrude(baffle_thickness)
            .translate((0, -baffle_radius / 2, 0))
        )
    return disc.intersect(cut_box).cut(pipes_for_cut)


baffle1 = make_baffle_disc(baffle_start)
baffle2 = make_baffle_disc(baffle_end)

semi_spacing = (baffle_end - baffle_start) / (num_semi_baffles + 1)
semi_list = [
    make_semi_baffle(i, baffle_start + (i + 1) * semi_spacing)
    for i in range(num_semi_baffles)
]

solid_baffles = baffle1.union(baffle2)
for sb in semi_list:
    solid_baffles = solid_baffles.union(sb)

###############################################################################
# 6.  FLUID 1 (tube side)
#     plenums + tube bores + axial nozzle bores (trimmed, outside main bore)
###############################################################################
interior = inner_bore

plenum_bot = interior.intersect(
    cq.Workplane("front")
    .workplane(offset=-inner_fluid_depth)
    .rect(2 * R, 2 * R)
    .extrude(inner_fluid_depth + baffle_start)
)
plenum_top = interior.intersect(
    cq.Workplane("front")
    .workplane(offset=baffle_end + baffle_thickness)
    .rect(2 * R, 2 * R)
    .extrude(H - (baffle_end + baffle_thickness) + inner_fluid_depth)
)

fluid1_tubes = tube_bore_parts[0]
for t in tube_bore_parts[1:]:
    fluid1_tubes = fluid1_tubes.union(t)

fluid_1 = (
    plenum_bot.union(plenum_top)
    .union(fluid1_tubes)
    .union(axial_bot_bore_trimmed)
    .union(axial_top_bore_trimmed)
)

###############################################################################
# 7.  FLUID 2 (shell side)
#     annular space + side nozzle bores (trimmed)
###############################################################################
fluid_2 = (
    interior.cut(plenum_bot)
    .cut(plenum_top)
    .cut(solid_baffles)
    .cut(pipes_for_cut)
    .union(inlet_bore_trimmed)
    .union(outlet_bore_trimmed)
)

solid_shell = solid_shell.cut(fluid_1).cut(fluid_2)
hx = cq.Assembly()
hx.add(solid_shell, name="solid_shell", color=cq.Color("gray"))
hx.add(solid_pipes, name="solid_pipes", color=cq.Color("blue"))
hx.add(solid_baffles, name="solid_baffles", color=cq.Color("red"))
hx.add(fluid_1, name="fluid_1", color=cq.Color("cyan", alpha=0.4))
hx.add(fluid_2, name="fluid_2", color=cq.Color("green", alpha=0.4))

###############################################################################
# 8.  EXPORT
###############################################################################

full_geometry = cq.Compound.makeCompound(
    [
        solid_shell.val(),
        solid_pipes.val(),
        solid_baffles.val(),
        fluid_1.val(),
        fluid_2.val(),
    ]
)

cq.exporters.export(full_geometry, "hx_fixed.brep")
