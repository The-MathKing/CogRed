# PyMOL script: overlay native (crystal) vs redocked best pose, for both
# 9L40 (failed redocking) and 9L4B (passed redocking), side by side.
bg_color white
set ray_opaque_background, 1
set orthoscopic, on

# --- Panel A: 9L40 (failed) ---
load 9l40_active_site_ligand_fixed.sdf, native_9l40
load 9l40_ATR_VE822_seed1000_redocked.sdf, docked_9l40
hide everything
show sticks, native_9l40
show sticks, docked_9l40
color grey30, native_9l40
color firebrick, docked_9l40
util.cnc native_9l40
util.cnc docked_9l40
color grey30, native_9l40 and elem C
color firebrick, docked_9l40 and elem C
orient native_9l40 or docked_9l40
zoom native_9l40 or docked_9l40, 2
set stick_radius, 0.15
ray 900, 700
png panel_9l40_overlay.png, dpi=300

# --- Panel B: 9L4B (passed) ---
delete all
load 9l4b_active_site_ligand_fixed.sdf, native_9l4b
load 9l4b_ATR_camonsertib_seed1000_redocked.sdf, docked_9l4b
hide everything
show sticks, native_9l4b
show sticks, docked_9l4b
color grey30, native_9l4b
color forest, docked_9l4b
util.cnc native_9l4b
util.cnc docked_9l4b
color grey30, native_9l4b and elem C
color forest, docked_9l4b and elem C
orient native_9l4b or docked_9l4b
zoom native_9l4b or docked_9l4b, 2
set stick_radius, 0.15
ray 900, 700
png panel_9l4b_overlay.png, dpi=300
