# Existing coordinates only; start PyMOL in this directory.
load apo_single_0.pdb, major
load apo_multi_3.pdb, alternative
load apo_multi_oracle.sdf, ligand
load ligand_ref_apo.pdb, reference
hide everything
show sticks, ligand or reference
show sticks, (major or alternative) within 5 of ligand
color marine, major
color red, alternative
color gray50, ligand
color green, reference
set stick_radius, 0.12
set transparency, 0.3, reference
zoom ligand, 8
bg_color white
