"""Imagens (PyMOL): visao geral com todas as poses e close do sitio de cada ligante.

Uso: pymol -cq figures.py -- OUT_DIR RECEPTOR_PDBQT COFACTOR_PDB|- NOME1=POSE1.pdbqt [NOME2=POSE2.pdbqt ...]
Um nome terminado em "_sitio2" (copia simetrica) usa a mesma cor do ligante de origem.
"""
import sys

from pymol import cmd, util

args = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
out, rec, cof = args[0], args[1], args[2]
ligs = [a.split("=", 1) for a in args[3:]]

CHAIN_COLORS = ["lightblue", "wheat", "palecyan", "lightpink", "palegreen", "lightorange"]
LIG_COLORS = ["lightorange", "violet", "tv_green", "hotpink", "slate", "deepteal", "olive", "salmon"]

cmd.bg_color("white")
for k, v in [("ray_opaque_background", 1), ("antialias", 2), ("ray_trace_mode", 1), ("ray_trace_color", "grey30"),
             ("cartoon_transparency", 0.15), ("stick_radius", 0.18), ("label_size", -0.85), ("label_font_id", 7),
             ("label_color", "black"), ("label_outline_color", "white"), ("label_position", (0, 0, 3)),
             ("dash_radius", 0.1), ("dash_gap", 0.25), ("depth_cue", 1), ("fog_start", 0.45), ("specular", 0.25)]:
    cmd.set(k, v)

cmd.load(rec, "rec")
cmd.remove("rec and hydro")
chains = cmd.get_chains("rec")
for i, c in enumerate(chains):
    cmd.color(CHAIN_COLORS[i % len(CHAIN_COLORS)], f"rec and chain {c}")
if cof != "-":
    cmd.load(cof, "cof")
    cmd.remove("cof and hydro")
bases = []
for name, path in ligs:
    cmd.load(path, name)
    base = name[:-len("_sitio2")] if name.endswith("_sitio2") else name
    if base not in bases:
        bases.append(base)
    cmd.color(LIG_COLORS[bases.index(base) % len(LIG_COLORS)], f"{name} and elem C")

# visao geral
cmd.hide("everything")
cmd.show("cartoon", "rec")
if cof != "-":
    cmd.show("sticks", "cof")
    cmd.set("stick_radius", 0.3, "cof")
    cmd.color("yelloworange", "cof and elem C")
for name, _ in ligs:
    cmd.show("spheres", f"{name} and not hydro")
    cmd.set("sphere_scale", 0.6, name)
util.cnc(" or ".join((["cof"] if cof != "-" else []) + [n for n, _ in ligs]))
cmd.orient("rec")
cmd.turn("x", -60)
cmd.zoom("rec", -6)
cmd.ray(1600, 1200)
cmd.png(f"{out}/00_visao_geral.png", dpi=300)

# sitio de cada ligante
for i, (name, _) in enumerate(ligs):
    cmd.hide("everything")
    for n2, _ in ligs:
        cmd.disable(n2)
    cmd.enable(name)
    lig = f"({name})"
    site = f"(byres ((rec) within 4.0 of {lig}))"
    cmd.show("cartoon", "rec")
    cmd.set("cartoon_transparency", 0.55)
    cmd.show("sticks", f"{site} and not name N+C+O")
    cmd.show("sticks", f"{site} and resn PRO and name N")
    cmd.color("grey70", f"{site} and elem C")
    cmd.show("sticks", f"{lig} and not (hydro and (neighbor elem C))")
    util.cnc(f"{site} or {lig}")
    cmd.distance(f"hb_{name}", f"{lig} and (elem N+O+Cl)", "rec and (elem N+O)", mode=0, cutoff=3.5)
    cmd.show("dashes", f"hb_{name}")
    cmd.hide("labels", f"hb_{name}")
    cmd.set("dash_color", "red", f"hb_{name}")
    cmd.label(f"({site} and name CB) or ({site} and resn GLY and name CA)", '"%s%s%s" % (resn.capitalize(), resi, chain)')
    cmd.orient(lig)
    cmd.zoom(f"{lig} or {site}", 3.0)
    cmd.ray(1600, 1200)
    cmd.png(f"{out}/{i + 1:02d}_{name}.png" if name.endswith("_sitio2") else f"{out}/{i + 1:02d}_{name}_sitio.png", dpi=300)
    cmd.delete(f"hb_{name}")
    cmd.label("all", "")
    cmd.set("cartoon_transparency", 0.15)
print("figuras em", out)
