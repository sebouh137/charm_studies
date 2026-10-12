#!/usr/bin/env python3
"""
convert_mvtx_gdml_to_gemc.py

Editable converter for the MVTX part of PHGeomUtility GDML -> GEMC TEXT input.

Version 3: adds an explicit MVTX envelope mother volume.  All generated
HalfStaves and StaveStructs are direct daughters of this envelope, and a single
user parameter translates the full MVTX assembly.

Default outputs:
  mvtx__geometry_original.txt
  mvtx__materials_original.txt
  mvtx_text.gcard

Dependencies:
  python3, numpy, scipy

Typical use:
  python3 convert_mvtx_gdml_to_gemc.py 'PHGeomUtility_geom (1).gdml' --outdir mvtx_text
  cd mvtx_text
  gemc mvtx_text.gcard

The most useful tuning knobs are in USER CONFIGURATION below.
"""

from __future__ import annotations

import argparse
import copy
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.spatial.transform import Rotation as SciRot

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Number of staves generated from the three GDML seed placements.
NSTAVES = [12, 16, 20]

# -----------------------------------------------------------------------------
# MVTX envelope placement
# -----------------------------------------------------------------------------
# The envelope uses the original GDML sol_MVTX_Wrapper polycone.  All HalfStaves
# and StaveStructs are placed as DIRECT daughters of this volume.  Their listed
# positions therefore remain local to the MVTX envelope.
#
# Change only this tuple to translate the complete MVTX rigidly in GEMC.
# Units are mm in the GEMC mother/root coordinate system.
MVTX_ENVELOPE_POSITION_MM = (0.0, 0.0, 85.0)

MVTX_ENVELOPE_NAME = "mvtx_envelope"
MVTX_ENVELOPE_MATERIAL = "ITS_AIR"
MVTX_ENVELOPE_VISIBLE = False

# Active silicon sensitivity.
MAKE_SENSORS_SENSITIVE = True
SENSITIVITY = "flux"
HIT_TYPE = "flux"

# Some GEMC builds do not define RECORD_PASSBY.  Leave it off by default.
# Turn this on only if your GEMC version supports that gcard option.
RECORD_PASSBY = False

# IMPORTANT GDML -> GEMC TEXT rotation convention.
# Geant4's GDML reader constructs physvol transforms with
# GetRotationMatrix(rotation).inverse(), whereas GEMC TEXT passes the parsed
# rotation matrix directly to G4PVPlacement.  Therefore ordinary GDML physvol
# rotations must be inverted when converted to GEMC TEXT.
INVERT_GDML_PHYSVOL_ROTATIONS = False

# For faster visualization, show only the 30 um active sensors.  All other
# volumes are still constructed and still participate in tracking/material.
SHOW_ONLY_SENSORS = False

# Make every HalfStave/chip/sensor explicit.  This gives unique manual IDs per
# sensor and avoids relying on CopyOf for sensitive logical volumes.
EXPLICIT_SENSORS = True
ADD_IDENTIFIERS = True

# Air logical envelopes are hidden by default.
HIDE_AIR_VOLUMES = True

# Short Boolean helper names avoid GEMC's historically fixed-length Operation:
# parsing.  Do not make these names long unless your GEMC parser is fixed.
BOOLEAN_HELPER_PREFIX = "b"
MAX_OPERATION_PAYLOAD = 190

# Sensor component ID written as the first GEMC identifier because the built-in
# flux hit process uses the first identifier as its component number.
def component_id(layer: int, stave: int, chip: int) -> int:
    return (layer + 1) * 10000 + (stave + 1) * 100 + (chip + 1)

# =============================================================================
# Small helpers
# =============================================================================

HEX_SUFFIX = re.compile(r"0x[0-9A-Fa-f]+$")


def strip_hex_suffix(name: str) -> str:
    name = HEX_SUFFIX.sub("", name)
    return re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_")


def unit_factor(unit: str) -> float:
    return {
        "mm": 1.0,
        "cm": 10.0,
        "m": 1000.0,
        "um": 0.001,
        "deg": 1.0,
        "rad": 180.0 / math.pi,
    }.get(unit, 1.0)


def attr_value(el, key: str, default_unit: str) -> float:
    if el is None:
        return 0.0
    return float(el.attrib.get(key, "0")) * unit_factor(el.attrib.get("unit", default_unit))


def gdml_position(parent) -> np.ndarray:
    p = parent.find("position")
    return np.array([attr_value(p, k, "mm") for k in ("x", "y", "z")], dtype=float)


def gdml_angles(parent) -> np.ndarray:
    r = parent.find("rotation")
    if r is None:
        return np.zeros(3)
    return np.array([attr_value(r, k, "deg") for k in ("x", "y", "z")], dtype=float)


def gdml_rotation(parent) -> SciRot:
    """Rotation matrix represented literally by a GDML <rotation>."""
    return SciRot.from_euler("xyz", gdml_angles(parent), degrees=True)


def gdml_physvol_rotation(parent) -> SciRot:
    """
    Rotation to write for an ordinary GDML physical-volume placement.

    IMPORTANT: use the GDML rotation LITERALLY.  Geant4's GDML reader first
    builds a G4Transform3D with GetRotationMatrix(rotation).inverse(), but the
    G4PVPlacement(G4Transform3D, ...) constructor inverts the transform rotation
    again before storing the physical-volume rotation.  GEMC TEXT constructs
    G4PVPlacement directly from the parsed G4RotationMatrix, so the same x/y/z
    angles as the GDML physvol give the same stored placement rotation.

    INVERT_GDML_PHYSVOL_ROTATIONS is retained only as a debugging switch.
    """
    r = gdml_rotation(parent)
    return r.inv() if INVERT_GDML_PHYSVOL_ROTATIONS else r


def rotate_placement_about_global_z(seed_position, seed_rotation, phi_deg):
    """Rigidly rotate an already-placed GDML physical volume about global +z.

    G4PVPlacement's stored rotation is a frame/passive rotation.  If P0 is the
    stored seed rotation and A0=P0^{-1} is the active local->mother rotation,
    then after an active global rotation Z(phi):

        A = Z(phi) A0
        P = A^{-1} = P0 Z(-phi)

    Thus the position rotates by +phi, while the stored placement matrix is
    right-multiplied by Z(-phi).  For a pure seed z rotation this reduces to
    z_text = z_seed - phi.
    """
    zplus = SciRot.from_euler("z", phi_deg, degrees=True)
    zminus = SciRot.from_euler("z", -phi_deg, degrees=True)
    new_position = zplus.apply(seed_position)
    new_rotation = seed_rotation * zminus
    return new_position, new_rotation


def compose_transform(t1, r1, t2, r2):
    """Apply transform 2, then transform 1."""
    return t1 + r1.apply(t2), r1 * r2


def fmt(x: float) -> str:
    if abs(x) < 5e-12:
        x = 0.0
    return f"{x:.12g}"


def pos_string(v) -> str:
    return " ".join(f"{fmt(x)}*mm" for x in v)


def rot_string(r: SciRot) -> str:
    a = r.as_euler("xyz", degrees=True)
    a[np.abs(a) < 5e-10] = 0.0
    return " ".join(f"{fmt(x)}*deg" for x in a)


# =============================================================================
# GEMC TEXT rows
# =============================================================================

class GemcRows:
    def __init__(self):
        self.rows: List[List[str]] = []
        self.names = set()

    def add(
        self, name, mother, description, position, rotation, color,
        type_, dimensions, material,
        magnetic_field="no", ncopy=1, pmany=1, exist=1,
        visible=1, style=1, sensitivity="no", hit_type="no",
        identifiers="no",
    ):
        if name in self.names:
            raise RuntimeError(f"Duplicate GEMC detector name: {name}")
        self.names.add(name)
        row = [
            name, mother, description, position, rotation, color,
            type_, dimensions, material, magnetic_field,
            str(ncopy), str(pmany), str(exist), str(visible), str(style),
            sensitivity, hit_type, identifiers,
        ]
        if len(row) != 18:
            raise RuntimeError("Internal error: geometry row does not have 18 fields")
        self.rows.append(row)

    def write(self, path: Path):
        path.write_text("\n".join(" | ".join(r) for r in self.rows) + "\n")


COLORS = {
    "ALPIDE_SI": "cc3333",
    "ALPIDE_METALSTACK": "cccccc",
    "ITS_ALUMINUM": "bbbbbb",
    "ITS_KAPTONPOLYCH2": "d4a017",
    "ITS_CarbonFleece": "444444",
    "ITS_K13D2U2k": "333333",
    "ITS_FGS003": "555555",
    "ITS_M55J6K": "333333",
    "ITS_M60J3K": "222222",
    "ITS_WATER": "3366cc",
    "ITS_INOX304": "999999",
    "ITS_PEEKCF30": "776655",
    "ITS_CERAMIC": "eee0c0",
    "ITS_GLUE": "c09060",
    "ITS_GLUE_IBFPC": "c09060",
    "ITS_AIR": "88aaff",
}


def color_for(material: str, component=False) -> str:
    return "aaaaaa" if component else COLORS.get(material, "77aa77")


# =============================================================================
# xtru -> G4GenericTrap helpers
# =============================================================================

def signed_area(poly):
    return 0.5 * sum(
        poly[i][0] * poly[(i + 1) % len(poly)][1]
        - poly[(i + 1) % len(poly)][0] * poly[i][1]
        for i in range(len(poly))
    )


def point_in_triangle(p, a, b, c):
    def cross(u, v, w):
        return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])
    c1, c2, c3 = cross(a, b, p), cross(b, c, p), cross(c, a, p)
    return ((c1 >= -1e-10 and c2 >= -1e-10 and c3 >= -1e-10)
            or (c1 <= 1e-10 and c2 <= 1e-10 and c3 <= 1e-10))


def triangulate_polygon(poly):
    poly = list(poly)
    indices = list(range(len(poly)))
    if signed_area(poly) < 0:
        indices.reverse()

    triangles = []
    guard = 0
    while len(indices) > 3 and guard < 1000:
        guard += 1
        clipped = False
        for ii in range(len(indices)):
            ia, ib, ic = indices[ii - 1], indices[ii], indices[(ii + 1) % len(indices)]
            a, b, c = poly[ia], poly[ib], poly[ic]
            z = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
            if z <= 1e-12:
                continue
            if any(j not in (ia, ib, ic) and point_in_triangle(poly[j], a, b, c)
                   for j in indices):
                continue
            triangles.append((a, b, c))
            indices.pop(ii)
            clipped = True
            break
        if not clipped:
            raise RuntimeError(f"Could not triangulate xtru polygon: {poly}")

    if len(indices) == 3:
        triangles.append(tuple(poly[j] for j in indices))
    return triangles


def geant4_clockwise_face(triangle):
    pts = [triangle[0], triangle[1], triangle[2], triangle[2]]

    def zorder(q):
        (x0, y0), (x1, y1), (x2, y2), (x3, y3) = q
        return (x2 - x0) * (y3 - y1) - (y2 - y0) * (x3 - x1)

    if zorder(pts) > 1e-12:
        pts[1], pts[3] = pts[3], pts[1]
    return pts


def generictrap_triangle_dimensions(triangle, dz):
    face = geant4_clockwise_face(triangle)
    vertices = face + list(face)
    vals = [f"{fmt(dz)}*mm"]
    for x, y in vertices:
        vals += [f"{fmt(x)}*mm", f"{fmt(y)}*mm"]
    return " ".join(vals)


# =============================================================================
# Converter
# =============================================================================

class Converter:
    def __init__(self, gdml_path: Path):
        import xml.etree.ElementTree as ET

        self.gdml_path = gdml_path
        self.root = ET.parse(gdml_path).getroot()
        self.structure = self.root.find("structure")
        self.solids_el = self.root.find("solids")
        self.materials_el = self.root.find("materials")

        self.volumes = {v.attrib["name"]: v for v in self.structure.findall("volume")}
        self.solids = {s.attrib["name"]: s for s in list(self.solids_el)}
        self.gdml_materials = {
            m.attrib["name"]: m for m in self.materials_el.findall("material")
        }

        self.half_name = next(n for n in self.volumes if n.startswith("MVTXHalfStave"))
        self.stave_name = next(n for n in self.volumes if n.startswith("MVTXStave_StaveStruct"))
        self.wrapper_name = next(n for n in self.volumes if n.startswith("log_MVTX_Wrapper"))
        self.wrapper = self.volumes[self.wrapper_name]
        self.seed_physvols = self.wrapper.findall("physvol")

        if len(self.seed_physvols) < 6:
            raise RuntimeError("Expected 6 seed placements in log_MVTX_Wrapper")

        self.rows = GemcRows()
        self.helper_counter = 0

    # ------------------------------- materials ------------------------------

    def material_name(self, gdml_ref):
        return strip_hex_suffix(gdml_ref)

    def source_material(self, volume_ref):
        return self.material_name(self.volumes[volume_ref].find("materialref").attrib["ref"])

    def mvtx_volume_closure(self):
        needed = set()
        stack = [self.half_name, self.stave_name]
        while stack:
            ref = stack.pop()
            if ref in needed:
                continue
            needed.add(ref)
            for pv in self.volumes[ref].findall("physvol"):
                stack.append(pv.find("volumeref").attrib["ref"])
        return needed

    def element_to_g4(self, ref):
        name = strip_hex_suffix(ref).upper()
        aliases = {
            "ALUMINIUM": "Al", "ALUMINUM": "Al", "OXYGEN": "O",
            "SILICON": "Si", "CARBON": "C", "HYDROGEN": "H",
            "NITROGEN": "N", "BARIUM": "Ba", "TITANIUM": "Ti",
            "ARGON": "Ar", "CHROMIUM": "Cr", "IRON": "Fe",
            "NICKEL": "Ni", "C": "C", "H": "H", "N": "N",
            "O": "O", "SI": "Si", "FE": "Fe", "CR": "Cr",
            "NI": "Ni", "AR": "Ar", "AL": "Al",
        }
        for key, symbol in aliases.items():
            if name == key or name.startswith(key + "_ELM") or name.startswith(key):
                return "G4_" + symbol
        raise RuntimeError(f"Unmapped GDML element: {ref}")

    def build_material_lines(self):
        refs = {
            self.volumes[v].find("materialref").attrib["ref"]
            for v in self.mvtx_volume_closure()
        }
        lines = []
        for ref in sorted(refs, key=self.material_name):
            clean = self.material_name(ref)
            m = self.gdml_materials[ref]
            density = float(m.find("D").attrib["value"])
            fractions = m.findall("fraction")

            if fractions:
                comps = [(self.element_to_g4(fr.attrib["ref"]), float(fr.attrib["n"]))
                         for fr in fractions]
            else:
                z = int(round(float(m.attrib["Z"])))
                z_to_symbol = {1: "H", 6: "C", 7: "N", 8: "O", 13: "Al", 14: "Si", 26: "Fe"}
                if z not in z_to_symbol:
                    raise RuntimeError(f"Unmapped elemental material Z={z}: {clean}")
                comps = [("G4_" + z_to_symbol[z], 1.0)]

            total = sum(f for _, f in comps)
            comps = [(name, f / total) for name, f in comps]
            comp_string = " ".join(f"{name} {f:.12g}" for name, f in comps)
            fields = [
                clean, f"from GDML {ref}", f"{density:.12g}", str(len(comps)),
                comp_string, "none", "none", "none", "none", "none",
            ]
            lines.append("|".join(fields))
        return lines

    # -------------------------------- solids --------------------------------

    def primitive_spec(self, solid_name):
        s = self.solids[solid_name]
        tag = s.tag
        lf = unit_factor(s.attrib.get("lunit", "mm"))
        af = unit_factor(s.attrib.get("aunit", "deg"))
        L = lambda k: float(s.attrib.get(k, "0")) * lf
        A = lambda k: float(s.attrib.get(k, "0")) * af

        if tag == "box":
            d = [L("x")/2, L("y")/2, L("z")/2]
            return "Box", " ".join(f"{fmt(x)}*mm" for x in d)
        if tag == "tube":
            d = [L("rmin"), L("rmax"), L("z")/2]
            a = [A("startphi"), A("deltaphi")]
            return "Tube", " ".join([f"{fmt(x)}*mm" for x in d] + [f"{fmt(x)}*deg" for x in a])
        if tag == "cone":
            d = [L("rmin1"), L("rmax1"), L("rmin2"), L("rmax2"), L("z")/2]
            a = [A("startphi"), A("deltaphi")]
            return "Cons", " ".join([f"{fmt(x)}*mm" for x in d] + [f"{fmt(x)}*deg" for x in a])
        if tag == "trd":
            d = [L("x1")/2, L("x2")/2, L("y1")/2, L("y2")/2, L("z")/2]
            return "Trd", " ".join(f"{fmt(x)}*mm" for x in d)
        if tag == "polycone":
            zp = s.findall("zplane")
            rmin = [float(z.attrib["rmin"])*lf for z in zp]
            rmax = [float(z.attrib["rmax"])*lf for z in zp]
            zz = [float(z.attrib["z"])*lf for z in zp]
            vals = [f'{fmt(A("startphi"))}*deg', f'{fmt(A("deltaphi"))}*deg', str(len(zp))]
            vals += [f"{fmt(x)}*mm" for x in rmin + rmax + zz]
            return "Polycone", " ".join(vals)
        raise KeyError(tag)

    def xtru_info(self, solid_name):
        s = self.solids[solid_name]
        poly = [(float(v.attrib["x"]), float(v.attrib["y"])) for v in s.findall("twoDimVertex")]
        sec = sorted(s.findall("section"), key=lambda e: int(e.attrib["zOrder"]))
        if len(sec) != 2:
            raise RuntimeError(f"{solid_name}: expected a two-section xtru")
        for q in sec:
            if (abs(float(q.attrib.get("scalingFactor", "1")) - 1) > 1e-12
                    or abs(float(q.attrib.get("xOffset", "0"))) > 1e-12
                    or abs(float(q.attrib.get("yOffset", "0"))) > 1e-12):
                raise RuntimeError(f"{solid_name}: nontrivial xtru scaling/offset not handled")
        z1, z2 = [float(q.attrib["zPosition"]) for q in sec]
        return poly, 0.5*(z1+z2), 0.5*abs(z2-z1)

    def intrinsic_offset(self, solid_name):
        s = self.solids[solid_name]
        if s.tag == "xtru":
            return np.array([0.0, 0.0, self.xtru_info(solid_name)[1]])
        if s.tag in ("union", "subtraction", "intersection"):
            return self.intrinsic_offset(s.find("first").attrib["ref"])
        return np.zeros(3)

    def new_helper_name(self):
        self.helper_counter += 1
        return f"{BOOLEAN_HELPER_PREFIX}{self.helper_counter:04d}"

    def check_operation(self, op):
        payload = op[10:] if op.startswith("Operation:") else op
        if len(payload) > MAX_OPERATION_PAYLOAD:
            raise RuntimeError(f"Operation payload too long ({len(payload)}): {op}")

    def compile_component(self, solid_name, t=None, r=None):
        if t is None:
            t = np.zeros(3)
        if r is None:
            r = SciRot.identity()
        s = self.solids[solid_name]

        if s.tag in ("box", "tube", "cone", "trd", "polycone"):
            name = self.new_helper_name()
            typ, dims = self.primitive_spec(solid_name)
            self.rows.add(name, "root", "Boolean helper", pos_string(t), rot_string(r),
                          color_for("Component", True), typ, dims, "Component", visible=0)
            return name

        if s.tag == "xtru":
            poly, zmid, dz = self.xtru_info(solid_name)
            tris = triangulate_polygon(poly)
            center_t, center_r = compose_transform(t, r, np.array([0.,0.,zmid]), SciRot.identity())
            pieces = []
            for tri in tris:
                name = self.new_helper_name()
                self.rows.add(name, "root", "xtru triangular prism", pos_string(center_t), rot_string(center_r),
                              color_for("Component", True), "G4GenericTrap",
                              generictrap_triangle_dimensions(tri, dz), "Component", visible=0)
                pieces.append(name)
            current = pieces[0]
            for other in pieces[1:]:
                name = self.new_helper_name()
                op = f"Operation:@ {current} + {other}"
                self.check_operation(op)
                self.rows.add(name, "root", "xtru polygon union", pos_string(center_t), rot_string(center_r),
                              color_for("Component", True), op, "0", "Component", visible=0)
                current = name
            return current

        if s.tag in ("union", "subtraction", "intersection"):
            first = s.find("first").attrib["ref"]
            second = s.find("second").attrib["ref"]
            op1 = self.compile_component(first, t, r)
            t2, r2 = compose_transform(t, r, gdml_position(s), gdml_rotation(s))
            op2 = self.compile_component(second, t2, r2)
            symbol = {"union": "+", "subtraction": "-", "intersection": "*"}[s.tag]
            result_t, result_r = compose_transform(t, r, self.intrinsic_offset(first), SciRot.identity())
            name = self.new_helper_name()
            op = f"Operation:@ {op1} {symbol} {op2}"
            self.check_operation(op)
            self.rows.add(name, "root", "GDML Boolean result", pos_string(result_t), rot_string(result_r),
                          color_for("Component", True), op, "0", "Component", visible=0)
            return name

        raise RuntimeError(f"Unsupported GDML solid: {s.tag} ({solid_name})")

    def actual_shape(self, solid_name):
        s = self.solids[solid_name]
        if s.tag in ("box", "tube", "cone", "trd", "polycone"):
            return self.primitive_spec(solid_name), np.zeros(3)

        if s.tag == "xtru":
            poly, zmid, dz = self.xtru_info(solid_name)
            tris = triangulate_polygon(poly)
            offset = np.array([0.,0.,zmid])
            if len(poly) <= 4:
                pts = list(poly)
                while len(pts) < 4:
                    pts.append(pts[-1])
                # Enforce clockwise G4GenericTrap faces.
                def zorder(q):
                    (x0,y0),(x1,y1),(x2,y2),(x3,y3)=q
                    return (x2-x0)*(y3-y1)-(y2-y0)*(x3-x1)
                if zorder(pts) > 1e-12:
                    pts[1], pts[3] = pts[3], pts[1]
                vals = [f"{fmt(dz)}*mm"]
                for x,y in pts + list(pts):
                    vals += [f"{fmt(x)}*mm", f"{fmt(y)}*mm"]
                return ("G4GenericTrap", " ".join(vals)), offset

            pieces = []
            for tri in tris:
                name = self.new_helper_name()
                self.rows.add(name, "root", "xtru triangular prism", pos_string(offset), "0*deg 0*deg 0*deg",
                              color_for("Component", True), "G4GenericTrap",
                              generictrap_triangle_dimensions(tri, dz), "Component", visible=0)
                pieces.append(name)
            current = pieces[0]
            for other in pieces[1:-1]:
                name = self.new_helper_name()
                op = f"Operation:@ {current} + {other}"
                self.check_operation(op)
                self.rows.add(name, "root", "xtru intermediate union", pos_string(offset), "0*deg 0*deg 0*deg",
                              color_for("Component", True), op, "0", "Component", visible=0)
                current = name
            op = f"Operation:@ {current} + {pieces[-1]}"
            self.check_operation(op)
            return (op, "0"), offset

        if s.tag in ("union", "subtraction", "intersection"):
            first = s.find("first").attrib["ref"]
            second = s.find("second").attrib["ref"]
            op1 = self.compile_component(first)
            op2 = self.compile_component(second, gdml_position(s), gdml_rotation(s))
            symbol = {"union": "+", "subtraction": "-", "intersection": "*"}[s.tag]
            op = f"Operation:@ {op1} {symbol} {op2}"
            self.check_operation(op)
            return (op, "0"), self.intrinsic_offset(first)

        raise RuntimeError(f"Unsupported solid type: {s.tag}")

    # --------------------------- recursive passive volume --------------------

    def define_volume(self, volume_ref, physical_name, mother, position, rotation,
                      copy_number=1, description=None, recurse=True):
        v = self.volumes[volume_ref]
        material = self.source_material(volume_ref)
        (typ, dims), offset = self.actual_shape(v.find("solidref").attrib["ref"])
        if np.linalg.norm(offset) > 1e-10 and v.findall("physvol"):
            raise RuntimeError(f"{volume_ref}: nonzero intrinsic solid offset with daughters")
        placed = position + rotation.apply(offset)
        self.rows.add(
            physical_name, mother, description or strip_hex_suffix(volume_ref),
            pos_string(placed), rot_string(rotation), color_for(material),
            typ, dims, material, ncopy=copy_number,
            visible=0 if (SHOW_ONLY_SENSORS or (HIDE_AIR_VOLUMES and material == "ITS_AIR")) else 1,
        )
        if recurse:
            for i, pv in enumerate(v.findall("physvol")):
                cref = pv.find("volumeref").attrib["ref"]
                cp = int(pv.attrib.get("copynumber", i))
                cname = f"{physical_name}_{strip_hex_suffix(cref)}_{cp:03d}_{i:03d}"
                self.define_volume(cref, cname, physical_name, gdml_position(pv), gdml_physvol_rotation(pv),
                                   copy_number=cp, description=f"{strip_hex_suffix(cref)} from GDML", recurse=True)

    # ------------------------------- geometry -------------------------------

    def add_mvtx_envelope(self):
        """Add the original GDML MVTX wrapper as one GEMC mother volume.

        The wrapper solid is a polycone in the source GDML.  The user-set
        MVTX_ENVELOPE_POSITION_MM translates only this mother; all HalfStave and
        StaveStruct coordinates remain local to it.
        """
        wrapper_solid_ref = self.wrapper.find("solidref").attrib["ref"]
        (etype, edims), eoff = self.actual_shape(wrapper_solid_ref)
        if np.linalg.norm(eoff) > 1e-10:
            raise RuntimeError(
                f"Unexpected intrinsic MVTX envelope solid offset: {eoff}"
            )

        epos = np.asarray(MVTX_ENVELOPE_POSITION_MM, dtype=float)
        if epos.shape != (3,):
            raise RuntimeError(
                "MVTX_ENVELOPE_POSITION_MM must contain exactly three values: "
                "(x_mm, y_mm, z_mm)"
            )

        self.rows.add(
            MVTX_ENVELOPE_NAME,
            "root",
            "MVTX envelope from source GDML wrapper",
            pos_string(epos),
            "0*deg 0*deg 0*deg",
            color_for(MVTX_ENVELOPE_MATERIAL),
            etype,
            edims,
            MVTX_ENVELOPE_MATERIAL,
            ncopy=1,
            visible=1 if MVTX_ENVELOPE_VISIBLE else 0,
        )

    def build_geometry(self):
        # Create the common mother first.  Every top-level HalfStave and
        # StaveStruct below is positioned in this envelope's local coordinates.
        self.add_mvtx_envelope()
        half_volume = self.volumes[self.half_name]
        module_pv = half_volume.find("physvol")
        module_ref = module_pv.find("volumeref").attrib["ref"]
        module_volume = self.volumes[module_ref]

        chip_pvs = [pv for pv in module_volume.findall("physvol")
                    if pv.find("volumeref").attrib["ref"].startswith("MVTXChip")]
        if len(chip_pvs) != 9:
            raise RuntimeError(f"Expected 9 MVTXChip placements, found {len(chip_pvs)}")
        chip_ref = chip_pvs[0].find("volumeref").attrib["ref"]
        chip_volume = self.volumes[chip_ref]

        sensor_pv = next(pv for pv in chip_volume.findall("physvol")
                         if pv.find("volumeref").attrib["ref"].startswith("MVTXSensor"))
        sensor_ref = sensor_pv.find("volumeref").attrib["ref"]
        metal_pv = next(pv for pv in chip_volume.findall("physvol")
                        if pv.find("volumeref").attrib["ref"].startswith("MetalStack"))
        metal_ref = metal_pv.find("volumeref").attrib["ref"]

        # Passive StaveStruct prototype; all other supports use CopyOf.
        stave_seed0 = self.seed_physvols[1]
        stave_proto = "mvtx_L0_S00_StaveStruct"
        self.define_volume(self.stave_name, stave_proto, MVTX_ENVELOPE_NAME, gdml_position(stave_seed0),
                           gdml_physvol_rotation(stave_seed0), copy_number=0,
                           description="MVTX stave support prototype", recurse=True)

        if not EXPLICIT_SENSORS:
            # Simple shared HalfStave version. Useful for debugging geometry only.
            half_seed0 = self.seed_physvols[0]
            half_proto = "mvtx_L0_S00_HalfStave"
            self.define_volume(self.half_name, half_proto, MVTX_ENVELOPE_NAME, gdml_position(half_seed0),
                               gdml_physvol_rotation(half_seed0), copy_number=0,
                               description="MVTX half-stave prototype", recurse=True)
            for layer, nstaves in enumerate(NSTAVES):
                hseed = self.seed_physvols[2*layer]
                sseed = self.seed_physvols[2*layer+1]
                for stave in range(nstaves):
                    if layer == 0 and stave == 0:
                        continue
                    phi = 360.0*stave/nstaves
                    hpos, hrot = rotate_placement_about_global_z(
                        gdml_position(hseed), gdml_physvol_rotation(hseed), phi
                    )
                    spos, srot = rotate_placement_about_global_z(
                        gdml_position(sseed), gdml_physvol_rotation(sseed), phi
                    )
                    hm = self.source_material(self.half_name)
                    self.rows.add(f"mvtx_L{layer}_S{stave:02d}_HalfStave", MVTX_ENVELOPE_NAME,
                                  f"Layer {layer} stave {stave} HalfStave",
                                  pos_string(hpos),
                                  rot_string(hrot), color_for(hm),
                                  f"CopyOf {half_proto}", "0", hm,
                                  ncopy=layer*100+stave, visible=0 if HIDE_AIR_VOLUMES else 1)
                    sm = self.source_material(self.stave_name)
                    self.rows.add(f"mvtx_L{layer}_S{stave:02d}_StaveStruct", MVTX_ENVELOPE_NAME,
                                  f"Layer {layer} stave {stave} StaveStruct",
                                  pos_string(spos),
                                  rot_string(srot), color_for(sm),
                                  f"CopyOf {stave_proto}", "0", sm,
                                  ncopy=layer*100+stave, visible=0 if HIDE_AIR_VOLUMES else 1)
            return

        # Shapes used by every explicit active hierarchy.
        (half_type, half_dims), half_off = self.actual_shape(half_volume.find("solidref").attrib["ref"])
        (module_type, module_dims), module_off = self.actual_shape(module_volume.find("solidref").attrib["ref"])
        (chip_type, chip_dims), chip_off = self.actual_shape(chip_volume.find("solidref").attrib["ref"])
        sensor_volume = self.volumes[sensor_ref]
        (sensor_type, sensor_dims), sensor_off = self.actual_shape(sensor_volume.find("solidref").attrib["ref"])
        metal_volume = self.volumes[metal_ref]
        (metal_type, metal_dims), metal_off = self.actual_shape(metal_volume.find("solidref").attrib["ref"])

        for label, off in [("HalfStave", half_off), ("Module", module_off), ("Chip", chip_off)]:
            if np.linalg.norm(off) > 1e-10:
                raise RuntimeError(f"Unexpected {label} intrinsic solid offset: {off}")

        passive_module_pvs = [pv for pv in module_volume.findall("physvol")
                              if not pv.find("volumeref").attrib["ref"].startswith("MVTXChip")]

        for layer, nstaves in enumerate(NSTAVES):
            hseed = self.seed_physvols[2*layer]
            sseed = self.seed_physvols[2*layer+1]
            for stave in range(nstaves):
                phi = 360.0*stave/nstaves
                hpos, hrot = rotate_placement_about_global_z(
                    gdml_position(hseed), gdml_physvol_rotation(hseed), phi
                )
                spos, srot = rotate_placement_about_global_z(
                    gdml_position(sseed), gdml_physvol_rotation(sseed), phi
                )

                hname = f"mvtx_L{layer}_S{stave:02d}_HalfStave"
                hmat = self.source_material(self.half_name)
                self.rows.add(hname, MVTX_ENVELOPE_NAME, f"Layer {layer} stave {stave} HalfStave",
                              pos_string(hpos),
                              rot_string(hrot), color_for(hmat),
                              half_type, half_dims, hmat, ncopy=layer*100+stave,
                              visible=0 if (HIDE_AIR_VOLUMES and hmat == "ITS_AIR") else 1)

                mname = f"{hname}_module"
                mmat = self.source_material(module_ref)
                self.rows.add(mname, hname, f"MVTX module L{layer} S{stave:02d}",
                              pos_string(gdml_position(module_pv)), rot_string(gdml_physvol_rotation(module_pv)),
                              color_for(mmat), module_type, module_dims, mmat,
                              visible=0 if (HIDE_AIR_VOLUMES and mmat == "ITS_AIR") else 1)

                # FPC, capacitors, resistors, etc.
                for i, pv in enumerate(passive_module_pvs):
                    cref = pv.find("volumeref").attrib["ref"]
                    cp = int(pv.attrib.get("copynumber", i))
                    cname = f"{mname}_{strip_hex_suffix(cref)}_{cp:03d}_{i:03d}"
                    self.define_volume(cref, cname, mname, gdml_position(pv), gdml_physvol_rotation(pv),
                                       copy_number=cp,
                                       description=f"{strip_hex_suffix(cref)} L{layer} S{stave:02d}", recurse=True)

                # 9 explicit ALPIDE chips and active sensors.
                for chip, cpv in enumerate(chip_pvs):
                    cname = f"{hname}_chip{chip}"
                    cmat = self.source_material(chip_ref)
                    self.rows.add(cname, mname, f"ALPIDE chip L{layer} S{stave:02d} C{chip}",
                                  pos_string(gdml_position(cpv)), rot_string(gdml_physvol_rotation(cpv)),
                                  color_for(cmat), chip_type, chip_dims, cmat,
                                  ncopy=chip+1, visible=0 if SHOW_ONLY_SENSORS else 1)

                    metname = f"{cname}_metal"
                    metmat = self.source_material(metal_ref)
                    self.rows.add(metname, cname, f"ALPIDE metal L{layer} S{stave:02d} C{chip}",
                                  pos_string(gdml_position(metal_pv)+metal_off),
                                  rot_string(gdml_physvol_rotation(metal_pv)), color_for(metmat),
                                  metal_type, metal_dims, metmat, visible=0 if SHOW_ONLY_SENSORS else 1)

                    sname = f"{cname}_sensor"
                    smat = self.source_material(sensor_ref)
                    sensitivity = SENSITIVITY if MAKE_SENSORS_SENSITIVE else "no"
                    hit_type = HIT_TYPE if MAKE_SENSORS_SENSITIVE else "no"
                    identifiers = "no"
                    if MAKE_SENSORS_SENSITIVE and ADD_IDENTIFIERS:
                        identifiers = (f"id manual {component_id(layer, stave, chip)} "
                                       f"layer manual {layer} stave manual {stave} chip manual {chip}")
                    self.rows.add(sname, cname, f"MVTX active sensor L{layer} S{stave:02d} C{chip}",
                                  pos_string(gdml_position(sensor_pv)+sensor_off),
                                  rot_string(gdml_physvol_rotation(sensor_pv)), color_for(smat),
                                  sensor_type, sensor_dims, smat, visible=1,
                                  sensitivity=sensitivity, hit_type=hit_type, identifiers=identifiers)

                # Support structure.
                if not (layer == 0 and stave == 0):
                    sm = self.source_material(self.stave_name)
                    self.rows.add(f"mvtx_L{layer}_S{stave:02d}_StaveStruct", MVTX_ENVELOPE_NAME,
                                  f"Layer {layer} stave {stave} StaveStruct",
                                  pos_string(spos),
                                  rot_string(srot), color_for(sm),
                                  f"CopyOf {stave_proto}", "0", sm,
                                  ncopy=layer*100+stave, visible=0 if HIDE_AIR_VOLUMES else 1)

    # ------------------------------- checks ---------------------------------

    def validate(self, material_lines):
        rows = self.rows.rows
        names = {r[0] for r in rows}
        if len(names) != len(rows):
            raise RuntimeError("Duplicate GEMC names")

        for r in rows:
            if r[1] != "root" and r[1] not in names:
                raise RuntimeError(f"{r[0]}: missing mother {r[1]}")
            if r[6].startswith("CopyOf "):
                target = r[6][7:].strip()
                if target not in names:
                    raise RuntimeError(f"{r[0]}: missing CopyOf target {target}")
            if r[6].startswith("Operation:"):
                if len(r[6][10:]) > MAX_OPERATION_PAYLOAD:
                    raise RuntimeError(f"{r[0]}: Operation too long")
                m = re.match(r"^Operation:(?:[@~])?\s*(.*?)\s+([+\-*])\s+(.*?)\s*$", r[6])
                if not m:
                    raise RuntimeError(f"{r[0]}: cannot parse operation {r[6]}")
                for operand in (m.group(1).strip(), m.group(3).strip()):
                    if operand not in names:
                        raise RuntimeError(f"{r[0]}: missing Boolean operand {operand}")

        mats = {line.split("|")[0].strip() for line in material_lines if line.strip()}
        for r in rows:
            if r[8] != "Component" and r[8] not in mats:
                raise RuntimeError(f"{r[0]}: undefined material {r[8]}")

        if MVTX_ENVELOPE_NAME not in names:
            raise RuntimeError(f"Missing MVTX envelope volume {MVTX_ENVELOPE_NAME}")

        top_level_re = re.compile(r"mvtx_L[012]_S\d\d_(HalfStave|StaveStruct)$")
        wrong_mother = [
            (r[0], r[1]) for r in rows
            if top_level_re.fullmatch(r[0]) and r[1] != MVTX_ENVELOPE_NAME
        ]
        if wrong_mother:
            raise RuntimeError(
                "Top-level MVTX staves are not daughters of the envelope: "
                + str(wrong_mother[:5])
            )

        if MAKE_SENSORS_SENSITIVE and EXPLICIT_SENSORS:
            sensors = [r for r in rows if r[0].endswith("_sensor") and r[15] == SENSITIVITY and r[16] == HIT_TYPE]
            expected = sum(NSTAVES) * 9
            if len(sensors) != expected:
                raise RuntimeError(f"Expected {expected} sensitive sensors, found {len(sensors)}")
            if ADD_IDENTIFIERS:
                ids = set()
                triples = set()
                for r in sensors:
                    m = re.fullmatch(r"id manual (\d+) layer manual (\d+) stave manual (\d+) chip manual (\d+)", r[17])
                    if not m:
                        raise RuntimeError(f"Bad identifier string: {r[0]} -> {r[17]}")
                    cid, L, S, C = map(int, m.groups())
                    if cid in ids or (L,S,C) in triples:
                        raise RuntimeError(f"Duplicate sensor identifier on {r[0]}")
                    ids.add(cid)
                    triples.add((L,S,C))

    # -------------------------------- write ---------------------------------

    def convert(self, outdir: Path):
        outdir.mkdir(parents=True, exist_ok=True)
        materials = self.build_material_lines()
        self.build_geometry()
        self.validate(materials)

        geom = outdir / "mvtx__geometry_original.txt"
        mats = outdir / "mvtx__materials_original.txt"
        gcard = outdir / "mvtx_text.gcard"

        self.rows.write(geom)
        mats.write_text("\n".join(materials) + "\n")
        passby = '  <option name="RECORD_PASSBY" value="1"/>\n' if RECORD_PASSBY else ""
        gcard.write_text(
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<gcard>\n'
            '  <detector name="mvtx" factory="TEXT" variation="original"/>\n'
            f'{passby}'
            '</gcard>\n'
        )

        sensors = [r for r in self.rows.rows if r[15] != "no"]
        operations = [r for r in self.rows.rows if r[6].startswith("Operation:")]
        print(f"Input GDML:            {self.gdml_path}")
        print(f"Geometry:              {geom}")
        print(f"Materials:             {mats}")
        print(f"Gcard:                 {gcard}")
        print(f"Geometry rows:         {len(self.rows.rows)}")
        print(f"Materials:             {len(materials)}")
        print(f"Sensitive sensors:     {len(sensors)}")
        print(f"Boolean operations:    {len(operations)}")
        print(f"Invert GDML physvol R: {INVERT_GDML_PHYSVOL_ROTATIONS}")
        print("Stave replication:      passive P0 * Rz(-phi)")
        print(f"Sensors-only display:  {SHOW_ONLY_SENSORS}")
        print(f"MVTX envelope mother:  {MVTX_ENVELOPE_NAME}")
        print(f"MVTX envelope position:{MVTX_ENVELOPE_POSITION_MM} mm")

        # Useful orientation sanity check for stave 0 in each layer.
        for layer in range(3):
            seed = self.seed_physvols[2*layer]
            gdml_z = gdml_angles(seed)[2]
            text_euler = gdml_physvol_rotation(seed).as_euler("xyz", degrees=True)
            print(f"L{layer} S00 GDML z={gdml_z:.6f} deg -> "
                  f"TEXT xyz=({text_euler[0]:.6f}, {text_euler[1]:.6f}, {text_euler[2]:.6f}) deg")
        if operations:
            print(f"Max operation payload: {max(len(r[6][10:]) for r in operations)}")
        print("Validation:            OK")


def main():
    p = argparse.ArgumentParser(description="Convert MVTX GDML to GEMC TEXT input")
    p.add_argument("gdml", type=Path, help="Original PHGeomUtility GDML file")
    p.add_argument("--outdir", type=Path, default=Path("mvtx_text"), help="Output directory")
    args = p.parse_args()
    if not args.gdml.exists():
        raise SystemExit(f"GDML file does not exist: {args.gdml}")
    Converter(args.gdml).convert(args.outdir)


if __name__ == "__main__":
    main()
