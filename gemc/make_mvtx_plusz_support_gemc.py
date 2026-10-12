#!/usr/bin/env python3
"""
make_mvtx_plusz_support_gemc.py

Generate the MVTX SUPPORT/ENCLOSURE MATERIAL on the +z side of the detector,
using the current sPHENIX PHG4MvtxSupport.cc definitions.

This is NOT the meter-long service barrel.  In the native sPHENIX coordinate
system that barrel is on -z.  The +z-side support material consists primarily of:

  * Three North aluminum end-wheel bases (layers 0,1,2)
  * Their stave-attachment steps
  * The North aluminum CYSS flange
  * The +z portion of the carbon-fiber CYSS cylinder

The script can model the drilled holes in the North end-wheel bases using GEMC
Boolean solids.  This is enabled by default.

Source:
  sPHENIX-Collaboration/coresoftware
  simulation/g4simulation/g4mvtx/PHG4MvtxSupport.cc
  PHG4MvtxDefs.h

Outputs:
  mvtx_plusz_support__geometry_original.txt
  mvtx_plusz_support__materials_original.txt
  mvtx_plusz_support.gcard
  MVTX_PLUSZ_SUPPORT_NOTES.txt

No third-party Python packages are required.
"""

from pathlib import Path
import argparse
import math


# =============================================================================
# USER SETTINGS
# =============================================================================

# Use "mvtx_envelope" if these rows are appended to the existing MVTX geometry.
# Use "root" to visualize/run this support geometry as a standalone detector.
MOTHER = "mvtx_envelope"

# Translation relative to the native MVTX coordinate system.
# Usually leave this at zero when MOTHER = "mvtx_envelope", because translating
# mvtx_envelope already translates the support together with the sensors.
PLUSZ_SUPPORT_POSITION_MM = (0.0, 0.0, 0.0)

VISIBLE = True

# Exact C++ end-wheel hole pattern via Boolean subtraction.
INCLUDE_ENDWHEEL_HOLES = True

# Exact C++ stave-attachment steps on the North end wheels.
INCLUDE_ENDWHEEL_STEPS = True

# The physical CYSS cylinder extends from negative z to positive z.
# For a deliberately "+z-only" material model, keep this True: the generated
# cylinder is clipped to z >= 0.  Set False to include the entire physical CYSS
# cylinder (-249.59-ish to +176.18-ish mm).
CLIP_CYSS_CYLINDER_AT_Z0 = False

COLORS = {
    "al": "b0b0b0",
    "cf": "555555",
    "component": "aaaaaa5",
}


# =============================================================================
# CONSTANTS FROM PHG4MvtxSupport.cc / PHG4MvtxDefs.h
# =============================================================================

# ServiceProperties
END_WHEEL_SN_HOLES_ZDIST = 308.0
END_W_STEP_HOLE_ZPOS = 4.0
END_W_STEP_HOLE_ZDIST = 4.0
END_WHEEL_N_LEN = 30.0

CYSS_FLG_SSF_FLG_NSF = 606.59
CYSS_FLG_SSF_CYLSF = 181.0

# PHG4MvtxDefs::mvtxdat
# layer: rMin, rMid, rMax, chips/stave, phi0 [rad], nStaves
MVTXDAT = [
    (24.610, 25.230, 27.930, 9, 0.2382, 12),
    (31.995, 33.360, 36.258, 9, 0.1937, 16),
    (39.930, 41.480, 44.260, 9, 0.1481, 20),
]

# North End Wheels
ENDW_N_RMAX = [30.57, 38.59, 46.30]
ENDW_N_RMIN = [22.83, 30.67, 38.70]
ENDW_N_WALL_R = [29.57, 37.59, 45.30]
ENDW_N_THICK = 2.8

ENDW_N_BASE_N_BIG_HOLES = 6
ENDW_N_BASE_N_SMALL_HOLES = 5
ENDW_N_BASE_BIG_HOLE_D = 4.0
ENDW_N_BASE_SMALL_HOLE_D = 1.6
ENDW_N_BASE_HOLES_R = [27.0, 35.0, 42.7]
ENDW_N_BASE_HOLES_PHI_DEG = 15.0

ENDW_N_WALL_NHOLES = [6, 8, 10]       # C++ loops over 2*N
ENDW_N_WALL_HOLE_D = 4.5

ENDW_N_STEP_XLEN = [9.88, 11.2, 8.38]
ENDW_N_STEP_YDISPL = [26.521, 33.891, 41.64]
ENDW_N_STEP_HOLE_PHI_DEG = [30.0, 22.5, 18.0]
ENDW_N_STEP_HOLE_TILT_RAD = [0.0, 0.295, 0.297]

# CYSS North flange
CYSS_FLG_EXT_THICK = 1.0
CYSS_FLG_INT_THICK = 1.5
CYSS_FLG_S_RING_LEN = 4.5
CYSS_FLG_L_RING_LEN = 7.5
CYSS_FLG_RMIN = 23.6
CYSS_FLG_RMAX = 52.8
CYSS_FLG_S_RING_RMIN = 46.4
CYSS_FLG_S_RING_RMAX = 47.4
CYSS_FLG_L_RING_RMIN = 50.6
CYSS_FLG_L_RING_RMAX = 51.5

# CYSS cylinder
CYSS_CYL_RMIN = 51.70
CYSS_CYL_RMAX = 52.82
CYSS_CYL_LEN = 425.77

# Materials from CreateMvtxSupportMaterials()
MVTX_AL_DENSITY = 2.7
MVTX_CF_DENSITY = 1.987
EPOXY_ATOMS = {"H": 32, "C": 2, "N": 4, "O": 15}
ATOMIC_MASS = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994}


# =============================================================================
# GEMC TEXT HELPERS
# =============================================================================

def fmt(x):
    if abs(x) < 5e-12:
        x = 0.0
    return f"{x:.12g}"


def pos(x, y, z):
    return f"{fmt(x)}*mm {fmt(y)}*mm {fmt(z)}*mm"


def rot(x=0.0, y=0.0, z=0.0):
    return f"{fmt(x)}*deg {fmt(y)}*deg {fmt(z)}*deg"


def tube_dims(rmin, rmax, halfz, phi0=0.0, dphi=360.0):
    return (
        f"{fmt(rmin)}*mm {fmt(rmax)}*mm {fmt(halfz)}*mm "
        f"{fmt(phi0)}*deg {fmt(dphi)}*deg"
    )


def box_dims(hx, hy, hz):
    return f"{fmt(hx)}*mm {fmt(hy)}*mm {fmt(hz)}*mm"


def polycone_dims(phi0, dphi, zplanes, rmins, rmaxs):
    n = len(zplanes)
    assert n == len(rmins) == len(rmaxs)
    out = [f"{fmt(phi0)}*deg", f"{fmt(dphi)}*deg", str(n)]
    out += [f"{fmt(v)}*mm" for v in rmins]
    out += [f"{fmt(v)}*mm" for v in rmaxs]
    out += [f"{fmt(v)}*mm" for v in zplanes]
    return " ".join(out)


def row(
    name, mother, desc, position, rotation, color,
    typ, dimensions, material, visible=True, ncopy=1
):
    fields = [
        name, mother, desc, position, rotation, color,
        typ, dimensions, material,
        "no",                 # mfield: inherit surrounding/global field
        str(ncopy), "1", "1",
        "1" if visible else "0",
        "1",
        "no", "no", "no",
    ]
    assert len(fields) == 18
    return " | ".join(fields)


def component(name, desc, position, rotation, typ, dimensions):
    return row(
        name, "root", desc, position, rotation,
        COLORS["component"], typ, dimensions, "Component", visible=False
    )


def operation_component(name, desc, left, op, right):
    return row(
        name, "root", desc, pos(0,0,0), rot(),
        COLORS["component"], f"Operation:@ {left} {op} {right}",
        "0", "Component", visible=False
    )


# =============================================================================
# MATERIALS
# =============================================================================

def epoxy_mass_fractions():
    masses = {el: EPOXY_ATOMS[el] * ATOMIC_MASS[el] for el in EPOXY_ATOMS}
    s = sum(masses.values())
    return {el: v/s for el, v in masses.items()}


def cf_mass_fractions():
    e = epoxy_mass_fractions()
    c = {
        "H": 0.30 * e["H"],
        "C": 0.70 + 0.30 * e["C"],
        "N": 0.30 * e["N"],
        "O": 0.30 * e["O"],
    }
    s = sum(c.values())
    return {el: v/s for el, v in c.items()}


def material_row(name, desc, density, components):
    comp = " ".join(f"G4_{el} {fmt(frac)}" for el, frac in components.items())
    fields = [
        name, desc, fmt(density), str(len(components)), comp,
        "none", "none", "none", "none", "none"
    ]
    assert len(fields) == 10
    return "|".join(fields)


# =============================================================================
# NATIVE Z POSITIONS
# =============================================================================

def north_endwheel_origin_z():
    # C++ line:
    # sEndWheelSNHolesZdist/2 - (sEndWStepHoleZpos+sEndWStepHoleZdist)
    return END_WHEEL_SN_HOLES_ZDIST/2.0 - (
        END_W_STEP_HOLE_ZPOS + END_W_STEP_HOLE_ZDIST
    )


def cyss_north_flange_origin_z():
    return (
        north_endwheel_origin_z()
        + END_WHEEL_N_LEN
        + CYSS_FLG_INT_THICK
    )


def cyss_cylinder_center_z():
    zpos = cyss_north_flange_origin_z()
    zpos -= (
        CYSS_FLG_INT_THICK
        - CYSS_CYL_LEN/2.0
        + CYSS_FLG_SSF_FLG_NSF
        - CYSS_FLG_SSF_CYLSF
    )
    return zpos


# =============================================================================
# NORTH END-WHEEL GEOMETRY
# =============================================================================

def make_north_endwheel(layer, helper_counter):
    rows = []
    dx, dy, dz = PLUSZ_SUPPORT_POSITION_MM

    z0 = north_endwheel_origin_z()

    # Base G4Polycone
    zplanes = [
        0.0,
        END_WHEEL_N_LEN - ENDW_N_THICK,
        END_WHEEL_N_LEN - ENDW_N_THICK,
        END_WHEEL_N_LEN,
    ]
    rmins = [
        ENDW_N_WALL_R[layer],
        ENDW_N_WALL_R[layer],
        ENDW_N_RMIN[layer],
        ENDW_N_RMIN[layer],
    ]
    rmaxs = [ENDW_N_RMAX[layer]] * 4

    base = f"cmp_ewn_L{layer}_base_{helper_counter:04d}"
    helper_counter += 1
    rows.append(component(
        base, f"L{layer} North EndWheel uncut basis",
        pos(0,0,0), rot(), "Polycone",
        polycone_dims(0,360,zplanes,rmins,rmaxs)
    ))

    subtractors = []

    if INCLUDE_ENDWHEEL_HOLES:
        # -------------------------------------------------------------
        # Holes through the vertical cylindrical wall.
        # -------------------------------------------------------------
        radial_mid = 0.5 * (ENDW_N_RMAX[layer] + ENDW_N_WALL_R[layer])
        dphi = 180.0 / ENDW_N_WALL_NHOLES[layer]
        phimin = dphi / 2.0

        for ihole in range(2 * ENDW_N_WALL_NHOLES[layer]):
            phi = phimin + ihole*dphi
            x = radial_mid * math.sin(math.radians(phi))
            y = radial_mid * math.cos(math.radians(phi))
            hname = f"cmp_ewn_L{layer}_wallhole_{ihole:02d}_{helper_counter:04d}"
            helper_counter += 1

            # C++: G4Tubs half-length = 1.5*(Rmax-WallR),
            # rotation Ra.set(-phi,90deg,0), center z=4 mm.
            rows.append(component(
                hname, f"L{layer} EndWheel wall hole {ihole}",
                pos(x,y,END_W_STEP_HOLE_ZPOS),
                rot(-phi,90.0,0.0),
                "Tube",
                tube_dims(
                    0.0,
                    ENDW_N_WALL_HOLE_D/2.0,
                    1.5*(ENDW_N_RMAX[layer]-ENDW_N_WALL_R[layer]),
                )
            ))
            subtractors.append(hname)

        # -------------------------------------------------------------
        # Holes through the end-wheel base plate.
        # C++ creates 11 azimuths and a diametrically-opposite copy.
        # -------------------------------------------------------------
        phi = 0.0
        plate_z = END_WHEEL_N_LEN - ENDW_N_THICK/2.0

        n_total = ENDW_N_BASE_N_BIG_HOLES + ENDW_N_BASE_N_SMALL_HOLES
        for ihole in range(n_total):
            phi += ENDW_N_BASE_HOLES_PHI_DEG
            x = ENDW_N_BASE_HOLES_R[layer]*math.cos(math.radians(phi))
            y = ENDW_N_BASE_HOLES_R[layer]*math.sin(math.radians(phi))

            # Reproduce the C++ small/big selection.
            is_small = (ihole < 2 or ihole == 5 or ihole > 8)
            diameter = (
                ENDW_N_BASE_SMALL_HOLE_D if is_small
                else ENDW_N_BASE_BIG_HOLE_D
            )

            for side, sx, sy in (("a",x,y),("b",-x,-y)):
                hname = (
                    f"cmp_ewn_L{layer}_basehole_{ihole:02d}{side}_"
                    f"{helper_counter:04d}"
                )
                helper_counter += 1
                rows.append(component(
                    hname, f"L{layer} EndWheel base hole {ihole}{side}",
                    pos(sx,sy,plate_z), rot(),
                    "Tube",
                    tube_dims(0.0,diameter/2.0,1.5*ENDW_N_THICK)
                ))
                subtractors.append(hname)

    # Chain all subtraction solids.  The last subtraction is made by the
    # actual physical-volume row so it acquires aluminum material.
    current = base

    if subtractors:
        for ih, hname in enumerate(subtractors[:-1]):
            opname = f"cmp_ewn_L{layer}_sub_{ih:03d}_{helper_counter:04d}"
            helper_counter += 1
            rows.append(operation_component(
                opname, f"L{layer} EndWheel subtraction chain",
                current, "-", hname
            ))
            current = opname

        final_type = f"Operation:@ {current} - {subtractors[-1]}"
        final_dims = "0"
    else:
        # No cuts requested: use a CopyOf-compatible physical type by creating
        # a directly identical Polycone row rather than trying to instantiate
        # a Component solid.
        final_type = "Polycone"
        final_dims = polycone_dims(0,360,zplanes,rmins,rmaxs)

    physical_name = f"mvtx_plusz_endwheel_L{layer}_basis"
    rows.append(row(
        physical_name, MOTHER,
        f"[EXACT_CPP] L{layer} North aluminum EndWheel basis"
        + (" with drilled holes" if INCLUDE_ENDWHEEL_HOLES else " (holes omitted)"),
        pos(dx,dy,dz+z0), rot(), COLORS["al"],
        final_type, final_dims, "MVTX_EW_Al", visible=VISIBLE
    ))

    # -----------------------------------------------------------------
    # Stave attachment step solid and copies
    # -----------------------------------------------------------------
    if INCLUDE_ENDWHEEL_STEPS:
        rmin = ENDW_N_WALL_R[layer]
        xlen = ENDW_N_STEP_XLEN[layer] - (0.01 if layer == 1 else 0.0)
        ydisp = ENDW_N_STEP_YDISPL[layer]

        xdispl = math.sqrt(rmin*rmin - ydisp*ydisp) - xlen
        ylen = math.sqrt(rmin*rmin - xdispl*xdispl) - ydisp

        step_zlen = END_WHEEL_N_LEN - ENDW_N_THICK
        box_hz = step_zlen/2.0 - 0.05

        xpos = xdispl + xlen/2.0
        ypos = ydisp + ylen/2.0
        rpos = math.hypot(xpos,ypos)

        phimin_deg = (
            90.0
            - math.degrees(math.acos(ydisp/rmin))
            - 5.0
        )
        dphi_deg = (
            90.0
            - math.degrees(math.asin(xdispl/rmin))
            - phimin_deg
            + 5.0
        )
        rmax = rmin + ylen

        boxname = f"cmp_ewn_L{layer}_step_box_{helper_counter:04d}"
        helper_counter += 1
        cutname = f"cmp_ewn_L{layer}_step_cut_{helper_counter:04d}"
        helper_counter += 1

        rows.append(component(
            boxname, f"L{layer} North EndWheel step box",
            pos(0,0,0), rot(), "Box",
            box_dims(xlen/2.0,ylen/2.0,box_hz)
        ))
        rows.append(component(
            cutname, f"L{layer} North EndWheel step curved cutter",
            pos(-xpos,-ypos,0), rot(), "Tube",
            tube_dims(
                rmin, rmax,
                1.5*step_zlen/2.0,
                phimin_deg,dphi_deg
            )
        ))

        step_type = f"Operation:@ {boxname} - {cutname}"
        dphi_local = math.degrees(math.atan2(ypos,xpos))
        phi0_rad = MVTXDAT[layer][4]
        nstaves = int(MVTXDAT[layer][5])

        template = None
        step_center_z = z0 + box_hz

        for j in range(nstaves):
            phi_rad = (
                phi0_rad
                + math.radians(j*ENDW_N_STEP_HOLE_PHI_DEG[layer])
                - ENDW_N_STEP_HOLE_TILT_RAD[layer]
            )
            x = rpos*math.cos(phi_rad)
            y = rpos*math.sin(phi_rad)
            zrot = dphi_local - math.degrees(phi_rad)

            name = f"mvtx_plusz_endwheel_L{layer}_step_S{j:02d}"
            if j == 0:
                rows.append(row(
                    name,MOTHER,
                    f"[EXACT_CPP] L{layer} North EndWheel stave step S{j:02d}",
                    pos(dx+x,dy+y,dz+step_center_z),
                    rot(0,0,zrot),
                    COLORS["al"],step_type,"0","MVTX_EW_Al",
                    visible=VISIBLE
                ))
                template = name
            else:
                rows.append(row(
                    name,MOTHER,
                    f"[EXACT_CPP] L{layer} North EndWheel stave step S{j:02d}",
                    pos(dx+x,dy+y,dz+step_center_z),
                    rot(0,0,zrot),
                    COLORS["al"],f"CopyOf {template}","0","MVTX_EW_Al",
                    visible=VISIBLE
                ))

    return rows, helper_counter


# =============================================================================
# CYSS +z GEOMETRY
# =============================================================================

def make_cyss_plusz():
    rows = []
    dx,dy,dz = PLUSZ_SUPPORT_POSITION_MM

    # -------------------------------------------------------------------------
    # North CYSS flange, exact C++ two-piece representation
    # -------------------------------------------------------------------------
    zf = cyss_north_flange_origin_z()

    z1 = [0.0,CYSS_FLG_INT_THICK,CYSS_FLG_INT_THICK,CYSS_FLG_S_RING_LEN]
    ri1 = [CYSS_FLG_RMIN,CYSS_FLG_RMIN,CYSS_FLG_S_RING_RMIN,CYSS_FLG_S_RING_RMIN]
    ro1 = [CYSS_FLG_L_RING_RMIN,CYSS_FLG_L_RING_RMIN,CYSS_FLG_S_RING_RMAX,CYSS_FLG_S_RING_RMAX]

    z2 = [0.0,CYSS_FLG_EXT_THICK,CYSS_FLG_EXT_THICK,CYSS_FLG_L_RING_LEN]
    ri2 = [CYSS_FLG_L_RING_RMIN]*4
    ro2 = [CYSS_FLG_RMAX,CYSS_FLG_RMAX,CYSS_FLG_L_RING_RMAX,CYSS_FLG_L_RING_RMAX]

    # C++ uses rotation Y = pi, so local +z points toward -global-z.
    rows.append(row(
        "mvtx_plusz_cyss_north_flange_inner",MOTHER,
        "[EXACT_CPP] North CYSS aluminum flange, inner polycone",
        pos(dx,dy,dz+zf),rot(0,180,0),COLORS["al"],
        "Polycone",polycone_dims(0,360,z1,ri1,ro1),
        "MVTX_EW_Al",visible=VISIBLE
    ))

    rows.append(row(
        "mvtx_plusz_cyss_north_flange_outer",MOTHER,
        "[EXACT_CPP] North CYSS aluminum flange, outer polycone",
        pos(dx,dy,dz+zf),rot(0,180,0),COLORS["al"],
        "Polycone",polycone_dims(0,360,z2,ri2,ro2),
        "MVTX_EW_Al",visible=VISIBLE
    ))

    # -------------------------------------------------------------------------
    # CYSS carbon-fiber cylinder.
    # Exact physical cylinder center is from the C++ arithmetic.
    # -------------------------------------------------------------------------
    zc = cyss_cylinder_center_z()
    zmin = zc - CYSS_CYL_LEN/2.0
    zmax = zc + CYSS_CYL_LEN/2.0

    if CLIP_CYSS_CYLINDER_AT_Z0:
        positive_min = max(0.0,zmin)
        positive_max = zmax
        if positive_max > positive_min:
            halfz = (positive_max-positive_min)/2.0
            center = 0.5*(positive_min+positive_max)
            rows.append(row(
                "mvtx_plusz_cyss_cf_cylinder",MOTHER,
                "[EXACT_CPP dimensions, +z-clipped] Carbon-fiber CYSS cylinder",
                pos(dx,dy,dz+center),rot(),COLORS["cf"],
                "Tube",tube_dims(CYSS_CYL_RMIN,CYSS_CYL_RMAX,halfz),
                "MVTX_CarbonFiber",visible=VISIBLE
            ))
    else:
        rows.append(row(
            "mvtx_plusz_cyss_cf_cylinder",MOTHER,
            "[EXACT_CPP] Full carbon-fiber CYSS cylinder",
            pos(dx,dy,dz+zc),rot(),COLORS["cf"],
            "Tube",tube_dims(CYSS_CYL_RMIN,CYSS_CYL_RMAX,CYSS_CYL_LEN/2.0),
            "MVTX_CarbonFiber",visible=VISIBLE
        ))

    return rows


# =============================================================================
# OUTPUT
# =============================================================================

def make_geometry():
    rows = []
    counter = 1
    for layer in range(3):
        r,counter = make_north_endwheel(layer,counter)
        rows += r
    rows += make_cyss_plusz()
    return rows


def make_materials():
    return [
        material_row(
            "MVTX_EW_Al",
            "sPHENIX MVTX support aluminum; PHG4MvtxSupport density 2.7 g/cm3",
            MVTX_AL_DENSITY,{"Al":1.0}
        ),
        material_row(
            "MVTX_CarbonFiber",
            "sPHENIX MVTX CF: 70 mass-% C + 30 mass-% H32C2N4O15 epoxy",
            MVTX_CF_DENSITY,cf_mass_fractions()
        ),
    ]


def validate(geometry,materials):
    for i,line in enumerate(geometry,1):
        if len(line.split("|")) != 18:
            raise RuntimeError(f"Geometry line {i} has wrong field count")
    for i,line in enumerate(materials,1):
        if len(line.split("|")) != 10:
            raise RuntimeError(f"Material line {i} has wrong field count")

    assert abs(north_endwheel_origin_z()-146.0) < 1e-12
    assert abs(cyss_north_flange_origin_z()-177.5) < 1e-12

    zc = cyss_cylinder_center_z()
    assert abs(zc - (-36.705)) < 1e-9

    zmax = zc+CYSS_CYL_LEN/2
    assert abs(zmax-176.18) < 1e-9


def notes():
    zc = cyss_cylinder_center_z()
    return f"""MVTX +Z SUPPORT / ENCLOSURE MODEL
================================

This is the geometry to use when the particle direction of interest is +z in
the native MVTX coordinate system.

It intentionally does NOT contain the meter-long Service Barrel generated by
PHG4MvtxSupport::CreateServiceBarrel(), because that object is located on -z.

Included
--------
North End Wheels, layers 0/1/2:
    base origin z = {north_endwheel_origin_z():.3f} mm
    bases extend to approximately +{north_endwheel_origin_z()+END_WHEEL_N_LEN:.3f} mm
    material = aluminum
    drilled holes included = {INCLUDE_ENDWHEEL_HOLES}
    stave attachment steps included = {INCLUDE_ENDWHEEL_STEPS}

North CYSS flange:
    placement origin z = {cyss_north_flange_origin_z():.3f} mm
    material = aluminum
    C++ rotation = 180 deg about y

CYSS carbon-fiber cylinder:
    physical full center z = {zc:.3f} mm
    physical full extent =
        {zc-CYSS_CYL_LEN/2:.3f} .. {zc+CYSS_CYL_LEN/2:.3f} mm
    radius =
        {CYSS_CYL_RMIN:.3f} .. {CYSS_CYL_RMAX:.3f} mm
    material = carbon fiber
    clipped at z=0 in output = {CLIP_CYSS_CYLINDER_AT_Z0}

Not included because it is on -z
--------------------------------
South End Wheels
Layer CF cones
CYSS cone
CYSS rib
CYSS south flange
meter-long Service Barrel
service-barrel/layer cable runs

Integration with the existing MVTX
----------------------------------
Default:
    MOTHER = "mvtx_envelope"

Therefore the coordinates are local to mvtx_envelope, and moving your existing
MVTX_ENVELOPE_POSITION_MM moves the sensors, stave supports, and this +z support
material together.

If run standalone, change:
    MOTHER = "root"
"""


def write_outputs(outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    geom=make_geometry()
    mats=make_materials()
    validate(geom,mats)

    (outdir/"mvtx_plusz_support__geometry_original.txt").write_text(
        "\n".join(geom)+"\n"
    )
    (outdir/"mvtx_plusz_support__materials_original.txt").write_text(
        "\n".join(mats)+"\n"
    )
    (outdir/"mvtx_plusz_support.gcard").write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<gcard>\n'
        '  <detector name="mvtx_plusz_support" factory="TEXT" variation="original"/>\n'
        '</gcard>\n'
    )
    (outdir/"MVTX_PLUSZ_SUPPORT_NOTES.txt").write_text(notes())

    physical = [
        l for l in geom
        if " | Component | " not in l and "|Component|" not in l
        and " | root | " not in l
    ]
    print("Created MVTX +z support model")
    print(f"  Mother:                         {MOTHER}")
    print(f"  North EndWheel z origin:        {north_endwheel_origin_z():.3f} mm")
    print(f"  North EndWheel max z:           {north_endwheel_origin_z()+END_WHEEL_N_LEN:.3f} mm")
    print(f"  North CYSS flange origin:       {cyss_north_flange_origin_z():.3f} mm")
    zc=cyss_cylinder_center_z()
    print(f"  CYSS cylinder physical extent:  {zc-CYSS_CYL_LEN/2:.3f} .. {zc+CYSS_CYL_LEN/2:.3f} mm")
    print(f"  CYSS clipped at z=0:            {CLIP_CYSS_CYLINDER_AT_Z0}")
    print(f"  EndWheel holes:                 {INCLUDE_ENDWHEEL_HOLES}")
    print(f"  EndWheel steps:                 {INCLUDE_ENDWHEEL_STEPS}")
    print(f"  Total geometry rows:            {len(geom)}")
    print("  Validation:                     OK")


def main():
    ap=argparse.ArgumentParser(
        description="Generate GEMC TEXT geometry for MVTX +z-side support material."
    )
    ap.add_argument("--outdir",type=Path,default=Path("mvtx_plusz_support_text"))
    args=ap.parse_args()
    write_outputs(args.outdir)


if __name__=="__main__":
    main()
