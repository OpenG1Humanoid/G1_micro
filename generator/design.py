"""OpenG1 micro: every dimension and every part.

World frame (zero pose, robot standing): X forward, Y left, Z up, mm, ground at z = 0.
Everything is designed for the LEFT side; right-side parts are mirror images (y -> -y).

Joint types
  * flange joint  : servo rides in the CHILD link, a standard HUB on the spline is bolted to
                    the parent; the hub flange is trapped in a groove in the child's housing
                    (hip pitch, shoulder pitch, waist yaw, hip yaw, shoulder yaw).
  * yoke joint    : servo sits in one link, the other link wraps it with two ears: a drive
                    ear screwed to a standard HORN PLATE on the spline, and an idler ear that
                    turns on a printed pin on the back of the servo housing
                    (hip/shoulder roll, knee, elbow, ankle pitch, ankle roll).

Moulding rule: every link is exactly two shells split by one flat plane. Each shell is a
straight-pull shape along the plane normal (every line along the pull hits it once), so each
shell can be cast in a simple two-part mould, and prints flat on its split face.
"""
import numpy as np
from sdf import *

# ------------------------------------------------------------------ servo envelope
# Servo frame S: origin at the top of the output spline, +z out of the spline,
# +x along the long side of the case (spline is SV_A from the near end), y across.
SV = dict(
    L=23.2,          # case length         (MG90S 22.8, MG92B 22.8, STS3032/SCS0009 23.2)
    W=12.2,          # case width          (12.0-12.2)
    A=6.0,           # spline centre to near end of case
    TOP=-6.3,        # case top face, below spline top
    BOT=-31.0,       # case bottom (MG92B 31; MG90S/STS3032 28.5 -> use the foam pad)
    EAR=(-15.5, -9.5),   # band the mounting ears can sit in (generous)
    EARLEN=32.6,     # tip to tip of the ears
)
CL = 0.15            # clearance per side around the case
WALL = 1.6           # 4 perimeters at 0.4 mm
LTOL = 0.6           # extra depth below the case (foam pad)
X0 = -SV["A"] - CL                     # -6.15
X1 = SV["L"] - SV["A"] + CL            # 17.35
HW = SV["W"] / 2 + CL                  # 6.25
ZB = SV["BOT"] - LTOL                  # -31.6
FLOOR = ZB - WALL                      # -33.2
TOPW = SV["TOP"] + WALL                # -4.7   (top wall of a yoke-type housing)

# horn: the 4-arm cross horn that comes with the servo, arms snipped to 13 mm tip-to-tip
HORN = dict(hub_r=3.9, reach=6.6, w_in=5.0, w_out=4.2, z0=-0.8, z1=1.4)
# hub (flange joints)
HUB = dict(fl_r=10.2, fl_z=(-6.0, -4.4), neck_r=8.0, top=2.4, bolt_r=5.6, recess=1.2)
# horn plate (yoke joints) and ears
HP = dict(r=9.5, z=(-0.6, 3.0), bolt_r=6.6)
EAR_D = (3.1, 5.5)        # drive ear, servo z
EAR_I = (-36.1, -33.5)    # idler ear
R_E = 10.5                # ear radius
BOSS_R = 3.0              # idler pin radius
SPAN_MID = (EAR_D[1] + EAR_I[0]) / 2   # -15.3: middle of a yoke along the servo axis

M2_CLR, M2_PILOT, M2_HEAD = 1.15, 0.8, 2.05   # radii
PIN_R, PIN_D = 0.95, 3.2   # 1.75 mm filament used as alignment dowels

# ------------------------------------------------------------------ kinematic layout (left)
KEEL = 4.5                                  # pelvis keel half thickness
HIP_Y0 = KEEL - HUB["recess"] + HUB["top"]  # 5.7  pitch servo spline top (y)
LEG_Y = HIP_Y0 + 20.8                       # 26.5 leg line
LEG_X = 12.5 + SPAN_MID                     # -2.8 leg line (x): roll spline top at x=12.5
ROLL_SPL_X = 12.5

Z_AR = 14.5                    # ankle roll axis
Z_A = Z_AR + 27.0              # ankle pitch axis
SHIN = 60.0
Z_K = Z_A + SHIN               # knee
# thigh stack, from knee up
Z_YS = Z_K + 52.05             # yaw servo spline top
Z_BF = Z_YS + HUB["top"] - HUB["recess"]   # bottom face of the roll-yoke spine (hub top sits 1.2 up in a recess)
Z_R = Z_BF + 17.1              # roll axis (spine below the roll servo is 7.5 mm deep)
Z_HP = Z_R + 26.6              # hip pitch axis (roll servo stands upright under the pitch servo)
Z_PT = Z_HP + 24.5             # pelvis top face
Z_SP = Z_HP + 124.0            # shoulder pitch axis
TORSO_Y = 38.0                 # torso half width (shoulder hub face)
SH_D = np.array([0.0, TORSO_Y - KEEL, Z_SP - Z_HP])   # hip -> shoulder offset of the 3-axis module
ARM_Y = LEG_Y + SH_D[1]
FOREARM = 56.0

# servo frames (left side, zero pose)
def F_(o, x, z):
    return Frame(o, x=x, z=z)

S_HP = F_((0, HIP_Y0, Z_HP), (0, 0, 1), (0, -1, 0))            # hip pitch, spline -> -y
S_HR = F_((ROLL_SPL_X, LEG_Y, Z_R), (0, 0, 1), (1, 0, 0))        # hip roll, spline -> +x, body up
S_HY = F_((LEG_X, LEG_Y, Z_YS), (0, -1, 0), (0, 0, 1))           # hip yaw, spline -> +z
KSPL = LEG_Y - SPAN_MID                                          # 41.8
S_KN = F_((LEG_X, KSPL, Z_K), (0, 0, 1), (0, 1, 0))              # knee, spline -> +y
S_AP = F_((LEG_X, KSPL, Z_A), (0, 0, 1), (0, 1, 0))              # ankle pitch
S_AR = F_((ROLL_SPL_X, LEG_Y, Z_AR), (0, -1, 0), (1, 0, 0))      # ankle roll (in foot)
S_WY = F_((0, 0, Z_PT - HUB["recess"] + HUB["top"]), (0, 1, 0), (0, 0, -1))  # waist, spline -> -z

def sh(fr, d=SH_D):
    return Frame(fr.o + d, x=fr.R[:, 0], z=fr.R[:, 2])

S_SP, S_SR, S_SY, S_EL = sh(S_HP), sh(S_HR), sh(S_HY), sh(S_KN)

# joints: name -> (parent, child, axis point, axis dir, lower, upper [deg])
AX = {
    "hip_pitch":   ((0, HIP_Y0, Z_HP), (0, 1, 0), -130, 60),
    "hip_roll":    ((0, LEG_Y, Z_R), (1, 0, 0), -30, 45),
    "hip_yaw":     ((LEG_X, LEG_Y, 0), (0, 0, 1), -60, 60),
    "knee":        ((LEG_X, 0, Z_K), (0, 1, 0), -5, 130),
    "ankle_pitch": ((LEG_X, 0, Z_A), (0, 1, 0), -50, 50),
    "ankle_roll":  ((0, LEG_Y, Z_AR), (1, 0, 0), -22, 22),
    "waist_yaw":   ((0, 0, 0), (0, 0, 1), -90, 90),
    "shoulder_pitch": ((0, HIP_Y0 + SH_D[1], Z_SP), (0, 1, 0), -180, 90),
    "shoulder_roll":  ((0, ARM_Y, Z_R + SH_D[2]), (1, 0, 0), -10, 100),
    "shoulder_yaw":   ((LEG_X, ARM_Y, 0), (0, 0, 1), -90, 90),
    "elbow":          ((LEG_X, 0, Z_K + SH_D[2]), (0, 1, 0), -135, 5),
}
SWEEP = {  # ranges used for clearance cuts (shared parts take the union of both uses)
    "roll": (-33, 102), "knee": (-6, 134), "elbow": (-137, 6), "ankle_pitch": (-52, 52),
    "ankle_roll": (-23, 23),
}


# ------------------------------------------------------------------ servo features
def horn_cut(S, z0=None, z1=None):
    """Pocket for the trimmed cross horn (arms at 45 deg in servo frame)."""
    z0 = HORN["z0"] if z0 is None else z0
    z1 = HORN["z1"] if z1 is None else z1
    parts = [cyl(S, HORN["hub_r"], z0 - 2.2, z1)]
    for k in range(4):
        a = np.radians(45 + 90 * k)
        Fa = Frame(S.o, x=S.wdir((np.cos(a), np.sin(a), 0)), z=S.R[:, 2])
        r, wi, wo = HORN["reach"], HORN["w_in"] / 2, HORN["w_out"] / 2
        parts.append(prism(Fa, [(0, -wi), (r - wo, -wo), (r - wo, wo), (0, wi)], z0, z1))
        parts.append(cyl(Fa, wo, z0, z1, c=(r - wo, 0)))
    return U(*parts)


def pocket(S, kind):
    """Space the servo needs (cut)."""
    e = (SV["EARLEN"] - SV["L"]) / 2
    cuts = [
        boxr(S, (X0, -HW, ZB), (X1, HW, SV["TOP"] + 0.01)),
        boxr(S, (X0 - e - 0.3, -HW, SV["EAR"][0] - 0.2), (X1 + e + 0.3, HW, SV["EAR"][1] + 0.2)),
        cyl(S, 3.6, ZB - 0.8, ZB + 0.1),                            # rear shaft (double-shaft servos)
        boxr(S, (X1 - 0.6, -2.2, ZB), (X1 + WALL + 1.5, 2.2, SV["EAR"][0])),   # cable slot
    ]
    if kind == "yoke":
        cuts += [cyl(S, 6.45, SV["TOP"] - 0.1, 8.0),
                 boxr(S, (0, -4.0, SV["TOP"] - 0.1), (8.5, 4.0, TOPW + 0.6))]  # gear bump
    else:  # flange
        cuts += [cyl(S, HUB["fl_r"] + 0.4, SV["TOP"] - 0.01, -4.1),   # hub flange groove
                 cyl(S, HUB["neck_r"] + 0.2, -4.2, 8.0)]                               # hub neck hole
    return U(*cuts)


def housing(S, kind, r=1.4):
    """Solid around the servo."""
    top = TOPW if kind == "yoke" else -1.8
    s = [rboxr(S, (X0 - WALL, -HW - WALL, FLOOR), (X1 + WALL, HW + WALL, top), r)]
    if kind == "flange":
        s.append(rcyl(S, 12.0, SV["TOP"] - 1.6, -1.8, 0.8))
    else:
        s.append(cyl(S, BOSS_R, EAR_I[0] + 0.1, FLOOR + 0.2))         # idler pin
    return U(*s)


def ear(S, which, toward, reach, R=R_E, extra_w=0.0):
    """Tombstone ear: disc of radius R round the servo axis plus a strap of width 2R
    running in direction `toward` (unit vector in servo x-y) out to `reach`."""
    z0, z1 = EAR_D if which == "drive" else EAR_I
    tx, ty = toward
    Fa = Frame(S.o, x=S.wdir((tx, ty, 0)), z=S.R[:, 2])
    w = R + extra_w
    return U(rcyl(S, R, z0, z1, 0.6), rboxr(Fa, (0, -w, z0), (reach, w, z1), 0.6))


def ear_cuts(S, which, bolt_ang=0.0):
    if which == "idler":
        return cyl(S, BOSS_R + 0.18, EAR_I[0] - 1, EAR_I[1] + 1)
    c = [cyl(S, 2.6, EAR_D[0] - 1, EAR_D[1] + 1)]                       # horn screw access
    for k in (0, 1):
        a = np.radians(bolt_ang + 180 * k)
        p = (HP["bolt_r"] * np.cos(a), HP["bolt_r"] * np.sin(a))
        c += [cyl(S, M2_CLR, EAR_D[0] - 1, EAR_D[1] + 1, c=p),
              cyl(S, M2_HEAD, EAR_D[1] - 0.9, EAR_D[1] + 3, c=p)]
    return U(*c)


def hub_seat(S, bolt_angles, head_from, head_len=20.0):
    """Cuts in a PARENT for a hub: recess for the neck + bolt holes from the inside.
    head_from: servo-z where the inside surface of the parent wall is (screw heads)."""
    top = HUB["top"]
    c = [cyl(S, HUB["neck_r"] + 0.2, top - HUB["recess"] - 8.0, top)]
    for a in bolt_angles:
        a = np.radians(a); p = (HUB["bolt_r"] * np.cos(a), HUB["bolt_r"] * np.sin(a))
        c += [cyl(S, M2_CLR, top, head_from + 0.5, c=p),
              cyl(S, M2_HEAD, head_from - 1.5, head_from - 1.5 + head_len, c=p)]
    return U(*c)


# ------------------------------------------------------------------ standard parts
def part_hub():
    S = Frame()
    solid = U(rcyl(S, HUB["fl_r"], *HUB["fl_z"], 0.4), rcyl(S, HUB["neck_r"], HUB["fl_z"][0] + 0.5, HUB["top"], 0.5))
    cuts = [cyl(S, 6.8, -8, HORN["z0"]), horn_cut(S, HORN["z0"] - 0.01, HORN["z1"]),
            cyl(S, 1.3, -8, 8), cyl(S, 2.45, HORN["z1"], 8)]
    for a in (0, 90, 180, 270):
        a = np.radians(a)
        cuts.append(cyl(S, M2_PILOT, HORN["z0"] - 0.5, 8, c=(HUB["bolt_r"] * np.cos(a), HUB["bolt_r"] * np.sin(a))))
    return dict(name="hub", solid=solid, cuts=cuts, bbox=((-11, -11, -7), (11, 11, 3.5)), plane=None,
                print_dir=(0, 0, 1), frame=S)


def part_hornplate():
    S = Frame()
    solid = rcyl(S, HP["r"], HP["z"][0], HP["z"][1], 0.6)
    cuts = [horn_cut(S, HP["z"][0] - 0.3, HORN["z1"]), cyl(S, 1.3, -8, 8), cyl(S, 2.45, HORN["z1"], 8)]
    for a in (0, 90, 180, 270):
        a = np.radians(a)
        cuts.append(cyl(S, M2_PILOT, -0.2, 8, c=(HP["bolt_r"] * np.cos(a), HP["bolt_r"] * np.sin(a))))
    return dict(name="horn_plate", solid=solid, cuts=cuts, bbox=((-10, -10, -1.5), (10, 10, 4)), plane=None,
                print_dir=(0, 0, -1), frame=S)


# ------------------------------------------------------------------ helpers for links
def screws_along(axis, plane_c, head_side, pts, head_depth, length, extra_cuts=None):
    """M2 self-tappers along `axis` (0/1/2) crossing the split plane.
    pts: list of (u, v) in the two other axes (cyclic order axis+1, axis+2).
    head_depth: distance of the screw head seat from the plane on the head side."""
    out = []
    for (u, v) in pts:
        def mk(u=u, v=v):
            ax = np.eye(3)[axis]; a1 = np.eye(3)[(axis + 1) % 3]; a2 = np.eye(3)[(axis + 2) % 3]
            o = ax * plane_c + a1 * u + a2 * v
            Fr = Frame(o, x=a1, z=ax * head_side)   # +z points to the head side
            return U(cyl(Fr, M2_CLR, -0.01, head_depth + 0.01),
                     cyl(Fr, M2_HEAD, head_depth, 60),
                     cyl(Fr, M2_PILOT, -(length - head_depth) - 0.5, 0.01))
        out.append(mk())
    return out


def pins_along(axis, plane_c, pts):
    out = []
    for (u, v) in pts:
        ax = np.eye(3)[axis]; a1 = np.eye(3)[(axis + 1) % 3]; a2 = np.eye(3)[(axis + 2) % 3]
        o = ax * plane_c + a1 * u + a2 * v
        Fr = Frame(o, x=a1, z=ax)
        out.append(cyl(Fr, PIN_R, -PIN_D, PIN_D))
    return out


def plane_pt(axis, c, u, v):
    p = np.zeros(3); p[axis] = c; p[(axis + 1) % 3] = u; p[(axis + 2) % 3] = v
    return p


# ================================================================== LINKS
def part_pelvis():
    W = WORLD
    keel = rboxr(W, (-13, -KEEL, Z_HP - 14), (13, KEEL, Z_PT), 2.5)
    plate = rboxr(W, (-15, -11, Z_HP + 21.2), (15, 11, Z_PT), 3.0)
    belt = rboxr(W, (-15, -KEEL, Z_HP - 14), (15, KEEL, Z_HP + 2), 3.0)
    solid = SU(SU(keel, plate, 3.0), belt, 2.0)
    cuts = []
    # hollow, open to the split plane y=0, leaving bosses round the hub bolts
    cav = boxr(W, (-10.5, -1.6, Z_HP - 11), (10.5, 1.6, Z_PT - 2.6))
    keep = [cyl(Frame((0, 0, Z_HP + s * HUB["bolt_r"]), x=(1, 0, 0), z=(0, 1, 0)), 3.0, -9, 9) for s in (1, -1)]
    cuts.append(D(cav, *keep))
    for side in (1, -1):
        S = S_HP if side == 1 else Frame(S_HP.o * np.array([1, -1, 1]), x=(0, 0, 1), z=(0, 1, 0))
        # recess for hub neck + 2 bolts (0/180 deg = along z); heads inside the hollow
        cuts.append(hub_seat(S, (0, 180), head_from=HIP_Y0 - 1.5 + 1.5, head_len=1.5 + 0.01))
    # waist hub on top: recess + 2 bolts on the plane (x = +-5.6), heads under the plate
    cuts.append(hub_seat(S_WY, (90, 270), head_from=HUB["top"] + HUB["recess"] + 2.4))
    for s in (1, -1):   # screwdriver access from below to the waist bolts
        cuts.append(cyl(Frame((s * HUB["bolt_r"], 0, 0), z=(0, 0, 1), x=(1, 0, 0)), 2.1, Z_HP - 20, Z_PT - 3.4))
    # halves joined by screws along y (heads on the right half outside face)
    pts = [(Z_HP - 9.5, -10.0), (Z_HP - 9.5, 10.0), (Z_PT - 6.5, 9.5), (Z_PT - 6.5, -9.5)]   # (z, x)
    cuts += screws_along(1, 0.0, -1, pts, head_depth=2.6, length=8)
    cuts += pins_along(1, 0.0, [(Z_HP + 8, 8.5), (Z_HP + 8, -8.5)])
    return dict(name="pelvis", solid=solid, cuts=cuts, plane=(1, 0.0),
                bbox=((-16, -12, Z_HP - 22), (16, 12, Z_PT + 1)), mirror=False)


def part_hip_pitch():
    """Hip / shoulder pitch link: pitch servo (flange joint to pelvis/torso) + roll servo."""
    hp = housing(S_HP, "flange"); hr = housing(S_HR, "yoke")
    # screw blocks: above the pitch servo (inner part) and at the outer end
    blk1 = rboxr(WORLD, (-7.85, 11.5, Z_HP + 15), (7.85, 20.5, Z_HP + 24.5), 2.5)
    blk2 = rboxr(WORLD, (-7.85, 7.5, Z_R + 2.0), (7.85, 19.0, Z_HP - 7.0), 2.5)
    solid = SU(SU(hp, hr, 4.0), U(blk1, blk2), 2.0)
    cuts = [pocket(S_HP, "flange"), pocket(S_HR, "yoke")]
    # pitch cable leaves at the top end of the pitch servo, roll cable at the inner end
    cuts += screws_along(0, 0.0, +1, [(15.8, Z_HP + 21.0), (12.6, Z_R + 9.0)], head_depth=4.0, length=10)
    cuts += pins_along(0, 0.0, [(17.8, Z_HP + 17.2), (15.6, Z_R + 15.0)])
    return dict(name="hip_pitch_link", solid=solid, cuts=cuts, plane=(0, 0.0),
                bbox=((-24, 4, Z_R - 9), (13, 45, Z_HP + 25.5)), mirror=True)


def part_hip_roll():
    """Roll yoke round the roll servo + seat for the yaw hub underneath."""
    S = S_HR
    down = (-1, 0)  # servo -x = world -z
    reach = -X0 + WALL + 0.8 + 5.0
    e1 = ear(S, "drive", down, reach); e2 = ear(S, "idler", down, reach)
    spine = rboxr(WORLD, (EAR_I[0] + ROLL_SPL_X, LEG_Y - 12.2, Z_BF), (EAR_D[1] + ROLL_SPL_X, LEG_Y + 12.2, Z_R - 8.55), 2.0)
    solid = U(e1, e2, spine)
    cuts = [ear_cuts(S, "drive", 0), ear_cuts(S, "idler")]
    cuts.append(hub_seat(S_HY, (0, 180), head_from=HUB["top"] + 2.0 + 1.5))
    # long M2x25 through the spine from the front ear face
    pts = [(LEG_Y - 9.7, Z_BF + 3.4), (LEG_Y + 9.7, Z_BF + 3.4)]
    cuts += screws_along(0, LEG_X, +1, pts, head_depth=EAR_D[1] + ROLL_SPL_X - LEG_X - 1.0, length=25)
    cuts += pins_along(0, LEG_X, [(LEG_Y - 10.6, Z_BF + 6.2), (LEG_Y + 10.6, Z_BF + 6.2)])
    return dict(name="hip_roll_link", solid=solid, cuts=cuts, plane=(0, LEG_X),
                bbox=((-25, 12, Z_BF - 1), (19.5, 41, Z_R + 11)), mirror=True,
                sweeps=[("hip_pitch_link", AX["hip_roll"][0], AX["hip_roll"][1], SWEEP["roll"], 0.6, "", 30.0)])


def part_thigh():
    """Thigh = upper arm: yaw servo on top (flange joint), knee/elbow servo at the bottom."""
    hy = housing(S_HY, "flange"); hk = housing(S_KN, "yoke")
    zlo, zhi = Z_K - 7.75, Z_YS - 1.8
    shellb = rboxr(WORLD, (LEG_X - 9.4, 9.3, zlo), (LEG_X + 9.4, KSPL + TOPW, zhi - 3.0), 5.5)
    rib = rboxr(WORLD, (LEG_X - 9.4, 30, Z_K + 21), (LEG_X + 9.4, 39.8, zhi - 6), 3.0)
    solid = SU(U(hy, hk), U(shellb, rib), 3.0)
    # keep the knee ears' swing free (shin ears turn round the knee axis)
    ears_free = [cyl(S_KN, R_E + 0.6, EAR_I[0] - 2, EAR_I[1] + 0.3), cyl(S_KN, R_E + 0.6, HP["z"][0] - 0.3, 12)]
    cuts = [pocket(S_HY, "flange"), pocket(S_KN, "yoke")] + ears_free
    pts = [(36.6, Z_YS - 12), (36.6, Z_YS - 25), (8.6, Z_K + 15.5)]
    cuts += screws_along(0, LEG_X, +1, pts, head_depth=4.5, length=10)
    cuts += pins_along(0, LEG_X, [(36.8, Z_YS - 18.5), (16, Z_K + 20.6)])
    return dict(name="thigh", solid=solid, cuts=cuts, plane=(0, LEG_X),
                bbox=((LEG_X - 12.5, 4, Z_K - 9), (LEG_X + 12.5, 41, Z_YS)), mirror=True)


def part_shin():
    S = S_KN
    down = (-1, 0)   # servo -x = world -z
    reach = 27.0
    e1 = ear(S, "drive", down, reach); e2 = ear(S, "idler", down, reach)
    ha = housing(S_AP, "yoke")
    ztop = Z_K - 7.75 - 1.0
    body = rboxr(WORLD, (LEG_X - 9.0, 8.1, Z_A + 8), (LEG_X + 9.6, 44.2, ztop), 4.5)
    bridge = rboxr(WORLD, (LEG_X - 9.0, KSPL + EAR_I[0], Z_K - 31), (LEG_X + 9.6, KSPL + EAR_D[1], ztop), 2.5)
    solid = SU(SU(body, ha, 3.0), U(e1, e2, bridge), 2.0)
    # hollow middle (open to the split plane)
    cav = rboxr(WORLD, (LEG_X - 7.4, 10.0, Z_A + 21), (LEG_X + 8.0, 41.8, Z_K - 24), 2.0)
    keep = [cyl(Frame((LEG_X, y, z), z=(1, 0, 0), x=(0, 1, 0)), 2.9, -12, 12) for (y, z) in ((14, Z_A + 31), (38, Z_A + 31), (26, Z_K - 27.5))]
    cuts = [pocket(S_AP, "yoke"), ear_cuts(S, "drive", 0), ear_cuts(S, "idler"), D(cav, *keep)]
    # ankle link ears swing round the ankle axis
    cuts += [D(cyl(S_AP, R_E + 0.7, EAR_I[0] - 2, EAR_I[1] + 0.3), boxr(S_AP, (-30, -30, EAR_I[1] + 0.3), (30, 30, 30))),
             cyl(S_AP, R_E + 0.7, HP["z"][0] - 0.3, 12)]
    pts = [(14, Z_A + 31), (38, Z_A + 31), (26, Z_K - 27.5)]
    cuts += screws_along(0, LEG_X, +1, pts, head_depth=5.0, length=10)
    cuts += pins_along(0, LEG_X, [(7.0, Z_K - 26.5), (45.6, Z_K - 26.5)])
    return dict(name="shin", solid=solid, cuts=cuts, plane=(0, LEG_X),
                bbox=((LEG_X - 12, 3, Z_A - 9), (LEG_X + 12, 49, Z_K + 12)), mirror=True,
                sweeps=[("thigh", AX["knee"][0], AX["knee"][1], SWEEP["knee"], 0.6, "", 21.5)])


def part_ankle():
    """Cross-shaped double yoke: pitch ears up round the shin's ankle servo, roll ears down
    round the foot's roll servo."""
    up = (1, 0)                 # S_AP: servo +x = world +z ... ears hang DOWN from the axis
    pdown = (-1, 0)
    z_sp0 = Z_AR + HW + WALL + 0.8        # spine bottom
    z_sp1 = Z_A - 7.75 - 0.8               # spine top
    reach_p = Z_A - z_sp0                  # pitch ears reach down to the spine bottom
    reach_r = z_sp1 - Z_AR
    e1 = ear(S_AP, "drive", pdown, reach_p); e2 = ear(S_AP, "idler", pdown, reach_p)
    e3 = ear(S_AR, "drive", (0, -1), reach_r); e4 = ear(S_AR, "idler", (0, -1), reach_r)   # S_AR +y = world -z
    armx = rboxr(WORLD, (EAR_I[0] + ROLL_SPL_X, LEG_Y - R_E, z_sp0), (EAR_D[1] + ROLL_SPL_X, LEG_Y + R_E, z_sp1), 1.5)
    army = rboxr(WORLD, (LEG_X - R_E, KSPL + EAR_I[0], z_sp0), (LEG_X + R_E, KSPL + EAR_D[1], z_sp1), 1.5)
    solid = U(e1, e2, e3, e4, SU(armx, army, 3.0))
    cuts = [ear_cuts(S_AP, "drive", 0), ear_cuts(S_AP, "idler"), ear_cuts(S_AR, "drive", 0), ear_cuts(S_AR, "idler")]
    zs = (z_sp0 + z_sp1) / 2
    cuts += screws_along(0, LEG_X, +1, [(LEG_Y - 4.5, zs), (LEG_Y + 4.5, zs)],
                         head_depth=EAR_D[1] + ROLL_SPL_X - LEG_X - 1.0, length=25)
    cuts += pins_along(0, LEG_X, [(KSPL + EAR_I[0] + 4.5, zs), (KSPL + EAR_D[1] - 4.5, zs)])
    return dict(name="ankle_link", solid=solid, cuts=cuts, plane=(0, LEG_X),
                bbox=((-25, 3, Z_AR - 12), (20, 49, Z_A + 12)), mirror=True,
                sweeps=[("shin", AX["ankle_pitch"][0], AX["ankle_pitch"][1], SWEEP["ankle_pitch"], 0.6, "", 22.0),
                        ("foot", AX["ankle_roll"][0], AX["ankle_roll"][1], SWEEP["ankle_roll"], 0.6, "parent_moves", 24.0)])


def part_foot():
    hr = housing(S_AR, "yoke")
    x0, x1 = LEG_X - 34, LEG_X + 58
    y0, y1 = LEG_Y - 20.0, LEG_Y + 19.5
    sole = rboxr(WORLD, (x0, y0, 0), (x1, y1, 5.0), 2.5)
    top = rboxr(WORLD, (x0 + 6, y0 + 0.2, -20), (x1 - 14, y1 - 2, Z_AR + 4), 3.5)
    solid = I(SU(U(sole, top), hr, 3.0), halfspace(WORLD, 2, 0.0))
    pts = [(LEG_Y + 15.0, LEG_X - 12), (LEG_Y + 15.0, LEG_X + 6), (LEG_Y - 15.8, LEG_X + 30), (LEG_Y + 15.0, LEG_X + 30)]   # (y, x)
    # hollow inside the raised shell (open to the split plane from both sides)
    hx0, hx1 = ROLL_SPL_X + EAR_I[0] - 1.5, ROLL_SPL_X + TOPW + 1.0
    keep = [boxr(WORLD, (hx0, LEG_Y - X1 - WALL - 1.0, -5), (hx1, LEG_Y - X0 + WALL + 1.0, 40)),
            boxr(WORLD, (ROLL_SPL_X - 1.27, LEG_Y - HP["r"] - 1.53, -5), (ROLL_SPL_X + EAR_D[1] + 2, LEG_Y + HP["r"] + 1.5, 40))]
    keep += [cyl(Frame((p[1], p[0], 0), z=(0, 0, 1), x=(1, 0, 0)), 3.0, -5, 40) for p in pts]
    keep += [cyl(Frame((x, y, 0), z=(0, 0, 1), x=(1, 0, 0)), 2.6, -5, 40) for (x, y) in ((LEG_X + 40, LEG_Y), (LEG_X - 28, LEG_Y + 5))]
    cav = D(I(offset(top, -1.6), halfspace(WORLD, 2, 1.63)), *keep)
    cuts = [pocket(S_AR, "yoke"), cav]
    # slots for the ankle link's roll ears
    for zz in (EAR_D, EAR_I):
        cuts.append(cyl(S_AR, R_E + 0.8, zz[0] - 0.4, zz[1] + 0.4))
        cuts.append(boxr(S_AR, (-R_E - 0.8, -40, zz[0] - 0.4), (R_E + 0.8, 0, zz[1] + 0.4)))
    cuts.append(cyl(S_AR, HP["r"] + 0.6, HP["z"][0] - 0.3, EAR_D[0]))
    # sole grip-pad recesses (1 mm), front and heel
    cuts.append(rboxr(WORLD, (x1 - 30, y0 + 3, -1), (x1 - 4, y1 - 3, 1.0), 2.0))
    cuts.append(rboxr(WORLD, (x0 + 4, y0 + 3, -1), (x0 + 24, y1 - 3, 1.0), 2.0))
    cuts += screws_along(2, Z_AR, +1, [(p[1], p[0]) for p in pts], head_depth=3.5, length=10)
    cuts += pins_along(2, Z_AR, [(LEG_X + 40, LEG_Y), (LEG_X - 28, LEG_Y + 5)])
    return dict(name="foot", solid=solid, cuts=cuts, plane=(2, Z_AR),
                bbox=((x0 - 1, y0 - 1, -1), (x1 + 1, y1 + 1, Z_AR + 9)), mirror=True)


def part_forearm():
    S = S_EL
    down = (-1, 0)
    reach = 27.0
    e1 = ear(S, "drive", down, reach); e2 = ear(S, "idler", down, reach)
    ze = S.o[2]; yc = ARM_Y; xc = LEG_X
    ztop = ze - 7.75 - 1.0
    bridge = rboxr(WORLD, (xc - 9.0, S.o[1] + EAR_I[0], ze - 31), (xc + 9.0, S.o[1] + EAR_D[1], ze - 20), 2.5)
    arm = rboxr(WORLD, (xc - 8.5, yc - 10, ze - 38), (xc + 8.5, yc + 10, ztop), 5.0)
    hand = SU(ellipsoid(WORLD, (9.0, 6.5, 11.5), (xc + 1.0, yc - 0.5, ze - 44.5)),
              ellipsoid(WORLD, (3.4, 3.0, 6.5), (xc + 7.2, yc - 5.0, ze - 39)), 2.0)
    solid = SU(SU(bridge, arm, 3.0), hand, 3.0)
    solid = U(solid, e1, e2)
    cav = rboxr(WORLD, (xc - 6.5, yc - 7.5, ze - 34), (xc + 6.5, yc + 7.5, ze - 21), 2.0)
    keep = [cyl(Frame((xc, yc, ze - 28), z=(1, 0, 0), x=(0, 1, 0)), 2.9, -12, 12)]
    cuts = [ear_cuts(S, "drive", 0), ear_cuts(S, "idler"), D(cav, *keep)]
    cuts += screws_along(0, xc, +1, [(yc, ze - 28), (yc - 1, ze - 46)], head_depth=4.5, length=10)
    cuts += pins_along(0, xc, [(S.o[1] + EAR_I[0] + 1.6, ze - 26), (S.o[1] + EAR_D[1] - 1.6, ze - 26)])
    return dict(name="forearm", solid=solid, cuts=cuts, plane=(0, xc),
                bbox=((xc - 12, yc - 23, ze - 58), (xc + 13, yc + 23, ze + 12)), mirror=True,
                sweeps=[("thigh@arm", AX["elbow"][0], AX["elbow"][1], SWEEP["elbow"], 0.6, "", 21.5)])


def part_torso():
    W = WORLD
    hw = housing(S_WY, "flange")
    zb, zt = Z_PT + 3.0, Z_SP + 14.0
    TY = TORSO_Y
    chest = rboxr(W, (-17.5, -TY, Z_SP - 24), (19.0, TY, zt), 7.0)
    belly = rboxr(W, (-16.5, -TY + 6, Z_PT + 22), (17.5, TY - 6, Z_SP - 10), 7.0)
    waist = rboxr(W, (-13.5, -20.5, zb), (14.5, 20.5, Z_PT + 30), 5.0)
    outer = SU(SU(chest, belly, 9.0), SU(waist, hw, 2.0), 7.0)
    # local flat faces for the shoulder hubs
    pads = [rboxr(W, (-11.5, s * TY - (5.0 if s > 0 else -5.0) * 0 - 5.0, Z_SP - 11.5), (11.5, s * TY + 0.0, Z_SP + 11.5), 3.0)
            if s > 0 else rboxr(W, (-11.5, -TY, Z_SP - 11.5), (11.5, -TY + 5.0, Z_SP + 11.5), 3.0) for s in (1, -1)]
    solid = U(outer, *pads)
    inner = I(offset(outer, -1.3), halfspace(W, 2, Z_PT + 35.5))
    keep = [rboxr(W, (-30, s * TY - (4.4 if s > 0 else -0.01), Z_SP - 10.5), (30, s * TY + (0.01 if s > 0 else 4.4), Z_SP + 10.5), 1.0) for s in (1, -1)]
    bosses = [cyl(Frame((0, yy, zz), z=(1, 0, 0), x=(0, 1, 0)), 3.0, -19, 20)
              for (yy, zz) in ((-17.0, Z_PT + 39), (17.0, Z_PT + 39), (-27.0, zt - 7.0), (27.0, zt - 7.0))]
    # 4 M2 bosses on the back wall for the controller board (40 x 20 mm hole pattern)
    board = [cyl(Frame((-12, yy, zz), z=(1, 0, 0), x=(0, 1, 0)), 2.6, -8, 0)
             for yy in (-20, 20) for zz in (Z_SP - 32, Z_SP - 12)]
    solid = U(solid, I(U(*bosses, *board), offset(inner, 0.5)))
    cav = D(inner, *keep, *bosses, *board)
    cuts = [pocket(S_WY, "flange"), cav]
    cuts += [cyl(Frame((-12, yy, zz), z=(1, 0, 0), x=(0, 1, 0)), M2_PILOT, -6, 1) for yy in (-20, 20) for zz in (Z_SP - 32, Z_SP - 12)]
    for s in (1, -1):
        S = S_SP if s == 1 else Frame(S_SP.o * np.array([1, -1, 1]), x=(0, 0, 1), z=(0, 1, 0))
        cuts.append(hub_seat(S, (0, 180), head_from=HUB["top"] + HUB["recess"] + 2.6 - 0.01))
    # head mount: 2 M2 from inside the top, on the split plane
    for yy in (-6, 6):
        Fr = Frame((0, yy, zt), z=(0, 0, 1), x=(1, 0, 0))
        cuts += [cyl(Fr, M2_CLR, -4, 1), cyl(Fr, M2_HEAD, -30, -2.4)]
    cuts.append(cyl(Frame((0, 0, zt), z=(0, 0, 1), x=(1, 0, 0)), 5.0, -0.8, 2))   # neck spigot recess
    # back: USB / switch openings, front: chest panel line (cosmetic, 0.6 mm groove)
    cuts.append(rboxr(W, (-25, -9, Z_SP - 28), (-14, 9, Z_SP - 20), 1.5))
    cuts.append(D(rboxr(W, (17.5, -TY + 7, Z_SP - 22), (30, TY - 7, zt - 4), 4.0),
                  rboxr(W, (17.5, -TY + 7.7, Z_SP - 21.3), (30, TY - 7.7, zt - 4.7), 3.3),
                  halfspace(W, 0, 18.4, keep_above=False)))
    pts = [(-17.0, Z_PT + 39), (17.0, Z_PT + 39), (-27.0, zt - 7.0), (27.0, zt - 7.0)]
    cuts += screws_along(0, 0.0, -1, pts, head_depth=12.0, length=18)
    cuts += pins_along(0, 0.0, [(-18.6, Z_PT + 18), (18.6, Z_PT + 18)])
    return dict(name="torso", solid=solid, cuts=cuts, plane=(0, 0.0),
                bbox=((-19, -TY - 1, Z_PT + 1.5), (20.5, TY + 1, zt + 1)), mirror=False)


def part_head():
    W = WORLD
    zt = Z_SP + 14.0
    neck = cyl(W, 4.8, zt - 0.6, zt + 6)
    face = rboxr(W, (-12, -13.5, zt + 3), (15, 13.5, zt + 24), 5.0)
    dome = I(ellipsoid(W, (33, 33, 24), (1.5, 0, zt + 19)), halfspace(W, 2, zt + 19))
    rim = rcyl(W, 33.0, zt + 17.6, zt + 19.6, 0.9)
    solid = U(SU(neck, face, 2.0), dome, rim)
    solid = I(solid, halfspace(W, 2, zt - 0.6))
    cuts = [I(ellipsoid(W, (31.8, 31.8, 22.8), (1.5, 0, zt + 19)), halfspace(W, 2, zt + 19.0)),
            rboxr(W, (-10.4, -11.9, zt + 6.5), (13.4, 11.9, zt + 21), 3.6),
            rboxr(W, (13.0, -10, zt + 8.5), (16.5, 10, zt + 16.5), 2.0),       # visor recess
            cyl(Frame((14, 0, zt + 12.5), z=(1, 0, 0), x=(0, 1, 0)), 3.6, -3, 5)]   # camera / LED hole
    for yy in (-6, 6):
        cuts.append(cyl(Frame((0, yy, zt), z=(0, 0, 1), x=(1, 0, 0)), M2_PILOT, -1, 6.5))
    return dict(name="head", solid=solid, cuts=cuts, plane=None, print_dir=(0, 0, 1),
                bbox=((-33, -35, zt - 1), (36, 35, zt + 44)), mirror=False, whole=True)


HEAD_SPLIT = Z_SP + 14.0 + 19.0


def part_head_dome():
    p = part_head()
    p.update(name="head_dome", solid=I(p["solid"], halfspace(WORLD, 2, HEAD_SPLIT)), print_dir=(0, 0, 1),
             bbox=((-35, -35, HEAD_SPLIT - 1), (36, 35, HEAD_SPLIT + 26)), mold_only=True)
    return p


def part_head_face():
    p = part_head()
    p.update(name="head_face", solid=I(p["solid"], halfspace(WORLD, 2, HEAD_SPLIT, keep_above=False)),
             plane=(0, 0.0), bbox=((-14, -15, Z_SP + 12), (17.5, 15, HEAD_SPLIT + 1)), mold_only=True)
    return p


LINKS = [part_pelvis, part_hip_pitch, part_hip_roll, part_thigh, part_shin, part_ankle, part_foot,
         part_forearm, part_torso, part_head, part_head_dome, part_head_face, part_hub, part_hornplate]


def part_fit_test():
    """Print first: one servo pocket (both shells) + an idler pin, to check your servo fits."""
    S = Frame()
    solid = U(housing(S, "yoke"))
    cuts = [pocket(S, "yoke")]
    cuts += screws_along(1, 0.0, +1, [], head_depth=3, length=8)
    return dict(name="fit_test_servo_pocket", solid=solid, cuts=cuts, plane=(1, 0.0),
                bbox=((-9, -9, -37), (20, 9, -3)), mirror=False)


LINKS.append(part_fit_test)
