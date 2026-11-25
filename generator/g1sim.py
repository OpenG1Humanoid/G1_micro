"""MuJoCo model of the G1 micro, made from Unitree's own g1_23dof_rev_1_0.xml:
   * every length scaled by S (body positions, inertial positions, meshes, contact spheres)
   * wrist roll joints removed (fixed), as on the micro
   * link masses replaced by the micro estimate (printed plastic + servos), inertias rescaled
   * actuators: position servos with micro-servo torque limits
Meshes are the original G1 STLs scaled in the XML, so the model looks exactly like the robot.
"""
import os, json, copy, shutil
import xml.etree.ElementTree as ET
import numpy as np
import g1

S = 0.35
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(g1.G1, "g1_23dof_rev_1_0.xml")
SERVO_G = {"STS3032": 20.0, "SCS0009": 13.2}
LEG_TORQUE, ARM_TORQUE = 0.44, 0.225     # N m stall (4.5 / 2.3 kg cm)


def scale_vec(s, k):
    return " ".join(f"{float(v) * k:.6g}" for v in s.split())


def main(link_mass_g):
    tree = ET.parse(SRC); root = tree.getroot()
    out_dir = os.path.join(ROOT, "sim"); shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(os.path.join(out_dir, "meshes"))
    comp = root.find("compiler"); comp.set("meshdir", "meshes")
    used = set()
    for m in root.iter("mesh"):
        m.set("scale", f"{S} {S} {S}"); used.add(m.get("file"))
    for f in used:
        shutil.copy(os.path.join(g1.G1, "meshes", f), os.path.join(out_dir, "meshes", f))
    total_g1 = 0.0
    for b in root.iter("body"):
        if b.get("pos"): b.set("pos", scale_vec(b.get("pos"), S))
        for g in b.findall("geom"):
            if g.get("size") and g.get("type") is None:
                g.set("size", scale_vec(g.get("size"), S))
            if g.get("pos"): g.set("pos", scale_vec(g.get("pos"), S))
        for st in b.findall("site"):
            if st.get("pos"): st.set("pos", scale_vec(st.get("pos"), S))
            if st.get("size"): st.set("size", scale_vec(st.get("size"), S))
        inn = b.find("inertial")
        if inn is not None:
            m_old = float(inn.get("mass")); total_g1 += m_old
            name = b.get("name")
            key = name if name in link_mass_g else name.replace("_link", "")
            m_new = link_mass_g.get(name, link_mass_g.get(key, None))
            if m_new is None: m_new = 1000 * m_old * 0.75 / 35.0
            m_new /= 1000.0
            inn.set("pos", scale_vec(inn.get("pos"), S)); inn.set("mass", f"{m_new:.4f}")
            if inn.get("diaginertia"):
                inn.set("diaginertia", scale_vec(inn.get("diaginertia"), (m_new / m_old) * S * S))
    # remove wrist roll joints (the hand is fixed to the forearm on the micro)
    for b in root.iter("body"):
        for j in list(b.findall("joint")):
            if "wrist" in j.get("name", ""):
                b.remove(j)
        for j in b.findall("joint"):
            if j.get("type") != "free":
                j.attrib.pop("actuatorfrcrange", None)
                j.set("damping", "0.03"); j.set("armature", "0.0015"); j.set("frictionloss", "0.01")
    act = root.find("actuator")
    if act is not None: root.remove(act)
    act = ET.SubElement(root, "actuator")
    for j in root.iter("joint"):
        n = j.get("name", "")
        if j.get("type") == "free" or not n: continue
        arm = any(t in n for t in ("shoulder", "elbow"))
        tq = ARM_TORQUE if arm and "shoulder_pitch" not in n else LEG_TORQUE
        lo, hi = j.get("range").split()
        ET.SubElement(act, "position", name=n.replace("_joint", ""), joint=n, kp=f"{2.5 if arm else 6.0}",
                      forcerange=f"-{tq} {tq}", ctrlrange=f"{lo} {hi}")
    for s in root.iter("sensor"):
        root.remove(s); break
    # floor
    wb = root.find("worldbody")
    if not any(g.get("type") == "plane" for g in wb.findall("geom")):
        fl = ET.Element("geom", name="floor", type="plane", size="2 2 0.05", rgba="0.85 0.85 0.85 1")
        wb.insert(0, fl); wb.insert(0, ET.Element("light", pos="0 0 2", dir="0 0 -1"))
    ET.indent(tree)
    tree.write(os.path.join(out_dir, "g1_micro.xml"))
    return os.path.join(out_dir, "g1_micro.xml")


if __name__ == "__main__":
    lm = json.load(open(os.path.join(ROOT, "report", "link_mass_g.json")))
    print(main(lm))
