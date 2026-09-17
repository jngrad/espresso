
import numpy as np

def rotate_vector(v, k, phi):
    """Rotates vector v around unit vector k by angle phi.
    Uses Rodrigues' rotation formula."""
    vrot = np.array(v) * np.cos(phi) + np.cross(k, v) * \
        np.sin(phi) + np.array(k) * np.dot(k, v) * (1.0 - np.cos(phi))
    return vrot


def dihedral_potential_and_forces(k, n, phase, p1, p2, p3, p4):
    """
    Calculate the potential and forces for a dihedral angle.
    """
    v12 = p2 - p1
    v23 = p3 - p2
    v34 = p4 - p3

    v12Xv23 = np.cross(v12, v23)
    l_v12Xv23 = np.linalg.norm(v12Xv23)
    v23Xv34 = np.cross(v23, v34)
    l_v23Xv34 = np.linalg.norm(v23Xv34)

    phi = np.arctan2(np.dot(v23, np.cross(v12Xv23, v23Xv34)),
                     np.dot(v23, v23) * np.dot(v12Xv23, v23Xv34))

    f1 = (v23Xv34 - np.cos(phi) * v12Xv23) / l_v12Xv23
    f4 = (v12Xv23 - np.cos(phi) * v23Xv34) / l_v23Xv34

    v23Xf1 = np.cross(v23, f1)
    v23Xf4 = np.cross(v23, f4)
    v34Xf4 = np.cross(v34, f4)
    v12Xf1 = np.cross(v12, f1)

    # handle singularity near TINY_SIN_VALUE
    if np.abs(np.sin(phi)) < 1e-10:
        coeff = -k * n**2 * np.cos(n * phi - phase) / np.cos(phi)
    else:
        coeff = -k * n * np.sin(n * phi - phase) / np.sin(phi)

    force1 = coeff * v23Xf1
    force2 = coeff * (v34Xf4 - v12Xf1 - v23Xf1)
    force3 = coeff * (v12Xf1 - v23Xf4 - v34Xf4)
    force4 = coeff * v23Xf4
    potential = k * (1 - np.cos(n * phi - phase))
    return (potential, (force1, force2, force3, force4))


class Particle:
    def __init__(self):
        self.pos = np.zeros(3)

def generate_traj(N=8):
    p0 = Particle()
    p1 = Particle()
    p2 = Particle()
    p3 = Particle()
    b0 = 0.16
    p1.pos = np.array([0., 0., 0.])
    p2.pos = p1.pos + [b0, 0., 0.]
    p0.pos = p1.pos + [0., b0, 0.]
    d_phi = 2 * np.pi / N
    axis = np.array([1., 0., 0.])
    dh_k = 2.
    dh_n, dh_phi0_div = (2, 3)
    dh_phi0 = np.pi / dh_phi0_div
    with open("dihedral.template.gro") as f:
        template = f.read()
    template = template.replace(" 5.000   5.000   5.000", "{: .3f}  {: .3f}  {: .3f}")
    template = template.replace(" 6.000   6.000   6.000", "{: .3f}  {: .3f}  {: .3f}")
    template = template.replace(" 7.000   7.000   7.000", "{: .3f}  {: .3f}  {: .3f}")
    template = template.replace(" 8.000   8.000   8.000", "{: .3f}  {: .3f}  {: .3f}")
    with open("dihedral.gro", "w") as f:
        for i in range(N):
            phi = i * d_phi + np.pi
            p3.pos = p2.pos + rotate_vector([0., b0, 0.], axis, phi)
            E_ref, forces_ref = dihedral_potential_and_forces(
                dh_k, dh_n, dh_phi0, p0.pos, p1.pos, p2.pos, p3.pos)
            f.write(template.format(
                p0.pos[0], p0.pos[1], p0.pos[2],
                p1.pos[0], p1.pos[1], p1.pos[2],
                p2.pos[0], p2.pos[1], p2.pos[2],
                p3.pos[0], p3.pos[1], p3.pos[2]))
    for i in range(N):
        with open(f"dihedral/frame{i:0>5}.gro", "w") as f:
            phi = i * d_phi
            p3.pos = p2.pos + rotate_vector([0., b0, 0.], axis, phi)
            E_ref, forces_ref = dihedral_potential_and_forces(
                dh_k, dh_n, dh_phi0, p0.pos, p1.pos, p2.pos, p3.pos)
            f.write(template.format(
                p0.pos[0], p0.pos[1], p0.pos[2],
                p1.pos[0], p1.pos[1], p1.pos[2],
                p2.pos[0], p2.pos[1], p2.pos[2],
                p3.pos[0], p3.pos[1], p3.pos[2]))

generate_traj(3200)
