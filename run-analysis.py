import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker

class Particle:
    def __init__(self):
        self.pos = np.zeros(3)


def rotate_vector(v, k, phi):
    """Rotates vector v around unit vector k by angle phi.
    Uses Rodrigues' rotation formula."""
    vrot = np.array(v) * np.cos(phi) + np.cross(k, v) * \
        np.sin(phi) + np.array(k) * np.dot(k, v) * (1.0 - np.cos(phi))
    return vrot


def dihedral_potential_and_forces(k, n, phi, phase, p1, p2, p3, p4):
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
    #phi = np.arctan2(np.dot(v23, np.cross(v12Xv23, v23Xv34)),
    #                 np.dot(v23, v23) * np.dot(v12Xv23, v23Xv34))
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
    return (force1, force2, force3, force4)


def subsample_min_dist(points, min_dist):
    n = len(points)
    kept = np.zeros(n, dtype=bool)
    kept[0] = True
    last = points[0]
    for i in range(1, n):
        if np.linalg.norm(points[i] - last) >= min_dist:
            kept[i] = True
            last = points[i]
    return (np.array(kept),)


def generate_traj(N):
    results = []
    p0 = Particle()
    p1 = Particle()
    p2 = Particle()
    p3 = Particle()
    b0 = 0.16
    p1.pos = np.array([0., 0., 0.])
    p2.pos = p1.pos + [b0, 0., 0.]
    p0.pos = p1.pos + [0., b0, 0.]
    axis = np.array([1., 0., 0.])
    # share/gromacs/top/amber99.ff/ffbonded.itp
    #  X   CA  CB  X     9     180.0     14.64400     2  ; intrpol.bsd.on C6H6
    dh_k = 14.64400 / 11.25 / 2.
    dh_n = 2 + 1
    dh_phi0 = np.pi
    d_phi = 2 * np.pi / N
    for i in range(N):
        phi = i * d_phi
        p3.pos = p2.pos + rotate_vector([0., b0, 0.], axis, phi)
        E_ref = potential = dh_k * (1 - np.cos(dh_n * phi - dh_phi0))
        forces_ref = dihedral_potential_and_forces(
            dh_k, dh_n, phi, dh_phi0, p0.pos, p1.pos, p2.pos, p3.pos)
        p_ref = (
            np.outer(forces_ref[0], p0.pos - p1.pos)
            + np.outer(forces_ref[2], p2.pos - p1.pos)
            + np.outer(forces_ref[3], p3.pos - p1.pos)
        ) / 10**3
        results.append([phi, E_ref] + [np.linalg.norm(x) for x in forces_ref] + list(p_ref.reshape((-1,))))
    return pd.DataFrame(results, columns=["Angle", "Potential", "f0", "f1", "f2", "f3", "Pres-XX", "Pres-XY", "Pres-XY", "Pres-YX", "Pres-YY", "Pres-YZ", "Pres-ZX", "Pres-ZY", "Pres-ZZ"])


df = pd.read_fwf("dihedral/results.dat")
df["Angle"] = np.arange(3200) * 2 * np.pi / 3200

df_ref = generate_traj(3200)
a=(df["Pres-YY"] / df_ref["Pres-YY"]).to_numpy()
np.median(a[~np.isnan(a)])


fig, axs = plt.subplots(2, 2)
ax_counter = 0
axes = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]

for x in df_ref.columns:
    if x.startswith("Pres-") and np.max(np.abs(df_ref[x])) > 1e-6:
        values = 648.63065 * df_ref[x].to_numpy() * 100.
        axes[ax_counter].plot(df_ref["Angle"], values, "-", label="$" + re.sub("Pres-([A-Z]+)", r"\\sigma_{\g<1>}", x) + "$ (ESPResSo)")
        ax_counter += 1

ax_counter = 0
correct_gmx_curve = False

for x in df.columns:
    if x.startswith("Pres-") and np.max(np.abs(df[x])) > 1e-4:
        #axes[ax_counter].set_prop_cycle(None)
        values = np.copy(df[x].to_numpy())
        original_values = np.copy(values)
        if correct_gmx_curve and x == "Pres-ZZ":
            factor = 1
            for i in range(2, len(values) - 2):
                values[i] *= factor
                if original_values[i] < original_values[i-2] and original_values[i]<original_values[i+2] and ((i < (len(values) // 2 - 6)) or (i > (len(values) // 2 + 6))):
                    factor *= -1
        xdata = df["Angle"].to_numpy()
        ydata = values * 100.
        min_dist = 4 * np.sqrt((2 * np.pi / 100)**2 + ((np.max(ydata) - np.min(ydata)) / 100)**2)
        if correct_gmx_curve:
            keep = subsample_min_dist(np.cumsum(np.sqrt(np.diff(xdata, prepend=0)**2 + np.diff(ydata, prepend=0)**2)), min_dist)
            axes[ax_counter].scatter(xdata[keep], ydata[keep], s=26, facecolors='none', edgecolors='k', label="$" + re.sub("Pres-([A-Z]+)", r"\\sigma_{\g<1>}", x) + "$ (GROMACS)")
        else:
            keep = slice(1, -1)
            axes[ax_counter].plot(xdata[keep], ydata[keep], label="$" + re.sub("Pres-([A-Z]+)", r"\\sigma_{\g<1>}", x) + "$ (GROMACS)")
        ax_counter += 1

for ax_counter in range(4):
    axes[ax_counter].legend()
    axes[ax_counter].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
       lambda val,pos: r'{:.0g}$\pi$'.format(val/np.pi) if val !=0 else '0'
    ))
    axes[ax_counter].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=np.pi))

plt.tight_layout()
plt.show()


print("ZZ component in Gromacs flips sign everytime the z component is zero, for some reason. but the total pressure is 0, so it must be an artifact")


#plt.plot(df["Angle"], df["Potential"], label="Gromacs potential")
#plt.plot(df_ref["Angle"], df_ref["Potential"], "--", label="ref")
#plt.legend()
#plt.show()



