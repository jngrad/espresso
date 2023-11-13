import numpy as np
import espressomd
import matplotlib.pyplot as plt
import tqdm

np.random.seed(seed=42)
n_int = 160
n_part = 1000
kTs = [1., 2., 3., 4., 5.]
multiverse = []
for _ in kTs:
    system = espressomd.System(box_l=[15., 15., 15.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.non_bonded_inter[0, 0].wca.set_params(epsilon=1., sigma=1.)
    system.part.add(pos=np.random.random((n_part//2, 3)) * [10., 10., 10.], q=(n_part//2)*[-1])
    system.part.add(pos=np.random.random((n_part//2, 3)) * [10., 10., 10.], q=(n_part//2)*[+1])
    system.integrator.set_steepest_descent(f_max=0, gamma=1e-3, max_displacement=1. / 100.)
    system.integrator.run(100)
    system.integrator.set_vv()
    multiverse.append(system)

times = np.zeros((len(kTs), n_int,))
temperature = np.zeros((len(kTs), n_int,))

for j, (system, kT) in enumerate(zip(multiverse, kTs)):
    system.thermostat.set_langevin(kT=kT, gamma=0.2, seed=42)
    system.electrostatics.solver = espressomd.electrostatics.P3M(prefactor=1., accuracy=0.1)
    n_steps = 1
    for i in tqdm.trange(n_int, desc=f"Simulate system #{j} at kT={kT:.0f}"):
        system.integrator.run(steps=n_steps)
        times[j][i] = (times[j][i - 1] if i >= 1 else 0) + n_steps
        temperature[j][i] = 2. / 3. * system.analysis.energy()["kinetic"] / n_part
        if n_steps < 20:
            n_steps += 1

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(1, 1, 1)
for j, kT in enumerate(kTs):
    ax.plot(times[j,:], temperature[j,:], label=f"system at kT={kT:.0f}")

handles, labels = ax.get_legend_handles_labels()
ax.legend(reversed(handles), reversed(labels))
plt.ylabel("Temperature")
plt.xlabel("Simulation time")
plt.show()
