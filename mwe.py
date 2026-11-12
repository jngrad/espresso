import numpy as np
import espressomd
system = espressomd.System(box_l=[10., 10., 10.])
system.cell_system.skin = 0.4
system.time_step = 1e-6
system.integrator.set_vv()  # no thermostat -> ROT_EULER
p = system.part.add(pos=[0., 0., 0.], rotation=3 * [True],
                    director=np.array([1., 0.5, 0.3]) / 1.34**0.5,
                    ext_torque=[1., -0.5, 2.])
for i in range(4):
    print("===")
    print(f"system.integrator.run(1, recalc_forces={(i == 2)})")
    system.integrator.run(1, recalc_forces=(i == 2))
    print(np.copy(p.omega_lab))
    np.testing.assert_allclose(np.copy(p.omega_lab),
                               np.copy(p.ext_torque) * system.time, rtol=1e-10)

system.time = 0.
system.part.clear()
system.integrator.set_symplectic_euler()
print("\n\nsystem.integrator.set_symplectic_euler()")
p = system.part.add(pos=[0., 0., 0.], rotation=3 * [True],
                    director=np.array([1., 0.5, 0.3]) / 1.34**0.5,
                    ext_torque=[1., -0.5, 2.])
for i in range(4):
    print("===")
    print(f"system.integrator.run(1, recalc_forces={(i == 2)})")
    system.integrator.run(1, recalc_forces=(i == 2))
    print(np.copy(p.omega_lab))
    np.testing.assert_allclose(np.copy(p.omega_lab),
                               np.copy(p.ext_torque) * system.time, rtol=1e-10)

