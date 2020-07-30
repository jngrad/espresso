.. _Running the simulation:

Running the simulation
======================

To run the integrator call the method
:meth:`espressomd.integrate.Integrator.run`::

    system.integrator.run(number_of_steps, recalc_forces=False, reuse_forces=False)

where ``number_of_steps`` is the number of time steps the integrator
should perform. The two main integration schemes of |es| are the Velocity Verlet algorithm
and an adaption of the Velocity Verlet algorithm to simulate an NpT ensemble.
A steepest descent implementation is also available.

.. _Velocity Verlet Algorithm:

Velocity Verlet algorithm
-------------------------

:meth:`espressomd.integrate.IntegratorHandle.set_vv`

The equations of motion for the trajectory of point-like particles read

.. math:: \dot v_i(t) = F_i(\{x_j\},v_i,t)/m_i \\ \dot x_i(t) = v_i(t),

where :math:`x_i`, :math:`v_i`, :math:`m_i` are position, velocity and mass of
particle :math:`i` and :math:`F_i(\{x_j\},v_i,t)` the forces acting on it.
These forces comprise all interactions with other particles and external fields
as well as non-deterministic contributions described in :ref:`Thermostats`.

For numerical integration, this equation is discretized to the following steps (:cite:`rapaport04` eqs. 3.5.8 - 3.5.10):

1. Calculate the velocity at the half step

   .. math:: v(t+dt/2) = v(t) + \frac{F(x(t),v(t-dt/2),t)}{m} dt/2

2. Calculate the new position

   .. math:: x(t+dt) = x(t) + v(t+dt/2) dt

3. Calculate the force based on the new position

   .. math:: F = F(x(t+dt), v(t+dt/2), t+dt)

4. Calculate the new velocity

   .. math:: v(t+dt) = v(t+dt/2) + \frac{F(x(t+dt),t+dt)}{m} dt/2

Note that this implementation of the Velocity Verlet algorithm reuses
forces in step 1. That is, they are computed once in step 3,
but used twice, in step 4 and in step 1 of the next iteration. In the first time
step after setting up, there are no forces present yet. Therefore, |es| has
to compute them before the first time step. That has two consequences:
first, random forces are redrawn, resulting in a narrower distribution
of the random forces, which we compensate by stretching. Second,
coupling forces of e.g. the lattice-Boltzmann fluid cannot be computed
and are therefore lacking in the first half time step. In order to
minimize these effects, |es| has a quite conservative heuristics to decide
whether a change makes it necessary to recompute forces before the first
time step. Therefore, calling 100 times
:meth:`espressomd.integrate.Integrator.run` with ``steps=1`` does the
same as with ``steps=100``, apart from some small calling overhead.

However, for checkpointing, there is no way for |es| to tell that the forces
that you read back in actually match the parameters that are set.
Therefore, |es| would recompute the forces before the first time step, which
makes it essentially impossible to checkpoint LB simulations, where it
is vital to keep the coupling forces. To work around this, there is
an additional parameter ``reuse_forces``, which tells integrate to not recalculate
the forces for the first time step, but use that the values still stored
with the particles. Use this only if you are absolutely sure that the
forces stored match your current setup!

The opposite problem occurs when timing interactions: In this case, one
would like to recompute the forces, despite the fact that they are
already correctly calculated. To this aim, the option ``recalc_forces`` can be used to
enforce force recalculation.

.. _Isotropic NPT integrator:

Isotropic NPT integrator
------------------------

:meth:`espressomd.integrate.IntegratorHandle.set_isotropic_npt`

As the NpT thermostat alters the way the equations of motion are integrated, it is
discussed here and only a brief summary is given in :ref:`Thermostats`.

To activate the NpT integrator, use :meth:`~espressomd.integrate.IntegratorHandle.set_isotropic_npt`
with parameters:

* ``ext_pressure``: The external pressure
* ``piston``: The mass of the applied piston
* ``direction``: Flags to enable/disable box dimensions to be subject to fluctuations. By default, all directions are enabled.

Additionally, a NpT thermostat has to be set by :meth:`~espressomd.thermostat.Thermostat.set_npt()`
with parameters:

* ``kT``: Thermal energy of the heat bath
* ``gamma0``: Friction coefficient of the bath
* ``gammav``: Artificial friction coefficient for the volume fluctuations.

A code snippet would look like::

    import espressomd

    system = espressomd.System()
    system.thermostat.set_npt(kT=1.0, gamma0=1.0, gammav=1.0, seed=42)
    system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)

The physical meaning of these parameters is described below:

The relaxation towards a desired pressure :math:`P` (parameter ``ext_pressure``)
is enabled by treating the box
volume :math:`V` as a degree of freedom with corresponding momentum :math:`\Pi = Q\dot{V}`,
where :math:`Q` (parameter ``piston``) is an artificial piston mass.
Which box dimensions are affected to change the volume can be controlled by a list of
boolean flags for parameter ``direction``.
An additional energy :math:`H_V = 1/(2Q)\Pi + PV`
associated with the volume is postulated. This results in a "force" on the box such that

.. math:: \dot{\Pi} = \mathcal{P} - P

where

.. math:: \mathcal{P} = \frac{1}{Vd} \sum_{i,j} f_{ij}x_{ij} + \frac{1}{Vd} \sum_i m_i v_i^2

Here :math:`\mathcal{P}` is the instantaneous pressure, :math:`d` the dimension
of the system (number of flags set by ``direction``), :math:`f_{ij}` the
short range interaction force between particles :math:`i` and :math:`j` and
:math:`x_{ij}= x_j - x_i`.

In addition to this deterministic force, a friction :math:`-\frac{\gamma^V}{Q}\Pi(t)`
and noise :math:`\sqrt{k_B T \gamma^V} \eta(t)` are added for the box
volume dynamics and the particle dynamics. This introduces three new parameters:
The friction coefficient for the box :math:`\gamma^V` (parameter ``gammav``),
the friction coefficient of the particles :math:`\gamma^0` (parameter ``gamma0``)
and the thermal energy :math:`k_BT` (parameter ``kT``).
For a discussion of these terms and their discretisation, see :ref:`Langevin thermostat`,
which uses the same approach, but only for particles.
As a result of box geometry changes, the particle positions and velocities have to be rescaled
during integration.

The discretisation consists of the following steps (see :cite:`kolb99a` for a full derivation of the algorithm):

1. Calculate the particle velocities at the half step

   .. math:: v'(t+dt/2) = v(t) + \frac{F(x(t),v(t-dt/2),t)}{m} dt/2

2. Calculate the instantaneous pressure and "volume momentum"

   .. math:: \mathcal{P} = \mathcal{P}(x(t),V(t),f(x(t)), v'(t+dt/2))
   .. math:: \Pi(t+dt/2) = \Pi(t) + (\mathcal{P}-P) dt/2 -\frac{\gamma^V}{Q}\Pi(t) dt/2  +  \sqrt{k_B T \gamma^V dt} \overline{\eta}

3. Calculate box volume and scaling parameter :math:`L` at half step and full step, scale the simulation box accordingly

   .. math:: V(t+dt/2) = V(t) + \frac{\Pi(t+dt/2)}{Q} dt/2
   .. math:: L(t+dt/2) = V(t+dt/2)^{1/d}
   .. math:: V(t+dt) = V(t+dt/2) + \frac{\Pi(t+dt/2)}{Q} dt/2
   .. math:: L(t+dt) = V(t+dt)^{1/d}

4. Update particle positions and scale velocities

   .. math:: x(t+dt) = \frac{L(t+dt)}{L(t)} \left[ x(t) + \frac{L^2(t)}{L^2(t+dt/2)} v(t+dt/2) dt \right]
   .. math:: v(t+dt/2) = \frac{L(t)}{L(t+dt)} v'(t+dt/2)

5. Calculate forces, instantaneous pressure and "volume momentum"

   .. math:: F = F(x(t+dt),v(t+dt/2),t)
   .. math:: \mathcal{P} = \mathcal{P}(x(t+dt),V(t+dt),f(x(t+dt)), v(t+dt/2))
   .. math:: \Pi(t+dt) = \Pi(t+dt/2) + (\mathcal{P}-P) dt/2 -\frac{\gamma^V}{Q}\Pi(t+dt/2) dt/2  +  \sqrt{k_B T \gamma^V dt} \overline{\eta}

   with uncorrelated numbers :math:`\overline{\eta}` drawn from a random uniform process :math:`\eta(t)`

6. Update the velocities

   .. math:: v(t+dt) = v(t+dt/2) + \frac{F(t+dt)}{m} dt/2

Notes:

* The NpT algorithm is only tested for all 3 directions enabled for scaling. Usage of ``direction`` is considered an experimental feature.
* In step 4, only those coordinates are scaled for which ``direction`` is set.
* For the instantaneous pressure, the same limitations of applicability hold as described in :ref:`Pressure`.
* The particle forces :math:`F` include interactions as well as a friction (:math:`\gamma^0`) and noise term (:math:`\sqrt{k_B T \gamma^0 dt} \overline{\eta}`) analogous to the terms in the :ref:`Langevin thermostat`.
* The particle forces are only calculated in step 5 and then reused in step 1 of the next iteration. See :ref:`Velocity Verlet Algorithm` for the implications of that.

.. _Rotational degrees of freedom and particle anisotropy:

Rotational degrees of freedom and particle anisotropy
-----------------------------------------------------

When the feature ``ROTATION`` is compiled in, particles not only have a position, but also an orientation that changes with an angular velocity. A torque on a particle leads to a change in angular velocity depending on the particles rotational inertia. The property :attr:`espressomd.particle_data.ParticleHandle.rinertia` has to be specified as the three eigenvalues of the particles rotational inertia tensor.

The rotational degrees of freedom are also integrated using a velocity Verlet scheme.
The implementation is based on a quaternion representation of the particle orientation and described in :cite:`omelyan98` with quaternion components indexing made according to the formalism :math:`q = a + b\mathbf{i} + c\mathbf{j} + d\mathbf{k}` :cite:`allen2017`.

When the Langevin thermostat is enabled, the rotational degrees of freedom are also thermalized.

Whether or not rotational degrees of freedom are propagated, is controlled on a per-particle and per-axis level, where the axes are the Cartesian axes of the particle in its body-fixed frame.
It is important to note that starting from version 4.0 and unlike in earlier versions of |es|, the particles' rotation is disabled by default.
In this way, just compiling in the ``ROTATION`` feature no longer changes the physics of the system.

The rotation of a particle is controlled via the :attr:`espressomd.particle_data.ParticleHandle.rotation` property. E.g., the following code adds a particle with rotation enabled on the x axis::

    import espressomd
    system = espressomd.System()
    system.part.add(pos=(0, 0, 0), rotation=(1, 0, 0))

Notes:

* The orientation of a particle is stored as a quaternion in the :attr:`espressomd.particle_data.ParticleHandle.quat` property. For a value of (1,0,0,0), the body and space frames coincide.
* The space-frame direction of the particle's z-axis in its body frame is accessible through the :attr:`espressomd.particle_data.ParticleHandle.director` property.
* Any other vector can be converted from body to space fixed frame using the :meth:`espressomd.particle_data.ParticleHandle.convert_vector_body_to_space` method.
* When ``DIPOLES`` are compiled in, the particles dipole moment is always co-aligned with the z-axis in the body-fixed frame.
* Changing the particles dipole moment and director will re-orient the particle such that its z-axis in space frame is aligned parallel to the given vector. No guarantees are made for the other two axes after setting the director or the dipole moment.


The following particle properties are related to rotation:

* :attr:`espressomd.particle_data.ParticleHandle.dip`
* :attr:`espressomd.particle_data.ParticleHandle.director`
* :attr:`espressomd.particle_data.ParticleHandle.ext_torque`
* :attr:`espressomd.particle_data.ParticleHandle.gamma_rot`
* :attr:`espressomd.particle_data.ParticleHandle.gamma_rot`
* :attr:`espressomd.particle_data.ParticleHandle.omega_body`
* :attr:`espressomd.particle_data.ParticleHandle.omega_lab`
* :attr:`espressomd.particle_data.ParticleHandle.quat`
* :attr:`espressomd.particle_data.ParticleHandle.rinertia`
* :attr:`espressomd.particle_data.ParticleHandle.rotation`
* :attr:`espressomd.particle_data.ParticleHandle.torque_lab`

.. _Steepest descent:

Steepest descent
----------------

:meth:`espressomd.integrate.IntegratorHandle.set_steepest_descent`

This feature is used to propagate each particle by a small distance parallel to the force acting on it.
When only conservative forces for which a potential exists are in use, this is equivalent to a steepest descent energy minimization.
A common application is removing overlap between randomly placed particles.

Please note that the behavior is undefined if a thermostat is activated.
It runs a simple steepest descent algorithm:

.. math:: \vec{r}_{i+1} = \vec{r}_i + \min(\gamma \vec{F}_i, \vec{r}_{\text{max_displacement}}),

while the maximal force is bigger than ``f_max`` or for at most ``max_steps`` times. The energy
is relaxed by ``gamma``, while the change per coordinate per step is limited to ``max_displacement``.
The combination of ``gamma`` and ``max_displacement`` can be used to get a poor man's adaptive update.
Rotational degrees of freedom are treated similarly: each particle is
rotated around an axis parallel to the torque acting on the particle.
Please be aware of the fact that this needs not to converge to a local
minimum in periodic boundary conditions. Translational and rotational
coordinates that are fixed using the ``fix`` and ``rotation`` attribute of particles are not altered.

Usage example::

    system.integrator.set_steepest_descent(
        f_max=0, gamma=0.1, max_displacement=0.1)
    system.integrator.run(20)
    system.integrator.set_vv()  # to switch back to velocity Verlet

A wrapper function :func:`~espressomd.minimize_energy.steepest_descent` is
available to set up the steepest descent integrator while preserving the
original integrator. The snippet above can be rewritten to::

    from espressomd.minimize_energy import steepest_descent
    steepest_descent(system, f_max=0, gamma=0.1, max_displacement=0.1, max_steps=20)

This convenience function only exists for historical reasons and remains available
for backward compatibility. New scripts should setup the steepest descent
integrator with the :meth:`~espressomd.integrate.IntegratorHandle.set_steepest_descent`
handle directly.

.. _Stokesian Dynamics:

Stokesian Dynamics
------------------

.. note::

    Requires ``STOKESIAN_DYNAMICS_LIB_*`` external features, enabled with
    ``-DWITH_STOKESIAN_DYNAMICS_CPU=ON`` or ``-DWITH_STOKESIAN_DYNAMICS_GPU=ON``,
    and features ``STOKESIAN_DYNAMICS`` or ``STOKESIAN_DYNAMICS_GPU``.

:meth:`espressomd.integrate.IntegratorHandle.set_stokesian_dynamics`

The Stokesian Dynamics method allows to study the behaviour of spherical
particles in a viscous fluid. It is targeted at systems with very low Reynolds
numbers. In such systems, particles come to a rest almost immediately as soon as
any force on them is removed. In other words, motion has no memory of the past.

The integration scheme is relatively simple. Only the particles' positions,
radii and forces (including torques) are needed to compute the momentary
velocities (including angular velocities). The particle positions are
integrated by the simple Euler scheme.

The computation of the velocities is an approximation with good results
in the far field.
The Stokesian Dynamics method is only available for open systems,
i.e. no periodic boundary conditions are supported. The box size has
no effect either.

The Stokesian Dynamics method is outlined in :cite:`durlofsky87a`.

The following minimal example illustrates how to use the SDM in |es|::

    import espressomd
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.part.add(pos=[0, 0, 0], rotation=[1, 0, 0])
    system.integrator.set_stokesian_dynamics(viscosity=1.0, radii={0: 1.0})

    system.integrator.run(100)

Because there is no force on the particle yet, nothing will move. You will need
to add your own actors to the system. The parameter ``radii`` is a dictionary
that maps particle types to different radii. ``viscosity`` is the dynamic
viscosity of the ambient infinite fluid. There are additional optional
parameters for ``set_stokesian_dynamics()``. For more information, see
:py:meth:`espressomd.integrate.IntegratorHandle.set_stokesian_dynamics()`.

Note that this setup represents a system at zero temperature. In order to
thermalize the system, the SD thermostat needs to be activated (see
:ref:`Stokesian thermostat`).

.. _Important_SD:

Important
~~~~~~~~~

The particles must be prevented from overlapping. It is mathematically allowed
for the particles to overlap to a certain degree. However, once the distance
of the sphere centers is less than 2/3 of the sphere diameter, the mobility
matrix is no longer positive definite and the Stokesian Dynamics integrator
will fail. Therefore, the particle centers must be kept apart from each
other by a strongly repulsive potential, for example the WCA potential
that is set to the appropriate particle radius (for more information about
the available interaction types see :ref:`Non-bonded interactions`).

The current implementation of SD only includes the far field approximation.
The near field (so-called lubrication) correction is planned. For now,
Stokesian Dynamics provides a good approximation of the hydrodynamics
in dilute systems where the average distance between particles is several
sphere diameters.
