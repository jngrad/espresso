#
# Copyright (C) 2023 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Rapid prototyping of new Monte Carlo methods in Python.

This sample provides a re-implementation of the core functionality
of :ref:`reaction methods <Reaction methods>` in Python,
with a focus on the :ref:`constant pH <Constant pH>` and
:ref:`reaction ensemble <Reaction Ensemble>` methods.
See :ref:`Writing new Monte Carlo methods` for more details.
The sample simulates the acid-based titration of polyelectrolyte chains.

The sample is designed to run with the :ref:`kernprof` profiler attached:

.. code-block:: bash

    pypresso --kernprof monte_carlo.py --mode=core
    pypresso --kernprof monte_carlo.py --mode=python

"""

import numpy as np
import itertools
import argparse
import math
import time
import tqdm
import os
import sys

import espressomd
import espressomd.polymer
import espressomd.electrostatics
import espressomd.reaction_methods
import espressomd.utils
import espressomd.code_features

required_features = ["P3M", "WCA"]
espressomd.assert_features(required_features)

parser = argparse.ArgumentParser(
    prog=f"pypresso --kernprof {os.path.basename(__file__)}",
    epilog=__doc__.lstrip().split("\n", 1)[0])
parser.add_argument("--mode", choices=["core", "python"], default="python",
                    help="use C++ (core) or Python (python) implementation")
parser.add_argument("--method", choices=["cph", "re","widom"], default="cph",
                    help="use constant pH (cph) or reaction ensemble (re)")
args = parser.parse_args()

if "line_profiler" not in dir():
    def profile(func):
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper



class MCMethods:
    _so_name = "MCMethods"
    def __init__(self, **kwargs):
        self.kT = kwargs["kT"]
        self.rng = np.random.default_rng(seed=kwargs["seed"])
        self.system = kwargs["system"]
        self.particle_inside_exclusion_range_touched = False
        self.exclusion = espressomd.reaction_methods.ExclusionRadius(**kwargs)
        if 'exclusion_radius' in kwargs:
                raise KeyError(
                    'the keyword `exclusion_radius` is obsolete. Currently, the equivalent keyword is `exclusion_range`')
        if not 'sip' in kwargs:
            espressomd.utils.check_valid_keys(self.valid_keys(), kwargs.keys())
        self.constraint_type="none"
    
    def required_keys(self):
        return {"kT", "seed"}

    def valid_keys(self):
        return {"kT", "exclusion_range", "seed",  "exclusion_radius_per_type"}

    def calculate_acceptance_probability(
            self, reaction, E_pot_diff, old_particle_numbers):
        raise NotImplementedError("Derived classes must implement this method")

    def check_exclusion_range(self, pid):
        self.particle_inside_exclusion_range_touched |= self.exclusion.check_exclusion_range(
            pid=pid)

    def get_random_position_in_box(self):
        """
        Returns a random position in the simulation box within the input boundaries

        Note:
            - method='slab' currently only supports a slab in the z-direction.
            - method='cylinder' currently only supports a cylinder aligned in the z-direction.
        """
        supported_constraint_types=["none","cylinder","slab"]
        if self.constraint_type in supported_constraint_types:
            raise ValueError(f"Constraint type {self.constraint_type="none"} is not currently supported, supported types are {supported_constraint_types}")
        box_l = system.box_l
        position = []
        if self.constraint_type == "none":
            for side in box_l:
                position.append(side*self.rng.uniform())
        elif self.constraint_type == "slab":
            for side in box_l[:2]:
                position.append(side*self.rng.uniform())
            coord_z=self.params_boundaries["slab_start_z"]+self.rng.uniform()*(self.params_boundaries["slab_end_z"]-self.params_boundaries["slab_start_z"])
            position.append(coord_z)
        elif self.constraint_type == "cylinder":
            radius=self.params_boundaries["radius"]*np.sqrt(self.rng.uniform())
            phi=2*np.pi*self.rng.uniform()
            position.append(box_l[0]+self.params_boundaries["center_x"]*np.cos(phi))
            position.append(box_l[1]+self.params_boundaries["center_y"]*np.sin(phi))
            position.append(box_l[2]*self.rng.uniform())
        return position

    def get_random_pids(self, ptype, size):
        pids = self.system.analysis.call_method(
            "get_pids_of_type", ptype=ptype)
        indices = self.rng.choice(len(pids), size=size, replace=False)
        return [pids[i] for i in indices]

    def set_cyl_constraint(self,center_x, center_y,  radius):
        """
        Constrains the MC sampling within a volume given by a cylender.

        NOTE:
            This function assumes that the cylender is along the z-axis

        """
        if center_x < 0 or center_x > self.system.box_l[0]:
            raise ValueError(f"center_x is outside the box")
        if center_y < 0 or center_y > self.system.box_l[1]:
            raise ValueError(f"center_y is outside the box")
        if radius < 0:
            raise ValueError(f"radius is invalid")

        self.constraint_type="cylinder"
        self.params_boundaries={"radius": radius,
                                "center_x": center_x,
                                "center_y": center_y}

    def set_slab_constraint(self, slab_start_z,slab_end_z):
        """
        Constrains the MC sampling within a volume given by a cylender.

        NOTE:
            This function assumes that the cylender is along the z-axis

        """
        if slab_start_z < 0 or slab_start_z > self.system.box_l[2]:
            raise ValueError("slab_start_z is outside the box")
        if slab_end_z < 0 or slab_end_z > self.system.box_l[2]:
            raise ValueError("slab_end_z is outside the box")
        if slab_end_z < slab_start_z:
            raise ValueError("slab_end_z must be >= slab_start_z")

        self.constraint_type="slab"
        self.params_boundaries={"slab_start_z": slab_start_z,
                                "slab_end_z": slab_end_z}
    

class SingleReaction:
    _so_name = "MCMethods::SingleReaction"
    def __init__(self, **kwargs):
        self.reactant_types = kwargs["reactant_types"]
        self.reactant_coefficients = kwargs["reactant_coefficients"]
        self.product_types = kwargs["product_types"]
        self.product_coefficients = kwargs["product_coefficients"]
        self.gamma = kwargs["gamma"]
        self.accepted_moves = 0
        self.trial_moves = 0
        self.accumulator_potential_energy_difference_exponential = []
        self.nu_bar = sum(self.product_coefficients) - \
            sum(self.reactant_coefficients)

    def get_acceptance_rate(self):
        return self.accepted_moves / self.trial_moves

    def make_backward_reaction(self):
        return SingleReaction(
            gamma=1. / self.gamma, reactant_types=self.product_types,
            reactant_coefficients=self.product_coefficients,
            product_types=self.reactant_types,
            product_coefficients=self.reactant_coefficients)


class ReactionAlgorithm(MCMethods):
    _so_name = "MCMethods::ReactionAlgorithm"
    def __init__(self, **kwargs):
        self.system = kwargs["system"]
        self.kT = kwargs["kT"]
        self.non_interacting_type = 100
        self.reactions = []
        self.particle_inside_exclusion_range_touched = False
        self.default_charges = {}
        self.m_empty_p_ids_smaller_than_max_seen_particle = []
        self.rng = np.random.default_rng(seed=kwargs["seed"])
        self.exclusion = espressomd.reaction_methods.ExclusionRadius(**kwargs)
        self.constraint_type="none"

        if self._so_name == ReactionAlgorithm._so_name:
            raise RuntimeError(
                    f"Base class '{self.__class__.__name__}' cannot be instantiated")
            if 'exclusion_radius' in kwargs:
                raise KeyError(
                    'the keyword `exclusion_radius` is obsolete. Currently, the equivalent keyword is `exclusion_range`')
            if not 'sip' in kwargs:
                espressomd.utils.check_valid_keys(self.valid_keys(), kwargs.keys())
        self.inicialize_particle_changes()
        
    def inicialize_particle_changes(self):
        self.particle_changes={"created":[],
                                "changed":[],
                                "hidden":[]}

    def valid_keys(self):
        return {"kT", "exclusion_range", "seed",
                "exclusion_radius_per_type", "search_algorithm"}

    def required_keys(self):
        return {"kT", "exclusion_range", "seed"}

    def set_non_interacting_type(self, type):
        self.non_interacting_type = type

    def add_reaction(self, **kwargs):
        """
        Set up a reaction in the forward and backward directions.
        """
        default_charges = kwargs.pop("default_charges")
        neutrality_check = kwargs.pop("check_for_electroneutrality", True)
        if not isinstance(default_charges, dict):
            raise TypeError("Argument 'default_charges' needs to be a dict")
        self.default_charges.update()
        forward_reaction = SingleReaction(**kwargs)
        backward_reaction = forward_reaction.make_backward_reaction()
        if neutrality_check:
            self._check_charge_neutrality(
                type2charge=default_charges,
                reaction=forward_reaction)
        self.default_charges=default_charges
        self.reactions.append(forward_reaction)
        self.reactions.append(backward_reaction)
        self.check_reaction_method()

    def delete_reaction(self, **kwargs):
        """
        Delete a reaction from the set of used reactions
        (the forward and backward reaction).
        The ``reaction_id`` which is assigned to a reaction
        depends on the order in which :meth:`add_reaction` was called.
        The 0th reaction has ``reaction_id=0``, the next added
        reaction needs to be addressed with ``reaction_id=1``, etc.
        After the deletion of a reaction subsequent reactions
        take the ``reaction_id`` of the deleted reaction.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Reaction id
        """
        del self.reactions[kwargs["reaction_id"]]

    @profile
    def reaction(self, steps):
        """
        Perform reaction steps. Chemical reactions are selected at random.
        """
        self.setup_bookkeeping_of_empty_pids()
        E_pot = self.system.analysis.potential_energy()
        random = self.rng.choice(len(self.reactions), size=steps, replace=True)
        for i in random:
            E_pot = self.generic_oneway_reaction(self.reactions[i], E_pot)

    @profile
    def generic_oneway_reaction(self, reaction, E_pot_old):
        """
        Carry out a generic one-way chemical reaction of the type
        `A + B + ... --> C + D + ...` and return the new potential
        energy if the trial move is accepted.
        """
        try:
            reaction.trial_moves += 1
            self.particle_inside_exclusion_range_touched = False
            if not self.all_reactant_particles_exist(reaction):
                return E_pot_old

            old_particle_numbers = self.save_old_particle_numbers(reaction)
            self.make_reaction_attempt(reaction)

            if self.particle_inside_exclusion_range_touched:  # reject trial move
                self.restore_system()
                self.particle_inside_exclusion_range_touched = False
                return E_pot_old

            E_pot_new = self.system.analysis.potential_energy()
            E_pot_diff = E_pot_new - E_pot_old
            bf = self.calculate_acceptance_probability(
                reaction, E_pot_diff, old_particle_numbers)
            reaction.accumulator_potential_energy_difference_exponential.append(
                math.exp(-E_pot_diff / self.kT))
            if self.rng.uniform() < bf:  # accept trial move
                self.delete_hidden_particles()
                reaction.accepted_moves += 1
                self.inicialize_particle_changes()
                return E_pot_new
            else:  # reject trial move
                self.restore_system()
                return E_pot_old
        except BaseException as err:
            tb = sys.exc_info()[2]
            raise RuntimeError(
                "An exception was raised by a chemical reaction; the particle "
                "state tracking is no longer guaranteed to be correct! -- "
                f"{err}").with_traceback(tb)

    @class_method
    def _factorial_Ni0_by_factorial_Ni0_plus_nu_i(cls,nu_i, N_i0):
        value = 1.
        if nu_i > 0:
            value /= math.factorial(N_i0 + nu_i) // math.factorial(N_i0)
        elif nu_i < 0:
            value *= math.factorial(N_i0) // math.factorial(N_i0 + nu_i)
        return value


    @profile
    def make_reaction_attempt(self, reaction):
        """
        Carry out a chemical reaction and save the old system state.
        """
        minimum_number_of_types = min(len(reaction.reactant_types),
                                      len(reaction.product_types))
        maximum_number_of_types = max(len(reaction.reactant_types),
                                      len(reaction.product_types))
        
        for index in range(minimum_number_of_types):
            r_type = reaction.reactant_types[index]
            p_type = reaction.product_types[index]
            r_charge = self.default_charges[r_type]
            p_charge = self.default_charges[p_type]

            # change reactant particles to product particles
            size = min(reaction.reactant_coefficients[index],
                       reaction.product_coefficients[index])
            for random_pid in self.get_random_pids(r_type, size):
                p = self.system.part.by_id(random_pid)
                p.type = p_type
                p.q = p_charge
                self.particle_changes["changed"].append(
                    {"pid": random_pid, "type": r_type, "charge": r_charge})

            # measure stoichiometric excess
            delta_n = reaction.product_coefficients[index] - \
                reaction.reactant_coefficients[index]

            if delta_n > 0:
                # create product particles
                for _ in range(delta_n):
                    pid = self.create_particle(p_type)
                    self.check_exclusion_range(pid)
                    self.particle_changes["created"].append(
                        {"pid": pid, "type": p_type, "charge": p_charge})
            elif delta_n < 0:
                # hide reactant particles
                for random_pid in self.get_random_pids(r_type, -delta_n):
                    self.particle_changes["hidden"].append(
                        {"pid": random_pid, "type": r_type, "charge": r_charge})
                    self.check_exclusion_range(random_pid)
                    self.hide_particle(random_pid)

        # create/hide particles with non-corresponding replacement types
        for index in range(minimum_number_of_types, maximum_number_of_types):
            if len(reaction.product_types) < len(reaction.reactant_types):
                r_type = reaction.reactant_types[index]
                r_charge = self.default_charges[r_type]
                size = reaction.reactant_coefficients[index]
                # hide superfluous reactant particles
                for random_pid in self.get_random_pids(r_type, size):
                    self.particle_changes["hidden"].append(
                        {"pid": random_pid, "type": r_type, "charge": r_charge})
                    self.check_exclusion_range(random_pid)
                    self.hide_particle(random_pid)
            else:
                p_type = reaction.product_types[index]
                p_charge = self.default_charges[p_type]
                # create additional product particles
                for _ in range(reaction.product_coefficients[index]):
                    pid = self.create_particle(p_type)
                    self.check_exclusion_range(pid)
                    self.particle_changes["created"].append(
                        {"pid": pid, "type": p_type, "charge": p_charge})



    def all_reactant_particles_exist(self, reaction):
        for r_type in reaction.reactant_types:
            r_index = reaction.reactant_types.index(r_type)
            r_coef = reaction.reactant_coefficients[r_index]
            if self.system.number_of_particles(type=r_type) < r_coef:
                return False
        return True

    def save_old_particle_numbers(self, reaction):
        old_particle_numbers = {}
        for r_type in reaction.reactant_types + reaction.product_types:
            old_particle_numbers[r_type] = self.system.number_of_particles(
                type=r_type)
        return old_particle_numbers

    def delete_created_particles(self):
        for particle_info in self.particle_changes["created"]:
            self.system.part.by_id(particle_info["pid"]).remove()

    def delete_hidden_particles(self):
        for particle_info in self.particle_changes["hidden"]:
            self.system.part.by_id(particle_info["pid"]).remove()

    def restore_system(self):
        # restore properties of changed and hidden particles
        for particle_info in self.particle_changes["changed"] + \
                self.particle_changes["hidden"]:
            p = self.system.part.by_id(particle_info["pid"])
            p.type = particle_info["type"]
            p.q = particle_info["charge"]
        # destroy created particles
        self.delete_created_particles()
        self.inicialize_particle_changes()

    def hide_particle(self, pid):
        p = self.system.part.by_id(pid)
        p.type = self.non_interacting_type
        p.q = 0.

    def create_particle(self, ptype):
        if len(self.m_empty_p_ids_smaller_than_max_seen_particle) == 0:
            pid = self.system.part.highest_particle_id + 1
        else:
            pid = min(self.m_empty_p_ids_smaller_than_max_seen_particle)
            self.m_empty_p_ids_smaller_than_max_seen_particle.remove(pid)
        self.system.part.add(id=pid, type=ptype, q=self.default_charges[ptype],
                             pos=self.rng.random((3,)) * self.system.box_l,
                             v=self.rng.normal(size=3) * math.sqrt(self.kT))
        return pid

    def setup_bookkeeping_of_empty_pids(self):
        particle_ids = self.system.part.all().id
        available_pids = self.find_missing_pids(pids_list=particle_ids)
        self.m_empty_p_ids_smaller_than_max_seen_particle = available_pids


    def find_missing_pids(self,pids_list):
        """
        Finds the missing particles ids in `pids_list`.
        NOTE:  `pids_list` must be a sorted list [0,1,3,5,7..]
        """
        return [i for x, y in zip(pids_list, pids_list[1:]) for i in range(x + 1, y) if y - x > 1]
    
    def check_reaction_method(self):
        if len(self.reactions) == 0:
            raise RuntimeError("Reaction system not initialized")

        # charges of all reactive types need to be known
        if espressomd.code_features.has_features("ELECTROSTATICS"):
            for reaction in self.reactions:
                for p_type in reaction.reactant_types:
                    if p_type not in self.default_charges:
                        raise RuntimeError(
                            f"Forgot to assign charge to type {p_type}")

    def _check_charge_neutrality(self, type2charge, reaction):
        charges = np.array(list(type2charge.values()))
        if np.count_nonzero(charges) == 0:
            # all particles have zero charge
            # no need to check electroneutrality
            return
        # calculate net change of electrical charge for the reaction
        net_charge_change = 0.0
        for coef, ptype in zip(
                reaction.reactant_coefficients, reaction.reactant_types):
            net_charge_change -= coef * type2charge[ptype]
        for coef, ptype in zip(
                reaction.product_coefficients, reaction.product_types):
            net_charge_change += coef * type2charge[ptype]
        min_abs_nonzero_charge = np.min(
            np.abs(charges[np.nonzero(charges)[0]]))
        if abs(net_charge_change) / min_abs_nonzero_charge > 1e-10:
            raise ValueError("Reaction system is not charge neutral")

    def _check_reaction_index(self, reaction_index):
        if reaction_index < 0 or reaction_index >= len(self.reactions):
            raise IndexError(f"No reaction with id {reaction_index}")

    def get_status(self):
        """
        Returns the status of the reaction ensemble in a dictionary containing
        the used reactions, the used kT and the used exclusion radius.

        """

        self.check_reaction_method()
        property_keys = {"reactant_coefficients", "reactant_types",
                         "product_coefficients", "product_types", "gamma"}
        reactions_list = [{key: getattr(reaction, key) for key in property_keys}
                          for reaction in self.reactions]

        return {"reactions": reactions_list, "kT": self.kT,
                "exclusion_range": self.exclusion_range,
                "exclusion_radius_per_type": self.exclusion_radius_per_type}

class ReactionEnsemble(ReactionAlgorithm):
    """
    This class implements the Reaction Ensemble.
    """
    _so_name = "MCMethods::ReactionEnsemble"
    def calculate_acceptance_probability(
            self, reaction, E_pot_diff, old_particle_numbers):
        """
        Calculate the acceptance probability of a Monte Carlo move.
        """

        volume = self.system.volume()
        expr = math.exp(-E_pot_diff / self.kT)
        expr *= volume**reaction.nu_bar * reaction.gamma

        # factorial contribution of reactants
        for i in range(len(reaction.reactant_types)):
            nu_i = -reaction.reactant_coefficients[i]
            N_i0 = old_particle_numbers[reaction.reactant_types[i]]
            expr *= self._factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0)

        # factorial contribution of products
        for i in range(len(reaction.product_types)):
            nu_i = reaction.product_coefficients[i]
            N_i0 = old_particle_numbers[reaction.product_types[i]]
            expr *= self._factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0)

        return expr


class ConstantpHEnsemble(ReactionAlgorithm):
    """
    This class implements the constant pH Ensemble.

    When adding an acid-base reaction, the acid and base particle types
    are always assumed to be at index 0 of the lists passed to arguments
    ``reactant_types`` and ``product_types``.

    Attributes
    ----------
    constant_pH : :obj:`float`
        Constant pH value.

    """
    _so_name = "MCMethods::ConstantpHEnsemble"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.constant_pH = kwargs["constant_pH"]

    def valid_keys(self):
        return {"kT", "exclusion_range", "seed",
                "constant_pH", "exclusion_radius_per_type", "search_algorithm"}

    def required_keys(self):
        return {"kT", "exclusion_range", "seed", "constant_pH"}

    def add_reaction(self, **kwargs):
        kwargs["reactant_coefficients"] = [1]
        kwargs["product_coefficients"] = [1, 1]
        super().add_reaction(**kwargs)

    def calculate_acceptance_probability(self, reaction, E_pot_diff, old_particle_numbers):
        """
        Calculate the acceptance probability of a Monte Carlo move.
        """

        ln_bf = E_pot_diff - reaction.nu_bar * self.kT * math.log(10.) * (
            self.constant_pH + reaction.nu_bar * math.log10(reaction.gamma))

        factorial_expr = math.exp(-ln_bf / self.kT)

        # factorial contribution of reactants
        nu_i = -reaction.reactant_coefficients[0]
        N_i0 = old_particle_numbers[reaction.reactant_types[0]]
        factorial_expr *= self._factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0)

        # factorial contribution of products
        nu_i = reaction.product_coefficients[0]
        N_i0 = old_particle_numbers[reaction.product_types[0]]
        factorial_expr *= self._factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0)

        return factorial_expr

class WidomInsertion(ReactionAlgorithm):
    """
    This class implements the Widom insertion method in the canonical ensemble
    for homogeneous systems, where the excess chemical potential is not
    depending on the location.

    """

    _so_name = "MCMethods::WidomInsertion"

    def required_keys(self):
        return {"kT", "seed"}

    def valid_keys(self):
        return {"kT", "seed"}

    def add_reaction(self, **kwargs):
        kwargs['gamma'] = 1.
        super().add_reaction(**kwargs)

    def calculate_particle_insertion_potential_energy(self, **kwargs):
        """
        Measures the potential energy when particles are inserted in the
        system following the reaction provided in ``reaction_id``. Please
        define the insertion moves by calling the method
        :meth:`~ReactionAlgorithm.add_reaction` (with only product types
        specified).

        Note that although this function does not provide directly
        the chemical potential, it can be used to calculate it.
        For an example of such an application please check
        :file:`/samples/widom_insertion.py`.

        Parameters
        ----------
        reaction_id : :obj:`int`
            Reaction identifier. Will be multiplied by 2 internally to
            skip reverse reactions, i.e. deletion reactions!

        Returns
        -------
        :obj:`float`
            The particle insertion potential energy.

        """
        reaction_id = kwargs.pop("reaction_id")
        reaction = self.reactions[reaction_id]
        if not self.all_reactant_particles_exist(reaction):
            raise RuntimeError("Trying to remove some non-existing particles "
                               "from the system via the inverse Widom scheme.")
        self.setup_bookkeeping_of_empty_pids()
        E_pot_old = self.system.analysis.potential_energy()
        self.make_reaction_attempt(reaction)
        E_pot_new = self.system.analysis.potential_energy()
        self.restore_system()
        return E_pot_new - E_pot_old

    def calculate_excess_chemical_potential(self, **kwargs):
        """
        Given a set of samples of the particle insertion potential energy,
        calculates the excess chemical potential and its statistical error.

        Parameters
        ----------
        particle_insertion_potential_energy_samples : array_like of :obj:`float`
            Samples of the particle insertion potential energy.
        N_blocks : :obj:`int`, optional
            Number of bins for binning analysis.

        Returns
        -------
        mean : :obj:`float`
            Mean excess chemical potential.
        error : :obj:`float`
            Standard error of the mean.

        """

        def do_block_analysis(samples, N_blocks):
            """
            Performs a binning analysis of samples.
            Divides the samples in ``N_blocks`` equispaced blocks
            and returns the mean and its uncertainty
            """
            size_block = int(len(samples) / N_blocks)
            block_list = []
            for block in range(N_blocks):
                block_list.append(
                    np.mean(samples[block * size_block:(block + 1) * size_block]))

            sample_mean = np.mean(block_list)
            sample_std = np.std(block_list, ddof=1)
            sample_uncertainty = sample_std / np.sqrt(N_blocks)

            return sample_mean, sample_uncertainty

        kT = self.kT

        gamma_samples = np.exp(-1.0 * np.array(
            kwargs["particle_insertion_potential_energy_samples"]) / kT)

        gamma_mean, gamma_std = do_block_analysis(
            samples=gamma_samples, N_blocks=kwargs.get("N_blocks", 16))

        mu_ex_mean = -kT * np.log(gamma_mean)

        # full propagation of error
        mu_ex_Delta = 0.5 * kT * abs(-np.log(gamma_mean + gamma_std) -
                                     (-np.log(gamma_mean - gamma_std)))

        return mu_ex_mean, mu_ex_Delta

class CanonicalEnsemble():
    _so_name = "MCMethods::CanonicalEnsemble"
    
    def move_particle_in_simulation_box(self,ptype,steps):
        """
        NOTE: Logic for the boundaries not implemented yet
        """
        accepted_moves=0
        for _ in range(steps):
            p_id = self.get_random_pids(ptype=ptype,size=1)
            old_position = self.system.part.by_id().pos
            E_pot_old = self.system.analysis.potential_energy()
            new_position = self.get_random_position_in_box()
            self.system.part.by_id.pos=new_position
            if check_exclusion_range:
                self.system.part.by_id().pos=old_position
            E_pot_new = self.system.analysis.potential_energy()
            bf = self.calculate_acceptance_probability(E_pot_new-E_pot_old)
            if self.rng.uniform() < bf:  # accept trial move
                accepted_moves+=1
            else:
                self.system.part.by_id().pos=old_position
        return

    def calculate_acceptance_probability(self,potential_energy_diff):
        return np.exp(-potential_energy_diff/self.kT)


if args.method == "widom":
    # System parameters
    cs_bulk = 0.1
    N0 = 70
    box_l = (N0 / cs_bulk)**(1.0 / 3.0)
    seed=23
    # Integration parameters
    system = espressomd.System(box_l=[box_l, box_l, box_l])
    np.random.seed(seed=42)
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    temperature = 1.0

    for i in range(N0):
        system.part.add(pos=np.random.random(3) * system.box_l, type=1, q=-1)
    for i in range(N0, 2 * N0):
        system.part.add(pos=np.random.random(3) * system.box_l, type=2, q=1)

    wca_eps = 1.0
    wca_sig = 1.0
    types = [0, 1, 2]
    for type_1 in types:
        for type_2 in types:
            system.non_bonded_inter[type_1, type_2].wca.set_params(
                epsilon=wca_eps, sigma=wca_sig)

    p3m = espressomd.electrostatics.P3M(prefactor=2.0, accuracy=1e-3)
    system.electrostatics.solver = p3m
    p3m_params = p3m.get_params()

    # Warmup
    #############################################################
    # warmup integration (steepest descent)
    warm_steps = 20
    warm_n_times = 20
    min_dist = 0.9 * wca_sig

    # minimize energy using min_dist as the convergence criterion
    system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                        max_displacement=0.01)
    i = 0
    while system.analysis.min_dist() < min_dist and i < warm_n_times:
        system.integrator.run(warm_steps)
        i += 1

    system.integrator.set_vv()

    # activate thermostat
    system.thermostat.set_langevin(kT=temperature, gamma=1.0, seed=42)

    if args.mode == "core":
        widom = espressomd.reaction_methods.WidomInsertion(
        kT=temperature, seed=seed)
    elif args.mode == "python":
        widom = WidomInsertion(kT=temperature, seed=seed, system=system)

    # add insertion reaction
    insertion_reaction_id = 0
    widom.add_reaction(reactant_types=[],
                    reactant_coefficients=[], product_types=[1, 2],
                    product_coefficients=[1, 1], default_charges={1: -1, 2: +1})
    system.setup_type_map(type_list=[0, 1, 2])


    # Set the hidden particle type to the lowest possible number to speed
    # up the simulation
    widom.set_non_interacting_type(type=max(types) + 1)

    particle_insertion_potential_energy_samples = []

    n_iterations = 50
    n_samples_per_iteration = 100

    for i in range(n_iterations):
        for _ in range(n_samples_per_iteration):
            particle_insertion_potential_energy_samples.append(
                widom.calculate_particle_insertion_potential_energy(reaction_id=insertion_reaction_id))
        system.integrator.run(steps=500)

    mu_ex_mean, mu_ex_Delta = widom.calculate_excess_chemical_potential(
        particle_insertion_potential_energy_samples=particle_insertion_potential_energy_samples)

    print(
        f"excess chemical potential for an ion pair {mu_ex_mean:.4g} +/- {mu_ex_Delta:.4g}")    
    exit()



# System parameters
#############################################################
box_l = 35
system = espressomd.System(box_l=[box_l] * 3)
SEED = 42
N_acid = 50
N_chain = 5
sigma = 1.
epsilon = 1.
system.time_step = 0.01
system.cell_system.skin = 1.0
N_steps_MD = 1000
N_steps_MC = 50

# Reaction parameters
#############################################################
pKa = 7.
pH = 7.25
kT = 1.
friction = 1.
types = {
    "HA": 0,
    "A-": 1,
    "H+": 2,
}
charges = {
    "HA": 0.,
    "A-": -1.,
    "H+": +1.,
}
params = {
    "kT": kT,
    "exclusion_range": sigma,
    "seed": SEED,
}
if args.method == "cph":
    params["constant_pH"] = pH

# Setup
#############################################################
np.random.seed(seed=SEED)

positions = espressomd.polymer.linear_polymer_positions(
    n_polymers=N_chain, beads_per_chain=N_acid // N_chain,
    bond_length=sigma + 0.1, seed=SEED)

bond = espressomd.interactions.HarmonicBond(k=10., r_0=sigma + 0.1)
system.bonded_inter.add(bond)
for polymer_pos in positions:
    bond_partner = None
    for pos in polymer_pos:
        p = system.part.add(pos=pos, type=types["A-"], q=charges["A-"])
        if bond_partner:
            p.add_bond((bond, bond_partner))
        bond_partner = p

for _ in range(N_acid):
    pos = np.random.random(3) * system.box_l
    system.part.add(pos=pos, type=types["H+"], q=charges["H+"])

for type_pair in itertools.combinations_with_replacement(types.values(), 2):
    system.non_bonded_inter[type_pair[0], type_pair[1]].wca.set_params(
        epsilon=epsilon, sigma=sigma)

p3m = espressomd.electrostatics.P3M(
    prefactor=2., accuracy=1e-2, mesh=8, cao=3, verbose=False)
dh = espressomd.electrostatics.DH(
    prefactor=2., kappa=0., r_cut=0.2 * box_l)

# energy minimize the system
system.integrator.set_steepest_descent(
    f_max=0., gamma=0.1, max_displacement=0.1 * sigma)
system.integrator.run(200)
system.electrostatics.solver = p3m
system.integrator.run(1000)

# thermalize the system
system.integrator.set_vv()
system.thermostat.set_langevin(kT=kT, gamma=friction, seed=SEED)
system.integrator.run(1000)
system.electrostatics.solver = dh

if args.mode == "core":
    if args.method == "cph":
        ReactionMethod = espressomd.reaction_methods.ConstantpHEnsemble
    elif args.method == "re":
        ReactionMethod = espressomd.reaction_methods.ReactionEnsemble
elif args.mode == "python":
    if args.method == "cph":
        ReactionMethod = ConstantpHEnsemble
    elif args.method == "re":
        ReactionMethod = ReactionEnsemble
    params["system"] = system

# set up reaction method
RE = ReactionMethod(**params)
RE.set_non_interacting_type(type=max(types.values()) + 1)
if args.method == "cph":
    RE.add_reaction(
        gamma=10**-pKa,
        reactant_types=[types["HA"]],
        product_types=[types["A-"], types["H+"]],
        default_charges={types[name]: charges[name] for name in types.keys()})
elif args.method == "re":
    RE.add_reaction(
        gamma=1e-3,
        reactant_types=[types["HA"]],
        reactant_coefficients=[1],
        product_types=[types["A-"], types["H+"]],
        product_coefficients=[1, 1],
        default_charges={types[name]: charges[name] for name in types.keys()})
reaction = RE.reactions[0]
system.setup_type_map(type_list=list(types.values()))

# equilibrate the polyelectrolyte chains
for i in range(5):
    RE.reaction(steps=10 * N_steps_MC)
    system.integrator.run(N_steps_MD)


@profile
def sample_alpha(length):
    alpha_list = []
    for _ in tqdm.tqdm(range(length)):
        system.integrator.run(steps=N_steps_MD)
        RE.reaction(steps=N_steps_MC)
        alpha = system.number_of_particles(type=types["A-"]) / N_acid
        alpha_list.append(alpha)
    return alpha_list


sample_size = 100
tick = time.time()
alphas = sample_alpha(sample_size)
tock = time.time()

alpha_avg = np.mean(alphas)
alpha_err = 1.96 * np.sqrt(np.var(alphas) / len(alphas))
acceptance_rate = reaction.get_acceptance_rate()
print(f"acceptance rate = {100. * acceptance_rate:.0f}%")
print(f"alpha = {alpha_avg:.2f} +/- {alpha_err:.2f}")
print(f"runtime = {tock - tick:.2f}s")
