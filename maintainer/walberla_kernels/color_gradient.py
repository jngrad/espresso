#
# Copyright (C) 2021-2023 The ESPResSo project
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

import sympy as sp
import pystencils as ps
from pystencils.sympyextensions import common_denominator

import lbmpy
import lbmpy.creationfunctions
import lbmpy.macroscopic_value_kernels
import lbmpy.forcemodels
import lbmpy.relaxationrates
import lbmpy.stencils
import lbmpy.moments

import numpy as np

import relaxation_rates
import pystencils_espresso

def generate_fields(stencil, data_type, field_layout='fzyx'):
    q = len(stencil)
    dim = len(stencil[0])

    fields = {}

    fields['velocity'] = ps.Field.create_generic(
        'velocity',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )
    fields['phasefield'] = ps.Field.create_generic(
        'phasefield',
        dim,
        data_type,
        index_dimensions=0,
        layout=field_layout,
    )

    # Symbols for component a
    fields['pdfs_a'] = ps.Field.create_generic(
        'pdfs_a',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['pdfs_a_tmp'] = ps.Field.create_generic(
        'pdfs_a_tmp',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['force_a'] = ps.Field.create_generic(
        'force_a',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )
    fields['rho_a'] = ps.Field.create_generic(
        'rho_a',
        dim,
        data_type,
        index_dimensions=0,
        layout=field_layout,
    )

    # Symbols for component b
    fields['pdfs_b'] = ps.Field.create_generic(
        'pdfs_b',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['pdfs_b_tmp'] = ps.Field.create_generic(
        'pdfs_b_tmp',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['force_b'] = ps.Field.create_generic(
        'force_b',
        dim,
        data_type,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )
    fields['rho_b'] = ps.Field.create_generic(
        'rho_b',
        dim,
        data_type,
        index_dimensions=0,
        layout=field_layout,
    )

    return fields

def create_methods(stencil, fields):
    # Both components collide with the same locally phase-interpolated relaxation rate
    omega_shear_a_sym = sp.Symbol("omega_shear_a")
    omega_shear_b_sym = sp.Symbol("omega_shear_b")
    omega_eff = get_interpolated_relaxation_rate(
        fields["phasefield"].center, omega_shear_a_sym, omega_shear_b_sym)
    omega_odd_eff = lbmpy.relaxationrates.relaxation_rate_from_magic_number(
        hydrodynamic_relaxation_rate=omega_eff)

    return (
        color_gradient_lb_method(stencil, fields["force_a"], omega_eff, omega_odd_eff),
        color_gradient_lb_method(stencil, fields["force_b"], omega_eff, omega_odd_eff),
    )
    
def create_configs(fields, methods):
    stencil = methods[0].stencil
    config_a = lbmpy.creationfunctions.LBMConfig(
        stencil=stencil,
        lb_method=methods[0],
        kernel_type="collide_only",
        streaming_pattern="pull",
        field_name=fields["pdfs_a"].name,
        density_input=fields["rho_a"],
        velocity_input=fields["velocity"],
    )
    config_b = lbmpy.creationfunctions.LBMConfig(
        stencil=stencil,
        lb_method=methods[1],
        kernel_type="collide_only",
        streaming_pattern="pull",
        field_name=fields["pdfs_b"].name,
        density_input=fields["rho_b"],
        velocity_input=fields["velocity"],
    )
    return config_a, config_b
    
def create_opts(fields):
    opt_a = lbmpy.LBMOptimisation(
        symbolic_field=fields["pdfs_a"],
        symbolic_temporary_field=fields["pdfs_a_tmp"],
        field_layout=fields["pdfs_a"].layout
        )
    opt_b = lbmpy.LBMOptimisation(
        symbolic_field=fields["pdfs_b"],
        symbolic_temporary_field=fields["pdfs_b_tmp"],
        field_layout=fields["pdfs_b"].layout
        )
    return opt_a, opt_b

def create_fluctuations(double_accuracy):
    """Noise parameters for the two per-color collision rules.

    Both components share a single ``seed``, one ``kT`` and the same Philox
    counter (time step and global lattice site). Their noise is kept
    independent by giving component b a disjoint range of Philox mode-index
    keys.
    """
    common = {
        "temperature": sp.symbols("kT"),
        "block_offsets": tuple(
            ps.TypedSymbol(f"block_offset_{i}", np.uint32) for i in range(3)),
        "seed": ps.TypedSymbol("seed", np.uint32),
    }
    fluctuating_a = {
        **common,
        "rng_node": pystencils_espresso.precision_rng[double_accuracy],
    }
    fluctuating_b = {
        **common,
        "rng_node": pystencils_espresso.precision_rng_shifted[double_accuracy],
    }
    return fluctuating_a, fluctuating_b

def get_rng_nodes(collection: ps.AssignmentCollection) -> list:
    return [s for s in collection.subexpressions
            if not isinstance(s, ps.Assignment)]

def sympy_cse_with_rng_nodes(ac: ps.AssignmentCollection) -> ps.AssignmentCollection:
    """ps.simp.sympy_cse() that tolerates RNG nodes in the subexpressions.
    """
    nodes = get_rng_nodes(ac)
    if not nodes:
        return ps.simp.sympy_cse(ac)
    stripped = ac.copy(
        ac.main_assignments,
        [s for s in ac.subexpressions if isinstance(s, ps.Assignment)],
    )
    simplified = ps.simp.sympy_cse(stripped)
    return simplified.copy(
        simplified.main_assignments, nodes + simplified.subexpressions)

def append_to_all_non_field_symbols(
    collection: ps.AssignmentCollection, appendix: str
    ) -> ps.AssignmentCollection:

    protected = {sym for node in get_rng_nodes(collection)
                 for sym in node.symbols_defined}

    substitutions = dict(
        [
            (val, sp.Symbol(f"{val.name}_{appendix}"))
            for val in collection.bound_symbols
            if not isinstance(val, ps.Field.Access) and val not in protected
            and not val.name.startswith("xi_")
        ]
    )
    xi_substitutions = dict(
        [
            (val, sp.Symbol(f"xi{appendix}_{val.name[len('xi_'):]}"))
            for val in collection.bound_symbols
            if not isinstance(val, ps.Field.Access) and val not in protected
            and val.name.startswith("xi_")
        ]
    )
    substitutions.update(xi_substitutions)
    return collection.new_with_substitutions(substitutions)

def move_main_assignment_to_subexpressions(ac, symbols):
    subexpressions = ac.subexpressions_dict
    main_assignments = ac.main_assignments_dict
    for sym in symbols:
        subexpressions[sym] = main_assignments.pop(sym)
    return ps.AssignmentCollection(
        subexpressions=subexpressions,
        main_assignments=main_assignments,
    )

def get_color_gradient(fields):
    # use second-order isotropic D3Q27 color gradient discretization from leclaire11a
    # where weights coincide with standard lattice-Boltzmann weights

    temp_stencil = lbmpy.stencils.LBStencil("D3Q27")

    # lattice speed of sound squared
    c_s_sq = sp.Rational(1, 3)

    # Convert weights to sympy Rationals to avoid float comparison issues
    raw_weights = lbmpy.get_weights(stencil=temp_stencil)
    fd_weights = list(
        map(
            lambda x: sp.nsimplify(x, rational=True) / c_s_sq,
            raw_weights,
        )
    )

    #compute colorgradient
    gradient = sum(
        (
            fields["phasefield"].center.get_shifted(*direction)
            * weight
            * sp.Matrix(direction)
            for direction, weight in zip(temp_stencil[1:], fd_weights[1:])
        ),
        fields["phasefield"].center.get_shifted(*temp_stencil[0])
        * fd_weights[0]
        * sp.Matrix(temp_stencil[0]),
    )
    return gradient

def get_interpolated_relaxation_rate(phasefield, omega_a, omega_b, delta=0.5):
    """
    Viscosity interpolation model taken from resi07a for the color-gradient model.
    The bulk value 'omega_a' is imposed for phi > 'delta' and 'omega_b' for phi < '-delta'.
    It is interpolating the harmonic mean for phi=0 piecewise with second-order polynomials towards bulk values at +/- delta.
    The constraints for the polynomaials are the bulk values at +/- delta, vanishing derivative at that point and the harmonic mean value
    at phi=0.
    """

    xi = 2 * omega_a * omega_b / (omega_a + omega_b)
    eta = 2 / delta * (omega_a - xi)
    kappa = -eta / (2 * delta)
    lam = 2 / delta * (xi - omega_b)
    nu = lam / (2 * delta)

    return sp.Piecewise(
        (omega_a, phasefield > delta),
        (omega_b, phasefield < -delta),
        (xi + eta * phasefield + kappa * phasefield**2, sp.And(phasefield > 0, phasefield <= delta)),
        (xi + lam * phasefield + nu * phasefield**2, True),
    )  # sp.And(phi >= -delta, phi < 0)

def get_weighted_orthogonal_moments(stencil, weights):
    """Reproduce the Gram-Schmidt orthogonalization that
    lbmpy.methods.creationfunctions.create_mrt_orthogonal
    applies to the raw default moment set
    """
    x, y, z = lbmpy.moments.MOMENT_SYMBOLS
    moments = lbmpy.moments.get_default_moment_set_for_stencil(stencil)
    diagonal_viscous_moments = [x**2 + y**2 + z**2, x**2, y**2 - z**2]
    for i, d in enumerate(lbmpy.moments.MOMENT_SYMBOLS[:stencil.D]):
        if d**2 in moments:
            moments[moments.index(d**2)] = diagonal_viscous_moments[i]
    orthogonal_moments = lbmpy.moments.gram_schmidt(moments, stencil, weights)
    return [e * common_denominator(e) for e in orthogonal_moments]

def color_gradient_lb_method(stencil, force_field: ps.Field, omega_eff, omega_odd_eff):

        def s(*args):
            for r in ps.sympyextensions.multidimensional_sum(*args, dim=len(stencil[0])):
                yield r

        def get_eq_value_offset(alpha):
            # for D3Q19
            values = {0: alpha, 1: (1 - alpha) / 12, sp.sqrt(2): (1 - alpha) / 24}
            return tuple(map(lambda x: values[sp.Matrix(x).norm()], stencil.stencil_entries))

        # lattice sound speed
        c_s = sp.sqrt(sp.Rational(1, 3))

        force_model = lbmpy.forcemodels.Schiller(force_field.center_vector)
        weights = lbmpy.maxwellian_equilibrium.get_weights(stencil)
        kd = ps.sympyextensions.kronecker_delta

        # offset values phi for the altered equilibrium
        # choose alpha = w_0 such that latice sound speed is c_s_sq = 1/3
        phi = get_eq_value_offset(alpha = weights[0])

        rho = sp.Symbol("rho")
        v = sp.symbols(f"u_:{stencil.D}")

        # calculate altered equilibrium
        equilibrium = []
        for d, w, phi_i in zip(stencil, weights, phi):
            result = sum(d[i] * v[i] for i, in s(1)) / (c_s**2)
            result += sum(
                v[i] * v[j] * (d[i] * d[j] - c_s**2 * kd(i, j)) for i, j in s(2)
            ) / (2 * c_s**4)
            result = rho * (phi_i + w * result)
            equilibrium.append(result)

        equilibrium = lbmpy.equilibrium.GenericDiscreteEquilibrium(
            stencil=stencil,
            equilibrium_pdfs=equilibrium,
            zeroth_order_moment_symbol=rho,
            first_order_moment_symbols=v,
            deviation_only=False,
        )

        # Use weighted-orthogonal moments -- required by lbmpy's fluctuating LB
        moments = get_weighted_orthogonal_moments(stencil, weights)

        # Build rr_dict from espresso's rr_getter function, but route every
        # non-conserved moment onto the effective relaxation rates
        # (omega_eff, omega_odd_eff) shared between both components 
        rr_dict = {}
        for m in moments:
            rr = relaxation_rates.rr_getter((m,))[0]
            if rr == 0:
                rr_dict[m] = 0
            elif rr.name == "omega_odd":
                rr_dict[m] = omega_odd_eff
            else:
                rr_dict[m] = omega_eff

        cqc = lbmpy.methods.DensityVelocityComputation(
            stencil=stencil, compressible=True, zero_centered=False, force_model=force_model
        )

        # create the altered lb method
        lb_method = lbmpy.methods.creationfunctions.create_from_equilibrium(
            stencil=stencil,
            equilibrium=equilibrium,
            conserved_quantity_computation=cqc,
            moment_to_relaxation_rate_dict=rr_dict,
            force_model=force_model,
        )
        lb_method.override_weights(weights)

        return lb_method

def create_collision_operator(configs, opts, fluctuations):
        
    # create collision operator for component a
    collision_rule_a = lbmpy.create_lb_collision_rule(
        lbm_config=configs[0], lbm_optimisation=opts[0],
        fluctuating=fluctuations[0]
    )
    accessor_a = lbmpy.fieldaccess.CollideOnlyInplaceAccessor
    collide_a = lbmpy.updatekernels.create_lbm_kernel(
        collision_rule=collision_rule_a,
        src_field=opts[0].symbolic_field,
        dst_field=opts[0].symbolic_field,
        accessor=accessor_a,
    )
    # create collision operator for component b
    collision_rule_b = lbmpy.create_lb_collision_rule(
        lbm_config=configs[1], lbm_optimisation=opts[1],
        fluctuating=fluctuations[1]
    )
    accessor_b = lbmpy.fieldaccess.CollideOnlyInplaceAccessor
    collide_b = lbmpy.updatekernels.create_lbm_kernel(
        collision_rule=collision_rule_b,
        src_field=opts[1].symbolic_field,
        dst_field=opts[1].symbolic_field,
        accessor=accessor_b,
    )

    #merge
    collide_a = append_to_all_non_field_symbols(collide_a, "a")
    collide_b = append_to_all_non_field_symbols(collide_b, "b")

    return ps.AssignmentCollection(
        main_assignments=collide_a.main_assignments + collide_b.main_assignments,
        subexpressions=collide_a.subexpressions + collide_b.subexpressions,
    )

def single_perturbation_operator(fields, method, config, opt, minimum_color_gradient):
    stencil = method.stencil

    pdfs_dst = lbmpy.macroscopic_value_kernels.get_field_accesses(
        lb_method=method,
        pdfs=opt.symbolic_field,
        streaming_pattern=config.streaming_pattern,
        previous_timestep=config.timestep,
        pre_collision_pdfs=False,
    )

    def get_b_value():
        """Implementation of eq. 17 in leclaire11a/leclaire17b."""
        values = {
            0: sp.Rational(-2, 9),
            1: sp.Rational(1, 54),
            sp.sqrt(2): sp.Rational(1, 27),
        }
        return tuple(map(lambda x: values[sp.Matrix(x).norm()], stencil.stencil_entries))

    omega_effective = get_interpolated_relaxation_rate(
        fields["phasefield"].center, sp.Symbol("omega_shear_a"), sp.Symbol("omega_shear_b"))

    A = sp.Rational(9, 4) * omega_effective * sp.Symbol("sigma")
    b = get_b_value()

    color_gradient = get_color_gradient(fields)
    f_norm = color_gradient.norm()

    # Create perturbation assignments
    result = []
    for dst, d, w, b_i in zip(pdfs_dst, method.stencil, method.weights, b):
        offset = sp.Piecewise(
            (
                A / 2
                * f_norm
                * (w * color_gradient.dot(sp.Matrix(d)) ** 2 / f_norm**2 - b_i),
                f_norm > minimum_color_gradient,
            ),
            (0, True),
        )

        result.append(
            ps.Assignment(
                dst,
                dst + offset,
            )
        )

    return ps.AssignmentCollection(result)

def perturbation_operator(fields, methods, configs, opts,
                          minimum_color_gradient: float = 0.0):
        perturbation_operator_a = single_perturbation_operator(
            fields, methods[0], configs[0], opts[0], minimum_color_gradient)
        perturbation_operator_b = single_perturbation_operator(
            fields, methods[1], configs[1], opts[1], minimum_color_gradient)
        perturbation = perturbation_operator_a.new_merged(perturbation_operator_b)

        return perturbation

def recoloring_operator(fields, methods, configs, opts, beta: float, minimum_color_gradient: float):
    stencil = methods[0].stencil
    from math import prod

    pdfs_dst_a = lbmpy.macroscopic_value_kernels.get_field_accesses(
        lb_method=methods[0],
        pdfs=opts[0].symbolic_field,
        streaming_pattern=configs[0].streaming_pattern,
        previous_timestep=configs[0].timestep,
        pre_collision_pdfs=False,
    )
    pdfs_dst_b = lbmpy.macroscopic_value_kernels.get_field_accesses(
        lb_method=methods[1],
        pdfs=opts[1].symbolic_field,
        streaming_pattern=configs[1].streaming_pattern,
        previous_timestep=configs[1].timestep,
        pre_collision_pdfs=False,
    )

    fluid_densities = (fields["rho_a"].center, fields["rho_b"].center)
    fluid_populations = (sp.Matrix(pdfs_dst_a), sp.Matrix(pdfs_dst_b))

    def get_population_equilibrium_for_zero_velocity(
        lb_method: lbmpy.methods.AbstractLbMethod, fluid_density: ps.Field.Access
    ):
        equilibrium = lb_method.get_equilibrium_terms()
        # set velocities to zero
        equilibrium = equilibrium.subs(
            dict([(x, 0) for x in lb_method.first_order_equilibrium_moment_symbols])
        )
        # replace density symbol with field
        equilibrium = equilibrium.subs(
            {lb_method.zeroth_order_equilibrium_moment_symbol: fluid_density}
        )

        return equilibrium

    # get equilibrium values for each population in each lb_method
    populations_equilibrium = tuple(
        map(
            lambda m, rho: get_population_equilibrium_for_zero_velocity(m, rho),
            [methods[0],methods[1]],
            fluid_densities,
        )
    )
    # calculate the sum of the equilibrium values between different lb_methods
    sum_populations_equilibrium = sum(
        populations_equilibrium[1:], populations_equilibrium[0]
    )

    # calculate the sum of the fluid population values between different lb_methods
    sum_populations = sum(fluid_populations[1:], fluid_populations[0])

    # calculate the total density
    sum_density = sum(fluid_densities)

    color_gradient = get_color_gradient(fields)
    color_gradient_norm = color_gradient.norm()

    # calculate the common term in eq. 10 containing the cos-value of the phasefield-projection
    common = [0]
    for d, f_i in zip(stencil[1:], sum_populations_equilibrium[1:]):
        d = sp.Matrix(d)
        cos_phi = sp.Piecewise(
            (
                color_gradient.dot(d) / (color_gradient_norm * d.norm()),
                color_gradient_norm > minimum_color_gradient,
            ),
            (0, True),
        )
        common.append(beta * prod(fluid_densities) / sum_density**2 * cos_phi * f_i)
    common = sp.Matrix(common)

    # calculate population update terms for each species
    terms = [
        rho_k / sum_density * sum_populations + (-1) ** k * common
        for k, rho_k in enumerate(fluid_densities)
    ]

    assignment_collections = tuple(
        map(
            lambda populations, value: ps.AssignmentCollection(
                [
                    ps.Assignment(populations, value),
                ]
            ),
            fluid_populations,
            terms,
        )
    )
    assert len(assignment_collections) == 2
    return ps.simp.sympy_cse(
        assignment_collections[0].new_merged(assignment_collections[1])
    )

def merged_collide_perturb_and_recoloring_operator(opts, collide, perturbation, recoloring):
    # append a suffix to the perturbation assignments to avoid name clashes
    collide = append_to_all_non_field_symbols(collide, "collide")
    perturbation = append_to_all_non_field_symbols(perturbation, "perturb")
    recoloring = append_to_all_non_field_symbols(recoloring, "recolor")

    tmp_populations_symbols_a = sp.symbols(
        f"tmp_a:{len(perturbation.main_assignments) // 2}"
    )
    symbols_a_iter = iter(tmp_populations_symbols_a)
    tmp_populations_symbols_b = sp.symbols(
        f"tmp_b:{len(perturbation.main_assignments) // 2}"
    )
    symbols_b_iter = iter(tmp_populations_symbols_b)

    # merge assignments
    merged = []
    for assignment in perturbation.main_assignments:
        assert (
            assignment.lhs in collide.defined_symbols
        ), f"could not find assignment in defined symbols of collision. {assignment.lhs}"

        assert (
            assignment.lhs in recoloring.defined_symbols
        ), f"could not find assignment in defined symbols of recoloring. {assignment.lhs}"

        if assignment.lhs.field == opts[0].symbolic_field:
            assign = next(symbols_a_iter)
        elif assignment.lhs.field == opts[1].symbolic_field:
            assign = next(symbols_b_iter)
        else:
            raise RuntimeError("Could not decide which color. No field matched")

        merged.append(
            ps.Assignment(
                assign,
                collide.main_assignments_dict[assignment.lhs]
                + assignment.rhs
                - assignment.lhs,
            )
        )

        # replace the reads in recoloring with the new assignments
        # so replace reads and writes and correct writes afterwards
        recoloring = recoloring.replace(assignment.lhs, assign)
        main_assignments_dict = recoloring.main_assignments_dict
        main_assignments_dict[assignment.lhs] = main_assignments_dict.pop(assign)
        recoloring = ps.AssignmentCollection(
            main_assignments=main_assignments_dict,
            subexpressions=recoloring.subexpressions,
        )

    merged = ps.AssignmentCollection(
        main_assignments=recoloring.main_assignments,
        subexpressions=collide.subexpressions
        + perturbation.subexpressions
        + merged
        + recoloring.subexpressions,
    )

    return sympy_cse_with_rng_nodes(merged)

def create_collide_perturb_recolor_operator(fields, methods, configs, opts, fluctuations,
                                           minimum_color_gradient: float):

    # create collide, perturbation and recoloring operator separately
    collide = create_collision_operator(configs, opts, fluctuations)
    perturbation = perturbation_operator(
        fields, methods, configs, opts, minimum_color_gradient=minimum_color_gradient)
    recoloring = recoloring_operator(
        fields, methods, configs, opts, beta=sp.Symbol("beta"),
        minimum_color_gradient=minimum_color_gradient)

    # merge collide, perturbation and recoloring operator and store it
    return merged_collide_perturb_and_recoloring_operator(opts, collide, perturbation, recoloring)

def create_init_operator(fields, methods):

        # Create init kernel for each component separately
        init_a = lbmpy.macroscopic_value_kernels.macroscopic_values_setter(
            methods[0],
            velocity=fields["velocity"].center_vector,
            pdfs=fields["pdfs_a"].center_vector,
            density=fields["rho_a"].center
            )
        init_b = lbmpy.macroscopic_value_kernels.macroscopic_values_setter(
            methods[1],
            velocity=fields["velocity"].center_vector,
            pdfs=fields["pdfs_b"].center_vector,
            density=fields["rho_b"].center
            )
        
        # Combine
        init_a = append_to_all_non_field_symbols(init_a, "a")
        init_b = append_to_all_non_field_symbols(init_b, "b")

        init = init_a.new_merged(init_b)
        
        # Add phasefield update
        rho_a = sp.Symbol(
            f"{methods[0].zeroth_order_equilibrium_moment_symbol.name}_a"
        )
        rho_b = sp.Symbol(
            f"{methods[1].zeroth_order_equilibrium_moment_symbol.name}_b"
        )
        sym_phasefield = sp.Symbol("phi")
        init = init.new_merged(
            ps.AssignmentCollection(
                subexpressions=[
                    ps.Assignment(sym_phasefield, (rho_a - rho_b) / (rho_a + rho_b))
                ],
                main_assignments=[
                    ps.Assignment(fields["phasefield"].center, sym_phasefield)
                ],
            )
        )

        return init

def create_stream_operator(fields, methods, configs, opts):
    stencil = methods[0].stencil
    density_symbol = sp.Symbol("rho_tmp")
    velocity_symbols = sp.symbols(f"vel_:{stencil.D}")
    
    accessor_a = lbmpy.advanced_streaming.utility.get_accessor(
    streaming_pattern=configs[0].streaming_pattern,
    timestep=configs[0].timestep,
    )
    accessor_b = lbmpy.advanced_streaming.utility.get_accessor(
        streaming_pattern=configs[1].streaming_pattern,
        timestep=configs[1].timestep,
    )

    stream_a = lbmpy.creationfunctions.create_stream_pull_with_output_kernel(
        lb_method=methods[0],
        src_field=opts[0].symbolic_field,
        dst_field=opts[0].symbolic_temporary_field,
        output={
            "density": density_symbol,
            "velocity": velocity_symbols,
        },
        accessor=accessor_a,
    )
    stream_b = lbmpy.creationfunctions.create_stream_pull_with_output_kernel(
        lb_method=methods[1],
        src_field=opts[1].symbolic_field,
        dst_field=opts[1].symbolic_temporary_field,
        output={
            "density": density_symbol,
            "velocity": velocity_symbols,
        },
        accessor=accessor_b,
    )

    stream_a = append_to_all_non_field_symbols(stream_a, "a")
    stream_b = append_to_all_non_field_symbols(stream_b, "b")

    density_symbol_a = sp.Symbol(f"{density_symbol.name}_a")
    density_symbol_b = sp.Symbol(f"{density_symbol.name}_b")
    velocity_symbols_a = sp.symbols(
        " ".join(map(lambda sym: f"{sym.name}_a", velocity_symbols))
    )
    velocity_symbols_b = sp.symbols(
        " ".join(map(lambda sym: f"{sym.name}_b", velocity_symbols))
    )

    stream_a_removed = move_main_assignment_to_subexpressions(
        stream_a, [density_symbol_a, *velocity_symbols_a]
    )
    stream_b_removed = move_main_assignment_to_subexpressions(
        stream_b, [density_symbol_b, *velocity_symbols_b]
    )

    stream = stream_a_removed.new_merged(stream_b_removed)

    subexpressions = stream.subexpressions_dict
    main_assignments = stream.main_assignments_dict

    # density-assignments
    main_assignments[fields["rho_a"].center] = density_symbol_a
    main_assignments[fields["rho_b"].center] = density_symbol_b

    # inverse density-assignments
    rho_total_inv = sp.Symbol("rho_total_inv")
    subexpressions[rho_total_inv] = 1 / (density_symbol_a + density_symbol_b)

    # velocity-assignments
    for vel_assignment, vel_sym_a, vel_sym_b in zip(
        fields["velocity"].center_vector,
        velocity_symbols_a,
        velocity_symbols_b,
    ):
        main_assignments[vel_assignment] = rho_total_inv * (
            density_symbol_a * vel_sym_a + density_symbol_b * vel_sym_b
        )

    # phasefield-assignment
    sym_phasefield = sp.Symbol("phi")
    subexpressions[sym_phasefield] = (
        density_symbol_a - density_symbol_b
    ) * rho_total_inv
    main_assignments[fields["phasefield"].center] = sym_phasefield

    stream = ps.AssignmentCollection(
        subexpressions=subexpressions,
        main_assignments=main_assignments,
    )
    return stream