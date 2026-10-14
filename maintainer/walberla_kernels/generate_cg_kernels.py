#
# Copyright (C) 2020-2024 The ESPResSo project
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

import re
import argparse
import packaging.specifiers

import pystencils as ps

import lbmpy
import lbmpy.boundaries

kernel_codes = "cg_init cg_stream cg_collide cg_boundary".split()
parser = argparse.ArgumentParser(description="Generate the waLBerla kernels.")
parser.add_argument("--single-precision", action="store_true", required=False,
                    help="Use single-precision")
parser.add_argument("--gpu", action="store_true")
parser.add_argument("--kernels", nargs="+", type=str, default="all",
                    choices=["all"] + kernel_codes,
                    help="Which kernels to generate")
args = parser.parse_args()

# Make sure we have the correct versions of the required dependencies
for module, requirement in [(ps, "==1.4.0"), (lbmpy, "==1.4.0")]:
    assert packaging.specifiers.SpecifierSet(requirement).contains(module.__version__), \
        f"{module.__name__} version {module.__version__} " \
        f"doesn't match requirement {requirement}"

import pystencils_walberla
import pystencils_espresso
import lbmpy.stencils
import lbmpy.enums

import color_gradient
import code_generation_context
import custom_additional_extensions
from kernel_generation_utils import paramlist, get_ext_header, get_ext_source, patch_openmp_kernels

if args.gpu:
    target = ps.Target.GPU
else:
    target = ps.Target.CPU
if args.kernels == "all" or args.kernels == ["all"]:
    args.kernels = kernel_codes

# vectorization parameters
parameters = {}
if target == ps.Target.GPU:
    default_key = "GPU"
    parameters["GPU"] = ({"target": target}, "CUDA")
else:
    default_key = "CPU"
    cpu_vectorize_info = {
        "instruction_set": "avx",
        "assume_inner_stride_one": True,
        "assume_aligned": True,
        "assume_sufficient_line_padding": False}
    parameters["CPU"] = ({"target": target,
                          "cpu_openmp": True}, "")
    parameters["CPU_linear"] = ({"target": target}, "")
    parameters["AVX"] = ({"target": target,
                          "cpu_openmp": True,
                          "cpu_vectorize_info": cpu_vectorize_info}, "AVX")

# global parameters
stencil = lbmpy.stencils.LBStencil(lbmpy.enums.Stencil.D3Q19)
np2cpp_t = pystencils_espresso.numpy_types_to_cpp_types


def generate_cg_init_kernels(ctx, cg_fields, cg_methods):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    init = color_gradient.create_init_operator(cg_fields, cg_methods)
    for params, target_suffix in paramlist(parameters, (default_key,)):
        stem = f"ColorGradientInitialPDFsSetter{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(ctx, stem, init, **params)
        ctx.patch_file(stem, get_ext_source(target_suffix), patch_openmp_kernels)


def generate_cg_collide_kernels(ctx, cg_fields, cg_methods, cg_configs, cg_opts, cg_fluctuations):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    # Recolouring depends only on the direction of the colour gradient, not on
    # its magnitude. Where the phase field is flat to round-off, the gradient
    # and its projections are round-off residues computed by different sums,
    # so the cosines are not bounded by 1 and can be arbitrarily large.
    # Recolouring and surface tension are therefore switched off below a
    # threshold well above round-off and far below any interface gradient.
    minimum_color_gradient = 1e-12 if ctx.double_accuracy else 1e-5
    collide = color_gradient.create_collide_perturb_recolor_operator(
        cg_fields, cg_methods, cg_configs, cg_opts, cg_fluctuations,
        minimum_color_gradient=minimum_color_gradient)
    block_offsets = cg_fluctuations[0]["block_offsets"]
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"ColorGradientCollideSweep{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(
            ctx, stem, collide, block_offset=block_offsets, **params)
        ctx.patch_file(stem, get_ext_source(target_suffix), patch_openmp_kernels)


def generate_cg_stream_kernels(ctx, cg_fields, cg_methods, cg_configs, cg_opts):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    stream = color_gradient.create_stream_operator(cg_fields, cg_methods, cg_configs, cg_opts)
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"ColorGradientStreamSweep{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(
            ctx, stem, stream,
            field_swaps=[
                (cg_fields["pdfs_a"], cg_fields["pdfs_a_tmp"]),
                (cg_fields["pdfs_b"], cg_fields["pdfs_b_tmp"]),
            ],
            **params)
        ctx.patch_file(stem, get_ext_source(target_suffix), patch_openmp_kernels)


def generate_cg_boundary_kernels(ctx, cg_methods, data_type):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    # Both color components share the same stencil and the same
    # (compressible, zero_centered=False) conserved-quantity convention --
    # they only differ in their force field and relaxation rate, neither of
    # which the UBB boundary condition depends on. One generated kernel is
    # therefore reused for both components' BoundaryModel instances, unlike
    # the single-component DynamicUBB kernel (generated with
    # zero_centered=True) which cannot be reused here: applying it to the
    # raw (non-zero-centered) CG PDF fields reconstructs the local density
    # as delta_rho + 1.0 instead of delta_rho, corrupting the velocity
    # correction term for any boundary with a nonzero slip velocity.
    method = cg_methods[0]
    ubb_dynamic = lbmpy.boundaries.UBB(
        lambda *args: None, dim=3, data_type=data_type, calculate_force_on_boundary=True)
    ubb_data_handler = custom_additional_extensions.BounceBackSlipVelocityUBB(
        method.stencil, ubb_dynamic)

    # pylint: disable=unused-argument
    def patch_boundary_header(content, target_suffix):
        # replace real_t by actual floating-point type
        return content.replace("real_t", f"{np2cpp_t[data_type]}") \
                      .replace("real_c", f"{np2cpp_t[data_type]}_c")

    def patch_boundary_kernel(content, target_suffix):
        if target_suffix in ["CUDA"]:
            # replace preprocessor macros and pragmas
            push, pop = custom_additional_extensions.generate_device_preprocessor(
                "ubb_boundary", defines=("RESTRICT",))
            content = re.sub(r"#ifdef __GNUC__[\s\S]+?#endif(?=\n\n|\n//)", "", content)  # nopep8
            content = re.sub(r"#ifdef __CUDACC__[\s\S]+?#endif(?=\n\n|\n//)", push, content, 1)  # nopep8
            content = re.sub(r"#ifdef __CUDACC__[\s\S]+?#endif(?=\n\n|\n//)", pop, content, 1)  # nopep8
            assert push in content
            assert pop in content
        return content

    for _, target_suffix in paramlist(parameters, ("CPU", "GPU")):
        class_name = f"DynamicUBBColorGradient{precision_prefix}{target_suffix}"
        custom_additional_extensions.generate_lb_boundary(
            ctx, class_name, ubb_dynamic, method,
            additional_data_handler=ubb_data_handler,
            streaming_pattern="pull", target=target,
            template_file="templates/Boundary.tmpl.h")
        ctx.patch_file(class_name, get_ext_header(target_suffix),
                       patch_boundary_header, target_suffix)
        ctx.patch_file(class_name, get_ext_source(target_suffix),
                       patch_boundary_kernel, target_suffix)
        ctx.patch_file(class_name, get_ext_source(target_suffix),
                       patch_openmp_kernels)


with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = not args.single_precision
    if target == ps.Target.CPU:
        ctx.openmp = True
    if target == ps.Target.GPU:
        ctx.gpu = True
        ctx.cuda = True

    # codegen configuration
    config = pystencils_espresso.generate_config(
        ctx, parameters[default_key][0])
    data_type = "float64" if ctx.double_accuracy else "float32"
    
    # CG fields    
    cg_fields = color_gradient.generate_fields(stencil, data_type)

    # CG Method definition
    cg_methods = color_gradient.create_methods(stencil, cg_fields)
    
    cg_configs = color_gradient.create_configs(cg_fields, cg_methods)
    cg_opts = color_gradient.create_opts(cg_fields)
    cg_fluctuations = color_gradient.create_fluctuations(ctx.double_accuracy)

    if "cg_init" in args.kernels:
        generate_cg_init_kernels(ctx, cg_fields, cg_methods)
    if "cg_stream" in args.kernels:
        generate_cg_stream_kernels(ctx, cg_fields, cg_methods, cg_configs, cg_opts)
    if "cg_collide" in args.kernels:
        generate_cg_collide_kernels(ctx, cg_fields, cg_methods, cg_configs, cg_opts, cg_fluctuations)
    if "cg_boundary" in args.kernels:
        generate_cg_boundary_kernels(ctx, cg_methods, data_type)