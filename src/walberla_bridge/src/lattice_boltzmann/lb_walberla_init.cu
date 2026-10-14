/*
 * Copyright (C) 2022-2026 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#if defined(__NVCC__)
#define RESTRICT __restrict__
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 554 // no implicit or explicit cast
#else
#pragma push
#pragma diag_suppress 554 // no implicit or explicit cast
#endif
#endif

#include "LBWalberlaImplColorGradient.hpp"
#include "LBWalberlaImplSingleComponent.hpp"

#if defined(__NVCC__)
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic pop
#else
#pragma pop
#endif
#endif

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <gpu/DeviceSelectMPI.h>
#include <waLBerlaDefinitions.h>

#include <memory>
#include <vector>
#include <stdexcept>

std::shared_ptr<LBWalberlaBase>
new_lb_walberla_gpu(std::shared_ptr<LatticeWalberla> const &lattice,
                    std::vector<double> viscosity, double density,
                    bool single_precision) {
#if not defined(WALBERLA_BUILD_WITH_CUDA)
  throw std::runtime_error("waLBerla was compiled without CUDA support");
#else
  if (viscosity.size() == 1u) {
    if (single_precision) {
      return std::make_shared<
          walberla::LBWalberlaImplSingleComponent<float, lbmpy::Arch::GPU>>(
          lattice, viscosity, density, false);
    }
    return std::make_shared<
        walberla::LBWalberlaImplSingleComponent<double, lbmpy::Arch::GPU>>(
        lattice, viscosity, density, false);
  }
  // single precision doesn't conserve mass and momentum of the two-component
  // model well enough; the script interface rejects it before this point
  if (single_precision) {
    throw std::runtime_error(
        "The two-component LB only supports double precision on the GPU");
  }
  return std::make_shared<
      walberla::LBWalberlaImplColorGradient<double, lbmpy::Arch::GPU>>(
      lattice, viscosity, density, true);
#endif
}

std::shared_ptr<LBWalberlaBase>
new_lb_walberla_gpu(std::shared_ptr<LatticeWalberla> const &lattice,
                    double viscosity, double density,
                    bool single_precision) {
  return new_lb_walberla_gpu(lattice, std::vector<double>{viscosity}, density,
                             single_precision);
}

void set_device_id_per_rank() {
#if not defined(WALBERLA_BUILD_WITH_CUDA)
  throw std::runtime_error("waLBerla was compiled without CUDA support");
#else
  walberla::gpu::selectDeviceBasedOnMpiRank();
#endif
}
