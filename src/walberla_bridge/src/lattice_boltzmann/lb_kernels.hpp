/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#pragma once

#include <walberla_bridge/Architecture.hpp>

#include "generated_kernels/DynamicUBBDoublePrecision.h"
#include "generated_kernels/DynamicUBBSinglePrecision.h"
#include "generated_kernels/DynamicUBBColorGradientDoublePrecision.h"
#include "generated_kernels/DynamicUBBColorGradientSinglePrecision.h"
#include "generated_kernels/FieldAccessorsDoublePrecision.h"
#include "generated_kernels/FieldAccessorsSinglePrecision.h"
#include "generated_kernels/InitialPDFsSetterDoublePrecision.h"
#include "generated_kernels/InitialPDFsSetterSinglePrecision.h"
#include "generated_kernels/PackInfoPdfDoublePrecision.h"
#include "generated_kernels/PackInfoPdfSinglePrecision.h"
#include "generated_kernels/PackInfoVecDoublePrecision.h"
#include "generated_kernels/PackInfoVecSinglePrecision.h"
#include "generated_kernels/UpdateVelFromPDFDoublePrecision.h"
#include "generated_kernels/UpdateVelFromPDFSinglePrecision.h"
#include "generated_kernels/ColorGradientInitialPDFsSetterDoublePrecision.h"
#include "generated_kernels/ColorGradientInitialPDFsSetterSinglePrecision.h"

#ifdef __AVX2__
#include "generated_kernels/StreamCollideSweepLeesEdwardsDoublePrecisionAVX.h"
#include "generated_kernels/StreamCollideSweepLeesEdwardsSinglePrecisionAVX.h"
#include "generated_kernels/StreamCollideSweepThermalizedDoublePrecisionAVX.h"
#include "generated_kernels/StreamCollideSweepThermalizedSinglePrecisionAVX.h"
#include "generated_kernels/ColorGradientCollideSweepDoublePrecisionAVX.h"
#include "generated_kernels/ColorGradientCollideSweepSinglePrecisionAVX.h"
#include "generated_kernels/ColorGradientStreamSweepDoublePrecisionAVX.h"
#include "generated_kernels/ColorGradientStreamSweepSinglePrecisionAVX.h"
#else
#include "generated_kernels/StreamCollideSweepLeesEdwardsDoublePrecision.h"
#include "generated_kernels/StreamCollideSweepLeesEdwardsSinglePrecision.h"
#include "generated_kernels/StreamCollideSweepThermalizedDoublePrecision.h"
#include "generated_kernels/StreamCollideSweepThermalizedSinglePrecision.h"
#include "generated_kernels/ColorGradientCollideSweepDoublePrecision.h"
#include "generated_kernels/ColorGradientCollideSweepSinglePrecision.h"
#include "generated_kernels/ColorGradientStreamSweepDoublePrecision.h"
#include "generated_kernels/ColorGradientStreamSweepSinglePrecision.h"
#endif

namespace walberla {
namespace detail {

template <typename FT = double, lbmpy::Arch AT = lbmpy::Arch::CPU>
struct KernelTrait {
#ifdef __AVX2__
  using StreamCollisionModelThermalized =
      pystencils::StreamCollideSweepThermalizedDoublePrecisionAVX;
  using StreamCollisionModelLeesEdwards =
      pystencils::StreamCollideSweepLeesEdwardsDoublePrecisionAVX;
  using StreamModelTwoComponent =
      pystencils::ColorGradientStreamSweepDoublePrecisionAVX;
  using CollisionModelTwoComponent =
      pystencils::ColorGradientCollideSweepDoublePrecisionAVX;
#else
  using StreamCollisionModelThermalized =
      pystencils::StreamCollideSweepThermalizedDoublePrecision;
  using StreamCollisionModelLeesEdwards =
      pystencils::StreamCollideSweepLeesEdwardsDoublePrecision;
  using StreamModelTwoComponent =
      pystencils::ColorGradientStreamSweepDoublePrecision;
  using CollisionModelTwoComponent =
      pystencils::ColorGradientCollideSweepDoublePrecision;
#endif
  using InitialPDFsSetter = pystencils::InitialPDFsSetterDoublePrecision;
  using UpdateVelFromPDF = pystencils::UpdateVelFromPDFDoublePrecision;
  using PackInfoPdf = pystencils::PackInfoPdfDoublePrecision;
  using PackInfoVec = pystencils::PackInfoVecDoublePrecision;
  using InitialPDFsSetterTwoComponent = pystencils::ColorGradientInitialPDFsSetterDoublePrecision;
  using DynamicUBB = lbm::DynamicUBBDoublePrecision;
  using DynamicUBBTwoComponent = lbm::DynamicUBBColorGradientDoublePrecision;
};

template <> struct KernelTrait<float, lbmpy::Arch::CPU> {
#ifdef __AVX2__
  using StreamCollisionModelThermalized =
      pystencils::StreamCollideSweepThermalizedSinglePrecisionAVX;
  using StreamCollisionModelLeesEdwards =
      pystencils::StreamCollideSweepLeesEdwardsSinglePrecisionAVX;
  using StreamModelTwoComponent =
      pystencils::ColorGradientStreamSweepSinglePrecisionAVX;
  using CollisionModelTwoComponent =
      pystencils::ColorGradientCollideSweepSinglePrecisionAVX;
#else
  using StreamCollisionModelThermalized =
      pystencils::StreamCollideSweepThermalizedSinglePrecision;
  using StreamCollisionModelLeesEdwards =
      pystencils::StreamCollideSweepLeesEdwardsSinglePrecision;
  using StreamModelTwoComponent =
      pystencils::ColorGradientStreamSweepSinglePrecision;
  using CollisionModelTwoComponent =
      pystencils::ColorGradientCollideSweepSinglePrecision;
#endif
  using InitialPDFsSetter = pystencils::InitialPDFsSetterSinglePrecision;
  using UpdateVelFromPDF = pystencils::UpdateVelFromPDFSinglePrecision;
  using PackInfoPdf = pystencils::PackInfoPdfSinglePrecision;
  using PackInfoVec = pystencils::PackInfoVecSinglePrecision;
  using InitialPDFsSetterTwoComponent = pystencils::ColorGradientInitialPDFsSetterSinglePrecision;
  using DynamicUBB = lbm::DynamicUBBSinglePrecision;
  using DynamicUBBTwoComponent = lbm::DynamicUBBColorGradientSinglePrecision;
};

} // namespace detail
} // namespace walberla
