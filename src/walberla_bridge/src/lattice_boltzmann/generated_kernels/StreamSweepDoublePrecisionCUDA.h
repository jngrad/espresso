//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file StreamSweepDoublePrecisionCUDA.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.2,
// lbmpy_walberla/pystencils_walberla from waLBerla commit
// 065ce5f311850371a97ac4766f47dbb5ca8424ba

#pragma once
#include "core/DataTypes.h"

#include "gpu/GPUField.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/SwapableCompare.h"

#include <functional>
#include <set>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace pystencils {

class StreamSweepDoublePrecisionCUDA {
public:
  StreamSweepDoublePrecisionCUDA(BlockDataID forceID_, BlockDataID pdfsID_,
                                 BlockDataID velocityID_)
      : forceID(forceID_), pdfsID(pdfsID_), velocityID(velocityID_){};

  ~StreamSweepDoublePrecisionCUDA() {
    for (auto p : cache_pdfs_) {
      delete p;
    }
  }

  void run(IBlock *block, gpuStream_t stream = nullptr);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block,
                         gpuStream_t stream = nullptr);

  void operator()(IBlock *block, gpuStream_t stream = nullptr) {
    run(block, stream);
  }

  static std::function<void(IBlock *)>
  getSweep(const shared_ptr<StreamSweepDoublePrecisionCUDA> &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *, gpuStream_t)> getSweepOnCellInterval(
      const shared_ptr<StreamSweepDoublePrecisionCUDA> &kernel,
      const shared_ptr<StructuredBlockStorage> &blocks,
      const CellInterval &globalCellInterval, cell_idx_t ghostLayers = 1) {
    return [kernel, blocks, globalCellInterval,
            ghostLayers](IBlock *b, gpuStream_t stream = nullptr) {
      kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b,
                                stream);
    };
  }

  std::function<void(IBlock *)> getSweep(gpuStream_t stream = nullptr) {
    return [this, stream](IBlock *b) { this->run(b, stream); };
  }

  std::function<void(IBlock *)>
  getSweepOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers = 1,
                         gpuStream_t stream = nullptr) {
    return [this, blocks, globalCellInterval, ghostLayers, stream](IBlock *b) {
      this->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b,
                              stream);
    };
  }

  BlockDataID forceID;
  BlockDataID pdfsID;
  BlockDataID velocityID;

private:
  std::set<gpu::GPUField<double> *,
           field::SwapableCompare<gpu::GPUField<double> *>>
      cache_pdfs_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif