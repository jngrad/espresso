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
//! \\author Martin Bauer <martin.bauer@fau.de>
//======================================================================================================================

#include <cmath>

#include "FluctuatingMRT_LatticeModel.h"
#include "core/DataTypes.h"
#include "core/Macros.h"
#include "lbm/field/PdfField.h"
#include "lbm/sweeps/Streaming.h"

#ifdef _MSC_VER
#pragma warning(disable : 4458)
#endif

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "philox_rand.h"

using namespace std;

namespace walberla {
namespace lbm {

namespace internal_kernel_streamCollide {
static FUNC_PREFIX void kernel_streamCollide(
    real_t *RESTRICT const _data_force, real_t *RESTRICT const _data_pdfs,
    real_t *RESTRICT _data_pdfs_tmp, int64_t const _size_force_0,
    int64_t const _size_force_1, int64_t const _size_force_2,
    int64_t const _stride_force_0, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t const _stride_pdfs_tmp_0, int64_t const _stride_pdfs_tmp_1,
    int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3,
    uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2,
    real_t omega_bulk, real_t omega_even, real_t omega_odd, real_t omega_shear,
    uint32_t seed, real_t temperature, uint32_t time_step) {
  const real_t xi_25 = -omega_bulk;
  const real_t xi_36 = -omega_shear;
  const real_t xi_37 = xi_36 + real_t(2.0);
  const real_t xi_38 = xi_37 * real_t(0.5);
  const real_t xi_43 = xi_37 * real_t(0.0833333333333333);
  const real_t xi_48 = xi_37 * real_t(0.166666666666667);
  const real_t xi_58 = xi_37 * real_t(0.25);
  const real_t xi_63 = xi_37 * real_t(0.0416666666666667);
  const real_t xi_90 = real_t(2.4494897427831779);
  const real_t xi_115 = omega_odd * real_t(0.25);
  const real_t xi_131 = omega_odd * real_t(0.0833333333333333);
  const real_t xi_196 = omega_shear * real_t(0.25);
  const real_t xi_211 = omega_odd * real_t(0.0416666666666667);
  const real_t xi_213 = omega_odd * real_t(0.125);
  const int64_t rr_0 = 0;
  const real_t xi_120 = real_c(rr_0) * real_t(0.166666666666667);
  const real_t xi_186 = real_c(rr_0) * real_t(0.0833333333333333);
  for (int ctr_2 = 1; ctr_2 < _size_force_2 - 1; ctr_2 += 1) {
    real_t *RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 14 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 18 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 11 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 15 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 12 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                         _stride_pdfs_2 + 5 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 16 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 17 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                        _stride_pdfs_2 + 6 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 13 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    real_t *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    real_t *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    real_t *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    real_t *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    real_t *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_30 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2;
    real_t *RESTRICT _data_pdfs_tmp_20_31 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_32 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_33 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_34 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_35 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_36 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_37 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_38 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_39 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_310 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_311 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_312 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_313 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_314 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_315 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_316 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_317 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_tmp_20_318 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3;
    for (int ctr_1 = 1; ctr_1 < _size_force_1 - 1; ctr_1 += 1) {
      real_t *RESTRICT _data_pdfs_2m1_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_314;
      real_t *RESTRICT _data_pdfs_21_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_318;
      real_t *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      real_t *RESTRICT _data_pdfs_2m1_311_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
      real_t *RESTRICT _data_pdfs_20_31_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
      real_t *RESTRICT _data_pdfs_21_315_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
      real_t *RESTRICT _data_pdfs_2m1_312_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
      real_t *RESTRICT _data_pdfs_2m1_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_35;
      real_t *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      real_t *RESTRICT _data_pdfs_20_39_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
      real_t *RESTRICT _data_pdfs_20_32_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
      real_t *RESTRICT _data_pdfs_21_316_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
      real_t *RESTRICT _data_pdfs_21_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_317;
      real_t *RESTRICT _data_pdfs_21_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_36;
      real_t *RESTRICT _data_pdfs_20_37_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
      real_t *RESTRICT _data_pdfs_2m1_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_313;
      real_t *RESTRICT _data_pdfs_20_310_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
      real_t *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      real_t *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      real_t *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      real_t *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      real_t *RESTRICT _data_pdfs_20_38_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
      real_t *RESTRICT _data_pdfs_tmp_20_30_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_30;
      real_t *RESTRICT _data_pdfs_tmp_20_31_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_31;
      real_t *RESTRICT _data_pdfs_tmp_20_32_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_32;
      real_t *RESTRICT _data_pdfs_tmp_20_33_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_33;
      real_t *RESTRICT _data_pdfs_tmp_20_34_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_34;
      real_t *RESTRICT _data_pdfs_tmp_20_35_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_35;
      real_t *RESTRICT _data_pdfs_tmp_20_36_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_36;
      real_t *RESTRICT _data_pdfs_tmp_20_37_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_37;
      real_t *RESTRICT _data_pdfs_tmp_20_38_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_38;
      real_t *RESTRICT _data_pdfs_tmp_20_39_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_39;
      real_t *RESTRICT _data_pdfs_tmp_20_310_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_310;
      real_t *RESTRICT _data_pdfs_tmp_20_311_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_311;
      real_t *RESTRICT _data_pdfs_tmp_20_312_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_312;
      real_t *RESTRICT _data_pdfs_tmp_20_313_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_313;
      real_t *RESTRICT _data_pdfs_tmp_20_314_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_314;
      real_t *RESTRICT _data_pdfs_tmp_20_315_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_315;
      real_t *RESTRICT _data_pdfs_tmp_20_316_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_316;
      real_t *RESTRICT _data_pdfs_tmp_20_317_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_317;
      real_t *RESTRICT _data_pdfs_tmp_20_318_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_318;
      for (int ctr_0 = 1; ctr_0 < _size_force_0 - 1; ctr_0 += 1) {

        float Dummy_299;
        float Dummy_300;
        float Dummy_301;
        float Dummy_302;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 3, seed, Dummy_299, Dummy_300,
                      Dummy_301, Dummy_302);

        float Dummy_295;
        float Dummy_296;
        float Dummy_297;
        float Dummy_298;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 2, seed, Dummy_295, Dummy_296,
                      Dummy_297, Dummy_298);

        float Dummy_291;
        float Dummy_292;
        float Dummy_293;
        float Dummy_294;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 1, seed, Dummy_291, Dummy_292,
                      Dummy_293, Dummy_294);

        float Dummy_287;
        float Dummy_288;
        float Dummy_289;
        float Dummy_290;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 0, seed, Dummy_287, Dummy_288,
                      Dummy_289, Dummy_290);

        const real_t xi_0 =
            _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] +
            _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_1 =
            xi_0 + _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_2 = _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] +
                            _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] +
                            _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        const real_t xi_3 = _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0] +
                            _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_4 =
            _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_5 = _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0] +
                            _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_6 =
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_8 =
            -_data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_9 =
            xi_8 -
            _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_10 =
            -_data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_11 =
            -_data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_12 =
            -_data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_13 = xi_10 + xi_11 + xi_12;
        const real_t xi_14 = -_data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_15 =
            -_data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_16 = xi_14 + xi_15;
        const real_t xi_17 = -_data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_18 = -_data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_19 = xi_17 + xi_18;
        const real_t xi_20 =
            -_data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_21 = xi_10 + xi_20;
        const real_t xi_22 = -_data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0];
        const real_t xi_23 = -_data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_24 = xi_17 + xi_22 + xi_23 +
                             _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        const real_t xi_42 =
            real_t(0.166666666666667) * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_50 =
            real_t(0.166666666666667) * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const real_t xi_54 =
            real_t(0.166666666666667) * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const real_t xi_57 =
            real_t(0.5) * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_61 =
            real_t(0.0833333333333333) * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const real_t xi_65 =
            real_t(0.0833333333333333) * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_75 =
            real_t(0.0833333333333333) * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const real_t xi_93 = -_data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_94 = xi_93 +
                             real_t(3.0) * _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] +
                             real_t(3.0) * _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_95 =
            omega_even *
            (xi_94 - real_t(3.0) * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] -
             real_t(3.0) * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] -
             real_t(3.0) * _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0] -
             real_t(3.0) * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0] +
             real_t(3.0) * _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] +
             real_t(3.0) * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0]);
        const real_t xi_96 =
            real_t(2.0) * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] +
            real_t(2.0) * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] +
            real_t(2.0) * _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0] +
            real_t(2.0) * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_97 =
            xi_96 +
            real_t(5.0) * _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            real_t(5.0) * _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_98 =
            omega_even *
            (xi_94 + xi_97 -
             real_t(2.0) * _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] -
             real_t(2.0) * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0] -
             real_t(5.0) *
                 _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] -
             real_t(5.0) *
                 _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] -
             real_t(5.0) * _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 +
                                         _stride_pdfs_0] -
             real_t(5.0) * _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 -
                                         _stride_pdfs_0]);
        const real_t xi_101 = -_data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        const real_t xi_102 = xi_101 + xi_18;
        const real_t xi_103 =
            -_data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_106 =
            -_data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_107 = xi_106 + xi_11 + xi_15 + xi_21;
        const real_t xi_109 =
            real_t(2.0) *
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_110 =
            real_t(2.0) *
            _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_111 =
            real_t(2.0) *
                _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            real_t(2.0) * _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_112 =
            omega_even *
            (xi_109 + xi_110 + xi_111 + xi_93 + xi_97 -
             real_t(4.0) * _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] -
             real_t(4.0) * _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0] -
             real_t(7.0) *
                 _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] -
             real_t(7.0) *
                 _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] -
             real_t(7.0) *
                 _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] -
             real_t(7.0) *
                 _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
             real_t(5.0) * _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] +
             real_t(5.0) * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0]);
        const real_t xi_113 =
            xi_101 + _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_114 = xi_113 + xi_14 + xi_22 +
                              _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] +
                              _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_116 = xi_114 * xi_115;
        const real_t xi_118 =
            xi_103 +
            _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_122 = Dummy_298 - real_t(0.5);
        const real_t xi_127 =
            real_t(2.0) * _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_128 =
            real_t(2.0) * _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_129 =
            -real_t(2.0) *
                _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] +
            real_t(2.0) * _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_130 = -xi_127 + xi_128 + xi_129 + xi_14 + xi_19 + xi_2;
        const real_t xi_132 = xi_130 * xi_131;
        const real_t xi_133 = Dummy_293 - real_t(0.5);
        const real_t xi_138 = Dummy_288 - real_t(0.5);
        const real_t xi_142 =
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_156 =
            xi_106 +
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_157 =
            xi_12 + xi_156 + xi_20 +
            _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] +
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_158 = xi_115 * xi_157;
        const real_t xi_159 = Dummy_296 - real_t(0.5);
        const real_t xi_161 = xi_1 + xi_127 - xi_128 + xi_129 + xi_13;
        const real_t xi_162 = xi_131 * xi_161;
        const real_t xi_163 = Dummy_295 - real_t(0.5);
        const real_t xi_168 = _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] +
                              _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_169 = xi_102 + xi_168 + xi_23 +
                              _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_170 = xi_115 * xi_169;
        const real_t xi_173 = Dummy_297 - real_t(0.5);
        const real_t xi_175 = -xi_109 - xi_110 + xi_111 + xi_24 + xi_3;
        const real_t xi_176 = xi_131 * xi_175;
        const real_t xi_177 = Dummy_294 - real_t(0.5);
        const real_t xi_184 = xi_112 * real_t(0.0138888888888889);
        const real_t xi_205 = xi_98 * -real_t(0.00714285714285714);
        const real_t xi_207 = xi_95 * real_t(0.025);
        const real_t xi_212 = xi_175 * xi_211;
        const real_t xi_214 = xi_169 * xi_213;
        const real_t xi_223 = xi_130 * xi_211;
        const real_t xi_224 = xi_114 * xi_213;
        const real_t xi_232 = xi_98 * real_t(0.0178571428571429);
        const real_t xi_238 = xi_157 * xi_213;
        const real_t xi_239 = xi_161 * xi_211;
        const real_t vel0Term =
            xi_1 +
            _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] +
            _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t vel1Term =
            xi_2 +
            _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t vel2Term =
            xi_3 +
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t rho = vel0Term + vel1Term + vel2Term + xi_4 + xi_5 + xi_6 +
                           _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_7 = 1 / (rho);
        const real_t xi_86 = rho * temperature;
        const real_t xi_87 =
            sqrt(xi_86 * (-((-omega_even + real_t(1.0)) * (-omega_even + real_t(1.0))) + real_t(1.0)));
        const real_t xi_88 = xi_87 * (Dummy_299 - real_t(0.5)) * real_t(3.7416573867739413);
        const real_t xi_89 = xi_87 * (Dummy_301 - real_t(0.5)) * real_t(5.4772255750516612);
        const real_t xi_91 =
            xi_90 * sqrt(xi_86 * (-((xi_25 + real_t(1.0)) * (xi_25 + real_t(1.0))) + real_t(1.0))) *
            (Dummy_292 - real_t(0.5));
        const real_t xi_92 = xi_87 * (Dummy_300 - real_t(0.5)) * real_t(8.3666002653407556);
        const real_t xi_123 =
            sqrt(xi_86 * (-((-omega_odd + real_t(1.0)) * (-omega_odd + real_t(1.0))) + real_t(1.0)));
        const real_t xi_124 = xi_123 * real_t(1.4142135623730951);
        const real_t xi_125 = xi_124 * real_t(0.5);
        const real_t xi_126 = xi_122 * xi_125;
        const real_t xi_134 = xi_123 * xi_90;
        const real_t xi_135 = xi_134 * real_t(0.166666666666667);
        const real_t xi_136 = xi_133 * xi_135;
        const real_t xi_137 = -xi_132 - xi_136;
        const real_t xi_139 =
            sqrt(xi_86 * (-((xi_36 + real_t(1.0)) * (xi_36 + real_t(1.0))) + real_t(1.0)));
        const real_t xi_140 = xi_139 * real_t(0.5);
        const real_t xi_141 = xi_138 * xi_140;
        const real_t xi_146 =
            xi_112 * -real_t(0.0198412698412698) + xi_88 * -real_t(0.119047619047619);
        const real_t xi_148 = xi_139 * (Dummy_287 - real_t(0.5)) * real_t(1.7320508075688772);
        const real_t xi_152 = xi_132 + xi_136;
        const real_t xi_160 = xi_125 * xi_159;
        const real_t xi_164 = xi_135 * xi_163;
        const real_t xi_165 = xi_162 + xi_164;
        const real_t xi_167 = -xi_162 - xi_164;
        const real_t xi_174 = xi_125 * xi_173;
        const real_t xi_178 = xi_135 * xi_177;
        const real_t xi_179 = -xi_176 - xi_178;
        const real_t xi_181 = xi_176 + xi_178;
        const real_t xi_182 = xi_138 * xi_139 * real_t(0.25);
        const real_t xi_185 = xi_88 * real_t(0.0833333333333333);
        const real_t xi_195 = xi_140 * (Dummy_289 - real_t(0.5));
        const real_t xi_204 = xi_140 * (Dummy_291 - real_t(0.5));
        const real_t xi_208 = xi_92 * -real_t(0.0142857142857143);
        const real_t xi_209 = xi_89 * real_t(0.05);
        const real_t xi_215 = xi_134 * real_t(0.0833333333333333);
        const real_t xi_216 = xi_177 * xi_215;
        const real_t xi_217 = xi_124 * real_t(0.25);
        const real_t xi_218 = xi_173 * xi_217;
        const real_t xi_220 =
            xi_112 * -real_t(0.00396825396825397) + xi_88 * -real_t(0.0238095238095238);
        const real_t xi_225 = xi_133 * xi_215;
        const real_t xi_226 = xi_122 * xi_217;
        const real_t xi_230 = -xi_182;
        const real_t xi_233 = xi_92 * real_t(0.0357142857142857);
        const real_t xi_235 = xi_140 * (Dummy_290 - real_t(0.5));
        const real_t xi_240 = xi_159 * xi_217;
        const real_t xi_241 = xi_163 * xi_215;
        const real_t u_0 = xi_7 * (vel0Term + xi_13 + xi_9);
        const real_t xi_26 =
            u_0 * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const real_t xi_27 = xi_26 * real_t(0.333333333333333);
        const real_t xi_33 = -xi_27;
        const real_t xi_99 = rho * (u_0 * u_0);
        const real_t xi_153 = rho * u_0;
        const real_t xi_154 =
            -vel0Term + xi_142 + xi_153 + xi_4 +
            _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_155 = xi_120 * xi_154;
        const real_t xi_191 = xi_154 * xi_186;
        const real_t u_1 =
            xi_7 *
            (vel1Term + xi_16 + xi_19 + xi_8 +
             _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]);
        const real_t xi_28 =
            u_1 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_29 = xi_28 * real_t(0.333333333333333);
        const real_t xi_34 = -xi_29;
        const real_t xi_56 = u_1 * real_t(0.5);
        const real_t xi_59 =
            xi_58 * (u_0 * xi_57 +
                     xi_56 * _data_force_20_30_10[_stride_force_0 * ctr_0]);
        const real_t xi_60 = -xi_59;
        const real_t xi_104 = rho * (u_1 * u_1);
        const real_t xi_105 = xi_103 + xi_104 + xi_9;
        const real_t xi_117 = rho * u_1;
        const real_t xi_119 =
            -vel1Term + xi_117 + xi_118 + xi_5 +
            _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_121 = xi_119 * xi_120;
        const real_t xi_187 = xi_119 * xi_186;
        const real_t xi_197 =
            xi_196 *
            (u_0 * xi_117 + xi_118 + xi_8 +
             _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]);
        const real_t xi_198 = -xi_195 - xi_197;
        const real_t xi_199 = xi_195 + xi_197;
        const real_t u_2 =
            xi_7 *
            (vel2Term + xi_21 + xi_24 +
             _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]);
        const real_t xi_30 =
            u_2 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const real_t xi_31 = xi_30 * real_t(0.333333333333333);
        const real_t xi_32 = (xi_25 +  real_t(2.0)) * (xi_27 + xi_29 + xi_31);
        const real_t xi_35 = xi_30 * real_t(0.666666666666667) + xi_33 + xi_34;
        const real_t xi_39 = -xi_31;
        const real_t xi_40 = xi_28 * real_t(0.666666666666667) + xi_33 + xi_39;
        const real_t xi_41 = xi_26 * real_t(0.666666666666667) + xi_34 + xi_39;
        const real_t xi_44 = xi_35 * xi_43;
        const real_t xi_45 = -xi_44;
        const real_t xi_46 = xi_41 * xi_43;
        const real_t xi_47 = -xi_46;
        const real_t xi_49 = xi_40 * xi_48 + xi_45 + xi_47;
        const real_t xi_51 = xi_40 * xi_43;
        const real_t xi_52 = -xi_51;
        const real_t xi_53 = xi_41 * xi_48 + xi_45 + xi_52;
        const real_t xi_55 = xi_35 * xi_48 + xi_47 + xi_52;
        const real_t xi_62 = xi_46 - xi_61;
        const real_t xi_64 = -xi_35 * xi_63;
        const real_t xi_66 = xi_32 * real_t(0.125);
        const real_t xi_67 = xi_51 + xi_66;
        const real_t xi_68 = xi_65 + xi_67;
        const real_t xi_69 = xi_64 + xi_68;
        const real_t xi_70 = xi_46 + xi_61;
        const real_t xi_71 = -xi_65 + xi_67;
        const real_t xi_72 = xi_64 + xi_71;
        const real_t xi_73 =
            xi_58 * (u_2 * xi_57 +
                     xi_56 * _data_force_20_32_10[_stride_force_0 * ctr_0]);
        const real_t xi_74 = -xi_41 * xi_63;
        const real_t xi_76 = xi_44 + xi_75;
        const real_t xi_77 = xi_74 + xi_76;
        const real_t xi_78 = -xi_73;
        const real_t xi_79 =
            xi_58 * (u_0 * real_t(0.5) * _data_force_20_32_10[_stride_force_0 * ctr_0] +
                     u_2 * real_t(0.5) * _data_force_20_30_10[_stride_force_0 * ctr_0]);
        const real_t xi_80 = -xi_79;
        const real_t xi_81 = -xi_40 * xi_63;
        const real_t xi_82 = xi_66 + xi_76 + xi_81;
        const real_t xi_83 = xi_44 - xi_75;
        const real_t xi_84 = xi_74 + xi_83;
        const real_t xi_85 = xi_66 + xi_81 + xi_83;
        const real_t xi_100 = rho * (u_2 * u_2);
        const real_t xi_108 =
            omega_bulk * (xi_100 + xi_102 + xi_105 + xi_107 + xi_17 + xi_22 +
                          xi_99 + _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0]);
        const real_t xi_143 = -xi_100 +
                              _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] +
                              _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_144 =
            omega_shear * (xi_0 + xi_105 + xi_142 + xi_143 + xi_16 -
                           _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0]);
        const real_t xi_145 = xi_144 * real_t(0.125);
        const real_t xi_147 =
            omega_shear *
            (xi_103 - xi_104 + xi_107 + xi_143 + xi_9 + xi_96 + xi_99 * real_t(2.0) -
             real_t(2.0) *
                 _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] -
             real_t(2.0) *
                 _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] +
             _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] +
             _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0]);
        const real_t xi_149 =
            xi_147 * -real_t(0.0416666666666667) + xi_148 * -real_t(0.166666666666667);
        const real_t xi_150 = xi_149 + xi_89 * -real_t(0.1) + xi_95 * -real_t(0.05);
        const real_t xi_151 = xi_141 + xi_145 + xi_146 + xi_150 +
                              xi_92 * real_t(0.0285714285714286) +
                              xi_98 * real_t(0.0142857142857143);
        const real_t xi_166 =
            xi_146 + xi_147 * real_t(0.0833333333333333) + xi_148 * real_t(0.333333333333333) +
            xi_92 * -real_t(0.0714285714285714) + xi_98 * -real_t(0.0357142857142857);
        const real_t xi_171 =
            rho * u_2 - vel2Term + xi_101 + xi_106 + xi_168 + xi_6 +
            _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_172 = xi_120 * xi_171;
        const real_t xi_180 = xi_112 * real_t(0.0158730158730159) - xi_141 - xi_145 +
                              xi_150 + xi_88 * real_t(0.0952380952380952) +
                              xi_92 * -real_t(0.0428571428571429) +
                              xi_98 * -real_t(0.0214285714285714);
        const real_t xi_183 = xi_144 * real_t(0.0625);
        const real_t xi_188 =
            xi_108 * real_t(0.0416666666666667) + xi_91 * real_t(0.0833333333333333);
        const real_t xi_189 = xi_187 + xi_188;
        const real_t xi_190 =
            xi_152 + xi_182 + xi_183 + xi_184 + xi_185 + xi_189;
        const real_t xi_192 =
            xi_147 * real_t(0.0208333333333333) + xi_148 * real_t(0.0833333333333333);
        const real_t xi_193 = -xi_191 + xi_192;
        const real_t xi_194 = xi_167 + xi_193;
        const real_t xi_200 = xi_191 + xi_192;
        const real_t xi_201 = xi_165 + xi_200;
        const real_t xi_202 = -xi_187 + xi_188;
        const real_t xi_203 =
            xi_137 + xi_182 + xi_183 + xi_184 + xi_185 + xi_202;
        const real_t xi_206 =
            xi_196 * (u_2 * xi_117 + xi_113 + xi_17 +
                      _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0]);
        const real_t xi_210 =
            xi_149 + xi_204 + xi_205 + xi_206 + xi_207 + xi_208 + xi_209;
        const real_t xi_219 = xi_171 * xi_186;
        const real_t xi_221 = xi_219 + xi_220;
        const real_t xi_222 = -xi_212 + xi_214 - xi_216 + xi_218 + xi_221;
        const real_t xi_227 = xi_189 - xi_223 + xi_224 - xi_225 + xi_226;
        const real_t xi_228 = xi_202 + xi_223 - xi_224 + xi_225 - xi_226;
        const real_t xi_229 =
            xi_149 - xi_204 + xi_205 - xi_206 + xi_207 + xi_208 + xi_209;
        const real_t xi_231 = -xi_183;
        const real_t xi_234 =
            xi_181 + xi_188 + xi_221 + xi_230 + xi_231 + xi_232 + xi_233;
        const real_t xi_236 =
            xi_196 *
            (u_2 * xi_153 + xi_10 + xi_156 +
             _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]);
        const real_t xi_237 = -xi_235 - xi_236;
        const real_t xi_242 = xi_193 - xi_238 + xi_239 - xi_240 + xi_241;
        const real_t xi_243 = xi_235 + xi_236;
        const real_t xi_244 = xi_200 + xi_238 - xi_239 + xi_240 - xi_241;
        const real_t xi_245 = -xi_219 + xi_220;
        const real_t xi_246 = xi_212 - xi_214 + xi_216 - xi_218 + xi_245;
        const real_t xi_247 =
            xi_179 + xi_188 + xi_230 + xi_231 + xi_232 + xi_233 + xi_245;
        const real_t forceTerm_0 =
            xi_32 * -real_t(1.5) - xi_35 * xi_38 - xi_38 * xi_40 - xi_38 * xi_41;
        const real_t forceTerm_1 = xi_42 + xi_49;
        const real_t forceTerm_2 = -xi_42 + xi_49;
        const real_t forceTerm_3 = -xi_50 + xi_53;
        const real_t forceTerm_4 = xi_50 + xi_53;
        const real_t forceTerm_5 = xi_54 + xi_55;
        const real_t forceTerm_6 = -xi_54 + xi_55;
        const real_t forceTerm_7 = xi_60 + xi_62 + xi_69;
        const real_t forceTerm_8 = xi_59 + xi_69 + xi_70;
        const real_t forceTerm_9 = xi_59 + xi_62 + xi_72;
        const real_t forceTerm_10 = xi_60 + xi_70 + xi_72;
        const real_t forceTerm_11 = xi_68 + xi_73 + xi_77;
        const real_t forceTerm_12 = xi_71 + xi_77 + xi_78;
        const real_t forceTerm_13 = xi_62 + xi_80 + xi_82;
        const real_t forceTerm_14 = xi_70 + xi_79 + xi_82;
        const real_t forceTerm_15 = xi_68 + xi_78 + xi_84;
        const real_t forceTerm_16 = xi_71 + xi_73 + xi_84;
        const real_t forceTerm_17 = xi_62 + xi_79 + xi_85;
        const real_t forceTerm_18 = xi_70 + xi_80 + xi_85;
        _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_0 + xi_108 * -real_t(0.5) + xi_112 * real_t(0.0238095238095238) +
            xi_88 * real_t(0.142857142857143) + xi_89 * real_t(0.2) - xi_91 +
            xi_92 * real_t(0.0857142857142857) + xi_95 * real_t(0.1) +
            xi_98 * 0.0428571428571429) +
            _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_1 - xi_116 + xi_121 - xi_126 + xi_137 + xi_151) +
            _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_2 + xi_116 - xi_121 + xi_126 + xi_151 + xi_152) +
            _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_3 - xi_155 + xi_158 + xi_160 + xi_165 + xi_166) +
            _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_4 + xi_155 - xi_158 - xi_160 + xi_166 + xi_167) +
            _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_5 - xi_170 + xi_172 - xi_174 + xi_179 + xi_180) +
            _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_6 + xi_170 - xi_172 + xi_174 + xi_180 + xi_181) +
            _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_7 + xi_190 + xi_194 + xi_198) +
            _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_8 + xi_190 + xi_199 + xi_201) +
            _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_9 + xi_194 + xi_199 + xi_203) +
            _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_10 + xi_198 + xi_201 + xi_203) +
            _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_11 + xi_210 + xi_222 + xi_227) +
            _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_12 + xi_222 + xi_228 + xi_229) +
            _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_13 + xi_234 + xi_237 + xi_242) +
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_14 + xi_234 + xi_243 + xi_244) +
            _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_15 + xi_227 + xi_229 + xi_246) +
            _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_16 + xi_210 + xi_228 + xi_246) +
            _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_17 + xi_242 + xi_243 + xi_247) +
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_18 + xi_237 + xi_244 + xi_247) +
            _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_kernel_streamCollide
namespace internal_kernel_collide {
static FUNC_PREFIX void
kernel_collide(real_t *RESTRICT const _data_force, real_t *RESTRICT _data_pdfs,
               int64_t const _size_force_0, int64_t const _size_force_1,
               int64_t const _size_force_2, int64_t const _stride_force_0,
               int64_t const _stride_force_1, int64_t const _stride_force_2,
               int64_t const _stride_force_3, int64_t const _stride_pdfs_0,
               int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
               int64_t const _stride_pdfs_3, uint32_t block_offset_0,
               uint32_t block_offset_1, uint32_t block_offset_2,
               real_t omega_bulk, real_t omega_even, real_t omega_odd,
               real_t omega_shear, uint32_t seed, real_t temperature,
               uint32_t time_step) {
  const real_t xi_25 = -omega_bulk;
  const real_t xi_36 = -omega_shear;
  const real_t xi_37 = xi_36 + real_t(2.0);
  const real_t xi_38 = xi_37 * real_t(0.5);
  const real_t xi_43 = xi_37 * real_t(0.0833333333333333);
  const real_t xi_48 = xi_37 * real_t(0.166666666666667);
  const real_t xi_58 = xi_37 * real_t(0.25);
  const real_t xi_63 = xi_37 * real_t(0.0416666666666667);
  const real_t xi_90 = real_t(2.4494897427831779);
  const real_t xi_115 = omega_odd * real_t(0.25);
  const real_t xi_131 = omega_odd * real_t(0.0833333333333333);
  const real_t xi_196 = omega_shear * real_t(0.25);
  const real_t xi_211 = omega_odd * real_t(0.0416666666666667);
  const real_t xi_213 = omega_odd * real_t(0.125);
  const int64_t rr_0 = 0;
  const real_t xi_120 = real_c(rr_0) * real_t(0.166666666666667);
  const real_t xi_186 = real_c(rr_0) * real_t(0.0833333333333333);
  for (int ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    real_t *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    real_t *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    real_t *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    real_t *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    real_t *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    real_t *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    real_t *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    real_t *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    for (int ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      real_t *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      real_t *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      real_t *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      real_t *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      real_t *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      real_t *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      real_t *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      real_t *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      real_t *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      real_t *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      real_t *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      real_t *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      real_t *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      real_t *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      real_t *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      real_t *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      real_t *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      real_t *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      real_t *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      real_t *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      real_t *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      real_t *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      for (int ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const real_t xi_248 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_249 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_250 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_251 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_252 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_253 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_254 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_255 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const real_t xi_256 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_257 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_258 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_259 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_260 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_261 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_262 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_263 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_264 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_265 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_266 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_267 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_268 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const real_t xi_269 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];

        float Dummy_299;
        float Dummy_300;
        float Dummy_301;
        float Dummy_302;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 3, seed, Dummy_299, Dummy_300,
                      Dummy_301, Dummy_302);

        float Dummy_295;
        float Dummy_296;
        float Dummy_297;
        float Dummy_298;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 2, seed, Dummy_295, Dummy_296,
                      Dummy_297, Dummy_298);

        float Dummy_291;
        float Dummy_292;
        float Dummy_293;
        float Dummy_294;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 1, seed, Dummy_291, Dummy_292,
                      Dummy_293, Dummy_294);

        float Dummy_287;
        float Dummy_288;
        float Dummy_289;
        float Dummy_290;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 0, seed, Dummy_287, Dummy_288,
                      Dummy_289, Dummy_290);

        const real_t xi_0 = xi_264 + xi_265;
        const real_t xi_1 = xi_0 + xi_254;
        const real_t xi_2 = xi_263 + xi_267 + xi_269;
        const real_t xi_3 = xi_250 + xi_262;
        const real_t xi_4 = xi_252 + xi_260;
        const real_t xi_5 = xi_248 + xi_266;
        const real_t xi_6 = xi_257 + xi_261;
        const real_t xi_8 = -xi_260;
        const real_t xi_9 = -xi_251 + xi_8;
        const real_t xi_10 = -xi_257;
        const real_t xi_11 = -xi_259;
        const real_t xi_12 = -xi_252;
        const real_t xi_13 = xi_10 + xi_11 + xi_12;
        const real_t xi_14 = -xi_266;
        const real_t xi_15 = -xi_249;
        const real_t xi_16 = xi_14 + xi_15;
        const real_t xi_17 = -xi_248;
        const real_t xi_18 = -xi_262;
        const real_t xi_19 = xi_17 + xi_18;
        const real_t xi_20 = -xi_265;
        const real_t xi_21 = xi_10 + xi_20;
        const real_t xi_22 = -xi_269;
        const real_t xi_23 = -xi_261;
        const real_t xi_24 = xi_17 + xi_22 + xi_23 + xi_267;
        const real_t xi_42 = xi_253 * real_t(0.166666666666667);
        const real_t xi_50 = xi_255 * real_t(0.166666666666667);
        const real_t xi_54 = xi_268 * real_t(0.166666666666667);
        const real_t xi_57 = xi_253 * real_t(0.5);
        const real_t xi_61 = xi_255 * real_t(0.0833333333333333);
        const real_t xi_65 = xi_253 * real_t(0.0833333333333333);
        const real_t xi_75 = xi_268 * real_t(0.0833333333333333);
        const real_t xi_93 = -xi_256;
        const real_t xi_94 = xi_250 * real_t(3.0) + xi_261 * real_t(3.0) + xi_93;
        const real_t xi_95 =
            omega_even * (xi_248 * -real_t(3.0) + xi_262 * -real_t(3.0) + xi_263 * real_t(3.0) +
                          xi_266 * real_t(3.0) + xi_267 * -real_t(3.0) + xi_269 * -real_t(3.0) + xi_94);
        const real_t xi_96 =
            xi_248 * real_t(2.0) + xi_262 * real_t(2.0) + xi_267 * real_t(2.0) + xi_269 * real_t(2.0);
        const real_t xi_97 = xi_252 * real_t(5.0) + xi_254 * real_t(5.0) + xi_96;
        const real_t xi_98 =
            omega_even *
            (xi_257 * -real_t(5.0) + xi_259 * -real_t(5.0) + xi_263 * -real_t(2.0) + xi_264 * -real_t(5.0) +
             xi_265 * -real_t(5.0) + xi_266 * -real_t(2.0) + xi_94 + xi_97);
        const real_t xi_101 = -xi_267;
        const real_t xi_102 = xi_101 + xi_18;
        const real_t xi_103 = -xi_258;
        const real_t xi_106 = -xi_264;
        const real_t xi_107 = xi_106 + xi_11 + xi_15 + xi_21;
        const real_t xi_109 = xi_259 * real_t(2.0);
        const real_t xi_110 = xi_264 * real_t(2.0);
        const real_t xi_111 = xi_257 * real_t(2.0) + xi_265 * real_t(2.0);
        const real_t xi_112 =
            omega_even *
            (xi_109 + xi_110 + xi_111 + xi_249 * -real_t(7.0) + xi_250 * -real_t(4.0) +
             xi_251 * -real_t(7.0) + xi_258 * -real_t(7.0) + xi_260 * -real_t(7.0) + xi_261 * -real_t(4.0) +
             xi_263 * real_t(5.0) + xi_266 * real_t(5.0) + xi_93 + xi_97);
        const real_t xi_113 = xi_101 + xi_262;
        const real_t xi_114 = xi_113 + xi_14 + xi_22 + xi_248 + xi_263;
        const real_t xi_116 = xi_114 * xi_115;
        const real_t xi_118 = xi_103 + xi_249;
        const real_t xi_122 = Dummy_298 - real_t(0.5);
        const real_t xi_127 = xi_251 * real_t(2.0);
        const real_t xi_128 = xi_249 * real_t(2.0);
        const real_t xi_129 = xi_258 * -real_t(2.0) + xi_260 * real_t(2.0);
        const real_t xi_130 = -xi_127 + xi_128 + xi_129 + xi_14 + xi_19 + xi_2;
        const real_t xi_132 = xi_130 * xi_131;
        const real_t xi_133 = Dummy_293 - real_t(0.5);
        const real_t xi_138 = Dummy_288 - real_t(0.5);
        const real_t xi_142 = xi_257 + xi_259;
        const real_t xi_156 = xi_106 + xi_259;
        const real_t xi_157 = xi_12 + xi_156 + xi_20 + xi_254 + xi_257;
        const real_t xi_158 = xi_115 * xi_157;
        const real_t xi_159 = Dummy_296 - real_t(0.5);
        const real_t xi_161 = xi_1 + xi_127 - xi_128 + xi_129 + xi_13;
        const real_t xi_162 = xi_131 * xi_161;
        const real_t xi_163 = Dummy_295 - real_t(0.5);
        const real_t xi_168 = xi_248 + xi_269;
        const real_t xi_169 = xi_102 + xi_168 + xi_23 + xi_250;
        const real_t xi_170 = xi_115 * xi_169;
        const real_t xi_173 = Dummy_297 - real_t(0.5);
        const real_t xi_175 = -xi_109 - xi_110 + xi_111 + xi_24 + xi_3;
        const real_t xi_176 = xi_131 * xi_175;
        const real_t xi_177 = Dummy_294 - real_t(0.5);
        const real_t xi_184 = xi_112 * real_t(0.0138888888888889);
        const real_t xi_205 = xi_98 * -real_t(0.00714285714285714);
        const real_t xi_207 = xi_95 * real_t(0.025);
        const real_t xi_212 = xi_175 * xi_211;
        const real_t xi_214 = xi_169 * xi_213;
        const real_t xi_223 = xi_130 * xi_211;
        const real_t xi_224 = xi_114 * xi_213;
        const real_t xi_232 = xi_98 * real_t(0.0178571428571429);
        const real_t xi_238 = xi_157 * xi_213;
        const real_t xi_239 = xi_161 * xi_211;
        const real_t vel0Term = xi_1 + xi_249 + xi_258;
        const real_t vel1Term = xi_2 + xi_251;
        const real_t vel2Term = xi_259 + xi_3;
        const real_t rho =
            vel0Term + vel1Term + vel2Term + xi_256 + xi_4 + xi_5 + xi_6;
        const real_t xi_7 = 1 / (rho);
        const real_t xi_86 = rho * temperature;
        const real_t xi_87 =
            sqrt(xi_86 * (-((-omega_even + real_t(1.0)) * (-omega_even +  real_t(1.0))) +  real_t(1.0)));
        const real_t xi_88 = xi_87 * (Dummy_299 - real_t(0.5)) * real_t(3.7416573867739413);
        const real_t xi_89 = xi_87 * (Dummy_301 - real_t(0.5)) * real_t(5.4772255750516612);
        const real_t xi_91 =
            xi_90 * sqrt(xi_86 * (-((xi_25 +  real_t(1.0)) * (xi_25 +  real_t(1.0))) +  real_t(1.0))) *
            (Dummy_292 - real_t(0.5));
        const real_t xi_92 = xi_87 * (Dummy_300 - real_t(0.5)) * real_t(8.3666002653407556);
        const real_t xi_123 =
            sqrt(xi_86 * (-((-omega_odd +  real_t(1.0)) * (-omega_odd +  real_t(1.0))) +  real_t(1.0)));
        const real_t xi_124 = xi_123 * real_t(1.4142135623730951);
        const real_t xi_125 = xi_124 * real_t(0.5);
        const real_t xi_126 = xi_122 * xi_125;
        const real_t xi_134 = xi_123 * xi_90;
        const real_t xi_135 = xi_134 * real_t(0.166666666666667);
        const real_t xi_136 = xi_133 * xi_135;
        const real_t xi_137 = -xi_132 - xi_136;
        const real_t xi_139 =
            sqrt(xi_86 * (-((xi_36 +  real_t(1.0)) * (xi_36 +  real_t(1.0))) +  real_t(1.0)));
        const real_t xi_140 = xi_139 * real_t(0.5);
        const real_t xi_141 = xi_138 * xi_140;
        const real_t xi_146 =
            xi_112 * -real_t(0.0198412698412698) + xi_88 * -real_t(0.119047619047619);
        const real_t xi_148 = xi_139 * (Dummy_287 - real_t(0.5)) * real_t(1.7320508075688772);
        const real_t xi_152 = xi_132 + xi_136;
        const real_t xi_160 = xi_125 * xi_159;
        const real_t xi_164 = xi_135 * xi_163;
        const real_t xi_165 = xi_162 + xi_164;
        const real_t xi_167 = -xi_162 - xi_164;
        const real_t xi_174 = xi_125 * xi_173;
        const real_t xi_178 = xi_135 * xi_177;
        const real_t xi_179 = -xi_176 - xi_178;
        const real_t xi_181 = xi_176 + xi_178;
        const real_t xi_182 = xi_138 * xi_139 * real_t(0.25);
        const real_t xi_185 = xi_88 * real_t(0.0833333333333333);
        const real_t xi_195 = xi_140 * (Dummy_289 - real_t(0.5));
        const real_t xi_204 = xi_140 * (Dummy_291 - real_t(0.5));
        const real_t xi_208 = xi_92 * -real_t(0.0142857142857143);
        const real_t xi_209 = xi_89 * real_t(0.05);
        const real_t xi_215 = xi_134 * real_t(0.0833333333333333);
        const real_t xi_216 = xi_177 * xi_215;
        const real_t xi_217 = xi_124 * real_t(0.25);
        const real_t xi_218 = xi_173 * xi_217;
        const real_t xi_220 =
            xi_112 * -real_t(0.00396825396825397) + xi_88 * -real_t(0.0238095238095238);
        const real_t xi_225 = xi_133 * xi_215;
        const real_t xi_226 = xi_122 * xi_217;
        const real_t xi_230 = -xi_182;
        const real_t xi_233 = xi_92 * real_t(0.0357142857142857);
        const real_t xi_235 = xi_140 * (Dummy_290 - real_t(0.5));
        const real_t xi_240 = xi_159 * xi_217;
        const real_t xi_241 = xi_163 * xi_215;
        const real_t u_0 = xi_7 * (vel0Term + xi_13 + xi_9);
        const real_t xi_26 = u_0 * xi_255;
        const real_t xi_27 = xi_26 * real_t(0.333333333333333);
        const real_t xi_33 = -xi_27;
        const real_t xi_99 = rho * (u_0 * u_0);
        const real_t xi_153 = rho * u_0;
        const real_t xi_154 = -vel0Term + xi_142 + xi_153 + xi_251 + xi_4;
        const real_t xi_155 = xi_120 * xi_154;
        const real_t xi_191 = xi_154 * xi_186;
        const real_t u_1 = xi_7 * (vel1Term + xi_16 + xi_19 + xi_258 + xi_8);
        const real_t xi_28 = u_1 * xi_253;
        const real_t xi_29 = xi_28 * real_t(0.333333333333333);
        const real_t xi_34 = -xi_29;
        const real_t xi_56 = u_1 * real_t(0.5);
        const real_t xi_59 = xi_58 * (u_0 * xi_57 + xi_255 * xi_56);
        const real_t xi_60 = -xi_59;
        const real_t xi_104 = rho * (u_1 * u_1);
        const real_t xi_105 = xi_103 + xi_104 + xi_9;
        const real_t xi_117 = rho * u_1;
        const real_t xi_119 =
            -vel1Term + xi_117 + xi_118 + xi_260 + xi_262 + xi_5;
        const real_t xi_121 = xi_119 * xi_120;
        const real_t xi_187 = xi_119 * xi_186;
        const real_t xi_197 = xi_196 * (u_0 * xi_117 + xi_118 + xi_251 + xi_8);
        const real_t xi_198 = -xi_195 - xi_197;
        const real_t xi_199 = xi_195 + xi_197;
        const real_t u_2 = xi_7 * (vel2Term + xi_21 + xi_24 + xi_264);
        const real_t xi_30 = u_2 * xi_268;
        const real_t xi_31 = xi_30 * real_t(0.333333333333333);
        const real_t xi_32 = (xi_25 +  real_t(2.0)) * (xi_27 + xi_29 + xi_31);
        const real_t xi_35 = xi_30 * real_t(0.666666666666667) + xi_33 + xi_34;
        const real_t xi_39 = -xi_31;
        const real_t xi_40 = xi_28 * real_t(0.666666666666667) + xi_33 + xi_39;
        const real_t xi_41 = xi_26 * real_t(0.666666666666667) + xi_34 + xi_39;
        const real_t xi_44 = xi_35 * xi_43;
        const real_t xi_45 = -xi_44;
        const real_t xi_46 = xi_41 * xi_43;
        const real_t xi_47 = -xi_46;
        const real_t xi_49 = xi_40 * xi_48 + xi_45 + xi_47;
        const real_t xi_51 = xi_40 * xi_43;
        const real_t xi_52 = -xi_51;
        const real_t xi_53 = xi_41 * xi_48 + xi_45 + xi_52;
        const real_t xi_55 = xi_35 * xi_48 + xi_47 + xi_52;
        const real_t xi_62 = xi_46 - xi_61;
        const real_t xi_64 = -xi_35 * xi_63;
        const real_t xi_66 = xi_32 * real_t(0.125);
        const real_t xi_67 = xi_51 + xi_66;
        const real_t xi_68 = xi_65 + xi_67;
        const real_t xi_69 = xi_64 + xi_68;
        const real_t xi_70 = xi_46 + xi_61;
        const real_t xi_71 = -xi_65 + xi_67;
        const real_t xi_72 = xi_64 + xi_71;
        const real_t xi_73 = xi_58 * (u_2 * xi_57 + xi_268 * xi_56);
        const real_t xi_74 = -xi_41 * xi_63;
        const real_t xi_76 = xi_44 + xi_75;
        const real_t xi_77 = xi_74 + xi_76;
        const real_t xi_78 = -xi_73;
        const real_t xi_79 = xi_58 * (u_0 * xi_268 * real_t(0.5) + u_2 * xi_255 * real_t(0.5));
        const real_t xi_80 = -xi_79;
        const real_t xi_81 = -xi_40 * xi_63;
        const real_t xi_82 = xi_66 + xi_76 + xi_81;
        const real_t xi_83 = xi_44 - xi_75;
        const real_t xi_84 = xi_74 + xi_83;
        const real_t xi_85 = xi_66 + xi_81 + xi_83;
        const real_t xi_100 = rho * (u_2 * u_2);
        const real_t xi_108 = omega_bulk * (xi_100 + xi_102 + xi_105 + xi_107 +
                                            xi_17 + xi_22 + xi_256 + xi_99);
        const real_t xi_143 = -xi_100 + xi_250 + xi_261;
        const real_t xi_144 =
            omega_shear * (xi_0 + xi_105 + xi_142 + xi_143 + xi_16 - xi_263);
        const real_t xi_145 = xi_144 * real_t(0.125);
        const real_t xi_147 =
            omega_shear *
            (xi_103 - xi_104 + xi_107 + xi_143 + xi_252 * -real_t(2.0) + xi_254 * -real_t(2.0) +
             xi_263 + xi_266 + xi_9 + xi_96 + xi_99 *  real_t(2.0));
        const real_t xi_149 =
            xi_147 * -real_t(0.0416666666666667) + xi_148 * -real_t(0.166666666666667);
        const real_t xi_150 = xi_149 + xi_89 * -real_t(0.1) + xi_95 * -real_t(0.05);
        const real_t xi_151 = xi_141 + xi_145 + xi_146 + xi_150 +
                              xi_92 * real_t(0.0285714285714286) +
                              xi_98 * real_t(0.0142857142857143);
        const real_t xi_166 =
            xi_146 + xi_147 * real_t(0.0833333333333333) + xi_148 * real_t(0.333333333333333) +
            xi_92 * -real_t(0.0714285714285714) + xi_98 * -real_t(0.0357142857142857);
        const real_t xi_171 =
            rho * u_2 - vel2Term + xi_101 + xi_106 + xi_168 + xi_265 + xi_6;
        const real_t xi_172 = xi_120 * xi_171;
        const real_t xi_180 = xi_112 * real_t(0.0158730158730159) - xi_141 - xi_145 +
                              xi_150 + xi_88 * real_t(0.0952380952380952) +
                              xi_92 * -real_t(0.0428571428571429) +
                              xi_98 * -real_t(0.0214285714285714);
        const real_t xi_183 = xi_144 * real_t(0.0625);
        const real_t xi_188 =
            xi_108 * real_t(0.0416666666666667) + xi_91 * real_t(0.0833333333333333);
        const real_t xi_189 = xi_187 + xi_188;
        const real_t xi_190 =
            xi_152 + xi_182 + xi_183 + xi_184 + xi_185 + xi_189;
        const real_t xi_192 =
            xi_147 * real_t(0.0208333333333333) + xi_148 * real_t(0.0833333333333333);
        const real_t xi_193 = -xi_191 + xi_192;
        const real_t xi_194 = xi_167 + xi_193;
        const real_t xi_200 = xi_191 + xi_192;
        const real_t xi_201 = xi_165 + xi_200;
        const real_t xi_202 = -xi_187 + xi_188;
        const real_t xi_203 =
            xi_137 + xi_182 + xi_183 + xi_184 + xi_185 + xi_202;
        const real_t xi_206 = xi_196 * (u_2 * xi_117 + xi_113 + xi_17 + xi_269);
        const real_t xi_210 =
            xi_149 + xi_204 + xi_205 + xi_206 + xi_207 + xi_208 + xi_209;
        const real_t xi_219 = xi_171 * xi_186;
        const real_t xi_221 = xi_219 + xi_220;
        const real_t xi_222 = -xi_212 + xi_214 - xi_216 + xi_218 + xi_221;
        const real_t xi_227 = xi_189 - xi_223 + xi_224 - xi_225 + xi_226;
        const real_t xi_228 = xi_202 + xi_223 - xi_224 + xi_225 - xi_226;
        const real_t xi_229 =
            xi_149 - xi_204 + xi_205 - xi_206 + xi_207 + xi_208 + xi_209;
        const real_t xi_231 = -xi_183;
        const real_t xi_234 =
            xi_181 + xi_188 + xi_221 + xi_230 + xi_231 + xi_232 + xi_233;
        const real_t xi_236 = xi_196 * (u_2 * xi_153 + xi_10 + xi_156 + xi_265);
        const real_t xi_237 = -xi_235 - xi_236;
        const real_t xi_242 = xi_193 - xi_238 + xi_239 - xi_240 + xi_241;
        const real_t xi_243 = xi_235 + xi_236;
        const real_t xi_244 = xi_200 + xi_238 - xi_239 + xi_240 - xi_241;
        const real_t xi_245 = -xi_219 + xi_220;
        const real_t xi_246 = xi_212 - xi_214 + xi_216 - xi_218 + xi_245;
        const real_t xi_247 =
            xi_179 + xi_188 + xi_230 + xi_231 + xi_232 + xi_233 + xi_245;
        const real_t forceTerm_0 =
            xi_32 * -real_t(1.5) - xi_35 * xi_38 - xi_38 * xi_40 - xi_38 * xi_41;
        const real_t forceTerm_1 = xi_42 + xi_49;
        const real_t forceTerm_2 = -xi_42 + xi_49;
        const real_t forceTerm_3 = -xi_50 + xi_53;
        const real_t forceTerm_4 = xi_50 + xi_53;
        const real_t forceTerm_5 = xi_54 + xi_55;
        const real_t forceTerm_6 = -xi_54 + xi_55;
        const real_t forceTerm_7 = xi_60 + xi_62 + xi_69;
        const real_t forceTerm_8 = xi_59 + xi_69 + xi_70;
        const real_t forceTerm_9 = xi_59 + xi_62 + xi_72;
        const real_t forceTerm_10 = xi_60 + xi_70 + xi_72;
        const real_t forceTerm_11 = xi_68 + xi_73 + xi_77;
        const real_t forceTerm_12 = xi_71 + xi_77 + xi_78;
        const real_t forceTerm_13 = xi_62 + xi_80 + xi_82;
        const real_t forceTerm_14 = xi_70 + xi_79 + xi_82;
        const real_t forceTerm_15 = xi_68 + xi_78 + xi_84;
        const real_t forceTerm_16 = xi_71 + xi_73 + xi_84;
        const real_t forceTerm_17 = xi_62 + xi_79 + xi_85;
        const real_t forceTerm_18 = xi_70 + xi_80 + xi_85;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_0 + xi_108 * -real_t(0.5) + xi_112 * real_t(0.0238095238095238) + xi_256 +
            xi_88 * real_t(0.142857142857143) + xi_89 * real_t(0.2) - xi_91 +
            xi_92 * real_t(0.0857142857142857) + xi_95 * real_t(0.1) +
            xi_98 * 0.0428571428571429);
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_1 - xi_116 + xi_121 - xi_126 + xi_137 + xi_151 + xi_263);
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_2 + xi_116 - xi_121 + xi_126 + xi_151 + xi_152 + xi_266);
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_3 - xi_155 + xi_158 + xi_160 + xi_165 + xi_166 + xi_252);
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_4 + xi_155 - xi_158 - xi_160 + xi_166 + xi_167 + xi_254);
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_5 - xi_170 + xi_172 - xi_174 + xi_179 + xi_180 + xi_250);
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_6 + xi_170 - xi_172 + xi_174 + xi_180 + xi_181 + xi_261);
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_7 + xi_190 + xi_194 + xi_198 + xi_251);
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_8 + xi_190 + xi_199 + xi_201 + xi_258);
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_9 + xi_194 + xi_199 + xi_203 + xi_260);
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_10 + xi_198 + xi_201 + xi_203 + xi_249);
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_11 + xi_210 + xi_222 + xi_227 + xi_267);
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_12 + xi_222 + xi_228 + xi_229 + xi_262);
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_13 + xi_234 + xi_237 + xi_242 + xi_259);
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_14 + xi_234 + xi_243 + xi_244 + xi_264);
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_15 + xi_227 + xi_229 + xi_246 + xi_269);
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_16 + xi_210 + xi_228 + xi_246 + xi_248);
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_17 + xi_242 + xi_243 + xi_247 + xi_257);
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_18 + xi_237 + xi_244 + xi_247 + xi_265);
      }
    }
  }
}
} // namespace internal_kernel_collide
namespace internal_kernel_stream {
static FUNC_PREFIX void kernel_stream(
    real_t *RESTRICT const _data_pdfs, real_t *RESTRICT _data_pdfs_tmp,
    int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
    int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_0,
    int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2,
    int64_t const _stride_pdfs_tmp_3) {
  for (int ctr_2 = 1; ctr_2 < _size_pdfs_2 - 1; ctr_2 += 1) {
    real_t *RESTRICT _data_pdfs_tmp_20_30 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2;
    real_t *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    real_t *RESTRICT _data_pdfs_tmp_20_31 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_32 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_33 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_34 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_35 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_2m1_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                         _stride_pdfs_2 + 5 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_36 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_21_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                        _stride_pdfs_2 + 6 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_37 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_38 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_39 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_310 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_311 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_2m1_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 11 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_312 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_2m1_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 12 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_313 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_2m1_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 13 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_314 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_2m1_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 14 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_315 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_21_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 15 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_316 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_21_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 16 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_317 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_21_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 17 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_tmp_20_318 =
        _data_pdfs_tmp + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3;
    real_t *RESTRICT _data_pdfs_21_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 18 * _stride_pdfs_3;
    for (int ctr_1 = 1; ctr_1 < _size_pdfs_1 - 1; ctr_1 += 1) {
      real_t *RESTRICT _data_pdfs_tmp_20_30_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_30;
      real_t *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      real_t *RESTRICT _data_pdfs_tmp_20_31_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_31;
      real_t *RESTRICT _data_pdfs_20_31_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_31;
      real_t *RESTRICT _data_pdfs_tmp_20_32_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_32;
      real_t *RESTRICT _data_pdfs_20_32_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_32;
      real_t *RESTRICT _data_pdfs_tmp_20_33_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_33;
      real_t *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      real_t *RESTRICT _data_pdfs_tmp_20_34_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_34;
      real_t *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      real_t *RESTRICT _data_pdfs_tmp_20_35_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_35;
      real_t *RESTRICT _data_pdfs_2m1_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_35;
      real_t *RESTRICT _data_pdfs_tmp_20_36_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_36;
      real_t *RESTRICT _data_pdfs_21_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_36;
      real_t *RESTRICT _data_pdfs_tmp_20_37_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_37;
      real_t *RESTRICT _data_pdfs_20_37_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_37;
      real_t *RESTRICT _data_pdfs_tmp_20_38_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_38;
      real_t *RESTRICT _data_pdfs_20_38_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_38;
      real_t *RESTRICT _data_pdfs_tmp_20_39_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_39;
      real_t *RESTRICT _data_pdfs_20_39_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_39;
      real_t *RESTRICT _data_pdfs_tmp_20_310_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_310;
      real_t *RESTRICT _data_pdfs_20_310_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_310;
      real_t *RESTRICT _data_pdfs_tmp_20_311_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_311;
      real_t *RESTRICT _data_pdfs_2m1_311_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_311;
      real_t *RESTRICT _data_pdfs_tmp_20_312_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_312;
      real_t *RESTRICT _data_pdfs_2m1_312_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_312;
      real_t *RESTRICT _data_pdfs_tmp_20_313_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_313;
      real_t *RESTRICT _data_pdfs_2m1_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_313;
      real_t *RESTRICT _data_pdfs_tmp_20_314_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_314;
      real_t *RESTRICT _data_pdfs_2m1_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_314;
      real_t *RESTRICT _data_pdfs_tmp_20_315_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_315;
      real_t *RESTRICT _data_pdfs_21_315_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_315;
      real_t *RESTRICT _data_pdfs_tmp_20_316_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_316;
      real_t *RESTRICT _data_pdfs_21_316_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_316;
      real_t *RESTRICT _data_pdfs_tmp_20_317_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_317;
      real_t *RESTRICT _data_pdfs_21_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_317;
      real_t *RESTRICT _data_pdfs_tmp_20_318_10 =
          _stride_pdfs_tmp_1 * ctr_1 + _data_pdfs_tmp_20_318;
      real_t *RESTRICT _data_pdfs_21_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_318;
      for (int ctr_0 = 1; ctr_0 < _size_pdfs_0 - 1; ctr_0 += 1) {
        _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0 * ctr_0] =
            _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_kernel_stream

const real_t FluctuatingMRT_LatticeModel::w[19] = {
    real_t(0.333333333333333),  real_t(0.0555555555555556), real_t(0.0555555555555556),
    real_t(0.0555555555555556), real_t(0.0555555555555556), real_t(0.0555555555555556),
    real_t(0.0555555555555556), real_t(0.0277777777777778), real_t(0.0277777777777778),
    real_t(0.0277777777777778), real_t(0.0277777777777778), real_t(0.0277777777777778),
    real_t(0.0277777777777778), real_t(0.0277777777777778), real_t(0.0277777777777778),
    real_t(0.0277777777777778), real_t(0.0277777777777778), real_t(0.0277777777777778),
    real_t(0.0277777777777778)};
const real_t FluctuatingMRT_LatticeModel::wInv[19] = {
    real_t(3.00000000000000), real_t(18.0000000000000), real_t(18.0000000000000), real_t(18.0000000000000),
    real_t(18.0000000000000), real_t(18.0000000000000), real_t(18.0000000000000), real_t(36.0000000000000),
    real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000),
    real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000),
    real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000)};

void FluctuatingMRT_LatticeModel::Sweep::streamCollide(
    IBlock *block, const uint_t numberOfGhostLayersToInclude) {
  auto pdfs = block->getData<field::GhostLayerField<real_t, 19>>(pdfsID);
  field::GhostLayerField<real_t, 19> *pdfs_tmp;
  {
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find(pdfs);
    if (it != cache_pdfs_.end()) {
      pdfs_tmp = *it;
    } else {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_.insert(pdfs_tmp);
    }
  }

  auto &lm = dynamic_cast<lbm::PdfField<FluctuatingMRT_LatticeModel> *>(pdfs)
                 ->latticeModel();
  WALBERLA_ASSERT_EQUAL(*(lm.blockId_), block->getId());

  auto &omega_even = lm.omega_even_;
  auto &temperature = lm.temperature_;
  auto &omega_odd = lm.omega_odd_;
  auto &seed = lm.seed_;
  auto &block_offset_0 = lm.block_offset_0_;
  auto &block_offset_2 = lm.block_offset_2_;
  auto &omega_bulk = lm.omega_bulk_;
  auto &time_step = lm.time_step_;
  auto &block_offset_1 = lm.block_offset_1_;
  auto &force = lm.force_;
  auto &omega_shear = lm.omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                                -int_c(force->nrOfGhostLayers()));
  real_t *RESTRICT const _data_force =
      force->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                    -cell_idx_c(numberOfGhostLayersToInclude) - 1,
                    -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                                -int_c(pdfs->nrOfGhostLayers()));
  real_t *RESTRICT const _data_pdfs =
      pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                   -cell_idx_c(numberOfGhostLayersToInclude) - 1,
                   -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                                -int_c(pdfs_tmp->nrOfGhostLayers()));
  real_t *RESTRICT _data_pdfs_tmp =
      pdfs_tmp->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                       -cell_idx_c(numberOfGhostLayersToInclude) - 1,
                       -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(
      force->xSizeWithGhostLayer(),
      int64_t(cell_idx_c(force->xSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2));
  const int64_t _size_force_0 =
      int64_t(cell_idx_c(force->xSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(
      force->ySizeWithGhostLayer(),
      int64_t(cell_idx_c(force->ySize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2));
  const int64_t _size_force_1 =
      int64_t(cell_idx_c(force->ySize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(
      force->zSizeWithGhostLayer(),
      int64_t(cell_idx_c(force->zSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2));
  const int64_t _size_force_2 =
      int64_t(cell_idx_c(force->zSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2);
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  internal_kernel_streamCollide::kernel_streamCollide(
      _data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1,
      _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
      _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
      _stride_pdfs_3, _stride_pdfs_tmp_0, _stride_pdfs_tmp_1,
      _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, block_offset_0, block_offset_1,
      block_offset_2, omega_bulk, omega_even, omega_odd, omega_shear, seed,
      temperature, time_step);
  pdfs->swapDataPointers(pdfs_tmp);
}

void FluctuatingMRT_LatticeModel::Sweep::collide(
    IBlock *block, const uint_t numberOfGhostLayersToInclude) {
  auto pdfs = block->getData<field::GhostLayerField<real_t, 19>>(pdfsID);

  auto &lm = dynamic_cast<lbm::PdfField<FluctuatingMRT_LatticeModel> *>(pdfs)
                 ->latticeModel();
  WALBERLA_ASSERT_EQUAL(*(lm.blockId_), block->getId());

  auto &omega_even = lm.omega_even_;
  auto &temperature = lm.temperature_;
  auto &omega_odd = lm.omega_odd_;
  auto &seed = lm.seed_;
  auto &block_offset_0 = lm.block_offset_0_;
  auto &block_offset_2 = lm.block_offset_2_;
  auto &omega_bulk = lm.omega_bulk_;
  auto &time_step = lm.time_step_;
  auto &block_offset_1 = lm.block_offset_1_;
  auto &force = lm.force_;
  auto &omega_shear = lm.omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude),
                                -int_c(force->nrOfGhostLayers()));
  real_t *RESTRICT const _data_force =
      force->dataAt(-cell_idx_c(numberOfGhostLayersToInclude),
                    -cell_idx_c(numberOfGhostLayersToInclude),
                    -cell_idx_c(numberOfGhostLayersToInclude), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude),
                                -int_c(pdfs->nrOfGhostLayers()));
  real_t *RESTRICT _data_pdfs =
      pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude),
                   -cell_idx_c(numberOfGhostLayersToInclude),
                   -cell_idx_c(numberOfGhostLayersToInclude), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(
      force->xSizeWithGhostLayer(),
      int64_t(cell_idx_c(force->xSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude)));
  const int64_t _size_force_0 =
      int64_t(cell_idx_c(force->xSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude));
  WALBERLA_ASSERT_GREATER_EQUAL(
      force->ySizeWithGhostLayer(),
      int64_t(cell_idx_c(force->ySize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude)));
  const int64_t _size_force_1 =
      int64_t(cell_idx_c(force->ySize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude));
  WALBERLA_ASSERT_GREATER_EQUAL(
      force->zSizeWithGhostLayer(),
      int64_t(cell_idx_c(force->zSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude)));
  const int64_t _size_force_2 =
      int64_t(cell_idx_c(force->zSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude));
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_kernel_collide::kernel_collide(
      _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
      _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
      _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
      block_offset_0, block_offset_1, block_offset_2, omega_bulk, omega_even,
      omega_odd, omega_shear, seed, temperature, time_step);
}

void FluctuatingMRT_LatticeModel::Sweep::stream(
    IBlock *block, const uint_t numberOfGhostLayersToInclude) {
  auto pdfs = block->getData<field::GhostLayerField<real_t, 19>>(pdfsID);
  field::GhostLayerField<real_t, 19> *pdfs_tmp;
  {
    // Getting temporary field pdfs_tmp
    auto it = cache_pdfs_.find(pdfs);
    if (it != cache_pdfs_.end()) {
      pdfs_tmp = *it;
    } else {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_.insert(pdfs_tmp);
    }
  }

  WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                                -int_c(pdfs->nrOfGhostLayers()));
  real_t *RESTRICT const _data_pdfs =
      pdfs->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                   -cell_idx_c(numberOfGhostLayersToInclude) - 1,
                   -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                                -int_c(pdfs_tmp->nrOfGhostLayers()));
  real_t *RESTRICT _data_pdfs_tmp =
      pdfs_tmp->dataAt(-cell_idx_c(numberOfGhostLayersToInclude) - 1,
                       -cell_idx_c(numberOfGhostLayersToInclude) - 1,
                       -cell_idx_c(numberOfGhostLayersToInclude) - 1, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(
      pdfs->xSizeWithGhostLayer(),
      int64_t(cell_idx_c(pdfs->xSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2));
  const int64_t _size_pdfs_0 =
      int64_t(cell_idx_c(pdfs->xSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(
      pdfs->ySizeWithGhostLayer(),
      int64_t(cell_idx_c(pdfs->ySize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2));
  const int64_t _size_pdfs_1 =
      int64_t(cell_idx_c(pdfs->ySize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2);
  WALBERLA_ASSERT_GREATER_EQUAL(
      pdfs->zSizeWithGhostLayer(),
      int64_t(cell_idx_c(pdfs->zSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2));
  const int64_t _size_pdfs_2 =
      int64_t(cell_idx_c(pdfs->zSize()) +
              2 * cell_idx_c(numberOfGhostLayersToInclude) + 2);
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_pdfs_tmp_0 = int64_t(pdfs_tmp->xStride());
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  internal_kernel_stream::kernel_stream(
      _data_pdfs, _data_pdfs_tmp, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
      _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
      _stride_pdfs_tmp_0, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2,
      _stride_pdfs_tmp_3);

  pdfs->swapDataPointers(pdfs_tmp);
}

} // namespace lbm
} // namespace walberla

// Buffer Packing

namespace walberla {
namespace mpi {

mpi::SendBuffer &
operator<<(mpi::SendBuffer &buf,
           const ::walberla::lbm::FluctuatingMRT_LatticeModel &lm) {
  buf << lm.currentLevel;
  return buf;
}

mpi::RecvBuffer &operator>>(mpi::RecvBuffer &buf,
                            ::walberla::lbm::FluctuatingMRT_LatticeModel &lm) {
  buf >> lm.currentLevel;
  return buf;
}

} // namespace mpi
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
