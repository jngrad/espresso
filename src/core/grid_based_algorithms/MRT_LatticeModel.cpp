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

#include "MRT_LatticeModel.h"
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
    real_t omega_bulk, real_t omega_even, real_t omega_odd,
    real_t omega_shear) {
  const real_t xi_35 = -omega_shear + real_t(2.0);
  const real_t xi_36 = xi_35 * real_t(0.5);
  const real_t xi_41 = xi_35 * real_t(0.0833333333333333);
  const real_t xi_46 = xi_35 * real_t(0.166666666666667);
  const real_t xi_56 = xi_35 * real_t(0.25);
  const real_t xi_61 = xi_35 * real_t(0.0416666666666667);
  const real_t xi_106 = omega_odd * real_t(0.25);
  const real_t xi_112 = omega_odd * real_t(0.0833333333333333);
  const real_t xi_148 = omega_shear * real_t(0.25);
  const real_t xi_170 = omega_odd * real_t(0.0416666666666667);
  const real_t xi_172 = omega_odd * real_t(0.125);
  const int64_t rr_0 = 0;
  const real_t xi_118 = real_c(rr_0) * real_t(0.166666666666667);
  const real_t xi_153 = real_c(rr_0) * real_t(0.0833333333333333);
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
        const real_t xi_40 =
            real_t(0.166666666666667) * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_48 =
            real_t(0.166666666666667) * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const real_t xi_52 =
            real_t(0.166666666666667) * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const real_t xi_55 =
            real_t(0.5) * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_59 =
            real_t(0.0833333333333333) * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const real_t xi_63 =
            real_t(0.0833333333333333) * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_73 =
            real_t(0.0833333333333333) * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const real_t xi_84 = -_data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_85 = xi_84 +
                             real_t(3.0) * _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] +
                             real_t(3.0) * _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_86 =
            omega_even *
            (xi_85 - real_t(3.0) * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] -
             real_t(3.0) * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] -
             real_t(3.0) * _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0] -
             real_t(3.0) * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0] +
             real_t(3.0) * _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] +
             real_t(3.0) * _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0]);
        const real_t xi_87 =
            real_t(2.0) * _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] +
            real_t(2.0) * _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0] +
            real_t(2.0) * _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0] +
            real_t(2.0) * _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_88 =
            xi_87 +
            real_t(5.0) * _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            real_t(5.0) * _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_89 =
            omega_even *
            (xi_85 + xi_88 -
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
        const real_t xi_92 = -_data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        const real_t xi_93 = xi_18 + xi_92;
        const real_t xi_94 =
            -_data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_97 =
            -_data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_98 = xi_11 + xi_15 + xi_21 + xi_97;
        const real_t xi_100 =
            real_t(2.0) *
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_101 =
            real_t(2.0) *
            _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_102 =
            real_t(2.0) *
                _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            real_t(2.0) * _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_103 =
            omega_even *
            (xi_100 + xi_101 + xi_102 + xi_84 + xi_88 -
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
        const real_t xi_104 =
            xi_92 + _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_105 = xi_104 + xi_14 + xi_22 +
                              _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] +
                              _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_107 = xi_105 * xi_106;
        const real_t xi_108 =
            real_t(2.0) * _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_109 =
            real_t(2.0) * _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_110 =
            -real_t(2.0) *
                _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] +
            real_t(2.0) * _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_111 = -xi_108 + xi_109 + xi_110 + xi_14 + xi_19 + xi_2;
        const real_t xi_113 = xi_111 * xi_112;
        const real_t xi_114 = -xi_113;
        const real_t xi_116 =
            xi_94 +
            _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_120 = xi_103 * -real_t(0.0198412698412698);
        const real_t xi_121 =
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_132 =
            xi_97 +
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_133 =
            xi_12 + xi_132 + xi_20 +
            _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] +
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_134 = xi_106 * xi_133;
        const real_t xi_135 = xi_1 + xi_108 - xi_109 + xi_110 + xi_13;
        const real_t xi_136 = xi_112 * xi_135;
        const real_t xi_138 = -xi_136;
        const real_t xi_139 = _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0] +
                              _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_140 = xi_139 + xi_23 + xi_93 +
                              _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_141 = xi_106 * xi_140;
        const real_t xi_142 = -xi_100 - xi_101 + xi_102 + xi_24 + xi_3;
        const real_t xi_143 = xi_112 * xi_142;
        const real_t xi_144 = -xi_143;
        const real_t xi_151 = xi_103 * real_t(0.0138888888888889);
        const real_t xi_166 = xi_89 * -real_t(0.00714285714285714);
        const real_t xi_168 = xi_86 * real_t(0.025);
        const real_t xi_171 = xi_142 * xi_170;
        const real_t xi_173 = xi_140 * xi_172;
        const real_t xi_174 = xi_103 * -real_t(0.00396825396825397);
        const real_t xi_178 = xi_111 * xi_170;
        const real_t xi_179 = xi_105 * xi_172;
        const real_t xi_185 = xi_89 * real_t(0.0178571428571429);
        const real_t xi_188 = xi_133 * xi_172;
        const real_t xi_189 = xi_135 * xi_170;
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
        const real_t u_0 = xi_7 * (vel0Term + xi_13 + xi_9);
        const real_t xi_25 =
            u_0 * _data_force_20_30_10[_stride_force_0 * ctr_0];
        const real_t xi_26 = xi_25 * real_t(0.333333333333333);
        const real_t xi_32 = -xi_26;
        const real_t xi_90 = rho * (u_0 * u_0);
        const real_t xi_129 = rho * u_0;
        const real_t xi_130 =
            -vel0Term + xi_121 + xi_129 + xi_4 +
            _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        const real_t xi_131 = xi_118 * xi_130;
        const real_t xi_156 = xi_130 * xi_153;
        const real_t u_1 =
            xi_7 *
            (vel1Term + xi_16 + xi_19 + xi_8 +
             _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]);
        const real_t xi_27 =
            u_1 * _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_28 = xi_27 * real_t(0.333333333333333);
        const real_t xi_33 = -xi_28;
        const real_t xi_54 = u_1 * real_t(0.5);
        const real_t xi_57 =
            xi_56 * (u_0 * xi_55 +
                     xi_54 * _data_force_20_30_10[_stride_force_0 * ctr_0]);
        const real_t xi_58 = -xi_57;
        const real_t xi_95 = rho * (u_1 * u_1);
        const real_t xi_96 = xi_9 + xi_94 + xi_95;
        const real_t xi_115 = rho * u_1;
        const real_t xi_117 =
            -vel1Term + xi_115 + xi_116 + xi_5 +
            _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] +
            _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        const real_t xi_119 = xi_117 * xi_118;
        const real_t xi_149 =
            xi_148 *
            (u_0 * xi_115 + xi_116 + xi_8 +
             _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0]);
        const real_t xi_154 = xi_117 * xi_153;
        const real_t xi_155 = xi_113 + xi_154;
        const real_t xi_164 = -xi_154;
        const real_t xi_165 = xi_114 + xi_164;
        const real_t xi_180 = xi_154 - xi_178 + xi_179;
        const real_t xi_181 = xi_164 + xi_178 - xi_179;
        const real_t u_2 =
            xi_7 *
            (vel2Term + xi_21 + xi_24 +
             _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]);
        const real_t xi_29 =
            u_2 * _data_force_20_32_10[_stride_force_0 * ctr_0];
        const real_t xi_30 = xi_29 * real_t(0.333333333333333);
        const real_t xi_31 = (-omega_bulk + real_t(2.0)) * (xi_26 + xi_28 + xi_30);
        const real_t xi_34 = xi_29 * real_t(0.666666666666667) + xi_32 + xi_33;
        const real_t xi_37 = -xi_30;
        const real_t xi_38 = xi_27 * real_t(0.666666666666667) + xi_32 + xi_37;
        const real_t xi_39 = xi_25 * real_t(0.666666666666667) + xi_33 + xi_37;
        const real_t xi_42 = xi_34 * xi_41;
        const real_t xi_43 = -xi_42;
        const real_t xi_44 = xi_39 * xi_41;
        const real_t xi_45 = -xi_44;
        const real_t xi_47 = xi_38 * xi_46 + xi_43 + xi_45;
        const real_t xi_49 = xi_38 * xi_41;
        const real_t xi_50 = -xi_49;
        const real_t xi_51 = xi_39 * xi_46 + xi_43 + xi_50;
        const real_t xi_53 = xi_34 * xi_46 + xi_45 + xi_50;
        const real_t xi_60 = xi_44 - xi_59;
        const real_t xi_62 = -xi_34 * xi_61;
        const real_t xi_64 = xi_31 * real_t(0.125);
        const real_t xi_65 = xi_49 + xi_64;
        const real_t xi_66 = xi_63 + xi_65;
        const real_t xi_67 = xi_62 + xi_66;
        const real_t xi_68 = xi_44 + xi_59;
        const real_t xi_69 = -xi_63 + xi_65;
        const real_t xi_70 = xi_62 + xi_69;
        const real_t xi_71 =
            xi_56 * (u_2 * xi_55 +
                     xi_54 * _data_force_20_32_10[_stride_force_0 * ctr_0]);
        const real_t xi_72 = -xi_39 * xi_61;
        const real_t xi_74 = xi_42 + xi_73;
        const real_t xi_75 = xi_72 + xi_74;
        const real_t xi_76 = -xi_71;
        const real_t xi_77 =
            xi_56 * (u_0 * real_t(0.5) * _data_force_20_32_10[_stride_force_0 * ctr_0] +
                     u_2 * real_t(0.5) * _data_force_20_30_10[_stride_force_0 * ctr_0]);
        const real_t xi_78 = -xi_77;
        const real_t xi_79 = -xi_38 * xi_61;
        const real_t xi_80 = xi_64 + xi_74 + xi_79;
        const real_t xi_81 = xi_42 - xi_73;
        const real_t xi_82 = xi_72 + xi_81;
        const real_t xi_83 = xi_64 + xi_79 + xi_81;
        const real_t xi_91 = rho * (u_2 * u_2);
        const real_t xi_99 =
            omega_bulk * (xi_17 + xi_22 + xi_90 + xi_91 + xi_93 + xi_96 +
                          xi_98 + _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0]);
        const real_t xi_122 = -xi_91 +
                              _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0] +
                              _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_123 =
            omega_shear * (xi_0 + xi_121 + xi_122 + xi_16 + xi_96 -
                           _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0]);
        const real_t xi_124 = xi_123 * real_t(0.125);
        const real_t xi_125 =
            omega_shear *
            (xi_122 + xi_87 + xi_9 + xi_90 * real_t(2.0) + xi_94 - xi_95 + xi_98 -
             real_t(2.0) *
                 _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] -
             real_t(2.0) *
                 _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] +
             _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0] +
             _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0]);
        const real_t xi_126 = xi_125 * -real_t(0.0416666666666667);
        const real_t xi_127 = xi_126 + xi_86 * -real_t(0.05);
        const real_t xi_128 =
            xi_120 + xi_124 + xi_127 + xi_89 * real_t(0.0142857142857143);
        const real_t xi_137 =
            xi_120 + xi_125 * real_t(0.0833333333333333) + xi_89 * -real_t(0.0357142857142857);
        const real_t xi_145 =
            rho * u_2 - vel2Term + xi_139 + xi_6 + xi_92 + xi_97 +
            _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        const real_t xi_146 = xi_118 * xi_145;
        const real_t xi_147 = xi_103 * real_t(0.0158730158730159) - xi_124 + xi_127 +
                              xi_89 * -real_t(0.0214285714285714);
        const real_t xi_150 = xi_123 * real_t(0.0625);
        const real_t xi_152 = -xi_149 + xi_150 + xi_151;
        const real_t xi_157 = xi_99 * real_t(0.0416666666666667);
        const real_t xi_158 = xi_125 * real_t(0.0208333333333333) + xi_157;
        const real_t xi_159 = -xi_156 + xi_158;
        const real_t xi_160 = xi_138 + xi_159;
        const real_t xi_161 = xi_149 + xi_150 + xi_151;
        const real_t xi_162 = xi_156 + xi_158;
        const real_t xi_163 = xi_136 + xi_162;
        const real_t xi_167 =
            xi_148 * (u_2 * xi_115 + xi_104 + xi_17 +
                      _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0]);
        const real_t xi_169 = xi_126 + xi_157 + xi_166 + xi_167 + xi_168;
        const real_t xi_175 = xi_145 * xi_153;
        const real_t xi_176 = xi_174 + xi_175;
        const real_t xi_177 = -xi_171 + xi_173 + xi_176;
        const real_t xi_182 = xi_126 + xi_157 + xi_166 - xi_167 + xi_168;
        const real_t xi_183 =
            xi_148 *
            (u_2 * xi_129 + xi_10 + xi_132 +
             _data_pdfs_21_318_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0]);
        const real_t xi_184 = -xi_150;
        const real_t xi_186 = -xi_183 + xi_184 + xi_185;
        const real_t xi_187 = xi_143 + xi_176;
        const real_t xi_190 = xi_159 - xi_188 + xi_189;
        const real_t xi_191 = xi_183 + xi_184 + xi_185;
        const real_t xi_192 = xi_162 + xi_188 - xi_189;
        const real_t xi_193 = xi_174 - xi_175;
        const real_t xi_194 = xi_171 - xi_173 + xi_193;
        const real_t xi_195 = xi_144 + xi_193;
        const real_t forceTerm_0 =
            xi_31 * -real_t(1.5) - xi_34 * xi_36 - xi_36 * xi_38 - xi_36 * xi_39;
        const real_t forceTerm_1 = xi_40 + xi_47;
        const real_t forceTerm_2 = -xi_40 + xi_47;
        const real_t forceTerm_3 = -xi_48 + xi_51;
        const real_t forceTerm_4 = xi_48 + xi_51;
        const real_t forceTerm_5 = xi_52 + xi_53;
        const real_t forceTerm_6 = -xi_52 + xi_53;
        const real_t forceTerm_7 = xi_58 + xi_60 + xi_67;
        const real_t forceTerm_8 = xi_57 + xi_67 + xi_68;
        const real_t forceTerm_9 = xi_57 + xi_60 + xi_70;
        const real_t forceTerm_10 = xi_58 + xi_68 + xi_70;
        const real_t forceTerm_11 = xi_66 + xi_71 + xi_75;
        const real_t forceTerm_12 = xi_69 + xi_75 + xi_76;
        const real_t forceTerm_13 = xi_60 + xi_78 + xi_80;
        const real_t forceTerm_14 = xi_68 + xi_77 + xi_80;
        const real_t forceTerm_15 = xi_66 + xi_76 + xi_82;
        const real_t forceTerm_16 = xi_69 + xi_71 + xi_82;
        const real_t forceTerm_17 = xi_60 + xi_77 + xi_83;
        const real_t forceTerm_18 = xi_68 + xi_78 + xi_83;
        _data_pdfs_tmp_20_30_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_0 + xi_103 * real_t(0.0238095238095238) + xi_86 * real_t(0.1) +
            xi_89 * real_t(0.0428571428571429) + xi_99 * -0.5) +
            _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_31_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_1 - xi_107 + xi_114 + xi_119 + xi_128) +
            _data_pdfs_20_31_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_32_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_2 + xi_107 + xi_113 - xi_119 + xi_128) +
            _data_pdfs_20_32_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_33_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_3 - xi_131 + xi_134 + xi_136 + xi_137) +
            _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_34_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_4 + xi_131 - xi_134 + xi_137 + xi_138) +
            _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_35_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_5 - xi_141 + xi_144 + xi_146 + xi_147) +
            _data_pdfs_2m1_35_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_36_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_6 + xi_141 + xi_143 - xi_146 + xi_147) +
            _data_pdfs_21_36_10[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_37_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_7 + xi_152 + xi_155 + xi_160) +
            _data_pdfs_20_37_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_38_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_8 + xi_155 + xi_161 + xi_163) +
            _data_pdfs_20_38_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_39_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_9 + xi_160 + xi_161 + xi_165) +
            _data_pdfs_20_39_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_310_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_10 + xi_152 + xi_163 + xi_165) +
            _data_pdfs_20_310_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_311_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_11 + xi_169 + xi_177 + xi_180) +
            _data_pdfs_2m1_311_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_312_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_12 + xi_177 + xi_181 + xi_182) +
            _data_pdfs_2m1_312_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_313_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_13 + xi_186 + xi_187 + xi_190) +
            _data_pdfs_2m1_313_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_314_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_14 + xi_187 + xi_191 + xi_192) +
            _data_pdfs_2m1_314_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_pdfs_tmp_20_315_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_15 + xi_180 + xi_182 + xi_194) +
            _data_pdfs_21_315_1m1[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_316_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_16 + xi_169 + xi_181 + xi_194) +
            _data_pdfs_21_316_11[_stride_pdfs_0 * ctr_0];
        _data_pdfs_tmp_20_317_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_17 + xi_190 + xi_191 + xi_195) +
            _data_pdfs_21_317_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_pdfs_tmp_20_318_10[_stride_pdfs_tmp_0 * ctr_0] =
            real_c(forceTerm_18 + xi_186 + xi_192 + xi_195) +
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
               int64_t const _stride_pdfs_3, real_t omega_bulk,
               real_t omega_even, real_t omega_odd, real_t omega_shear) {
  const real_t xi_35 = -real_c(omega_shear) + real_t(2.0);
  const real_t xi_36 = xi_35 * real_t(0.5);
  const real_t xi_41 = xi_35 * real_t(0.0833333333333333);
  const real_t xi_46 = xi_35 * real_t(0.166666666666667);
  const real_t xi_56 = xi_35 * real_t(0.25);
  const real_t xi_61 = xi_35 * real_t(0.0416666666666667);
  const real_t xi_106 = omega_odd * real_t(0.25);
  const real_t xi_112 = omega_odd * real_t(0.0833333333333333);
  const real_t xi_148 = omega_shear * real_t(0.25);
  const real_t xi_170 = omega_odd * real_t(0.0416666666666667);
  const real_t xi_172 = omega_odd * real_t(0.125);
  const int64_t rr_0 = 0;
  const real_t xi_118 = real_c(rr_0) * real_t(0.166666666666667);
  const real_t xi_153 = real_c(rr_0) * real_t(0.0833333333333333);
  for (int ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    real_t *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    real_t *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    real_t *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    real_t *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    real_t *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    real_t *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    real_t *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    real_t *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    for (int ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      real_t *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      real_t *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      real_t *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      real_t *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      real_t *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      real_t *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      real_t *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      real_t *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      real_t *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      real_t *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      real_t *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      real_t *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      real_t *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      real_t *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      real_t *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      real_t *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      real_t *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      real_t *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      real_t *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      real_t *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      real_t *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      real_t *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      for (int ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const real_t xi_196 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_197 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_198 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_199 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const real_t xi_200 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_201 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_202 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_203 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_204 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const real_t xi_205 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_206 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_207 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_208 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_209 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_210 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_211 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_212 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_213 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const real_t xi_214 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_215 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_216 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_217 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const real_t xi_0 = xi_212 + xi_215;
        const real_t xi_1 = xi_0 + xi_196;
        const real_t xi_2 = xi_198 + xi_205 + xi_210;
        const real_t xi_3 = xi_201 + xi_206;
        const real_t xi_4 = xi_197 + xi_207;
        const real_t xi_5 = xi_200 + xi_214;
        const real_t xi_6 = xi_202 + xi_211;
        const real_t xi_8 = -xi_207;
        const real_t xi_9 = -xi_208 + xi_8;
        const real_t xi_10 = -xi_202;
        const real_t xi_11 = -xi_203;
        const real_t xi_12 = -xi_197;
        const real_t xi_13 = xi_10 + xi_11 + xi_12;
        const real_t xi_14 = -xi_200;
        const real_t xi_15 = -xi_216;
        const real_t xi_16 = xi_14 + xi_15;
        const real_t xi_17 = -xi_214;
        const real_t xi_18 = -xi_201;
        const real_t xi_19 = xi_17 + xi_18;
        const real_t xi_20 = -xi_212;
        const real_t xi_21 = xi_10 + xi_20;
        const real_t xi_22 = -xi_210;
        const real_t xi_23 = -xi_211;
        const real_t xi_24 = xi_17 + xi_198 + xi_22 + xi_23;
        const real_t xi_40 = xi_213 * real_t(0.166666666666667);
        const real_t xi_48 = xi_199 * real_t(0.166666666666667);
        const real_t xi_52 = xi_204 * real_t(0.166666666666667);
        const real_t xi_55 = xi_213 * real_t(0.5);
        const real_t xi_59 = xi_199 * real_t(0.0833333333333333);
        const real_t xi_63 = xi_213 * real_t(0.0833333333333333);
        const real_t xi_73 = xi_204 * real_t(0.0833333333333333);
        const real_t xi_84 = -xi_217;
        const real_t xi_85 = xi_206 * real_t(3.0) + xi_211 * real_t(3.0) + xi_84;
        const real_t xi_86 =
            omega_even * (xi_198 * -real_t(3.0) + xi_200 * real_t(3.0) + xi_201 * -real_t(3.0) +
                          xi_205 * real_t(3.0) + xi_210 * -real_t(3.0) + xi_214 * -real_t(3.0) + xi_85);
        const real_t xi_87 =
            xi_198 * real_t(2.0) + xi_201 * real_t(2.0) + xi_210 * real_t(2.0) + xi_214 * real_t(2.0);
        const real_t xi_88 = xi_196 * real_t(5.0) + xi_197 * real_t(5.0) + xi_87;
        const real_t xi_89 =
            omega_even *
            (xi_200 * -real_t(2.0) + xi_202 * -real_t(5.0) + xi_203 * -real_t(5.0) + xi_205 * -real_t(2.0) +
             xi_212 * -real_t(5.0) + xi_215 * -real_t(5.0) + xi_85 + xi_88);
        const real_t xi_92 = -xi_198;
        const real_t xi_93 = xi_18 + xi_92;
        const real_t xi_94 = -xi_209;
        const real_t xi_97 = -xi_215;
        const real_t xi_98 = xi_11 + xi_15 + xi_21 + xi_97;
        const real_t xi_100 = xi_203 * real_t(2.0);
        const real_t xi_101 = xi_215 * real_t(2.0);
        const real_t xi_102 = xi_202 * real_t(2.0) + xi_212 * real_t(2.0);
        const real_t xi_103 =
            omega_even *
            (xi_100 + xi_101 + xi_102 + xi_200 * real_t(5.0) + xi_205 * real_t(5.0) +
             xi_206 * -real_t(4.0) + xi_207 * -real_t(7.0) + xi_208 * -real_t(7.0) + xi_209 * -real_t(7.0) +
             xi_211 * -real_t(4.0) + xi_216 * -real_t(7.0) + xi_84 + xi_88);
        const real_t xi_104 = xi_201 + xi_92;
        const real_t xi_105 = xi_104 + xi_14 + xi_205 + xi_214 + xi_22;
        const real_t xi_107 = xi_105 * xi_106;
        const real_t xi_108 = xi_208 * real_t(2.0);
        const real_t xi_109 = xi_216 * real_t(2.0);
        const real_t xi_110 = xi_207 * real_t(2.0) + xi_209 * -real_t(2.0);
        const real_t xi_111 = -xi_108 + xi_109 + xi_110 + xi_14 + xi_19 + xi_2;
        const real_t xi_113 = xi_111 * xi_112;
        const real_t xi_114 = -xi_113;
        const real_t xi_116 = xi_216 + xi_94;
        const real_t xi_120 = xi_103 * -real_t(0.0198412698412698);
        const real_t xi_121 = xi_202 + xi_203;
        const real_t xi_132 = xi_203 + xi_97;
        const real_t xi_133 = xi_12 + xi_132 + xi_196 + xi_20 + xi_202;
        const real_t xi_134 = xi_106 * xi_133;
        const real_t xi_135 = xi_1 + xi_108 - xi_109 + xi_110 + xi_13;
        const real_t xi_136 = xi_112 * xi_135;
        const real_t xi_138 = -xi_136;
        const real_t xi_139 = xi_210 + xi_214;
        const real_t xi_140 = xi_139 + xi_206 + xi_23 + xi_93;
        const real_t xi_141 = xi_106 * xi_140;
        const real_t xi_142 = -xi_100 - xi_101 + xi_102 + xi_24 + xi_3;
        const real_t xi_143 = xi_112 * xi_142;
        const real_t xi_144 = -xi_143;
        const real_t xi_151 = xi_103 * real_t(0.0138888888888889);
        const real_t xi_166 = xi_89 * -real_t(0.00714285714285714);
        const real_t xi_168 = xi_86 * real_t(0.025);
        const real_t xi_171 = xi_142 * xi_170;
        const real_t xi_173 = xi_140 * xi_172;
        const real_t xi_174 = xi_103 * -real_t(0.00396825396825397);
        const real_t xi_178 = xi_111 * xi_170;
        const real_t xi_179 = xi_105 * xi_172;
        const real_t xi_185 = xi_89 * real_t(0.0178571428571429);
        const real_t xi_188 = xi_133 * xi_172;
        const real_t xi_189 = xi_135 * xi_170;
        const real_t vel0Term = xi_1 + xi_209 + xi_216;
        const real_t vel1Term = xi_2 + xi_208;
        const real_t vel2Term = xi_203 + xi_3;
        const real_t rho =
            vel0Term + vel1Term + vel2Term + xi_217 + xi_4 + xi_5 + xi_6;
        const real_t xi_7 = 1 / (rho);
        const real_t u_0 = xi_7 * (vel0Term + xi_13 + xi_9);
        const real_t xi_25 = u_0 * xi_199;
        const real_t xi_26 = xi_25 * real_t(0.333333333333333);
        const real_t xi_32 = -xi_26;
        const real_t xi_90 = rho * (u_0 * u_0);
        const real_t xi_129 = rho * u_0;
        const real_t xi_130 = -vel0Term + xi_121 + xi_129 + xi_208 + xi_4;
        const real_t xi_131 = xi_118 * xi_130;
        const real_t xi_156 = xi_130 * xi_153;
        const real_t u_1 = xi_7 * (vel1Term + xi_16 + xi_19 + xi_209 + xi_8);
        const real_t xi_27 = u_1 * xi_213;
        const real_t xi_28 = xi_27 * real_t(0.333333333333333);
        const real_t xi_33 = -xi_28;
        const real_t xi_54 = u_1 * real_t(0.5);
        const real_t xi_57 = xi_56 * (u_0 * xi_55 + xi_199 * xi_54);
        const real_t xi_58 = -xi_57;
        const real_t xi_95 = rho * (u_1 * u_1);
        const real_t xi_96 = xi_9 + xi_94 + xi_95;
        const real_t xi_115 = rho * u_1;
        const real_t xi_117 =
            -vel1Term + xi_115 + xi_116 + xi_201 + xi_207 + xi_5;
        const real_t xi_119 = xi_117 * xi_118;
        const real_t xi_149 = xi_148 * (u_0 * xi_115 + xi_116 + xi_208 + xi_8);
        const real_t xi_154 = xi_117 * xi_153;
        const real_t xi_155 = xi_113 + xi_154;
        const real_t xi_164 = -xi_154;
        const real_t xi_165 = xi_114 + xi_164;
        const real_t xi_180 = xi_154 - xi_178 + xi_179;
        const real_t xi_181 = xi_164 + xi_178 - xi_179;
        const real_t u_2 = xi_7 * (vel2Term + xi_21 + xi_215 + xi_24);
        const real_t xi_29 = u_2 * xi_204;
        const real_t xi_30 = xi_29 * real_t(0.333333333333333);
        const real_t xi_31 = (-omega_bulk + real_t(2.0)) * (xi_26 + xi_28 + xi_30);
        const real_t xi_34 = xi_29 * real_t(0.666666666666667) + xi_32 + xi_33;
        const real_t xi_37 = -xi_30;
        const real_t xi_38 = xi_27 * real_t(0.666666666666667) + xi_32 + xi_37;
        const real_t xi_39 = xi_25 * real_t(0.666666666666667) + xi_33 + xi_37;
        const real_t xi_42 = xi_34 * xi_41;
        const real_t xi_43 = -xi_42;
        const real_t xi_44 = xi_39 * xi_41;
        const real_t xi_45 = -xi_44;
        const real_t xi_47 = xi_38 * xi_46 + xi_43 + xi_45;
        const real_t xi_49 = xi_38 * xi_41;
        const real_t xi_50 = -xi_49;
        const real_t xi_51 = xi_39 * xi_46 + xi_43 + xi_50;
        const real_t xi_53 = xi_34 * xi_46 + xi_45 + xi_50;
        const real_t xi_60 = xi_44 - xi_59;
        const real_t xi_62 = -xi_34 * xi_61;
        const real_t xi_64 = xi_31 * real_t(0.125);
        const real_t xi_65 = xi_49 + xi_64;
        const real_t xi_66 = xi_63 + xi_65;
        const real_t xi_67 = xi_62 + xi_66;
        const real_t xi_68 = xi_44 + xi_59;
        const real_t xi_69 = -xi_63 + xi_65;
        const real_t xi_70 = xi_62 + xi_69;
        const real_t xi_71 = xi_56 * (u_2 * xi_55 + xi_204 * xi_54);
        const real_t xi_72 = -xi_39 * xi_61;
        const real_t xi_74 = xi_42 + xi_73;
        const real_t xi_75 = xi_72 + xi_74;
        const real_t xi_76 = -xi_71;
        const real_t xi_77 = xi_56 * (u_0 * xi_204 * real_t(0.5) + u_2 * xi_199 * real_t(0.5));
        const real_t xi_78 = -xi_77;
        const real_t xi_79 = -xi_38 * xi_61;
        const real_t xi_80 = xi_64 + xi_74 + xi_79;
        const real_t xi_81 = xi_42 - xi_73;
        const real_t xi_82 = xi_72 + xi_81;
        const real_t xi_83 = xi_64 + xi_79 + xi_81;
        const real_t xi_91 = rho * (u_2 * u_2);
        const real_t xi_99 = omega_bulk * (xi_17 + xi_217 + xi_22 + xi_90 +
                                           xi_91 + xi_93 + xi_96 + xi_98);
        const real_t xi_122 = xi_206 + xi_211 - xi_91;
        const real_t xi_123 =
            omega_shear * (xi_0 + xi_121 + xi_122 + xi_16 - xi_205 + xi_96);
        const real_t xi_124 = xi_123 * real_t(0.125);
        const real_t xi_125 =
            omega_shear *
            (xi_122 + xi_196 * -real_t(2.0) + xi_197 * -real_t(2.0) + xi_200 + xi_205 + xi_87 +
             xi_9 + xi_90 * real_t(2.0) + xi_94 - xi_95 + xi_98);
        const real_t xi_126 = xi_125 * -real_t(0.0416666666666667);
        const real_t xi_127 = xi_126 + xi_86 * -real_t(0.05);
        const real_t xi_128 =
            xi_120 + xi_124 + xi_127 + xi_89 * real_t(0.0142857142857143);
        const real_t xi_137 =
            xi_120 + xi_125 * real_t(0.0833333333333333) + xi_89 * -real_t(0.0357142857142857);
        const real_t xi_145 =
            rho * u_2 - vel2Term + xi_139 + xi_212 + xi_6 + xi_92 + xi_97;
        const real_t xi_146 = xi_118 * xi_145;
        const real_t xi_147 = xi_103 * real_t(0.0158730158730159) - xi_124 + xi_127 +
                              xi_89 * -real_t(0.0214285714285714);
        const real_t xi_150 = xi_123 * real_t(0.0625);
        const real_t xi_152 = -xi_149 + xi_150 + xi_151;
        const real_t xi_157 = xi_99 * real_t(0.0416666666666667);
        const real_t xi_158 = xi_125 * real_t(0.0208333333333333) + xi_157;
        const real_t xi_159 = -xi_156 + xi_158;
        const real_t xi_160 = xi_138 + xi_159;
        const real_t xi_161 = xi_149 + xi_150 + xi_151;
        const real_t xi_162 = xi_156 + xi_158;
        const real_t xi_163 = xi_136 + xi_162;
        const real_t xi_167 = xi_148 * (u_2 * xi_115 + xi_104 + xi_17 + xi_210);
        const real_t xi_169 = xi_126 + xi_157 + xi_166 + xi_167 + xi_168;
        const real_t xi_175 = xi_145 * xi_153;
        const real_t xi_176 = xi_174 + xi_175;
        const real_t xi_177 = -xi_171 + xi_173 + xi_176;
        const real_t xi_182 = xi_126 + xi_157 + xi_166 - xi_167 + xi_168;
        const real_t xi_183 = xi_148 * (u_2 * xi_129 + xi_10 + xi_132 + xi_212);
        const real_t xi_184 = -xi_150;
        const real_t xi_186 = -xi_183 + xi_184 + xi_185;
        const real_t xi_187 = xi_143 + xi_176;
        const real_t xi_190 = xi_159 - xi_188 + xi_189;
        const real_t xi_191 = xi_183 + xi_184 + xi_185;
        const real_t xi_192 = xi_162 + xi_188 - xi_189;
        const real_t xi_193 = xi_174 - xi_175;
        const real_t xi_194 = xi_171 - xi_173 + xi_193;
        const real_t xi_195 = xi_144 + xi_193;
        const real_t forceTerm_0 =
            xi_31 * -real_t(1.5) - xi_34 * xi_36 - xi_36 * xi_38 - xi_36 * xi_39;
        const real_t forceTerm_1 = xi_40 + xi_47;
        const real_t forceTerm_2 = -xi_40 + xi_47;
        const real_t forceTerm_3 = -xi_48 + xi_51;
        const real_t forceTerm_4 = xi_48 + xi_51;
        const real_t forceTerm_5 = xi_52 + xi_53;
        const real_t forceTerm_6 = -xi_52 + xi_53;
        const real_t forceTerm_7 = xi_58 + xi_60 + xi_67;
        const real_t forceTerm_8 = xi_57 + xi_67 + xi_68;
        const real_t forceTerm_9 = xi_57 + xi_60 + xi_70;
        const real_t forceTerm_10 = xi_58 + xi_68 + xi_70;
        const real_t forceTerm_11 = xi_66 + xi_71 + xi_75;
        const real_t forceTerm_12 = xi_69 + xi_75 + xi_76;
        const real_t forceTerm_13 = xi_60 + xi_78 + xi_80;
        const real_t forceTerm_14 = xi_68 + xi_77 + xi_80;
        const real_t forceTerm_15 = xi_66 + xi_76 + xi_82;
        const real_t forceTerm_16 = xi_69 + xi_71 + xi_82;
        const real_t forceTerm_17 = xi_60 + xi_77 + xi_83;
        const real_t forceTerm_18 = xi_68 + xi_78 + xi_83;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_0 + xi_103 * real_t(0.0238095238095238) + xi_217 + xi_86 * real_t(0.1) +
            xi_89 * real_t(0.0428571428571429) + xi_99 * -0.5);
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_1 - xi_107 + xi_114 + xi_119 + xi_128 + xi_205);
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_2 + xi_107 + xi_113 - xi_119 + xi_128 + xi_200);
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_3 - xi_131 + xi_134 + xi_136 + xi_137 + xi_197);
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_4 + xi_131 - xi_134 + xi_137 + xi_138 + xi_196);
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_5 - xi_141 + xi_144 + xi_146 + xi_147 + xi_206);
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_6 + xi_141 + xi_143 - xi_146 + xi_147 + xi_211);
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_7 + xi_152 + xi_155 + xi_160 + xi_208);
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_8 + xi_155 + xi_161 + xi_163 + xi_209);
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_9 + xi_160 + xi_161 + xi_165 + xi_207);
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_10 + xi_152 + xi_163 + xi_165 + xi_216);
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_11 + xi_169 + xi_177 + xi_180 + xi_198);
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_12 + xi_177 + xi_181 + xi_182 + xi_201);
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_13 + xi_186 + xi_187 + xi_190 + xi_203);
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_14 + xi_187 + xi_191 + xi_192 + xi_215);
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_15 + xi_180 + xi_182 + xi_194 + xi_210);
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_16 + xi_169 + xi_181 + xi_194 + xi_214);
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_17 + xi_190 + xi_191 + xi_195 + xi_202);
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            real_c(forceTerm_18 + xi_186 + xi_192 + xi_195 + xi_212);
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

const real_t MRT_LatticeModel::w[19] = {
    real_t(0.333333333333333),  real_t(0.0555555555555556), real_t(0.0555555555555556),
    real_t(0.0555555555555556), real_t(0.0555555555555556), real_t(0.0555555555555556),
    real_t(0.0555555555555556), real_t(0.0277777777777778), real_t(0.0277777777777778),
    real_t(0.0277777777777778), real_t(0.0277777777777778), real_t(0.0277777777777778),
    real_t(0.0277777777777778), real_t(0.0277777777777778), real_t(0.0277777777777778),
    real_t(0.0277777777777778), real_t(0.0277777777777778), real_t(0.0277777777777778),
    real_t(0.0277777777777778)};
const real_t MRT_LatticeModel::wInv[19] = {
    real_t(3.00000000000000), real_t(18.0000000000000), real_t(18.0000000000000), real_t(18.0000000000000),
    real_t(18.0000000000000), real_t(18.0000000000000), real_t(18.0000000000000), real_t(36.0000000000000),
    real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000),
    real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000),
    real_t(36.0000000000000), real_t(36.0000000000000), real_t(36.0000000000000)};

void MRT_LatticeModel::Sweep::streamCollide(
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

  auto &lm =
      dynamic_cast<lbm::PdfField<MRT_LatticeModel> *>(pdfs)->latticeModel();
  WALBERLA_ASSERT_EQUAL(*(lm.blockId_), block->getId());

  auto &omega_even = lm.omega_even_;
  auto &omega_odd = lm.omega_odd_;
  auto &omega_bulk = lm.omega_bulk_;
  auto &omega_shear = lm.omega_shear_;
  auto &force = lm.force_;
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
      _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, omega_bulk, omega_even, omega_odd,
      omega_shear);
  pdfs->swapDataPointers(pdfs_tmp);
}

void MRT_LatticeModel::Sweep::collide(
    IBlock *block, const uint_t numberOfGhostLayersToInclude) {
  auto pdfs = block->getData<field::GhostLayerField<real_t, 19>>(pdfsID);

  auto &lm =
      dynamic_cast<lbm::PdfField<MRT_LatticeModel> *>(pdfs)->latticeModel();
  WALBERLA_ASSERT_EQUAL(*(lm.blockId_), block->getId());

  auto &omega_even = lm.omega_even_;
  auto &omega_odd = lm.omega_odd_;
  auto &omega_bulk = lm.omega_bulk_;
  auto &omega_shear = lm.omega_shear_;
  auto &force = lm.force_;
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
      omega_bulk, omega_even, omega_odd, omega_shear);
}

void MRT_LatticeModel::Sweep::stream(
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

mpi::SendBuffer &operator<<(mpi::SendBuffer &buf,
                            const ::walberla::lbm::MRT_LatticeModel &lm) {
  buf << lm.currentLevel;
  return buf;
}

mpi::RecvBuffer &operator>>(mpi::RecvBuffer &buf,
                            ::walberla::lbm::MRT_LatticeModel &lm) {
  buf >> lm.currentLevel;
  return buf;
}

} // namespace mpi
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
