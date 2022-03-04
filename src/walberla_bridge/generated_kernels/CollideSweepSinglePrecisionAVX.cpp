// kernel generated with pystencils v0.4.4, lbmpy v0.4.4,
// lbmpy_walberla/pystencils_walberla from commit
// 2527a5f799da52b4e5293a05ed691fc35bf7122b

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
//! \\file CollideSweepSinglePrecisionAVX.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecisionAVX.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include <immintrin.h>

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning push
#pragma warning(disable : 1599)
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_dfe735f2f0357dcc08d993f108917b07 {
static FUNC_PREFIX void
collidesweepsingleprecisionavx_collidesweepsingleprecisionavx(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, float omega_bulk, float omega_even,
    float omega_odd, float omega_shear) {
  const float xi_29 = omega_bulk * 0.50000000000000000f;
  const float xi_56 = omega_shear * 0.041666666666666667f;
  const float xi_61 = omega_bulk * 0.041666666666666667f;
  const float xi_72 = omega_shear * 0.12500000000000000f;
  const float xi_128 = omega_odd * 0.25000000000000000f;
  const float xi_134 = omega_odd * 0.083333333333333333f;
  const float xi_171 = omega_shear * 0.25000000000000000f;
  const float xi_194 = omega_odd * 0.041666666666666667f;
  const float xi_196 = omega_odd * 0.12500000000000000f;
  const float rr_0 = 0.0f;
  const float xi_54 = rr_0 * 0.041666666666666667f;
  const float xi_140 = rr_0 * 0.16666666666666667f;
  const float xi_176 = rr_0 * 0.083333333333333333f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      {
        for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_force_0) / (8)) * (8);
             ctr_0 += 8) {
          const __m256 xi_220 = _mm256_load_ps(&_data_pdfs_20_33_10[ctr_0]);
          const __m256 xi_221 = _mm256_load_ps(&_data_pdfs_20_37_10[ctr_0]);
          const __m256 xi_222 = _mm256_load_ps(&_data_pdfs_20_317_10[ctr_0]);
          const __m256 xi_223 = _mm256_load_ps(&_data_pdfs_20_34_10[ctr_0]);
          const __m256 xi_224 = _mm256_load_ps(&_data_pdfs_20_310_10[ctr_0]);
          const __m256 xi_225 = _mm256_load_ps(&_data_pdfs_20_39_10[ctr_0]);
          const __m256 xi_226 = _mm256_load_ps(&_data_force_20_30_10[ctr_0]);
          const __m256 xi_227 = _mm256_load_ps(&_data_pdfs_20_30_10[ctr_0]);
          const __m256 xi_228 = _mm256_load_ps(&_data_pdfs_20_312_10[ctr_0]);
          const __m256 xi_229 = _mm256_load_ps(&_data_pdfs_20_36_10[ctr_0]);
          const __m256 xi_230 = _mm256_load_ps(&_data_pdfs_20_318_10[ctr_0]);
          const __m256 xi_231 = _mm256_load_ps(&_data_force_20_31_10[ctr_0]);
          const __m256 xi_232 = _mm256_load_ps(&_data_pdfs_20_313_10[ctr_0]);
          const __m256 xi_233 = _mm256_load_ps(&_data_pdfs_20_35_10[ctr_0]);
          const __m256 xi_234 = _mm256_load_ps(&_data_pdfs_20_38_10[ctr_0]);
          const __m256 xi_235 = _mm256_load_ps(&_data_pdfs_20_315_10[ctr_0]);
          const __m256 xi_236 = _mm256_load_ps(&_data_pdfs_20_314_10[ctr_0]);
          const __m256 xi_237 = _mm256_load_ps(&_data_pdfs_20_311_10[ctr_0]);
          const __m256 xi_238 = _mm256_load_ps(&_data_pdfs_20_316_10[ctr_0]);
          const __m256 xi_239 = _mm256_load_ps(&_data_force_20_32_10[ctr_0]);
          const __m256 xi_240 = _mm256_load_ps(&_data_pdfs_20_32_10[ctr_0]);
          const __m256 xi_241 = _mm256_load_ps(&_data_pdfs_20_31_10[ctr_0]);
          const __m256 xi_0 = _mm256_add_ps(xi_230, xi_236);
          const __m256 xi_1 = _mm256_add_ps(xi_0, xi_223);
          const __m256 xi_2 =
              _mm256_add_ps(_mm256_add_ps(xi_235, xi_237), xi_241);
          const __m256 xi_3 = _mm256_add_ps(xi_228, xi_233);
          const __m256 xi_4 = _mm256_add_ps(xi_220, xi_225);
          const __m256 xi_5 = _mm256_add_ps(xi_238, xi_240);
          const __m256 xi_6 = _mm256_add_ps(xi_222, xi_229);
          const __m256 xi_9 =
              _mm256_mul_ps(xi_225, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_10 = _mm256_add_ps(
              _mm256_mul_ps(xi_221, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
              xi_9);
          const __m256 xi_11 =
              _mm256_mul_ps(xi_222, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_12 =
              _mm256_mul_ps(xi_232, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_13 =
              _mm256_mul_ps(xi_220, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_14 =
              _mm256_add_ps(_mm256_add_ps(xi_11, xi_12), xi_13);
          const __m256 xi_15 =
              _mm256_mul_ps(xi_240, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_16 =
              _mm256_mul_ps(xi_224, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_17 = _mm256_add_ps(xi_15, xi_16);
          const __m256 xi_18 =
              _mm256_mul_ps(xi_238, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_19 =
              _mm256_mul_ps(xi_228, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_20 = _mm256_add_ps(xi_18, xi_19);
          const __m256 xi_21 =
              _mm256_mul_ps(xi_230, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_22 = _mm256_add_ps(xi_11, xi_21);
          const __m256 xi_23 =
              _mm256_mul_ps(xi_235, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_24 =
              _mm256_mul_ps(xi_229, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_25 = _mm256_add_ps(
              _mm256_add_ps(_mm256_add_ps(xi_18, xi_23), xi_237), xi_24);
          const __m256 xi_30 = _mm256_mul_ps(
              xi_231,
              _mm256_set_ps(0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f));
          const __m256 xi_31 = _mm256_mul_ps(
              xi_231,
              _mm256_set_ps(0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f));
          const __m256 xi_43 = _mm256_mul_ps(
              xi_226,
              _mm256_set_ps(0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f));
          const __m256 xi_44 = _mm256_mul_ps(
              xi_226,
              _mm256_set_ps(0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f));
          const __m256 xi_50 = _mm256_mul_ps(
              xi_239,
              _mm256_set_ps(0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f,
                            0.16666666666666667f, 0.16666666666666667f));
          const __m256 xi_51 = _mm256_mul_ps(
              xi_239,
              _mm256_set_ps(0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f));
          const __m256 xi_68 = _mm256_mul_ps(
              xi_231,
              _mm256_set_ps(0.25000000000000000f, 0.25000000000000000f,
                            0.25000000000000000f, 0.25000000000000000f,
                            0.25000000000000000f, 0.25000000000000000f,
                            0.25000000000000000f, 0.25000000000000000f));
          const __m256 xi_73 =
              _mm256_mul_ps(xi_231, _mm256_set_ps(xi_72, xi_72, xi_72, xi_72,
                                                  xi_72, xi_72, xi_72, xi_72));
          const __m256 xi_106 =
              _mm256_mul_ps(xi_227, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_107 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_229, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f,
                                                      3.0f, 3.0f, 3.0f, 3.0f)),
                  _mm256_mul_ps(xi_233, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f,
                                                      3.0f, 3.0f, 3.0f, 3.0f))),
              xi_106);
          const __m256 xi_108 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_mul_ps(
                                          xi_228,
                                          _mm256_set_ps(-3.0f, -3.0f, -3.0f,
                                                        -3.0f, -3.0f, -3.0f,
                                                        -3.0f, -3.0f)),
                                      _mm256_mul_ps(
                                          xi_235,
                                          _mm256_set_ps(-3.0f, -3.0f, -3.0f,
                                                        -3.0f, -3.0f, -3.0f,
                                                        -3.0f, -3.0f))),
                                  _mm256_mul_ps(
                                      xi_237, _mm256_set_ps(-3.0f, -3.0f, -3.0f,
                                                            -3.0f, -3.0f, -3.0f,
                                                            -3.0f, -3.0f))),
                              _mm256_mul_ps(xi_238,
                                            _mm256_set_ps(-3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f))),
                          _mm256_mul_ps(xi_240,
                                        _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f,
                                                      3.0f, 3.0f, 3.0f, 3.0f))),
                      _mm256_mul_ps(xi_241,
                                    _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                                                  3.0f, 3.0f, 3.0f))),
                  xi_107),
              _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                            omega_even, omega_even, omega_even, omega_even));
          const __m256 xi_109 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_228,
                                    _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f)),
                      _mm256_mul_ps(xi_235,
                                    _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f))),
                  _mm256_mul_ps(xi_237, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f))),
              _mm256_mul_ps(xi_238, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f)));
          const __m256 xi_110 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_220, _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f,
                                                      5.0f, 5.0f, 5.0f, 5.0f)),
                  _mm256_mul_ps(xi_223, _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f,
                                                      5.0f, 5.0f, 5.0f, 5.0f))),
              xi_109);
          const __m256 xi_111 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              xi_222,
                                              _mm256_set_ps(-5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f)),
                                          _mm256_mul_ps(
                                              xi_230,
                                              _mm256_set_ps(-5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f))),
                                      _mm256_mul_ps(
                                          xi_232,
                                          _mm256_set_ps(-5.0f, -5.0f, -5.0f,
                                                        -5.0f, -5.0f, -5.0f,
                                                        -5.0f, -5.0f))),
                                  _mm256_mul_ps(
                                      xi_236, _mm256_set_ps(-5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f))),
                              _mm256_mul_ps(xi_240,
                                            _mm256_set_ps(-2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f))),
                          _mm256_mul_ps(xi_241,
                                        _mm256_set_ps(-2.0f, -2.0f, -2.0f,
                                                      -2.0f, -2.0f, -2.0f,
                                                      -2.0f, -2.0f))),
                      xi_107),
                  xi_110),
              _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                            omega_even, omega_even, omega_even, omega_even));
          const __m256 xi_114 =
              _mm256_mul_ps(xi_237, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_115 = _mm256_add_ps(xi_114, xi_19);
          const __m256 xi_116 =
              _mm256_mul_ps(xi_234, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_119 =
              _mm256_mul_ps(xi_236, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_120 = _mm256_add_ps(
              _mm256_add_ps(_mm256_add_ps(xi_119, xi_12), xi_16), xi_22);
          const __m256 xi_122 =
              _mm256_mul_ps(xi_232, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f));
          const __m256 xi_123 =
              _mm256_mul_ps(xi_236, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f));
          const __m256 xi_124 = _mm256_add_ps(
              _mm256_mul_ps(xi_222, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f)),
              _mm256_mul_ps(xi_230, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f)));
          const __m256 xi_125 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_mul_ps(
                                                                  xi_221,
                                                                  _mm256_set_ps(
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f)),
                                                              _mm256_mul_ps(
                                                                  xi_224,
                                                                  _mm256_set_ps(
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f,
                                                                      -7.0f))),
                                                          _mm256_mul_ps(
                                                              xi_225,
                                                              _mm256_set_ps(
                                                                  -7.0f, -7.0f,
                                                                  -7.0f, -7.0f,
                                                                  -7.0f, -7.0f,
                                                                  -7.0f,
                                                                  -7.0f))),
                                                      _mm256_mul_ps(
                                                          xi_229,
                                                          _mm256_set_ps(
                                                              -4.0f, -4.0f,
                                                              -4.0f, -4.0f,
                                                              -4.0f, -4.0f,
                                                              -4.0f, -4.0f))),
                                                  _mm256_mul_ps(
                                                      xi_233,
                                                      _mm256_set_ps(
                                                          -4.0f, -4.0f, -4.0f,
                                                          -4.0f, -4.0f, -4.0f,
                                                          -4.0f, -4.0f))),
                                              _mm256_mul_ps(
                                                  xi_234,
                                                  _mm256_set_ps(-7.0f, -7.0f,
                                                                -7.0f, -7.0f,
                                                                -7.0f, -7.0f,
                                                                -7.0f, -7.0f))),
                                          _mm256_mul_ps(
                                              xi_240,
                                              _mm256_set_ps(5.0f, 5.0f, 5.0f,
                                                            5.0f, 5.0f, 5.0f,
                                                            5.0f, 5.0f))),
                                      _mm256_mul_ps(
                                          xi_241, _mm256_set_ps(
                                                      5.0f, 5.0f, 5.0f, 5.0f,
                                                      5.0f, 5.0f, 5.0f, 5.0f))),
                                  xi_106),
                              xi_110),
                          xi_122),
                      xi_123),
                  xi_124),
              _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                            omega_even, omega_even, omega_even, omega_even));
          const __m256 xi_126 = _mm256_add_ps(xi_114, xi_228);
          const __m256 xi_127 = _mm256_add_ps(
              _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_126, xi_15), xi_23),
                            xi_238),
              xi_241);
          const __m256 xi_129 = _mm256_mul_ps(
              xi_127, _mm256_set_ps(xi_128, xi_128, xi_128, xi_128, xi_128,
                                    xi_128, xi_128, xi_128));
          const __m256 xi_130 =
              _mm256_mul_ps(xi_221, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f));
          const __m256 xi_131 =
              _mm256_mul_ps(xi_224, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f));
          const __m256 xi_132 = _mm256_add_ps(
              _mm256_mul_ps(xi_225, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                  2.0f, 2.0f, 2.0f)),
              _mm256_mul_ps(xi_234, _mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f,
                                                  -2.0f, -2.0f, -2.0f, -2.0f)));
          const __m256 xi_133 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(xi_130,
                                            _mm256_set_ps(-1.0, -1.0, -1.0,
                                                          -1.0, -1.0, -1.0,
                                                          -1.0, -1.0)),
                              xi_131),
                          xi_132),
                      xi_15),
                  xi_2),
              xi_20);
          const __m256 xi_135 = _mm256_mul_ps(
              xi_133, _mm256_set_ps(xi_134, xi_134, xi_134, xi_134, xi_134,
                                    xi_134, xi_134, xi_134));
          const __m256 xi_136 =
              _mm256_mul_ps(xi_135, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_138 = _mm256_add_ps(xi_116, xi_224);
          const __m256 xi_142 = _mm256_add_ps(xi_222, xi_232);
          const __m256 xi_146 = _mm256_mul_ps(
              xi_125,
              _mm256_set_ps(-0.019841269841269841f, -0.019841269841269841f,
                            -0.019841269841269841f, -0.019841269841269841f,
                            -0.019841269841269841f, -0.019841269841269841f,
                            -0.019841269841269841f, -0.019841269841269841f));
          const __m256 xi_154 = _mm256_add_ps(xi_119, xi_232);
          const __m256 xi_155 = _mm256_add_ps(
              _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_13, xi_154), xi_21),
                            xi_222),
              xi_223);
          const __m256 xi_156 = _mm256_mul_ps(
              xi_155, _mm256_set_ps(xi_128, xi_128, xi_128, xi_128, xi_128,
                                    xi_128, xi_128, xi_128));
          const __m256 xi_157 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(xi_131,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                          xi_1),
                      xi_130),
                  xi_132),
              xi_14);
          const __m256 xi_158 = _mm256_mul_ps(
              xi_157, _mm256_set_ps(xi_134, xi_134, xi_134, xi_134, xi_134,
                                    xi_134, xi_134, xi_134));
          const __m256 xi_160 =
              _mm256_mul_ps(xi_158, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_161 = _mm256_add_ps(xi_235, xi_238);
          const __m256 xi_162 = _mm256_add_ps(
              _mm256_add_ps(_mm256_add_ps(xi_115, xi_161), xi_233), xi_24);
          const __m256 xi_163 = _mm256_mul_ps(
              xi_162, _mm256_set_ps(xi_128, xi_128, xi_128, xi_128, xi_128,
                                    xi_128, xi_128, xi_128));
          const __m256 xi_166 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(xi_122,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                          _mm256_mul_ps(xi_123,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                      xi_124),
                  xi_25),
              xi_3);
          const __m256 xi_167 = _mm256_mul_ps(
              xi_166, _mm256_set_ps(xi_134, xi_134, xi_134, xi_134, xi_134,
                                    xi_134, xi_134, xi_134));
          const __m256 xi_168 =
              _mm256_mul_ps(xi_167, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_170 = xi_167;
          const __m256 xi_174 = _mm256_mul_ps(
              xi_125,
              _mm256_set_ps(0.013888888888888889f, 0.013888888888888889f,
                            0.013888888888888889f, 0.013888888888888889f,
                            0.013888888888888889f, 0.013888888888888889f,
                            0.013888888888888889f, 0.013888888888888889f));
          const __m256 xi_190 = _mm256_mul_ps(
              xi_111,
              _mm256_set_ps(-0.0071428571428571429f, -0.0071428571428571429f,
                            -0.0071428571428571429f, -0.0071428571428571429f,
                            -0.0071428571428571429f, -0.0071428571428571429f,
                            -0.0071428571428571429f, -0.0071428571428571429f));
          const __m256 xi_192 = _mm256_mul_ps(
              xi_108,
              _mm256_set_ps(0.025000000000000000f, 0.025000000000000000f,
                            0.025000000000000000f, 0.025000000000000000f,
                            0.025000000000000000f, 0.025000000000000000f,
                            0.025000000000000000f, 0.025000000000000000f));
          const __m256 xi_195 = _mm256_mul_ps(
              xi_166, _mm256_set_ps(xi_194, xi_194, xi_194, xi_194, xi_194,
                                    xi_194, xi_194, xi_194));
          const __m256 xi_197 = _mm256_mul_ps(
              xi_162, _mm256_set_ps(xi_196, xi_196, xi_196, xi_196, xi_196,
                                    xi_196, xi_196, xi_196));
          const __m256 xi_198 = _mm256_mul_ps(
              xi_125,
              _mm256_set_ps(-0.0039682539682539683f, -0.0039682539682539683f,
                            -0.0039682539682539683f, -0.0039682539682539683f,
                            -0.0039682539682539683f, -0.0039682539682539683f,
                            -0.0039682539682539683f, -0.0039682539682539683f));
          const __m256 xi_202 = _mm256_mul_ps(
              xi_133, _mm256_set_ps(xi_194, xi_194, xi_194, xi_194, xi_194,
                                    xi_194, xi_194, xi_194));
          const __m256 xi_203 = _mm256_mul_ps(
              xi_127, _mm256_set_ps(xi_196, xi_196, xi_196, xi_196, xi_196,
                                    xi_196, xi_196, xi_196));
          const __m256 xi_209 = _mm256_mul_ps(
              xi_111,
              _mm256_set_ps(0.017857142857142857f, 0.017857142857142857f,
                            0.017857142857142857f, 0.017857142857142857f,
                            0.017857142857142857f, 0.017857142857142857f,
                            0.017857142857142857f, 0.017857142857142857f));
          const __m256 xi_212 = _mm256_mul_ps(
              xi_155, _mm256_set_ps(xi_196, xi_196, xi_196, xi_196, xi_196,
                                    xi_196, xi_196, xi_196));
          const __m256 xi_213 = _mm256_mul_ps(
              xi_157, _mm256_set_ps(xi_194, xi_194, xi_194, xi_194, xi_194,
                                    xi_194, xi_194, xi_194));
          const __m256 xi_32 =
              _mm256_mul_ps(xi_31, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                 rr_0, rr_0, rr_0));
          const __m256 xi_45 =
              _mm256_mul_ps(xi_44, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                 rr_0, rr_0, rr_0));
          const __m256 xi_52 =
              _mm256_mul_ps(xi_51, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                                 rr_0, rr_0, rr_0));
          const __m256 xi_55 =
              _mm256_mul_ps(xi_226, _mm256_set_ps(xi_54, xi_54, xi_54, xi_54,
                                                  xi_54, xi_54, xi_54, xi_54));
          const __m256 xi_60 =
              _mm256_mul_ps(xi_231, _mm256_set_ps(xi_54, xi_54, xi_54, xi_54,
                                                  xi_54, xi_54, xi_54, xi_54));
          const __m256 xi_82 =
              _mm256_mul_ps(xi_239, _mm256_set_ps(xi_54, xi_54, xi_54, xi_54,
                                                  xi_54, xi_54, xi_54, xi_54));
          const __m256 vel0Term =
              _mm256_add_ps(_mm256_add_ps(xi_1, xi_224), xi_234);
          const __m256 vel1Term = _mm256_add_ps(xi_2, xi_221);
          const __m256 vel2Term = _mm256_add_ps(xi_232, xi_3);
          const __m256 rho = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(_mm256_add_ps(vel0Term, vel1Term),
                                        vel2Term),
                          xi_227),
                      xi_4),
                  xi_5),
              xi_6);
          const __m256 xi_7 = _mm256_div_ps(
              _mm256_set_ps(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), rho);
          const __m256 xi_8 = _mm256_mul_ps(
              xi_7, _mm256_set_ps(0.50000000000000000f, 0.50000000000000000f,
                                  0.50000000000000000f, 0.50000000000000000f,
                                  0.50000000000000000f, 0.50000000000000000f,
                                  0.50000000000000000f, 0.50000000000000000f));
          const __m256 u_0 = _mm256_add_ps(
              _mm256_mul_ps(
                  xi_7, _mm256_add_ps(_mm256_add_ps(vel0Term, xi_10), xi_14)),
              _mm256_mul_ps(xi_226, xi_8));
          const __m256 xi_26 = _mm256_mul_ps(u_0, xi_226);
          const __m256 xi_38 = _mm256_mul_ps(
              xi_26, _mm256_set_ps(0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f));
          const __m256 xi_39 = _mm256_mul_ps(
              xi_26,
              _mm256_set_ps(0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f));
          const __m256 xi_40 = _mm256_mul_ps(
              xi_39, _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                   omega_shear, omega_shear, omega_shear,
                                   omega_shear, omega_shear));
          const __m256 xi_41 = _mm256_add_ps(
              _mm256_mul_ps(xi_38, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              xi_40);
          const __m256 xi_57 = _mm256_add_ps(
              _mm256_mul_ps(
                  _mm256_mul_ps(xi_26, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  _mm256_set_ps(xi_56, xi_56, xi_56, xi_56, xi_56, xi_56, xi_56,
                                xi_56)),
              xi_38);
          const __m256 xi_58 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_44, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_55),
              xi_57);
          const __m256 xi_62 = _mm256_mul_ps(
              _mm256_mul_ps(xi_26, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              _mm256_set_ps(xi_61, xi_61, xi_61, xi_61, xi_61, xi_61, xi_61,
                            xi_61));
          const __m256 xi_69 = _mm256_mul_ps(u_0, xi_68);
          const __m256 xi_74 = _mm256_mul_ps(u_0, xi_73);
          const __m256 xi_78 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_55, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_44),
              xi_57);
          const __m256 xi_85 =
              _mm256_mul_ps(xi_39, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0));
          const __m256 xi_96 = _mm256_mul_ps(u_0, xi_239);
          const __m256 xi_97 = _mm256_mul_ps(
              xi_96, _mm256_set_ps(0.25000000000000000f, 0.25000000000000000f,
                                   0.25000000000000000f, 0.25000000000000000f,
                                   0.25000000000000000f, 0.25000000000000000f,
                                   0.25000000000000000f, 0.25000000000000000f));
          const __m256 xi_100 =
              _mm256_mul_ps(xi_96, _mm256_set_ps(xi_72, xi_72, xi_72, xi_72,
                                                 xi_72, xi_72, xi_72, xi_72));
          const __m256 xi_112 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_0, u_0)));
          const __m256 xi_151 = _mm256_mul_ps(rho, u_0);
          const __m256 xi_152 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(vel0Term,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                          xi_142),
                      xi_151),
                  xi_221),
              xi_4);
          const __m256 xi_153 = _mm256_mul_ps(
              xi_152, _mm256_set_ps(xi_140, xi_140, xi_140, xi_140, xi_140,
                                    xi_140, xi_140, xi_140));
          const __m256 xi_180 = _mm256_mul_ps(
              xi_152, _mm256_set_ps(xi_176, xi_176, xi_176, xi_176, xi_176,
                                    xi_176, xi_176, xi_176));
          const __m256 u_1 = _mm256_add_ps(
              _mm256_mul_ps(
                  xi_7,
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(_mm256_add_ps(vel1Term, xi_17), xi_20),
                          xi_234),
                      xi_9)),
              _mm256_mul_ps(xi_231, xi_8));
          const __m256 xi_27 = _mm256_mul_ps(u_1, xi_231);
          const __m256 xi_33 = _mm256_mul_ps(
              xi_27, _mm256_set_ps(0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f));
          const __m256 xi_46 = _mm256_mul_ps(
              xi_27,
              _mm256_set_ps(0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f));
          const __m256 xi_47 = _mm256_mul_ps(
              xi_46, _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                   omega_shear, omega_shear, omega_shear,
                                   omega_shear, omega_shear));
          const __m256 xi_48 = _mm256_add_ps(
              _mm256_mul_ps(xi_33, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              xi_47);
          const __m256 xi_63 = _mm256_mul_ps(
              _mm256_mul_ps(xi_27, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              _mm256_set_ps(xi_61, xi_61, xi_61, xi_61, xi_61, xi_61, xi_61,
                            xi_61));
          const __m256 xi_70 = _mm256_mul_ps(
              u_1, _mm256_set_ps(0.25000000000000000f, 0.25000000000000000f,
                                 0.25000000000000000f, 0.25000000000000000f,
                                 0.25000000000000000f, 0.25000000000000000f,
                                 0.25000000000000000f, 0.25000000000000000f));
          const __m256 xi_71 = _mm256_mul_ps(xi_226, xi_70);
          const __m256 xi_75 =
              _mm256_mul_ps(u_1, _mm256_set_ps(xi_72, xi_72, xi_72, xi_72,
                                               xi_72, xi_72, xi_72, xi_72));
          const __m256 xi_76 = _mm256_mul_ps(xi_226, xi_75);
          const __m256 xi_77 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_69,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
                      _mm256_mul_ps(xi_71,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0))),
                  xi_74),
              xi_76);
          const __m256 xi_79 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_74,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
                      _mm256_mul_ps(xi_76,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0))),
                  xi_69),
              xi_71);
          const __m256 xi_87 = _mm256_mul_ps(xi_239, xi_70);
          const __m256 xi_89 = _mm256_mul_ps(xi_239, xi_75);
          const __m256 xi_94 =
              _mm256_mul_ps(xi_46, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0));
          const __m256 xi_117 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_1, u_1)));
          const __m256 xi_118 =
              _mm256_add_ps(_mm256_add_ps(xi_10, xi_116), xi_117);
          const __m256 xi_137 = _mm256_mul_ps(rho, u_1);
          const __m256 xi_139 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(vel1Term,
                                            _mm256_set_ps(-1.0, -1.0, -1.0,
                                                          -1.0, -1.0, -1.0,
                                                          -1.0, -1.0)),
                              xi_137),
                          xi_138),
                      xi_225),
                  xi_228),
              xi_5);
          const __m256 xi_141 = _mm256_mul_ps(
              xi_139, _mm256_set_ps(xi_140, xi_140, xi_140, xi_140, xi_140,
                                    xi_140, xi_140, xi_140));
          const __m256 xi_172 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_mul_ps(u_0, xi_137), xi_138),
                      xi_221),
                  xi_9),
              _mm256_set_ps(xi_171, xi_171, xi_171, xi_171, xi_171, xi_171,
                            xi_171, xi_171));
          const __m256 xi_177 = _mm256_mul_ps(
              xi_139, _mm256_set_ps(xi_176, xi_176, xi_176, xi_176, xi_176,
                                    xi_176, xi_176, xi_176));
          const __m256 xi_178 = xi_177;
          const __m256 xi_179 = _mm256_add_ps(xi_135, xi_178);
          const __m256 xi_188 =
              _mm256_mul_ps(xi_177, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_189 = _mm256_add_ps(xi_136, xi_188);
          const __m256 xi_204 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_202, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_178),
              xi_203);
          const __m256 xi_205 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_203, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_188),
              xi_202);
          const __m256 u_2 = _mm256_add_ps(
              _mm256_mul_ps(
                  xi_7,
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(vel2Term, xi_22), xi_236),
                      xi_25)),
              _mm256_mul_ps(xi_239, xi_8));
          const __m256 xi_28 = _mm256_mul_ps(u_2, xi_239);
          const __m256 xi_34 = _mm256_mul_ps(
              xi_28, _mm256_set_ps(0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f,
                                   0.16666666666666667f, 0.16666666666666667f));
          const __m256 xi_35 = _mm256_mul_ps(
              xi_28,
              _mm256_set_ps(0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f,
                            0.083333333333333333f, 0.083333333333333333f));
          const __m256 xi_36 = _mm256_mul_ps(
              xi_35, _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                   omega_shear, omega_shear, omega_shear,
                                   omega_shear, omega_shear));
          const __m256 xi_37 = _mm256_add_ps(
              _mm256_mul_ps(xi_34, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              xi_36);
          const __m256 xi_42 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          xi_27,
                          _mm256_set_ps(
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_33,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_37),
              xi_41);
          const __m256 xi_49 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          xi_26,
                          _mm256_set_ps(
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_38,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_37),
              xi_48);
          const __m256 xi_53 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          xi_28,
                          _mm256_set_ps(
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f,
                              0.33333333333333333f, 0.33333333333333333f)),
                      _mm256_mul_ps(
                          _mm256_mul_ps(xi_34,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_41),
              xi_48);
          const __m256 xi_59 =
              _mm256_mul_ps(xi_35, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0));
          const __m256 xi_64 = _mm256_mul_ps(
              _mm256_mul_ps(xi_28, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                 -1.0, -1.0, -1.0)),
              _mm256_set_ps(xi_61, xi_61, xi_61, xi_61, xi_61, xi_61, xi_61,
                            xi_61));
          const __m256 xi_65 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(
                              _mm256_mul_ps(
                                  xi_27, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                              _mm256_set_ps(xi_56, xi_56, xi_56, xi_56, xi_56,
                                            xi_56, xi_56, xi_56)),
                          xi_33),
                      xi_62),
                  xi_63),
              xi_64);
          const __m256 xi_66 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_60, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_31),
              xi_65);
          const __m256 xi_67 =
              _mm256_add_ps(_mm256_add_ps(xi_36, xi_59), xi_66);
          const __m256 xi_80 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_31, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_60),
              xi_65);
          const __m256 xi_81 =
              _mm256_add_ps(_mm256_add_ps(xi_36, xi_59), xi_80);
          const __m256 xi_83 = _mm256_add_ps(
              _mm256_mul_ps(
                  _mm256_mul_ps(xi_28, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  _mm256_set_ps(xi_56, xi_56, xi_56, xi_56, xi_56, xi_56, xi_56,
                                xi_56)),
              xi_34);
          const __m256 xi_84 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_82, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_51),
              xi_83);
          const __m256 xi_86 =
              _mm256_add_ps(_mm256_add_ps(xi_40, xi_66), xi_85);
          const __m256 xi_88 = _mm256_mul_ps(u_2, xi_68);
          const __m256 xi_90 = _mm256_mul_ps(u_2, xi_73);
          const __m256 xi_91 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_89,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
                      _mm256_mul_ps(xi_90,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0))),
                  xi_87),
              xi_88);
          const __m256 xi_92 =
              _mm256_add_ps(_mm256_add_ps(xi_40, xi_80), xi_85);
          const __m256 xi_93 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_87,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
                      _mm256_mul_ps(xi_88,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0))),
                  xi_89),
              xi_90);
          const __m256 xi_95 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(xi_47, xi_62), xi_63), xi_64),
                  xi_84),
              xi_94);
          const __m256 xi_98 = _mm256_mul_ps(u_2, xi_226);
          const __m256 xi_99 = _mm256_mul_ps(
              xi_98, _mm256_set_ps(0.25000000000000000f, 0.25000000000000000f,
                                   0.25000000000000000f, 0.25000000000000000f,
                                   0.25000000000000000f, 0.25000000000000000f,
                                   0.25000000000000000f, 0.25000000000000000f));
          const __m256 xi_101 =
              _mm256_mul_ps(xi_98, _mm256_set_ps(xi_72, xi_72, xi_72, xi_72,
                                                 xi_72, xi_72, xi_72, xi_72));
          const __m256 xi_102 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_97,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
                      _mm256_mul_ps(xi_99,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0))),
                  xi_100),
              xi_101);
          const __m256 xi_103 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(xi_100,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
                      _mm256_mul_ps(xi_101,
                                    _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0))),
                  xi_97),
              xi_99);
          const __m256 xi_104 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_51, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_82),
              xi_83);
          const __m256 xi_105 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(xi_104, xi_47), xi_62),
                      xi_63),
                  xi_64),
              xi_94);
          const __m256 xi_113 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_2, u_2)));
          const __m256 xi_121 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(_mm256_add_ps(xi_112, xi_113),
                                                xi_115),
                                  xi_118),
                              xi_120),
                          xi_18),
                      xi_227),
                  xi_23),
              _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk,
                            omega_bulk, omega_bulk, omega_bulk, omega_bulk));
          const __m256 xi_143 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_113, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_229),
              xi_233);
          const __m256 xi_144 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(xi_241,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0)),
                                  xi_0),
                              xi_118),
                          xi_142),
                      xi_143),
                  xi_17),
              _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                            omega_shear, omega_shear, omega_shear,
                            omega_shear));
          const __m256 xi_145 = _mm256_mul_ps(
              xi_144,
              _mm256_set_ps(0.12500000000000000f, 0.12500000000000000f,
                            0.12500000000000000f, 0.12500000000000000f,
                            0.12500000000000000f, 0.12500000000000000f,
                            0.12500000000000000f, 0.12500000000000000f));
          const __m256 xi_147 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          xi_112,
                                                          _mm256_set_ps(
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f)),
                                                      _mm256_mul_ps(
                                                          xi_117,
                                                          _mm256_set_ps(
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                                                  _mm256_mul_ps(
                                                      xi_220,
                                                      _mm256_set_ps(
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f))),
                                              _mm256_mul_ps(
                                                  xi_223,
                                                  _mm256_set_ps(-2.0f, -2.0f,
                                                                -2.0f, -2.0f,
                                                                -2.0f, -2.0f,
                                                                -2.0f, -2.0f))),
                                          xi_10),
                                      xi_109),
                                  xi_116),
                              xi_120),
                          xi_143),
                      xi_240),
                  xi_241),
              _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                            omega_shear, omega_shear, omega_shear,
                            omega_shear));
          const __m256 xi_148 = _mm256_mul_ps(
              xi_147,
              _mm256_set_ps(-0.041666666666666667f, -0.041666666666666667f,
                            -0.041666666666666667f, -0.041666666666666667f,
                            -0.041666666666666667f, -0.041666666666666667f,
                            -0.041666666666666667f, -0.041666666666666667f));
          const __m256 xi_149 = _mm256_add_ps(
              _mm256_mul_ps(
                  xi_108, _mm256_set_ps(
                              -0.050000000000000000f, -0.050000000000000000f,
                              -0.050000000000000000f, -0.050000000000000000f,
                              -0.050000000000000000f, -0.050000000000000000f,
                              -0.050000000000000000f, -0.050000000000000000f)),
              xi_148);
          const __m256 xi_150 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          xi_111,
                          _mm256_set_ps(
                              0.014285714285714286f, 0.014285714285714286f,
                              0.014285714285714286f, 0.014285714285714286f,
                              0.014285714285714286f, 0.014285714285714286f,
                              0.014285714285714286f, 0.014285714285714286f)),
                      xi_145),
                  xi_146),
              xi_149);
          const __m256 xi_159 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(
                      xi_111,
                      _mm256_set_ps(
                          -0.035714285714285714f, -0.035714285714285714f,
                          -0.035714285714285714f, -0.035714285714285714f,
                          -0.035714285714285714f, -0.035714285714285714f,
                          -0.035714285714285714f, -0.035714285714285714f)),
                  _mm256_mul_ps(
                      xi_147,
                      _mm256_set_ps(
                          0.083333333333333333f, 0.083333333333333333f,
                          0.083333333333333333f, 0.083333333333333333f,
                          0.083333333333333333f, 0.083333333333333333f,
                          0.083333333333333333f, 0.083333333333333333f))),
              xi_146);
          const __m256 xi_164 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(rho, u_2),
                                  _mm256_mul_ps(vel2Term,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              xi_114),
                          xi_119),
                      xi_161),
                  xi_230),
              xi_6);
          const __m256 xi_165 = _mm256_mul_ps(
              xi_164, _mm256_set_ps(xi_140, xi_140, xi_140, xi_140, xi_140,
                                    xi_140, xi_140, xi_140));
          const __m256 xi_169 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          xi_111,
                          _mm256_set_ps(
                              -0.021428571428571429f, -0.021428571428571429f,
                              -0.021428571428571429f, -0.021428571428571429f,
                              -0.021428571428571429f, -0.021428571428571429f,
                              -0.021428571428571429f, -0.021428571428571429f)),
                      _mm256_mul_ps(
                          xi_125,
                          _mm256_set_ps(
                              0.015873015873015873f, 0.015873015873015873f,
                              0.015873015873015873f, 0.015873015873015873f,
                              0.015873015873015873f, 0.015873015873015873f,
                              0.015873015873015873f, 0.015873015873015873f))),
                  _mm256_mul_ps(xi_145, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
              xi_149);
          const __m256 xi_173 = _mm256_mul_ps(
              xi_144,
              _mm256_set_ps(0.062500000000000000f, 0.062500000000000000f,
                            0.062500000000000000f, 0.062500000000000000f,
                            0.062500000000000000f, 0.062500000000000000f,
                            0.062500000000000000f, 0.062500000000000000f));
          const __m256 xi_175 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_172, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_173),
              xi_174);
          const __m256 xi_181 = _mm256_mul_ps(
              xi_121,
              _mm256_set_ps(0.041666666666666667f, 0.041666666666666667f,
                            0.041666666666666667f, 0.041666666666666667f,
                            0.041666666666666667f, 0.041666666666666667f,
                            0.041666666666666667f, 0.041666666666666667f));
          const __m256 xi_182 = _mm256_add_ps(
              _mm256_mul_ps(
                  xi_147,
                  _mm256_set_ps(0.020833333333333333f, 0.020833333333333333f,
                                0.020833333333333333f, 0.020833333333333333f,
                                0.020833333333333333f, 0.020833333333333333f,
                                0.020833333333333333f, 0.020833333333333333f)),
              xi_181);
          const __m256 xi_183 = _mm256_add_ps(
              _mm256_mul_ps(xi_180, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
              xi_182);
          const __m256 xi_184 = _mm256_add_ps(xi_160, xi_183);
          const __m256 xi_185 =
              _mm256_add_ps(_mm256_add_ps(xi_172, xi_173), xi_174);
          const __m256 xi_186 = _mm256_add_ps(xi_180, xi_182);
          const __m256 xi_187 = _mm256_add_ps(xi_158, xi_186);
          const __m256 xi_191 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_mul_ps(u_2, xi_137), xi_126), xi_18),
                  xi_235),
              _mm256_set_ps(xi_171, xi_171, xi_171, xi_171, xi_171, xi_171,
                            xi_171, xi_171));
          const __m256 xi_193 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(_mm256_add_ps(xi_148, xi_181), xi_190), xi_191),
              xi_192);
          const __m256 xi_199 = _mm256_mul_ps(
              xi_164, _mm256_set_ps(xi_176, xi_176, xi_176, xi_176, xi_176,
                                    xi_176, xi_176, xi_176));
          const __m256 xi_200 = _mm256_add_ps(xi_198, xi_199);
          const __m256 xi_201 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_195, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_197),
              xi_200);
          const __m256 xi_206 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_mul_ps(xi_191,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                          xi_148),
                      xi_181),
                  xi_190),
              xi_192);
          const __m256 xi_207 = _mm256_mul_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_mul_ps(u_2, xi_151), xi_11), xi_154),
                  xi_230),
              _mm256_set_ps(xi_171, xi_171, xi_171, xi_171, xi_171, xi_171,
                            xi_171, xi_171));
          const __m256 xi_208 =
              _mm256_mul_ps(xi_173, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0));
          const __m256 xi_210 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_207, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_208),
              xi_209);
          const __m256 xi_211 = _mm256_add_ps(xi_170, xi_200);
          const __m256 xi_214 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_212, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_183),
              xi_213);
          const __m256 xi_215 =
              _mm256_add_ps(_mm256_add_ps(xi_207, xi_208), xi_209);
          const __m256 xi_216 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_213, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_186),
              xi_212);
          const __m256 xi_217 = _mm256_add_ps(
              _mm256_mul_ps(xi_199, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                  -1.0, -1.0, -1.0)),
              xi_198);
          const __m256 xi_218 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_197, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                  xi_195),
              xi_217);
          const __m256 xi_219 = _mm256_add_ps(xi_168, xi_217);
          const __m256 forceTerm_0 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_mul_ps(
                                  xi_26, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                              _mm256_mul_ps(xi_26,
                                            _mm256_set_ps(xi_29, xi_29, xi_29,
                                                          xi_29, xi_29, xi_29,
                                                          xi_29, xi_29))),
                          _mm256_mul_ps(xi_27,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                      _mm256_mul_ps(xi_27,
                                    _mm256_set_ps(xi_29, xi_29, xi_29, xi_29,
                                                  xi_29, xi_29, xi_29, xi_29))),
                  _mm256_mul_ps(xi_28, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0))),
              _mm256_mul_ps(xi_28, _mm256_set_ps(xi_29, xi_29, xi_29, xi_29,
                                                 xi_29, xi_29, xi_29, xi_29)));
          const __m256 forceTerm_1 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_32, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_30),
              xi_42);
          const __m256 forceTerm_2 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_30, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_32),
              xi_42);
          const __m256 forceTerm_3 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_43, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_45),
              xi_49);
          const __m256 forceTerm_4 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_45, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_43),
              xi_49);
          const __m256 forceTerm_5 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_52, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_50),
              xi_53);
          const __m256 forceTerm_6 = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                     -1.0, -1.0, -1.0, -1.0)),
                  xi_52),
              xi_53);
          const __m256 forceTerm_7 =
              _mm256_add_ps(_mm256_add_ps(xi_58, xi_67), xi_77);
          const __m256 forceTerm_8 =
              _mm256_add_ps(_mm256_add_ps(xi_67, xi_78), xi_79);
          const __m256 forceTerm_9 =
              _mm256_add_ps(_mm256_add_ps(xi_58, xi_79), xi_81);
          const __m256 forceTerm_10 =
              _mm256_add_ps(_mm256_add_ps(xi_77, xi_78), xi_81);
          const __m256 forceTerm_11 =
              _mm256_add_ps(_mm256_add_ps(xi_84, xi_86), xi_91);
          const __m256 forceTerm_12 =
              _mm256_add_ps(_mm256_add_ps(xi_84, xi_92), xi_93);
          const __m256 forceTerm_13 =
              _mm256_add_ps(_mm256_add_ps(xi_102, xi_58), xi_95);
          const __m256 forceTerm_14 =
              _mm256_add_ps(_mm256_add_ps(xi_103, xi_78), xi_95);
          const __m256 forceTerm_15 =
              _mm256_add_ps(_mm256_add_ps(xi_104, xi_86), xi_93);
          const __m256 forceTerm_16 =
              _mm256_add_ps(_mm256_add_ps(xi_104, xi_91), xi_92);
          const __m256 forceTerm_17 =
              _mm256_add_ps(_mm256_add_ps(xi_103, xi_105), xi_58);
          const __m256 forceTerm_18 =
              _mm256_add_ps(_mm256_add_ps(xi_102, xi_105), xi_78);
          _mm256_store_ps(
              &_data_pdfs_20_30_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      xi_108,
                                      _mm256_set_ps(0.10000000000000000f,
                                                    0.10000000000000000f,
                                                    0.10000000000000000f,
                                                    0.10000000000000000f,
                                                    0.10000000000000000f,
                                                    0.10000000000000000f,
                                                    0.10000000000000000f,
                                                    0.10000000000000000f)),
                                  _mm256_mul_ps(
                                      xi_111,
                                      _mm256_set_ps(0.042857142857142857f,
                                                    0.042857142857142857f,
                                                    0.042857142857142857f,
                                                    0.042857142857142857f,
                                                    0.042857142857142857f,
                                                    0.042857142857142857f,
                                                    0.042857142857142857f,
                                                    0.042857142857142857f))),
                              _mm256_mul_ps(
                                  xi_121,
                                  _mm256_set_ps(-0.50000000000000000f,
                                                -0.50000000000000000f,
                                                -0.50000000000000000f,
                                                -0.50000000000000000f,
                                                -0.50000000000000000f,
                                                -0.50000000000000000f,
                                                -0.50000000000000000f,
                                                -0.50000000000000000f))),
                          _mm256_mul_ps(
                              xi_125,
                              _mm256_set_ps(
                                  0.023809523809523810f, 0.023809523809523810f,
                                  0.023809523809523810f, 0.023809523809523810f,
                                  0.023809523809523810f, 0.023809523809523810f,
                                  0.023809523809523810f,
                                  0.023809523809523810f))),
                      forceTerm_0),
                  xi_227));
          _mm256_store_ps(
              &_data_pdfs_20_31_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(xi_129,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0)),
                                  forceTerm_1),
                              xi_136),
                          xi_141),
                      xi_150),
                  xi_241));
          _mm256_store_ps(
              &_data_pdfs_20_32_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(xi_141,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0)),
                                  forceTerm_2),
                              xi_129),
                          xi_135),
                      xi_150),
                  xi_240));
          _mm256_store_ps(
              &_data_pdfs_20_33_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(xi_153,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0)),
                                  forceTerm_3),
                              xi_156),
                          xi_158),
                      xi_159),
                  xi_220));
          _mm256_store_ps(
              &_data_pdfs_20_34_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(xi_156,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0)),
                                  forceTerm_4),
                              xi_153),
                          xi_159),
                      xi_160),
                  xi_223));
          _mm256_store_ps(
              &_data_pdfs_20_35_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(xi_163,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0)),
                                  forceTerm_5),
                              xi_165),
                          xi_168),
                      xi_169),
                  xi_233));
          _mm256_store_ps(
              &_data_pdfs_20_36_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(xi_165,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0)),
                                  forceTerm_6),
                              xi_163),
                          xi_169),
                      xi_170),
                  xi_229));
          _mm256_store_ps(
              &_data_pdfs_20_37_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_7, xi_175), xi_179),
                      xi_184),
                  xi_221));
          _mm256_store_ps(
              &_data_pdfs_20_38_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_8, xi_179), xi_185),
                      xi_187),
                  xi_234));
          _mm256_store_ps(
              &_data_pdfs_20_39_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_9, xi_184), xi_185),
                      xi_189),
                  xi_225));
          _mm256_store_ps(
              &_data_pdfs_20_310_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_10, xi_175),
                                    xi_187),
                      xi_189),
                  xi_224));
          _mm256_store_ps(
              &_data_pdfs_20_311_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_11, xi_193),
                                    xi_201),
                      xi_204),
                  xi_237));
          _mm256_store_ps(
              &_data_pdfs_20_312_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_12, xi_201),
                                    xi_205),
                      xi_206),
                  xi_228));
          _mm256_store_ps(
              &_data_pdfs_20_313_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_13, xi_210),
                                    xi_211),
                      xi_214),
                  xi_232));
          _mm256_store_ps(
              &_data_pdfs_20_314_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_14, xi_211),
                                    xi_215),
                      xi_216),
                  xi_236));
          _mm256_store_ps(
              &_data_pdfs_20_315_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_15, xi_204),
                                    xi_206),
                      xi_218),
                  xi_235));
          _mm256_store_ps(
              &_data_pdfs_20_316_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_16, xi_193),
                                    xi_205),
                      xi_218),
                  xi_238));
          _mm256_store_ps(
              &_data_pdfs_20_317_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_17, xi_214),
                                    xi_215),
                      xi_219),
                  xi_222));
          _mm256_store_ps(
              &_data_pdfs_20_318_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(_mm256_add_ps(forceTerm_18, xi_210),
                                    xi_216),
                      xi_219),
                  xi_230));
        }
        for (int64_t ctr_0 = (int64_t)((_size_force_0) / (8)) * (8);
             ctr_0 < _size_force_0; ctr_0 += 1) {
          const float xi_220 = _data_pdfs_20_33_10[ctr_0];
          const float xi_221 = _data_pdfs_20_37_10[ctr_0];
          const float xi_222 = _data_pdfs_20_317_10[ctr_0];
          const float xi_223 = _data_pdfs_20_34_10[ctr_0];
          const float xi_224 = _data_pdfs_20_310_10[ctr_0];
          const float xi_225 = _data_pdfs_20_39_10[ctr_0];
          const float xi_226 = _data_force_20_30_10[ctr_0];
          const float xi_227 = _data_pdfs_20_30_10[ctr_0];
          const float xi_228 = _data_pdfs_20_312_10[ctr_0];
          const float xi_229 = _data_pdfs_20_36_10[ctr_0];
          const float xi_230 = _data_pdfs_20_318_10[ctr_0];
          const float xi_231 = _data_force_20_31_10[ctr_0];
          const float xi_232 = _data_pdfs_20_313_10[ctr_0];
          const float xi_233 = _data_pdfs_20_35_10[ctr_0];
          const float xi_234 = _data_pdfs_20_38_10[ctr_0];
          const float xi_235 = _data_pdfs_20_315_10[ctr_0];
          const float xi_236 = _data_pdfs_20_314_10[ctr_0];
          const float xi_237 = _data_pdfs_20_311_10[ctr_0];
          const float xi_238 = _data_pdfs_20_316_10[ctr_0];
          const float xi_239 = _data_force_20_32_10[ctr_0];
          const float xi_240 = _data_pdfs_20_32_10[ctr_0];
          const float xi_241 = _data_pdfs_20_31_10[ctr_0];
          const float xi_0 = xi_230 + xi_236;
          const float xi_1 = xi_0 + xi_223;
          const float xi_2 = xi_235 + xi_237 + xi_241;
          const float xi_3 = xi_228 + xi_233;
          const float xi_4 = xi_220 + xi_225;
          const float xi_5 = xi_238 + xi_240;
          const float xi_6 = xi_222 + xi_229;
          const float xi_9 = -xi_225;
          const float xi_10 = -xi_221 + xi_9;
          const float xi_11 = -xi_222;
          const float xi_12 = -xi_232;
          const float xi_13 = -xi_220;
          const float xi_14 = xi_11 + xi_12 + xi_13;
          const float xi_15 = -xi_240;
          const float xi_16 = -xi_224;
          const float xi_17 = xi_15 + xi_16;
          const float xi_18 = -xi_238;
          const float xi_19 = -xi_228;
          const float xi_20 = xi_18 + xi_19;
          const float xi_21 = -xi_230;
          const float xi_22 = xi_11 + xi_21;
          const float xi_23 = -xi_235;
          const float xi_24 = -xi_229;
          const float xi_25 = xi_18 + xi_23 + xi_237 + xi_24;
          const float xi_30 = xi_231 * 0.16666666666666667f;
          const float xi_31 = xi_231 * 0.083333333333333333f;
          const float xi_43 = xi_226 * 0.16666666666666667f;
          const float xi_44 = xi_226 * 0.083333333333333333f;
          const float xi_50 = xi_239 * 0.16666666666666667f;
          const float xi_51 = xi_239 * 0.083333333333333333f;
          const float xi_68 = xi_231 * 0.25000000000000000f;
          const float xi_73 = xi_231 * xi_72;
          const float xi_106 = -xi_227;
          const float xi_107 = xi_106 + xi_229 * 3.0f + xi_233 * 3.0f;
          const float xi_108 =
              omega_even *
              (xi_107 + xi_228 * -3.0f + xi_235 * -3.0f + xi_237 * -3.0f +
               xi_238 * -3.0f + xi_240 * 3.0f + xi_241 * 3.0f);
          const float xi_109 =
              xi_228 * 2.0f + xi_235 * 2.0f + xi_237 * 2.0f + xi_238 * 2.0f;
          const float xi_110 = xi_109 + xi_220 * 5.0f + xi_223 * 5.0f;
          const float xi_111 =
              omega_even * (xi_107 + xi_110 + xi_222 * -5.0f + xi_230 * -5.0f +
                            xi_232 * -5.0f + xi_236 * -5.0f + xi_240 * -2.0f +
                            xi_241 * -2.0f);
          const float xi_114 = -xi_237;
          const float xi_115 = xi_114 + xi_19;
          const float xi_116 = -xi_234;
          const float xi_119 = -xi_236;
          const float xi_120 = xi_119 + xi_12 + xi_16 + xi_22;
          const float xi_122 = xi_232 * 2.0f;
          const float xi_123 = xi_236 * 2.0f;
          const float xi_124 = xi_222 * 2.0f + xi_230 * 2.0f;
          const float xi_125 =
              omega_even *
              (xi_106 + xi_110 + xi_122 + xi_123 + xi_124 + xi_221 * -7.0f +
               xi_224 * -7.0f + xi_225 * -7.0f + xi_229 * -4.0f +
               xi_233 * -4.0f + xi_234 * -7.0f + xi_240 * 5.0f + xi_241 * 5.0f);
          const float xi_126 = xi_114 + xi_228;
          const float xi_127 = xi_126 + xi_15 + xi_23 + xi_238 + xi_241;
          const float xi_129 = xi_127 * xi_128;
          const float xi_130 = xi_221 * 2.0f;
          const float xi_131 = xi_224 * 2.0f;
          const float xi_132 = xi_225 * 2.0f + xi_234 * -2.0f;
          const float xi_133 = -xi_130 + xi_131 + xi_132 + xi_15 + xi_2 + xi_20;
          const float xi_135 = xi_133 * xi_134;
          const float xi_136 = -xi_135;
          const float xi_138 = xi_116 + xi_224;
          const float xi_142 = xi_222 + xi_232;
          const float xi_146 = xi_125 * -0.019841269841269841f;
          const float xi_154 = xi_119 + xi_232;
          const float xi_155 = xi_13 + xi_154 + xi_21 + xi_222 + xi_223;
          const float xi_156 = xi_128 * xi_155;
          const float xi_157 = xi_1 + xi_130 - xi_131 + xi_132 + xi_14;
          const float xi_158 = xi_134 * xi_157;
          const float xi_160 = -xi_158;
          const float xi_161 = xi_235 + xi_238;
          const float xi_162 = xi_115 + xi_161 + xi_233 + xi_24;
          const float xi_163 = xi_128 * xi_162;
          const float xi_166 = -xi_122 - xi_123 + xi_124 + xi_25 + xi_3;
          const float xi_167 = xi_134 * xi_166;
          const float xi_168 = -xi_167;
          const float xi_170 = xi_167;
          const float xi_174 = xi_125 * 0.013888888888888889f;
          const float xi_190 = xi_111 * -0.0071428571428571429f;
          const float xi_192 = xi_108 * 0.025000000000000000f;
          const float xi_195 = xi_166 * xi_194;
          const float xi_197 = xi_162 * xi_196;
          const float xi_198 = xi_125 * -0.0039682539682539683f;
          const float xi_202 = xi_133 * xi_194;
          const float xi_203 = xi_127 * xi_196;
          const float xi_209 = xi_111 * 0.017857142857142857f;
          const float xi_212 = xi_155 * xi_196;
          const float xi_213 = xi_157 * xi_194;
          const float xi_32 = rr_0 * xi_31;
          const float xi_45 = rr_0 * xi_44;
          const float xi_52 = rr_0 * xi_51;
          const float xi_55 = xi_226 * xi_54;
          const float xi_60 = xi_231 * xi_54;
          const float xi_82 = xi_239 * xi_54;
          const float vel0Term = xi_1 + xi_224 + xi_234;
          const float vel1Term = xi_2 + xi_221;
          const float vel2Term = xi_232 + xi_3;
          const float rho =
              vel0Term + vel1Term + vel2Term + xi_227 + xi_4 + xi_5 + xi_6;
          const float xi_7 = 1 / (rho);
          const float xi_8 = xi_7 * 0.50000000000000000f;
          const float u_0 = xi_226 * xi_8 + xi_7 * (vel0Term + xi_10 + xi_14);
          const float xi_26 = u_0 * xi_226;
          const float xi_38 = xi_26 * 0.16666666666666667f;
          const float xi_39 = xi_26 * 0.083333333333333333f;
          const float xi_40 = omega_shear * xi_39;
          const float xi_41 = -xi_38 + xi_40;
          const float xi_57 = -xi_26 * xi_56 + xi_38;
          const float xi_58 = -xi_44 + xi_55 + xi_57;
          const float xi_62 = -xi_26 * xi_61;
          const float xi_69 = u_0 * xi_68;
          const float xi_74 = u_0 * xi_73;
          const float xi_78 = xi_44 - xi_55 + xi_57;
          const float xi_85 = -xi_39;
          const float xi_96 = u_0 * xi_239;
          const float xi_97 = xi_96 * 0.25000000000000000f;
          const float xi_100 = xi_72 * xi_96;
          const float xi_112 = rho * (u_0 * u_0);
          const float xi_151 = rho * u_0;
          const float xi_152 = -vel0Term + xi_142 + xi_151 + xi_221 + xi_4;
          const float xi_153 = xi_140 * xi_152;
          const float xi_180 = xi_152 * xi_176;
          const float u_1 =
              xi_231 * xi_8 + xi_7 * (vel1Term + xi_17 + xi_20 + xi_234 + xi_9);
          const float xi_27 = u_1 * xi_231;
          const float xi_33 = xi_27 * 0.16666666666666667f;
          const float xi_46 = xi_27 * 0.083333333333333333f;
          const float xi_47 = omega_shear * xi_46;
          const float xi_48 = -xi_33 + xi_47;
          const float xi_63 = -xi_27 * xi_61;
          const float xi_70 = u_1 * 0.25000000000000000f;
          const float xi_71 = xi_226 * xi_70;
          const float xi_75 = u_1 * xi_72;
          const float xi_76 = xi_226 * xi_75;
          const float xi_77 = -xi_69 - xi_71 + xi_74 + xi_76;
          const float xi_79 = xi_69 + xi_71 - xi_74 - xi_76;
          const float xi_87 = xi_239 * xi_70;
          const float xi_89 = xi_239 * xi_75;
          const float xi_94 = -xi_46;
          const float xi_117 = rho * (u_1 * u_1);
          const float xi_118 = xi_10 + xi_116 + xi_117;
          const float xi_137 = rho * u_1;
          const float xi_139 =
              -vel1Term + xi_137 + xi_138 + xi_225 + xi_228 + xi_5;
          const float xi_141 = xi_139 * xi_140;
          const float xi_172 = xi_171 * (u_0 * xi_137 + xi_138 + xi_221 + xi_9);
          const float xi_177 = xi_139 * xi_176;
          const float xi_178 = xi_177;
          const float xi_179 = xi_135 + xi_178;
          const float xi_188 = -xi_177;
          const float xi_189 = xi_136 + xi_188;
          const float xi_204 = xi_178 - xi_202 + xi_203;
          const float xi_205 = xi_188 + xi_202 - xi_203;
          const float u_2 =
              xi_239 * xi_8 + xi_7 * (vel2Term + xi_22 + xi_236 + xi_25);
          const float xi_28 = u_2 * xi_239;
          const float xi_34 = xi_28 * 0.16666666666666667f;
          const float xi_35 = xi_28 * 0.083333333333333333f;
          const float xi_36 = omega_shear * xi_35;
          const float xi_37 = -xi_34 + xi_36;
          const float xi_42 = -omega_shear * xi_33 +
                              xi_27 * 0.33333333333333333f + xi_37 + xi_41;
          const float xi_49 = -omega_shear * xi_38 +
                              xi_26 * 0.33333333333333333f + xi_37 + xi_48;
          const float xi_53 = -omega_shear * xi_34 +
                              xi_28 * 0.33333333333333333f + xi_41 + xi_48;
          const float xi_59 = -xi_35;
          const float xi_64 = -xi_28 * xi_61;
          const float xi_65 = -xi_27 * xi_56 + xi_33 + xi_62 + xi_63 + xi_64;
          const float xi_66 = xi_31 - xi_60 + xi_65;
          const float xi_67 = xi_36 + xi_59 + xi_66;
          const float xi_80 = -xi_31 + xi_60 + xi_65;
          const float xi_81 = xi_36 + xi_59 + xi_80;
          const float xi_83 = -xi_28 * xi_56 + xi_34;
          const float xi_84 = xi_51 - xi_82 + xi_83;
          const float xi_86 = xi_40 + xi_66 + xi_85;
          const float xi_88 = u_2 * xi_68;
          const float xi_90 = u_2 * xi_73;
          const float xi_91 = xi_87 + xi_88 - xi_89 - xi_90;
          const float xi_92 = xi_40 + xi_80 + xi_85;
          const float xi_93 = -xi_87 - xi_88 + xi_89 + xi_90;
          const float xi_95 = xi_47 + xi_62 + xi_63 + xi_64 + xi_84 + xi_94;
          const float xi_98 = u_2 * xi_226;
          const float xi_99 = xi_98 * 0.25000000000000000f;
          const float xi_101 = xi_72 * xi_98;
          const float xi_102 = xi_100 + xi_101 - xi_97 - xi_99;
          const float xi_103 = -xi_100 - xi_101 + xi_97 + xi_99;
          const float xi_104 = -xi_51 + xi_82 + xi_83;
          const float xi_105 = xi_104 + xi_47 + xi_62 + xi_63 + xi_64 + xi_94;
          const float xi_113 = rho * (u_2 * u_2);
          const float xi_121 = omega_bulk * (xi_112 + xi_113 + xi_115 + xi_118 +
                                             xi_120 + xi_18 + xi_227 + xi_23);
          const float xi_143 = -xi_113 + xi_229 + xi_233;
          const float xi_144 =
              omega_shear * (xi_0 + xi_118 + xi_142 + xi_143 + xi_17 - xi_241);
          const float xi_145 = xi_144 * 0.12500000000000000f;
          const float xi_147 =
              omega_shear *
              (xi_10 + xi_109 + xi_112 * 2.0f + xi_116 - xi_117 + xi_120 +
               xi_143 + xi_220 * -2.0f + xi_223 * -2.0f + xi_240 + xi_241);
          const float xi_148 = xi_147 * -0.041666666666666667f;
          const float xi_149 = xi_108 * -0.050000000000000000f + xi_148;
          const float xi_150 =
              xi_111 * 0.014285714285714286f + xi_145 + xi_146 + xi_149;
          const float xi_159 = xi_111 * -0.035714285714285714f + xi_146 +
                               xi_147 * 0.083333333333333333f;
          const float xi_164 =
              rho * u_2 - vel2Term + xi_114 + xi_119 + xi_161 + xi_230 + xi_6;
          const float xi_165 = xi_140 * xi_164;
          const float xi_169 = xi_111 * -0.021428571428571429f +
                               xi_125 * 0.015873015873015873f - xi_145 + xi_149;
          const float xi_173 = xi_144 * 0.062500000000000000f;
          const float xi_175 = -xi_172 + xi_173 + xi_174;
          const float xi_181 = xi_121 * 0.041666666666666667f;
          const float xi_182 = xi_147 * 0.020833333333333333f + xi_181;
          const float xi_183 = -xi_180 + xi_182;
          const float xi_184 = xi_160 + xi_183;
          const float xi_185 = xi_172 + xi_173 + xi_174;
          const float xi_186 = xi_180 + xi_182;
          const float xi_187 = xi_158 + xi_186;
          const float xi_191 =
              xi_171 * (u_2 * xi_137 + xi_126 + xi_18 + xi_235);
          const float xi_193 = xi_148 + xi_181 + xi_190 + xi_191 + xi_192;
          const float xi_199 = xi_164 * xi_176;
          const float xi_200 = xi_198 + xi_199;
          const float xi_201 = -xi_195 + xi_197 + xi_200;
          const float xi_206 = xi_148 + xi_181 + xi_190 - xi_191 + xi_192;
          const float xi_207 =
              xi_171 * (u_2 * xi_151 + xi_11 + xi_154 + xi_230);
          const float xi_208 = -xi_173;
          const float xi_210 = -xi_207 + xi_208 + xi_209;
          const float xi_211 = xi_170 + xi_200;
          const float xi_214 = xi_183 - xi_212 + xi_213;
          const float xi_215 = xi_207 + xi_208 + xi_209;
          const float xi_216 = xi_186 + xi_212 - xi_213;
          const float xi_217 = xi_198 - xi_199;
          const float xi_218 = xi_195 - xi_197 + xi_217;
          const float xi_219 = xi_168 + xi_217;
          const float forceTerm_0 = xi_26 * xi_29 - xi_26 + xi_27 * xi_29 -
                                    xi_27 + xi_28 * xi_29 - xi_28;
          const float forceTerm_1 = xi_30 - xi_32 + xi_42;
          const float forceTerm_2 = -xi_30 + xi_32 + xi_42;
          const float forceTerm_3 = -xi_43 + xi_45 + xi_49;
          const float forceTerm_4 = xi_43 - xi_45 + xi_49;
          const float forceTerm_5 = xi_50 - xi_52 + xi_53;
          const float forceTerm_6 = -xi_50 + xi_52 + xi_53;
          const float forceTerm_7 = xi_58 + xi_67 + xi_77;
          const float forceTerm_8 = xi_67 + xi_78 + xi_79;
          const float forceTerm_9 = xi_58 + xi_79 + xi_81;
          const float forceTerm_10 = xi_77 + xi_78 + xi_81;
          const float forceTerm_11 = xi_84 + xi_86 + xi_91;
          const float forceTerm_12 = xi_84 + xi_92 + xi_93;
          const float forceTerm_13 = xi_102 + xi_58 + xi_95;
          const float forceTerm_14 = xi_103 + xi_78 + xi_95;
          const float forceTerm_15 = xi_104 + xi_86 + xi_93;
          const float forceTerm_16 = xi_104 + xi_91 + xi_92;
          const float forceTerm_17 = xi_103 + xi_105 + xi_58;
          const float forceTerm_18 = xi_102 + xi_105 + xi_78;
          _data_pdfs_20_30_10[ctr_0] =
              forceTerm_0 + xi_108 * 0.10000000000000000f +
              xi_111 * 0.042857142857142857f + xi_121 * -0.50000000000000000f +
              xi_125 * 0.023809523809523810f + xi_227;
          _data_pdfs_20_31_10[ctr_0] =
              forceTerm_1 - xi_129 + xi_136 + xi_141 + xi_150 + xi_241;
          _data_pdfs_20_32_10[ctr_0] =
              forceTerm_2 + xi_129 + xi_135 - xi_141 + xi_150 + xi_240;
          _data_pdfs_20_33_10[ctr_0] =
              forceTerm_3 - xi_153 + xi_156 + xi_158 + xi_159 + xi_220;
          _data_pdfs_20_34_10[ctr_0] =
              forceTerm_4 + xi_153 - xi_156 + xi_159 + xi_160 + xi_223;
          _data_pdfs_20_35_10[ctr_0] =
              forceTerm_5 - xi_163 + xi_165 + xi_168 + xi_169 + xi_233;
          _data_pdfs_20_36_10[ctr_0] =
              forceTerm_6 + xi_163 - xi_165 + xi_169 + xi_170 + xi_229;
          _data_pdfs_20_37_10[ctr_0] =
              forceTerm_7 + xi_175 + xi_179 + xi_184 + xi_221;
          _data_pdfs_20_38_10[ctr_0] =
              forceTerm_8 + xi_179 + xi_185 + xi_187 + xi_234;
          _data_pdfs_20_39_10[ctr_0] =
              forceTerm_9 + xi_184 + xi_185 + xi_189 + xi_225;
          _data_pdfs_20_310_10[ctr_0] =
              forceTerm_10 + xi_175 + xi_187 + xi_189 + xi_224;
          _data_pdfs_20_311_10[ctr_0] =
              forceTerm_11 + xi_193 + xi_201 + xi_204 + xi_237;
          _data_pdfs_20_312_10[ctr_0] =
              forceTerm_12 + xi_201 + xi_205 + xi_206 + xi_228;
          _data_pdfs_20_313_10[ctr_0] =
              forceTerm_13 + xi_210 + xi_211 + xi_214 + xi_232;
          _data_pdfs_20_314_10[ctr_0] =
              forceTerm_14 + xi_211 + xi_215 + xi_216 + xi_236;
          _data_pdfs_20_315_10[ctr_0] =
              forceTerm_15 + xi_204 + xi_206 + xi_218 + xi_235;
          _data_pdfs_20_316_10[ctr_0] =
              forceTerm_16 + xi_193 + xi_205 + xi_218 + xi_238;
          _data_pdfs_20_317_10[ctr_0] =
              forceTerm_17 + xi_214 + xi_215 + xi_219 + xi_222;
          _data_pdfs_20_318_10[ctr_0] =
              forceTerm_18 + xi_210 + xi_216 + xi_219 + xi_230;
        }
      }
    }
  }
}
} // namespace internal_dfe735f2f0357dcc08d993f108917b07

void CollideSweepSinglePrecisionAVX::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &omega_shear = this->omega_shear_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_dfe735f2f0357dcc08d993f108917b07::
      collidesweepsingleprecisionavx_collidesweepsingleprecisionavx(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, omega_bulk, omega_even, omega_odd,
          omega_shear);
}

void CollideSweepSinglePrecisionAVX::runOnCellInterval(
    const shared_ptr<StructuredBlockStorage> &blocks,
    const CellInterval &globalCellInterval, cell_idx_t ghostLayers,
    IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &omega_shear = this->omega_shear_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_dfe735f2f0357dcc08d993f108917b07::
      collidesweepsingleprecisionavx_collidesweepsingleprecisionavx(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, omega_bulk, omega_even, omega_odd,
          omega_shear);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif