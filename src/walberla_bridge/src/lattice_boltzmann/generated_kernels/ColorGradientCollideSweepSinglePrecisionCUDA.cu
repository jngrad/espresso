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
//! \\file ColorGradientCollideSweepSinglePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepSinglePrecisionCUDA.h"


#include "philox_rand.h"



#define FUNC_PREFIX __global__

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning push
#pragma warning( disable :  1599 )
#endif

using namespace std;

namespace walberla {
namespace pystencils {


namespace internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda(float * RESTRICT const _data_force_a, float * RESTRICT const _data_force_b, float * RESTRICT  _data_pdfs_a, float * RESTRICT  _data_pdfs_b, float * RESTRICT const _data_phasefield, float * RESTRICT const _data_rho_a, float * RESTRICT const _data_rho_b, float * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_0, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_0, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_0, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_0, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_0, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_0, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_0, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, float beta, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, float kT, float omega_shear_a, float omega_shear_b, uint32_t seed, float sigma, uint32_t time_step)
{
   if (blockDim.x*blockIdx.x + threadIdx.x + 1 < _size_force_a_0 - 1 && blockDim.y*blockIdx.y + threadIdx.y + 1 < _size_force_a_1 - 1 && blockDim.z*blockIdx.z + threadIdx.z + 1 < _size_force_a_2 - 1)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x + 1;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y + 1;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z + 1;
      
      float random_3_0{};
      float random_3_1{};
      float random_3_2{};
      float random_3_3{};
      if (kT > 0.) {  
      philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);
      }
      
      float random_2_0{};
      float random_2_1{};
      float random_2_2{};
      float random_2_3{};
      if (kT > 0.) {  
      philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);
      }
      
      float random_1_0{};
      float random_1_1{};
      float random_1_2{};
      float random_1_3{};
      if (kT > 0.) {  
      philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);
      }
      
      float random_0_0{};
      float random_0_1{};
      float random_0_2{};
      float random_0_3{};
      if (kT > 0.) {  
      philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);
      }
      
      float random_7_0_shifted{};
      float random_7_1_shifted{};
      float random_7_2_shifted{};
      float random_7_3_shifted{};
      if (kT > 0.) {  
      philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 19, seed, random_7_0_shifted, random_7_1_shifted, random_7_2_shifted, random_7_3_shifted);
      }
      
      float random_6_0_shifted{};
      float random_6_1_shifted{};
      float random_6_2_shifted{};
      float random_6_3_shifted{};
      if (kT > 0.) {  
      philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 18, seed, random_6_0_shifted, random_6_1_shifted, random_6_2_shifted, random_6_3_shifted);
      }
      
      float random_5_0_shifted{};
      float random_5_1_shifted{};
      float random_5_2_shifted{};
      float random_5_3_shifted{};
      if (kT > 0.) {  
      philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 17, seed, random_5_0_shifted, random_5_1_shifted, random_5_2_shifted, random_5_3_shifted);
      }
      
      float random_4_0_shifted{};
      float random_4_1_shifted{};
      float random_4_2_shifted{};
      float random_4_3_shifted{};
      if (kT > 0.) {  
      philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 16, seed, random_4_0_shifted, random_4_1_shifted, random_4_2_shifted, random_4_3_shifted);
      }
      const float xi_2 = ((1.0f) / (omega_shear_a + omega_shear_b));
      const float xi_3 = omega_shear_a*omega_shear_b;
      const float xi_4 = xi_2*xi_3;
      const float xi_5 = xi_4*2.0f;
      const float xi_6 = omega_shear_a*-4.0f;
      const float xi_7 = xi_4*8.0f;
      const float xi_8 = xi_6 + xi_7;
      const float xi_9 = -xi_8;
      const float xi_12 = omega_shear_b*-4.0f;
      const float xi_13 = xi_12 + xi_7;
      const float xi_212 = 3.7416573867739413f;
      const float xi_213 = xi_212*0.14285714285714285f;
      const float xi_216 = random_3_0 - 0.5f;
      const float xi_218 = 5.4772255750516612f;
      const float xi_219 = xi_218*0.20000000000000001f;
      const float xi_220 = random_3_2 - 0.5f;
      const float xi_222 = random_0_0 - 0.5f;
      const float xi_223 = 2.4494897427831779f;
      const float xi_225 = 8.3666002653407556f;
      const float xi_226 = xi_225*0.085714285714285715f;
      const float xi_227 = random_3_1 - 0.5f;
      const float xi_229 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_230 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_231 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_232 = -0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_233 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_234 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_235 = -xi_231 + xi_232 + xi_233 + xi_234;
      const float xi_236 = -0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_237 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_238 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_239 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_240 = xi_236 - xi_237 + xi_238 + xi_239;
      const float xi_241 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_242 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_243 = -xi_241 + xi_242;
      const float xi_244 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_245 = -xi_244;
      const float xi_246 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_247 = xi_245 + xi_246;
      const float xi_248 = -0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_249 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_250 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_251 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_252 = -xi_250 + xi_251;
      const float xi_253 = xi_248 + xi_249 + xi_252;
      const float xi_254 = xi_247 + xi_253;
      const float xi_255 = -xi_229 + xi_230 + xi_235 + xi_240 + xi_243 + xi_254;
      const float xi_256 = xi_255*xi_255;
      const float xi_257 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_258 = -xi_257;
      const float xi_259 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_260 = xi_245 + xi_246 + xi_258 + xi_259;
      const float xi_261 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_262 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_263 = -xi_261 + xi_262;
      const float xi_264 = xi_229 - xi_230;
      const float xi_265 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_266 = -xi_265;
      const float xi_267 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_268 = xi_250 - xi_251 + xi_266 + xi_267;
      const float xi_269 = xi_235 + xi_237 - xi_239 + xi_260 + xi_263 + xi_264 + xi_268;
      const float xi_270 = xi_269*xi_269;
      const float xi_271 = xi_261 - xi_262;
      const float xi_272 = xi_258 + xi_259;
      const float xi_273 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const float xi_274 = -xi_273;
      const float xi_275 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const float xi_276 = xi_241 - xi_242 + xi_274 + xi_275;
      const float xi_277 = xi_231 + xi_232 + xi_233 - xi_234 + xi_240 + xi_264 + xi_271 + xi_272 + xi_276;
      const float xi_278 = xi_277*xi_277;
      const float xi_279 = xi_256 + xi_270 + xi_278;
      const float xi_280 = powf(xi_279, 0.5f);
      const float xi_281 = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2]*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const float xi_282 = sigma*xi_280*((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] > 0.5f) ? (omega_shear_a): ((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] < -0.5f) ? (omega_shear_b): ((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] > 0.0f) ? (xi_281*xi_8 + xi_5 + xi_9*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2]): (xi_13*xi_281 + xi_13*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] + xi_5))));
      const bool xi_283 = xi_280 > 1.0000000000000001e-5f;
      const float xi_284 = ((xi_283) ? (xi_282*0.25f): (0.0f));
      const float xi_286 = ((1.0f) / (xi_279));
      const float xi_287 = xi_282*1.125f;
      const float xi_288 = -((xi_283) ? (xi_287*(xi_270*xi_286*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_293 = random_0_2 - 0.5f;
      const float xi_308 = random_2_3 - 0.5f;
      const float xi_310 = 1.4142135623730951f;
      const float xi_311 = xi_310*0.5f;
      const float xi_318 = random_1_2 - 0.5f;
      const float xi_332 = 1.7320508075688772f;
      const float xi_333 = xi_332*0.16666666666666666f;
      const float xi_334 = random_0_1 - 0.5f;
      const float xi_340 = ((xi_283) ? (xi_287*(xi_286*0.055555555555555552f*(xi_269*xi_269) - 0.018518518518518517f)): (0.0f));
      const float xi_348 = ((xi_283) ? (xi_287*(xi_286*0.055555555555555552f*(xi_277*xi_277) - 0.018518518518518517f)): (0.0f));
      const float xi_352 = random_2_1 - 0.5f;
      const float xi_356 = random_2_0 - 0.5f;
      const float xi_360 = -((xi_283) ? (xi_287*(xi_278*xi_286*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_361 = -((xi_283) ? (xi_287*(xi_286*0.055555555555555552f*(xi_255*xi_255) - 0.018518518518518517f)): (0.0f));
      const float xi_367 = random_2_2 - 0.5f;
      const float xi_370 = random_1_3 - 0.5f;
      const float xi_374 = ((xi_283) ? (xi_287*(xi_256*xi_286*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_375 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_376 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_377 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_378 = xi_236 + xi_238;
      const float xi_379 = xi_276 + xi_378;
      const float xi_380 = -0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_381 = xi_265 - xi_267 + xi_380;
      const float xi_382 = xi_244 - xi_246 + xi_252 - xi_375 + xi_376 + xi_377 + xi_379 + xi_381;
      const float xi_383 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_382*xi_382) - 0.037037037037037035f)): (0.0f));
      const float xi_397 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_398 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_399 = -0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_400 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_401 = xi_247 + xi_268 + xi_379 - xi_397 + xi_398 + xi_399 + xi_400;
      const float xi_402 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_401*xi_401) - 0.037037037037037035f)): (0.0f));
      const float xi_404 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_401*xi_401) - 0.037037037037037035f)): (0.0f));
      const float xi_406 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_382*xi_382) - 0.037037037037037035f)): (0.0f));
      const float xi_407 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_408 = xi_243 + xi_248 + xi_249 + xi_378;
      const float xi_409 = xi_257 - xi_259 + xi_397 - xi_398;
      const float xi_410 = xi_271 + xi_381 + xi_407 + xi_408 + xi_409;
      const float xi_411 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_410*xi_410) - 0.037037037037037035f)): (0.0f));
      const float xi_420 = xi_310*0.25f;
      const float xi_435 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_436 = xi_263 + xi_375 - xi_376;
      const float xi_437 = xi_266 + xi_267 + xi_272 + xi_399 + xi_408 + xi_435 + xi_436;
      const float xi_438 = -((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_437*xi_437) - 0.037037037037037035f)): (0.0f));
      const float xi_441 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_442 = xi_253 + xi_260 + xi_271 + xi_274 + xi_275 + xi_380 + xi_399 + xi_441;
      const float xi_443 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_442*xi_442) - 0.037037037037037035f)): (0.0f));
      const float xi_457 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_458 = xi_254 + xi_273 - xi_275 + xi_409 + xi_436 + xi_457;
      const float xi_459 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_458*xi_458) - 0.037037037037037035f)): (0.0f));
      const float xi_462 = -((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_437*xi_437) - 0.037037037037037035f)): (0.0f));
      const float xi_464 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_410*xi_410) - 0.037037037037037035f)): (0.0f));
      const float xi_465 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_458*xi_458) - 0.037037037037037035f)): (0.0f));
      const float xi_467 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776f*(xi_442*xi_442) - 0.037037037037037035f)): (0.0f));
      const float xi_488 = random_7_0_shifted - 0.5f;
      const float xi_490 = random_7_2_shifted - 0.5f;
      const float xi_492 = random_4_0_shifted - 0.5f;
      const float xi_494 = random_7_1_shifted - 0.5f;
      const float xi_501 = random_4_2_shifted - 0.5f;
      const float xi_515 = random_6_3_shifted - 0.5f;
      const float xi_524 = random_5_2_shifted - 0.5f;
      const float xi_538 = random_4_1_shifted - 0.5f;
      const float xi_554 = random_6_1_shifted - 0.5f;
      const float xi_558 = random_6_0_shifted - 0.5f;
      const float xi_567 = random_6_2_shifted - 0.5f;
      const float xi_570 = random_5_3_shifted - 0.5f;
      const float xia_3_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3];
      const float xi_303 = -xia_3_collide;
      const float xia_4_collide = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const bool xi_0 = xia_4_collide > 0.5f;
      const bool xi_1 = xia_4_collide < -0.5f;
      const float xi_10 = xia_4_collide*xia_4_collide;
      const bool xi_11 = xia_4_collide > 0.0f;
      const float xi_14 = ((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_8 + xi_5 + xi_9*xia_4_collide): (xi_10*xi_13 + xi_13*xia_4_collide + xi_5))));
      const float xi_15 = ((1.0f) / (xi_14*-0.25f + 2.0f));
      const float xia_5_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2];
      const float xi_43 = xia_5_collide*0.083333333333333329f;
      const float xia_6_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const float xia_7_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3];
      const float xi_302 = -xia_7_collide;
      const float xia_8_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3];
      const float xia_9_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3];
      const float xi_285 = -xia_9_collide;
      const float xia_10_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3];
      const float xia_11_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3];
      const float xi_200 = xia_11_collide + xia_7_collide;
      const float xia_12_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3];
      const float xia_13_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const float xia_14_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3];
      const float xi_353 = -xia_14_collide;
      const float xia_15_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3];
      const float xi_207 = xia_15_collide*2.0f;
      const float xia_16_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3];
      const float xi_304 = xi_303 + xia_16_collide;
      const float xi_305 = xi_302 + xi_304 + xia_11_collide;
      const float xi_364 = -xia_16_collide;
      const float xi_365 = xi_200 + xi_303 + xi_364;
      const float xia_17_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3];
      const float xia_18_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3];
      const float xi_296 = -xia_17_collide - xia_18_collide;
      const float xi_313 = xia_18_collide*2.0f;
      const float xia_19_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3];
      const float xi_28 = xia_19_collide*0.083333333333333329f;
      const float xi_62 = xia_19_collide*0.25f;
      const float xia_20_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3];
      const float xi_52 = xia_20_collide*0.083333333333333329f;
      const float xia_21_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const float xia_22_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3];
      const float xi_203 = xia_10_collide + xia_22_collide;
      const float xi_204 = xi_203 + xia_17_collide + xia_18_collide;
      const float xi_314 = xia_10_collide*2.0f + xia_22_collide*-2.0f;
      const float xia_23_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3];
      const float xi_327 = -xia_23_collide;
      const float xia_24_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2];
      const float xi_194 = -xia_24_collide;
      const float xi_195 = xi_194 + xia_11_collide*2.0f + xia_14_collide*5.0f + xia_16_collide*2.0f + xia_3_collide*2.0f + xia_7_collide*2.0f + xia_8_collide*5.0f;
      const float xia_25_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3];
      const float xi_208 = xia_25_collide*2.0f;
      const float xi_349 = xia_15_collide - xia_25_collide;
      const float xia_26_collide = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const float xia_27_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3];
      const float xi_297 = -xia_27_collide;
      const float xi_306 = -xi_285 - xi_305 - xia_27_collide;
      const float xi_315 = xi_297 + xi_305 - xi_313 + xi_314 + xia_17_collide*2.0f + xia_9_collide;
      const float xia_28_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3];
      const float xi_192 = xia_11_collide*-3.0f + xia_16_collide*-3.0f + xia_23_collide*3.0f - xia_24_collide + xia_27_collide*3.0f + xia_28_collide*3.0f + xia_3_collide*-3.0f + xia_7_collide*-3.0f + xia_9_collide*3.0f;
      const float xi_328 = -xia_28_collide;
      const float xi_366 = -xi_327 - xi_365 - xia_28_collide;
      const float xia_29_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3];
      const float xi_196 = xi_195 + xia_12_collide*-5.0f + xia_15_collide*-5.0f + xia_23_collide*3.0f + xia_25_collide*-5.0f + xia_27_collide*-2.0f + xia_28_collide*3.0f + xia_29_collide*-5.0f + xia_9_collide*-2.0f;
      const float xi_209 = xia_12_collide*2.0f + xia_29_collide*2.0f;
      const float xi_210 = xi_195 + xi_207 + xi_208 + xi_209 + xia_10_collide*-7.0f + xia_17_collide*-7.0f + xia_18_collide*-7.0f + xia_22_collide*-7.0f + xia_23_collide*-4.0f + xia_27_collide*5.0f + xia_28_collide*-4.0f + xia_9_collide*5.0f;
      const float xi_350 = xi_349 + xia_12_collide - xia_29_collide;
      const float xi_351 = xi_350 + xia_14_collide - xia_8_collide;
      const float xi_354 = -xi_313 - xi_314 - xi_350 - xi_353 + xia_17_collide*2.0f - xia_8_collide;
      const float xi_368 = -xi_207 - xi_208 + xi_209 + xi_328 + xi_365 + xia_23_collide;
      const float xia_0_collide = xi_15;
      const float xia_1_collide = xi_2;
      const float xi_16 = xi_3*xia_1_collide;
      const float xi_17 = xi_16*2.0f;
      const float xi_18 = xi_16*8.0f;
      const float xi_19 = xi_18 + xi_6;
      const float xi_20 = xi_12 + xi_18;
      const float xi_21 = ((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_19 + xi_17 - xi_19*xia_4_collide): (xi_10*xi_20 + xi_17 + xi_20*xia_4_collide))));
      const float rr_0_a_collide = 0.0f;
      const float xi_29 = rr_0_a_collide*xi_28 + xia_19_collide*-0.16666666666666666f;
      const float xi_44 = rr_0_a_collide*xi_43 + xia_5_collide*-0.16666666666666666f;
      const float xi_53 = rr_0_a_collide*xi_52 + xia_20_collide*-0.16666666666666666f;
      const float xi_56 = rr_0_a_collide*0.041666666666666664f;
      const float xi_57 = -xi_28 + xi_56*xia_19_collide;
      const float xi_58 = xi_56*xia_5_collide;
      const float xi_59 = xi_43 - xi_58;
      const float xi_60 = xi_57 + xi_59;
      const float xi_72 = -xi_43 + xi_58;
      const float xi_73 = xi_57 + xi_72;
      const float xi_77 = xi_56*xia_20_collide;
      const float xi_78 = -xi_52 + xi_77;
      const float xi_79 = xi_57 + xi_78;
      const float xi_86 = xi_52 + xi_57 - xi_77;
      const float xi_88 = xi_59 + xi_78;
      const float xi_97 = xi_72 + xi_78;
      const float rr_1_a_collide = xi_21;
      const float xi_25 = rr_1_a_collide*0.5f;
      const float xi_66 = rr_1_a_collide*0.125f;
      const float xi_67 = xi_66*xia_19_collide;
      const float xi_75 = rr_1_a_collide*0.041666666666666664f;
      const float xi_193 = rr_1_a_collide*xi_192;
      const float xi_197 = rr_1_a_collide*xi_196;
      const float xi_211 = rr_1_a_collide*xi_210;
      const float xi_289 = xi_211*0.01984126984126984f;
      const float xi_324 = xi_193*0.050000000000000003f;
      const float xi_345 = xi_197*0.035714285714285712f;
      const float xi_362 = xi_197*0.021428571428571429f;
      const float xi_386 = xi_211*0.013888888888888888f;
      const float xi_394 = rr_1_a_collide*0.25f;
      const float xi_429 = xi_197*0.0071428571428571426f;
      const float xi_431 = xi_211*0.003968253968253968f;
      const float xi_446 = xi_197*0.017857142857142856f;
      const float rr_2_a_collide = xia_0_collide*(xi_21*-2.0f + 4.0f);
      const float xi_307 = rr_2_a_collide*0.25f;
      const float xi_316 = rr_2_a_collide*0.083333333333333329f;
      const float xi_317 = xi_315*xi_316;
      const float xi_355 = xi_316*xi_354;
      const float xi_369 = xi_316*xi_368;
      const float xi_413 = rr_2_a_collide*0.041666666666666664f;
      const float xi_414 = xi_315*xi_413;
      const float xi_416 = rr_2_a_collide*0.125f;
      const float xi_417 = xi_306*xi_416;
      const float xi_424 = xi_368*xi_413;
      const float xi_425 = xi_366*xi_416;
      const float xi_452 = xi_351*xi_416;
      const float xi_453 = xi_354*xi_413;
      const float rho_a_collide = xia_26_collide;
      const float xi_214 = kT*rho_a_collide;
      const float xi_215 = powf(xi_214*(1.0f - ((-xi_14 + 1.0f)*(-xi_14 + 1.0f))), 0.5f);
      const float xi_217 = xi_215*xi_216;
      const float xi_221 = xi_215*xi_220;
      const float xi_224 = xi_215*xi_222*xi_223;
      const float xi_228 = xi_215*xi_227;
      const float xi_290 = xi_212*xi_217;
      const float xi_291 = xi_290*0.11904761904761904f;
      const float xi_292 = xi_289 + xi_291;
      const float xi_294 = xi_215*0.5f;
      const float xi_295 = xi_293*xi_294;
      const float xi_309 = powf(xi_214*(1.0f - ((-xi_15*(xi_14*-2.0f + 4.0f) + 1.0f)*(-xi_15*(xi_14*-2.0f + 4.0f) + 1.0f))), 0.5f);
      const float xi_312 = xi_309*xi_311;
      const float xi_319 = xi_223*xi_309;
      const float xi_320 = xi_319*0.16666666666666666f;
      const float xi_321 = xi_318*xi_320;
      const float xi_322 = xi_317 + xi_321;
      const float xi_323 = xi_306*xi_307 + xi_308*xi_312 + xi_322;
      const float xi_325 = xi_218*xi_221;
      const float xi_326 = xi_325*0.10000000000000001f;
      const float xi_335 = xi_215*xi_334;
      const float xi_336 = xi_333*xi_335;
      const float xi_339 = xi_225*xi_228;
      const float xi_341 = -xi_289 - xi_291;
      const float xi_346 = xi_339*0.071428571428571425f;
      const float xi_347 = xi_332*xi_335;
      const float xi_357 = xi_320*xi_356;
      const float xi_358 = xi_355 + xi_357;
      const float xi_359 = xi_307*xi_351 + xi_312*xi_352 + xi_358;
      const float xi_363 = xi_339*0.042857142857142858f;
      const float xi_371 = xi_320*xi_370;
      const float xi_372 = xi_369 + xi_371;
      const float xi_373 = xi_307*xi_366 + xi_312*xi_367 + xi_372;
      const float xi_384 = xi_215*xi_293*0.25f;
      const float xi_387 = xi_290*0.083333333333333329f;
      const float xi_391 = -xi_355 - xi_357;
      const float xi_392 = xi_294*(random_0_3 - 0.5f);
      const float xi_412 = xi_294*(random_1_1 - 0.5f);
      const float xi_418 = xi_319*0.083333333333333329f;
      const float xi_419 = xi_318*xi_418;
      const float xi_421 = xi_309*xi_420;
      const float xi_422 = xi_308*xi_421;
      const float xi_426 = xi_370*xi_418;
      const float xi_427 = xi_367*xi_421;
      const float xi_428 = -xi_424 + xi_425 - xi_426 + xi_427;
      const float xi_430 = xi_339*0.014285714285714285f;
      const float xi_432 = xi_290*0.023809523809523808f;
      const float xi_433 = -xi_431 - xi_432;
      const float xi_439 = xi_424 - xi_425 + xi_426 - xi_427;
      const float xi_444 = -xi_384;
      const float xi_447 = xi_339*0.035714285714285712f;
      const float xi_449 = xi_294*(random_1_0 - 0.5f);
      const float xi_454 = xi_352*xi_421;
      const float xi_455 = xi_356*xi_418;
      const float xi_456 = -xi_452 + xi_453 - xi_454 + xi_455;
      const float xi_461 = xi_452 - xi_453 + xi_454 - xi_455;
      const float u_0_a_collide = xia_6_collide;
      const float xi_22 = u_0_a_collide*xia_5_collide;
      const float xi_34 = xi_22*0.16666666666666666f;
      const float xi_35 = xi_22*0.083333333333333329f;
      const float xi_36 = rr_1_a_collide*xi_35;
      const float xi_37 = xi_34 - xi_36;
      const float xi_40 = -xi_34 + xi_36;
      const float xi_42 = rr_1_a_collide*xi_34;
      const float xi_63 = u_0_a_collide*xi_62;
      const float xi_68 = u_0_a_collide*xi_67;
      const float xi_90 = u_0_a_collide*xia_20_collide;
      const float xi_91 = xi_90*0.25f;
      const float xi_94 = xi_66*xi_90;
      const float xi_198 = u_0_a_collide*u_0_a_collide;
      const float u_1_a_collide = xia_13_collide;
      const float xi_23 = u_1_a_collide*xia_19_collide;
      const float xi_26 = xi_23*0.16666666666666666f;
      const float xi_27 = rr_1_a_collide*xi_26;
      const float xi_45 = xi_23*0.083333333333333329f;
      const float xi_46 = rr_1_a_collide*xi_45;
      const float xi_47 = -xi_26 + xi_46;
      const float xi_49 = xi_26 - xi_46;
      const float xi_54 = xi_37 + xi_49;
      const float xi_55 = xi_40 + xi_47;
      const float xi_64 = u_1_a_collide*0.25f;
      const float xi_65 = xi_64*xia_5_collide;
      const float xi_69 = u_1_a_collide*xi_66;
      const float xi_70 = xi_69*xia_5_collide;
      const float xi_71 = xi_63 + xi_65 - xi_68 - xi_70;
      const float xi_74 = -xi_63 - xi_65 + xi_68 + xi_70;
      const float xi_81 = xi_64*xia_20_collide;
      const float xi_83 = xi_69*xia_20_collide;
      const float xi_199 = rho_a_collide*(u_1_a_collide*u_1_a_collide);
      const float xi_298 = xi_199 + xi_285 + xi_297;
      const float xi_393 = rho_a_collide*u_1_a_collide;
      const float xi_395 = xi_394*(u_0_a_collide*xi_393 + xi_203 + xi_296);
      const float xi_396 = -xi_392 - xi_395;
      const float xi_403 = xi_392 + xi_395;
      const float u_2_a_collide = xia_21_collide;
      const float xi_24 = u_2_a_collide*xia_20_collide;
      const float xi_30 = xi_24*0.16666666666666666f;
      const float xi_31 = xi_24*0.083333333333333329f;
      const float xi_32 = rr_1_a_collide*xi_31;
      const float xi_33 = xi_30 - xi_32;
      const float xi_38 = xi_33 + xi_37;
      const float xi_39 = -xi_30 + xi_32;
      const float xi_41 = xi_39 + xi_40;
      const float xi_48 = xi_39 + xi_47;
      const float xi_50 = xi_33 + xi_49;
      const float xi_51 = rr_1_a_collide*xi_30;
      const float xi_61 = rr_1_a_collide*u_2_a_collide*xia_20_collide*-0.041666666666666664f + xi_31 + xi_55;
      const float xi_76 = xi_24*xi_75 - xi_31 + xi_54;
      const float xi_80 = rr_1_a_collide*u_0_a_collide*xia_5_collide*-0.041666666666666664f + xi_35 + xi_48;
      const float xi_82 = u_2_a_collide*xi_62;
      const float xi_84 = u_2_a_collide*xi_67;
      const float xi_85 = -xi_81 - xi_82 + xi_83 + xi_84;
      const float xi_87 = xi_22*xi_75 - xi_35 + xi_50;
      const float xi_89 = rr_1_a_collide*u_1_a_collide*xia_19_collide*-0.041666666666666664f + xi_41 + xi_45;
      const float xi_92 = u_2_a_collide*xia_5_collide;
      const float xi_93 = xi_92*0.25f;
      const float xi_95 = xi_66*xi_92;
      const float xi_96 = xi_91 + xi_93 - xi_94 - xi_95;
      const float xi_98 = -xi_91 - xi_93 + xi_94 + xi_95;
      const float xi_99 = xi_81 + xi_82 - xi_83 - xi_84;
      const float xi_100 = xi_23*xi_75 + xi_38 - xi_45;
      const float xi_201 = rho_a_collide*(u_2_a_collide*u_2_a_collide);
      const float xi_202 = -xi_201 + xia_12_collide + xia_15_collide + xia_25_collide + xia_29_collide;
      const float xi_205 = rho_a_collide*xi_198 - xi_194 + xi_199 - xi_200 - xi_202 - xi_204 - xia_16_collide - xia_3_collide;
      const float xi_206 = rr_1_a_collide*xi_205;
      const float xi_299 = rr_1_a_collide*(xi_202 + xi_296 + xi_298 - xia_10_collide - xia_22_collide + xia_23_collide + xia_28_collide);
      const float xi_300 = xi_299*0.125f;
      const float xi_301 = -xi_295 - xi_300;
      const float xi_329 = rho_a_collide*xi_198*2.0f - xi_201 - xi_204 - xi_298 - xi_327 - xi_328 + xia_11_collide*2.0f - xia_12_collide + xia_14_collide*-2.0f - xia_15_collide + xia_16_collide*2.0f - xia_25_collide - xia_29_collide + xia_3_collide*2.0f + xia_7_collide*2.0f + xia_8_collide*-2.0f;
      const float xi_330 = rr_1_a_collide*xi_329;
      const float xi_331 = xi_330*0.041666666666666664f;
      const float xi_337 = xi_331 + xi_336;
      const float xi_338 = xi_324 + xi_326 + xi_337;
      const float xi_342 = xi_295 + xi_300;
      const float xi_343 = -xi_331 - xi_336;
      const float xi_344 = -xi_324 - xi_326 + xi_343;
      const float xi_385 = xi_299*0.0625f;
      const float xi_388 = xi_206*0.041666666666666664f + xi_224*0.083333333333333329f;
      const float xi_389 = xi_330*0.020833333333333332f + xi_347*0.083333333333333329f + xi_388;
      const float xi_390 = xi_322 + xi_384 + xi_385 + xi_386 + xi_387 + xi_389;
      const float xi_405 = -xi_317 - xi_321 + xi_384 + xi_385 + xi_386 + xi_387 + xi_389;
      const float xi_415 = xi_394*(u_2_a_collide*xi_393 + xi_304 - xia_11_collide + xia_7_collide);
      const float xi_423 = xi_412 - xi_414 + xi_415 + xi_417 - xi_419 + xi_422;
      const float xi_434 = xi_193*0.025000000000000001f + xi_325*0.050000000000000003f + xi_343 + xi_388 - xi_429 - xi_430 + xi_433;
      const float xi_440 = rr_1_a_collide*xi_192*-0.025000000000000001f + rr_1_a_collide*xi_205*-0.041666666666666664f + xi_215*xi_218*xi_220*-0.050000000000000003f + xi_215*xi_222*xi_223*-0.083333333333333329f + xi_337 + xi_429 + xi_430 + xi_431 + xi_432;
      const float xi_445 = -xi_385;
      const float xi_448 = xi_372 + xi_389 + xi_433 + xi_444 + xi_445 + xi_446 + xi_447;
      const float xi_450 = xi_394*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_349 - xia_12_collide + xia_29_collide);
      const float xi_451 = -xi_449 - xi_450;
      const float xi_460 = xi_449 + xi_450;
      const float xi_463 = xi_412 + xi_414 + xi_415 - xi_417 + xi_419 - xi_422;
      const float xi_466 = -xi_369 - xi_371 + xi_389 + xi_433 + xi_444 + xi_445 + xi_446 + xi_447;
      const float forceTerm_0_a_collide = xi_22*xi_25 - xi_22 + xi_23*xi_25 - xi_23 + xi_24*xi_25 - xi_24;
      const float forceTerm_1_a_collide = u_1_a_collide*xia_19_collide*0.33333333333333331f - xi_27 - xi_29 - xi_38;
      const float forceTerm_2_a_collide = xi_23*0.33333333333333331f - xi_27 + xi_29 + xi_41;
      const float forceTerm_3_a_collide = xi_22*0.33333333333333331f - xi_42 + xi_44 + xi_48;
      const float forceTerm_4_a_collide = u_0_a_collide*xia_5_collide*0.33333333333333331f - xi_42 - xi_44 - xi_50;
      const float forceTerm_5_a_collide = u_2_a_collide*xia_20_collide*0.33333333333333331f - xi_51 - xi_53 - xi_54;
      const float forceTerm_6_a_collide = xi_24*0.33333333333333331f - xi_51 + xi_53 + xi_55;
      const float forceTerm_7_a_collide = -xi_60 - xi_61 - xi_71;
      const float forceTerm_8_a_collide = -xi_61 - xi_73 - xi_74;
      const float forceTerm_9_a_collide = xi_71 + xi_73 + xi_76;
      const float forceTerm_10_a_collide = xi_60 + xi_74 + xi_76;
      const float forceTerm_11_a_collide = -xi_79 - xi_80 - xi_85;
      const float forceTerm_12_a_collide = xi_85 + xi_86 + xi_87;
      const float forceTerm_13_a_collide = -xi_88 - xi_89 - xi_96;
      const float forceTerm_14_a_collide = -xi_89 - xi_97 - xi_98;
      const float forceTerm_15_a_collide = -xi_80 - xi_86 - xi_99;
      const float forceTerm_16_a_collide = xi_79 + xi_87 + xi_99;
      const float forceTerm_17_a_collide = xi_100 + xi_96 + xi_97;
      const float forceTerm_18_a_collide = xi_100 + xi_88 + xi_98;
      const float xib_3_collide = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const bool xi_101 = xib_3_collide > 0.5f;
      const bool xi_102 = xib_3_collide < -0.5f;
      const float xi_103 = xib_3_collide*xib_3_collide;
      const bool xi_104 = xib_3_collide > 0.0f;
      const float xi_105 = ((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_8 + xi_5 + xi_9*xib_3_collide): (xi_103*xi_13 + xi_13*xib_3_collide + xi_5))));
      const float xi_106 = ((1.0f) / (xi_105*-0.25f + 2.0f));
      const float xib_4_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3];
      const float xi_564 = -xib_4_collide;
      const float xib_5_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3];
      const float xi_505 = -xib_5_collide;
      const float xib_6_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3];
      const float xib_7_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const float xib_8_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3];
      const float xi_533 = -xib_8_collide;
      const float xib_9_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3];
      const float xi_496 = -xib_9_collide;
      const float xib_10_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3];
      const float xib_11_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const float xib_12_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2];
      const float xi_470 = -xib_12_collide;
      const float xib_13_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3];
      const float xib_14_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3];
      const float xib_15_collide = _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xib_16_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3];
      const float xi_504 = -xib_10_collide - xib_16_collide;
      const float xi_520 = xib_10_collide*2.0f + xib_16_collide*-2.0f;
      const float xib_17_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const float xib_18_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3];
      const float xi_476 = xib_14_collide + xib_18_collide;
      const float xi_511 = xib_14_collide - xib_18_collide;
      const float xib_19_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2];
      const float xi_134 = xib_19_collide*0.083333333333333329f;
      const float xib_20_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3];
      const float xib_21_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3];
      const float xi_552 = xib_21_collide - xib_6_collide;
      const float xib_22_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3];
      const float xi_519 = xib_22_collide*2.0f;
      const float xib_23_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3];
      const float xi_483 = xib_13_collide*2.0f + xib_23_collide*2.0f;
      const float xi_551 = xib_13_collide - xib_23_collide;
      const float xib_24_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3];
      const float xi_534 = -xib_24_collide;
      const float xib_25_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3];
      const float xi_553 = xi_551 + xi_552 - xib_20_collide + xib_25_collide;
      const float xi_555 = xib_20_collide - xib_25_collide;
      const float xib_26_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3];
      const float xi_479 = xib_22_collide + xib_26_collide;
      const float xi_480 = xi_479 + xib_10_collide + xib_16_collide;
      const float xi_518 = xib_26_collide*2.0f;
      const float xi_556 = xi_518 - xi_519 + xi_520 + xi_552 + xi_555 - xib_13_collide + xib_23_collide;
      const float xib_27_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3];
      const float xi_143 = xib_27_collide*0.083333333333333329f;
      const float xib_28_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3];
      const float xi_119 = xib_28_collide*0.083333333333333329f;
      const float xi_154 = xib_28_collide*0.25f;
      const float xib_29_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3];
      const float xi_468 = -xib_12_collide + xib_14_collide*-3.0f + xib_18_collide*-3.0f + xib_24_collide*3.0f + xib_29_collide*-3.0f + xib_4_collide*-3.0f + xib_5_collide*3.0f + xib_8_collide*3.0f + xib_9_collide*3.0f;
      const float xi_471 = xi_470 + xib_14_collide*2.0f + xib_18_collide*2.0f + xib_21_collide*5.0f + xib_29_collide*2.0f + xib_4_collide*2.0f + xib_6_collide*5.0f;
      const float xi_472 = xi_471 + xib_13_collide*-5.0f + xib_20_collide*-5.0f + xib_23_collide*-5.0f + xib_24_collide*3.0f + xib_25_collide*-5.0f + xib_5_collide*-2.0f + xib_8_collide*3.0f + xib_9_collide*-2.0f;
      const float xi_484 = xi_471 + xi_483 + xib_10_collide*-7.0f + xib_16_collide*-7.0f + xib_20_collide*2.0f + xib_22_collide*-7.0f + xib_24_collide*-4.0f + xib_25_collide*2.0f + xib_26_collide*-7.0f + xib_5_collide*5.0f + xib_8_collide*-4.0f + xib_9_collide*5.0f;
      const float xi_510 = -xib_29_collide;
      const float xi_512 = xi_510 + xi_511 + xib_4_collide;
      const float xi_513 = -xi_496 - xi_512 - xib_5_collide;
      const float xi_521 = xi_505 + xi_512 - xi_518 + xi_519 + xi_520 + xib_9_collide;
      const float xi_565 = xi_476 + xi_510 + xi_564;
      const float xi_566 = xi_533 + xi_565 + xib_24_collide;
      const float xi_568 = -xi_483 - xi_534 - xi_565 + xib_20_collide*2.0f + xib_25_collide*2.0f - xib_8_collide;
      const float xib_0_collide = xi_106;
      const float xib_1_collide = xi_2;
      const float xi_107 = xi_3*xib_1_collide;
      const float xi_108 = xi_107*2.0f;
      const float xi_109 = xi_107*8.0f;
      const float xi_110 = xi_109 + xi_6;
      const float xi_111 = xi_109 + xi_12;
      const float xi_112 = ((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_110 + xi_108 - xi_110*xib_3_collide): (xi_103*xi_111 + xi_108 + xi_111*xib_3_collide))));
      const float rr_0_b_collide = 0.0f;
      const float xi_120 = rr_0_b_collide*xi_119 + xib_28_collide*-0.16666666666666666f;
      const float xi_135 = rr_0_b_collide*xi_134 + xib_19_collide*-0.16666666666666666f;
      const float xi_144 = rr_0_b_collide*xi_143 + xib_27_collide*-0.16666666666666666f;
      const float xi_147 = rr_0_b_collide*0.041666666666666664f;
      const float xi_148 = xi_147*xib_28_collide;
      const float xi_149 = xi_119 - xi_148;
      const float xi_150 = -xi_134 + xi_147*xib_19_collide;
      const float xi_151 = xi_149 + xi_150;
      const float xi_164 = -xi_119 + xi_148;
      const float xi_165 = xi_150 + xi_164;
      const float xi_168 = xi_147*xib_27_collide;
      const float xi_169 = -xi_143 + xi_168;
      const float xi_170 = xi_164 + xi_169;
      const float xi_177 = xi_149 + xi_169;
      const float xi_179 = xi_143 + xi_150 - xi_168;
      const float xi_188 = xi_150 + xi_169;
      const float rr_1_b_collide = xi_112;
      const float xi_116 = rr_1_b_collide*0.5f;
      const float xi_152 = rr_1_b_collide*0.041666666666666664f;
      const float xi_158 = rr_1_b_collide*0.125f;
      const float xi_159 = xi_158*xib_28_collide;
      const float xi_469 = rr_1_b_collide*xi_468;
      const float xi_473 = rr_1_b_collide*xi_472;
      const float xi_485 = rr_1_b_collide*xi_484;
      const float xi_497 = xi_485*0.01984126984126984f;
      const float xi_530 = xi_469*0.050000000000000003f;
      const float xi_548 = xi_473*0.035714285714285712f;
      const float xi_562 = xi_473*0.021428571428571429f;
      const float xi_576 = xi_485*0.013888888888888888f;
      const float xi_584 = rr_1_b_collide*0.25f;
      const float xi_605 = xi_473*0.0071428571428571426f;
      const float xi_607 = xi_485*0.003968253968253968f;
      const float xi_615 = xi_473*0.017857142857142856f;
      const float rr_2_b_collide = xib_0_collide*(xi_112*-2.0f + 4.0f);
      const float xi_514 = rr_2_b_collide*0.25f;
      const float xi_522 = rr_2_b_collide*0.083333333333333329f;
      const float xi_523 = xi_521*xi_522;
      const float xi_557 = xi_522*xi_556;
      const float xi_569 = xi_522*xi_568;
      const float xi_590 = rr_2_b_collide*0.041666666666666664f;
      const float xi_591 = xi_521*xi_590;
      const float xi_593 = rr_2_b_collide*0.125f;
      const float xi_594 = xi_513*xi_593;
      const float xi_600 = xi_568*xi_590;
      const float xi_601 = xi_566*xi_593;
      const float xi_621 = xi_553*xi_593;
      const float xi_622 = xi_556*xi_590;
      const float rho_b_collide = xib_15_collide;
      const float xi_486 = kT*rho_b_collide;
      const float xi_487 = powf(xi_486*(1.0f - ((-xi_105 + 1.0f)*(-xi_105 + 1.0f))), 0.5f);
      const float xi_489 = xi_487*xi_488;
      const float xi_491 = xi_487*xi_490;
      const float xi_493 = xi_223*xi_487*xi_492;
      const float xi_495 = xi_487*xi_494;
      const float xi_498 = xi_212*xi_489;
      const float xi_499 = xi_498*0.11904761904761904f;
      const float xi_500 = xi_497 + xi_499;
      const float xi_502 = xi_487*0.5f;
      const float xi_503 = xi_501*xi_502;
      const float xi_516 = powf(xi_486*(1.0f - ((-xi_106*(xi_105*-2.0f + 4.0f) + 1.0f)*(-xi_106*(xi_105*-2.0f + 4.0f) + 1.0f))), 0.5f);
      const float xi_517 = xi_311*xi_516;
      const float xi_525 = xi_223*xi_516;
      const float xi_526 = xi_525*0.16666666666666666f;
      const float xi_527 = xi_524*xi_526;
      const float xi_528 = xi_523 + xi_527;
      const float xi_529 = xi_513*xi_514 + xi_515*xi_517 + xi_528;
      const float xi_531 = xi_218*xi_491;
      const float xi_532 = xi_531*0.10000000000000001f;
      const float xi_539 = xi_487*xi_538;
      const float xi_540 = xi_333*xi_539;
      const float xi_543 = xi_225*xi_495;
      const float xi_544 = -xi_497 - xi_499;
      const float xi_549 = xi_543*0.071428571428571425f;
      const float xi_550 = xi_332*xi_539;
      const float xi_559 = xi_526*xi_558;
      const float xi_560 = xi_557 + xi_559;
      const float xi_561 = xi_514*xi_553 + xi_517*xi_554 + xi_560;
      const float xi_563 = xi_543*0.042857142857142858f;
      const float xi_571 = xi_526*xi_570;
      const float xi_572 = xi_569 + xi_571;
      const float xi_573 = xi_514*xi_566 + xi_517*xi_567 + xi_572;
      const float xi_574 = xi_487*xi_501*0.25f;
      const float xi_577 = xi_498*0.083333333333333329f;
      const float xi_581 = -xi_557 - xi_559;
      const float xi_582 = xi_502*(random_4_3_shifted - 0.5f);
      const float xi_589 = xi_502*(random_5_1_shifted - 0.5f);
      const float xi_595 = xi_525*0.083333333333333329f;
      const float xi_596 = xi_524*xi_595;
      const float xi_597 = xi_420*xi_516;
      const float xi_598 = xi_515*xi_597;
      const float xi_602 = xi_570*xi_595;
      const float xi_603 = xi_567*xi_597;
      const float xi_604 = -xi_600 + xi_601 - xi_602 + xi_603;
      const float xi_606 = xi_543*0.014285714285714285f;
      const float xi_608 = xi_498*0.023809523809523808f;
      const float xi_609 = -xi_607 - xi_608;
      const float xi_611 = xi_600 - xi_601 + xi_602 - xi_603;
      const float xi_613 = -xi_574;
      const float xi_616 = xi_543*0.035714285714285712f;
      const float xi_618 = xi_502*(random_5_0_shifted - 0.5f);
      const float xi_623 = xi_554*xi_597;
      const float xi_624 = xi_558*xi_595;
      const float xi_625 = -xi_621 + xi_622 - xi_623 + xi_624;
      const float xi_627 = xi_621 - xi_622 + xi_623 - xi_624;
      const float u_0_b_collide = xib_7_collide;
      const float xi_113 = u_0_b_collide*xib_19_collide;
      const float xi_125 = xi_113*0.16666666666666666f;
      const float xi_126 = xi_113*0.083333333333333329f;
      const float xi_127 = rr_1_b_collide*xi_126;
      const float xi_128 = xi_125 - xi_127;
      const float xi_131 = -xi_125 + xi_127;
      const float xi_133 = rr_1_b_collide*xi_125;
      const float xi_155 = u_0_b_collide*xi_154;
      const float xi_160 = u_0_b_collide*xi_159;
      const float xi_181 = u_0_b_collide*xib_27_collide;
      const float xi_182 = xi_181*0.25f;
      const float xi_185 = xi_158*xi_181;
      const float xi_474 = u_0_b_collide*u_0_b_collide;
      const float u_1_b_collide = xib_11_collide;
      const float xi_114 = u_1_b_collide*xib_28_collide;
      const float xi_117 = xi_114*0.16666666666666666f;
      const float xi_118 = rr_1_b_collide*xi_117;
      const float xi_136 = xi_114*0.083333333333333329f;
      const float xi_137 = rr_1_b_collide*xi_136;
      const float xi_138 = -xi_117 + xi_137;
      const float xi_140 = xi_117 - xi_137;
      const float xi_145 = xi_128 + xi_140;
      const float xi_146 = xi_131 + xi_138;
      const float xi_156 = u_1_b_collide*0.25f;
      const float xi_157 = xi_156*xib_19_collide;
      const float xi_161 = u_1_b_collide*xi_158;
      const float xi_162 = xi_161*xib_19_collide;
      const float xi_163 = -xi_155 - xi_157 + xi_160 + xi_162;
      const float xi_167 = xi_155 + xi_157 - xi_160 - xi_162;
      const float xi_172 = xi_156*xib_27_collide;
      const float xi_174 = xi_161*xib_27_collide;
      const float xi_475 = rho_b_collide*(u_1_b_collide*u_1_b_collide);
      const float xi_506 = xi_475 + xi_496 + xi_505;
      const float xi_583 = rho_b_collide*u_1_b_collide;
      const float xi_585 = xi_584*(u_0_b_collide*xi_583 + xi_479 + xi_504);
      const float xi_586 = -xi_582 - xi_585;
      const float xi_587 = xi_582 + xi_585;
      const float u_2_b_collide = xib_17_collide;
      const float xi_115 = u_2_b_collide*xib_27_collide;
      const float xi_121 = xi_115*0.16666666666666666f;
      const float xi_122 = xi_115*0.083333333333333329f;
      const float xi_123 = rr_1_b_collide*xi_122;
      const float xi_124 = xi_121 - xi_123;
      const float xi_129 = xi_124 + xi_128;
      const float xi_130 = -xi_121 + xi_123;
      const float xi_132 = xi_130 + xi_131;
      const float xi_139 = xi_130 + xi_138;
      const float xi_141 = xi_124 + xi_140;
      const float xi_142 = rr_1_b_collide*xi_121;
      const float xi_153 = xi_115*xi_152 - xi_122 + xi_145;
      const float xi_166 = rr_1_b_collide*u_2_b_collide*xib_27_collide*-0.041666666666666664f + xi_122 + xi_146;
      const float xi_171 = rr_1_b_collide*u_0_b_collide*xib_19_collide*-0.041666666666666664f + xi_126 + xi_139;
      const float xi_173 = u_2_b_collide*xi_154;
      const float xi_175 = u_2_b_collide*xi_159;
      const float xi_176 = -xi_172 - xi_173 + xi_174 + xi_175;
      const float xi_178 = xi_172 + xi_173 - xi_174 - xi_175;
      const float xi_180 = xi_114*xi_152 + xi_129 - xi_136;
      const float xi_183 = u_2_b_collide*xib_19_collide;
      const float xi_184 = xi_183*0.25f;
      const float xi_186 = xi_158*xi_183;
      const float xi_187 = -xi_182 - xi_184 + xi_185 + xi_186;
      const float xi_189 = rr_1_b_collide*u_1_b_collide*xib_28_collide*-0.041666666666666664f + xi_132 + xi_136;
      const float xi_190 = xi_113*xi_152 - xi_126 + xi_141;
      const float xi_191 = xi_182 + xi_184 - xi_185 - xi_186;
      const float xi_477 = rho_b_collide*(u_2_b_collide*u_2_b_collide);
      const float xi_478 = -xi_477 + xib_13_collide + xib_20_collide + xib_23_collide + xib_25_collide;
      const float xi_481 = rho_b_collide*xi_474 - xi_470 + xi_475 - xi_476 - xi_478 - xi_480 - xib_29_collide - xib_4_collide;
      const float xi_482 = rr_1_b_collide*xi_481;
      const float xi_507 = rr_1_b_collide*(xi_478 + xi_504 + xi_506 - xib_22_collide + xib_24_collide - xib_26_collide + xib_8_collide);
      const float xi_508 = xi_507*0.125f;
      const float xi_509 = -xi_503 - xi_508;
      const float xi_535 = rho_b_collide*xi_474*2.0f - xi_477 - xi_480 - xi_506 - xi_533 - xi_534 - xib_13_collide + xib_14_collide*2.0f + xib_18_collide*2.0f - xib_20_collide + xib_21_collide*-2.0f - xib_23_collide - xib_25_collide + xib_29_collide*2.0f + xib_4_collide*2.0f + xib_6_collide*-2.0f;
      const float xi_536 = rr_1_b_collide*xi_535;
      const float xi_537 = xi_536*0.041666666666666664f;
      const float xi_541 = xi_537 + xi_540;
      const float xi_542 = xi_530 + xi_532 + xi_541;
      const float xi_545 = xi_503 + xi_508;
      const float xi_546 = -xi_537 - xi_540;
      const float xi_547 = -xi_530 - xi_532 + xi_546;
      const float xi_575 = xi_507*0.0625f;
      const float xi_578 = xi_482*0.041666666666666664f + xi_493*0.083333333333333329f;
      const float xi_579 = xi_536*0.020833333333333332f + xi_550*0.083333333333333329f + xi_578;
      const float xi_580 = xi_528 + xi_574 + xi_575 + xi_576 + xi_577 + xi_579;
      const float xi_588 = -xi_523 - xi_527 + xi_574 + xi_575 + xi_576 + xi_577 + xi_579;
      const float xi_592 = xi_584*(u_2_b_collide*xi_583 + xi_511 + xi_564 + xib_29_collide);
      const float xi_599 = xi_589 - xi_591 + xi_592 + xi_594 - xi_596 + xi_598;
      const float xi_610 = xi_469*0.025000000000000001f + xi_531*0.050000000000000003f + xi_546 + xi_578 - xi_605 - xi_606 + xi_609;
      const float xi_612 = rr_1_b_collide*xi_468*-0.025000000000000001f + rr_1_b_collide*xi_481*-0.041666666666666664f + xi_218*xi_487*xi_490*-0.050000000000000003f + xi_223*xi_487*xi_492*-0.083333333333333329f + xi_541 + xi_605 + xi_606 + xi_607 + xi_608;
      const float xi_614 = -xi_575;
      const float xi_617 = xi_572 + xi_579 + xi_609 + xi_613 + xi_614 + xi_615 + xi_616;
      const float xi_619 = xi_584*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_551 + xi_555);
      const float xi_620 = -xi_618 - xi_619;
      const float xi_626 = xi_618 + xi_619;
      const float xi_628 = xi_589 + xi_591 + xi_592 - xi_594 + xi_596 - xi_598;
      const float xi_629 = -xi_569 - xi_571 + xi_579 + xi_609 + xi_613 + xi_614 + xi_615 + xi_616;
      const float forceTerm_0_b_collide = xi_113*xi_116 - xi_113 + xi_114*xi_116 - xi_114 + xi_115*xi_116 - xi_115;
      const float forceTerm_1_b_collide = u_1_b_collide*xib_28_collide*0.33333333333333331f - xi_118 - xi_120 - xi_129;
      const float forceTerm_2_b_collide = xi_114*0.33333333333333331f - xi_118 + xi_120 + xi_132;
      const float forceTerm_3_b_collide = xi_113*0.33333333333333331f - xi_133 + xi_135 + xi_139;
      const float forceTerm_4_b_collide = u_0_b_collide*xib_19_collide*0.33333333333333331f - xi_133 - xi_135 - xi_141;
      const float forceTerm_5_b_collide = u_2_b_collide*xib_27_collide*0.33333333333333331f - xi_142 - xi_144 - xi_145;
      const float forceTerm_6_b_collide = xi_115*0.33333333333333331f - xi_142 + xi_144 + xi_146;
      const float forceTerm_7_b_collide = xi_151 + xi_153 + xi_163;
      const float forceTerm_8_b_collide = -xi_163 - xi_165 - xi_166;
      const float forceTerm_9_b_collide = xi_153 + xi_165 + xi_167;
      const float forceTerm_10_b_collide = -xi_151 - xi_166 - xi_167;
      const float forceTerm_11_b_collide = -xi_170 - xi_171 - xi_176;
      const float forceTerm_12_b_collide = -xi_171 - xi_177 - xi_178;
      const float forceTerm_13_b_collide = xi_179 + xi_180 + xi_187;
      const float forceTerm_14_b_collide = -xi_187 - xi_188 - xi_189;
      const float forceTerm_15_b_collide = xi_176 + xi_177 + xi_190;
      const float forceTerm_16_b_collide = xi_170 + xi_178 + xi_190;
      const float forceTerm_17_b_collide = xi_180 + xi_188 + xi_191;
      const float forceTerm_18_b_collide = -xi_179 - xi_189 - xi_191;
      const float tmp_a0 = forceTerm_0_a_collide + xi_193*0.10000000000000001f + xi_197*0.042857142857142858f + xi_206*-0.5f + xi_211*0.023809523809523808f + xi_213*xi_217 + xi_219*xi_221 - xi_224 + xi_226*xi_228 + xi_284 + xia_24_collide;
      const float tmp_a1 = forceTerm_1_a_collide + rr_1_a_collide*xi_196*0.014285714285714285f + xi_215*xi_225*xi_227*0.028571428571428571f - xi_285 - xi_288 - xi_292 - xi_301 - xi_323 - xi_338;
      const float tmp_a2 = forceTerm_2_a_collide + xi_197*0.014285714285714285f + xi_323 + xi_339*0.028571428571428571f + xi_340 + xi_341 + xi_342 + xi_344 + xia_27_collide;
      const float tmp_a3 = forceTerm_3_a_collide + xi_330*0.083333333333333329f + xi_341 - xi_345 - xi_346 + xi_347*0.33333333333333331f + xi_348 + xi_359 + xia_8_collide;
      const float tmp_a4 = forceTerm_4_a_collide + rr_1_a_collide*xi_329*0.083333333333333329f + xi_215*xi_332*xi_334*0.33333333333333331f - xi_292 - xi_345 - xi_346 - xi_353 - xi_359 - xi_360;
      const float tmp_a5 = forceTerm_5_a_collide + rr_1_a_collide*xi_210*0.015873015873015872f + xi_212*xi_215*xi_216*0.095238095238095233f - xi_327 - xi_338 - xi_342 - xi_361 - xi_362 - xi_363 - xi_373;
      const float tmp_a6 = forceTerm_6_a_collide + xi_211*0.015873015873015872f + xi_290*0.095238095238095233f + xi_301 + xi_344 - xi_362 - xi_363 + xi_373 + xi_374 + xia_28_collide;
      const float tmp_a7 = forceTerm_7_a_collide + xi_383 + xi_390 + xi_391 + xi_396 + xia_22_collide;
      const float tmp_a8 = forceTerm_8_a_collide + xi_358 + xi_390 + xi_402 + xi_403 + xia_18_collide;
      const float tmp_a9 = forceTerm_9_a_collide + xi_391 + xi_403 + xi_404 + xi_405 + xia_17_collide;
      const float tmp_a10 = forceTerm_10_a_collide + xi_358 + xi_396 + xi_405 + xi_406 + xia_10_collide;
      const float tmp_a11 = forceTerm_11_a_collide + xi_411 + xi_423 + xi_428 + xi_434 + xia_11_collide;
      const float tmp_a12 = forceTerm_12_a_collide - xi_302 - xi_423 - xi_438 - xi_439 - xi_440;
      const float tmp_a13 = forceTerm_13_a_collide + xi_443 + xi_448 + xi_451 + xi_456 + xia_15_collide;
      const float tmp_a14 = forceTerm_14_a_collide + xi_448 + xi_459 + xi_460 + xi_461 + xia_25_collide;
      const float tmp_a15 = forceTerm_15_a_collide - xi_364 - xi_428 - xi_440 - xi_462 - xi_463;
      const float tmp_a16 = forceTerm_16_a_collide + xi_434 + xi_439 + xi_463 + xi_464 + xia_3_collide;
      const float tmp_a17 = forceTerm_17_a_collide + xi_456 + xi_460 + xi_465 + xi_466 + xia_12_collide;
      const float tmp_a18 = forceTerm_18_a_collide + xi_451 + xi_461 + xi_466 + xi_467 + xia_29_collide;
      const float tmp_b0 = forceTerm_0_b_collide + xi_213*xi_489 + xi_219*xi_491 + xi_226*xi_495 + xi_284 + xi_469*0.10000000000000001f + xi_473*0.042857142857142858f + xi_482*-0.5f + xi_485*0.023809523809523808f - xi_493 + xib_12_collide;
      const float tmp_b1 = forceTerm_1_b_collide + rr_1_b_collide*xi_472*0.014285714285714285f + xi_225*xi_487*xi_494*0.028571428571428571f - xi_288 - xi_496 - xi_500 - xi_509 - xi_529 - xi_542;
      const float tmp_b2 = forceTerm_2_b_collide + xi_340 + xi_473*0.014285714285714285f + xi_529 + xi_543*0.028571428571428571f + xi_544 + xi_545 + xi_547 + xib_5_collide;
      const float tmp_b3 = forceTerm_3_b_collide + xi_348 + xi_536*0.083333333333333329f + xi_544 - xi_548 - xi_549 + xi_550*0.33333333333333331f + xi_561 + xib_6_collide;
      const float tmp_b4 = forceTerm_4_b_collide + rr_1_b_collide*xi_535*0.083333333333333329f + xi_332*xi_487*xi_538*0.33333333333333331f - xi_360 - xi_500 - xi_548 - xi_549 - xi_561 + xib_21_collide;
      const float tmp_b5 = forceTerm_5_b_collide + rr_1_b_collide*xi_484*0.015873015873015872f + xi_212*xi_487*xi_488*0.095238095238095233f - xi_361 - xi_534 - xi_542 - xi_545 - xi_562 - xi_563 - xi_573;
      const float tmp_b6 = forceTerm_6_b_collide + xi_374 + xi_485*0.015873015873015872f + xi_498*0.095238095238095233f + xi_509 + xi_547 - xi_562 - xi_563 + xi_573 + xib_8_collide;
      const float tmp_b7 = forceTerm_7_b_collide + xi_383 + xi_580 + xi_581 + xi_586 + xib_26_collide;
      const float tmp_b8 = forceTerm_8_b_collide + xi_402 + xi_560 + xi_580 + xi_587 + xib_16_collide;
      const float tmp_b9 = forceTerm_9_b_collide + xi_404 + xi_581 + xi_587 + xi_588 + xib_10_collide;
      const float tmp_b10 = forceTerm_10_b_collide + xi_406 + xi_560 + xi_586 + xi_588 + xib_22_collide;
      const float tmp_b11 = forceTerm_11_b_collide + xi_411 + xi_599 + xi_604 + xi_610 + xib_4_collide;
      const float tmp_b12 = forceTerm_12_b_collide - xi_438 - xi_510 - xi_599 - xi_611 - xi_612;
      const float tmp_b13 = forceTerm_13_b_collide + xi_443 + xi_617 + xi_620 + xi_625 + xib_13_collide;
      const float tmp_b14 = forceTerm_14_b_collide + xi_459 + xi_617 + xi_626 + xi_627 + xib_23_collide;
      const float tmp_b15 = forceTerm_15_b_collide - xi_462 - xi_604 - xi_612 - xi_628 + xib_14_collide;
      const float tmp_b16 = forceTerm_16_b_collide + xi_464 + xi_610 + xi_611 + xi_628 + xib_18_collide;
      const float tmp_b17 = forceTerm_17_b_collide + xi_465 + xi_625 + xi_626 + xi_629 + xib_25_collide;
      const float tmp_b18 = forceTerm_18_b_collide + xi_467 + xi_620 + xi_627 + xi_629 + xib_20_collide;
      const float xirecolor_0 = tmp_a0 + tmp_b0;
      const float xirecolor_1 = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_2 = ((1.0f) / (xirecolor_1));
      const float xi_651 = xirecolor_2*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const float xirecolor_4 = tmp_a1 + tmp_b1;
      const float xirecolor_5 = xi_239;
      const float xirecolor_6 = xi_237;
      const float xirecolor_7 = xi_257;
      const float xi_630 = -xirecolor_7;
      const float xirecolor_8 = xi_630;
      const float xirecolor_9 = xi_244;
      const float xirecolor_10 = -xirecolor_9;
      const float xirecolor_11 = xi_246;
      const float xi_631 = xirecolor_10 + xirecolor_11;
      const float xirecolor_12 = xi_259;
      const float xi_632 = xirecolor_12 + xirecolor_8;
      const float xirecolor_13 = xi_631 + xi_632;
      const float xirecolor_14 = xi_261;
      const float xirecolor_15 = xi_262;
      const float xi_633 = xirecolor_14 - xirecolor_15;
      const float xirecolor_16 = -xi_633;
      const float xirecolor_17 = xi_231;
      const float xirecolor_18 = xi_232;
      const float xirecolor_19 = xi_233;
      const float xi_634 = xirecolor_18 + xirecolor_19;
      const float xirecolor_20 = xi_234;
      const float xirecolor_21 = xi_634 - xirecolor_17 + xirecolor_20;
      const float xirecolor_22 = xi_230;
      const float xirecolor_23 = xi_229;
      const float xi_635 = xirecolor_22 - xirecolor_23;
      const float xirecolor_24 = -xi_635;
      const float xirecolor_25 = xi_265;
      const float xirecolor_26 = -xirecolor_25;
      const float xirecolor_27 = xi_251;
      const float xirecolor_28 = xi_250;
      const float xirecolor_29 = xi_267;
      const float xi_636 = xirecolor_26 + xirecolor_29;
      const float xirecolor_30 = xi_636 - xirecolor_27 + xirecolor_28;
      const float xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
      const float xirecolor_32 = xi_236;
      const float xirecolor_33 = xi_238;
      const float xi_637 = xirecolor_32 + xirecolor_33;
      const float xirecolor_34 = xi_637 + xirecolor_5 - xirecolor_6;
      const float xirecolor_35 = xi_241;
      const float xirecolor_36 = xi_242;
      const float xi_638 = xirecolor_35 - xirecolor_36;
      const float xirecolor_37 = -xi_638;
      const float xirecolor_38 = xi_631;
      const float xirecolor_39 = xi_248;
      const float xirecolor_40 = xi_249;
      const float xi_639 = xirecolor_39 + xirecolor_40;
      const float xirecolor_41 = xirecolor_27 - xirecolor_28;
      const float xirecolor_42 = xi_639 + xirecolor_41;
      const float xirecolor_43 = xirecolor_38 + xirecolor_42;
      const float xirecolor_44 = xi_635 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
      const float xirecolor_45 = xi_633;
      const float xirecolor_46 = xi_632;
      const float xirecolor_47 = xi_273;
      const float xirecolor_48 = -xirecolor_47;
      const float xirecolor_49 = xi_275;
      const float xi_640 = xirecolor_48 + xirecolor_49;
      const float xirecolor_50 = xi_638 + xi_640;
      const float xirecolor_51 = xi_634 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
      const float xirecolor_52 = powf(xirecolor_31*xirecolor_31 + xirecolor_44*xirecolor_44 + xirecolor_51*xirecolor_51, 0.5f);
      const float xirecolor_53 = ((1.0f) / (xirecolor_52));
      const float xi_641 = xirecolor_31*xirecolor_53;
      const float xi_642 = xirecolor_51*xirecolor_53;
      const float xi_643 = xirecolor_44*xirecolor_53;
      const bool xirecolor_54 = xirecolor_52 > 1.0000000000000001e-5f;
      const float xirecolor_55 = beta*((1.0f) / (xirecolor_1*xirecolor_1))*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2]*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_56 = xirecolor_55*(0.055555555555555552f*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.055555555555555552f*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const float xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_641): (0.0f));
      const float xirecolor_58 = tmp_a2 + tmp_b2;
      const float xirecolor_59 = xirecolor_56*((xirecolor_54) ? (-xi_641): (0.0f));
      const float xirecolor_60 = tmp_a3 + tmp_b3;
      const float xirecolor_61 = xirecolor_56*((xirecolor_54) ? (-xi_642): (0.0f));
      const float xirecolor_62 = tmp_a4 + tmp_b4;
      const float xirecolor_63 = xirecolor_56*((xirecolor_54) ? (xi_642): (0.0f));
      const float xirecolor_64 = tmp_a5 + tmp_b5;
      const float xirecolor_65 = xirecolor_56*((xirecolor_54) ? (-xi_643): (0.0f));
      const float xirecolor_66 = tmp_a6 + tmp_b6;
      const float xirecolor_67 = xirecolor_56*((xirecolor_54) ? (xi_643): (0.0f));
      const float xirecolor_68 = tmp_a7 + tmp_b7;
      const float xirecolor_69 = xi_375;
      const float xirecolor_70 = xi_376;
      const float xirecolor_71 = xi_637;
      const float xirecolor_72 = xirecolor_50 + xirecolor_71;
      const float xirecolor_73 = xi_380;
      const float xirecolor_74 = xirecolor_25 - xirecolor_29 + xirecolor_73;
      const float xirecolor_75 = xi_377 - xirecolor_11 + xirecolor_41 - xirecolor_69 + xirecolor_70 + xirecolor_72 + xirecolor_74 + xirecolor_9;
      const float xirecolor_76 = xi_311*xirecolor_53;
      const float xi_644 = xirecolor_75*xirecolor_76;
      const float xirecolor_77 = xirecolor_55*(0.027777777777777776f*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.027777777777777776f*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const float xirecolor_78 = xirecolor_77*((xirecolor_54) ? (-xi_644): (0.0f));
      const float xirecolor_79 = tmp_a8 + tmp_b8;
      const float xirecolor_80 = xi_397;
      const float xirecolor_81 = xi_398;
      const float xi_645 = -xirecolor_80 + xirecolor_81;
      const float xirecolor_82 = xi_399;
      const float xirecolor_83 = xi_400 + xi_645 + xirecolor_30 + xirecolor_38 + xirecolor_72 + xirecolor_82;
      const float xi_646 = xirecolor_76*xirecolor_83;
      const float xirecolor_84 = xirecolor_77*((xirecolor_54) ? (xi_646): (0.0f));
      const float xirecolor_85 = tmp_a9 + tmp_b9;
      const float xirecolor_86 = xirecolor_77*((xirecolor_54) ? (-xi_646): (0.0f));
      const float xirecolor_87 = tmp_a10 + tmp_b10;
      const float xirecolor_88 = xirecolor_77*((xirecolor_54) ? (xi_644): (0.0f));
      const float xirecolor_89 = tmp_a11 + tmp_b11;
      const float xirecolor_90 = xi_639 + xirecolor_37 + xirecolor_71;
      const float xirecolor_91 = -xi_630 - xi_645 - xirecolor_12;
      const float xirecolor_92 = xi_407 + xirecolor_45 + xirecolor_74 + xirecolor_90 + xirecolor_91;
      const float xi_647 = xirecolor_76*xirecolor_92;
      const float xirecolor_93 = xirecolor_77*((xirecolor_54) ? (-xi_647): (0.0f));
      const float xirecolor_94 = tmp_a12 + tmp_b12;
      const float xirecolor_95 = xirecolor_16 + xirecolor_69 - xirecolor_70;
      const float xirecolor_96 = xi_435 + xi_636 + xirecolor_46 + xirecolor_82 + xirecolor_90 + xirecolor_95;
      const float xi_648 = xirecolor_76*xirecolor_96;
      const float xirecolor_97 = xirecolor_77*((xirecolor_54) ? (-xi_648): (0.0f));
      const float xirecolor_98 = tmp_a13 + tmp_b13;
      const float xirecolor_99 = xi_441 + xi_640 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_73 + xirecolor_82;
      const float xi_649 = xirecolor_76*xirecolor_99;
      const float xirecolor_100 = xirecolor_77*((xirecolor_54) ? (-xi_649): (0.0f));
      const float xirecolor_101 = tmp_a14 + tmp_b14;
      const float xirecolor_102 = xi_457 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_91 + xirecolor_95;
      const float xi_650 = xirecolor_102*xirecolor_76;
      const float xirecolor_103 = xirecolor_77*((xirecolor_54) ? (-xi_650): (0.0f));
      const float xirecolor_104 = tmp_a15 + tmp_b15;
      const float xirecolor_105 = xirecolor_77*((xirecolor_54) ? (xi_648): (0.0f));
      const float xirecolor_106 = tmp_a16 + tmp_b16;
      const float xirecolor_107 = xirecolor_77*((xirecolor_54) ? (xi_647): (0.0f));
      const float xirecolor_108 = tmp_a17 + tmp_b17;
      const float xirecolor_109 = xirecolor_77*((xirecolor_54) ? (xi_650): (0.0f));
      const float xirecolor_110 = tmp_a18 + tmp_b18;
      const float xirecolor_111 = xirecolor_77*((xirecolor_54) ? (xi_649): (0.0f));
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2] = xirecolor_0*xirecolor_3;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3] = xirecolor_3*xirecolor_4 + xirecolor_57;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3] = xirecolor_3*xirecolor_58 + xirecolor_59;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3] = xirecolor_3*xirecolor_60 + xirecolor_61;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3] = xirecolor_3*xirecolor_62 + xirecolor_63;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3] = xirecolor_3*xirecolor_64 + xirecolor_65;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3] = xirecolor_3*xirecolor_66 + xirecolor_67;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3] = xirecolor_3*xirecolor_68 + xirecolor_78;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3] = xirecolor_3*xirecolor_79 + xirecolor_84;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3] = xirecolor_3*xirecolor_85 + xirecolor_86;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3] = xirecolor_3*xirecolor_87 + xirecolor_88;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3] = xirecolor_3*xirecolor_89 + xirecolor_93;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3] = xirecolor_3*xirecolor_94 + xirecolor_97;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3] = xirecolor_100 + xirecolor_3*xirecolor_98;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3] = xirecolor_101*xirecolor_3 + xirecolor_103;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3] = xirecolor_104*xirecolor_3 + xirecolor_105;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3] = xirecolor_106*xirecolor_3 + xirecolor_107;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3] = xirecolor_108*xirecolor_3 + xirecolor_109;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3] = xirecolor_110*xirecolor_3 + xirecolor_111;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2] = xi_651*xirecolor_0;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3] = xi_651*xirecolor_4 - xirecolor_57;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3] = xi_651*xirecolor_58 - xirecolor_59;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3] = xi_651*xirecolor_60 - xirecolor_61;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3] = xi_651*xirecolor_62 - xirecolor_63;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3] = xi_651*xirecolor_64 - xirecolor_65;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3] = xi_651*xirecolor_66 - xirecolor_67;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3] = xi_651*xirecolor_68 - xirecolor_78;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3] = xi_651*xirecolor_79 - xirecolor_84;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3] = xi_651*xirecolor_85 - xirecolor_86;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3] = xi_651*xirecolor_87 - xirecolor_88;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3] = xi_651*xirecolor_89 - xirecolor_93;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3] = xi_651*xirecolor_94 - xirecolor_97;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3] = xi_651*xirecolor_98 - xirecolor_100;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3] = xi_651*xirecolor_101 - xirecolor_103;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3] = xi_651*xirecolor_104 - xirecolor_105;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3] = xi_651*xirecolor_106 - xirecolor_107;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3] = xi_651*xirecolor_108 - xirecolor_109;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3] = xi_651*xirecolor_110 - xirecolor_111;
   } 
}
}


void ColorGradientCollideSweepSinglePrecisionCUDA::run(IBlock * block, gpuStream_t stream)
{
   if (!this->configured_)
      WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")
         
    auto velocity = block->getData< gpu::GPUField<float> >(velocityID);
    auto force_a = block->getData< gpu::GPUField<float> >(force_aID);
    auto pdfs_b = block->getData< gpu::GPUField<float> >(pdfs_bID);
    auto rho_a = block->getData< gpu::GPUField<float> >(rho_aID);
    auto pdfs_a = block->getData< gpu::GPUField<float> >(pdfs_aID);
    auto force_b = block->getData< gpu::GPUField<float> >(force_bID);
    auto phasefield = block->getData< gpu::GPUField<float> >(phasefieldID);
    auto rho_b = block->getData< gpu::GPUField<float> >(rho_bID);

    auto & block_offset_2 = this->block_offset_2_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & beta = this->beta_;
    auto & seed = this->seed_;
    auto & time_step = this->time_step_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & block_offset_0 = this->block_offset_0_;
    auto & sigma = this->sigma_;
    auto & block_offset_1 = this->block_offset_1_;
    auto & kT = this->kT_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(force_a->xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(force_a->xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(force_a->ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(force_a->ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(force_a->zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(force_a->zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    const int64_t _stride_force_a_0 = int64_t(force_a->xStride());
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_0 = int64_t(force_b->xStride());
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_0 = int64_t(pdfs_a->xStride());
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_0 = int64_t(pdfs_b->xStride());
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_0 = int64_t(phasefield->xStride());
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_0 = int64_t(rho_a->xStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_0 = int64_t(rho_b->xStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    dim3 _block(uint32_c(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)), uint32_c(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))), uint32_c(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))));
    dim3 _grid(uint32_c(( (_size_force_a_0 - 2) % (((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) == 0 ? (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) : ( (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) ) +1 )), uint32_c(( (_size_force_a_1 - 2) % (((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) == 0 ? (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) : ( (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) ) +1 )), uint32_c(( (_size_force_a_2 - 2) % (((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) == 0 ? (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) : ( (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) ) +1 )));
    internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda::colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, block_offset_0, block_offset_1, block_offset_2, kT, omega_shear_a, omega_shear_b, seed, sigma, time_step);
    
}


void ColorGradientCollideSweepSinglePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block, gpuStream_t stream)
{
   if (!this->configured_)
      WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")
         
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto velocity = block->getData< gpu::GPUField<float> >(velocityID);
    auto force_a = block->getData< gpu::GPUField<float> >(force_aID);
    auto pdfs_b = block->getData< gpu::GPUField<float> >(pdfs_bID);
    auto rho_a = block->getData< gpu::GPUField<float> >(rho_aID);
    auto pdfs_a = block->getData< gpu::GPUField<float> >(pdfs_aID);
    auto force_b = block->getData< gpu::GPUField<float> >(force_bID);
    auto phasefield = block->getData< gpu::GPUField<float> >(phasefieldID);
    auto rho_b = block->getData< gpu::GPUField<float> >(rho_bID);

    auto & block_offset_2 = this->block_offset_2_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & beta = this->beta_;
    auto & seed = this->seed_;
    auto & time_step = this->time_step_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & block_offset_0 = this->block_offset_0_;
    auto & sigma = this->sigma_;
    auto & block_offset_1 = this->block_offset_1_;
    auto & kT = this->kT_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(ci.xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(ci.ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(ci.zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    const int64_t _stride_force_a_0 = int64_t(force_a->xStride());
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_0 = int64_t(force_b->xStride());
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_0 = int64_t(pdfs_a->xStride());
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_0 = int64_t(pdfs_b->xStride());
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_0 = int64_t(phasefield->xStride());
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_0 = int64_t(rho_a->xStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_0 = int64_t(rho_b->xStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    dim3 _block(uint32_c(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)), uint32_c(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))), uint32_c(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))));
    dim3 _grid(uint32_c(( (_size_force_a_0 - 2) % (((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) == 0 ? (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) : ( (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) ) +1 )), uint32_c(( (_size_force_a_1 - 2) % (((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) == 0 ? (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) : ( (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) ) +1 )), uint32_c(( (_size_force_a_2 - 2) % (((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) == 0 ? (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) : ( (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) ) +1 )));
    internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda::colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, block_offset_0, block_offset_1, block_offset_2, kT, omega_shear_a, omega_shear_b, seed, sigma, time_step);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
