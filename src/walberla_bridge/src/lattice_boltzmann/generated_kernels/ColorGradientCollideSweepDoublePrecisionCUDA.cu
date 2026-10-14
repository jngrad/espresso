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
//! \\file ColorGradientCollideSweepDoublePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepDoublePrecisionCUDA.h"


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


namespace internal_colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda(double * RESTRICT const _data_force_a, double * RESTRICT const _data_force_b, double * RESTRICT  _data_pdfs_a, double * RESTRICT  _data_pdfs_b, double * RESTRICT const _data_phasefield, double * RESTRICT const _data_rho_a, double * RESTRICT const _data_rho_b, double * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_0, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_0, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_0, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_0, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_0, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_0, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_0, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, double beta, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double kT, double omega_shear_a, double omega_shear_b, uint32_t seed, double sigma, uint32_t time_step)
{
   if (blockDim.x*blockIdx.x + threadIdx.x + 1 < _size_force_a_0 - 1 && blockDim.y*blockIdx.y + threadIdx.y + 1 < _size_force_a_1 - 1 && blockDim.z*blockIdx.z + threadIdx.z + 1 < _size_force_a_2 - 1)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x + 1;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y + 1;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z + 1;
      
      double random_7_0{};
      double random_7_1{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 7, seed, random_7_0, random_7_1);
      }
      
      double random_6_0{};
      double random_6_1{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 6, seed, random_6_0, random_6_1);
      }
      
      double random_5_0{};
      double random_5_1{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 5, seed, random_5_0, random_5_1);
      }
      
      double random_4_0{};
      double random_4_1{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 4, seed, random_4_0, random_4_1);
      }
      
      double random_3_0{};
      double random_3_1{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, random_3_0, random_3_1);
      }
      
      double random_2_0{};
      double random_2_1{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, random_2_0, random_2_1);
      }
      
      double random_1_0{};
      double random_1_1{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, random_1_0, random_1_1);
      }
      
      double random_0_0{};
      double random_0_1{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, random_0_0, random_0_1);
      }
      
      double random_15_0_shifted{};
      double random_15_1_shifted{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 23, seed, random_15_0_shifted, random_15_1_shifted);
      }
      
      double random_14_0_shifted{};
      double random_14_1_shifted{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 22, seed, random_14_0_shifted, random_14_1_shifted);
      }
      
      double random_13_0_shifted{};
      double random_13_1_shifted{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 21, seed, random_13_0_shifted, random_13_1_shifted);
      }
      
      double random_12_0_shifted{};
      double random_12_1_shifted{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 20, seed, random_12_0_shifted, random_12_1_shifted);
      }
      
      double random_11_0_shifted{};
      double random_11_1_shifted{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 19, seed, random_11_0_shifted, random_11_1_shifted);
      }
      
      double random_10_0_shifted{};
      double random_10_1_shifted{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 18, seed, random_10_0_shifted, random_10_1_shifted);
      }
      
      double random_9_0_shifted{};
      double random_9_1_shifted{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 17, seed, random_9_0_shifted, random_9_1_shifted);
      }
      
      double random_8_0_shifted{};
      double random_8_1_shifted{};
      if (kT > 0.) {  
      philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 16, seed, random_8_0_shifted, random_8_1_shifted);
      }
      const double xi_2 = ((1.0) / (omega_shear_a + omega_shear_b));
      const double xi_3 = omega_shear_a*omega_shear_b;
      const double xi_4 = xi_2*xi_3;
      const double xi_5 = xi_4*2.0;
      const double xi_6 = omega_shear_a*-4.0;
      const double xi_7 = xi_4*8.0;
      const double xi_8 = xi_6 + xi_7;
      const double xi_9 = -xi_8;
      const double xi_12 = omega_shear_b*-4.0;
      const double xi_13 = xi_12 + xi_7;
      const double xi_212 = 3.7416573867739413;
      const double xi_213 = xi_212*0.14285714285714285;
      const double xi_216 = random_6_0 - 0.5;
      const double xi_218 = 5.4772255750516612;
      const double xi_219 = xi_218*0.20000000000000001;
      const double xi_220 = random_7_0 - 0.5;
      const double xi_222 = random_0_0 - 0.5;
      const double xi_223 = 2.4494897427831779;
      const double xi_225 = 8.3666002653407556;
      const double xi_226 = xi_225*0.085714285714285715;
      const double xi_227 = random_6_1 - 0.5;
      const double xi_229 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_230 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_231 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_232 = -0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_233 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_234 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_235 = -xi_231 + xi_232 + xi_233 + xi_234;
      const double xi_236 = -0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_237 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_238 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_239 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_240 = xi_236 - xi_237 + xi_238 + xi_239;
      const double xi_241 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_242 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_243 = -xi_241 + xi_242;
      const double xi_244 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_245 = -xi_244;
      const double xi_246 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_247 = xi_245 + xi_246;
      const double xi_248 = -0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_249 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_250 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_251 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_252 = -xi_250 + xi_251;
      const double xi_253 = xi_248 + xi_249 + xi_252;
      const double xi_254 = xi_247 + xi_253;
      const double xi_255 = -xi_229 + xi_230 + xi_235 + xi_240 + xi_243 + xi_254;
      const double xi_256 = xi_255*xi_255;
      const double xi_257 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_258 = -xi_257;
      const double xi_259 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_260 = xi_245 + xi_246 + xi_258 + xi_259;
      const double xi_261 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_262 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_263 = -xi_261 + xi_262;
      const double xi_264 = xi_229 - xi_230;
      const double xi_265 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_266 = -xi_265;
      const double xi_267 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_268 = xi_250 - xi_251 + xi_266 + xi_267;
      const double xi_269 = xi_235 + xi_237 - xi_239 + xi_260 + xi_263 + xi_264 + xi_268;
      const double xi_270 = xi_269*xi_269;
      const double xi_271 = xi_261 - xi_262;
      const double xi_272 = xi_258 + xi_259;
      const double xi_273 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const double xi_274 = -xi_273;
      const double xi_275 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const double xi_276 = xi_241 - xi_242 + xi_274 + xi_275;
      const double xi_277 = xi_231 + xi_232 + xi_233 - xi_234 + xi_240 + xi_264 + xi_271 + xi_272 + xi_276;
      const double xi_278 = xi_277*xi_277;
      const double xi_279 = xi_256 + xi_270 + xi_278;
      const double xi_280 = pow(xi_279, 0.5);
      const double xi_281 = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2]*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const double xi_282 = sigma*xi_280*((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] > 0.5) ? (omega_shear_a): ((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] < -0.5) ? (omega_shear_b): ((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] > 0.0) ? (xi_281*xi_8 + xi_5 + xi_9*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2]): (xi_13*xi_281 + xi_13*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] + xi_5))));
      const bool xi_283 = xi_280 > 9.9999999999999998e-13;
      const double xi_284 = ((xi_283) ? (xi_282*0.25): (0.0));
      const double xi_286 = ((1.0) / (xi_279));
      const double xi_287 = xi_282*1.125;
      const double xi_288 = -((xi_283) ? (xi_287*(xi_270*xi_286*0.055555555555555552 - 0.018518518518518517)): (0.0));
      const double xi_293 = random_1_0 - 0.5;
      const double xi_307 = random_5_1 - 0.5;
      const double xi_309 = 1.4142135623730951;
      const double xi_310 = xi_309*0.5;
      const double xi_317 = random_3_0 - 0.5;
      const double xi_331 = 1.7320508075688772;
      const double xi_332 = xi_331*0.16666666666666666;
      const double xi_333 = random_0_1 - 0.5;
      const double xi_339 = ((xi_283) ? (xi_287*(xi_286*0.055555555555555552*(xi_269*xi_269) - 0.018518518518518517)): (0.0));
      const double xi_347 = ((xi_283) ? (xi_287*(xi_286*0.055555555555555552*(xi_277*xi_277) - 0.018518518518518517)): (0.0));
      const double xi_352 = random_4_1 - 0.5;
      const double xi_355 = random_4_0 - 0.5;
      const double xi_359 = -((xi_283) ? (xi_287*(xi_278*xi_286*0.055555555555555552 - 0.018518518518518517)): (0.0));
      const double xi_360 = -((xi_283) ? (xi_287*(xi_286*0.055555555555555552*(xi_255*xi_255) - 0.018518518518518517)): (0.0));
      const double xi_366 = random_5_0 - 0.5;
      const double xi_369 = random_3_1 - 0.5;
      const double xi_373 = ((xi_283) ? (xi_287*(xi_256*xi_286*0.055555555555555552 - 0.018518518518518517)): (0.0));
      const double xi_374 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_375 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_376 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_377 = xi_236 + xi_238;
      const double xi_378 = xi_276 + xi_377;
      const double xi_379 = -0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_380 = xi_265 - xi_267 + xi_379;
      const double xi_381 = xi_244 - xi_246 + xi_252 - xi_374 + xi_375 + xi_376 + xi_378 + xi_380;
      const double xi_382 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_381*xi_381) - 0.037037037037037035)): (0.0));
      const double xi_396 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_397 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_398 = -0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_399 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_400 = xi_247 + xi_268 + xi_378 - xi_396 + xi_397 + xi_398 + xi_399;
      const double xi_401 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_400*xi_400) - 0.037037037037037035)): (0.0));
      const double xi_403 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_400*xi_400) - 0.037037037037037035)): (0.0));
      const double xi_405 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_381*xi_381) - 0.037037037037037035)): (0.0));
      const double xi_406 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_407 = xi_243 + xi_248 + xi_249 + xi_377;
      const double xi_408 = xi_257 - xi_259 + xi_396 - xi_397;
      const double xi_409 = xi_271 + xi_380 + xi_406 + xi_407 + xi_408;
      const double xi_410 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_409*xi_409) - 0.037037037037037035)): (0.0));
      const double xi_419 = xi_309*0.25;
      const double xi_434 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_435 = xi_263 + xi_374 - xi_375;
      const double xi_436 = xi_266 + xi_267 + xi_272 + xi_398 + xi_407 + xi_434 + xi_435;
      const double xi_437 = -((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_436*xi_436) - 0.037037037037037035)): (0.0));
      const double xi_440 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_441 = xi_253 + xi_260 + xi_271 + xi_274 + xi_275 + xi_379 + xi_398 + xi_440;
      const double xi_442 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_441*xi_441) - 0.037037037037037035)): (0.0));
      const double xi_456 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_457 = xi_254 + xi_273 - xi_275 + xi_408 + xi_435 + xi_456;
      const double xi_458 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_457*xi_457) - 0.037037037037037035)): (0.0));
      const double xi_461 = -((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_436*xi_436) - 0.037037037037037035)): (0.0));
      const double xi_463 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_409*xi_409) - 0.037037037037037035)): (0.0));
      const double xi_464 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_457*xi_457) - 0.037037037037037035)): (0.0));
      const double xi_466 = ((xi_283) ? (xi_287*(xi_286*0.027777777777777776*(xi_441*xi_441) - 0.037037037037037035)): (0.0));
      const double xi_489 = random_14_0_shifted - 0.5;
      const double xi_491 = random_15_0_shifted - 0.5;
      const double xi_493 = random_8_0_shifted - 0.5;
      const double xi_495 = random_14_1_shifted - 0.5;
      const double xi_502 = random_9_0_shifted - 0.5;
      const double xi_517 = random_13_1_shifted - 0.5;
      const double xi_526 = random_11_0_shifted - 0.5;
      const double xi_540 = random_8_1_shifted - 0.5;
      const double xi_556 = random_12_1_shifted - 0.5;
      const double xi_560 = random_12_0_shifted - 0.5;
      const double xi_569 = random_13_0_shifted - 0.5;
      const double xi_572 = random_11_1_shifted - 0.5;
      const double xia_3_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3];
      const double xi_363 = -xia_3_collide;
      const double xia_4_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3];
      const double xia_5_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3];
      const double xia_6_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3];
      const double xi_203 = xia_4_collide + xia_6_collide;
      const double xi_313 = xia_4_collide*2.0 + xia_6_collide*-2.0;
      const double xia_7_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3];
      const double xia_8_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3];
      const double xi_348 = -xia_8_collide;
      const double xia_9_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3];
      const double xi_28 = xia_9_collide*0.083333333333333329;
      const double xi_64 = xia_9_collide*0.25;
      const double xia_10_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3];
      const double xi_52 = xia_10_collide*0.083333333333333329;
      const double xia_11_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2];
      const double xi_194 = -xia_11_collide;
      const double xia_12_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2];
      const double xi_43 = xia_12_collide*0.083333333333333329;
      const double xia_13_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3];
      const double xia_14_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const double xia_15_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3];
      const double xi_200 = xia_15_collide + xia_7_collide;
      const double xi_303 = xia_15_collide - xia_7_collide;
      const double xia_16_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3];
      const double xi_195 = xi_194 + xia_15_collide*2.0 + xia_16_collide*2.0 + xia_3_collide*2.0 + xia_5_collide*5.0 + xia_7_collide*2.0 + xia_8_collide*5.0;
      const double xi_302 = -xia_16_collide;
      const double xi_304 = xi_302 + xi_303 + xia_3_collide;
      const double xi_364 = xi_200 + xi_302 + xi_363;
      const double xia_17_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3];
      const double xi_297 = -xia_17_collide;
      const double xia_18_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const double xia_19_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3];
      const double xi_207 = xia_19_collide*2.0;
      const double xia_20_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3];
      const double xi_285 = -xia_20_collide;
      const double xi_305 = xi_297 + xi_304 + xia_20_collide;
      const double xia_21_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3];
      const double xi_208 = xia_21_collide*2.0;
      const double xia_22_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3];
      const double xi_312 = xia_22_collide*2.0;
      const double xia_23_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const double xia_24_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3];
      const double xi_327 = -xia_24_collide;
      const double xia_25_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3];
      const double xi_192 = -xia_11_collide + xia_15_collide*-3.0 + xia_16_collide*-3.0 + xia_17_collide*3.0 + xia_20_collide*3.0 + xia_24_collide*3.0 + xia_25_collide*3.0 + xia_3_collide*-3.0 + xia_7_collide*-3.0;
      const double xi_326 = -xia_25_collide;
      const double xi_365 = -xi_326 - xi_364 - xia_24_collide;
      const double xia_26_collide = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const bool xi_0 = xia_26_collide > 0.5;
      const bool xi_1 = xia_26_collide < -0.5;
      const double xi_10 = xia_26_collide*xia_26_collide;
      const bool xi_11 = xia_26_collide > 0.0;
      const double xi_14 = ((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_8 + xi_5 + xi_9*xia_26_collide): (xi_10*xi_13 + xi_13*xia_26_collide + xi_5))));
      const double xi_15 = ((1.0) / (xi_14*-0.25 + 2.0));
      const double xia_27_collide = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const double xia_28_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3];
      const double xi_196 = xi_195 + xia_13_collide*-5.0 + xia_17_collide*-2.0 + xia_19_collide*-5.0 + xia_20_collide*-2.0 + xia_21_collide*-5.0 + xia_24_collide*3.0 + xia_25_collide*3.0 + xia_28_collide*-5.0;
      const double xi_209 = xia_13_collide*2.0 + xia_28_collide*2.0;
      const double xi_349 = xia_13_collide - xia_28_collide;
      const double xi_350 = xi_349 - xia_19_collide + xia_21_collide;
      const double xi_351 = -xi_348 - xi_350 - xia_5_collide;
      const double xi_367 = -xi_207 - xi_208 + xi_209 + xi_327 + xi_364 + xia_25_collide;
      const double xia_29_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3];
      const double xi_204 = xi_203 + xia_22_collide + xia_29_collide;
      const double xi_210 = xi_195 + xi_207 + xi_208 + xi_209 + xia_17_collide*5.0 + xia_20_collide*5.0 + xia_22_collide*-7.0 + xia_24_collide*-4.0 + xia_25_collide*-4.0 + xia_29_collide*-7.0 + xia_4_collide*-7.0 + xia_6_collide*-7.0;
      const double xi_296 = -xia_22_collide - xia_29_collide;
      const double xi_314 = -xi_285 - xi_304 - xi_312 - xi_313 - xia_17_collide + xia_29_collide*2.0;
      const double xi_353 = -xi_312 + xi_313 + xi_350 + xia_29_collide*2.0 - xia_5_collide + xia_8_collide;
      const double xia_0_collide = xi_15;
      const double xia_1_collide = xi_2;
      const double xi_16 = xi_3*xia_1_collide;
      const double xi_17 = xi_16*2.0;
      const double xi_18 = xi_16*8.0;
      const double xi_19 = xi_18 + xi_6;
      const double xi_20 = xi_12 + xi_18;
      const double xi_21 = ((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_19 + xi_17 - xi_19*xia_26_collide): (xi_10*xi_20 + xi_17 + xi_20*xia_26_collide))));
      const double rr_0_a_collide = 0.0;
      const double xi_29 = rr_0_a_collide*xi_28 + xia_9_collide*-0.16666666666666666;
      const double xi_44 = rr_0_a_collide*xi_43 + xia_12_collide*-0.16666666666666666;
      const double xi_53 = rr_0_a_collide*xi_52 + xia_10_collide*-0.16666666666666666;
      const double xi_56 = rr_0_a_collide*0.041666666666666664;
      const double xi_57 = xi_56*xia_9_collide;
      const double xi_58 = xi_28 - xi_57;
      const double xi_59 = xi_56*xia_12_collide;
      const double xi_60 = -xi_43 + xi_59;
      const double xi_61 = xi_58 + xi_60;
      const double xi_74 = -xi_28 + xi_57;
      const double xi_75 = xi_60 + xi_74;
      const double xi_78 = -xi_52 + xi_56*xia_10_collide;
      const double xi_79 = xi_74 + xi_78;
      const double xi_86 = xi_58 + xi_78;
      const double xi_88 = xi_43 - xi_59 + xi_78;
      const double xi_97 = xi_60 + xi_78;
      const double rr_1_a_collide = xi_21;
      const double xi_25 = rr_1_a_collide*0.5;
      const double xi_62 = rr_1_a_collide*0.041666666666666664;
      const double xi_68 = rr_1_a_collide*0.125;
      const double xi_69 = xi_68*xia_9_collide;
      const double xi_193 = rr_1_a_collide*xi_192;
      const double xi_197 = rr_1_a_collide*xi_196;
      const double xi_211 = rr_1_a_collide*xi_210;
      const double xi_289 = xi_211*0.01984126984126984;
      const double xi_323 = xi_193*0.050000000000000003;
      const double xi_344 = xi_197*0.035714285714285712;
      const double xi_361 = xi_197*0.021428571428571429;
      const double xi_385 = xi_211*0.013888888888888888;
      const double xi_393 = rr_1_a_collide*0.25;
      const double xi_428 = xi_197*0.0071428571428571426;
      const double xi_430 = xi_211*0.003968253968253968;
      const double xi_445 = xi_197*0.017857142857142856;
      const double rr_2_a_collide = xia_0_collide*(xi_21*-2.0 + 4.0);
      const double xi_306 = rr_2_a_collide*0.25;
      const double xi_315 = rr_2_a_collide*0.083333333333333329;
      const double xi_316 = xi_314*xi_315;
      const double xi_354 = xi_315*xi_353;
      const double xi_368 = xi_315*xi_367;
      const double xi_412 = rr_2_a_collide*0.041666666666666664;
      const double xi_413 = xi_314*xi_412;
      const double xi_415 = rr_2_a_collide*0.125;
      const double xi_416 = xi_305*xi_415;
      const double xi_423 = xi_367*xi_412;
      const double xi_424 = xi_365*xi_415;
      const double xi_451 = xi_351*xi_415;
      const double xi_452 = xi_353*xi_412;
      const double rho_a_collide = xia_27_collide;
      const double xi_214 = kT*rho_a_collide;
      const double xi_215 = pow(xi_214*(1.0 - ((-xi_14 + 1.0)*(-xi_14 + 1.0))), 0.5);
      const double xi_217 = xi_215*xi_216;
      const double xi_221 = xi_215*xi_220;
      const double xi_224 = xi_215*xi_222*xi_223;
      const double xi_228 = xi_215*xi_227;
      const double xi_290 = xi_212*xi_217;
      const double xi_291 = xi_290*0.11904761904761904;
      const double xi_292 = xi_289 + xi_291;
      const double xi_294 = xi_215*0.5;
      const double xi_295 = xi_293*xi_294;
      const double xi_308 = pow(xi_214*(1.0 - ((-xi_15*(xi_14*-2.0 + 4.0) + 1.0)*(-xi_15*(xi_14*-2.0 + 4.0) + 1.0))), 0.5);
      const double xi_311 = xi_308*xi_310;
      const double xi_318 = xi_223*xi_308;
      const double xi_319 = xi_318*0.16666666666666666;
      const double xi_320 = xi_317*xi_319;
      const double xi_321 = xi_316 + xi_320;
      const double xi_322 = xi_305*xi_306 + xi_307*xi_311 + xi_321;
      const double xi_324 = xi_218*xi_221;
      const double xi_325 = xi_324*0.10000000000000001;
      const double xi_334 = xi_215*xi_333;
      const double xi_335 = xi_332*xi_334;
      const double xi_338 = xi_225*xi_228;
      const double xi_340 = -xi_289 - xi_291;
      const double xi_345 = xi_338*0.071428571428571425;
      const double xi_346 = xi_331*xi_334;
      const double xi_356 = xi_319*xi_355;
      const double xi_357 = xi_354 + xi_356;
      const double xi_358 = xi_306*xi_351 + xi_311*xi_352 + xi_357;
      const double xi_362 = xi_338*0.042857142857142858;
      const double xi_370 = xi_319*xi_369;
      const double xi_371 = xi_368 + xi_370;
      const double xi_372 = xi_306*xi_365 + xi_311*xi_366 + xi_371;
      const double xi_383 = xi_215*xi_293*0.25;
      const double xi_386 = xi_290*0.083333333333333329;
      const double xi_390 = -xi_354 - xi_356;
      const double xi_391 = xi_294*(random_1_1 - 0.5);
      const double xi_411 = xi_294*(random_2_1 - 0.5);
      const double xi_417 = xi_318*0.083333333333333329;
      const double xi_418 = xi_317*xi_417;
      const double xi_420 = xi_308*xi_419;
      const double xi_421 = xi_307*xi_420;
      const double xi_425 = xi_369*xi_417;
      const double xi_426 = xi_366*xi_420;
      const double xi_427 = -xi_423 + xi_424 - xi_425 + xi_426;
      const double xi_429 = xi_338*0.014285714285714285;
      const double xi_431 = xi_290*0.023809523809523808;
      const double xi_432 = -xi_430 - xi_431;
      const double xi_438 = xi_423 - xi_424 + xi_425 - xi_426;
      const double xi_443 = -xi_383;
      const double xi_446 = xi_338*0.035714285714285712;
      const double xi_448 = xi_294*(random_2_0 - 0.5);
      const double xi_453 = xi_352*xi_420;
      const double xi_454 = xi_355*xi_417;
      const double xi_455 = -xi_451 + xi_452 - xi_453 + xi_454;
      const double xi_460 = xi_451 - xi_452 + xi_453 - xi_454;
      const double u_0_a_collide = xia_23_collide;
      const double xi_22 = u_0_a_collide*xia_12_collide;
      const double xi_34 = xi_22*0.16666666666666666;
      const double xi_35 = xi_22*0.083333333333333329;
      const double xi_36 = rr_1_a_collide*xi_35;
      const double xi_37 = xi_34 - xi_36;
      const double xi_40 = -xi_34 + xi_36;
      const double xi_42 = rr_1_a_collide*xi_34;
      const double xi_65 = u_0_a_collide*xi_64;
      const double xi_70 = u_0_a_collide*xi_69;
      const double xi_90 = u_0_a_collide*xia_10_collide;
      const double xi_91 = xi_90*0.25;
      const double xi_94 = xi_68*xi_90;
      const double xi_198 = u_0_a_collide*u_0_a_collide;
      const double u_1_a_collide = xia_14_collide;
      const double xi_23 = u_1_a_collide*xia_9_collide;
      const double xi_26 = xi_23*0.16666666666666666;
      const double xi_27 = rr_1_a_collide*xi_26;
      const double xi_45 = xi_23*0.083333333333333329;
      const double xi_46 = rr_1_a_collide*xi_45;
      const double xi_47 = -xi_26 + xi_46;
      const double xi_49 = xi_26 - xi_46;
      const double xi_54 = xi_37 + xi_49;
      const double xi_55 = xi_40 + xi_47;
      const double xi_66 = u_1_a_collide*0.25;
      const double xi_67 = xi_66*xia_12_collide;
      const double xi_71 = u_1_a_collide*xi_68;
      const double xi_72 = xi_71*xia_12_collide;
      const double xi_73 = -xi_65 - xi_67 + xi_70 + xi_72;
      const double xi_77 = xi_65 + xi_67 - xi_70 - xi_72;
      const double xi_81 = xi_66*xia_10_collide;
      const double xi_83 = xi_71*xia_10_collide;
      const double xi_199 = rho_a_collide*(u_1_a_collide*u_1_a_collide);
      const double xi_298 = xi_199 + xi_285 + xi_297;
      const double xi_392 = rho_a_collide*u_1_a_collide;
      const double xi_394 = xi_393*(u_0_a_collide*xi_392 + xi_203 + xi_296);
      const double xi_395 = -xi_391 - xi_394;
      const double xi_402 = xi_391 + xi_394;
      const double u_2_a_collide = xia_18_collide;
      const double xi_24 = u_2_a_collide*xia_10_collide;
      const double xi_30 = xi_24*0.16666666666666666;
      const double xi_31 = xi_24*0.083333333333333329;
      const double xi_32 = rr_1_a_collide*xi_31;
      const double xi_33 = xi_30 - xi_32;
      const double xi_38 = xi_33 + xi_37;
      const double xi_39 = -xi_30 + xi_32;
      const double xi_41 = xi_39 + xi_40;
      const double xi_48 = xi_39 + xi_47;
      const double xi_50 = xi_33 + xi_49;
      const double xi_51 = rr_1_a_collide*xi_30;
      const double xi_63 = xi_24*xi_62 - xi_31 + xi_54;
      const double xi_76 = rr_1_a_collide*u_2_a_collide*xia_10_collide*-0.041666666666666664 + xi_31 + xi_55;
      const double xi_80 = rr_1_a_collide*u_0_a_collide*xia_12_collide*-0.041666666666666664 + xi_35 + xi_48;
      const double xi_82 = u_2_a_collide*xi_64;
      const double xi_84 = u_2_a_collide*xi_69;
      const double xi_85 = -xi_81 - xi_82 + xi_83 + xi_84;
      const double xi_87 = xi_81 + xi_82 - xi_83 - xi_84;
      const double xi_89 = rr_1_a_collide*u_1_a_collide*xia_9_collide*-0.041666666666666664 + xi_41 + xi_45;
      const double xi_92 = u_2_a_collide*xia_12_collide;
      const double xi_93 = xi_92*0.25;
      const double xi_95 = xi_68*xi_92;
      const double xi_96 = xi_91 + xi_93 - xi_94 - xi_95;
      const double xi_98 = -xi_91 - xi_93 + xi_94 + xi_95;
      const double xi_99 = xi_22*xi_62 - xi_35 + xi_50;
      const double xi_100 = xi_23*xi_62 + xi_38 - xi_45;
      const double xi_201 = rho_a_collide*(u_2_a_collide*u_2_a_collide);
      const double xi_202 = -xi_201 + xia_13_collide + xia_19_collide + xia_21_collide + xia_28_collide;
      const double xi_205 = rho_a_collide*xi_198 - xi_194 + xi_199 - xi_200 - xi_202 - xi_204 - xia_16_collide - xia_3_collide;
      const double xi_206 = rr_1_a_collide*xi_205;
      const double xi_299 = rr_1_a_collide*(xi_202 + xi_296 + xi_298 + xia_24_collide + xia_25_collide - xia_4_collide - xia_6_collide);
      const double xi_300 = xi_299*0.125;
      const double xi_301 = -xi_295 - xi_300;
      const double xi_328 = rho_a_collide*xi_198*2.0 - xi_201 - xi_204 - xi_298 - xi_326 - xi_327 - xia_13_collide + xia_15_collide*2.0 + xia_16_collide*2.0 - xia_19_collide - xia_21_collide - xia_28_collide + xia_3_collide*2.0 + xia_5_collide*-2.0 + xia_7_collide*2.0 + xia_8_collide*-2.0;
      const double xi_329 = rr_1_a_collide*xi_328;
      const double xi_330 = xi_329*0.041666666666666664;
      const double xi_336 = xi_330 + xi_335;
      const double xi_337 = xi_323 + xi_325 + xi_336;
      const double xi_341 = xi_295 + xi_300;
      const double xi_342 = -xi_330 - xi_335;
      const double xi_343 = -xi_323 - xi_325 + xi_342;
      const double xi_384 = xi_299*0.0625;
      const double xi_387 = xi_206*0.041666666666666664 + xi_224*0.083333333333333329;
      const double xi_388 = xi_329*0.020833333333333332 + xi_346*0.083333333333333329 + xi_387;
      const double xi_389 = xi_321 + xi_383 + xi_384 + xi_385 + xi_386 + xi_388;
      const double xi_404 = -xi_316 - xi_320 + xi_383 + xi_384 + xi_385 + xi_386 + xi_388;
      const double xi_414 = xi_393*(u_2_a_collide*xi_392 + xi_303 + xi_363 + xia_16_collide);
      const double xi_422 = xi_411 - xi_413 + xi_414 + xi_416 - xi_418 + xi_421;
      const double xi_433 = xi_193*0.025000000000000001 + xi_324*0.050000000000000003 + xi_342 + xi_387 - xi_428 - xi_429 + xi_432;
      const double xi_439 = rr_1_a_collide*xi_192*-0.025000000000000001 + rr_1_a_collide*xi_205*-0.041666666666666664 + xi_215*xi_218*xi_220*-0.050000000000000003 + xi_215*xi_222*xi_223*-0.083333333333333329 + xi_336 + xi_428 + xi_429 + xi_430 + xi_431;
      const double xi_444 = -xi_384;
      const double xi_447 = xi_371 + xi_388 + xi_432 + xi_443 + xi_444 + xi_445 + xi_446;
      const double xi_449 = xi_393*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_349 + xia_19_collide - xia_21_collide);
      const double xi_450 = -xi_448 - xi_449;
      const double xi_459 = xi_448 + xi_449;
      const double xi_462 = xi_411 + xi_413 + xi_414 - xi_416 + xi_418 - xi_421;
      const double xi_465 = -xi_368 - xi_370 + xi_388 + xi_432 + xi_443 + xi_444 + xi_445 + xi_446;
      const double forceTerm_0_a_collide = xi_22*xi_25 - xi_22 + xi_23*xi_25 - xi_23 + xi_24*xi_25 - xi_24;
      const double forceTerm_1_a_collide = u_1_a_collide*xia_9_collide*0.33333333333333331 - xi_27 - xi_29 - xi_38;
      const double forceTerm_2_a_collide = xi_23*0.33333333333333331 - xi_27 + xi_29 + xi_41;
      const double forceTerm_3_a_collide = xi_22*0.33333333333333331 - xi_42 + xi_44 + xi_48;
      const double forceTerm_4_a_collide = u_0_a_collide*xia_12_collide*0.33333333333333331 - xi_42 - xi_44 - xi_50;
      const double forceTerm_5_a_collide = u_2_a_collide*xia_10_collide*0.33333333333333331 - xi_51 - xi_53 - xi_54;
      const double forceTerm_6_a_collide = xi_24*0.33333333333333331 - xi_51 + xi_53 + xi_55;
      const double forceTerm_7_a_collide = xi_61 + xi_63 + xi_73;
      const double forceTerm_8_a_collide = -xi_73 - xi_75 - xi_76;
      const double forceTerm_9_a_collide = xi_63 + xi_75 + xi_77;
      const double forceTerm_10_a_collide = -xi_61 - xi_76 - xi_77;
      const double forceTerm_11_a_collide = -xi_79 - xi_80 - xi_85;
      const double forceTerm_12_a_collide = -xi_80 - xi_86 - xi_87;
      const double forceTerm_13_a_collide = -xi_88 - xi_89 - xi_96;
      const double forceTerm_14_a_collide = -xi_89 - xi_97 - xi_98;
      const double forceTerm_15_a_collide = xi_85 + xi_86 + xi_99;
      const double forceTerm_16_a_collide = xi_79 + xi_87 + xi_99;
      const double forceTerm_17_a_collide = xi_100 + xi_96 + xi_97;
      const double forceTerm_18_a_collide = xi_100 + xi_88 + xi_98;
      const double xib_3_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3];
      const double xi_506 = -xib_3_collide;
      const double xib_4_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3];
      const double xi_511 = -xib_4_collide;
      const double xib_5_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2];
      const double xi_134 = xib_5_collide*0.083333333333333329;
      const double xib_6_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2];
      const double xi_469 = -xib_6_collide;
      const double xib_7_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3];
      const double xi_512 = -xib_7_collide;
      const double xib_8_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const double xib_9_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3];
      const double xi_143 = xib_9_collide*0.083333333333333329;
      const double xib_10_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3];
      const double xib_11_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3];
      const double xi_475 = xib_11_collide + xib_4_collide;
      const double xib_12_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3];
      const double xi_497 = -xib_12_collide;
      const double xib_13_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const double xib_14_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3];
      const double xi_513 = xi_512 + xib_14_collide;
      const double xi_514 = xi_511 + xi_513 + xib_11_collide;
      const double xi_515 = -xi_497 - xi_514 - xib_3_collide;
      const double xi_566 = -xib_14_collide;
      const double xi_567 = xi_475 + xi_512 + xi_566;
      const double xib_15_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3];
      const double xi_482 = xib_15_collide*2.0;
      const double xib_16_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3];
      const double xi_521 = xib_16_collide*2.0;
      const double xib_17_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3];
      const double xi_483 = xib_17_collide*2.0;
      const double xi_553 = -xib_15_collide + xib_17_collide;
      const double xib_18_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const double xib_19_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3];
      const double xi_470 = xi_469 + xib_10_collide*5.0 + xib_11_collide*2.0 + xib_14_collide*2.0 + xib_19_collide*5.0 + xib_4_collide*2.0 + xib_7_collide*2.0;
      const double xi_554 = xib_10_collide - xib_19_collide;
      const double xib_20_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3];
      const double xi_478 = xib_16_collide + xib_20_collide;
      const double xi_520 = xib_20_collide*2.0;
      const double xib_21_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3];
      const double xi_536 = -xib_21_collide;
      const double xib_22_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3];
      const double xib_23_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3];
      const double xi_484 = xib_22_collide*2.0 + xib_23_collide*2.0;
      const double xi_555 = xi_553 + xi_554 - xib_22_collide + xib_23_collide;
      const double xi_557 = xib_22_collide - xib_23_collide;
      const double xib_24_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3];
      const double xi_467 = xib_11_collide*-3.0 + xib_12_collide*3.0 + xib_14_collide*-3.0 + xib_21_collide*3.0 + xib_24_collide*3.0 + xib_3_collide*3.0 + xib_4_collide*-3.0 - xib_6_collide + xib_7_collide*-3.0;
      const double xi_471 = xi_470 + xib_12_collide*-2.0 + xib_15_collide*-5.0 + xib_17_collide*-5.0 + xib_21_collide*3.0 + xib_22_collide*-5.0 + xib_23_collide*-5.0 + xib_24_collide*3.0 + xib_3_collide*-2.0;
      const double xi_535 = -xib_24_collide;
      const double xi_568 = -xi_535 - xi_567 - xib_21_collide;
      const double xi_570 = -xi_482 - xi_483 + xi_484 + xi_536 + xi_567 + xib_24_collide;
      const double xib_25_collide = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const bool xi_101 = xib_25_collide > 0.5;
      const bool xi_102 = xib_25_collide < -0.5;
      const double xi_103 = xib_25_collide*xib_25_collide;
      const bool xi_104 = xib_25_collide > 0.0;
      const double xi_105 = ((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_8 + xi_5 + xi_9*xib_25_collide): (xi_103*xi_13 + xi_13*xib_25_collide + xi_5))));
      const double xi_106 = ((1.0) / (xi_105*-0.25 + 2.0));
      const double xib_26_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3];
      const double xib_27_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3];
      const double xi_479 = xi_478 + xib_26_collide + xib_27_collide;
      const double xi_485 = xi_470 + xi_482 + xi_483 + xi_484 + xib_12_collide*5.0 + xib_16_collide*-7.0 + xib_20_collide*-7.0 + xib_21_collide*-4.0 + xib_24_collide*-4.0 + xib_26_collide*-7.0 + xib_27_collide*-7.0 + xib_3_collide*5.0;
      const double xi_505 = -xib_26_collide - xib_27_collide;
      const double xi_522 = xib_26_collide*2.0 + xib_27_collide*-2.0;
      const double xi_523 = xi_506 + xi_514 - xi_520 + xi_521 + xi_522 + xib_12_collide;
      const double xi_558 = xi_520 - xi_521 + xi_522 + xi_554 + xi_557 + xib_15_collide - xib_17_collide;
      const double xib_28_collide = _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const double xib_29_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3];
      const double xi_119 = xib_29_collide*0.083333333333333329;
      const double xi_152 = xib_29_collide*0.25;
      const double xib_0_collide = xi_106;
      const double xib_1_collide = xi_2;
      const double xi_107 = xi_3*xib_1_collide;
      const double xi_108 = xi_107*2.0;
      const double xi_109 = xi_107*8.0;
      const double xi_110 = xi_109 + xi_6;
      const double xi_111 = xi_109 + xi_12;
      const double xi_112 = ((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_110 + xi_108 - xi_110*xib_25_collide): (xi_103*xi_111 + xi_108 + xi_111*xib_25_collide))));
      const double rr_0_b_collide = 0.0;
      const double xi_120 = rr_0_b_collide*xi_119 + xib_29_collide*-0.16666666666666666;
      const double xi_135 = rr_0_b_collide*xi_134 + xib_5_collide*-0.16666666666666666;
      const double xi_144 = rr_0_b_collide*xi_143 + xib_9_collide*-0.16666666666666666;
      const double xi_147 = rr_0_b_collide*0.041666666666666664;
      const double xi_148 = xi_147*xib_5_collide;
      const double xi_149 = -xi_119 + xi_147*xib_29_collide;
      const double xi_150 = xi_134 - xi_148 + xi_149;
      const double xi_162 = -xi_134 + xi_148;
      const double xi_163 = xi_149 + xi_162;
      const double xi_167 = xi_147*xib_9_collide;
      const double xi_168 = -xi_143 + xi_167;
      const double xi_169 = xi_149 + xi_168;
      const double xi_176 = xi_143 - xi_167;
      const double xi_177 = xi_149 + xi_176;
      const double xi_179 = xi_162 + xi_176;
      const double xi_188 = xi_162 + xi_168;
      const double rr_1_b_collide = xi_112;
      const double xi_116 = rr_1_b_collide*0.5;
      const double xi_156 = rr_1_b_collide*0.125;
      const double xi_157 = xi_156*xib_29_collide;
      const double xi_165 = rr_1_b_collide*0.041666666666666664;
      const double xi_468 = rr_1_b_collide*xi_467;
      const double xi_472 = rr_1_b_collide*xi_471;
      const double xi_486 = rr_1_b_collide*xi_485;
      const double xi_498 = xi_486*0.01984126984126984;
      const double xi_532 = xi_468*0.050000000000000003;
      const double xi_550 = xi_472*0.035714285714285712;
      const double xi_564 = xi_472*0.021428571428571429;
      const double xi_578 = xi_486*0.013888888888888888;
      const double xi_586 = rr_1_b_collide*0.25;
      const double xi_607 = xi_472*0.0071428571428571426;
      const double xi_609 = xi_486*0.003968253968253968;
      const double xi_617 = xi_472*0.017857142857142856;
      const double rr_2_b_collide = xib_0_collide*(xi_112*-2.0 + 4.0);
      const double xi_516 = rr_2_b_collide*0.25;
      const double xi_524 = rr_2_b_collide*0.083333333333333329;
      const double xi_525 = xi_523*xi_524;
      const double xi_559 = xi_524*xi_558;
      const double xi_571 = xi_524*xi_570;
      const double xi_592 = rr_2_b_collide*0.041666666666666664;
      const double xi_593 = xi_523*xi_592;
      const double xi_595 = rr_2_b_collide*0.125;
      const double xi_596 = xi_515*xi_595;
      const double xi_602 = xi_570*xi_592;
      const double xi_603 = xi_568*xi_595;
      const double xi_623 = xi_555*xi_595;
      const double xi_624 = xi_558*xi_592;
      const double rho_b_collide = xib_28_collide;
      const double xi_487 = kT*rho_b_collide;
      const double xi_488 = pow(xi_487*(1.0 - ((-xi_105 + 1.0)*(-xi_105 + 1.0))), 0.5);
      const double xi_490 = xi_488*xi_489;
      const double xi_492 = xi_488*xi_491;
      const double xi_494 = xi_223*xi_488*xi_493;
      const double xi_496 = xi_488*xi_495;
      const double xi_499 = xi_212*xi_490;
      const double xi_500 = xi_499*0.11904761904761904;
      const double xi_501 = xi_498 + xi_500;
      const double xi_503 = xi_488*0.5;
      const double xi_504 = xi_502*xi_503;
      const double xi_518 = pow(xi_487*(1.0 - ((-xi_106*(xi_105*-2.0 + 4.0) + 1.0)*(-xi_106*(xi_105*-2.0 + 4.0) + 1.0))), 0.5);
      const double xi_519 = xi_310*xi_518;
      const double xi_527 = xi_223*xi_518;
      const double xi_528 = xi_527*0.16666666666666666;
      const double xi_529 = xi_526*xi_528;
      const double xi_530 = xi_525 + xi_529;
      const double xi_531 = xi_515*xi_516 + xi_517*xi_519 + xi_530;
      const double xi_533 = xi_218*xi_492;
      const double xi_534 = xi_533*0.10000000000000001;
      const double xi_541 = xi_488*xi_540;
      const double xi_542 = xi_332*xi_541;
      const double xi_545 = xi_225*xi_496;
      const double xi_546 = -xi_498 - xi_500;
      const double xi_551 = xi_545*0.071428571428571425;
      const double xi_552 = xi_331*xi_541;
      const double xi_561 = xi_528*xi_560;
      const double xi_562 = xi_559 + xi_561;
      const double xi_563 = xi_516*xi_555 + xi_519*xi_556 + xi_562;
      const double xi_565 = xi_545*0.042857142857142858;
      const double xi_573 = xi_528*xi_572;
      const double xi_574 = xi_571 + xi_573;
      const double xi_575 = xi_516*xi_568 + xi_519*xi_569 + xi_574;
      const double xi_576 = xi_488*xi_502*0.25;
      const double xi_579 = xi_499*0.083333333333333329;
      const double xi_583 = -xi_559 - xi_561;
      const double xi_584 = xi_503*(random_9_1_shifted - 0.5);
      const double xi_591 = xi_503*(random_10_1_shifted - 0.5);
      const double xi_597 = xi_527*0.083333333333333329;
      const double xi_598 = xi_526*xi_597;
      const double xi_599 = xi_419*xi_518;
      const double xi_600 = xi_517*xi_599;
      const double xi_604 = xi_572*xi_597;
      const double xi_605 = xi_569*xi_599;
      const double xi_606 = -xi_602 + xi_603 - xi_604 + xi_605;
      const double xi_608 = xi_545*0.014285714285714285;
      const double xi_610 = xi_499*0.023809523809523808;
      const double xi_611 = -xi_609 - xi_610;
      const double xi_613 = xi_602 - xi_603 + xi_604 - xi_605;
      const double xi_615 = -xi_576;
      const double xi_618 = xi_545*0.035714285714285712;
      const double xi_620 = xi_503*(random_10_0_shifted - 0.5);
      const double xi_625 = xi_556*xi_599;
      const double xi_626 = xi_560*xi_597;
      const double xi_627 = -xi_623 + xi_624 - xi_625 + xi_626;
      const double xi_629 = xi_623 - xi_624 + xi_625 - xi_626;
      const double u_0_b_collide = xib_18_collide;
      const double xi_113 = u_0_b_collide*xib_5_collide;
      const double xi_125 = xi_113*0.16666666666666666;
      const double xi_126 = xi_113*0.083333333333333329;
      const double xi_127 = rr_1_b_collide*xi_126;
      const double xi_128 = xi_125 - xi_127;
      const double xi_131 = -xi_125 + xi_127;
      const double xi_133 = rr_1_b_collide*xi_125;
      const double xi_153 = u_0_b_collide*xi_152;
      const double xi_158 = u_0_b_collide*xi_157;
      const double xi_181 = u_0_b_collide*xib_9_collide;
      const double xi_182 = xi_181*0.25;
      const double xi_185 = xi_156*xi_181;
      const double xi_473 = u_0_b_collide*u_0_b_collide;
      const double u_1_b_collide = xib_8_collide;
      const double xi_114 = u_1_b_collide*xib_29_collide;
      const double xi_117 = xi_114*0.16666666666666666;
      const double xi_118 = rr_1_b_collide*xi_117;
      const double xi_136 = xi_114*0.083333333333333329;
      const double xi_137 = rr_1_b_collide*xi_136;
      const double xi_138 = -xi_117 + xi_137;
      const double xi_140 = xi_117 - xi_137;
      const double xi_145 = xi_128 + xi_140;
      const double xi_146 = xi_131 + xi_138;
      const double xi_154 = u_1_b_collide*0.25;
      const double xi_155 = xi_154*xib_5_collide;
      const double xi_159 = u_1_b_collide*xi_156;
      const double xi_160 = xi_159*xib_5_collide;
      const double xi_161 = xi_153 + xi_155 - xi_158 - xi_160;
      const double xi_164 = -xi_153 - xi_155 + xi_158 + xi_160;
      const double xi_171 = xi_154*xib_9_collide;
      const double xi_173 = xi_159*xib_9_collide;
      const double xi_474 = rho_b_collide*(u_1_b_collide*u_1_b_collide);
      const double xi_507 = xi_474 + xi_497 + xi_506;
      const double xi_585 = rho_b_collide*u_1_b_collide;
      const double xi_587 = xi_586*(u_0_b_collide*xi_585 + xi_478 + xi_505);
      const double xi_588 = -xi_584 - xi_587;
      const double xi_589 = xi_584 + xi_587;
      const double u_2_b_collide = xib_13_collide;
      const double xi_115 = u_2_b_collide*xib_9_collide;
      const double xi_121 = xi_115*0.16666666666666666;
      const double xi_122 = xi_115*0.083333333333333329;
      const double xi_123 = rr_1_b_collide*xi_122;
      const double xi_124 = xi_121 - xi_123;
      const double xi_129 = xi_124 + xi_128;
      const double xi_130 = -xi_121 + xi_123;
      const double xi_132 = xi_130 + xi_131;
      const double xi_139 = xi_130 + xi_138;
      const double xi_141 = xi_124 + xi_140;
      const double xi_142 = rr_1_b_collide*xi_121;
      const double xi_151 = rr_1_b_collide*u_2_b_collide*xib_9_collide*-0.041666666666666664 + xi_122 + xi_146;
      const double xi_166 = xi_115*xi_165 - xi_122 + xi_145;
      const double xi_170 = rr_1_b_collide*u_0_b_collide*xib_5_collide*-0.041666666666666664 + xi_126 + xi_139;
      const double xi_172 = u_2_b_collide*xi_152;
      const double xi_174 = u_2_b_collide*xi_157;
      const double xi_175 = -xi_171 - xi_172 + xi_173 + xi_174;
      const double xi_178 = xi_113*xi_165 - xi_126 + xi_141;
      const double xi_180 = xi_114*xi_165 + xi_129 - xi_136;
      const double xi_183 = u_2_b_collide*xib_5_collide;
      const double xi_184 = xi_183*0.25;
      const double xi_186 = xi_156*xi_183;
      const double xi_187 = -xi_182 - xi_184 + xi_185 + xi_186;
      const double xi_189 = rr_1_b_collide*u_1_b_collide*xib_29_collide*-0.041666666666666664 + xi_132 + xi_136;
      const double xi_190 = xi_171 + xi_172 - xi_173 - xi_174;
      const double xi_191 = xi_182 + xi_184 - xi_185 - xi_186;
      const double xi_476 = rho_b_collide*(u_2_b_collide*u_2_b_collide);
      const double xi_477 = -xi_476 + xib_15_collide + xib_17_collide + xib_22_collide + xib_23_collide;
      const double xi_480 = rho_b_collide*xi_473 - xi_469 + xi_474 - xi_475 - xi_477 - xi_479 - xib_14_collide - xib_7_collide;
      const double xi_481 = rr_1_b_collide*xi_480;
      const double xi_508 = rr_1_b_collide*(xi_477 + xi_505 + xi_507 - xib_16_collide - xib_20_collide + xib_21_collide + xib_24_collide);
      const double xi_509 = xi_508*0.125;
      const double xi_510 = -xi_504 - xi_509;
      const double xi_537 = rho_b_collide*xi_473*2.0 - xi_476 - xi_479 - xi_507 - xi_535 - xi_536 + xib_10_collide*-2.0 + xib_11_collide*2.0 + xib_14_collide*2.0 - xib_15_collide - xib_17_collide + xib_19_collide*-2.0 - xib_22_collide - xib_23_collide + xib_4_collide*2.0 + xib_7_collide*2.0;
      const double xi_538 = rr_1_b_collide*xi_537;
      const double xi_539 = xi_538*0.041666666666666664;
      const double xi_543 = xi_539 + xi_542;
      const double xi_544 = xi_532 + xi_534 + xi_543;
      const double xi_547 = xi_504 + xi_509;
      const double xi_548 = -xi_539 - xi_542;
      const double xi_549 = -xi_532 - xi_534 + xi_548;
      const double xi_577 = xi_508*0.0625;
      const double xi_580 = xi_481*0.041666666666666664 + xi_494*0.083333333333333329;
      const double xi_581 = xi_538*0.020833333333333332 + xi_552*0.083333333333333329 + xi_580;
      const double xi_582 = xi_530 + xi_576 + xi_577 + xi_578 + xi_579 + xi_581;
      const double xi_590 = -xi_525 - xi_529 + xi_576 + xi_577 + xi_578 + xi_579 + xi_581;
      const double xi_594 = xi_586*(u_2_b_collide*xi_585 + xi_513 - xib_11_collide + xib_4_collide);
      const double xi_601 = xi_591 - xi_593 + xi_594 + xi_596 - xi_598 + xi_600;
      const double xi_612 = xi_468*0.025000000000000001 + xi_533*0.050000000000000003 + xi_548 + xi_580 - xi_607 - xi_608 + xi_611;
      const double xi_614 = rr_1_b_collide*xi_467*-0.025000000000000001 + rr_1_b_collide*xi_480*-0.041666666666666664 + xi_218*xi_488*xi_491*-0.050000000000000003 + xi_223*xi_488*xi_493*-0.083333333333333329 + xi_543 + xi_607 + xi_608 + xi_609 + xi_610;
      const double xi_616 = -xi_577;
      const double xi_619 = xi_574 + xi_581 + xi_611 + xi_615 + xi_616 + xi_617 + xi_618;
      const double xi_621 = xi_586*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_553 + xi_557);
      const double xi_622 = -xi_620 - xi_621;
      const double xi_628 = xi_620 + xi_621;
      const double xi_630 = xi_591 + xi_593 + xi_594 - xi_596 + xi_598 - xi_600;
      const double xi_631 = -xi_571 - xi_573 + xi_581 + xi_611 + xi_615 + xi_616 + xi_617 + xi_618;
      const double forceTerm_0_b_collide = xi_113*xi_116 - xi_113 + xi_114*xi_116 - xi_114 + xi_115*xi_116 - xi_115;
      const double forceTerm_1_b_collide = u_1_b_collide*xib_29_collide*0.33333333333333331 - xi_118 - xi_120 - xi_129;
      const double forceTerm_2_b_collide = xi_114*0.33333333333333331 - xi_118 + xi_120 + xi_132;
      const double forceTerm_3_b_collide = xi_113*0.33333333333333331 - xi_133 + xi_135 + xi_139;
      const double forceTerm_4_b_collide = u_0_b_collide*xib_5_collide*0.33333333333333331 - xi_133 - xi_135 - xi_141;
      const double forceTerm_5_b_collide = u_2_b_collide*xib_9_collide*0.33333333333333331 - xi_142 - xi_144 - xi_145;
      const double forceTerm_6_b_collide = xi_115*0.33333333333333331 - xi_142 + xi_144 + xi_146;
      const double forceTerm_7_b_collide = -xi_150 - xi_151 - xi_161;
      const double forceTerm_8_b_collide = -xi_151 - xi_163 - xi_164;
      const double forceTerm_9_b_collide = xi_161 + xi_163 + xi_166;
      const double forceTerm_10_b_collide = xi_150 + xi_164 + xi_166;
      const double forceTerm_11_b_collide = -xi_169 - xi_170 - xi_175;
      const double forceTerm_12_b_collide = xi_175 + xi_177 + xi_178;
      const double forceTerm_13_b_collide = xi_179 + xi_180 + xi_187;
      const double forceTerm_14_b_collide = -xi_187 - xi_188 - xi_189;
      const double forceTerm_15_b_collide = -xi_170 - xi_177 - xi_190;
      const double forceTerm_16_b_collide = xi_169 + xi_178 + xi_190;
      const double forceTerm_17_b_collide = xi_180 + xi_188 + xi_191;
      const double forceTerm_18_b_collide = -xi_179 - xi_189 - xi_191;
      const double tmp_a0 = forceTerm_0_a_collide + xi_193*0.10000000000000001 + xi_197*0.042857142857142858 + xi_206*-0.5 + xi_211*0.023809523809523808 + xi_213*xi_217 + xi_219*xi_221 - xi_224 + xi_226*xi_228 + xi_284 + xia_11_collide;
      const double tmp_a1 = forceTerm_1_a_collide + rr_1_a_collide*xi_196*0.014285714285714285 + xi_215*xi_225*xi_227*0.028571428571428571 - xi_285 - xi_288 - xi_292 - xi_301 - xi_322 - xi_337;
      const double tmp_a2 = forceTerm_2_a_collide + xi_197*0.014285714285714285 + xi_322 + xi_338*0.028571428571428571 + xi_339 + xi_340 + xi_341 + xi_343 + xia_17_collide;
      const double tmp_a3 = forceTerm_3_a_collide + xi_329*0.083333333333333329 + xi_340 - xi_344 - xi_345 + xi_346*0.33333333333333331 + xi_347 + xi_358 + xia_5_collide;
      const double tmp_a4 = forceTerm_4_a_collide + rr_1_a_collide*xi_328*0.083333333333333329 + xi_215*xi_331*xi_333*0.33333333333333331 - xi_292 - xi_344 - xi_345 - xi_348 - xi_358 - xi_359;
      const double tmp_a5 = forceTerm_5_a_collide + rr_1_a_collide*xi_210*0.015873015873015872 + xi_212*xi_215*xi_216*0.095238095238095233 - xi_326 - xi_337 - xi_341 - xi_360 - xi_361 - xi_362 - xi_372;
      const double tmp_a6 = forceTerm_6_a_collide + xi_211*0.015873015873015872 + xi_290*0.095238095238095233 + xi_301 + xi_343 - xi_361 - xi_362 + xi_372 + xi_373 + xia_24_collide;
      const double tmp_a7 = forceTerm_7_a_collide + xi_382 + xi_389 + xi_390 + xi_395 + xia_4_collide;
      const double tmp_a8 = forceTerm_8_a_collide + xi_357 + xi_389 + xi_401 + xi_402 + xia_22_collide;
      const double tmp_a9 = forceTerm_9_a_collide + xi_390 + xi_402 + xi_403 + xi_404 + xia_29_collide;
      const double tmp_a10 = forceTerm_10_a_collide + xi_357 + xi_395 + xi_404 + xi_405 + xia_6_collide;
      const double tmp_a11 = forceTerm_11_a_collide + xi_410 + xi_422 + xi_427 + xi_433 + xia_7_collide;
      const double tmp_a12 = forceTerm_12_a_collide - xi_422 - xi_437 - xi_438 - xi_439 + xia_15_collide;
      const double tmp_a13 = forceTerm_13_a_collide + xi_442 + xi_447 + xi_450 + xi_455 + xia_19_collide;
      const double tmp_a14 = forceTerm_14_a_collide + xi_447 + xi_458 + xi_459 + xi_460 + xia_21_collide;
      const double tmp_a15 = forceTerm_15_a_collide - xi_302 - xi_427 - xi_439 - xi_461 - xi_462;
      const double tmp_a16 = forceTerm_16_a_collide + xi_433 + xi_438 + xi_462 + xi_463 + xia_3_collide;
      const double tmp_a17 = forceTerm_17_a_collide + xi_455 + xi_459 + xi_464 + xi_465 + xia_28_collide;
      const double tmp_a18 = forceTerm_18_a_collide + xi_450 + xi_460 + xi_465 + xi_466 + xia_13_collide;
      const double tmp_b0 = forceTerm_0_b_collide + xi_213*xi_490 + xi_219*xi_492 + xi_226*xi_496 + xi_284 + xi_468*0.10000000000000001 + xi_472*0.042857142857142858 + xi_481*-0.5 + xi_486*0.023809523809523808 - xi_494 + xib_6_collide;
      const double tmp_b1 = forceTerm_1_b_collide + rr_1_b_collide*xi_471*0.014285714285714285 + xi_225*xi_488*xi_495*0.028571428571428571 - xi_288 - xi_497 - xi_501 - xi_510 - xi_531 - xi_544;
      const double tmp_b2 = forceTerm_2_b_collide + xi_339 + xi_472*0.014285714285714285 + xi_531 + xi_545*0.028571428571428571 + xi_546 + xi_547 + xi_549 + xib_3_collide;
      const double tmp_b3 = forceTerm_3_b_collide + xi_347 + xi_538*0.083333333333333329 + xi_546 - xi_550 - xi_551 + xi_552*0.33333333333333331 + xi_563 + xib_19_collide;
      const double tmp_b4 = forceTerm_4_b_collide + rr_1_b_collide*xi_537*0.083333333333333329 + xi_331*xi_488*xi_540*0.33333333333333331 - xi_359 - xi_501 - xi_550 - xi_551 - xi_563 + xib_10_collide;
      const double tmp_b5 = forceTerm_5_b_collide + rr_1_b_collide*xi_485*0.015873015873015872 + xi_212*xi_488*xi_489*0.095238095238095233 - xi_360 - xi_535 - xi_544 - xi_547 - xi_564 - xi_565 - xi_575;
      const double tmp_b6 = forceTerm_6_b_collide + xi_373 + xi_486*0.015873015873015872 + xi_499*0.095238095238095233 + xi_510 + xi_549 - xi_564 - xi_565 + xi_575 + xib_21_collide;
      const double tmp_b7 = forceTerm_7_b_collide + xi_382 + xi_582 + xi_583 + xi_588 + xib_20_collide;
      const double tmp_b8 = forceTerm_8_b_collide + xi_401 + xi_562 + xi_582 + xi_589 + xib_27_collide;
      const double tmp_b9 = forceTerm_9_b_collide + xi_403 + xi_583 + xi_589 + xi_590 + xib_26_collide;
      const double tmp_b10 = forceTerm_10_b_collide + xi_405 + xi_562 + xi_588 + xi_590 + xib_16_collide;
      const double tmp_b11 = forceTerm_11_b_collide + xi_410 + xi_601 + xi_606 + xi_612 + xib_11_collide;
      const double tmp_b12 = forceTerm_12_b_collide - xi_437 - xi_511 - xi_601 - xi_613 - xi_614;
      const double tmp_b13 = forceTerm_13_b_collide + xi_442 + xi_619 + xi_622 + xi_627 + xib_17_collide;
      const double tmp_b14 = forceTerm_14_b_collide + xi_458 + xi_619 + xi_628 + xi_629 + xib_15_collide;
      const double tmp_b15 = forceTerm_15_b_collide - xi_461 - xi_566 - xi_606 - xi_614 - xi_630;
      const double tmp_b16 = forceTerm_16_b_collide + xi_463 + xi_612 + xi_613 + xi_630 + xib_7_collide;
      const double tmp_b17 = forceTerm_17_b_collide + xi_464 + xi_627 + xi_628 + xi_631 + xib_23_collide;
      const double tmp_b18 = forceTerm_18_b_collide + xi_466 + xi_622 + xi_629 + xi_631 + xib_22_collide;
      const double xirecolor_0 = tmp_a0 + tmp_b0;
      const double xirecolor_1 = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const double xirecolor_2 = ((1.0) / (xirecolor_1));
      const double xi_653 = xirecolor_2*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const double xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const double xirecolor_4 = tmp_a1 + tmp_b1;
      const double xirecolor_5 = xi_239;
      const double xirecolor_6 = xi_237;
      const double xirecolor_7 = xi_257;
      const double xi_632 = -xirecolor_7;
      const double xirecolor_8 = xi_632;
      const double xirecolor_9 = xi_244;
      const double xirecolor_10 = -xirecolor_9;
      const double xirecolor_11 = xi_246;
      const double xi_633 = xirecolor_10 + xirecolor_11;
      const double xirecolor_12 = xi_259;
      const double xi_634 = xirecolor_12 + xirecolor_8;
      const double xirecolor_13 = xi_633 + xi_634;
      const double xirecolor_14 = xi_261;
      const double xirecolor_15 = xi_262;
      const double xi_635 = xirecolor_14 - xirecolor_15;
      const double xirecolor_16 = -xi_635;
      const double xirecolor_17 = xi_231;
      const double xirecolor_18 = xi_232;
      const double xirecolor_19 = xi_233;
      const double xi_636 = xirecolor_18 + xirecolor_19;
      const double xirecolor_20 = xi_234;
      const double xirecolor_21 = xi_636 - xirecolor_17 + xirecolor_20;
      const double xirecolor_22 = xi_230;
      const double xirecolor_23 = xi_229;
      const double xi_637 = xirecolor_22 - xirecolor_23;
      const double xirecolor_24 = -xi_637;
      const double xirecolor_25 = xi_265;
      const double xirecolor_26 = -xirecolor_25;
      const double xirecolor_27 = xi_251;
      const double xirecolor_28 = xi_250;
      const double xirecolor_29 = xi_267;
      const double xi_638 = xirecolor_26 + xirecolor_29;
      const double xirecolor_30 = xi_638 - xirecolor_27 + xirecolor_28;
      const double xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
      const double xirecolor_32 = xi_236;
      const double xirecolor_33 = xi_238;
      const double xi_639 = xirecolor_32 + xirecolor_33;
      const double xirecolor_34 = xi_639 + xirecolor_5 - xirecolor_6;
      const double xirecolor_35 = xi_241;
      const double xirecolor_36 = xi_242;
      const double xi_640 = xirecolor_35 - xirecolor_36;
      const double xirecolor_37 = -xi_640;
      const double xirecolor_38 = xi_633;
      const double xirecolor_39 = xi_248;
      const double xirecolor_40 = xi_249;
      const double xi_641 = xirecolor_39 + xirecolor_40;
      const double xirecolor_41 = xirecolor_27 - xirecolor_28;
      const double xirecolor_42 = xi_641 + xirecolor_41;
      const double xirecolor_43 = xirecolor_38 + xirecolor_42;
      const double xirecolor_44 = xi_637 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
      const double xirecolor_45 = xi_635;
      const double xirecolor_46 = xi_634;
      const double xirecolor_47 = xi_273;
      const double xirecolor_48 = -xirecolor_47;
      const double xirecolor_49 = xi_275;
      const double xi_642 = xirecolor_48 + xirecolor_49;
      const double xirecolor_50 = xi_640 + xi_642;
      const double xirecolor_51 = xi_636 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
      const double xirecolor_52 = pow(xirecolor_31*xirecolor_31 + xirecolor_44*xirecolor_44 + xirecolor_51*xirecolor_51, 0.5);
      const double xirecolor_53 = ((1.0) / (xirecolor_52));
      const double xi_643 = xirecolor_31*xirecolor_53;
      const double xi_644 = xirecolor_51*xirecolor_53;
      const double xi_645 = xirecolor_44*xirecolor_53;
      const bool xirecolor_54 = xirecolor_52 > 9.9999999999999998e-13;
      const double xirecolor_55 = beta*((1.0) / (xirecolor_1*xirecolor_1))*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2]*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const double xirecolor_56 = xirecolor_55*(0.055555555555555552*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.055555555555555552*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const double xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_643): (0.0));
      const double xirecolor_58 = tmp_a2 + tmp_b2;
      const double xirecolor_59 = xirecolor_56*((xirecolor_54) ? (-xi_643): (0.0));
      const double xirecolor_60 = tmp_a3 + tmp_b3;
      const double xirecolor_61 = xirecolor_56*((xirecolor_54) ? (-xi_644): (0.0));
      const double xirecolor_62 = tmp_a4 + tmp_b4;
      const double xirecolor_63 = xirecolor_56*((xirecolor_54) ? (xi_644): (0.0));
      const double xirecolor_64 = tmp_a5 + tmp_b5;
      const double xirecolor_65 = xirecolor_56*((xirecolor_54) ? (-xi_645): (0.0));
      const double xirecolor_66 = tmp_a6 + tmp_b6;
      const double xirecolor_67 = xirecolor_56*((xirecolor_54) ? (xi_645): (0.0));
      const double xirecolor_68 = tmp_a7 + tmp_b7;
      const double xirecolor_69 = xi_374;
      const double xirecolor_70 = xi_375;
      const double xirecolor_71 = xi_639;
      const double xirecolor_72 = xirecolor_50 + xirecolor_71;
      const double xirecolor_73 = xi_379;
      const double xirecolor_74 = xirecolor_25 - xirecolor_29 + xirecolor_73;
      const double xirecolor_75 = xi_376 - xirecolor_11 + xirecolor_41 - xirecolor_69 + xirecolor_70 + xirecolor_72 + xirecolor_74 + xirecolor_9;
      const double xirecolor_76 = xi_310*xirecolor_53;
      const double xi_646 = xirecolor_75*xirecolor_76;
      const double xirecolor_77 = xirecolor_55*(0.027777777777777776*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.027777777777777776*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const double xirecolor_78 = xirecolor_77*((xirecolor_54) ? (-xi_646): (0.0));
      const double xirecolor_79 = tmp_a8 + tmp_b8;
      const double xirecolor_80 = xi_396;
      const double xirecolor_81 = xi_397;
      const double xi_647 = -xirecolor_80 + xirecolor_81;
      const double xirecolor_82 = xi_398;
      const double xirecolor_83 = xi_399 + xi_647 + xirecolor_30 + xirecolor_38 + xirecolor_72 + xirecolor_82;
      const double xi_648 = xirecolor_76*xirecolor_83;
      const double xirecolor_84 = xirecolor_77*((xirecolor_54) ? (xi_648): (0.0));
      const double xirecolor_85 = tmp_a9 + tmp_b9;
      const double xirecolor_86 = xirecolor_77*((xirecolor_54) ? (-xi_648): (0.0));
      const double xirecolor_87 = tmp_a10 + tmp_b10;
      const double xirecolor_88 = xirecolor_77*((xirecolor_54) ? (xi_646): (0.0));
      const double xirecolor_89 = tmp_a11 + tmp_b11;
      const double xirecolor_90 = xi_641 + xirecolor_37 + xirecolor_71;
      const double xirecolor_91 = -xi_632 - xi_647 - xirecolor_12;
      const double xirecolor_92 = xi_406 + xirecolor_45 + xirecolor_74 + xirecolor_90 + xirecolor_91;
      const double xi_649 = xirecolor_76*xirecolor_92;
      const double xirecolor_93 = xirecolor_77*((xirecolor_54) ? (-xi_649): (0.0));
      const double xirecolor_94 = tmp_a12 + tmp_b12;
      const double xirecolor_95 = xirecolor_16 + xirecolor_69 - xirecolor_70;
      const double xirecolor_96 = xi_434 + xi_638 + xirecolor_46 + xirecolor_82 + xirecolor_90 + xirecolor_95;
      const double xi_650 = xirecolor_76*xirecolor_96;
      const double xirecolor_97 = xirecolor_77*((xirecolor_54) ? (-xi_650): (0.0));
      const double xirecolor_98 = tmp_a13 + tmp_b13;
      const double xirecolor_99 = xi_440 + xi_642 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_73 + xirecolor_82;
      const double xi_651 = xirecolor_76*xirecolor_99;
      const double xirecolor_100 = xirecolor_77*((xirecolor_54) ? (-xi_651): (0.0));
      const double xirecolor_101 = tmp_a14 + tmp_b14;
      const double xirecolor_102 = xi_456 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_91 + xirecolor_95;
      const double xi_652 = xirecolor_102*xirecolor_76;
      const double xirecolor_103 = xirecolor_77*((xirecolor_54) ? (-xi_652): (0.0));
      const double xirecolor_104 = tmp_a15 + tmp_b15;
      const double xirecolor_105 = xirecolor_77*((xirecolor_54) ? (xi_650): (0.0));
      const double xirecolor_106 = tmp_a16 + tmp_b16;
      const double xirecolor_107 = xirecolor_77*((xirecolor_54) ? (xi_649): (0.0));
      const double xirecolor_108 = tmp_a17 + tmp_b17;
      const double xirecolor_109 = xirecolor_77*((xirecolor_54) ? (xi_652): (0.0));
      const double xirecolor_110 = tmp_a18 + tmp_b18;
      const double xirecolor_111 = xirecolor_77*((xirecolor_54) ? (xi_651): (0.0));
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
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2] = xi_653*xirecolor_0;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3] = xi_653*xirecolor_4 - xirecolor_57;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3] = xi_653*xirecolor_58 - xirecolor_59;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3] = xi_653*xirecolor_60 - xirecolor_61;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3] = xi_653*xirecolor_62 - xirecolor_63;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3] = xi_653*xirecolor_64 - xirecolor_65;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3] = xi_653*xirecolor_66 - xirecolor_67;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3] = xi_653*xirecolor_68 - xirecolor_78;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3] = xi_653*xirecolor_79 - xirecolor_84;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3] = xi_653*xirecolor_85 - xirecolor_86;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3] = xi_653*xirecolor_87 - xirecolor_88;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3] = xi_653*xirecolor_89 - xirecolor_93;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3] = xi_653*xirecolor_94 - xirecolor_97;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3] = xi_653*xirecolor_98 - xirecolor_100;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3] = xi_653*xirecolor_101 - xirecolor_103;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3] = xi_653*xirecolor_104 - xirecolor_105;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3] = xi_653*xirecolor_106 - xirecolor_107;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3] = xi_653*xirecolor_108 - xirecolor_109;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3] = xi_653*xirecolor_110 - xirecolor_111;
   } 
}
}


void ColorGradientCollideSweepDoublePrecisionCUDA::run(IBlock * block, gpuStream_t stream)
{
   if (!this->configured_)
      WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")
         
    auto rho_a = block->getData< gpu::GPUField<double> >(rho_aID);
    auto pdfs_b = block->getData< gpu::GPUField<double> >(pdfs_bID);
    auto force_b = block->getData< gpu::GPUField<double> >(force_bID);
    auto pdfs_a = block->getData< gpu::GPUField<double> >(pdfs_aID);
    auto phasefield = block->getData< gpu::GPUField<double> >(phasefieldID);
    auto rho_b = block->getData< gpu::GPUField<double> >(rho_bID);
    auto velocity = block->getData< gpu::GPUField<double> >(velocityID);
    auto force_a = block->getData< gpu::GPUField<double> >(force_aID);

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
    double * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
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
    internal_colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda::colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, block_offset_0, block_offset_1, block_offset_2, kT, omega_shear_a, omega_shear_b, seed, sigma, time_step);
    
}


void ColorGradientCollideSweepDoublePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block, gpuStream_t stream)
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

    auto rho_a = block->getData< gpu::GPUField<double> >(rho_aID);
    auto pdfs_b = block->getData< gpu::GPUField<double> >(pdfs_bID);
    auto force_b = block->getData< gpu::GPUField<double> >(force_bID);
    auto pdfs_a = block->getData< gpu::GPUField<double> >(pdfs_aID);
    auto phasefield = block->getData< gpu::GPUField<double> >(phasefieldID);
    auto rho_b = block->getData< gpu::GPUField<double> >(rho_bID);
    auto velocity = block->getData< gpu::GPUField<double> >(velocityID);
    auto force_a = block->getData< gpu::GPUField<double> >(force_aID);

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
    double * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
    internal_colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda::colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, block_offset_0, block_offset_1, block_offset_2, kT, omega_shear_a, omega_shear_b, seed, sigma, time_step);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
