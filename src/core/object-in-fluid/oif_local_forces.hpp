/*
  Copyright (C) 2012-2018 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _OBJECT_IN_FLUID_OIF_LOCAL_FORCES_H
#define _OBJECT_IN_FLUID_OIF_LOCAL_FORCES_H

/** \file
 *  Routines to calculate the OIF_LOCAL_FORCES
 *  for a particle quadruple (two neighboring triangles with common edge).
 * (Dupin2007) \ref forces.cpp
 */

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include <utils/math/triangle_functions.hpp>

inline int cimo_vector_product(double &res0, double &res1, double &res2, double a0, double a1, double a2, double b0, double b1, double b2) {
    res0 = a1 * b2 - a2 * b1;
    res1 = a2 * b0 - a0 * b2;
    res2 = a0 * b1 - a1 * b0;
    return 1;
}

// set parameters for local forces
int oif_local_forces_set_params(int bond_type, double r0, double ks,
                                double kslin, double phi0, double kb,
                                double A01, double A02, double kal,
                                double kvisc);

inline double KS(double lambda) { // Defined by (19) from Dupin2007
  double res;
  res = (pow(lambda, 0.5) + pow(lambda, -2.5)) / (lambda + pow(lambda, -3.));
  return res;
}

/** Computes the local forces (Dupin2007) and adds them
 *  to the particle forces.
 *  @param p1           %Particle of triangle 1.
 *  @param p2 , p3      Particles common to triangle 1 and triangle 2.
 *  @param p4           %Particle of triangle 2.
 *  @param iaparams     Bonded parameters for the OIF interaction.
 *  @param force        Force on @p p1.
 *  @param force2       Force on @p p2.
 *  @param force3       Force on @p p3.
 *  @param force4       Force on @p p4.
 *  @return 0
 */
inline int calc_oif_local(Particle *p2, Particle *p1, Particle *p3,
                          Particle *p4, Bonded_ia_parameters *iaparams,
                          double force[3], double force2[3], double force3[3],
                          double force4[3]) // first-fold-then-the-same approach
{

  auto const fp2 = unfolded_position(*p2);
  auto const fp1 = fp2 + get_mi_vector(p1->r.p, fp2);
  auto const fp3 = fp2 + get_mi_vector(p3->r.p, fp2);
  auto const fp4 = fp2 + get_mi_vector(p4->r.p, fp2);

  for (int i = 0; i < 3; i++) {
    force[i] = 0;
    force2[i] = 0;
    force3[i] = 0;
    force4[i] = 0;
  }

  // non-linear stretching
  if (iaparams->p.oif_local_forces.ks > TINY_OIF_ELASTICITY_COEFFICIENT) {
    auto const dx = fp2 - fp3;
    auto const len = dx.norm();
    auto const dr = len - iaparams->p.oif_local_forces.r0;
    auto const lambda = 1.0 * len / iaparams->p.oif_local_forces.r0;
    auto const fac =
        -iaparams->p.oif_local_forces.ks * KS(lambda) * dr; // no normalization
    for (int i = 0; i < 3; i++) {
      force2[i] += fac * dx[i] / len;
      force3[i] += -fac * dx[i] / len;
    }
  }

  // linear stretching
  if (iaparams->p.oif_local_forces.kslin > TINY_OIF_ELASTICITY_COEFFICIENT) {
    auto const dx = fp2 - fp3;
    auto const len = dx.norm();
    auto const dr = len - iaparams->p.oif_local_forces.r0;
    auto const fac =
        -iaparams->p.oif_local_forces.kslin * dr; // no normalization

    for (int i = 0; i < 3; i++) {
      force2[i] += fac * dx[i] / len;
      force3[i] += -fac * dx[i] / len;
    }
  }

  // viscous force
  if (iaparams->p.oif_local_forces.kvisc >
      TINY_OIF_ELASTICITY_COEFFICIENT) { // to be implemented....
    auto const dx = fp2 - fp3;
    auto const len2 = dx.norm2();
    auto const v_ij = p3->m.v - p2->m.v;

    // Variant A
    // Here the force is in the direction of relative velocity btw points

    // Code:
    // for(int i=0;i<3;i++) {
    // force2[i] += iaparams->p.oif_local_forces.kvisc*v[i];
    // force3[i] -= iaparams->p.oif_local_forces.kvisc*v[i];
    //}

    // Variant B
    // Here the force is the projection of relative velocity btw points onto
    // line btw the points

    // denote p vector between p2 and p3
    // denote v the velocity difference between the points p2 and p3
    // denote alpha the angle between p and v
    // denote x the projevted v onto p
    // cos alpha = |x|/|v|
    // cos alpha = (v,p)/(|v||p|)
    // together we get |x|=(v,p)/|p|
    // also, x is along p, so x = |x|.p/|p|
    // so x = p/|p| . (v,p)/|p|
    // altogether x = p . (v,p)/(|p|)^2
    // |p|^2 is stored in len2

    // Code:
    auto const fac = iaparams->p.oif_local_forces.kvisc * (dx * v_ij) / len2;

    for (int i = 0; i < 3; i++) {
      force2[i] += fac * dx[i];
      force3[i] -= fac * dx[i];
    }
  }

  /* bending
     implemented according to updated expressions in book Computational Blood Cell Mechanics, by I.Cimrak and I.Jancigova, see errata...*/
  if (iaparams->p.oif_local_forces.kb > TINY_OIF_ELASTICITY_COEFFICIENT) {
    // how fp1 - fp4 correspond to points A,B,C,D from the book, Figure A.1:
     //    fp1 -> C
     //    fp2 -> A
     //    fp3 -> B
     //    fp4 -> D
    
    auto const Nc = Utils::get_n_triangle(fp1, fp2, fp3);    // returns (fp2 - fp1)x(fp3 - fp1), thus Nc = (A - C)x(B - C)  
    auto const Nd = Utils::get_n_triangle(fp4, fp3, fp2);    // returns (fp3 - fp4)x(fp2 - fp4), thus Nd = (B - D)x(A - D)  

    auto const phi = Utils::angle_btw_triangles(fp1, fp2, fp3, fp4);
    auto const aa = (phi - iaparams->p.oif_local_forces
                               .phi0); // no renormalization by phi0, to be
                                       // consistent with Krueger and Fedosov
    auto const BminA = fp3 - fp2;
    auto const fac = iaparams->p.oif_local_forces.kb * aa;
    auto const factorFaNc = (fp2 - fp3) * (fp1 - fp3) / BminA.norm() / Nc.norm2();
    auto const factorFaNd = (fp2 - fp3) * (fp4 - fp3) / BminA.norm() / Nd.norm2();
    auto const factorFbNc = (fp2 - fp3) * (fp2 - fp1) / BminA.norm() / Nc.norm2();
    auto const factorFbNd = (fp2 - fp3) * (fp2 - fp4) / BminA.norm() / Nd.norm2();
    
    for (int i = 0; i < 3; i++) {
      force[i] -= fac * BminA.norm()/Nc.norm2() * Nc[i];                // Fc
      force2[i] +=  fac * (factorFaNc * Nc[i] + factorFaNd * Nd[i]);    // Fa
      force3[i] +=  fac * (factorFbNc * Nc[i] + factorFbNd * Nd[i]);    // Fb
      force4[i] -= fac * BminA.norm()/Nd.norm2() * Nd[i];               // Fd
    }
    
    
    double newforce[3],newforce2[3],newforce3[3],newforce4[3];
    for (int i = 0; i < 3; i++) {
      newforce[i] = - fac * BminA.norm()/Nc.norm2() * Nc[i];                // Fc
      newforce2[i] =  fac * (factorFaNc * Nc[i] + factorFaNd * Nd[i]);   // Fa
      newforce3[i] =  fac * (factorFbNc * Nc[i] + factorFbNd * Nd[i]);   // Fb
      newforce4[i] = - fac * BminA.norm()/Nd.norm2() * Nd[i];               // Fd
    }
    
    
    // comparison with force- and torque-free implementation from Krueger, see book Computational Blood Cell mechanics, appendix A.2
    auto const nc = Utils::get_n_triangle(fp1, fp2, fp3).normalize();
    auto const nd = Utils::get_n_triangle(fp4, fp3, fp2).normalize(); 
    
    auto const denominator = sqrt(1 - (nc * nd) * (nc * nd));
    double inv_twiceSABC = 1.0 / Utils::get_n_triangle(fp1, fp2, fp3).norm();
    double inv_twiceSABD = 1.0 / Utils::get_n_triangle(fp4, fp3, fp2).norm();
    double scal = nc * nd;
    auto const nd_min_scal_nc = nd - scal * nc;
    auto const nc_min_scal_nd = nc - scal * nd;
    
    double cross1Fa0,cross1Fa1,cross1Fa2;
    cimo_vector_product(cross1Fa0,cross1Fa1,cross1Fa2,(fp1 - fp3)[0],(fp1 - fp3)[1],(fp1 - fp3)[2],nd_min_scal_nc[0],nd_min_scal_nc[1],nd_min_scal_nc[2]);
    double cross2Fa0,cross2Fa1,cross2Fa2;
    cimo_vector_product(cross2Fa0,cross2Fa1,cross2Fa2,(fp3 - fp4)[0],(fp3 - fp4)[1],(fp3 - fp4)[2],nc_min_scal_nd[0],nc_min_scal_nd[1],nc_min_scal_nd[2]);

    double cross1Fb0,cross1Fb1,cross1Fb2;
    cimo_vector_product(cross1Fb0,cross1Fb1,cross1Fb2,(fp2 - fp1)[0],(fp2 - fp1)[1],(fp2 - fp1)[2],nd_min_scal_nc[0],nd_min_scal_nc[1],nd_min_scal_nc[2]);
    double cross2Fb0,cross2Fb1,cross2Fb2;
    cimo_vector_product(cross2Fb0,cross2Fb1,cross2Fb2,(fp4 - fp2)[0],(fp4 - fp2)[1],(fp4 - fp2)[2],nc_min_scal_nd[0],nc_min_scal_nd[1],nc_min_scal_nd[2]);

    double crossFc0,crossFc1,crossFc2;
    cimo_vector_product(crossFc0,crossFc1,crossFc2,(fp3 - fp2)[0],(fp3 - fp2)[1],(fp3 - fp2)[2],nd_min_scal_nc[0],nd_min_scal_nc[1],nd_min_scal_nc[2]);

    double crossFd0,crossFd1,crossFd2;
    cimo_vector_product(crossFd0,crossFd1,crossFd2,(fp2 - fp3)[0],(fp2 - fp3)[1],(fp2 - fp3)[2],nc_min_scal_nd[0],nc_min_scal_nd[1],nc_min_scal_nd[2]);

    // Krueger
    double FA[3],FB[3],FC[3],FD[3];
    FA[0] = fac * 1.0 / denominator * (inv_twiceSABC*cross1Fa0 + inv_twiceSABD*cross2Fa0);    
    FA[1] = fac * 1.0 / denominator * (inv_twiceSABC*cross1Fa1 + inv_twiceSABD*cross2Fa1);    
    FA[2] = fac * 1.0 / denominator * (inv_twiceSABC*cross1Fa2 + inv_twiceSABD*cross2Fa2);    

    FB[0] = fac * 1.0 / denominator * (inv_twiceSABC*cross1Fb0 + inv_twiceSABD*cross2Fb0);    
    FB[1] = fac * 1.0 / denominator * (inv_twiceSABC*cross1Fb1 + inv_twiceSABD*cross2Fb1);    
    FB[2] = fac * 1.0 / denominator * (inv_twiceSABC*cross1Fb2 + inv_twiceSABD*cross2Fb2);    

    FC[0] = fac * 1.0 / denominator * inv_twiceSABC*crossFc0;
    FC[1] = fac * 1.0 / denominator * inv_twiceSABC*crossFc1;    
    FC[2] = fac * 1.0 / denominator * inv_twiceSABC*crossFc2;    

    FD[0] = fac * 1.0 / denominator * inv_twiceSABD*crossFd0;
    FD[1] = fac * 1.0 / denominator * inv_twiceSABD*crossFd1;    
    FD[2] = fac * 1.0 / denominator * inv_twiceSABD*crossFd2;    

    //for (int i = 0; i < 3; i++) {
      //force[i] += FC[i];                // Fc
      //force2[i] +=  FA[i];    // Fa
      //force3[i] +=  FB[i];    // Fb
      //force4[i] += FD[i];               // Fc
    //}
 



    printf("%2.10lf %2.10lf %2.10lf \n",FA[0] - newforce2[0], FA[1] - newforce2[1], FA[2] - newforce2[2]);
    printf("%2.10lf %2.10lf %2.10lf \n",FB[0] - newforce3[0], FB[1] - newforce3[1], FB[2] - newforce3[2]);
    printf("%2.10lf %2.10lf %2.10lf \n",FC[0] - newforce[0], FC[1] - newforce[1], FC[2] - newforce[2]);
    printf("%2.10lf %2.10lf %2.10lf \n\n",FD[0] - newforce4[0], FD[1] - newforce4[1], FD[2] - newforce4[2]);
    printf("%2.10lf %2.10lf %2.10lf %2.10lf %2.10lf %2.10lf \n",FA[0], FA[1], FA[2], newforce2[0],  newforce2[1], newforce2[2]);


    auto const tt = 0.25*(fp1 + fp2 + fp3 + fp4); //tt is centroid
    double torqueA0,torqueA1,torqueA2;
    double torqueB0,torqueB1,torqueB2;
    double torqueC0,torqueC1,torqueC2;
    double torqueD0,torqueD1,torqueD2;
    cimo_vector_product(torqueC0,torqueC1,torqueC2,(fp1 - tt)[0],(fp1 - tt)[1],(fp1 - tt)[2],FC[0],FC[1],FC[2]);
    cimo_vector_product(torqueA0,torqueA1,torqueA2,(fp2 - tt)[0],(fp2 - tt)[1],(fp2 - tt)[2],FA[0],FA[1],FA[2]);
    cimo_vector_product(torqueB0,torqueB1,torqueB2,(fp3 - tt)[0],(fp3 - tt)[1],(fp3 - tt)[2],FB[0],FB[1],FB[2]);
    cimo_vector_product(torqueD0,torqueD1,torqueD2,(fp4 - tt)[0],(fp4 - tt)[1],(fp4 - tt)[2],FD[0],FD[1],FD[2]);
    double total_torque0,total_torque1,total_torque2;
    total_torque0 = torqueA0 + torqueB0 + torqueC0 + torqueD0;
    total_torque1 = torqueA1 + torqueB1 + torqueC1 + torqueD1;
    total_torque2 = torqueA2 + torqueB2 + torqueC2 + torqueD2;
    printf("    torque free old implementation (should be non-torque-free): %lf %lf %lf\n\n", total_torque0, total_torque1, total_torque2);






  }

  /* local area
     for both triangles
     only 1/3 of calculated forces are added, because each triangle will enter
     this calculation 3 times (one time per edge)

              Proportional distribution of forces, implemented according to the
     article I.Jancigova, I.Cimrak, Non-uniform force allocation for area
     preservation in spring network models, International Journal for Numerical
     Methods in Biomedical Engineering, DOI: 10.1002/cnm.2757

  */
  if (iaparams->p.oif_local_forces.kal > TINY_OIF_ELASTICITY_COEFFICIENT) {

    auto handle_triangle = [](double kal, double A0, Utils::Vector3d const &fp1,
                              Utils::Vector3d const &fp2,
                              Utils::Vector3d const &fp3, double force1[3],
                              double force2[3], double force3[3]) {
      auto const h = (1. / 3.) * (fp1 + fp2 + fp3);
      auto const A = Utils::area_triangle(fp1, fp2, fp3);
      auto const t = sqrt(A / A0) - 1.0;

      auto const m1 = h - fp1;
      auto const m2 = h - fp2;
      auto const m3 = h - fp3;

      auto const m1_length = m1.norm();
      auto const m2_length = m2.norm();
      auto const m3_length = m3.norm();

      auto const fac = kal * A0 * (2 * t + t * t) /
                       (m1_length * m1_length + m2_length * m2_length +
                        m3_length * m3_length);

      for (int i = 0; i < 3; i++) { // local area force for p1
        force1[i] += fac * m1[i] / 3.0;
        force2[i] += fac * m2[i] / 3.0;
        force3[i] += fac * m3[i] / 3.0;
      }
    };

    handle_triangle(iaparams->p.oif_local_forces.kal,
                    iaparams->p.oif_local_forces.A01, fp1, fp2, fp3, force,
                    force2, force3);
    handle_triangle(iaparams->p.oif_local_forces.kal,
                    iaparams->p.oif_local_forces.A02, fp2, fp3, fp4, force2,
                    force3, force4);
  }
  return 0;
}

#endif
