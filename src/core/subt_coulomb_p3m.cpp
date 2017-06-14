/*
   Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
   Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
   Max-Planck-Institute for Polymer Research, Theory Group

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
/** \file subt_coulomb_p3m.cpp
 *
 *  Implementation of \ref subt_coublom_p3m.hpp
 */
#include "subt_coulomb_p3m.hpp"

#ifdef P3M 
#include <cmath>
#include "utils.hpp"
#include "communication.hpp"
#include "p3m-common.hpp"
#include "grid.hpp"

int subt_coulomb_p3m_set_params(int bond_type, int res, double r_max)
{
    if (bond_type < 0)
        return ES_ERROR;

    if (p3m.params.alpha == 0)
    {
        fprintf(stderr, "ERROR: P3M alpha is not set. Long-Range contribution cannot be precomputed. Tune P3M before creating this bond!\n");
        return ES_ERROR;
    }

    make_bond_type_exist(bond_type);

    bonded_ia_params[bond_type].type = BONDED_IA_SUBT_COULOMB_P3M;  
    bonded_ia_params[bond_type].num = 1;
    
    /*Precompute long-range forces

    Issues:
    -Box len changes
    -Alpha changes
    -Proper resolution ? User ?
    -Bond deleted ?
    
    */ 

    bonded_ia_params[bond_type].p.subt_coulomb_p3m.res = res;
    bonded_ia_params[bond_type].p.subt_coulomb_p3m.r_max = r_max;
    
    bonded_ia_params[bond_type].p.subt_coulomb_p3m.long_range_forces = (double*)Utils::malloc(res*sizeof(double));

    double pre1 = 1.0/(4.0*p3m.params.alpha*p3m.params.alpha);

    int sn = 10;
    double kvec[3];
    double ksq;

    for (int di = 0; di < res; ++di) {
        double d = r_max*di/res;
        double F = 0;
        for (int m=-sn; m<sn; ++m) {
            for (int n=-sn; n<sn; ++n) {
                for (int o=-sn; o<sn; ++o) {
                    if (m==0 && n==0 && o==0) 
                        continue;
                    kvec[0] = 2.0 * PI/box_l[0] * m;
                    kvec[1] = 2.0 * PI/box_l[1] * n;
                    kvec[2] = 2.0 * PI/box_l[2] * o;
                    ksq = kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];

                    F += exp(-ksq*pre1)*4.0*PI/ksq*kvec[0]*sin(kvec[0]*d);
                }
            }
        }
        //for(i=0;i<3;i++)
        //    F[i] *= q1*q2/V
        bonded_ia_params[bond_type].p.subt_coulomb_p3m.long_range_forces[di] = F;
    }


    mpi_bcast_ia_params(bond_type, -1); 

    return ES_OK;
}

#endif

