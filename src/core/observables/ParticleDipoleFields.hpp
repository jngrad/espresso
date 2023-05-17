/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef OBSERVABLES_PARTICLEDIPOLEFIELDS_HPP
#define OBSERVABLES_PARTICLEDIPOLEFIELDS_HPP

#include "PidObservable.hpp"

namespace Observables {

/** Extract particle directors.
 *  For \f$n\f$ particles, return \f$3 n\f$ directors ordered as
 *  \f$(d1_x, d1_y, d1_z, \dots, dn_x, dn_y, dn_z)\f$.
 */
using ParticleDipoleFields = ParticleObservable<ParticleObservables::DipoleFields>;

} // namespace Observables
#endif
