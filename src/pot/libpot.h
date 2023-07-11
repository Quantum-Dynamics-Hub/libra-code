/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef LIBPOT_H
#define LIBPOT_H


#include "Switching_functions.h"
#include "Potentials_bonds.h"
#include "Potentials_angles.h"
#include "Potentials_stretch_bend.h"
#include "Potentials_dihedrals.h"
#include "Potentials_oop.h"
#include "Potentials_vdw.h"
#include "Potentials_elec.h"
#include "Potentials_frag.h"

#include "Potentials_mb_vdw.h"
#include "Potentials_mb_elec.h"

/// liblibra namespace
namespace liblibra{

namespace libpot{

void export_Pot_objects();

}// namespace libpot

}// liblibra

#endif// LIBPOT_H
