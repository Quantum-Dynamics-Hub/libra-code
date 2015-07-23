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


namespace libpot{

void export_Pot_objects();

}// namespace libpot

#endif// LIBPOT_H
