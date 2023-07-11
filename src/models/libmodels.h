/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libmodels.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_MODELS_H
#define LIB_MODELS_H


#include "../math_linalg/liblinalg.h"
#include "../math_meigen/libmeigen.h"

#include "Models_1_state.h"
#include "Models_2_state.h"
#include "Model_cubic.h"
#include "Model_DAC.h"
#include "Model_double_well.h"
#include "Model_ECWR.h"
#include "Model_Marcus.h"
#include "Model_Rabi2.h"
#include "Model_SAC.h"
#include "Model_SEXCH.h"
#include "Model_sin.h"
#include "Model_sin_2D.h"


/// liblibra namespace
namespace liblibra{

/// libmodels
namespace libmodels{

  void export_models_objects();


}// namespace libmodels
}// liblibra

#endif// LIB_MODELS_H
