/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef IVR_H
#define IVR_H

#include "../math_linalg/liblinalg.h"
#include "../math_random/librandom.h"
#include "../math_meigen/libmeigen.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace librandom;


namespace libivr{


MATRIX ivr_Husimi(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd);
vector<MATRIX> ivr_Husimi(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size);

MATRIX ivr_LSC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd);
vector<MATRIX> ivr_LSC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size);

MATRIX ivr_DHK(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd);
vector<MATRIX> ivr_DHK(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size);

MATRIX ivr_FB_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, int flag, Random& rnd);
vector<MATRIX> ivr_FB_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, int flag, Random& rnd, int sample_size);

MATRIX ivr_FF_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, Random& rnd);
vector<MATRIX> ivr_FF_MQC(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, Random& rnd, int sample_size);

}// namespace libivr
}// liblibra

#endif // IVR_H
