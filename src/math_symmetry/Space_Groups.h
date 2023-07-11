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
  \file Space_Groups.h
  \brief The file describes and implements the SPACE_GROUP class and describes related functions
    
*/

#ifndef SPACE_GROUPS_H
#define SPACE_GROUPS_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <string>
#include <vector>
#endif

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libsymmetry namespace
namespace libsymmetry{


//----------- Some functions ----------------
void Apply_Symmetry(std::string space_group_name,VECTOR r,std::vector<VECTOR>& r_equiv);

//-------------------------------------------


class SPACE_GROUP{
/**
  This class holds and generates the set of symmetry operations that correspond to a specific
  space group of symmetry

  Internet source: http://img.chem.ucl.ac.uk/sgp/large/sgp.htm
*/

public:
      //------- Data ------------
      vector<MATRIX> operators; ///< The list of matrices representing the symmetry operations. These matrices are 3x4 in size
                                ///< The first 3 coloumns form 3x3 matrix of rotations/reflections operations. The last coloumn
                                ///< gives an additional shift operation 


      //------- Methods ---------
      SPACE_GROUP(){;;}  ///< Default constructor
      SPACE_GROUP(std::string space_group_name){
      /**
        This constructor creates all the symmetry operation for given space group name

        \param[in] space_group_name Self-explanatory. Examples: P_1, P_-1, C1_2_1, I_a_-3_d, etc.
      */
 
      int sz,i;
      MATRIX m(3,4);  // This is actually a temporary object = R(3,3) + T(3,1) 
      //================== TRICLINIC ========================
      if(space_group_name == "P_1"){
         m = 0.0;
         m.M[0] =  1.0;  
                         m.M[5] = 1.0; 
                                        m.M[10] = 1.0;     operators.push_back(m);  // x,  y,  z      

      }// P1
      else if(space_group_name == "P_-1"){
         m = 0.0;
         m.M[0] =  1.0; 
                         m.M[5] =  1.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;  
                         m.M[5] = -1.0; 
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -y,  -z

      }

      //================== MONOCLINIC =======================
      else if(space_group_name == "P1_2_1"){
         m = 0.0;
         m.M[0] =  1.0;  
                         m.M[5] =  1.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;  
                         m.M[5] =  1.0; 
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y,  -z

      }
      else if(space_group_name == "P1_21_1"){
         m = 0.0;
         m.M[0] =  1.0;  
                         m.M[5] =  1.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;  
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y+0.5,  -z

      }
      else if(space_group_name == "C1_2_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y,  -z
         // + (0.5, 0.5, 0.0)
         sz = 2;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "P1_m_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x, -y,  z

      }
      else if(space_group_name == "P1_c_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // x, -y,  z+0.5

      }
      else if(space_group_name == "C1_m_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x, -y,  z

         // + (0.5, 0.5, 0.0)
         sz = 2;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }


      }
      else if(space_group_name == "C1_c_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // x, -y,  z+0.5
         
         // + (0.5, 0.5, 0.0)
         sz = 2;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "P1_2/m_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  //  x, -y,  z

  
      }
      else if(space_group_name == "P1_21/m_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y+0.5,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] =  1.0;    operators.push_back(m);  //  x, 0.5-y,  z

      }
      else if(space_group_name == "C1_2/m_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  //  x, -y,  z
         
         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "P1_2/c_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x, y,  0.5-z
         m = 0.0;
         m.M[0] = -1.0;                                
                         m.M[5] = -1.0;                 
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                
                         m.M[5] = -1.0;                 
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // x, -y,  0.5+z

      }
      else if(space_group_name == "P1_21/c_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x, 0.5+y,  0.5-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // x, 0.5-y,  0.5+z

      }
      else if(space_group_name == "C1_2/c_1"){

      //std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Here are two wariants : 
      // a) - with origin at -1 on c-glide plane
        std::cout<<"Using variant with origin at -1 on c-glide plane"<<std::endl;
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;                 
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x, y,  0.5-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // x, -y,  0.5+z

         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }



      
      // b) - with origin at -1 on n-glide plane
        
      }

      //================== ORTHORHOMBIC =====================
 
      else if(space_group_name == "P2_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y,  z

      }
      else if(space_group_name == "P2_2_21"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "P21_21_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 0.5+x, 0.5-y,  -z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5; 
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 0.5-x, 0.5+y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y,  z

      }
      else if(space_group_name == "P21_2_21"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "C2_2_21"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x, y, 0.5-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x, -y, 0.5+z

          // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "C2_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z

         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "F2_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z

         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

         // + (0.5, 0.0, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

         // + (0.0, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }


      }
      else if(space_group_name == "I2_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z

         // + (0.5, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "I21_21_21"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // x, -y,  0.5-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 0.5-x, y,  -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, 0.5-y, z

         // + (0.5, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Pm_m_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z

      }
      else if(space_group_name == "Pm_c_21"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //  x, -y, 0.5+z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x, -y, 0.5+z

      }
      else if(space_group_name == "Pc_c_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x, y,  0.5+z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //  x, -y, 0.5+z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z

      }
      else if(space_group_name == "Pm_a_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // 0.5-x, y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // 0.5+x, -y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z

      }
      else if(space_group_name == "Pc_a_21"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // 0.5-x, y,  0.5+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // 0.5+x, -y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x, -y, 0.5+z

      }
      else if(space_group_name == "Pn_c_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x, 0.5+y,  0.5+z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //  x, 0.5-y, 0.5+z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z

      }
      else if(space_group_name == "Pm_n_21"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //  0.5+x, -y, 0.5+z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // 0.5-x, -y, 0.5+z

      }
      else if(space_group_name == "Pb_a_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;    operators.push_back(m);  // 0.5-x, 0.5+y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;    operators.push_back(m);  // 0.5+x, 0.5-y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z

      }
      else if(space_group_name == "Pn_a_21"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // 0.5-x, 0.5+y,  0.5+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;    operators.push_back(m);  //  0.5+x, 0.5-y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x, -y, 0.5+z

      }
      else if(space_group_name == "Pn_n_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // 0.5-x, 0.5+y,  0.5+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //  0.5+x, 0.5-y, 0.5+z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z

      }
      else if(space_group_name == "Cm_m_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z

         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

          
      }
      else if(space_group_name == "Cm_c_21"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //  x, -y, 0.5+z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x, -y, 0.5+z

         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         } 


      }
      else if(space_group_name == "Cc_c_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x, y,  0.5+z
         m = 0.0; 
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //  x, -y, 0.5+z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Am_m_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.0, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Ab_m_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, 0.5+y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, 0.5-y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.0, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Am_a_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 0.5-x, y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 0.5+x, -y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.0, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Ab_a_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 0.5-x, 0.5+y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 0.5+x, 0.5-y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.0, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Fm_m_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.0, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

         // + (0.5, 0.0, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Fd_d_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.25;
                         m.M[5] = 1.0;                  m.M[7] = 0.25;
                                        m.M[10] = 1.0;  m.M[11]= 0.25;   operators.push_back(m);  // 1/4-x, 1/4+y,  1/4+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.25;
                         m.M[5] = -1.0;                 m.M[7] = 0.25;
                                        m.M[10] = 1.0;  m.M[11]= 0.25;   operators.push_back(m);  //  1/4+x, 1/4-y, 1/4+z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.0, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }
         // + (0.5, 0.0, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }
         // + (0.5, 0.5, 0.0)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;   operators.push_back(m);
         }


      }
      else if(space_group_name == "Im_m_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.5, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Ib_a_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  // -x, y,  0.5+z
         m = 0.0;  
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  //  x, -y, 0.5+z
         m = 0.0;        
         m.M[0] = -1.0;  
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.5, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Im_a_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 0.5-x, y,  z
         m = 0.0;  
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  0.5+x, -y, z
         m = 0.0;        
         m.M[0] = -1.0;  
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.5, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;   operators.push_back(m);
         }

      }
      else if(space_group_name == "Pm_m_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                               
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0; 
         m.M[0] =  1.0;                                
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, -y, z
         m = 0.0; 
         m.M[0] =  1.0; 
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z
      
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;     operators.push_back(m);  // x,-y, -z
         m = 0.0; 
         m.M[0] = -1.0; 
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;     operators.push_back(m);  // -x, y,-z
         m = 0.0; 
         m.M[0] = -1.0; 
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z


      }
      else if(space_group_name == "Pn_n_n"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "Pc_c_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                               
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  // -x, y,  0.5+z
         m = 0.0; 
         m.M[0] =  1.0;                                
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  //  x, -y, 0.5+z
         m = 0.0; 
         m.M[0] =  1.0; 
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z
      

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;    operators.push_back(m);  // x,-y, 0.5-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;    operators.push_back(m);  // -x, y,0.5-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z


      }
      else if(space_group_name == "Pb_a_n"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 0.5-x, 0.5+y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 0.5+x, 0.5-y, z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 0.5+x, 0.5+y, -z

        m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 0.5-x, 0.5-y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;     operators.push_back(m);  // x,-y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;     operators.push_back(m);  // -x, y,-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z


      }
      else if(space_group_name == "Pm_m_a"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 1/2-x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;     operators.push_back(m);  // 1/2+x,-y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;     operators.push_back(m);  // -x, y,-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // 1/2-x, -y, z


      }
      else if(space_group_name == "Pn_n_a"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  // -x, 1/2+y,  1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  1/2+x, 1/2-y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5; 
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // x, 1/2-y, 1/2-z
         m = 0.0;                                                   
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // 1/2-x,1/2+y,1/2-z
         m = 0.0;                                                   
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // 1/2-x, -y, z


      }
      else if(space_group_name == "Pm_n_a"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  // 1/2-x, y,  1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  1/2+x, -y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2+x, y, 1/2-z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // 1/2-x,y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2-x, -y, 1/2+z


      }
      else if(space_group_name == "Pc_c_a"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  // 1/2-x, y,  1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  x, -y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/2+x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2+x, -y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // -x,y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // 1/2-x, -y, z


      }
      else if(space_group_name == "Pb_a_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 1/2-x, 1/2+y, z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  1/2+x, 1/2-y, z
         m = 0.0;
         m.M[0] =  1.0;                                 
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2-x,1/2+y,-z
         m = 0.0;
         m.M[0] = -1.0;                                
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z


      }
      else if(space_group_name == "Pc_c_n"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  // 1/2-x, y,  1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  1/2+x, 1/2-y, z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/2+x, 1/2+y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2+x, -y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // -x,1/2+y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] =  1.0;    operators.push_back(m);  // 1/2-x, 1/2-y, z


      }
      else if(space_group_name == "Pb_c_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, 1/2+y, z
         m = 0.0;
         m.M[0] =  1.0;                                
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  x, 1/2-y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // x, y, 1/2-z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;                                 
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // -x,1/2y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11]= 0.5;   operators.push_back(m);  // -x, -y, 1/2+z


      }
      else if(space_group_name == "Pn_n_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  // 1/2-x, 1/2+y,  1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  1/2+x, 1/2-y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2+x, 1/2-y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // 1/2-x,1/2+y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z


      }
      else if(space_group_name == "Pm_m_n"){

         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  y,  z

         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  -y,  z

         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] =-1.0;  m.M[11]= 0.0;   operators.push_back(m);  // 1/2+x, 1/2+y,  -z

         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] =-1.0;                  m.M[7] = 0.5;
                                        m.M[10] =-1.0;  m.M[11]= 0.0;   operators.push_back(m);  // 1/2-x, 1/2-y,  -z

         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =-1.0;                  m.M[7] = 0.5;
                                        m.M[10] =-1.0;  m.M[11]= 0.0;   operators.push_back(m);  // 1/2+x, 1/2-y,  -z

         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] =-1.0;  m.M[11]= 0.0;   operators.push_back(m);  // 1/2-x, 1/2+y,  -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z



      //std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "Pb_c_n"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 1/2-x, 1/2+y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 
                         m.M[5] = -1.0;                 
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  x, -y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2+x, 1/2+y, 1/2-z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;                                 
                         m.M[5] =  1.0;                 
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // -x,y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // 1/2-x, 1/2-y, 1/2+z


      }
      else if(space_group_name == "Pb_c_a"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  // 1/2-x, 1/2+y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  x, 1/2-y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                 
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2+x, y, 1/2-z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // -x,1/2+y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // 1/2-x, -y, 1/2+z


      }
      else if(space_group_name == "Pn_m_a"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  // 1/2-x, 1/2+y,  1/2+z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, 1/2-y, z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] =  1.0;                
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2+x, y, 1/2-z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // 1/2+x, 1/2-y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,1/2+y,-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // 1/2-x, -y, 1/2+z


      }
      else if(space_group_name == "Cm_c_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;                                 
                         m.M[5] = 1.0;                
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  x, -y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;                                 
                         m.M[5] =  1.0;               
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // x, y, 1/2-z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;                               
                         m.M[5] = -1.0;                
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11]= 0.5;     operators.push_back(m);  // -x,y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 
                         m.M[5] = -1.0;                 
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x, -y, 1/2+z
         // + (0.5, 0.5, 0.0)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; operators.push_back(m);
         }
      }
      else if(space_group_name == "Cm_c_a"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11]= 0.5;   operators.push_back(m);  //  x, 1/2-y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // x, 1/2+y, 1/2-z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  // -x,1/2+y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] =  1.0; m.M[11]= 0.5;   operators.push_back(m);  // -x, 1/2-y, 1/2+z
         // + (0.5, 0.5, 0.0)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; operators.push_back(m);
         }
      }
      else if(space_group_name == "Cm_m_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;                
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;                 
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;                 
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,y,-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;                
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.5, 0.5, 0.0)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; operators.push_back(m);
         }
      }
      else if(space_group_name == "Cc_c_m"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "Cm_m_a"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "Cc_c_a"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Three wariants
      }
      else if(space_group_name == "Fm_m_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,y,-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.5, 0.5, 0.0)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; operators.push_back(m);
         }
         // + (0.5, 0.0, 0.5)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }
         // + (0.0, 0.5, 0.5)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

      }
      else if(space_group_name == "Fd_d_d"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "Im_m_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  // -x, y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;     operators.push_back(m);  //  x, -y, z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,y,-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x, -y, z
         // + (0.5, 0.5, 0.5)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
      }
      else if(space_group_name == "Ib_a_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  // -x, y,  1/2+z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  //  x, -y, 1/2+z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x,y,1/2-z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x, -y, 1/2+z
         // + (0.5, 0.5, 0.5)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
      }
      else if(space_group_name == "Ib_c_a"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "Im_m_a"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }



      //================== TETRAGONAL =======================
      else if(space_group_name == "P_4"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z


      }
      else if(space_group_name == "P_41"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x,  -y,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.25;   operators.push_back(m);  // -y,  x,  1/4+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.75;   operators.push_back(m);  // y,  -x,  3/4+z


      }
      else if(space_group_name == "P_42"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -y,  x,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // y,  -x,  1/2+z


      }
      else if(space_group_name == "P_43"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x,  -y,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.75;   operators.push_back(m);  // -y,  x,  3/4+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.25;   operators.push_back(m);  // y,  -x,  1/4+z


      }
      else if(space_group_name == "I_4"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         // + (0.5, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }

      }
      else if(space_group_name == "I_41"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "P_-4"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y, -x, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-y,  x, -z


      }
      else if(space_group_name == "I_-4"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y, -x, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-y,  x, -z

         // + (0.5, 0.5, 0.5)
         sz = 4;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }


      }
      else if(space_group_name == "P_4/m"){
         m = 0.0; 
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0; 
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,-y, -z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,  y, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y, -x, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-y, x, -z 

      }
      else if(space_group_name == "P_42/m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -y,  x,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // y,  -x,  1/2+z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,-y, -z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,  y, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // y, -x, 1/2-z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  //-y, x, 1/2-z

      }
      else if(space_group_name == "P_4/n"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Four wariants
      }
      else if(space_group_name == "P_42/n"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "I_4/m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,-y, -z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,  y, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y, -x, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-y, x, -z

         // + (0.5, 0.5, 0.5)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }

      }
      else if(space_group_name == "I_41/a"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Ten wariants
      }
      else if(space_group_name == "P4_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,-y, -z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  y, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y, -x, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y, x, -z

      }
      else if(space_group_name == "P4_21_2"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "P41_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x,  -y,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.25;   operators.push_back(m);  // -y,  x,  1/4+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.75;   operators.push_back(m);  // y,  -x,  3/4+z

         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // x,-y, 1/2-z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  y, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0; m.M[11] = 0.25;   operators.push_back(m);  // -y, -x, 1/4-z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0; m.M[11] = 0.75;   operators.push_back(m);  // y, x, 3/4-z

      }
      else if(space_group_name == "P41_21_2"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Three wariants
      }
      else if(space_group_name == "P42_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -y,  x,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // y,  -x,  1/2+z

         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,-y, -z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  y, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // -y, -x, 1/2-z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // y, x, 1/2-z

      }
      else if(space_group_name == "P42_21_2"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "P43_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x,  -y,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.75;   operators.push_back(m);  // -y,  x,  3/4+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.25;   operators.push_back(m);  // y,  -x,  1/4+z

         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = 0.5;   operators.push_back(m);  // x,-y, 1/2-z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  y, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0; m.M[11] = 0.75;   operators.push_back(m);  // -y, -x, 3/4-z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0; m.M[11] = 0.25;   operators.push_back(m);  // y, x, 1/4-z

      }
      else if(space_group_name == "P43_21_2"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Three wariants
      }
      else if(space_group_name == "I4_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,-y, -z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  y, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y, -x, -z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y, x, -z

         // + (0.5, 0.5, 0.5)
         sz = 8;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }

      }
      else if(space_group_name == "I41_2_2"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two variants
      }
      else if(space_group_name == "P4_m_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         m = 0.0;
         m.M[0] =-1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x,y, z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] =-1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  -y, z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y, x, z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y, -x, z

      }
      else if(space_group_name == "P4_b_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         m = 0.0;
         m.M[0] =-1.0;                                  m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;    operators.push_back(m);  //1/2-x,1/2+y, z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;                  m.M[3] = 0.5;
         m.M[4] = 0.0;   m.M[5] =-1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;    operators.push_back(m);  // 1/2+x,  1/2-y, z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;                   m.M[3] = 0.5;
         m.M[4] = 1.0;  m.M[5] = 0.0;                   m.M[7] = 0.5;
                                        m.M[10] =  1.0;    operators.push_back(m);  // 1/2+y, 1/2+x, z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;                   m.M[3] = 0.5;
         m.M[4] =-1.0;  m.M[5] = 0.0;                   m.M[7] = 0.5;
                                        m.M[10] =  1.0;    operators.push_back(m);  // 1/2-y, 1/2-x, z

      }
      else if(space_group_name == "P42_c_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -y,  x,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // y,  -x,  1/2+z

         m = 0.0;
         m.M[0] =-1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //-x,y, 1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] =-1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // x,  -y, 1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y, x, z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y, -x, z

      }
      else if(space_group_name == "P42_n_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;                   m.M[3] = 0.5;
         m.M[4] = 1.0;  m.M[5] = 0.0;                   m.M[7] = 0.5;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // 1/2-y, 1/2+x,  1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;                   m.M[3] = 0.5;
         m.M[4] =-1.0;  m.M[5] = 0.0;                   m.M[7] = 0.5;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // 1/2+y, 1/2-x,  1/2+z

         m = 0.0;
         m.M[0] =-1.0;                                  m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  //1/2-x,1/2+y, 1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;                  m.M[3] = 0.5;
         m.M[4] = 0.0;   m.M[5] =-1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;  operators.push_back(m);  // 1/2+x, 1/2-y, 1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y, x, z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y, -x, z

      }
      else if(space_group_name == "P4_c_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         m = 0.0;
         m.M[0] =-1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  //-x,y, 1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] =-1.0;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  // x,  -y, 1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // y, x, 1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -y, -x, 1/2+z

      }
      else if(space_group_name == "P4_n_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 0.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x,  z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] = 1.0;
         m.M[4] =-1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x,  z

         m = 0.0;
         m.M[0] =-1.0;                                  m.M[3] = 0.5;
                         m.M[5] = 1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  //1/2-x,1/2+y, 1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = 0.0;                  m.M[3] = 0.5;
         m.M[4] = 0.0;   m.M[5] =-1.0;                  m.M[7] = 0.5;
                                        m.M[10] = 1.0;  m.M[11] = 0.5;   operators.push_back(m);  // 1/2+x, 1/2-y, 1/2+z
         m = 0.0; 
         m.M[0] = 0.0;  m.M[1] = 1.0;                   m.M[3] = 0.5;
         m.M[4] = 1.0;  m.M[5] = 0.0;                   m.M[7] = 0.5;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // 1/2+y, 1/2+x, 1/2+z
         m = 0.0;
         m.M[0] = 0.0;  m.M[1] =-1.0;                   m.M[3] = 0.5;
         m.M[4] =-1.0;  m.M[5] = 0.0;                   m.M[7] = 0.5;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // 1/2-y, 1/2-x, 1/2+z

      }











      //================== TRIGONAL =========================

      else if(space_group_name == "P_3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;  
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


      }
      else if(space_group_name == "P_31"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -y,  x-y,  1/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -x+y,  -x,  2/3+z


      }
      else if(space_group_name == "P_32"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -y,  x-y,  2/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -x+y,  -x,  1/3+z


      }
      else if(space_group_name == "R_3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;   
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x


      }
      else if(space_group_name == "P_-3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  -x+y, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y, x, -z

      }
      else if(space_group_name == "R_-3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -y, -z
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] = -1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // -z, -x, -y
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -y, -z, -x

      }
      else if(space_group_name == "P3_1_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
                         m.M[1] = -1.0;  
        m.M[4] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y, -x, -z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] =  0.0;   m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y,  y, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, x-y, -z

      }
      else if(space_group_name == "P3_2_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  x, -z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] =  -1.0;
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  -y, -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x, -x+y, -z

      }
      else if(space_group_name == "P31_1_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -y,  x-y,  1/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -x+y,  -x,  2/3+z

         m = 0.0;
                         m.M[1] = -1.0;
        m.M[4] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -y, -x, 2/3-z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] =  0.0;   m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -x+y,  y, 1/3-z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, x-y, -z

      }
      else if(space_group_name == "P31_2_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -y,  x-y,  1/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;               
                                        m.M[10] = 1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // -x+y, -x,  2/3+z
         m = 0.0;
                         m.M[1] =  1.0;
         m.M[4] = 1.0;                
                                        m.M[10] = -1.0;   operators.push_back(m);  //  y, x, -z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // x-y, -y, 2/3-z
         m = 0.0;
         m.M[0] = -1.0; 
         m.M[4] = -1.0;  m.M[5] = 1.0; 
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -x, -x+y, 1/3-z


      }
      else if(space_group_name == "P32_1_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -y,  x-y,  2/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -x+y,  -x,  1/3+z

         m = 0.0;
                         m.M[1] = -1.0;
        m.M[4] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -y, -x, 1/3-z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] =  0.0;   m.M[5] =  1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -x+y,  y, 2/3-z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, x-y, -z

      }// P32_1_2     
      else if(space_group_name == "P32_2_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -y,  x-y,  2/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;               
                                        m.M[10] = 1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // -x+y, -x,  1/3+z
         m = 0.0;
                         m.M[1] =  1.0;
         m.M[4] = 1.0;                
                                        m.M[10] = -1.0;   operators.push_back(m);  //  y, x, -z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // x-y, -y, 1/3-z
         m = 0.0;
         m.M[0] = -1.0; 
         m.M[4] = -1.0;  m.M[5] = 1.0; 
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -x, -x+y, 2/3-z


      }// P32_2_1
      else if(space_group_name == "R_32"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;   
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -y, -z
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] = -1.0;   
                        m.M[9] =  -1.0;                     operators.push_back(m);  // -z, -x, -y
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -y, -z, -x


      }// R_32
      else if(space_group_name == "P_3_m_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;  
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] =  -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  -x,  z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] = 1.0;  
         m.M[4] = 0.0;   m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  y,  z
         m = 0.0;
         m.M[0] =  1.0;  m.M[1] = 0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  x-y,  z


      }// P_3_m_1
      else if(space_group_name == "P_3_1_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;  
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  x,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;  
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x-y, -y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -x+y,  z


      }// P_3_1_m
      else if(space_group_name == "P_3_c_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;  
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] =  -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -y,  -x,  1/2+z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] = 1.0;  
         m.M[4] = 0.0;   m.M[5] =  1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x+y,  y,  1/2+z
         m = 0.0;
         m.M[0] =  1.0;  m.M[1] = 0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0; 
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // x,  x-y,  1/2+z


      }// P_3_c_1
      else if(space_group_name == "P_3_1_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;  
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] =  1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // y,  x,  1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;  
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // x-y, -y, 1/2+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0; 
                                        m.M[10] =  1.0; m.M[11] = 0.5;   operators.push_back(m);  // -x,  -x+y, 1/2+z


      }// P_3_1_c
      else if(space_group_name == "R_3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;   
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x

         m = 0.0;
         m.M[0] =  1.0;
                                        m.M[6] =  1.0;
                       m.M[9] =  1.0;                      operators.push_back(m);  // x,  z,  y
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;   
         m.M[8] =  1.0;                                    operators.push_back(m);  // z, y, x
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y, x, z


      }// R_3_m
      else if(space_group_name == "R_3_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;   
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x

         m = 0.0;  m.M[3] = 0.5;  m.M[7] = 0.5; m.M[11] = 0.5;
         m.M[0] =  1.0;
                                        m.M[6] =  1.0;
                       m.M[9] =  1.0;                      operators.push_back(m);  // 1/2+x,  1/2+z,  1/2+y
         m = 0.0;  m.M[3] = 0.5;  m.M[7] = 0.5; m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;   
         m.M[8] =  1.0;                                    operators.push_back(m);  // 1/2+z, 1/2+y, 1/2+x
         m = 0.0;  m.M[3] = 0.5;  m.M[7] = 0.5; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // 1/2+y, 1/2+x, 1/2+z


      }// R_3_c
      else if(space_group_name == "P_-3_1_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;  
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  x,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;  
         m.M[4] = 0.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x-y, -y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0; 
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -x+y,  z


	 // Add inversion
         sz = 6;
         for(i=0;i<sz;i++){
         m = operators[i];
         operators.push_back((-m));
         }


      }// P_-3_1_m
      else if(space_group_name == "P_-3_1_c"){
      // Two variants
      }// P_-3_1_c
      else if(space_group_name == "P_-3_m_1"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] =  -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  -x,  z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] = 1.0;
         m.M[4] = 0.0;   m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  y,  z
         m = 0.0;
         m.M[0] =  1.0;  m.M[1] = 0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  x-y,  z

          // Add inversion
         sz = 6;
         for(i=0;i<sz;i++){
         m = operators[i];
         operators.push_back((-m));
         }

      }// P_-3_m_1
      else if(space_group_name == "P_-3_c_1"){
      // Two variants
      }// P_-3_c_1
      else if(space_group_name == "R_-3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x

         m = 0.0;
         m.M[0] =  1.0;
                                        m.M[6] =  1.0;
                       m.M[9] =  1.0;                      operators.push_back(m);  // x,  z,  y
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // z, y, x
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y, x, z

         // Add inversion
         sz = 6;
         for(i=0;i<sz;i++){
         m = operators[i];
         operators.push_back((-m));
         }


      }// R_-3_m
      else if(space_group_name == "R_-3_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x

         m = 0.0;  m.M[3] = 0.5;  m.M[7] = 0.5; m.M[11] = 0.5;
         m.M[0] =  1.0;
                                        m.M[6] =  1.0;
                       m.M[9] =  1.0;                      operators.push_back(m);  // 1/2+x,  1/2+z,  1/2+y
         m = 0.0;  m.M[3] = 0.5;  m.M[7] = 0.5; m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // 1/2+z, 1/2+y, 1/2+x
         m = 0.0;  m.M[3] = 0.5;  m.M[7] = 0.5; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // 1/2+y, 1/2+x, 1/2+z


         // Add inversion
         sz = 6;
         for(i=0;i<sz;i++){
         m = operators[i];
         operators.push_back((-m));
         }


      }// R_-3_c


      //================== HEXAGONAL ========================
      else if(space_group_name == "P_6"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x-y,  x,  z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // y,  -x+y,  z

      }// P_6

      else if(space_group_name == "P_61"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // -y,  x-y,  1/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // -x+y,  -x,  2/3+z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x,  -y,  1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] =  1.0;  m.M[11] = (1.0/6.0);  operators.push_back(m);  // x-y,  x,  1/6+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] =  1.0;  m.M[11] = (5.0/6.0);  operators.push_back(m);  // y,  -x+y,  5/6+z

      }// P_61

      else if(space_group_name == "P_65"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // -y,  x-y,  2/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // -x+y,  -x,  1/3+z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;  m.M[11] = 0.5;  operators.push_back(m);  // -x,  -y,  1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] =  1.0;  m.M[11] = (5.0/6.0);  operators.push_back(m);  // x-y,  x,  5/6+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] =  1.0;  m.M[11] = (1.0/6.0);  operators.push_back(m);  // y,  -x+y,  1/6+z

      }// P_65

      else if(space_group_name == "P_62"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // -y,  x-y,  2/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // -x+y,  -x,  1/3+z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] =  1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // x-y,  x,  1/3+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] =  1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // y,  -x+y,  2/3+z

      }// P_62

      if(space_group_name == "P_64"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // -y,  x-y,  1/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // -x+y,  -x,  2/3+z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x,  -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] =  1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // x-y,  x,  2/3+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] =  1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // y,  -x+y,  1/3+z

      }// P_64

      else if(space_group_name == "P_63"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;  operators.push_back(m);  // -x,  -y,  1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;  operators.push_back(m);  // x-y,  x,  1/2+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] =  1.0; m.M[11] = 0.5;  operators.push_back(m);  // y,  -x+y,  1/2+z

      }// P_63

      else if(space_group_name == "P_-6"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,  y,  -z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = -1.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, -z

      }// P_-6

      else if(space_group_name == "P_6/m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  -x+y,  -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  x,  -z


         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,  y,  -z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = -1.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, -z

      }// P_6/m

      else if(space_group_name == "P_63/m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


  	 m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, 1/2+z

         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  -x+y,  -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  x,  -z


         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,  y,  1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] = -1.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, 1/2-z

      }// P_63/m

      else if(space_group_name == "P_6_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, z


         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  x, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 0.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  -y,  -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -x+y,  -z


         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  -x,  -z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] =  1.0;  m.M[1] =  0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,  x-y, -z

      }// P_6_2_2

      else if(space_group_name == "P_61_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -y,  x-y,  1/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -x+y,  -x,  2/3+z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = (1.0/2.0);  operators.push_back(m);  // -x, -y,  1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;  m.M[11] = (1.0/6.0);  operators.push_back(m);  // x-y, x, 1/6+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = (5.0/6.0);  operators.push_back(m);  // y,  -x+y, 5/6+z


         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // y,  x, 1/3-z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 0.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  -y,  -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -x,  -x+y,  2/3-z


         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (5.0/6.0);   operators.push_back(m);  // -y,  -x,  5/6-z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/2.0);   operators.push_back(m);  // -x+y, -x, 1/2-z
         m = 0.0;
         m.M[0] =  1.0;  m.M[1] =  0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/6.0);   operators.push_back(m);  // x,  x-y, 1/6-z

      }// P_61_2_2

      else if(space_group_name == "P_65_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -y,  x-y,  2/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -x+y,  -x,  1/3+z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;  m.M[11] = (1.0/2.0);  operators.push_back(m);  // -x, -y,  1/2+z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;  m.M[11] = (5.0/6.0);  operators.push_back(m);  // x-y, x, 5/6+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = (1.0/6.0);  operators.push_back(m);  // y,  -x+y, 1/6+z


         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // y,  x, 2/3-z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 0.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  -y,  -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -x,  -x+y,  1/3-z


         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/6.0);   operators.push_back(m);  // -y,  -x,  1/6-z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/2.0);   operators.push_back(m);  // -x+y, -x, 1/2-z
         m = 0.0;
         m.M[0] =  1.0;  m.M[1] =  0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (5.0/6.0);   operators.push_back(m);  // x,  x-y, 5/6-z

      }// P_65_2_2

      else if(space_group_name == "P_62_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -y,  x-y,  2/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -x+y,  -x,  1/3+z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -x, -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // x-y, x, 1/3+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // y,  -x+y, 2/3+z


         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // y,  x, 2/3-z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 0.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  -y,  -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -x,  -x+y,  1/3-z


         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -y,  -x,  2/3-z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] =  1.0;  m.M[1] =  0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // x,  x-y, 1/3-z

      }// P_62_2_2

      else if(space_group_name == "P_64_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -y,  x-y,  1/3+z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -x+y,  -x,  2/3+z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -x, -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;  m.M[11] = (2.0/3.0);  operators.push_back(m);  // x-y, x, 2/3+z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;  m.M[11] = (1.0/3.0);  operators.push_back(m);  // y,  -x+y, 1/3+z


         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // y,  x, 1/3-z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 0.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  -y,  -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // -x,  -x+y,  2/3-z


         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (1.0/3.0);   operators.push_back(m);  // -y,  -x,  1/3-z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] =  1.0;  m.M[1] =  0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0; m.M[11] = (2.0/3.0);   operators.push_back(m);  // x,  x-y, 2/3-z

      }// P_64_2_2

      else if(space_group_name == "P_63_2_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


	 m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, 1/2+z


         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  x, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 0.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  -y,  -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 0.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -x+y,  -z


         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  -x,  1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;   m.M[1] =  1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  1.0;  m.M[1] =  0.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x,  x-y, 1/2-z

      }// P_63_2_2

      else if(space_group_name == "P_6_m_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, z



         m = 0.0;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, z
         m = 0.0;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  z


         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, z

      }// P_6_m_m

      else if(space_group_name == "P_6_c_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, z



 	 m = 0.0; m.M[11] = 0.5;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  1/2+z


         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, 1/2+z

      }// P_6_c_c

      else if(space_group_name == "P_63_c_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


	 m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, 1/2+z



         m = 0.0; m.M[11] = 0.5;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  1/2+z


         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, z

      }// P_63_c_m

      else if(space_group_name == "P_63_m_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


	 m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, 1/2+z



         m = 0.0;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, z
         m = 0.0;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  z


         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, 1/2+z

      }// P_63_m_c

      else if(space_group_name == "P_-6_m_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y, -z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] =  1.0;
         m.M[4] =-1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x+y, -x, -z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] =-1.0;
         m.M[4] =  1.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, -z



         m = 0.0;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, z
         m = 0.0;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  z


         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-y, -x, -z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 1.0;
         m.M[4] =-1.0;   m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y,-x,-z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] =  0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, x-y, -z

      }// P_-6_m_2

      else if(space_group_name == "P_-6_c_2"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two variants

      }// P_-6_c_2

      else if(space_group_name == "P_-6_2_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z


         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y,  -z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 1.0;
         m.M[4] =-1.0;   m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x+y,-x, -z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = -1.0;
         m.M[4] =  1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, -z



         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  x, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] =-1.0;
         m.M[4] = 0.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y, -y, -z
         m = 0.0;
         m.M[0] =-1.0;  m.M[1] = 0.0;
         m.M[4] =-1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -x+y, -z


         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, z

      }// P_-6_2_m

      else if(space_group_name == "P_-6_2_c"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two variants

      }// P_-6_2_c

      else if(space_group_name == "P_6/mmm"){
//1
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

// 4
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, z

//-------------------
// 7
         m = 0.0;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, z
         m = 0.0;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  z

// 10
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, z


//-------------------
//13
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  -x+y,  -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  x, -z

//16
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y,  -z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 1.0;
         m.M[4] =-1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] =-1.0;
         m.M[4] =  1.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, -z

//----------------
// 19
         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  x, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] =-1.0;
         m.M[4] = 0.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y, -y,  -z
         m = 0.0;
         m.M[0] =-1.0;  m.M[1] = 0.0;
         m.M[4] =-1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -x+y,  -z

// 22
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  -x,  -z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] =  0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, x-y, -z


      }// P_6/mmm

      else if(space_group_name == "P_6/mcc"){
//1
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

// 4
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, z

//-------------------
// 7
         m = 0.0; m.M[11] = 0.5;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  1/2+z

// 10
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, 1/2+z


//-------------------
//13
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;   m.M[5] = 1.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // y,  -x+y,  -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // x-y,  x,  -z

//16
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y,  -z
         m = 0.0;
         m.M[0] =-1.0;   m.M[1] = 1.0;
         m.M[4] =-1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] =  0.0;  m.M[1] =-1.0;
         m.M[4] =  1.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, -z

//----------------
// 19
         m = 0.0; m.M[11] = 0.5;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;  m.M[1] =-1.0;
         m.M[4] = 0.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y, -y,  1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =-1.0;  m.M[1] = 0.0;
         m.M[4] =-1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -x+y,  1/2-z

// 22
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  -x,  1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;  m.M[1] =  0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, x-y, 1/2-z


      }// P_6/mcc

      else if(space_group_name == "P_63/mcm"){
//1
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

// 4
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, 1/2+z

//-------------------
// 7
         m = 0.0; m.M[11] = 0.5;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  1/2+z

// 10
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  z
         m = 0.0;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, z


//-------------------
//13
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  -x+y,  -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  x, -z

//16
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y,  1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =-1.0;   m.M[1] = 1.0;
         m.M[4] =-1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] =-1.0;
         m.M[4] =  1.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, 1/2-z

//----------------
// 19
         m = 0.0; m.M[11] = 0.5;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;  m.M[1] =-1.0;
         m.M[4] = 0.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y, -y,  1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =-1.0;  m.M[1] = 0.0;
         m.M[4] =-1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -x+y,  1/2-z

// 22
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  -x,  -z
         m = 0.0;
         m.M[0] = -1.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] =  0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, x-y, -z


      }// P63/mcm

      else if(space_group_name == "P_63/mmc"){
//1
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] = -1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -y,  x-y,  z
         m = 0.0;
         m.M[0] = -1.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 0.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // -x+y,  -x,  z

// 4
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -y,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] = 1.0;
         m.M[4] = -1.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  -x+y, 1/2+z

//-------------------
// 7
         m = 0.0;
                        m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -y,  -x, z
         m = 0.0;
         m.M[0] = -1.0; m.M[1] = 1.0;
         m.M[4] = 0.0;  m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x+y,  y,  z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = 0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x,  x-y,  z


// 10
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // y,  x,  1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;   m.M[1] = -1.0;
         m.M[4] = 1.0;   m.M[5] =  0.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x-y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;  m.M[1] =  0.0;
         m.M[4] = -1.0;  m.M[5] =  1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, -x+y, 1/2+z


//-------------------
//13
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] = 0.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;   m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  -x+y,  -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] = -1.0;
         m.M[4] = 1.0;  m.M[5] = 0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y,  x, -z

//16
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, y,  1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =-1.0;   m.M[1] = 1.0;
         m.M[4] =-1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] =  0.0;  m.M[1] =-1.0;
         m.M[4] =  1.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  x-y, 1/2-z
//----------------
// 19
         m = 0.0;
                        m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // y,  x, -z
         m = 0.0;
         m.M[0] = 1.0;  m.M[1] =-1.0;
         m.M[4] = 0.0;  m.M[5] =-1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x-y, -y,  -z
         m = 0.0;
         m.M[0] =-1.0;  m.M[1] = 0.0;
         m.M[4] =-1.0;  m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x,  -x+y,  -z

// 22
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -y,  -x,  1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;   m.M[1] = 1.0;
         m.M[4] = -1.0;   m.M[5] =  0.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // -x+y, -x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;  m.M[1] =  0.0;
         m.M[4] = 1.0;  m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, x-y, 1/2-z


      }// P_63/mmc










     //===================== CUBIC ============================

      else if(space_group_name == "P2_3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x


      }
      else if(space_group_name == "F2_3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x


         // + (0.0, 0.5, 0.5)
         sz = 12;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
         // + (0.5, 0.0, 0.5)
         sz = 12;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
         // + (0.5, 0.5, 0.0)
         sz = 12;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;  operators.push_back(m);
         }

      }// end of F2_3
      else if(space_group_name == "F2_3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x


         // + (0.5, 0.5, 0.5)
         sz = 12;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }

      }// end of I2_3
      else if(space_group_name == "P21_3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                      	 m.M[5] = -1.0;                 m.M[7] = 0.5; 
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;
                   	 m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  //-x,  1/2+y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] =  1.0; m.M[11]= 0.5;   operators.push_back(m);  //1/2-x, -y,  1/2+z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;  m.M[3] = 0.5;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                 m.M[11]= 0.5;    operators.push_back(m);  //1/2-z,-x, 1/2+y
         m = 0.0;
                                        m.M[2] = 1.0;  m.M[3] = 0.5;
	 m.M[4] =-1.0;                                 m.M[7] = 0.5;
                        m.M[9] = -1.0;                     operators.push_back(m);  // 1/2+z,1/2-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
	 m.M[4] = 1.0;                                 m.M[7] = 0.5;
                        m.M[9] = -1.0;                 m.M[11]= 0.5;    operators.push_back(m);  //-z, 1/2+x,1/2-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
           			        m.M[6] = 1.0;  m.M[7] = 0.5;
         m.M[8] = -1.0;                                m.M[11]= 0.5;    operators.push_back(m);  //-y, 1/2+z,1/2-x
         m = 0.0;
                	 m.M[1] =-1.0;                 m.M[3] = 0.5;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                m.M[11]= 0.5;    operators.push_back(m);  //1/2-y, -z, 1/2+x
         m = 0.0;
                	 m.M[1] = 1.0;                 m.M[3] = 0.5;
                		        m.M[6] =-1.0;  m.M[7] = 0.5;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2+y, 1/2-z, -x


      }// end of P21_3
      else if(space_group_name == "I21_3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                      	 m.M[5] = -1.0;                 m.M[7] = 0.5; 
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;
                   	 m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  //-x,  1/2+y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0; m.M[11]= 0.5;   operators.push_back(m);  //1/2-x, -y,  1/2+z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;  m.M[3] = 0.5;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                 m.M[11]= 0.5;    operators.push_back(m);  //1/2-z,-x, 1/2+y
         m = 0.0;
                                        m.M[2] = 1.0;  m.M[3] = 0.5;
	 m.M[4] =-1.0;                                 m.M[7] = 0.5;
                        m.M[9] = -1.0;                     operators.push_back(m);  // 1/2+z,1/2-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
	 m.M[4] = 1.0;                                 m.M[7] = 0.5;
                        m.M[9] = -1.0;                 m.M[11]= 0.5;    operators.push_back(m);  //-z, 1/2+x,1/2-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
           			        m.M[6] = 1.0;  m.M[7] = 0.5;
         m.M[8] = -1.0;                                m.M[11]= 0.5;    operators.push_back(m);  //-y, 1/2+z,1/2-x
         m = 0.0;
                	 m.M[1] =-1.0;                 m.M[3] = 0.5;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                m.M[11]= 0.5;    operators.push_back(m);  //1/2-y, -z, 1/2+x
         m = 0.0;
                	 m.M[1] = 1.0;                 m.M[3] = 0.5;
                		        m.M[6] =-1.0;  m.M[7] = 0.5;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2+y, 1/2-z, -x

      
         // + (0.5, 0.5, 0.5)
         sz = 12;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }

      }// end of I21_3
      else if(space_group_name == "P_m_-3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  //-x, y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x, -y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, y,  z



         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = -1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // -z,-x,-y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z, x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z,-x,y



         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //y, z,-x
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  //-y, z, x


      }// end of P_m_-3
      else if(space_group_name == "P_n_-3"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two variants
      }
      else if(space_group_name == "F_m_-3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  //-x, y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x, -y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, y,  z



         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = -1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // -z,-x,-y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z, x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z,-x,y



         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //y, z,-x
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  //-y, z, x

 
         // + (0.0, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
         // + (0.5, 0.0, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
         // + (0.5, 0.5, 0.0)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;  operators.push_back(m);
         }

      }// end of F_m_-3
      else if(space_group_name == "F_d_-3"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two variants
      }
      else if(space_group_name == "I_m_-3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  //-x, y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x, -y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, y, -z



         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = -1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // -z,-x,-y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z, x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z,-x,y



         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //y, z,-x
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  //-y, z, x

 
         // + (0.5, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }


      }// end of I_m_-3
      else if(space_group_name == "P_a_-3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;                                
                	 m.M[5] =  1.0;                 m.M[7] = 0.5;
                                        m.M[10] = -1.0; m.M[11]= 0.5;   operators.push_back(m);  //-x, 1/2+y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;                                 m.M[3] = 0.5;
                         m.M[5] = -1.0;               
                                        m.M[10] = 1.0;  m.M[11]= 0.5;    operators.push_back(m);  //1/2-x, -y,  1/2+z


         m = 0.0;
                               	        m.M[2] = 1.0;  
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                 	                m.M[2] =-1.0;  m.M[3] = 0.5;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                 m.M[11]= 0.5;     operators.push_back(m);  //1/2-z,-x, 1/2+y
         m = 0.0;
                           	        m.M[2] = 1.0;  m.M[3] = 0.5;
	 m.M[4] =-1.0;                                 m.M[7] = 0.5;
                        m.M[9] = -1.0;                     operators.push_back(m);  // 1/2+z,1/2-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
 	 m.M[4] = 1.0;                                 m.M[7] = 0.5;
                        m.M[9] = -1.0;                 m.M[11]= 0.5;    operators.push_back(m);  //-z, 1/2+x,1/2-y



         m = 0.0;
                         m.M[1] = 1.0;
          	         	        m.M[6] = 1.0;  
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                 	  	        m.M[6] = 1.0;  m.M[7] = 0.5;
         m.M[8] = -1.0;                                m.M[11]= 0.5;     operators.push_back(m);  //-y, 1/2+z,1/2-x
         m = 0.0;
                     	 m.M[1] =-1.0;                 m.M[3] = 0.5;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                m.M[11]= 0.5;    operators.push_back(m);  //1/2-y, -z, 1/2+x
         m = 0.0;
                	 m.M[1] = 1.0;                 m.M[3] = 0.5;
                		        m.M[6] =-1.0;  m.M[7] = 0.5;
         m.M[8] = -1.0;                                      operators.push_back(m);  // 1/2+y, 1/2-z, -x





         m = 0.0;  
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
	 m = 0.0;   m.M[3] = m.M[7] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  //1/2-x, 1/2+y, z
  	 m = 0.0;   m.M[7] = m.M[11]= 0.5;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // x, 1/2-y,1/2+z
	 m = 0.0;   m.M[3] = m.M[11] = 0.5;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //1/2+x, y,  1/2-z



         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = -1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // -z,-x,-y
         m = 0.0;   m.M[3] = m.M[11] = 0.5;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // 1/2+z, x,1/2-y
         m = 0.0;   m.M[3] = m.M[7] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // 1/2-z,1/2+x,y
         m = 0.0;   m.M[7] = m.M[11] = 0.5;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z,1/2-x,1/2+y



         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;   m.M[7] = m.M[11] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, 1/2-z, 1/2+x
         m = 0.0;   m.M[3] = m.M[11] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //1/2+y, z,1/2-x
         m = 0.0;   m.M[3] = m.M[7] = 0.5;
                         m.M[1] = -1.0;
                                        m.M[6] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  //1/2-y, 1/2+z, x

 


      }// end of P_a_-3
      else if(space_group_name == "I_a_-3"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
	 m = 0.0;  m.M[11] = 0.5;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, 1/2-z
         m = 0.0;  m.M[3] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //1/2-x,  y, -z
         m = 0.0;  m.M[7] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, 1/2-y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;  m.M[7] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,1/2-x, y
         m = 0.0;  m.M[11] = 0.5;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,1/2-y
         m = 0.0;  m.M[3] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //1/2-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;  m.M[3] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //1/2-y, z,-x
         m = 0.0;  m.M[7] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, 1/2-z, x
         m = 0.0;  m.M[11] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, 1/2-x





         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x, -y, -z
         m = 0.0;  m.M[11] = 0.5;
         m.M[0] = -1.0; 
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  //-x, y, 1/2+z
         m = 0.0;  m.M[3] = 0.5;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // 1/2+x, -y, z
         m = 0.0;  m.M[7] = 0.5;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, 1/2+y, -z



         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = -1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // -z,-x,-y
         m = 0.0;  m.M[7] = 0.5;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z, 1/2+x,-y
         m = 0.0;  m.M[11] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,1/2+y
         m = 0.0;  m.M[3] = 0.5;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //1/2+z,-x,y



         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;  m.M[3] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2+y, -z, x
         m = 0.0;  m.M[7] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //y, 1/2+z,-x
         m = 0.0;  m.M[11] = 0.5;
                         m.M[1] = -1.0;
                                        m.M[6] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  //-y, z, 1/2+x

 
         // + (0.5, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }


      }// end of I_a_-3  
      else if(space_group_name == "P_4_3_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, z, -y
         m = 0.0;
         m.M[0] =-1.0;
                                        m.M[6] = -1.0;
                         m.M[9] =-1.0;                     operators.push_back(m);  //-x, -z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, z, y



         m = 0.0;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, y, -x
         m = 0.0;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, y, x
         m = 0.0;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, -y, -x
         m = 0.0;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, -y, x



         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, x, z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] =-1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, -x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // -y, -x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // y, x, -z

 

      }// end of P_4_3_2
      else if(space_group_name == "P_42_3_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/2+x, 1/2-z, 1/2+y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/2+x, 1/2+z, 1/2-y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
         m.M[0] =-1.0;
                                        m.M[6] = -1.0;
                         m.M[9] =-1.0;                     operators.push_back(m);  //1/2-x, 1/2-z, 1/2-y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/2+x, 1/2+z, 1/2+y



         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2+z, 1/2+y, 1/2-x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2-z, 1/2+y, 1/2+x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2-z, 1/2-y, 1/2-x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2+z, 1/2-y, 1/2+x



         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/2-y, 1/2+x, 1/2+z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] =-1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/2+y, 1/2-x, 1/2+z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 1/2-y, 1/2-x, 1/2-z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 1/2+y, 1/2+x, 1/2-z

 

      }// end of P_42_3_2
      else if(space_group_name == "F_4_3_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, z, -y
         m = 0.0;
         m.M[0] =-1.0;
                                        m.M[6] = -1.0;
                         m.M[9] =-1.0;                     operators.push_back(m);  //-x, -z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, z, y



         m = 0.0;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, y, -x
         m = 0.0;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, y, x
         m = 0.0;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, -y, -x
         m = 0.0;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, -y, x



         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, x, z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] =-1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, -x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // -y, -x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // y, x, -z

 
         // + (0.0, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
         // + (0.5, 0.0, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
         // + (0.5, 0.5, 0.0)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;  operators.push_back(m);
         }


      }// end of F_4_3_2
      else if(space_group_name == "F_41_3_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/4+x, 1/4-z, 1/4+y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/4+x, 1/4+z, 1/4-y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
         m.M[0] =-1.0;
                                        m.M[6] = -1.0;
                         m.M[9] =-1.0;                     operators.push_back(m);  //1/4-x, 1/4-z, 1/4-y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/4-x, 1/4+z, 1/4+y



         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4+z, 1/4+y, 1/4-x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4-z, 1/4+y, 1/4+x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4-z, 1/4-y, 1/4-x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4+z, 1/4-y, 1/4+x



         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/4-y, 1/4+x, 1/4+z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] =-1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/4+y, 1/4-x, 1/4+z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 1/4-y, 1/4-x, 1/4-z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 1/4+y, 1/4+x, 1/4-z


         // + (0.0, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
         // + (0.5, 0.0, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }
         // + (0.5, 0.5, 0.0)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;  operators.push_back(m);
         }


      }// end of F_41_3_2
      else if(space_group_name == "I_4_3_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, z, -y
         m = 0.0;
         m.M[0] =-1.0;
                                        m.M[6] = -1.0;
                         m.M[9] =-1.0;                     operators.push_back(m);  //-x, -z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, z, y



         m = 0.0;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, y, -x
         m = 0.0;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, y, x
         m = 0.0;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, -y, -x
         m = 0.0;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, -y, x



         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, x, z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] =-1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, -x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // -y, -x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // y, x, -z

 
         // + (0.5, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5;  operators.push_back(m);
         }


      }// end of I_4_3_2
      else if(space_group_name == "P_43_3_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;  m.M[3] = m.M[7] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[7] = m.M[11] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  1/2+y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;  m.M[3] = m.M[11] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //1/2-x, -y,  1/2+z

//5
         m = 0.0;        
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;       m.M[3] = m.M[11] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //1/2-z,-x, 1/2+y
         m = 0.0;       m.M[3] = m.M[7] = 0.5;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // 1/2+z,1/2-x,-y
         m = 0.0;       m.M[7] = m.M[11] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, 1/2+x,1/2-y


//9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;        m.M[7] = m.M[11] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, 1/2+z,1/2-x
         m = 0.0;        m.M[3] = m.M[11] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //1/2-y, -z, 1/2+x
         m = 0.0;        m.M[3] = m.M[7] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2+y, 1/2-z, -x





         m = 0.0;  m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.25;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //3/4+x, 3/4-z, 1/4+y
         m = 0.0;  m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75; 
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/4+x, 3/4+z, 3/4-y
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
         m.M[0] =-1.0;
                                        m.M[6] = -1.0;
                         m.M[9] =-1.0;                     operators.push_back(m);  //1/4-x, 1/4-z, 1/4-y
         m = 0.0;  m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.75;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //3/4+x, 1/4+z, 3/4+y


// 17
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4+z, 3/4+y, 3/4-x
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.75;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 3/4-z, 1/4+y, 3/4+x
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4-z, 1/4-y, 1/4-x
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.25;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 3/4+z, 3/4-y, 1/4+x


// 21
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.75; 
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 3/4-y, 1/4+x, 3/4+z
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] =-1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 3/4+y, 3/4-x, 1/4+z
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 1/4-y, 1/4-x, 1/4-z
         m = 0.0;    m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 1/4+y, 3/4+x, 3/4-z


      }// end of P_43_3_2
      else if(space_group_name == "P_41_3_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;  m.M[3] = m.M[7] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // 1/2+x, 1/2-y, -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[7] = m.M[11] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  1/2+y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;  m.M[3] = m.M[11] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //1/2-x, -y,  1/2+z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;       m.M[3] = m.M[11] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //1/2-z,-x, 1/2+y
         m = 0.0;       m.M[3] = m.M[7] = 0.5;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // 1/2+z,1/2-x,-y
         m = 0.0;       m.M[7] = m.M[11] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, 1/2+x,1/2-y

// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;        m.M[7] = m.M[11] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, 1/2+z,1/2-x
         m = 0.0;        m.M[3] = m.M[11] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //1/2-y, -z, 1/2+x
         m = 0.0;        m.M[3] = m.M[7] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2+y, 1/2-z, -x




//13
         m = 0.0;  m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.75;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/4+x, 1/4-z, 3/4+y
         m = 0.0;  m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.25;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //3/4+x, 1/4+z, 1/4-y
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.75;
         m.M[0] =-1.0;
                                        m.M[6] = -1.0;
                         m.M[9] =-1.0;                     operators.push_back(m);  //3/4-x, 3/4-z, 3/4-y
         m = 0.0;  m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.25;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/4+x, 3/4+z, 1/4+y

// 17
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.25;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 3/4+z, 1/4+y, 1/4-x
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.25;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4-z, 3/4+y, 1/4+x
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.75;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 3/4-z, 3/4-y, 3/4-x
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.75;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4+z, 1/4-y, 3/4+x


// 21
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/4-y, 3/4+x, 1/4+z
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.75;
                         m.M[1] = 1.0;
         m.M[4] =-1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/4+y, 1/4-x, 3/4+z
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.75;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 3/4-y, 3/4-x, 3/4-z
         m = 0.0;    m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 3/4+y, 1/4+x, 1/4-z


      }// end of P_41_3_2

      else if(space_group_name == "I_41_3_2"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;  m.M[11] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;  m.M[3] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //1/2-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[7] =  0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, 1/2-y,  z

//5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;       m.M[7] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z, 1/2-x, y
         m = 0.0;       m.M[11]  = 0.5;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x, 1/2-y
         m = 0.0;       m.M[3] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //1/2-z, x,-y

//9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;        m.M[3] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //1/2-y, z,-x
         m = 0.0;        m.M[7] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, 1/2-z, x
         m = 0.0;        m.M[11] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, 1/2-x





         m = 0.0;  m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.25;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //3/4+x, 3/4-z, 1/4+y
         m = 0.0;  m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/4+x, 3/4+z, 3/4-y
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
         m.M[0] =-1.0;
                                        m.M[6] = -1.0;
                         m.M[9] =-1.0;                     operators.push_back(m);  //1/4-x, 1/4-z, 1/4-y
         m = 0.0;  m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.75;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //3/4+x, 1/4+z, 3/4+y


// 17
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4+z, 3/4+y, 3/4-x
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.75;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 3/4-z, 1/4+y, 3/4+x
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4-z, 1/4-y, 1/4-x
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.25;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 3/4+z, 3/4-y, 1/4+x


// 21
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.25;  m.M[11] = 0.75;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 3/4-y, 1/4+x, 3/4+z
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] =-1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 3/4+y, 3/4-x, 1/4+z
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] =-1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 1/4-y, 1/4-x, 1/4-z
         m = 0.0;    m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] =-1.0;   operators.push_back(m);  // 1/4+y, 3/4+x, 3/4-z

         // + (0.5, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }



      }// end of I_41_3_2

      else if(space_group_name == "P_-4_3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x




// 13
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //-x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // x, z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, -z, -y

// 17
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, -y, x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, y, -x


// 21
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, -x, -z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, -x, z




      }// end of P_-4_3_m

      else if(space_group_name == "F_-4_3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x




// 13
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //-x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // x, z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, -z, -y

// 17
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, -y, x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, y, -x


// 21
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, -x, -z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, -x, z


         // + (0.0, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

         // + (0.5, 0.0, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

         // + (0.5, 0.5, 0.0)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;  operators.push_back(m);
         }


      }// end of F_-4_3_m

      else if(space_group_name == "I_-4_3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x




// 13
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //-x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // x, z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, -z, -y

// 17
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, -y, x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, y, -x


// 21
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, -x, -z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, -x, z


         // + (0.5, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }


      }// end of I_-4_3_m

      else if(space_group_name == "P_-4_3_n"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z


         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y



         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x



// 13
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/2-x, 1/2+z, 1/2-y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/2-x, 1/2-z, 1/2+y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/2+x, 1/2+z, 1/2+y
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/2+x, 1/2-z, 1/2-y


// 17
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2-z, 1/2-y, 1/2+x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2+z, 1/2-y, 1/2-x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // 1/2+z, 1/2+y, 1/2+x
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2-z, 1/2+y, 1/2-x


// 21
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/2+y, 1/2-x, 1/2-z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/2-y, 1/2+x, 1/2-z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/2+y, 1/2+x, 1/2+z
         m = 0.0;  m.M[3] = m.M[7] = m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/2-y, 1/2-x, 1/2+z



      }// end of P_-4_3_n
      else if(space_group_name == "F_-4_3_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x




// 13
         m = 0.0;  m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //-x, z, 1/2-y
         m = 0.0;  m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, -z, 1/2+y
         m = 0.0;  m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // x, z, 1/2+y
         m = 0.0;  m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, -z, 1/2-y

// 17
         m = 0.0;  m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, -y, 1/2+x
         m = 0.0;  m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, -y, 1/2-x
         m = 0.0;  m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, y, 1/2+x
         m = 0.0;  m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, y, 1/2-x


// 21
         m = 0.0;  m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, -x, 1/2-z
         m = 0.0;  m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, x, 1/2-z
         m = 0.0;  m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, x, 1/2+z
	 m = 0.0;  m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, -x, 1/2+z


         // + (0.0, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

         // + (0.5, 0.0, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

         // + (0.5, 0.5, 0.0)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;  operators.push_back(m);
         }


      }// end of F_-4_3_c

      else if(space_group_name == "I_-4_3_d"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;  m.M[11] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, 1/2-z
         m = 0.0;
         m.M[0] = -1.0;  m.M[3] = 0.5;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //1/2-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;  m.M[7] = 0.5;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, 1/2-y,  z

//5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;       m.M[7] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z, 1/2-x, y
         m = 0.0;       m.M[11]  = 0.5;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x, 1/2-y
         m = 0.0;       m.M[3] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //1/2-z, x,-y

//9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;        m.M[3] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //1/2-y, z,-x
         m = 0.0;        m.M[7] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, 1/2-z, x
         m = 0.0;        m.M[11] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, 1/2-x




//13 
         m = 0.0;  m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.75;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/4-x, 1/4+z, 3/4-y
         m = 0.0;  m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/4-x, 3/4-z, 3/4+y
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
         m.M[0] =1.0;
                                        m.M[6] = 1.0;
                         m.M[9] =1.0;                     operators.push_back(m);  //1/4+x, 1/4+z, 1/4+y
         m = 0.0;  m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.25;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/4+x, 3/4-z, 1/4-y


// 17
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75;
                                        m.M[2] =-1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4-z, 3/4-y, 3/4+x
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.25;
                                        m.M[2] =1.0;
                        m.M[5] =  -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4+z, 3/4-y, 1/4-x
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
                                        m.M[2] =1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4+z, 1/4+y, 1/4+x
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.75;
                                        m.M[2] =-1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4-z, 1/4+y, 3/4-x


// 21
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/4+y, 3/4-x, 1/4-z
         m = 0.0;   m.M[3] = 0.75; m.M[7] = 0.75;  m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] =1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 3/4-y, 3/4+x, 1/4-z
         m = 0.0;   m.M[3] = 0.25; m.M[7] = 0.25;  m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] =1.0;
                                        m.M[10] =1.0;   operators.push_back(m);  // 1/4+y, 1/4+x, 1/4+z
         m = 0.0;    m.M[3] = 0.25; m.M[7] = 0.75;  m.M[11] = 0.75;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] =1.0;   operators.push_back(m);  // 1/4-y, 3/4-x, 3/4+z

         // + (0.5, 0.5, 0.5)
         sz = 24;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }



      }// end of I_-4_3_d

      else if(space_group_name == "P_m_-3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





// 13
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  // -x, -z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, z, y

// 17
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, y, -x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, -y, x


// 21
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, x, z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, -x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, -x,-z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, x, -z






// 25
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] =  -1.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //x,  -y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, y,  -z

// 29
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] = -1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  // -z, -x, -y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  //z,x, -y
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] =1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z, -x,y


// 33
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, -z, x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] =1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  //y, z, -x
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] =1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -y, z, x





// 37
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //-x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // x, z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, -z, -y

// 41
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, -y, x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, y, -x


// 45
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, -x, -z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, -x, z




      }// end of P_m_-3_m
      else if(space_group_name == "P_n_-3_n"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "P_m_-3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





// 13
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/2+x, 1/2-z, 1/2+y
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/2+x, 1/2+z, 1/2-y
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  // 1/2-x, 1/2-z, 1/2-y
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/2-x, 1/2+z, 1/2+y

// 17
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2+z, 1/2+y, 1/2-x
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2-z, 1/2+y, 1/2+x
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2-z, 1/2-y, 1/2-x
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2+z, 1/2-y, 1/2+x


// 21 
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/2-y, 1/2+x, 1/2+z
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/2+y, 1/2-x, 1/2+z
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/2-y, 1/2-x,1/2-z
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/2+y, 1/2+x, 1/2-z






// 25
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] =  -1.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //x,  -y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, y,  -z

// 29
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] = -1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  // -z, -x, -y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  //z,x, -y
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] =1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z, -x,y


// 33
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, -z, x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] =1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  //y, z, -x
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] =1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -y, z, x





// 37
         m = 0.0;  m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/2-x, 1/2+z, 1/2-y
         m = 0.0;  m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/2-x, 1/2-z, 1/2+y
         m = 0.0;  m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // 1/2+x, 1/2+z, 1/2+y
         m = 0.0;  m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/2+x, 1/2-z, 1/2-y

// 41
         m = 0.0;  m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2-z, 1/2-y, 1/2+x
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2+z, 1/2-y, 1/2-x
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2+z, 1/2+y, 1/2+x
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/2-z, 1/2+y, 1/2-x


// 45
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/2+y, 1/2-x, 1/2-z
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/2-y, 1/2+x, 1/2-z
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/2+y, 1/2+x, 1/2+z
         m = 0.0; m.M[3] = 0.5; m.M[7] = 0.5;  m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/2-y, 1/2-x, 1/2+z




      }// end of P_m_-3_n
      else if(space_group_name == "P_n_-3_m"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }
      else if(space_group_name == "F_m_-3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





// 13
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  // -x, -z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, z, y

// 17
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, y, -x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, -y, x


// 21
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, x, z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, -x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, -x,-z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, x, -z






// 25
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] =  -1.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //x,  -y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, y,  -z

// 29
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] = -1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  // -z, -x, -y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  //z,x, -y
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] =1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z, -x,y


// 33
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, -z, x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] =1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  //y, z, -x
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] =1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -y, z, x





// 37
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //-x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // x, z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, -z, -y

// 41
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, -y, x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, y, -x


// 45
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, -x, -z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, -x, z


         // + (0.0, 0.5, 0.5)
         sz = 48;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

         // + (0.5, 0.0, 0.5)
         sz = 48;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

         // + (0.5, 0.5, 0.0)
         sz = 48;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;  operators.push_back(m);
         }



      }// end of F_m_-3_m
      else if(space_group_name == "F_m_-3_c"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





// 13
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, -z, 1/2+y
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, z, 1/2-y
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  // -x, -z, 1/2-y
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, z, 1/2+y

// 17
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, y, 1/2-x
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, y, 1/2+x
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, -y, 1/2-x
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, -y, 1/2+x


// 21
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, -x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, -x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, x, 1/2-z






// 25
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] =  -1.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //x,  -y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, y,  -z

// 29
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] = -1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  // -z, -x, -y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  //z,x, -y
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] =1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z, -x,y


// 33
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, -z, x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] =1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  //y, z, -x
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] =1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -y, z, x





// 37
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //-x, z, 1/2-y
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, -z, 1/2+y
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // x, z, 1/2+y
         m = 0.0; m.M[11] = 0.5;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, -z, 1/2-y

// 41
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, -y, 1/2+x
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, -y, 1/2-x
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, y, 1/2+x
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, y, 1/2-x


// 45
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, -x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, x, 1/2-z
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, x, 1/2+z
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, -x, 1/2+z


         // + (0.0, 0.5, 0.5)
         sz = 48;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

         // + (0.5, 0.0, 0.5)
         sz = 48;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }

         // + (0.5, 0.5, 0.0)
         sz = 48;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5;  operators.push_back(m);
         }



      }// end of F_m_-3_c

      else if(space_group_name == "F_d_-3_m"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }

      else if(space_group_name == "F_d_-3_c"){

      std::cout<<"This space group is not programmed yet!"<<std::endl;
      // Two wariants
      }

      else if(space_group_name == "I_m_-3_m"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //-x,  y, -z
         m = 0.0;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, -y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,-x, y
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,-y
         m = 0.0;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //-y, z,-x
         m = 0.0;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, -z, x
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, -x





// 13
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  // -x, -z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, z, y

// 17
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, y, -x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, -y, x


// 21
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, x, z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, -x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, -x,-z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, x, -z






// 25
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] =  -1.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // -x,  -y,  -z
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //x,  -y, z
         m = 0.0;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, y,  -z

// 29
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] = -1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  // -z, -x, -y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  //z,x, -y
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] =1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,y
         m = 0.0;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //z, -x,y


// 33
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // y, -z, x
         m = 0.0;
                         m.M[1] =1.0;
                                        m.M[6] =1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  //y, z, -x
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] =1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -y, z, x





// 37
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //-x, z, -y
         m = 0.0;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //-x, -z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // x, z, y
         m = 0.0;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //x, -z, -y

// 41
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -z, -y, x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // z, -y, -x
         m = 0.0;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // z, y, x
         m = 0.0;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // -z, y, -x


// 45
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // y, -x, -z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // -y, x, -z
         m = 0.0;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // y, x, z
         m = 0.0;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // -y, -x, z


         // + (0.5, 0.5, 0.5)
         sz = 48;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }
        


      }// end of I_m_-3_m

      else if(space_group_name == "I_a_-3_d"){
         m = 0.0;
         m.M[0] =  1.0;
                         m.M[5] =  1.0;
                                        m.M[10] =  1.0;    operators.push_back(m);  // x,  y,  z
  	 m = 0.0;  m.M[11] = 0.5;
         m.M[0] =  1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  // x, -y, 1/2-z
         m = 0.0;  m.M[3] = 0.5;
         m.M[0] = -1.0;                                 
                         m.M[5] =  1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //1/2-x,  y, -z
         m = 0.0;  m.M[7] = 0.5;
         m.M[0] = -1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //-x, 1/2-y,  z

// 5
         m = 0.0;
                                        m.M[2] = 1.0;
         m.M[4] = 1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  // z, x, y
         m = 0.0;  m.M[7] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] =-1.0;
                        m.M[9] =  1.0;                     operators.push_back(m);  //-z,1/2-x, y
         m = 0.0;  m.M[11] = 0.5;
                                        m.M[2] = 1.0;
         m.M[4] =-1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  // z,-x,1/2-y
         m = 0.0;  m.M[3] = 0.5;
                                        m.M[2] =-1.0;
         m.M[4] = 1.0;
                        m.M[9] = -1.0;                     operators.push_back(m);  //1/2-z, x,-y


// 9
         m = 0.0;
                         m.M[1] = 1.0;
                                        m.M[6] = 1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  // y, z, x
         m = 0.0;  m.M[3] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  //1/2-y, z,-x
         m = 0.0;  m.M[7] = 0.5;
                         m.M[1] =-1.0;
                                        m.M[6] =-1.0;
         m.M[8] =  1.0;                                    operators.push_back(m);  //-y, 1/2-z, x
         m = 0.0;  m.M[11] = 0.5;
                         m.M[1] = 1.0;
                                        m.M[6] =-1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // y, -z, 1/2-x





// 13
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.25; m.M[11] = 0.75;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/4+x, 1/4-z, 3/4+y
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.75; m.M[11] = 0.75;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/4+x, 3/4+z, 3/4-y
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.25; m.M[11] = 0.25;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  // 1/4-x, 1/4-z, 1/4-y
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.25; m.M[11] = 0.75;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //3/4-x, 1/4+z, 3/4+y

// 17
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.75; m.M[11] = 0.75;
                                        m.M[2] = 1.0;
                        m.M[5] = 1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4+z, 3/4+y, 3/4-x
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.75; m.M[11] = 0.25;
                                        m.M[2] = -1.0;
                        m.M[5] = 1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4-z, 3/4+y, 1/4+x
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.25; m.M[11] = 0.25;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4-z, 1/4-y, 1/4-x
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.75; m.M[11] = 0.25;
                                        m.M[2] = 1.0;
                        m.M[5] =  -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 3/4+z, 3/4-y, 1/4+x


// 21
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.25; m.M[11] = 0.75;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 3/4-y, 1/4+x, 3/4+z
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.75; m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 3/4+y, 3/4-x, 1/4+z
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.25; m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 1/4-y, 1/4-x,1/4-z
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.25; m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 3/4+y, 1/4+x, 1/4-z






// 25
         m = 0.0;
         m.M[0] =  -1.0;
                         m.M[5] =  -1.0;
                                        m.M[10] =  -1.0;    operators.push_back(m);  // -x,  -y,  -z
 	 m = 0.0; m.M[11] = 0.5;
         m.M[0] =  -1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  // -x, y, 1/2+z
         m = 0.0; m.M[3] = 0.5;
         m.M[0] = 1.0;
                         m.M[5] = -1.0;
                                        m.M[10] = 1.0;    operators.push_back(m);  //1/2+x,  -y, z
         m = 0.0; m.M[7] = 0.5;
         m.M[0] = 1.0;
                         m.M[5] = 1.0;
                                        m.M[10] = -1.0;    operators.push_back(m);  //x, 1/2+y,  -z

// 29
         m = 0.0;
                                        m.M[2] = -1.0;
         m.M[4] = -1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  // -z, -x, -y
         m = 0.0; m.M[7] = 0.5;
                                        m.M[2] =1.0;
         m.M[4] =1.0;
                        m.M[9] =  -1.0;                     operators.push_back(m);  //z,1/2+x, -y
         m = 0.0; m.M[11] = 0.5;
                                        m.M[2] = -1.0;
         m.M[4] =1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  // -z,x,1/2+y
         m = 0.0; m.M[3] = 0.5;
                                        m.M[2] =1.0;
         m.M[4] = -1.0;
                        m.M[9] = 1.0;                     operators.push_back(m);  //1/2+z, -x,y


// 33
         m = 0.0;
                         m.M[1] = -1.0;
                                        m.M[6] = -1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  // -y, -z, -x
         m = 0.0; m.M[3] = 0.5;
                         m.M[1] =1.0;
                                        m.M[6] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/2+y, -z, x
         m = 0.0; m.M[7] = 0.5;
                         m.M[1] =1.0;
                                        m.M[6] =1.0;
         m.M[8] =  -1.0;                                    operators.push_back(m);  //y, 1/2+z, -x
         m = 0.0; m.M[11] = 0.5;
                         m.M[1] = -1.0;
                                        m.M[6] =1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // -y, z, 1/2+x





// 37
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.25; m.M[11] = 0.75;
         m.M[0] = -1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //1/4-x, 1/4+z, 3/4-y
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.75; m.M[11] = 0.75;
         m.M[0] = -1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  //1/4-x, 3/4-z, 3/4+y
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.25; m.M[11] = 0.25;
         m.M[0] = 1.0;
                                        m.M[6] = 1.0;
                         m.M[9] = 1.0;                     operators.push_back(m);  // 1/4+x, 1/4+z, 1/4+y
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.25; m.M[11] = 0.75;
         m.M[0] = 1.0;
                                        m.M[6] = -1.0;
                         m.M[9] = -1.0;                     operators.push_back(m);  //3/4+x, 1/4-z, 3/4-y

// 41
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.75; m.M[11] = 0.75;
                                        m.M[2] = -1.0;
                        m.M[5] = -1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4-z, 3/4-y, 3/4+x
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.75; m.M[11] = 0.25;
                                        m.M[2] = 1.0;
                        m.M[5] = -1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 1/4+z, 3/4-y, 1/4-x
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.25; m.M[11] = 0.25;
                                        m.M[2] = 1.0;
                        m.M[5] =  1.0;
         m.M[8] = 1.0;                                    operators.push_back(m);  // 1/4+z, 1/4+y, 1/4+x
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.75; m.M[11] = 0.25;
                                        m.M[2] = -1.0;
                        m.M[5] =  1.0;
         m.M[8] = -1.0;                                    operators.push_back(m);  // 3/4-z, 3/4+y, 1/4-x


// 45
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.25; m.M[11] = 0.75;
                         m.M[1] = 1.0;
         m.M[4] = -1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 3/4+y, 1/4-x, 3/4-z
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.75; m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] = 1.0;
                                        m.M[10] = -1.0;   operators.push_back(m);  // 3/4-y, 3/4+x, 1/4-z
         m = 0.0; m.M[3] = 0.25; m.M[7] = 0.25; m.M[11] = 0.25;
                         m.M[1] = 1.0;
         m.M[4] = 1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 1/4+y, 1/4+x, 1/4+z
         m = 0.0; m.M[3] = 0.75; m.M[7] = 0.25; m.M[11] = 0.25;
                         m.M[1] = -1.0;
         m.M[4] = -1.0;
                                        m.M[10] = 1.0;   operators.push_back(m);  // 3/4-y, 1/4-x, 1/4+z


         // + (0.5, 0.5, 0.5)
         sz = 48;
         for(i=0;i<sz;i++){
         m = operators[i];
         m.M[3] = 0.5; m.M[7] = 0.5; m.M[11] = 0.5; operators.push_back(m);
         }
        


      }// end of I_a_-3_d
   
    

      }// end of constructor
};


}//namespace libsymmetry
}//namespace liblibra




#endif // SPACE_GROUPS_H

