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

#include "ForceField.h"

/// liblibra namespace
namespace liblibra{


namespace libforcefield{


string ForceField::mmff94_type(string elt,int geometry,string func_grp,int min_ring,int& coordination){

  //------- Transform element symbol to common format  --------------
  string type = "";
  int sz = elt.size();
  coordination = 2;

  if(sz==1){   type  = toupper(elt.c_str()[0]);  type += "_";    }
  else if(sz>=2){  type  = toupper(elt.c_str()[0]); type += tolower(elt.c_str()[1]);  }

  //-----------------Type determination ------------------------------

  if(type=="C_"){ if(geometry==1){type = "C%-"; coordination = 2; } // isonitrile carbon
                  else if(geometry==2){type = "CSP"; coordination = 1;} // linear, sp carbon
                  else if(geometry==3){ coordination = 3;  type = "C=C"; // default: sp2 carbon (vinilyc)
                    if(func_grp=="Carbonyl"||func_grp=="Aldehyde"||func_grp=="Haloformyl"||
                       func_grp=="Carbonate_ester"||func_grp=="Ester"||func_grp=="Carboxamide"||
                       func_grp=="Carboxyl"){ type = "C=O"; }
                    else if(func_grp=="Carboxylate"){ type = "CO2M"; }

                    if(min_ring==4){ type = "CE4R"; } // olefinic carbon in 4-membered rings
                    else if(min_ring==5){ type = "C5"; } // general carbon in 5-membered heteroaromatic ring
                    else if(min_ring==6){ type = "CB"; } // carbon as in benzen, pyrrole
                  }
                  else if((geometry>=4)||(geometry==0)){ type = "CR"; // default: alkyl carbon
                    if(min_ring==3){ type = "CR3R"; }
                    else if(min_ring==4){ type = "CR4R"; }
                  }                              
                } // Uncovered types are: CNN+(57), C5A(63), C5B(64), CIM+(80)
                
  else if(type=="N_"){ if(geometry==1){ type = "NSP"; coordination = 1; }   // nitrogen, triple bonded , linear, sp nitrogen
                       else if(geometry==2){ coordination = 3; type = "N=C"; // default: nitrogen in imines
                         if(func_grp=="Nitroso"){ type = "N=O"; }
                         else if(func_grp=="Isonitrile"){ type = "NR%"; }

                         if(min_ring==6){ type = "NPYD"; } // aromatic N, pyridine
                         else if(min_ring==5){ type = "N5"; }// general aromatic 5-ring nitrogen
                         else if(min_ring==4){ type = "CR4R"; }
                       } // trigonal-planar, sp2 with 2 connected atoms
                       else if(geometry==3){ type = "NR";  // default: amine N
                         if(func_grp=="Amide"){ type = "NC=O"; }
                         else if(func_grp=="Nitro"){ type = "NO2"; }

                         if(min_ring==5){ type = "NPYL"; } // aromatic N, pyrrole
                         else if(min_ring==6){ type = "NPD+";} // pyridinium ion
                       }
                       else if(geometry==4){ type = "NR+";} //quaternary nitrogen,sp3,positively charged, Sp3 N with four connected atoms
                       else if(geometry==0){   }
                     }  // Uncovered types are: N3OX(68),NAZT(47),NSO(48),=N=(53),NM(62),N5A(65),N5B(66),NC=C(40),NSO2(43),N+=C(54),NCN+(55)
                        // NGD+(56),N2OX(67),NPOX(69),NIM+(81),N5OX(82), N5M 

  else if(type=="O_"){ if(geometry==1){ coordination = 1; type = "O=C";} // default: O=C generic, Oxygen with one connected atom
                       else if(geometry==2){ type = "OR"; // default: Ether and ester oxygen
                         if(func_grp=="Water"){ type = "OH2"; }
                         else if(func_grp=="Hydroxyl") { type = "oh"; }// Oxygen in hydroxyl group                                   

                         if(min_ring==5){ type = "OFUR"; } // Aromatic O, Furan
                       }
                       else if(geometry==3){ type = "O+";} // Oxonium
                     } // Uncovered types are: O2CM(32),OM(35),O+=(51)

  else if(type=="S_"){ if(geometry==1){ type = "S=C"; } // doubly bonded to carbon
                       else if(geometry==2){ type = "S"; // default: thiol, sulfide
                         if(func_grp=="Sulfinyl"){ type = "=S=O"; }  // S in =S=O (sulfinyl) group
                         if(min_ring==5){ type = "STHI";}   // S in thiophene
                       }
                       else if(geometry==3){
                         if(func_grp=="Sulfinyl"){ type = "S=O";}   // S in -S=O (sulfoxide) group
                       }
                       else if(geometry==4){ type = "SO2";}  // S in sulfones
                     } // Uncovered types are: SM(72),

  else if(type=="P_"){ if(geometry==2){ type = "-P=C"; } // P doubly bonded to C
                       else if(geometry==3){ type = "P";}    // tricoordinate P
                       else if(geometry==4){ type = "PO4";}  // phosphodiester
                     }

  else if(type=="F_"){ if(geometry==1) {type = "F";}
                       else{            type = "F-"; }
                     }
  else if(type=="Cl"){ if(geometry==1) {type = "CL";}
                       else{            type = "CL-";}
                     }
  else if(type=="Br"){ if(geometry==1) {type = "BR";}
                       else{            type = "BR-";}
                     }
  else if(type=="I_"){ type = "I";}
  else if(type=="Si"){ type = "SI";}
  else if(type=="H_"){ }  // Skip on this step

  return type;
}


}// namespace libforcefield
}// namespace liblibra

