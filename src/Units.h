/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Units.h
  \brief The file defines the fundamental physical and conversion constants
    
*/

#ifndef UNITS_H
#define UNITS_H

/// liblibra namespace
namespace liblibra{


//================= Mathematical and numeric constants ====================
//    radian   conversion factor from radians to degrees
//    M_PI     numerical value of the geometric constant
#define radian 57.29577951308232088
#define radians_to_degrees     57.29577951308232088                  //   = 180/pi
#define rad_to_deg             57.29577951308232088
#define deg_to_rad             0.017453292519943295768261569865802


#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288419716939937510 // 50 digits
#endif

/*
 The rest of the digits
 Total number of digits = 1120
 Taken from: http://en.wikipedia.org/wiki/Pi

 58209749445923078164062862089986280348253421170679821480865132823066470938446095
 50582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610
 45432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036
 57595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798
 60943702770539217176293176752384674818467669405132000568127145263560827785771342757789609173637178721468440901224953430146549585371050792279
 68925892354201995611212902196086403441815981362977477130996051870721134999999837297804995105973173281609631859502445945534690830264252230825
 33446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066
 130019278766111959092164201989380952572010654858632788659361533818279682303019520353018529689957736225994138912497217752834791315155748572424541506959
*/


//======================== Conversion factors for Physical constants ==============

/*     avogadro    Avogadro's number (N) in particles/mole
*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/K/mole
*     gasconst    ideal gas constant (R) in kcal/mole/K
*     lightspd    speed of light in vacuum (c) in cm/ps
*     bohr        conversion from Bohrs to Angstroms
*     joule       conversion from calories to joules
*     evolt       conversion from Hartree to electron-volts
*     hartree     conversion from Hartree to kcal/mole
*     electric    conversion from electron**2/Ang to kcal/mole
*     debye       conversion from electron-Ang to Debyes
*     prescon     conversion from kcal/mole/Ang**3 to Atm
*/

// Length units
#define Angst 1.889725989                    // 1 Angstrom in atomic units
#define bohr       0.5291772083              // 1 Bohr = bohr Angstrom
#define nanometers_to_angstrem 10.0
#define bohrs_to_angstrem      0.5291772083
#define meters_to_angstrem     1.0E+10


// Time units
#define FS 41.34145          // 1 fs in atomic units of time
#define FS_TO_TAU 0.020454828280872953835407068780373 
#define ps_to_tau 20.45482828087295383540706878037      //   = 10*sqrt(4.184)
#define tau_to_ps 48.88821290839616117449108217105E-3    
#define tau_to_fs 48.88821290839616117449108217105      //   = 100/sqrt(4.184)
#define tau_to_s  48.88821290839616117449108217105E-15  


// Energy units
#define Ha 1.0                 // 1 Ha = 1 Ha 
#define Rydberg 0.5            // 1 Ry = 0.5 Ha;
#define eV 0.036749309         // 1 eV = 0.036.. Ha
#define joule      4.184       // 1 kcal = 4.184 kJ
#define evolt      27.2113834  // 1 Ha = 27.21.. eV
#define hartree    627.5094709 // 1 Ha = 627.5.. kcal/mol
#define cal_to_J   4.184       // 1 cal = 4.184 J
#define J_to_cal   0.23901     // 1 J = 0.239.. cal


// Pressure units --- Conversion from kcal/(mol*Angstrom^3) to Pa
#define pressure_to_Pa         1.6605387280149467216398197213547E+3  //  = 10000/avogadro
#define pressure_unit          694.76940380145370833410057141479E+2  //  1 kcal/(mol*Angstrom^3) = pressure_unit    bar
#define pressure_to_atm        68568.4                               //  1 kcal/(mol*Angstrom^3) = pressure_to_atm  atmospheres
#define atm_to_int             0.000014583977459004439362738520951342 // 1 atm = atm_to_int internal units


// Basic phisical constants
#define avogadro   6.02214199E+23
#define lightspd   2.99792458E-2
#define electric   332.05382
#define debye      4.8033324
#define prescon    6.85695E+4
#define boltzmann  1.9872065E-3        // in kcal/(mol*K)
//#define boltzmann  0.83143435
//#define gasconst   1.9872065E-3


// Electric dipole moments
#define Debye 0.393430307   // 1 Debye in atomic units [ e * Bohr]




}// namespace liblibra


#endif // UNITS_H
