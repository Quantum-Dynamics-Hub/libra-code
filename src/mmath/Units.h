#define PI  3.1415926

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
*     convert     conversion from kcal to g*Ang**2/ps**2
*/

#define avogadro   6.02214199E+23
//#define boltzmann  0.83143435
//#define gasconst   1.9872065E-3
#define lightspd   2.99792458E-2
#define bohr       0.5291772083
#define joule      4.184
#define evolt      27.2113834
#define hartree    627.5094709
#define electric   332.05382
#define debye      4.8033324
#define prescon    6.85695E+4
//#define convert    4.184E+2


#define boltzmann  1.9872065E-3        // in kcal/(mol*K)

#define nanometers_to_angstrem 10
#define bohrs_to_angstrem      0.5291772083
#define meters_to_angstrem     1.0E+10

#define ps_to_tau              20.45482828087295383540706878037      //   = 10*sqrt(4.184)
#define tau_to_ps              48.88821290839616117449108217105E-3    
#define tau_to_fs              48.88821290839616117449108217105      //   = 100/sqrt(4.184)
#define tau_to_s               48.88821290839616117449108217105E-15  

#define radians_to_degrees     57.29577951308232088                  //   = 180/pi

#define cal_to_J               4.184
#define J_to_cal               0.23900573613766730401529636711281    //   =  1/4.184

#define deg_to_rad             0.017453292519943295768261569865802
#define rad_to_deg             57.29577951308232088

//--- Conversion from kcal/(mol*Angstrem^3) to Pa
#define pressure_to_Pa         1.6605387280149467216398197213547E+3  //  = 10000/avogadro
#define pressure_unit          694.76940380145370833410057141479E+2  //  1 kcal/(mol*Angstrem^3) = pressure_unit    bar
#define pressure_to_atm        68568.4                               //  1 kcal/(mol*Angstrom^3) = pressure_to_atm  atmospheres
#define atm_to_int             0.000014583977459004439362738520951342 // 1 atm = atm_to_int internal units

#define fs_to_tau              0.020454828280872953835407068780373 

